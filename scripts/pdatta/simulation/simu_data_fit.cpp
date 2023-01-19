/* 
   This macro fits simulation dxHCAL plot to data.
   E.g. Config. File: sbs14-sbs70p-simu/conf_simu_data_fit.json
   -----
   P. Datta  Created  11-21-2022 [Older version: simu_data_fit_old.cpp]
*/

#include "TH1F.h"
#include "TFile.h"
#include "TLatex.h"

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include "../../../include/gmn-ana.h"
#include "../../../dflay/src/JSONManager.cxx"

void CompareHisto (TH1F*, TH1F*);
TH1F* MakeHisto (std::string, double, int, double, double);
TH1F* MakeHistoRB (std::string, double, double, int, double, double);

double fit_parabola (double *x, double *par) {
  return par[0] + par[1] * pow((x[0] - par[2]), 2); 
}

double fit_paraboloid (double *x, double *par) {
  return par[0] + par[1] * pow((x[0] - par[2]), 2) + par[3] * pow((x[1] - par[4]), 2); 
}

int simu_data_fit (const char *configfilename, std::string filebase="siout/simu_data_fit")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // reading in proper data and simu output files
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  int model = jmgr->GetValueFromKey<int>("model");
  int pass = jmgr->GetValueFromKey<int>("pass");
  TFile *fdata = new TFile(Form("../pdout/qelas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root", conf, sbsmag, model, pass));
  TFile *fsimu = new TFile(Form("siout/qelas_ana_simu_sbs%d_sbs%dp_model%d_simu.root", conf, sbsmag, model));

  // reading in mode of analysis and then defining the outputfile names
  bool fit_RB_simul = jmgr->GetValueFromKey<int>("fit_RB_simul"); // 1 => Simulataneous fit, 0 => Independent fit (i.e. fit B given Rmin)
  TString outFile, outtxtFile;
  if (!fit_RB_simul) {
    outFile = Form("%s_sbs%d_sbs%dp_model%d_RBindep.root", filebase.c_str(), conf, sbsmag, model);
    outtxtFile = Form("%s_sbs%d_sbs%dp_model%d_RBindep.csv", filebase.c_str(), conf, sbsmag, model);
  } else {
    outFile = Form("%s_sbs%d_sbs%dp_model%d_RBsimul.root", filebase.c_str(), conf, sbsmag, model);
    outtxtFile = Form("%s_sbs%d_sbs%dp_model%d_RBsimul.csv", filebase.c_str(), conf, sbsmag, model);
  }
  ofstream fit_data; fit_data.open(outtxtFile);
  fit_data << "R,B,chi2" << std::endl;  
  TFile *fout = new TFile(outFile,"RECREATE");

  // read in important histograms
  TH1F *h_dxHCAL_data; fdata->GetObject("h_dxHCAL_data",h_dxHCAL_data);
  TH1F *h_dxHCAL_bg; fdata->GetObject("h_dxHCAL_bg",h_dxHCAL_bg);
  TH1F *h_dxHCAL_simu; fsimu->GetObject("h_dxHCAL_simu",h_dxHCAL_simu);
  TH1F *h_dxHCAL_simu_p; fsimu->GetObject("h_dxHCAL_simu_p",h_dxHCAL_simu_p); 
  TH1F *h_dxHCAL_simu_n; fsimu->GetObject("h_dxHCAL_simu_n",h_dxHCAL_simu_n); 

  // sanity checks (nbin & ranges of histos need to be same for comparison)
  int nbin = h_dxHCAL_data->GetNbinsX();
  double hmin = h_dxHCAL_data->GetXaxis()->GetXmin();
  double hmax = h_dxHCAL_data->GetXaxis()->GetXmax();
  CompareHisto(h_dxHCAL_data, h_dxHCAL_bg);
  CompareHisto(h_dxHCAL_data, h_dxHCAL_simu);
  CompareHisto(h_dxHCAL_data, h_dxHCAL_simu_p);
  CompareHisto(h_dxHCAL_data, h_dxHCAL_simu_n);
  
  // normalize relevant histograms
  h_dxHCAL_data->Scale(1./h_dxHCAL_data->Integral()); h_dxHCAL_data->Write();
  h_dxHCAL_bg->Scale(1./h_dxHCAL_bg->Integral()); h_dxHCAL_bg->Write();
  h_dxHCAL_simu->Scale(1./h_dxHCAL_simu->Integral()); h_dxHCAL_simu->Write();
  // let's not scale the following two histograms yet
  h_dxHCAL_simu_p->Write();
  h_dxHCAL_simu_n->Write();

  // define fit range
  // (E.g. leave the tails out. Need to implement radiative correction first.)
  vector<double> fit_range; jmgr->GetVectorFromKey<double>("fit_range", fit_range);
  int lbin = h_dxHCAL_data->FindBin(fit_range[0]);
  int hbin = h_dxHCAL_data->FindBin(fit_range[1]);
  std::cout << Form("\nFit range: (%.2f,%.2f) | Corr. bins: (%d,%d)",fit_range[0],fit_range[1],lbin,hbin) << std::endl; 

  if (!fit_RB_simul) {

    // read in ranges to vary parameter R
    double R = 0;
    vector<double> R_lims; jmgr->GetVectorFromKey<double>("R_lims", R_lims);
    TH1F *h_n_R[int(R_lims[0])];
    TH1F *h_p_R[int(R_lims[0])];
    TH1F *h_comb_MC[int(R_lims[0])];

    std::cout << Form("\nVarying parameter R... [Range: (%.2f,%.2f)]",R_lims[1],R_lims[2]) << std::endl;
    // varying R
    double arrayR[int(R_lims[0])], arrayChi2[int(R_lims[0])];
    for (int iR=0; iR<int(R_lims[0]); iR++) {
      // construct R
      double Rwidth = (R_lims[2] - R_lims[1]) / R_lims[0];
      R = R_lims[1] + iR*Rwidth;
    
      // create histos
      h_comb_MC[iR] = MakeHisto("h_comb_MC_R", R, nbin, hmin, hmax);
      h_n_R[iR] = MakeHisto("h_n_R", R, nbin, hmin, hmax);
      h_p_R[iR] = MakeHisto("h_p_R", R, nbin, hmin, hmax);
    
      double norm = 1. / (h_dxHCAL_simu_p->Integral() + R*h_dxHCAL_simu_n->Integral());
      // Looping over bins to create combined simulation histo using AJRP method
      for (int ibin=lbin; ibin<hbin; ibin++) {
	double simu = norm*(h_dxHCAL_simu_p->GetBinContent(ibin) + R*h_dxHCAL_simu_n->GetBinContent(ibin));
	h_comb_MC[iR]->SetBinContent(ibin, simu);
	h_n_R[iR]->SetBinContent(ibin, norm*(R*h_dxHCAL_simu_n->GetBinContent(ibin)));
	h_p_R[iR]->SetBinContent(ibin, norm*(h_dxHCAL_simu_p->GetBinContent(ibin)));
      }

      // looping over bins again to calculate chi2
      double chi2 = 0.;
      for (int ibin=lbin; ibin<hbin; ibin++) {
	double simu = h_comb_MC[iR]->GetBinContent(ibin);
	double data = h_dxHCAL_data->GetBinContent(ibin);
	if (data>0) { 
	  double dataErr = sqrt(data)/sqrt(h_dxHCAL_data->GetEntries());
	  chi2 += (data-simu)*(data-simu) / ((dataErr)*(dataErr));
	}
      }

      // report R vs chi2
      arrayR[iR] = R; arrayChi2[iR] = chi2;
      // std::cout << Form("R = %.2f,  chi2 = %.2f",R,chi2) << std::endl;
      std::cout << Form("%.2f,%.2f",R,chi2) << std::endl;
    }
    std::cout << std::endl;

    // let's plot chi2 vs R
    TGraph *gchi2R = new TGraph(int(R_lims[0]), arrayR, arrayChi2);
    gchi2R->SetMarkerStyle(20); gchi2R->SetMarkerColor(4);
    gchi2R->SetTitle("#chi^{2} vs R");
    gchi2R->GetXaxis()->SetTitle("R");
    gchi2R->GetYaxis()->SetTitle("#chi^{2}");  gchi2R->GetYaxis()->SetMaxDigits(3);

    // let's fit chi2 vs R
    TF1 *chi2Rfn = new TF1("chi2Rfn",fit_parabola,R_lims[1],R_lims[2],3);
    chi2Rfn->SetNpx(500);
    // gchi2R->Fit(chi2Rfn, "R");
    // first try
    chi2Rfn->SetParameters(1,1,1);
    gchi2R->Fit("chi2Rfn","QR0");
    // second try
    chi2Rfn->SetParameters(16591,10113,1.3);
    gchi2R->Fit("chi2Rfn","QRV+","ep");

    // getting Rmin
    double Rmin = chi2Rfn->GetParameter(2);
    std::cout << " Rmin = " << Rmin << std::endl;

    // canvas
    TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
    c1->Divide(2,1);
    c1->cd(1);
    gchi2R->Draw("AP");
    // h_comb_MC[0]->Draw("p"); h_comb_MC[0]->SetMarkerStyle(20); h_comb_MC[0]->SetMarkerColor(2);
    // // h_n_R[0]->Draw("same p"); h_n_R->SetMarkerStyle(20); h_n_R->SetMarkerColor(4);
    // // h_p_R[0]->Draw("same p"); h_p_R->SetMarkerStyle(20); h_p_R->SetMarkerColor(8);
    // h_dxHCAL_data->Draw("same");

    // c1->cd(2);
    // h_dxHCAL_simu->Draw(); h_dxHCAL_simu->SetMarkerStyle(20); h_dxHCAL_simu->SetMarkerColor(2);
    // h_dxHCAL_data->Draw("same");


    // ofstream datafile1("chi2_vs_R_sbs50_.82T.csv", ios_base::app | ios_base::out);
    // datafile1 << R << "," << chi2 << endl;
    // ofstream datafile2("chi2_vs_R_sbs50_.82T_ap.csv", ios_base::app | ios_base::out);
    // datafile2 << R << "," << chi2_ap << endl;

    // ============================ Varying B ==============================
    double B = 0;
    vector<double> B_lims; jmgr->GetVectorFromKey<double>("B_lims", B_lims);
    TH1F *h_n_RB[int(B_lims[0])];
    TH1F *h_p_RB[int(B_lims[0])];
    TH1F *h_bg[int(B_lims[0])];
    TH1F *h_comb_MC_RB[int(B_lims[0])];

    std::cout << Form("\nVarying parameter B... [Range: (%.2f,%.2f)]",B_lims[1],B_lims[2]) << std::endl;
    // varying B
    double arrayB[int(B_lims[0])], arrayRBChi2[int(B_lims[0])];
    for (int iB=0; iB<int(B_lims[0]); iB++) {
      // construct B
      double Bwidth = (B_lims[2] - B_lims[1]) / B_lims[0];
      B = B_lims[1] + iB*Bwidth;
    
      // create histos
      h_comb_MC_RB[iB] = MakeHistoRB("h_comb_MC_RB", Rmin, B, nbin, hmin, hmax);
      h_n_RB[iB] = MakeHistoRB("h_n_RB", Rmin, B, nbin, hmin, hmax);
      h_p_RB[iB] = MakeHistoRB("h_p_RB", Rmin, B, nbin, hmin, hmax);
      h_bg[iB] = MakeHistoRB("h_bg", B, Rmin, nbin, hmin, hmax);
    
      double norm = 1. / (h_dxHCAL_simu_p->Integral() + Rmin*h_dxHCAL_simu_n->Integral() + B*h_dxHCAL_bg->Integral());
      // Looping over bins to create combined simulation histo using AJRP method
      for (int ibin=lbin; ibin<hbin; ibin++) {
	double simu = norm*(h_dxHCAL_simu_p->GetBinContent(ibin) + Rmin*h_dxHCAL_simu_n->GetBinContent(ibin) + B*h_dxHCAL_bg->GetBinContent(ibin));
	h_comb_MC_RB[iB]->SetBinContent(ibin, simu);
	h_n_RB[iB]->SetBinContent(ibin, norm*(Rmin*h_dxHCAL_simu_n->GetBinContent(ibin)));
	h_p_RB[iB]->SetBinContent(ibin, norm*(h_dxHCAL_simu_p->GetBinContent(ibin)));
	h_bg[iB]->SetBinContent(ibin, norm*(B*h_dxHCAL_bg->GetBinContent(ibin)));
      }

      // looping over bins again to calculate chi2
      double chi2 = 0.;
      for (int ibin=lbin; ibin<hbin; ibin++) {
	double simu = h_comb_MC_RB[iB]->GetBinContent(ibin);
	double data = h_dxHCAL_data->GetBinContent(ibin);
	if (data>0) { 
	  double dataErr = sqrt(data)/sqrt(h_dxHCAL_data->GetEntries());
	  chi2 += (data-simu)*(data-simu) / ((dataErr)*(dataErr));
	}
      }

      // report B vs chi2
      arrayB[iB] = B; arrayRBChi2[iB] = chi2;
      // std::cout << Form("%.4f,%.2f",B,chi2) << std::endl;
      fit_data << Form("%.4f,%.2f",B,chi2) << std::endl;
    }
    std::cout << std::endl;

    // let's plot chi2 vs B
    TGraph *gchi2B = new TGraph(int(B_lims[0]), arrayB, arrayRBChi2);
    gchi2B->SetMarkerStyle(20); gchi2B->SetMarkerColor(4);
    gchi2B->SetTitle("#chi^{2} vs B");
    gchi2B->GetXaxis()->SetTitle("B");
    gchi2B->GetYaxis()->SetTitle("#chi^{2}");  gchi2B->GetYaxis()->SetMaxDigits(3);

    // let's fit chi2 vs B
    TF1 *chi2Bfn = new TF1("chi2Bfn",fit_parabola,B_lims[1],B_lims[2],3);
    chi2Bfn->SetNpx(500);
    // gchi2B->Fit(chi2Bfn, "R");
    // first try
    chi2Bfn->SetParameters(1,1,1);
    gchi2B->Fit("chi2Bfn","QR0");
    // second try
    chi2Bfn->SetParameters(884,5000,0.032);
    gchi2B->Fit("chi2Bfn","QRV+","ep");

    c1->cd(2);
    gchi2B->Draw("AP");

    // getting Bmin
    double Bmin = chi2Bfn->GetParameter(2);
    std::cout << " Bmin = " << Bmin << std::endl;

    // **** Now that we have both Rmin and Bmin, it's time to draw the best fit histogram ****
    TH1F *h_comb_MC_RB_bfit = new TH1F("h_comb_MC_RB_bfit", Form("Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), nbin, hmin, hmax);
    TH1F *h_n_RB_bfit = new TH1F("h_n_RB_bfit", Form("n | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), nbin, hmin, hmax);
    TH1F *h_p_RB_bfit = new TH1F("h_p_RB_bfit", Form("p | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), nbin, hmin, hmax);
    TH1F *h_bg_bfit = new TH1F("h_bg_bfit", Form("bg | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), nbin, hmin, hmax);

    double norm = 1. / (h_dxHCAL_simu_p->Integral() + Rmin*h_dxHCAL_simu_n->Integral() + Bmin*h_dxHCAL_bg->Integral());
    // Looping over bins to create combined simulation histo using AJRP method
    for (int ibin=lbin; ibin<hbin; ibin++) {
      double simu = norm*(h_dxHCAL_simu_p->GetBinContent(ibin) + Rmin*h_dxHCAL_simu_n->GetBinContent(ibin) + Bmin*h_dxHCAL_bg->GetBinContent(ibin));
      h_comb_MC_RB_bfit->SetBinContent(ibin, simu);
      h_n_RB_bfit->SetBinContent(ibin, norm*(Rmin*h_dxHCAL_simu_n->GetBinContent(ibin)));
      h_p_RB_bfit->SetBinContent(ibin, norm*(h_dxHCAL_simu_p->GetBinContent(ibin)));
      h_bg_bfit->SetBinContent(ibin, norm*(Bmin*h_dxHCAL_bg->GetBinContent(ibin)));
    }
  
    // setting histogram styles
    h_comb_MC_RB_bfit->SetLineColor(kRed);
    h_n_RB_bfit->SetLineColor(kGreen);
    h_p_RB_bfit->SetLineColor(kBlue);
    h_bg_bfit->SetLineColor(kMagenta);
  }
  
  else {
    // ============================ Varying R & B simulatenously ==============================
    double R = 0;
    vector<double> R_lims; jmgr->GetVectorFromKey<double>("R_lims", R_lims);
    double B = 0;
    vector<double> B_lims; jmgr->GetVectorFromKey<double>("B_lims", B_lims);
    int Niter = int(R_lims[0])*int(B_lims[0]), iter = 0;
    TH1F *h_n_RB_simul[Niter];
    TH1F *h_p_RB_simul[Niter];
    TH1F *h_bg_simul[Niter];
    TH1F *h_comb_MC_RB_simul[Niter];

    // varying R
    double arrayR[Niter], arrayB[Niter], arrayChi2[Niter];
    for (int iR=0; iR<int(R_lims[0]); iR++) {
      // construct R
      double Rwidth = (R_lims[2] - R_lims[1]) / R_lims[0];
      R = R_lims[1] + iR*Rwidth;
      
      // varying B
      for (int iB=0; iB<int(B_lims[0]); iB++) {
	// construct B
	double Bwidth = (B_lims[2] - B_lims[1]) / B_lims[0];
	B = B_lims[1] + iB*Bwidth;
    
	// create histos
	h_comb_MC_RB_simul[iter] = MakeHistoRB("h_comb_MC_RB_simul", R, B, nbin, hmin, hmax);
	h_n_RB_simul[iter] = MakeHistoRB("h_n_RB_simul", R, B, nbin, hmin, hmax);
	h_p_RB_simul[iter] = MakeHistoRB("h_p_RB_simul", R, B, nbin, hmin, hmax);
	h_bg_simul[iter] = MakeHistoRB("h_bg_simul", R, B, nbin, hmin, hmax);
    
	double norm = 1. / (h_dxHCAL_simu_p->Integral() + R*h_dxHCAL_simu_n->Integral() + B*h_dxHCAL_bg->Integral());
	// Looping over bins to create combined simulation histo using AJRP method
	for (int ibin=lbin; ibin<hbin; ibin++) {
	  double simu = norm*(h_dxHCAL_simu_p->GetBinContent(ibin) + R*h_dxHCAL_simu_n->GetBinContent(ibin) + B*h_dxHCAL_bg->GetBinContent(ibin));
	  h_comb_MC_RB_simul[iter]->SetBinContent(ibin, simu);
	  h_n_RB_simul[iter]->SetBinContent(ibin, norm*(R*h_dxHCAL_simu_n->GetBinContent(ibin)));
	  h_p_RB_simul[iter]->SetBinContent(ibin, norm*(h_dxHCAL_simu_p->GetBinContent(ibin)));
	  h_bg_simul[iter]->SetBinContent(ibin, norm*(B*h_dxHCAL_bg->GetBinContent(ibin)));
	}

	// looping over bins again to calculate chi2
	double chi2 = 0.;
	for (int ibin=lbin; ibin<hbin; ibin++) {
	  double simu = h_comb_MC_RB_simul[iter]->GetBinContent(ibin);
	  double data = h_dxHCAL_data->GetBinContent(ibin);
	  if (data>0) { 
	    double dataErr = sqrt(data)/sqrt(h_dxHCAL_data->GetEntries());
	    chi2 += (data-simu)*(data-simu) / ((dataErr)*(dataErr));
	  }
	}

	// report R, B vs chi2
	arrayR[iter] = R; arrayB[iter] = B; arrayChi2[iter] = chi2;
	// std::cout << Form("%.2f,%.4f,%.2f",R,B,chi2) << std::endl;
	fit_data << Form("%.2f,%.4f,%.2f",R,B,chi2) << std::endl;

	iter++;
      } // for B

    } // for R
    std::cout << std::endl;

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800); c1->Divide(2,1);
    c1->cd(1);

    // let's plot chi2 vs R & B
    TGraph2D *gchi2RB = new TGraph2D(Niter, arrayR, arrayB, arrayChi2);
    gchi2RB->SetMarkerStyle(20); gchi2RB->SetMarkerColor(4);
    gchi2RB->SetTitle("#chi^{2} vs R & B");
    gchi2RB->GetXaxis()->SetTitle("R");
    gchi2RB->GetYaxis()->SetTitle("B");
    gchi2RB->GetZaxis()->SetTitle("#chi^{2}"); // gchi2RB->GetYaxis()->SetMaxDigits(3);
    gchi2RB->Draw("AP");

    // let's minimize chi2
    c1->cd(2);    
    TF2 *chi2RBfn = new TF2("chi2RBfn",fit_paraboloid,R_lims[1],R_lims[2],B_lims[1],B_lims[2],5);
    chi2RBfn->SetNpx(500);
    // first try
    chi2RBfn->SetParameters(1,1,1,1,1);
    gchi2RB->Fit("chi2RBfn","QR0");
    // second try
    chi2RBfn->SetParameters(800,1600,0.98,100000,0.03);
    gchi2RB->Fit("chi2RBfn","QRV+","ep");
    chi2RBfn->Draw();

    // getting Rmin & Bmin
    double Rmin = chi2RBfn->GetParameter(2), Bmin = chi2RBfn->GetParameter(4);
    std::cout << Form(" Rmin = %.2f, Bmin = %.4f", Rmin , Bmin)<< std::endl;

    // **** Now that we have both Rmin and Bmin, it's time to draw the best fit histogram ****
    TH1F *h_comb_MC_RB_bfit = new TH1F("h_comb_MC_RB_bfit", Form("Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), nbin, hmin, hmax);
    TH1F *h_n_RB_bfit = new TH1F("h_n_RB_bfit", Form("n | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), nbin, hmin, hmax);
    TH1F *h_p_RB_bfit = new TH1F("h_p_RB_bfit", Form("p | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), nbin, hmin, hmax);
    TH1F *h_bg_bfit = new TH1F("h_bg_bfit", Form("bg | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), nbin, hmin, hmax);

    double norm = 1. / (h_dxHCAL_simu_p->Integral() + Rmin*h_dxHCAL_simu_n->Integral() + Bmin*h_dxHCAL_bg->Integral());
    // Looping over bins to create combined simulation histo using AJRP method
    for (int ibin=lbin; ibin<hbin; ibin++) {
      double simu = norm*(h_dxHCAL_simu_p->GetBinContent(ibin) + Rmin*h_dxHCAL_simu_n->GetBinContent(ibin) + Bmin*h_dxHCAL_bg->GetBinContent(ibin));
      h_comb_MC_RB_bfit->SetBinContent(ibin, simu);
      h_n_RB_bfit->SetBinContent(ibin, norm*(Rmin*h_dxHCAL_simu_n->GetBinContent(ibin)));
      h_p_RB_bfit->SetBinContent(ibin, norm*(h_dxHCAL_simu_p->GetBinContent(ibin)));
      h_bg_bfit->SetBinContent(ibin, norm*(Bmin*h_dxHCAL_bg->GetBinContent(ibin)));
    }
  
    // setting histogram styles
    h_comb_MC_RB_bfit->SetLineColor(kRed);
    h_n_RB_bfit->SetLineColor(kGreen);
    h_p_RB_bfit->SetLineColor(kBlue);
    h_bg_bfit->SetLineColor(kMagenta);
    c1->Write();

  } // else

  fout->Write(); //fout->Close(); 
  delete jmgr;
  return 0;
}


// ---------------- Compare Histos for Sanity Checks ----------------
void CompareHisto (TH1F *h_base, TH1F *h_test) 
{
  if (h_test->GetNbinsX() != h_base->GetNbinsX()) {
    std::cout << "*!* nbin mismatch - " << h_test->GetName() << std::endl;
    throw;
  } 
  if (h_test->GetXaxis()->GetXmin() != h_base->GetXaxis()->GetXmin()) {
    std::cout << "*!* hmin mismatch - " << h_test->GetName() << std::endl;
    throw;
  }
  if (h_test->GetXaxis()->GetXmax() != h_base->GetXaxis()->GetXmax()) {
    std::cout << "*!* hmax mismatch - " << h_test->GetName() << std::endl;
    throw;
  }
}

// ---------------- Create generic histogram function ----------------
TH1F* MakeHisto (std::string pref, double param, int nbin, double min, double max)
{
  TH1F *h = new TH1F(TString::Format("%s_%.2f", pref.c_str(), param),
		     TString::Format("%s_%.2f", pref.c_str(), param), nbin, min, max);
  // h->SetStats(0);
  // h->SetLineWidth(2);
  // h->GetYaxis()->SetLabelSize(0.1);
  // h->GetYaxis()->SetLabelOffset(-0.17);
  // h->GetYaxis()->SetNdivisions(5);
  return h;
}
TH1F* MakeHistoRB (std::string pref, double R, double B, int nbin, double min, double max)
{
  TH1F *h = new TH1F(TString::Format("%s_%.2f_%.4f", pref.c_str(), R, B),
		     TString::Format("%s_%.2f_%.4f", pref.c_str(), R, B), nbin, min, max);
  return h;
}
