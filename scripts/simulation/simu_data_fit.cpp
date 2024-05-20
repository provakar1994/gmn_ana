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

#include "../../include/gmn_ana.h"
#include "../../dflay/src/JSONManager.cxx"

void CompareHisto (TH1F*, TH1F*);
TH1F* MakeHisto (std::string, double, int, double, double);
TH1F* MakeHistoRB (std::string, double, double, int, double, double);

double fit_parabola (double *x, double *par) {
  return par[0] + par[1] * pow((x[0] - par[2]), 2); 
}

double fit_paraboloid (double *x, double *par) {
  return par[0] + par[1] * pow((x[0] - par[2]), 2) + par[3] * pow((x[1] - par[4]), 2); 
}

int simu_data_fit (const char *configfilename, std::string filebase="siout/test_simu_data_fit")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // reading in proper data and simu output files
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  int model = jmgr->GetValueFromKey<int>("model");
  int pass = jmgr->GetValueFromKey<int>("pass");

  //############################//
  // Defining output file names //
  //############################//
  bool fit_RB_simul = jmgr->GetValueFromKey<int>("fit_RB_simul"); // 1 => Simulataneous fit, 0 => Independent fit (i.e. fit B given Rmin)
  bool fit_only_R = jmgr->GetValueFromKey<int>("fit_only_R"); // 1 => No bg included in the fit
  if (fit_only_R) fit_RB_simul = 0;
  TString outFile, outtxtFile, outtxtFileR;
  if (!fit_RB_simul) {
    if (!fit_only_R) {
      outFile = Form("%s_sbs%d_sbs%dp_model%d_RBindep.root", filebase.c_str(), conf, sbsmag, model);
      outtxtFileR = Form("%s_sbs%d_sbs%dp_model%d_RBindep_RvsChi2.csv", filebase.c_str(), conf, sbsmag, model);
      outtxtFile = Form("%s_sbs%d_sbs%dp_model%d_RBindep_BvsChi2.csv", filebase.c_str(), conf, sbsmag, model);
    } else {
      outFile = Form("%s_sbs%d_sbs%dp_model%d_R_no_bg.root", filebase.c_str(), conf, sbsmag, model);
      outtxtFileR = Form("%s_sbs%d_sbs%dp_model%d_R_no_bg.csv", filebase.c_str(), conf, sbsmag, model);
    }  
  } else {
    outFile = Form("%s_sbs%d_sbs%dp_model%d_RBsimul.root", filebase.c_str(), conf, sbsmag, model);
    outtxtFile = Form("%s_sbs%d_sbs%dp_model%d_RBsimul.csv", filebase.c_str(), conf, sbsmag, model);
  }
  ofstream fit_dataR; fit_dataR.open(outtxtFileR);
  fit_dataR << "R,chi2" << std::endl;
  ofstream fit_data; fit_data.open(outtxtFile);
  if (!fit_RB_simul)
    fit_data << "B,chi2" << std::endl; 
  else
    fit_data << "R,B,chi2" << std::endl; 
  TFile *fout = new TFile(outFile,"RECREATE");

  //####################################################//
  // Creating all the important histograms for analysis //
  //####################################################//
  // // read in important histograms (Old method, doesn't use RDataFrame)
  // TFile *fdata = new TFile(Form("../pdout/qelas_ana_data_sbs%d_sbs%dp_model%d_data.root", conf, sbsmag, model));
  // // TFile *fdata = new TFile(Form("../pdout/qelas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root", conf, sbsmag, model, pass));
  // TFile *fsimu = new TFile(Form("siout/qelas_ana_simu_sbs%d_sbs%dp_model%d.root", conf, sbsmag, model));
  // TH1F *h_dxHCAL_data; fdata->GetObject("h_dxHCAL_data",h_dxHCAL_data);
  // TH1F *h_dxHCAL_bg; fdata->GetObject("h_dxHCAL_bg",h_dxHCAL_bg);
  // TH1F *h_dxHCAL_simu; fsimu->GetObject("h_dxHCAL_simu",h_dxHCAL_simu);
  // TH1F *h_dxHCAL_simu_p; fsimu->GetObject("h_dxHCAL_simu_p",h_dxHCAL_simu_p); 
  // TH1F *h_dxH%CAL_simu_n; fsimu->GetObject("h_dxHCAL_simu_n",h_dxHCAL_simu_n); 

  // // sanity checks (nbin & ranges of histos need to be same for comparison)
  // int nbin = h_dxHCAL_data->GetNbinsX();
  // double hmin = h_dxHCAL_data->GetXaxis()->GetXmin();
  // double hmax = h_dxHCAL_data->GetXaxis()->GetXmax();
  // CompareHisto(h_dxHCAL_data, h_dxHCAL_bg);
  // CompareHisto(h_dxHCAL_data, h_dxHCAL_simu);
  // CompareHisto(h_dxHCAL_data, h_dxHCAL_simu_p);
  // CompareHisto(h_dxHCAL_data, h_dxHCAL_simu_n);

  /* 
     Switching to RDataFrame. Should fit perfectly for such analysis. Instead of reading the histograms from general data
     MC QE analysis ROOT files, I'll create them here. This will ensure the use of exactly the same cuts and histogram ranges
     between data and MC.
  */
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame data_rdf("Tout", Form("../pdout/qelas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",conf,sbsmag,model,pass));
  ROOT::RDataFrame simu_rdf("Tout", Form("siout/qelas_ana_simu_sbs%d_sbs%dp_model%d.root",conf,sbsmag,model));

  // Applying cuts
  std::string cuts_for_signal = jmgr->GetValueFromKey_str("cuts_for_signal");
  std::string cuts_for_bg = jmgr->GetValueFromKey_str("cuts_for_bg");
  auto data_rdf_filtered = data_rdf.Filter(cuts_for_signal);
  auto simu_rdf_filtered = simu_rdf.Filter(cuts_for_signal);
  auto bg_rdf_filtered = data_rdf.Filter(cuts_for_bg);

  // Creating important histograms
  vector<double> h_dx; jmgr->GetVectorFromKey<double>("h_dx", h_dx);
  auto h_dxHCAL_data = data_rdf_filtered.Histo1D({"h_dxHCAL_data", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx");
  auto h_dxHCAL_simu = simu_rdf_filtered.Histo1D({"h_dxHCAL_simu", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx", "weight");
  auto h_dxHCAL_simu_p = simu_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_simu_p", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx", "weight");
  auto h_dxHCAL_simu_n = simu_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_simu_n", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx", "weight");
  auto h_dxHCAL_bg = bg_rdf_filtered.Histo1D({"h_dxHCAL_bg", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx");
  
  // customize data histogram
  h_dxHCAL_data->SetMarkerStyle(20); h_dxHCAL_data->SetMarkerStyle(kBlack); 
  h_dxHCAL_data->SetMarkerSize(0.6); h_dxHCAL_data->SetLineColor(kBlack); 
  // normalize relevant histograms
  h_dxHCAL_data->Scale(1./h_dxHCAL_data->Integral()); h_dxHCAL_data->Write();
  h_dxHCAL_bg->Scale(1./h_dxHCAL_bg->Integral()); if(!fit_only_R) h_dxHCAL_bg->Write();
  h_dxHCAL_simu->Scale(1./h_dxHCAL_simu->Integral()); h_dxHCAL_simu->Write();
  // let's not scale the following two histograms yet
  h_dxHCAL_simu_p->Write();
  h_dxHCAL_simu_n->Write();

  // define fit range
  // (E.g. leave the tails out. Need to implement radiative correction first.)
  bool write_all_fit_histos = jmgr->GetValueFromKey<int>("write_all_fit_histos");
  vector<double> fit_range; jmgr->GetVectorFromKey<double>("fit_range", fit_range);
  int lbin = h_dxHCAL_data->FindBin(fit_range[0]);
  int hbin = h_dxHCAL_data->FindBin(fit_range[1]);
  std::cout << Form("\nFit range: (%.2f,%.2f) | Corr. bins: (%d,%d)",fit_range[0],fit_range[1],lbin,hbin) << std::endl; 

  //###################//
  // Starting analysis //
  //###################//
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
      h_comb_MC[iR] = MakeHisto("h_comb_MC_R", R, int(h_dx[0]), h_dx[1], h_dx[2]);
      h_n_R[iR] = MakeHisto("h_n_R", R, int(h_dx[0]), h_dx[1], h_dx[2]);
      h_p_R[iR] = MakeHisto("h_p_R", R, int(h_dx[0]), h_dx[1], h_dx[2]);
    
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
	  // The calculation of data error is a little tricky:
	  // 1. The data histogram is normalized here by 1/N, so cosidering Poission dist.
	  //    the stat error of each bin is sqrt(n_i) x 1/N, Here n_i is the content of ith bin
	  // 2. However, the variable "data" here is actually n_i/N due to scaling.
	  // 3. Hence, the correct stat error of ith bin should be sqrt(data)/sqrt(N) = sqrt(n_i)/N
	  double dataErr = sqrt(data)/sqrt(h_dxHCAL_data->GetEntries());
	  chi2 += (data-simu)*(data-simu) / ((dataErr)*(dataErr));
	}
      }

      // report R vs chi2
      arrayR[iR] = R; arrayChi2[iR] = chi2;
      fit_dataR << Form("%.2f,%.2f",R,chi2) << std::endl;
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
    vector<double> fpR; jmgr->GetVectorFromKey<double>("fit_params_R", fpR);
    TF1 *chi2Rfn = new TF1("chi2Rfn",fit_parabola,R_lims[1],R_lims[2],3);
    chi2Rfn->SetNpx(500);
    // gchi2R->Fit(chi2Rfn, "R");
    // first try
    chi2Rfn->SetParameters(1,1,1);
    gchi2R->Fit("chi2Rfn","QR0");
    // second try
    chi2Rfn->SetParameters(fpR[0],fpR[1],fpR[2]);
    gchi2R->Fit("chi2Rfn","QRV+","ep");

    // getting Rmin
    double Rmin = chi2Rfn->GetParameter(2);
    std::cout << " Rmin = " << Rmin << std::endl;

    // canvas
    TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
    c1->Divide(2,2);
    c1->cd(1);
    gchi2R->Draw("AP");

    if (fit_only_R) {
      // **** Fitting only with R (No background) 
      // **** Now that we have Rmin, it's time to draw the best fit histogram ****
      TH1F *h_comb_MC_R_bfit = new TH1F("h_comb_MC_R_bfit", Form("Rmin = %0.2f", Rmin), int(h_dx[0]), h_dx[1], h_dx[2]);
      TH1F *h_n_R_bfit = new TH1F("h_n_R_bfit", Form("n | Rmin = %0.2f", Rmin), int(h_dx[0]), h_dx[1], h_dx[2]);
      TH1F *h_p_R_bfit = new TH1F("h_p_R_bfit", Form("p | Rmin = %0.2f", Rmin), int(h_dx[0]), h_dx[1], h_dx[2]);

      double norm = 1. / (h_dxHCAL_simu_p->Integral() + Rmin*h_dxHCAL_simu_n->Integral());
      // Looping over bins to create combined simulation histo using AJRP method
      for (int ibin=lbin; ibin<hbin; ibin++) {
	double simu = norm*(h_dxHCAL_simu_p->GetBinContent(ibin) + Rmin*h_dxHCAL_simu_n->GetBinContent(ibin));
	h_comb_MC_R_bfit->SetBinContent(ibin, simu);
	h_n_R_bfit->SetBinContent(ibin, norm*(Rmin*h_dxHCAL_simu_n->GetBinContent(ibin)));
	h_p_R_bfit->SetBinContent(ibin, norm*(h_dxHCAL_simu_p->GetBinContent(ibin)));
      }
  
      // setting histogram styles
      h_comb_MC_R_bfit->SetLineColor(kRed);
      h_n_R_bfit->SetLineColor(kGreen);
      h_p_R_bfit->SetLineColor(kBlue);

      c1->cd(2);
      h_comb_MC_R_bfit->Draw();
      h_dxHCAL_data->DrawClone("same");

      c1->cd(3);
      h_comb_MC_R_bfit->Draw();
      h_n_R_bfit->Draw("same");
      h_p_R_bfit->Draw("same");
      h_dxHCAL_data->DrawClone("same");

      if (!write_all_fit_histos) {
	h_comb_MC_R_bfit->Write();
	h_n_R_bfit->Write();
	h_p_R_bfit->Write();
      }
      c1->Write();
 
    } else {

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
	h_comb_MC_RB[iB] = MakeHistoRB("h_comb_MC_RB", Rmin, B, int(h_dx[0]), h_dx[1], h_dx[2]);
	h_n_RB[iB] = MakeHistoRB("h_n_RB", Rmin, B, int(h_dx[0]), h_dx[1], h_dx[2]);
	h_p_RB[iB] = MakeHistoRB("h_p_RB", Rmin, B, int(h_dx[0]), h_dx[1], h_dx[2]);
	h_bg[iB] = MakeHistoRB("h_bg", B, Rmin, int(h_dx[0]), h_dx[1], h_dx[2]);
    
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
	    // The calculation of data error is a little tricky:
	    // 1. The data histogram is normalized here by 1/N, so cosidering Poission dist.
	    //    the stat error of each bin is sqrt(n_i) x 1/N, Here n_i is the content of ith bin
	    // 2. However, the variable "data" here is actually n_i/N due to scaling.
	    // 3. Hence, the correct stat error of ith bin should be sqrt(data)/sqrt(N) = sqrt(n_i)/N
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
      vector<double> fpB; jmgr->GetVectorFromKey<double>("fit_params_B", fpB);
      TF1 *chi2Bfn = new TF1("chi2Bfn",fit_parabola,B_lims[1],B_lims[2],3);
      chi2Bfn->SetNpx(500);
      // gchi2B->Fit(chi2Bfn, "R");
      // first try
      chi2Bfn->SetParameters(1,1,1);
      gchi2B->Fit("chi2Bfn","QR0");
      // second try
      chi2Bfn->SetParameters(fpB[0],fpB[1],fpB[2]);
      gchi2B->Fit("chi2Bfn","QRV+","ep");

      c1->cd(2);
      gchi2B->Draw("AP");

      // getting Bmin
      double Bmin = chi2Bfn->GetParameter(2);
      std::cout << " Bmin = " << Bmin << std::endl;

      // **** Now that we have both Rmin and Bmin, it's time to draw the best fit histogram ****
      TH1F *h_comb_MC_RB_bfit = new TH1F("h_comb_MC_RB_bfit", Form("Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), int(h_dx[0]), h_dx[1], h_dx[2]);
      TH1F *h_n_RB_bfit = new TH1F("h_n_RB_bfit", Form("n | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), int(h_dx[0]), h_dx[1], h_dx[2]);
      TH1F *h_p_RB_bfit = new TH1F("h_p_RB_bfit", Form("p | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), int(h_dx[0]), h_dx[1], h_dx[2]);
      TH1F *h_bg_bfit = new TH1F("h_bg_bfit", Form("bg | Rmin = %0.2f, Bmin = %0.4f", Rmin, Bmin), int(h_dx[0]), h_dx[1], h_dx[2]);

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

      c1->cd(3);
      h_comb_MC_RB_bfit->Draw();
      h_dxHCAL_data->DrawClone("same");

      c1->cd(4);
      h_comb_MC_RB_bfit->Draw();
      h_n_RB_bfit->Draw("same");
      h_p_RB_bfit->Draw("same");
      h_bg_bfit->Draw("same");
      h_dxHCAL_data->DrawClone("same");

      if (!write_all_fit_histos) {
	h_comb_MC_RB_bfit->Write();
	h_n_RB_bfit->Write();
	h_p_RB_bfit->Write();
	h_bg_bfit->Write();
      }
      c1->Write();

    } // fit_only_R = False
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
	h_comb_MC_RB_simul[iter] = MakeHistoRB("h_comb_MC_RB_simul", R, B, int(h_dx[0]), h_dx[1], h_dx[2]);
	h_n_RB_simul[iter] = MakeHistoRB("h_n_RB_simul", R, B, int(h_dx[0]), h_dx[1], h_dx[2]);
	h_p_RB_simul[iter] = MakeHistoRB("h_p_RB_simul", R, B, int(h_dx[0]), h_dx[1], h_dx[2]);
	h_bg_simul[iter] = MakeHistoRB("h_bg_simul", R, B, int(h_dx[0]), h_dx[1], h_dx[2]);
    
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
	    // The calculation of data error is a little tricky:
	    // 1. The data histogram is normalized here by 1/N, so cosidering Poission dist.
	    //    the stat error of each bin is sqrt(n_i) x 1/N, Here n_i is the content of ith bin
	    // 2. However, the variable "data" here is actually n_i/N due to scaling.
	    // 3. Hence, the correct stat error of ith bin should be sqrt(data)/sqrt(N) = sqrt(n_i)/N
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

    TCanvas *c1 = new TCanvas("c1", "c1", 1200, 800); c1->Divide(2,2);
    c1->cd(1);

    // let's plot chi2 vs R & B
    TGraph2D *gchi2RB = new TGraph2D(Niter, arrayR, arrayB, arrayChi2);
    gchi2RB->SetMarkerStyle(20); gchi2RB->SetMarkerColor(4);
    gchi2RB->SetTitle("#chi^{2} vs R & B");
    gchi2RB->GetXaxis()->SetTitle("R");
    gchi2RB->GetYaxis()->SetTitle("B");
    gchi2RB->GetZaxis()->SetTitle("#chi^{2}"); // gchi2RB->GetYaxis()->SetMaxDigits(3);
    gchi2RB->Draw("AP");

    c1->cd(2);    
    // let's minimize chi2
    vector<double> fpRB; jmgr->GetVectorFromKey<double>("fit_params_RB", fpRB);
    TF2 *chi2RBfn = new TF2("chi2RBfn",fit_paraboloid,R_lims[1],R_lims[2],B_lims[1],B_lims[2],5);
    chi2RBfn->SetNpx(500);
    // first try
    chi2RBfn->SetParameters(1,1,1,1,1);
    gchi2RB->Fit("chi2RBfn","QR0");
    // second try
    chi2RBfn->SetParameters(fpRB[0],fpRB[1],fpRB[2],fpRB[3],fpRB[4]);
    gchi2RB->Fit("chi2RBfn","QRV+","ep");
    chi2RBfn->Draw();

    // getting Rmin & Bmin
    double Rmin = chi2RBfn->GetParameter(2), Bmin = chi2RBfn->GetParameter(4), Chi2 = chi2RBfn->Eval(Rmin,Bmin);
    std::cout << Form(" Rmin = %.2f, Bmin = %.4f, Chi2 = %.2f", Rmin , Bmin, Chi2) << std::endl;

    // **** Now that we have both Rmin and Bmin, it's time to draw the best fit histogram ****
    TH1F *h_comb_MC_RB_bfit = new TH1F("h_comb_MC_RB_bfit", Form("Rmin = %0.2f, Bmin = %0.4f, #chi^{2} = %0.2f",Rmin,Bmin,Chi2), int(h_dx[0]), h_dx[1], h_dx[2]);
    TH1F *h_n_RB_bfit = new TH1F("h_n_RB_bfit", Form("n | Rmin = %0.2f, Bmin = %0.4f, #chi^{2} = %0.2f",Rmin,Bmin,Chi2), int(h_dx[0]), h_dx[1], h_dx[2]);
    TH1F *h_p_RB_bfit = new TH1F("h_p_RB_bfit", Form("p | Rmin = %0.2f, Bmin = %0.4f, #chi^{2} = %0.2f",Rmin,Bmin,Chi2), int(h_dx[0]), h_dx[1], h_dx[2]);
    TH1F *h_bg_bfit = new TH1F("h_bg_bfit", Form("bg | Rmin = %0.2f, Bmin = %0.4f, #chi^{2} = %0.2f",Rmin,Bmin,Chi2), int(h_dx[0]), h_dx[1], h_dx[2]);

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

    c1->cd(3);
    h_comb_MC_RB_bfit->Draw();
    h_dxHCAL_data->DrawClone("same");

    c1->cd(4);
    h_comb_MC_RB_bfit->Draw();
    h_n_RB_bfit->Draw("same");
    h_p_RB_bfit->Draw("same");
    h_bg_bfit->Draw("same");
    h_dxHCAL_data->DrawClone("same");

    if (!write_all_fit_histos) {
      h_comb_MC_RB_bfit->Write();
      h_n_RB_bfit->Write();
      h_p_RB_bfit->Write();
      h_bg_bfit->Write();
    }
    c1->Write();

  } // else

  if (write_all_fit_histos) fout->Write(); //fout->Close(); 
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
