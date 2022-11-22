/* 
   This macro fits simulation dxHCAL plot to data.
   E.g. Config. File: sbs14-sbs70p-simu/conf_simu_data_fit.json
   -----
   P. Datta  Created  11-21-2022 [Older version: simu_data_fit_old.cpp]
*/

#include "TH1F.h"
#include "TFile.h"

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include "../../../include/gmn-ana.h"
#include "../../../dflay/src/JSONManager.cxx"

void CompareHisto (TH1F*, TH1F*);
TH1F* MakeHisto (std::string, double, int, double, double);

int simu_data_fit (const char *configfilename, std::string filebase="siout/simu_data_fit")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // reading in proper data and simu output files
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  int model = jmgr->GetValueFromKey<int>("model");
  TFile *fdata = new TFile(Form("../pdout/qelas_ana_data_sbs%d_sbs%dp_model%d_data.root", conf, sbsmag, model));
  TFile *fsimu = new TFile(Form("siout/qelas_ana_simu_sbs%d_sbs%dp_model%d_simu.root", conf, sbsmag, model));

  // defining the outputfile
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d.root", filebase.c_str(), conf, sbsmag, model);
  TFile *fout = new TFile(outFile,"RECREATE");

  // read in important histograms
  TH1F *h_dxHCAL_data; fdata->GetObject("h_dxHCAL_data",h_dxHCAL_data);
  TH1F *h_dxHCAL_bg; fdata->GetObject("h_dxHCAL_bg",h_dxHCAL_bg);
  TH1F *h_dxHCAL_simu; fsimu->GetObject("h_dxHCAL_simu",h_dxHCAL_simu);
  TH1F *h_dxHCAL_simu_p; fsimu->GetObject("h_dxHCAL_simu_p",h_dxHCAL_simu_p); h_dxHCAL_simu_p->Write();
  TH1F *h_dxHCAL_simu_n; fsimu->GetObject("h_dxHCAL_simu_n",h_dxHCAL_simu_n); h_dxHCAL_simu_n->Write();

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

  // define fit range
  // (E.g. leave the tails out. Need to implement radiative correction first.)
  vector<double> fit_range; jmgr->GetVectorFromKey<double>("fit_range", fit_range);
  int lbin = h_dxHCAL_data->FindBin(fit_range[0]);
  int hbin = h_dxHCAL_data->FindBin(fit_range[1]);

  // read in ranges to vary parameter R
  double R = 0;
  vector<double> R_lims; jmgr->GetVectorFromKey<double>("R_lims", R_lims);
  TH1F *h_n[int(R_lims[0])];
  TH1F *h_p[int(R_lims[0])];
  TH1F *h_comb_MC[int(R_lims[0])];

  vector<pair<double,double>> Rchi2;
  for (int iR=0; iR<int(R_lims[0]); iR++) {
    // construct R
    double Rwidth = (R_lims[2] - R_lims[1]) / R_lims[0];
    R = R_lims[1] + iR*Rwidth;
    
    // create histos
    h_comb_MC[iR] = MakeHisto("h_comb_MC_R", R, nbin, hmin, hmax);
    h_n[iR] = MakeHisto("h_n_R", R, nbin, hmin, hmax);
    h_p[iR] = MakeHisto("h_p_R", R, nbin, hmin, hmax);
    
    double norm = 1. / (h_dxHCAL_simu_p->Integral() + R*h_dxHCAL_simu_n->Integral());
    // Looping over bins to create combined simulation histo using AJRP method
    for (int ibin=lbin; ibin<hbin; ibin++) {
      double simu = norm*(h_dxHCAL_simu_p->GetBinContent(ibin) + R*h_dxHCAL_simu_n->GetBinContent(ibin));
      h_comb_MC[iR]->SetBinContent(ibin, simu);
      h_n[iR]->SetBinContent(ibin, norm*(R*h_dxHCAL_simu_n->GetBinContent(ibin)));
      h_p[iR]->SetBinContent(ibin, norm*(h_dxHCAL_simu_p->GetBinContent(ibin)));
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

    Rchi2.push_back(make_pair(R, chi2));
    // print out useful info
    cout << "R = " << R << " chi2 = " << chi2 << " lbin = " << lbin << " hbin = " << hbin << endl;
  }

  // canvas
  TCanvas *c1 = new TCanvas("c1","c1",1600,1200);
  c1->Divide(2,1);
  c1->cd(1);
  h_comb_MC[0]->Draw("p"); h_comb_MC[0]->SetMarkerStyle(20); h_comb_MC[0]->SetMarkerColor(2);
  // h_n->Draw("same p"); h_n->SetMarkerStyle(20); h_n->SetMarkerColor(4);
  // h_p->Draw("same p"); h_p->SetMarkerStyle(20); h_p->SetMarkerColor(8);
  h_dxHCAL_data->Draw("same");

  c1->cd(2);
  h_dxHCAL_simu->Draw(); h_dxHCAL_simu->SetMarkerStyle(20); h_dxHCAL_simu->SetMarkerColor(2);
  h_dxHCAL_data->Draw("same");


  // ofstream datafile1("chi2_vs_R_sbs50_.82T.csv", ios_base::app | ios_base::out);
  // datafile1 << R << "," << chi2 << endl;
  // ofstream datafile2("chi2_vs_R_sbs50_.82T_ap.csv", ios_base::app | ios_base::out);
  // datafile2 << R << "," << chi2_ap << endl;

  // h_dxHCAL_n->Scale(1./h_dxHCAL_n->Integral());
  // h_dxHCAL_p->Scale(1./h_dxHCAL_p->Integral());
  // c1->cd(2);
  // h_dxHCAL_data->Draw();
  // h_comb_MC_ap->Draw("same");
  // h_n->Draw("same");
  // h_p->Draw("same");

  // TString plotsfilename = outputfilename;
  // plotsfilename.ReplaceAll(".root",".pdf");
  
  // c1->Print(plotsfilename.Data(),"pdf");
  // plotsfilename.ReplaceAll(".pdf",".png");
  // c1->Print(plotsfilename.Data(),"png");

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
