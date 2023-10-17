/*
  This macro will fit W2 distribution in various different ways.
  Plan:
  1. MC signal, no bg
  2. MC signal + poly bg
  3. Signal from data + poly bg
  -------
  P. Datta Created 05-02-2023
*/

#include "TH1F.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include "../include/gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

void customize_W2(TH1F* h)
{
  h->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
}

void customize_ht(TH1F* h) 
{
  h->SetLineColor(kBlue);
  h->SetLineWidth(2);
}

void customize_hs(TH1F* h) 
{
  h->SetMarkerColor(kBlack);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(20);
  h->SetLineColor(kBlack);
}

void customize_hbg(TH1F* h) 
{
  h->SetMarkerColor(kRed);
  h->SetMarkerStyle(29);
  h->SetLineColor(kRed);
}

int fit_W2 (const char *configfilename, 
	    bool is_elastic = 1, // 1=>Yes, 0=>QE
	    std::string filebase="pdout/test_fit_W2") 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  char const * key = is_elastic ? "elas" : "qelas";

  // reading in proper data and simu output files
  int conf = jmgr->GetValueFromSubKey<int>(key,"SBS_config");
  int sbsmag = jmgr->GetValueFromSubKey<int>(key,"SBS_magnet_percent");
  int model = jmgr->GetValueFromSubKey<int>(key,"model");
  int pass = jmgr->GetValueFromSubKey<int>(key,"pass");
  char const * process = is_elastic ? "elas" : "qelas";
  std::string gen = jmgr->GetValueFromSubKey_str(key,"generator");
  std::string dfprefix = jmgr->GetValueFromSubKey_str(key,"data_file_prefix");
  std::string sfprefix = jmgr->GetValueFromSubKey_str(key,"simu_file_prefix");
  dfprefix = dfprefix.empty() ? "" : dfprefix + "_";
  sfprefix = sfprefix.empty() ? "" : sfprefix + "_";
  // char const * dfp = dfprefix.empty() ? "" : (dfprefix + "_").c_str();
  // char const * sfp = sfprefix.empty() ? "" : (sfprefix + "_").c_str();

  // defining output files
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d.root",filebase.c_str(),conf,sbsmag,model);

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/%s%s_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfprefix.c_str(),process,conf,sbsmag,model,pass));
  ROOT::RDataFrame simu_rdf("Tout",Form("simulation/siout/%s%s_ana_%s_sbs%d_sbs%dp_model%d.root",sfprefix.c_str(),process,gen.c_str(),conf,sbsmag,model));

  // Applying cuts
  std::string cuts_for_signal_data = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_data");
  std::string cuts_for_signal_simu = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_simu");
  std::string cuts_for_bg = jmgr->GetValueFromSubKey_str(key,"cuts_for_bg");
  std::string coinT_cut = jmgr->GetValueFromSubKey_str(key,"coinT_cut");
  int Opoly = jmgr->GetValueFromSubKey<int>(key,"Order_of_poly_bg");
  auto data_rdf_filtered = data_rdf.Filter(cuts_for_signal_data);
  auto simu_rdf_filtered_1 = simu_rdf.Filter(cuts_for_signal_simu);
  double offset = jmgr->GetValueFromSubKey<double>(key,"W2_peak_offset_for_MC");
  std::string W2_shifted = "W2+" + std::to_string(offset);
  auto simu_rdf_filtered = simu_rdf_filtered_1.Define("W2_shifted",W2_shifted.c_str());
  auto bg_rdf_filtered = data_rdf.Filter(cuts_for_bg);

  // Creating important histograms
  vector<double> h_W2; jmgr->GetVectorFromSubKey<double>(key,"h_W2",h_W2);
  TH1F *h_W2_data_raw = (TH1F*)data_rdf.Histo1D({"h_W2_data_raw","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2")->Clone();
  TH1F *h_W2_data = (TH1F*)data_rdf_filtered.Histo1D({"h_W2_data","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2")->Clone();
  TH1F *h_W2_data_CT = (TH1F*)data_rdf_filtered.Filter(coinT_cut.c_str()).Histo1D({"h_W2_data_CT","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2")->Clone();
  TH1F *h_W2_simu = (TH1F*)simu_rdf_filtered.Histo1D({"h_W2_simu","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2_shifted","weight")->Clone();
  TH1F *h_W2_simu_p = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_W2_simu_p","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2_shifted","weight")->Clone();
  TH1F *h_W2_simu_n = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_W2_simu_n","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2_shifted","weight")->Clone();
  TH1F *h_W2_bg = (TH1F*)bg_rdf_filtered.Histo1D({"h_W2_bg","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2")->Clone();

  // Fits
  vector<double> W2_fit_range; jmgr->GetVectorFromSubKey<double>(key,"W2_fit_range",W2_fit_range);

  // Canvas 0 : Fitting data/MC w/o any background
  TCanvas *c0 = util_pd::TC("c0",1,1);
  c0->cd(); gStyle->SetOptFit(1);
  vector<TH1F*> ho;
  TF1 *f0 = fit::fit_1hs_nbg_THI(W2_fit_range,
  				 h_W2_data_CT,h_W2_simu,
  				 ho);
  ho[0]->Draw(); customize_ht(ho[0]); customize_W2(ho[0]);
  ho[1]->Draw("same"); customize_hs(ho[1]); customize_W2(ho[1]);
  TLegend *l0=new TLegend(0.10,0.77,0.38,0.9);
  l0->SetTextFont(42);
  l0->AddEntry(ho[0],"Data","l");
  l0->AddEntry(f0,"Global Fit","l");
  l0->AddEntry(ho[1],"Signal (from MC)","lep");
  l0->Draw();
  // --- 

  // Canvas 1 : Fitting data/MC w/ polynomial background
  TCanvas *c1 = util_pd::TC("c1",1,1);
  c1->cd(); gStyle->SetOptFit(1);
  vector<TH1F*> ho1;
  TF1 *f1 = fit::fit_1hs_1pbg_THI(W2_fit_range,
  				  h_W2_data_raw,h_W2_simu,Opoly,
  				  ho1);
  ho1[0]->Draw(); customize_ht(ho1[0]); customize_W2(ho1[0]);
  ho1[1]->Draw("same"); customize_hs(ho1[1]); customize_W2(ho1[1]);
  ho1[2]->Draw("same"); customize_hbg(ho1[2]); 
  // drawing the polynomial background as well
  FitFn *ffn = new FitFn(Opoly);
  TF1* bg1 = new TF1("bg1",ffn,&FitFn::ffn_poly,W2_fit_range[0],W2_fit_range[1],Opoly+1);
  bg1->SetNpx(500);
  bg1->SetParameters(&fit::GetFitParams(f1)[1]);
  bg1->SetLineColor(kGreen+2);
  bg1->Draw("same");
  TLegend *l1 = new TLegend(0.10,0.77,0.38,0.9);
  l1->SetTextFont(42);
  l1->AddEntry(ho1[0],"Data","l");
  l1->AddEntry(f1,"Global Fit (MC + Data bg.)","l");
  l1->AddEntry(ho1[1],"Signal (from MC)","lep");
  l1->AddEntry(ho1[2],Form("Bg. (%d^{th} order poly.)",Opoly),"lep");
  l1->Draw();
  std::cout << Form("\nTotal # elastics (from fitted MC signal): %d\n",
		    int(ho1[1]->Integral(ho1[1]->FindBin(W2_fit_range[0]),ho1[1]->FindBin(W2_fit_range[1]))));

  double bgcount = (int)bg1->Integral(W2_fit_range[0],W2_fit_range[1])/h_W2_data->GetBinWidth(1);
  double totcount = h_W2_data_raw->Integral(h_W2_data_raw->FindBin(W2_fit_range[0]),h_W2_data_raw->FindBin(W2_fit_range[1]));
  int sigcount = totcount - bgcount;
  std::cout << " ***** \n";
  std::cout << " Total count: " << totcount << "\n";
  std::cout << " Background count: " << bgcount << "\n";
  std::cout << " Signal count: " << sigcount << "\n\n";
  // --- 

  // Canvas 2 : Fitting w/ signal from data & polynomial background
  TCanvas *c2 = util_pd::TC("c2",1,1);
  c2->cd(); gStyle->SetOptFit(1);
  vector<TH1F*> ho2;
  TF1 *f2 = fit::fit_1hs_1pbg_THI(W2_fit_range,
  				  h_W2_data_raw,h_W2_data,Opoly,
  				  ho2);
  ho2[0]->Draw(); customize_ht(ho2[0]); customize_W2(ho2[0]);
  ho2[1]->Draw("same"); customize_hs(ho2[1]); customize_W2(ho2[1]);
  ho2[2]->Draw("same"); customize_hbg(ho2[2]);
  TLegend *l2=new TLegend(0.10,0.77,0.39,0.9);
  l2->SetTextFont(42);
  l2->AddEntry(ho2[0],"Data","l");
  l2->AddEntry(f2,"Global Fit (MC + poly. bg.)","l");
  l2->AddEntry(ho2[1],"Signal (from data)","lep");
  l2->AddEntry(ho2[2],Form("Bg. (%d^{th} order poly.)",Opoly),"lep");
  l2->Draw();
  std::cout << Form("\nTotal # elastics (from fitted data signal): %d\n",
		    int(ho2[1]->Integral(ho2[1]->FindBin(W2_fit_range[0]),ho2[1]->FindBin(W2_fit_range[1]))));

  return 0;
}
