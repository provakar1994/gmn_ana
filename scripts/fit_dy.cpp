/*
  This macro will fit dx in various different ways.
  Plan:
  1. MC signal, no bg
  2. MC signal + bg from data
  3. MC signal + poly bg
  -------
  P. Datta Created 05-02-2023 (Based on a script written by AJRP)
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

// TLegend* TL_fithistos(std::vector<TH1F*> const & ho){
//   TLegend *l = new TLegend;
// }

void customize_dy(TH1F* h)
{
  h->GetXaxis()->SetTitle("xHCAL_{obs} - xHCAL_{exp} (m)");
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

int fit_dy (const char *configfilename, std::string filebase="pdout/test_fit_dy") 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // reading in proper data and simu output files
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  int model = jmgr->GetValueFromKey<int>("model");
  int pass = jmgr->GetValueFromKey<int>("pass");
  std::string gen = jmgr->GetValueFromKey_str("generator");
  std::string dfprefix = jmgr->GetValueFromKey_str("data_file_prefix");
  std::string sfprefix = jmgr->GetValueFromKey_str("simu_file_prefix");
  char const * dfp = dfprefix.empty() ? "" : (dfprefix + "_").c_str();
  char const * sfp = sfprefix.empty() ? "" : (sfprefix + "_").c_str();

  // defining output files
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d.root",filebase.c_str(),conf,sbsmag,model);

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/%sqelas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfp,conf,sbsmag,model,pass));
  ROOT::RDataFrame simu_rdf("Tout",Form("simulation/siout/%sqelas_ana_%s_sbs%d_sbs%dp_model%d.root",sfp,gen.c_str(),conf,sbsmag,model));

  // Applying cuts
  std::string cuts_for_signal_data = jmgr->GetValueFromKey_str("cuts_for_signal_data");
  std::string cuts_for_signal_simu = jmgr->GetValueFromKey_str("cuts_for_signal_simu");
  std::string cuts_for_bg = jmgr->GetValueFromKey_str("cuts_for_bg");
  std::string coinT_cut = jmgr->GetValueFromKey_str("coinT_cut");
  int Opoly = jmgr->GetValueFromKey<int>("Order_of_poly_bg");
  auto data_rdf_filtered = data_rdf.Filter(cuts_for_signal_data);
  auto simu_rdf_filtered = simu_rdf.Filter(cuts_for_signal_simu);
  auto bg_rdf_filtered = data_rdf.Filter(cuts_for_bg);

  // Creating important histograms
  vector<double> h_dy; jmgr->GetVectorFromKey<double>("h_dy",h_dy);
  TH1F *h_dyHCAL_data = (TH1F*)data_rdf_filtered.Histo1D({"h_dyHCAL_data","",int(h_dy[0]),h_dy[1],h_dy[2]},"dy")->Clone();
  TH1F *h_dyHCAL_data_CT = (TH1F*)data_rdf_filtered.Filter(coinT_cut.c_str()).Histo1D({"h_dyHCAL_data_CT","",int(h_dy[0]),h_dy[1],h_dy[2]},"dy")->Clone();
  TH1F *h_dyHCAL_simu = (TH1F*)simu_rdf_filtered.Histo1D({"h_dyHCAL_simu","",int(h_dy[0]),h_dy[1],h_dy[2]},"dy","weight")->Clone();
  TH1F *h_dyHCAL_simu_p = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dyHCAL_simu_p","",int(h_dy[0]),h_dy[1],h_dy[2]},"dy","weight")->Clone();
  TH1F *h_dyHCAL_simu_n = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dyHCAL_simu_n","",int(h_dy[0]),h_dy[1],h_dy[2]},"dy","weight")->Clone();
  TH1F *h_dyHCAL_bg = (TH1F*)bg_rdf_filtered.Histo1D({"h_dyHCAL_bg","",int(h_dy[0]),h_dy[1],h_dy[2]},"dy")->Clone();

  // Fits
  vector<double> dy_fit_range; jmgr->GetVectorFromKey<double>("dy_fit_range",dy_fit_range);

  // --
  TCanvas *c0 = util_pd::TC("c0",1,1);
  c0->cd(); gStyle->SetOptFit(1);

  vector<TH1F*> ho;
  TF1 *f1 = fit::fit_1hs_nbg_THI(dy_fit_range,
  				 h_dyHCAL_data_CT,h_dyHCAL_simu,
  				 ho);
  ho[0]->Draw(); customize_ht(ho[0]); customize_dy(ho[0]);
  // hsp0_c->Draw("same");
  // hsn0_c->Draw("same");
  ho[1]->Draw("same"); customize_hs(ho[1]); customize_dy(ho[1]);
  // l0
  TLegend *l0=new TLegend(0.10,0.77,0.38,0.9);
  l0->SetTextFont(42);
  //l0->SetTextSize(0.02);
  l0->AddEntry(ho[0],"Data","l");
  l0->AddEntry(f1,"Global Fit","l");
  l0->AddEntry(ho[1],"Signal","lep");
  // l0->AddEntry(psignal,"p Signal Fit","l");
  // l0->AddEntry(nsignal,"n Signal Fit","l");
  // l0->AddEntry(bg,"Background fit","l");
  l0->Draw();

  // TCanvas *c1 = util_pd::TC("c1",1,1);
  // c1->cd(); gStyle->SetOptFit(1);
  // vector<TH1F*> ho1;
  // TF1 *f2 = fit::fit_2hs_1hbg_THI(dy_fit_range,
  // 				  h_dyHCAL_data,h_dyHCAL_simu_p,h_dyHCAL_simu_n,h_dyHCAL_bg,
  // 				  ho1);
  // ho1[0]->Draw(); customize_ht(ho1[0]); customize_dy(ho1[0]);
  // ho1[1]->Draw("same"); customize_hs(ho1[1]); customize_dy(ho1[1]);
  // // hsp0_c->Draw("same");
  // // hsn0_c->Draw("same");
  // ho1[2]->Draw("same"); customize_hbg(ho1[2]); customize_dy(ho1[2]);
  // // l2
  // TLegend *l2=new TLegend(0.10,0.77,0.38,0.9);
  // l2->SetTextFont(42);
  // //l2->SetTextSize(0.02);
  // l2->AddEntry(ho1[0],"Data","l");
  // l2->AddEntry(f2,"Global Fit (MC + Data bg.)","l");
  // l2->AddEntry(ho1[1],"Signal (from MC)","lep");
  // // l2->AddEntry(psignal,"p Signal Fit","l");
  // // l2->AddEntry(nsignal,"n Signal Fit","l");
  // l2->AddEntry(ho1[2],"Bg. (from Data)","lep");
  // l2->Draw();

  TCanvas *c2 = util_pd::TC("c2",1,1);
  c2->cd(); gStyle->SetOptFit(1);
  vector<TH1F*> ho2;
  TF1 *f3 = fit::fit_1hs_1pbg_THI(dy_fit_range,
  				  h_dyHCAL_data_CT,h_dyHCAL_simu,Opoly,
  				  ho2);
  ho2[0]->Draw(); customize_ht(ho2[0]); customize_dy(ho2[0]);
  ho2[1]->Draw("same"); customize_hs(ho2[1]); customize_dy(ho2[1]);
  // hsp0_c->Draw("same");
  // hsn0_c->Draw("same");
  ho2[2]->Draw("same"); customize_hbg(ho2[2]);
  // l2
  TLegend *l2=new TLegend(0.10,0.77,0.39,0.9);
  l2->SetTextFont(42);
  //l2->SetTextSize(0.02);
  l2->AddEntry(ho2[0],"Data","l");
  l2->AddEntry(f3,"Global Fit (MC + poly. bg.)","l");
  l2->AddEntry(ho2[1],"Signal (from MC)","lep");
  // l2->AddEntry(psignal,"p Signal Fit","l");
  // l2->AddEntry(nsignal,"n Signal Fit","l");
  l2->AddEntry(ho2[2],Form("Bg. (%d^{th} order poly.)",Opoly),"lep");
  l2->Draw();

  return 0;
}
