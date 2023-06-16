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

void customize_dx(TH1F* h)
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

int fit_dx (const char *configfilename, std::string filebase="pdout/test_fit_dx") 
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

  // defining output files
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d.root",filebase.c_str(),conf,sbsmag,model);

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/qelas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",conf,sbsmag,model,pass));
  //ROOT::RDataFrame simu_rdf("Tout",Form("simulation/siout/qelas_ana_%s_sbs%d_sbs%dp_model%d.root",gen.c_str(),conf,sbsmag,model));
  //ROOT::RDataFrame simu_rdf("Tout",Form("simulation/siout/65sf_qelas_ana_%s_sbs%d_sbs%dp_model%d.root",gen.c_str(),conf,sbsmag,model));
  ROOT::RDataFrame simu_rdf("Tout",Form("simulation/siout/tmp_hoff75cm_40sf_qelas_ana_%s_sbs%d_sbs%dp_model%d.root",gen.c_str(),conf,sbsmag,model));

  // Applying cuts
  std::string cuts_for_signal = jmgr->GetValueFromKey_str("cuts_for_signal");
  std::string cuts_for_bg = jmgr->GetValueFromKey_str("cuts_for_bg");
  auto data_rdf_filtered = data_rdf.Filter(cuts_for_signal);
  auto simu_rdf_filtered = simu_rdf.Filter(cuts_for_signal);
  auto bg_rdf_filtered = data_rdf.Filter(cuts_for_bg);

  // Creating important histograms
  vector<double> h_dx; jmgr->GetVectorFromKey<double>("h_dx", h_dx);
  auto h_dxHCAL_data = data_rdf_filtered.Histo1D({"h_dxHCAL_data", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx");
  auto h_dxHCAL_data_CT = data_rdf_filtered.Filter("abs((atimeHCAL-cltmeanHODO)-48.)<12.").Histo1D({"h_dxHCAL_data_CT", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx");
  auto h_dxHCAL_simu = simu_rdf_filtered.Histo1D({"h_dxHCAL_simu", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx", "weight");
  auto h_dxHCAL_simu_p = simu_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_simu_p", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx", "weight");
  auto h_dxHCAL_simu_n = simu_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_simu_n", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx", "weight");
  auto h_dxHCAL_bg = bg_rdf_filtered.Histo1D({"h_dxHCAL_bg", "", int(h_dx[0]), h_dx[1], h_dx[2]}, "dx");
  
  // converting to TH1F
  TH1F *h_dx_data = (TH1F*)h_dxHCAL_data->Clone();
  TH1F *h_dx_data_CT = (TH1F*)h_dxHCAL_data_CT->Clone();
  TH1F *h_dx_simu = (TH1F*)h_dxHCAL_simu->Clone();
  TH1F *h_dx_simu_p = (TH1F*)h_dxHCAL_simu_p->Clone();
  TH1F *h_dx_simu_n = (TH1F*)h_dxHCAL_simu_n->Clone();
  TH1F *h_dx_bg = (TH1F*)h_dxHCAL_bg->Clone();

  // TH1F *h_dx_data = new TH1F(*h_dxHCAL_data);
  // TH1F *h_dx_data_CT = new TH1F(*h_dxHCAL_data_CT);
  // TH1F *h_dx_simu = new TH1F(*h_dxHCAL_simu);
  // TH1F *h_dx_simu_p = new TH1F(*h_dxHCAL_simu_p);
  // TH1F *h_dx_simu_n = new TH1F(*h_dxHCAL_simu_n);
  // TH1F *h_dx_bg = new TH1F(*h_dxHCAL_bg);

  // Fits
  vector<double> fit_range; jmgr->GetVectorFromKey<double>("fit_range",fit_range);

  // Method 1: 2hs + nbg
  // TCanvas *c0 = util_pd::TC("c0",1,1);
  // c0->cd(); gStyle->SetOptFit(1);

  // --
  TCanvas *c0 = util_pd::TC("c0",1,1);
  c0->cd(); gStyle->SetOptFit(1);

  vector<TH1F*> ho;
  TF1 *f1 = fit::fit_2hs_nbg_THI(fit_range,
  				 h_dx_data_CT,h_dx_simu_p,h_dx_simu_n,
  				 ho);
  ho[0]->Draw(); customize_ht(ho[0]); customize_dx(ho[0]);
  // hsp0_c->Draw("same");
  // hsn0_c->Draw("same");
  ho[1]->Draw("same"); customize_hs(ho[1]); customize_dx(ho[1]);
  // legend
  TLegend *legend=new TLegend(0.10,0.77,0.38,0.9);
  legend->SetTextFont(42);
  //legend->SetTextSize(0.02);
  legend->AddEntry(ho[0],"Data","l");
  legend->AddEntry(f1,"Global Fit","l");
  legend->AddEntry(ho[1],"Signal","lep");
  // legend->AddEntry(psignal,"p Signal Fit","l");
  // legend->AddEntry(nsignal,"n Signal Fit","l");
  // legend->AddEntry(bg,"Background fit","l");
  legend->Draw();


  TCanvas *c3 = util_pd::TC("c3",1,1);
  c3->cd(); gStyle->SetOptFit(1);
  vector<TH1F*> ho1;
  TF1 *f2 = fit::fit_2hs_1hbg_THI(fit_range,
  				  h_dx_data_CT,h_dx_simu_p,h_dx_simu_n,h_dx_bg,
  				  ho1);
  ho1[0]->Draw(); customize_ht(ho1[0]); customize_dx(ho1[0]);
  ho1[1]->Draw("same"); customize_hs(ho1[1]); customize_dx(ho1[1]);
  // hsp0_c->Draw("same");
  // hsn0_c->Draw("same");
  ho1[2]->Draw("same"); customize_hbg(ho1[2]); customize_dx(ho1[2]);
  // l3
  TLegend *l3=new TLegend(0.10,0.77,0.38,0.9);
  l3->SetTextFont(42);
  //l3->SetTextSize(0.02);
  l3->AddEntry(ho1[0],"Data","l");
  l3->AddEntry(f2,"Global Fit (MC + Data bg.)","l");
  l3->AddEntry(ho1[1],"Signal (from MC)","lep");
  // l3->AddEntry(psignal,"p Signal Fit","l");
  // l3->AddEntry(nsignal,"n Signal Fit","l");
  l3->AddEntry(ho1[2],"Bg. (from Data)","lep");
  l3->Draw();

  int Opoly = 6;
  TCanvas *c4 = util_pd::TC("c4",1,1);
  c4->cd(); gStyle->SetOptFit(1);
  vector<TH1F*> ho2;
  TF1 *f3 = fit::fit_2hs_1pbg_THI(fit_range,
  				  h_dx_data_CT,h_dx_simu_p,h_dx_simu_n,Opoly,
  				  ho2);
  ho2[0]->Draw(); customize_ht(ho2[0]); customize_dx(ho2[0]);
  ho2[1]->Draw("same"); customize_hs(ho2[1]); customize_dx(ho2[1]);
  // hsp0_c->Draw("same");
  // hsn0_c->Draw("same");
  ho2[2]->Draw("same"); customize_hbg(ho2[2]);
  // l4
  TLegend *l4=new TLegend(0.10,0.77,0.39,0.9);
  l4->SetTextFont(42);
  //l4->SetTextSize(0.02);
  l4->AddEntry(ho2[0],"Data","l");
  l4->AddEntry(f3,"Global Fit (MC + poly. bg.)","l");
  l4->AddEntry(ho2[1],"Signal (from MC)","lep");
  // l4->AddEntry(psignal,"p Signal Fit","l");
  // l4->AddEntry(nsignal,"n Signal Fit","l");
  l4->AddEntry(ho2[2],Form("Bg. (%d^{th} order poly.)",Opoly),"lep");
  l4->Draw();

  return 0;
}
