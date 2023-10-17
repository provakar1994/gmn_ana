/*
  This macro will fit dx distribution in various different ways.
  Plan:
  1. MC signal, no bg
  2. MC signal + bg from data
  3. MC signal + poly bg
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

void customize_dx(TH1F* h)
{
  h->GetXaxis()->SetTitle("x_{HCAL}^{obs} - x_{HCAL}^{exp} (m)");
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

double total_fit (double * x, double * par) {
  FitFn *ffn = new FitFn(6);
  return ffn->ffn_gaus(x,&par[0]) + ffn->ffn_poly(x,&par[3]);
}

int fit_dx (const char *configfilename, 
	    bool is_elastic = 1, // 1=>Yes, 0=>QE
	    std::string filebase="pdout/test_fit_dx") 
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
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/%s%s_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfprefix.c_str(),key,conf,sbsmag,model,pass));
  ROOT::RDataFrame simu_rdf("Tout",Form("simulation/siout/%s%s_ana_%s_sbs%d_sbs%dp_model%d.root",sfprefix.c_str(),key,gen.c_str(),conf,sbsmag,model));

  // Applying cuts
  std::string cuts_for_signal_data = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_data");
  std::string cuts_for_signal_simu = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_simu");
  std::string cuts_for_bg = jmgr->GetValueFromSubKey_str(key,"cuts_for_bg");
  std::string coinT_cut = jmgr->GetValueFromSubKey_str(key,"coinT_cut");
  int Opoly = jmgr->GetValueFromSubKey<int>(key,"Order_of_poly_bg");
  auto data_rdf_filtered = data_rdf.Filter(cuts_for_signal_data);
  auto simu_rdf_filtered_1 = simu_rdf.Filter(cuts_for_signal_simu);
  double offset = jmgr->GetValueFromSubKey<double>(key,"dx_peak_offset_for_MC");
  std::string dx_shifted = "dx+" + std::to_string(offset);
  auto simu_rdf_filtered = simu_rdf_filtered_1.Define("dx_shifted",dx_shifted.c_str());
  auto bg_rdf_filtered = data_rdf.Filter(cuts_for_bg);

  // Creating important histograms
  vector<double> h_dx; jmgr->GetVectorFromSubKey<double>(key,"h_dx",h_dx);
  TH1F *h_dxHCAL_data = (TH1F*)data_rdf_filtered.Histo1D({"h_dxHCAL_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
  TH1F *h_dxHCAL_data_CT = (TH1F*)data_rdf_filtered.Filter(coinT_cut.c_str()).Histo1D({"h_dxHCAL_data_CT","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
  TH1F *h_dxHCAL_simu = (TH1F*)simu_rdf_filtered.Histo1D({"h_dxHCAL_simu","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted","weight")->Clone();
  TH1F *h_dxHCAL_simu_p = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_simu_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted","weight")->Clone();
  TH1F *h_dxHCAL_simu_n;
  if (!is_elastic) h_dxHCAL_simu_n = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_simu_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted","weight")->Clone();
  TH1F *h_dxHCAL_bg = (TH1F*)bg_rdf_filtered.Histo1D({"h_dxHCAL_bg","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();

  // Fits
  vector<double> dx_fit_range; jmgr->GetVectorFromSubKey<double>(key,"dx_fit_range",dx_fit_range);

  /*######################################
    ## Fitting elastic dx distributions ##
    ###################################### */
  if (is_elastic) {
    // Canvas 0 : Fitting data/MC w/o any background
    TCanvas *c0 = util_pd::TC("c0",1,1);
    c0->cd(); gStyle->SetOptFit(1);
    vector<TH1F*> ho;
    TF1 *f0 = fit::fit_1hs_nbg_THI(dx_fit_range,
				   h_dxHCAL_data_CT,h_dxHCAL_simu_p,
				   ho);
    ho[0]->Draw(); customize_ht(ho[0]); customize_dx(ho[0]);
    ho[1]->Draw("same"); customize_hs(ho[1]); customize_dx(ho[1]);
    TLegend *l0=new TLegend(0.10,0.77,0.38,0.9);
    l0->SetTextFont(42);
    l0->AddEntry(ho[0],"Data","l");
    l0->AddEntry(f0,"Global Fit","l");
    l0->AddEntry(ho[1],"Signal","lep");
    l0->Draw();
    // --- 

    // // Canvas 1 : Fitting data/MC w/ background from data
    // TCanvas *c1 = util_pd::TC("c1",1,1);
    // c1->cd(); gStyle->SetOptFit(1);
    // vector<TH1F*> ho1;
    // TF1 *f1 = fit::fit_1hs_1hbg_THI(dx_fit_range,
    // 				    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_bg,
    // 				    ho1);
    // ho1[0]->Draw(); customize_ht(ho1[0]); customize_dx(ho1[0]);
    // ho1[1]->Draw("same"); customize_hs(ho1[1]); customize_dx(ho1[1]);
    // TLegend *l1=new TLegend(0.10,0.77,0.38,0.9);
    // l1->SetTextFont(42);
    // l1->AddEntry(ho1[0],"Data","l");
    // l1->AddEntry(f1,"Global Fit (MC + Data bg.)","l");
    // l1->AddEntry(ho1[1],"Signal (from MC)","lep");
    // l1->AddEntry(ho1[2],"Bg. (from Data)","lep");
    // l1->Draw();
    // // --- 

    // Canvas 2 : Fitting data/MC w/ polynomial background
    TCanvas *c2 = util_pd::TC("c2",1,1);
    c2->cd(); gStyle->SetOptFit(1);
    vector<TH1F*> ho2;
    TF1 *f2 = fit::fit_1hs_1pbg_THI(dx_fit_range,
				    h_dxHCAL_data,h_dxHCAL_simu_p,Opoly,
				    ho2);
    ho2[0]->Draw(); customize_ht(ho2[0]); customize_dx(ho2[0]);
    ho2[1]->Draw("same"); customize_hs(ho2[1]); customize_dx(ho2[1]);
    //cout << " *** " << ho2[1]->Integral() << "\n";
    ho2[2]->Draw("same"); customize_hbg(ho2[2]);
    //cout << " *** " << ho2[2]->Integral() << "\n";
    // drawing the polynomial background as well
    FitFn *ffn = new FitFn(Opoly);
    TF1* bg2 = new TF1("bg2",ffn,&FitFn::ffn_poly,dx_fit_range[0],dx_fit_range[1],Opoly+1);
    bg2->SetNpx(500);
    bg2->SetParameters(&fit::GetFitParams(f2)[1]);
    bg2->SetLineColor(kGreen+2);
    bg2->Draw("same");
    TLegend *l2=new TLegend(0.10,0.77,0.39,0.9);
    l2->SetTextFont(42);
    l2->AddEntry(ho2[0],"Data","l");
    l2->AddEntry(f2,"Global Fit (MC + poly. bg.)","l");
    l2->AddEntry(ho2[1],"Signal (from MC)","lep");
    l2->AddEntry(ho2[2],Form("Bg. (%d^{th} order poly.)",Opoly),"lep");
    l2->Draw();    
    // ---

    // Canvas 3 : Sideband fit (distribution from data)
    TCanvas *c3 = util_pd::TC("c3",1,1);
    c3->cd(); gStyle->SetOptFit(1);
    vector<double> reject_points{-1.3,-0.15};
    TF1* bg3 = fit::fit_1pbg_SB(dx_fit_range,
				reject_points,
				Opoly,
				fit::GetFitParams(f2),
				h_dxHCAL_data);
    bg3->Draw("same");

    double bgcount = (int)bg3->Integral(dx_fit_range[0],dx_fit_range[1])/h_dxHCAL_data->GetBinWidth(1);
    double totcount = h_dxHCAL_data->Integral(h_dxHCAL_data->FindBin(dx_fit_range[0]),h_dxHCAL_data->FindBin(dx_fit_range[1]));
    int sigcount = totcount - bgcount;
    std::cout << " ***** Reporting # elastics from side band fit ***** \n";
    std::cout << " Total count: " << totcount << "\n";
    std::cout << " Background count: " << bgcount << "\n";
    std::cout << " Signal count: " << sigcount << "\n\n";
    // ---

    // Canvas 4 : Fitting using a Gaussian and fixed polynomial (got the params from f2)
    // TF1* fpfit = new TF1("fpfit",total_fit,-2,1,10);
    // fpfit->SetNpx(500);
    // fpfit->SetLineColor(kRed);
    // vector<double> parGuess{31587,-0.65,0.12,4219,-1892,-672,1204,-181,-517,-98};
    // fpfit->SetParameters(&parGuess[0]);
    // // fpfit->SetParameter(0,31587);
    // // fpfit->SetParameter(1,-0.65);
    // // fpfit->SetParameter(2,0.12);
    // // fpfit->FixParameter(3,4219);
    // // fpfit->FixParameter(4,-1892);
    // // fpfit->FixParameter(5,-672);
    // // fpfit->FixParameter(6,1204);
    // // fpfit->FixParameter(7,-181);
    // // fpfit->FixParameter(8,-517);
    // // fpfit->FixParameter(9,-98);
    // h_dxHCAL_data->Fit("fpfit","RV+","ep");

  } // elastic

  /*#################################
    ## Fitting QE dx distributions ##
    ################################# */
  else {
    // Canvas 0 : Fitting data/MC w/o any background
    TCanvas *c0 = util_pd::TC("c0",1,1);
    c0->cd(); gStyle->SetOptFit(1);
    vector<TH1F*> ho;
    TF1 *f0 = fit::fit_2hs_nbg_THI(dx_fit_range,
				   h_dxHCAL_data_CT,h_dxHCAL_simu_p,h_dxHCAL_simu_n,
				   ho);
    ho[0]->Draw(); customize_ht(ho[0]); customize_dx(ho[0]);
    ho[1]->Draw("same"); customize_hs(ho[1]); customize_dx(ho[1]);
    TLegend *l0=new TLegend(0.10,0.77,0.38,0.9);
    l0->SetTextFont(42);
    l0->AddEntry(ho[0],"Data","l");
    l0->AddEntry(f0,"Global Fit","l");
    l0->AddEntry(ho[1],"Signal","lep");
    l0->Draw();
    // --- 

    // Canvas 1 : Fitting data/MC w/ background from data
    TCanvas *c1 = util_pd::TC("c1",1,1);
    c1->cd(); gStyle->SetOptFit(1);
    vector<TH1F*> ho1;
    TF1 *f1 = fit::fit_2hs_1hbg_THI(dx_fit_range,
				    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg,
				    ho1);
    ho1[0]->Draw(); customize_ht(ho1[0]); customize_dx(ho1[0]);
    ho1[1]->Draw("same"); customize_hs(ho1[1]); customize_dx(ho1[1]);
    //cout << " *** " << ho1[1]->Integral() << "\n";
    ho1[2]->Draw("same"); customize_hbg(ho1[2]); customize_dx(ho1[2]);
    //cout << " *** " << ho1[2]->Integral() << "\n";
    TLegend *l1=new TLegend(0.10,0.77,0.38,0.9);
    l1->SetTextFont(42);
    l1->AddEntry(ho1[0],"Data","l");
    l1->AddEntry(f1,"Global Fit (MC + Data bg.)","l");
    l1->AddEntry(ho1[1],"Signal (from MC)","lep");
    l1->AddEntry(ho1[2],"Bg. (from Data)","lep");
    l1->Draw();
    // --- 

    // Canvas 2 : Fitting data/MC w/ polynomial background
    TCanvas *c2 = util_pd::TC("c2",1,1);
    c2->cd(); gStyle->SetOptFit(1);
    vector<TH1F*> ho2;
    TF1 *f2 = fit::fit_2hs_1pbg_THI(dx_fit_range,
				    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,Opoly,
				    ho2);
    ho2[0]->Draw(); customize_ht(ho2[0]); customize_dx(ho2[0]);
    ho2[1]->Draw("same"); customize_hs(ho2[1]); customize_dx(ho2[1]);
    ho2[2]->Draw("same"); customize_hbg(ho2[2]);
    // drawing the polynomial background as well
    FitFn *ffn = new FitFn(Opoly);
    TF1* bg2 = new TF1("bg2",ffn,&FitFn::ffn_poly,dx_fit_range[0],dx_fit_range[1],Opoly+1);
    bg2->SetNpx(500);
    bg2->SetParameters(&fit::GetFitParams(f2)[2]);
    bg2->SetLineColor(kGreen+2);
    bg2->Draw("same");
    TLegend *l2=new TLegend(0.10,0.77,0.39,0.9);
    l2->SetTextFont(42);
    l2->AddEntry(ho2[0],"Data","l");
    l2->AddEntry(f2,"Global Fit (MC + poly. bg.)","l");
    l2->AddEntry(ho2[1],"Signal (from MC)","lep");
    l2->AddEntry(ho2[2],Form("Bg. (%d^{th} order poly.)",Opoly),"lep");
    l2->Draw();
  } // QE

  return 0;
}
