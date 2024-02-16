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
#include "TPad.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include "../include/gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

void gStyleFitCanvas() 
{
  gStyle->SetOptStat("e"); gStyle->SetOptFit(1); 
  gStyle->SetErrorX(0);
}

void customize_residual(TH1F* h)
{
  h->GetXaxis()->SetLabelOffset(0.03);
  h->GetXaxis()->SetLabelSize(0.12);
  h->GetYaxis()->SetLabelOffset(0.005);
  h->GetYaxis()->SetLabelSize(0.11);
  h->SetMarkerStyle(22);
  h->SetMarkerColor(46);
  h->SetLineColor(46);
  h->SetStats(0);
}

void customize_data(TH1F* h)
{
  h->SetMarkerStyle(21);
  h->SetMarkerSize(0.8);
  h->SetMarkerColor(kBlack);
  h->SetLineColor(kBlack);
}

void customize_gfit(TH1F* h)
{
  h->SetLineColor(kRed);
  h->SetLineStyle(7);
  h->SetLineWidth(3);
  //***
  h->SetFillColor(kRed);
  h->SetFillColorAlpha(kRed,0.2);
}

void customize_psig(TH1F* h)
{
  h->SetLineColor(kBlue);
  h->SetLineStyle(4);
  h->SetLineWidth(3);
  //***
  h->SetFillColor(kBlue);
  h->SetFillColorAlpha(kBlue,0.3);
}

void customize_nsig(TH1F* h)
{
  h->SetLineColor(kGreen+2);
  h->SetLineStyle(8);
  h->SetLineWidth(3);
  //***
  h->SetFillColor(kGreen+2);
  h->SetFillColorAlpha(kGreen+2,0.3);
}

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
  h->SetMarkerStyle(22);
  h->SetLineColor(kBlack);
}

void customize_hbg(TH1F* h) 
{
  h->SetMarkerColor(6);
  h->SetMarkerStyle(29);
  h->SetLineColor(6);
  h->SetLineStyle(9);
  h->SetLineWidth(3);
  //***
  h->SetFillColor(6);
  h->SetFillColorAlpha(6,0.3);
}

double total_fit (double * x, double * par) {
  FitFn *ffn = new FitFn(6);
  return ffn->ffn_gaus(x,&par[0]) + ffn->ffn_poly(x,&par[3]);
}

int fit_dx (const char *configfilename, 
	    bool is_elastic = 1, // 1=>Yes, 0=>QE
	    std::string filebase="pdout/fit_dx") 
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
  std::string infprefix = jmgr->GetValueFromSubKey_str(key,"inel_file_prefix");
  dfprefix = dfprefix.empty() ? "" : dfprefix + "_";
  sfprefix = sfprefix.empty() ? "" : sfprefix + "_";
  infprefix = infprefix.empty() ? "" : infprefix + "_";

  // defining output files
  TString outFile = Form("%s_%s_pass%d_%s_sbs%d_sbs%dp_model%d.root",filebase.c_str(),key,pass,gen.c_str(),conf,sbsmag,model);
  TString outPlot = outFile; outPlot.ReplaceAll(".root",".pdf");
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/%s%s_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfprefix.c_str(),key,conf,sbsmag,model,pass));
  ROOT::RDataFrame simu_rdf("Tout",Form("simulation/siout/%s%s_ana_%s_sbs%d_sbs%dp_model%d.root",sfprefix.c_str(),key,gen.c_str(),conf,sbsmag,model));
  ROOT::RDataFrame inel_rdf("Tout",Form("simulation/siout/%sinel_%s_ana_%s_sbs%d_sbs%dp_model%d.root",infprefix.c_str(),key,gen.c_str(),conf,sbsmag,model));

  // ***********
  // Inelastic generator. make this part user configuration later. Add a flag, etc.
  // ***********
  // ** sbs4
  //ROOT::RDataFrame inel_rdf("Tout",Form("simulation/siout/inel_elas_ana_g4sbs_sbs4_sbs0p_model2.root"));
  //ROOT::RDataFrame inel_rdf("Tout",Form("simulation/siout/inel_qelas_ana_g4sbs_sbs4_sbs50p_model2.root"));
  //auto inel_rdf_filtered = inel_rdf.Filter("W2>0.89&&trP>1.16&&eHCAL>0");
  // ** --
  ////ROOT::RDataFrame inel_rdf("Tout",Form("simulation/siout/0p815sf_inel_elas_ana_g4sbs_sbs7_sbs85p_model2.root"));
  //ROOT::RDataFrame inel_rdf("Tout",Form("simulation/siout/0p815sf_inel_qelas_ana_g4sbs_sbs7_sbs85p_model2.root"));
  ////auto inel_rdf_filtered = inel_rdf.Filter("W2>0.89&&trP>1.2&&eHCAL>0&&fiduCut");
  ////vector<double> h_dx; jmgr->GetVectorFromSubKey<double>(key,"h_dx",h_dx); //temporary replacement
  ////TH1F *h_dxHCAL_bg_inel = (TH1F*)inel_rdf_filtered.Histo1D({"h_dxHCAL_bg_inel","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight")->Clone();
  // ***********

  // Applying cuts
  std::string cuts_for_signal_data = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_data");
  std::string cuts_for_signal_simu = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_simu");
  std::string cuts_for_bg_data = jmgr->GetValueFromSubKey_str(key,"cuts_for_bg_data");
  std::string cuts_for_bg_simu = jmgr->GetValueFromSubKey_str(key,"cuts_for_bg_simu");
  std::string coinT_cut = jmgr->GetValueFromSubKey_str(key,"coinT_cut");
  int Opoly = jmgr->GetValueFromSubKey<int>(key,"Order_of_poly_bg");
  auto data_rdf_filtered = data_rdf.Filter(cuts_for_signal_data);
  auto simu_rdf_filtered_1 = simu_rdf.Filter(cuts_for_signal_simu);
  auto inel_rdf_filtered = inel_rdf.Filter("W2>0.89&&trP>1.2&&eHCAL>0&&fiduCut");
  // Reading in offsets for dx peak position in MC for both p and n
  double dx_offset_p = jmgr->GetValueFromSubKey<double>(key,"dx_offset_MC_for_p");
  double dx_offset_n = !is_elastic ? jmgr->GetValueFromSubKey<double>(key,"dx_offset_MC_for_n") : 0;
  std::string dx_shifted_p = "dx+" + std::to_string(dx_offset_p);
  std::string dx_shifted_n = "dx+" + std::to_string(dx_offset_n);
  auto simu_rdf_filtered = simu_rdf_filtered_1.Define("dx_shifted_p",dx_shifted_p.c_str())
    .Define("dx_shifted_n",dx_shifted_n.c_str());
  auto bg_data_rdf_filtered = data_rdf.Filter(cuts_for_bg_data);

  // Creating important histograms
  vector<double> h_dx; jmgr->GetVectorFromSubKey<double>(key,"h_dx",h_dx);
  TH1F *h_dxHCAL_data = (TH1F*)data_rdf_filtered.Histo1D({"h_dxHCAL_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
  TH1F *h_dxHCAL_data_CT = (TH1F*)data_rdf_filtered.Filter(coinT_cut.c_str()).Histo1D({"h_dxHCAL_data_CT","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
  //TH1F *h_dxHCAL_simu = (TH1F*)simu_rdf_filtered.Histo1D({"h_dxHCAL_simu","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted","weight")->Clone();
  TH1F *h_dxHCAL_simu_p = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_simu_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight")->Clone();
  TH1F *h_dxHCAL_simu_n;
  if (!is_elastic) h_dxHCAL_simu_n = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_simu_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight")->Clone();
  TH1F *h_dxHCAL_bg_data = (TH1F*)bg_data_rdf_filtered.Histo1D({"h_dxHCAL_bg_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
  TH1F *h_dxHCAL_bg_inel = (TH1F*)inel_rdf_filtered.Histo1D({"h_dxHCAL_bg_inel","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight")->Clone();

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
    l0->AddEntry(f0,"Fit","l");
    l0->AddEntry(ho[1],"Signal","lep");
    l0->Draw();
    // --- 

    // // Canvas 1 : Fitting data/MC w/ background from data
    // TCanvas *c1 = util_pd::TC("c1",1,1);
    // c1->cd(); gStyle->SetOptFit(1);
    // vector<TH1F*> ho1;
    // TF1 *f1 = fit::fit_1hs_1hbg_THI(dx_fit_range,
    // 				    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_bg_data,
    // 				    ho1);
    // ho1[0]->Draw(); customize_ht(ho1[0]); customize_dx(ho1[0]);
    // ho1[1]->Draw("same"); customize_hs(ho1[1]); customize_dx(ho1[1]);
    // TLegend *l1=new TLegend(0.10,0.77,0.38,0.9);
    // l1->SetTextFont(42);
    // l1->AddEntry(ho1[0],"Data","l");
    // l1->AddEntry(f1,"Fit (MC + Data bg.)","l");
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
    l2->AddEntry(f2,"Fit (MC + poly. bg.)","l");
    l2->AddEntry(ho2[1],"Signal (from MC)","lep");
    l2->AddEntry(ho2[2],Form("Bg. (%d^{th} order poly.)",Opoly),"lep");
    l2->Draw();    
    // ---

    // Canvas 3 : Sideband fit (distribution from data)
    TCanvas *c3 = util_pd::TC("c3",1,1);
    c3->cd(); gStyle->SetOptFit(1);
    vector<double> reject_points{-1.3,-0.15};
    vector<TH1F*> ho3;
    TF1* bg3 = fit::fit_1pbg_SB(dx_fit_range,
				reject_points,
				Opoly,
				fit::GetFitParams(f2),
				h_dxHCAL_data,
				ho3);
    //bg3->Draw("same");
    ho3[0]->Draw(); customize_ht(ho3[0]); customize_dx(ho3[0]);
    ho3[1]->Draw("same ep"); customize_hs(ho3[1]); customize_dx(ho3[1]);
    ho3[2]->Draw("same"); customize_hbg(ho3[2]); customize_dx(ho3[2]);    

    //double bgcount = (int)bg3->Integral(dx_fit_range[0],dx_fit_range[1])/h_dxHCAL_data->GetBinWidth(1);
    double bgcount = ho3[2]->Integral(ho3[2]->FindBin(dx_fit_range[0]),ho3[2]->FindBin(dx_fit_range[1]));;
    double totcount = h_dxHCAL_data->Integral(h_dxHCAL_data->FindBin(dx_fit_range[0]),h_dxHCAL_data->FindBin(dx_fit_range[1]));
    int sigcount = totcount - bgcount;
    int sigcount_2ndmethod = ho3[1]->Integral(ho3[1]->FindBin(dx_fit_range[0]),ho3[1]->FindBin(dx_fit_range[1]));;
    std::cout << " ***** Reporting # elastics from side band fit ***** \n";
    std::cout << " Total count: " << totcount << "\n";
    std::cout << " Background count: " << bgcount << "\n";
    std::cout << " Signal count: " << sigcount << "\n";
    std::cout << " Signal count (2nd method): " << sigcount_2ndmethod << "\n\n";
    // ---

    // Canvas 4 : Fitting w/ signal and background from MC
    TCanvas *c4 = util_pd::TC("c4",1,1);
    c4->cd(); gStyle->SetOptFit(1);
    vector<TH1F*> ho4;
    TF1 *f4 = fit::fit_1hs_1hbg_THI(dx_fit_range,
				    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_bg_inel,
				    ho4);
    ho4[0]->Draw(); customize_ht(ho4[0]); customize_dx(ho4[0]);
    ho4[1]->Draw("same"); customize_hs(ho4[1]); customize_dx(ho4[1]);
    ho4[2]->Draw("same"); customize_hbg(ho4[2]); customize_dx(ho4[2]);
    TLegend *l4=new TLegend(0.10,0.77,0.38,0.9);
    l4->SetTextFont(42);
    l4->AddEntry(ho4[0],"Data","l");
    l4->AddEntry(f4,"Fit","l");
    l4->AddEntry(ho4[1],"Signal (from MC)","lep");
    l4->AddEntry(ho4[2],"Bg. (from MC)","lep");
    l4->Draw();
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
    std::vector<double> R_vals, Rerr_vals;

    // Canvas 0 : Fitting data/MC w/o any background
    TCanvas *c0 = util_pd::TC("c0",1,1);
    c0->cd(); gStyle->SetOptFit(1);
    vector<TH1F*> ho;
    TF1 *f0 = fit::fit_2hs_nbg_THI(dx_fit_range,
				   h_dxHCAL_data_CT,h_dxHCAL_simu_p,h_dxHCAL_simu_n,
				   ho);
    R_vals.push_back(f0->GetParameter(1)); Rerr_vals.push_back(f0->GetParError(1));
    ho[0]->Draw(); customize_ht(ho[0]); customize_dx(ho[0]);
    ho[1]->Draw("same"); customize_hs(ho[1]); customize_dx(ho[1]);
    TLegend *l0=new TLegend(0.10,0.77,0.38,0.9);
    l0->SetTextFont(42);
    l0->AddEntry(ho[0],"Data","l");
    l0->AddEntry(f0,"Fit","l");
    l0->AddEntry(ho[1],"Signal","lep");
    l0->Draw();
    // --- 

    // Canvas 1 : Fitting data/MC w/ background from data
    TCanvas *c1 = util_pd::TC("c1",1,1); gStyleFitCanvas();
    c1->cd();
    // performing the fit
    vector<TH1F*> ho1;
    TF1 *f1 = fit::fit_2hs_1hbg_THI(dx_fit_range,
    				    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_data,
    				    ho1);
    R_vals.push_back(f1->GetParameter(1)); Rerr_vals.push_back(f1->GetParError(1));
    // converting fit fn to a hostogram
    TH1F *hf1 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f1,hf1);
    ho1[0]->Draw(); c1->Update(); 
    // grabbing statbox of the fitted histo
    TPaveStats *st1 = (TPaveStats*)ho1[0]->FindObject("stats");
    // getting pads for pull plot
    std::vector<TPad*> p1 = util_pd::GetPadsForPullPlot(c1);
    //
    // preparing the pad for data/MC fit
    //
    p1[0]->cd();
    // drawing all the histograms
    //ho1[0]->Draw(); customize_ht(ho1[0]); customize_dx(ho1[0]);
    h_dxHCAL_data->Draw("E"); customize_data(h_dxHCAL_data);
    hf1->Draw("same HIST"); customize_gfit(hf1);
    ho1[4]->Draw("same HIST"); customize_psig(ho1[4]);
    ho1[5]->Draw("same HIST"); customize_nsig(ho1[5]);
    ho1[2]->Draw("same HIST"); customize_hbg(ho1[2]);
    // redrawing the stat box
    st1->SetX1NDC(0.6); st1->SetX2NDC(0.9); st1->SetY2NDC(0.9);
    st1->Draw("same");    
    // drawing a legend
    TLegend *l1=new TLegend(0.10,0.64,0.35,0.9);
    l1->SetTextFont(42);
    l1->AddEntry(h_dxHCAL_data,"Data","p");
    l1->AddEntry(hf1,"Fit (QE MC + bg.)","lf");
    l1->AddEntry(ho1[4],"p signal (from MC)","lf");
    l1->AddEntry(ho1[5],"n signal (from MC)","lf");
    l1->AddEntry(ho1[2],"Bg. (from Data)","lf");
    l1->AddEntry(ho1[3],"Residual","p");
    l1->Draw();
    //
    // preparing the pad for residual
    //  
    p1[1]->cd();
    ho1[3]->Draw(); customize_residual(ho1[3]);
    // drawing a horizontal line at y = 0
    util_pd::DrawZeroLine(p1[1],dx_fit_range[0],dx_fit_range[1]);
    // ** ----- ***

    // Canvas 2 : Fitting data/MC w/ background from data
    TCanvas *c2 = util_pd::TC("c2",1,1); gStyleFitCanvas();
    c2->cd();
    // performing the fit
    vector<TH1F*> ho2;
    TF1 *f2 = fit::fit_2hs_1hbg_THI(dx_fit_range,
    				    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_inel,
    				    ho2);
    R_vals.push_back(f2->GetParameter(1)); Rerr_vals.push_back(f2->GetParError(1));
    // converting fit fn to a hostogram
    TH1F *hf2 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f2,hf2);
    ho2[0]->Draw(); c2->Update(); 
    // grabbing statbox of the fitted histo
    TPaveStats *st2 = (TPaveStats*)ho2[0]->FindObject("stats");
    // getting pads for pull plot
    std::vector<TPad*> p2 = util_pd::GetPadsForPullPlot(c2);
    //
    // preparing the pad for data/MC fit
    //
    p2[0]->cd();
    // drawing all the histograms
    //ho2[0]->Draw(); customize_ht(ho2[0]); customize_dx(ho2[0]);
    h_dxHCAL_data->Draw("E"); customize_data(h_dxHCAL_data);
    hf2->Draw("same HIST"); customize_gfit(hf2);
    ho2[4]->Draw("same HIST"); customize_psig(ho2[4]);
    ho2[5]->Draw("same HIST"); customize_nsig(ho2[5]);
    ho2[2]->Draw("same HIST"); customize_hbg(ho2[2]);
    // redrawing the stat box
    st2->SetX1NDC(0.6); st2->SetX2NDC(0.9); st2->SetY2NDC(0.9);
    st2->Draw("same");    
    // drawing a legend
    TLegend *l2=new TLegend(0.10,0.64,0.35,0.9);
    l2->SetTextFont(42);
    l2->AddEntry(h_dxHCAL_data,"Data","p");
    l2->AddEntry(hf2,"Fit (QE MC + bg.)","lf");
    l2->AddEntry(ho2[4],"p signal (from MC)","lf");
    l2->AddEntry(ho2[5],"n signal (from MC)","lf");
    l2->AddEntry(ho2[2],"Bg. (from MC)","lf");
    l2->AddEntry(ho2[3],"Residual","p");
    l2->Draw();
    //
    // preparing the pad for residual
    //  
    p2[1]->cd();
    ho2[3]->Draw(); customize_residual(ho2[3]);
    // drawing a horizontal line at y = 0
    util_pd::DrawZeroLine(p2[1],dx_fit_range[0],dx_fit_range[1]);
    // ** ----- ***
    
    // Canvas 3 : Fitting data/MC w/ polynomial background
    TCanvas *c3 = util_pd::TC("c3",1,1); gStyleFitCanvas();
    c3->cd(); 
    // performing the fit
    vector<TH1F*> ho3;
    TF1 *f3 = fit::fit_2hs_1pbg_THI(dx_fit_range,
    				    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,Opoly,
    				    ho3);
    R_vals.push_back(f3->GetParameter(1)); Rerr_vals.push_back(f3->GetParError(1));
    // converting fit fn to a hostogram
    TH1F *hf3 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f3,hf3);
    ho3[0]->Draw(); c3->Update(); 
    // grabbing statbox of the fitted histo
    TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
    // getting pads for pull plot
    std::vector<TPad*> p3 = util_pd::GetPadsForPullPlot(c3);
    //
    // preparing the pad for data/MC fit
    //
    p3[0]->cd();
    // drawing all the histograms
    //ho3[0]->Draw(); customize_ht(ho3[0]); customize_dx(ho3[0]);
    h_dxHCAL_data->Draw("E"); customize_data(h_dxHCAL_data);
    hf3->Draw("same HIST"); customize_gfit(hf3);
    ho3[4]->Draw("same HIST"); customize_psig(ho3[4]);
    ho3[5]->Draw("same HIST"); customize_nsig(ho3[5]);
    ho3[2]->Draw("same HIST"); customize_hbg(ho3[2]);
    // redrawing the stat box
    st3->SetX1NDC(0.62); st3->SetX2NDC(0.9); st3->SetY2NDC(0.9);
    st3->Draw("same");
    // drawing a legend
    TLegend *l3=new TLegend(0.10,0.6,0.36,0.9);
    l3->SetTextFont(42);
    l3->AddEntry(h_dxHCAL_data,"Data","p");
    //l3->AddEntry(f3,"Fit (MC + poly. bg.)","l");
    l3->AddEntry(hf3,"Fit (QE MC + bg.)","lf");
    l3->AddEntry(ho3[4],"p signal (from MC)","lf");
    l3->AddEntry(ho3[5],"n signal (from MC)","lf");
    l3->AddEntry(ho3[2],Form("Bg. (poly. of order %d)",Opoly),"lf");
    l3->AddEntry(ho3[3],"Residual","p");
    //l3->SetFillStyle(0); // makes legend box transparent
    l3->Draw();
    //
    // preparing the pad for residual
    //  
    p3[1]->cd();
    ho3[3]->Draw(); customize_residual(ho3[3]);
    // drawing a horizontal line at y = 0
    util_pd::DrawZeroLine(p3[1],dx_fit_range[0],dx_fit_range[1]);

    // Writing out fit parameters
    std::cout << "\n--- Reporting fit params ---\n";
    std::cout << "R0,R0err,R1,R1err,R2,R2err,R3,R3err\n";
    for(size_t i=0; i < R_vals.size(); i++){
      std::cout << R_vals[i] << "," << Rerr_vals[i] << ",";
    }
    std::cout << "\n------\n";

    // writing out the canvases
    c0->Update(); c0->Write(); c0->SaveAs(Form("%s[",outPlot.Data())); c0->SaveAs(Form("%s",outPlot.Data())); 
    // turning off the stats of data histo [IMPORTANT!]
    h_dxHCAL_data->SetStats(0);
    c1->Update(); c1->Write(); c1->SaveAs(Form("%s",outPlot.Data())); 
    c2->Update(); c2->Write(); c2->SaveAs(Form("%s",outPlot.Data())); 
    c3->Update(); c3->Write(); c3->SaveAs(Form("%s",outPlot.Data())); c3->SaveAs(Form("%s]",outPlot.Data())); 
  } // QE

  return 0;
}
