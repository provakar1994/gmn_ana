/*
  This macro will fit dx distribution in various different ways.
  Plan:
  1. MC signal, no bg
  2. MC signal + bg from data
  3. MC signal + poly bg
  -------
  P. Datta Created 05-02-2023
*/

/*
  To-Do:
  1. A gif for the fits
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
  h->SetMarkerColor(kRed);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(22);
  h->SetLineColor(kRed);
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

void customize_hcut(TH1F* h)
{
  h->SetLineWidth(2);
  h->SetLineColor(kBlack);
  h->SetStats(0);
  h->GetXaxis()->CenterTitle(true);
}

void customize_hsummary(TH1F* h, std::vector<std::string> const & lcuts)
{
  gStyle->SetErrorX(0);
  gPad->SetGridy();
  h->SetStats(0);
  h->SetMarkerStyle(20);
  h->SetMarkerColor(2);
  h->SetLineColor(2);
  for (int i=0; i<lcuts.size(); i++) h->GetXaxis()->SetBinLabel(i+1,lcuts[i].c_str());
  h->LabelsOption("v","X");
}

void AddCutToLegend(TLegend *leg, std::string cut) {
  TLegendEntry* legE = leg->AddEntry((TObject*)0,Form("%s",cut.c_str()),"");
  legE->SetTextColor(2);
}

double total_fit (double * x, double * par) {
  FitFn *ffn = new FitFn(6);
  return ffn->ffn_gaus(x,&par[0]) + ffn->ffn_poly(x,&par[3]);
}

int fit_dx (const char *configfilename, 
	    bool is_elastic = 1) // 1=>Yes, 0=>QE 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

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

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/%s%s_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfprefix.c_str(),key,conf,sbsmag,model,pass));
  ROOT::RDataFrame simu_rdf("Tout",Form("simulation/siout/%s%s_ana_%s_sbs%d_sbs%dp_model%d.root",sfprefix.c_str(),key,gen.c_str(),conf,sbsmag,model));
  ROOT::RDataFrame inel_rdf("Tout",Form("simulation/siout/%sinel_%s_ana_g4sbs_sbs%d_sbs%dp_model%d.root",infprefix.c_str(),key,conf,sbsmag,model));

  // Applying cuts
  std::string cuts_for_signal_data = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_data");
  std::string cuts_for_signal_simu = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_simu");
  std::string cuts_for_bg_data = jmgr->GetValueFromSubKey_str(key,"cuts_for_bg_data");
  std::string cuts_for_bg_simu = jmgr->GetValueFromSubKey_str(key,"cuts_for_bg_simu");
  std::string coinT_cut = jmgr->GetValueFromSubKey_str(key,"coinT_cut");
  bool is_vary_cut = jmgr->GetValueFromSubKey<int>(key,"is_vary_cut");
  int Opoly = jmgr->GetValueFromSubKey<int>(key,"Order_of_poly_bg");
  auto data_rdf_filtered = data_rdf.Filter(cuts_for_signal_data);
  auto simu_rdf_filtered_1 = simu_rdf.Filter(cuts_for_signal_simu);
  auto inel_rdf_filtered_1 = inel_rdf.Filter(cuts_for_bg_simu);
  // Reading in offsets for dx peak position in MC for both p and n
  double dx_offset_p = jmgr->GetValueFromSubKey<double>(key,"dx_offset_MC_for_p");
  double dx_offset_n = !is_elastic ? jmgr->GetValueFromSubKey<double>(key,"dx_offset_MC_for_n") : 0;
  std::string dx_shifted_p = "dx+" + std::to_string(dx_offset_p);
  std::string dx_shifted_n = "dx+" + std::to_string(dx_offset_n);
  auto simu_rdf_filtered = simu_rdf_filtered_1.Define("dx_shifted_p",dx_shifted_p.c_str()).Define("dx_shifted_n",dx_shifted_n.c_str());
  auto inel_rdf_filtered = inel_rdf_filtered_1.Define("dx_shifted_p",dx_shifted_p.c_str()).Define("dx_shifted_n",dx_shifted_n.c_str());
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
  //TH1F *h_dxHCAL_bg_inel = (TH1F*)inel_rdf_filtered.Histo1D({"h_dxHCAL_bg_inel","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();  
  TH1F *h_dxHCAL_bg_inel_p = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_bg_inel_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight")->Clone();
  TH1F *h_dxHCAL_bg_inel_n;
  if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_bg_inel_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight")->Clone();
  TH1F *h_dxHCAL_bg_inel = (TH1F*)h_dxHCAL_bg_inel_p->Clone(); h_dxHCAL_bg_inel->Add(h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n);

  // defining output files
  std::string filebase = jmgr->GetValueFromSubKey_str(key,"output_filebase");
  char const * outdir = !is_vary_cut ? "pdout/" : "pdout/sysstdy/";
  TString outFile = Form("%s%s_%s_pass%d_%s_sbs%d_sbs%dp_model%d.root",outdir,filebase.c_str(),key,pass,gen.c_str(),conf,sbsmag,model);
  TString outPlot = outFile; outPlot.ReplaceAll(".root",".pdf");
  TString outData = outFile; outData.ReplaceAll(".root",".csv"); ofstream outdata; outdata.open(outData);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // Fits
  vector<double> dx_fit_range; jmgr->GetVectorFromSubKey<double>(key,"dx_fit_range",dx_fit_range);
  vector<double> reject_points; jmgr->GetVectorFromSubKey<double>(key,"SB_reject_points",reject_points);

  // Cut variation **************
  // Although not necessary, keeping this loop separate gives more control and clarity
  // forming the cuts
  std::string param_to_vary = jmgr->GetValueFromSubKey_str(key,"param_to_vary");
  vector<double> cut_range; jmgr->GetVectorFromSubKey<double>(key,"cut_iter_min_width",cut_range);
  vector<double> h_cut_param; jmgr->GetVectorFromSubKey<double>(key,"h_cut_param",h_cut_param);
  double min = cut_range[1], width = cut_range[2]; 
  int iter = is_vary_cut ? (int)cut_range[0] : 1;
  // choosing cut variation style
  int cut_vary_style = jmgr->GetValueFromSubKey<int>(key,"cut_vary_style");
  // 0 -> scan from low to high with fixed width
  // 1 -> scan around a fixed mean with increasing width. NOTE: cut_range[1] = mean, in this case
  // 2 -> increase threshold by a fixed amount
  double low = cut_vary_style==1 ? min-width : min;
  double high = cut_vary_style==2 ? h_cut_param[2] : min+width; 
  std::vector<double> minval, maxval;
  std::vector<std::string> cuts, cuts_p, cuts_n, cuts_2;
  for (int i=0; i<iter; i++) {
    minval.push_back(low); maxval.push_back(high);
    std::string cut, cut_2;
    char low_buff[20]; std::snprintf(low_buff,20,"%.2f",low); std::string low_str(low_buff);
    char high_buff[20]; std::snprintf(high_buff,20,"%.2f",high); std::string high_str(high_buff);
    if (cut_vary_style==0 || cut_vary_style==1) {
      cut = param_to_vary+">"+low_str+"&&"+param_to_vary+"<="+high_str; 
      cut_2 = low_str+"<"+param_to_vary+"<="+high_str;
    } 
    else if (cut_vary_style==2) { cut = param_to_vary+">"+low_str; cut_2 = cut; }
    cuts.push_back(cut); cuts_2.push_back(cut_2);
    //std::cout << cut << "\n";
    std::string cut_p = cut + "&&mc_fnucl==1"; cuts_p.push_back(cut_p);
    std::string cut_n = cut + "&&mc_fnucl==0"; cuts_n.push_back(cut_n);
    // update high and low
    if (cut_vary_style==0) { low = high; high += width; }
    else if (cut_vary_style==1) { low -= width; high += width; }
    else if (cut_vary_style==2) { low += width; high = h_cut_param[2]; }
  }

  // various outputs
  outdata << "cut,min,max,chi20,NDF0,R0,R0err,B0,B0err,chi21,NDF1,R1,R1err,B1,B1err,chi22,NDF2,R2,R2err,B2,B2err,"
	  << "chi23,NDF3,R3,R3err,B3,B3err,chi24,NDF4,R4,R4err,B4,B4err\n";
  TString outGIF = outFile; outGIF.ReplaceAll(".root",".gif");

  // summary histo
  TH1F *htemp = new TH1F("htemp","",iter,-0.5,iter-0.5);

  // draawing cut histo
  TH1F *hcut = (TH1F*)data_rdf_filtered.Histo1D({"hcut","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]},param_to_vary)->Clone();
  customize_hcut(hcut); hcut->SetTitle(Form("%s {%s}",param_to_vary.c_str(),cuts_for_signal_data.c_str()));
  // Canvas to plot cut region
  TCanvas *cCut = util_pd::TC("cCut",1,1);
  if (!is_vary_cut) cCut->Delete();
  for (int i=0; i<iter; i++) {
    if (is_vary_cut) {
      h_dxHCAL_data = (TH1F*)data_rdf_filtered.Filter(cuts[i]).Histo1D({"h_dxHCAL_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      h_dxHCAL_data_CT = (TH1F*)data_rdf_filtered.Filter(coinT_cut.c_str()).Histo1D({"h_dxHCAL_data_CT","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      h_dxHCAL_simu_p = (TH1F*)simu_rdf_filtered.Filter(cuts_p[i]).Histo1D({"h_dxHCAL_simu_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight")->Clone();
      if (!is_elastic) h_dxHCAL_simu_n = (TH1F*)simu_rdf_filtered.Filter(cuts_n[i]).Histo1D({"h_dxHCAL_simu_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight")->Clone();

      // drawing cut histos
      cCut->SetTickx(); cCut->SetTicky(); cCut->cd();
      hcut->Draw();
      util_pd::PlotCutRegion(minval[i],maxval[i]);
      TLegend *lCut = new TLegend(0.10,0.80,0.30,0.9);
      lCut->SetTextFont(42);
      lCut->AddEntry(hcut,Form("%s",param_to_vary.c_str()),"l");
      AddCutToLegend(lCut,cuts_2[i].c_str());
      lCut->Draw();
      cCut->Update(); cCut->Write(); if (is_vary_cut&&i==0) cCut->SaveAs(Form("%s[",outPlot.Data())); 
      cCut->SaveAs(Form("%s",outPlot.Data())); 
      cCut->SaveAs(Form("%s+150",outGIF.Data()));
    }
    // ***

    /*######################################
      ## Fitting elastic dx distributions ##
      ###################################### */
    std::vector<double> R_vals,Rerr_vals,chi2,NDF,B_vals,Berr_vals;
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
      l0->AddEntry(ho[1],"MC Signal","lep");
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
      //vector<double> reject_points{-1.3,-0.15};
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
      // Canvas 0 : Fitting data/MC w/o any background
      TCanvas *c0 = util_pd::TC("c0",1,1);
      c0->cd(); //gStyle->SetOptFit(1);
      vector<TH1F*> ho;
      TF1 *f0 = fit::fit_2hs_nbg_THI(dx_fit_range,
				     h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,
				     ho);
      // grabbing fit params for future use
      chi2.push_back(f0->GetChisquare()); NDF.push_back(f0->GetNDF());
      R_vals.push_back(f0->GetParameter(1)); Rerr_vals.push_back(f0->GetParError(1));
      B_vals.push_back(0); Berr_vals.push_back(0);
      //ho[0]->Draw("E"); customize_data(ho[0]); customize_dx(ho[0]);
      h_dxHCAL_data->Draw("E"); customize_data(h_dxHCAL_data);
      ho[1]->Draw("same"); customize_hs(ho[1]); //customize_dx(ho[1]);
      //f0->Draw("same");
      TLegend *l0 = new TLegend(0.10,0.78,0.30,0.9);
      l0->SetTextFont(42);
      l0->AddEntry(h_dxHCAL_data,"Data","p");
      //l0->AddEntry(f0,"Fit","l");
      l0->AddEntry(ho[1],"MC Signal","p");
      if (is_vary_cut) AddCutToLegend(l0,cuts_2[i].c_str());
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
      // grabbing fit params for future use
      chi2.push_back(f1->GetChisquare()); NDF.push_back(f1->GetNDF());
      R_vals.push_back(f1->GetParameter(1)); Rerr_vals.push_back(f1->GetParError(1));
      B_vals.push_back(f1->GetParameter(2)); Berr_vals.push_back(f1->GetParError(2));
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
      if (is_vary_cut) AddCutToLegend(l1,cuts_2[i].c_str());
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
      if (is_vary_cut) {
	htemp->SetBinContent(i+1,f2->GetParameter(1));
	htemp->SetBinError(i+1,f2->GetParError(1));
      }
      // grabbing fit params for future use
      chi2.push_back(f2->GetChisquare()); NDF.push_back(f2->GetNDF());
      R_vals.push_back(f2->GetParameter(1)); Rerr_vals.push_back(f2->GetParError(1));
      B_vals.push_back(f2->GetParameter(2)); Berr_vals.push_back(f2->GetParError(2));
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
      if (is_vary_cut) AddCutToLegend(l2,cuts_2[i].c_str());
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
      // grabbing fit params for future use
      chi2.push_back(f3->GetChisquare()); NDF.push_back(f3->GetNDF());
      R_vals.push_back(f3->GetParameter(1)); Rerr_vals.push_back(f3->GetParError(1));
      B_vals.push_back(f3->GetParameter(2)); Berr_vals.push_back(f3->GetParError(2));
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
      if (is_vary_cut) AddCutToLegend(l3,cuts_2[i].c_str());
      //l3->SetFillStyle(0); // makes legend box transparent
      l3->Draw();
      //
      // preparing the pad for residual
      //  
      p3[1]->cd();
      ho3[3]->Draw(); customize_residual(ho3[3]);
      // drawing a horizontal line at y = 0
      util_pd::DrawZeroLine(p3[1],dx_fit_range[0],dx_fit_range[1]);

      // ** --
      // Canvas 4 : Fitting w/ polynomial background using side band method
      // Steps: 1. Subrtact bg using sideband method, 2. Perform data/MC fit w/o bg
      TCanvas *c4 = util_pd::TC("c4",1,1); gStyleFitCanvas();
      c4->cd(); 
      // Let's perform the sideband fit first
      //vector<double> reject_points{-1.2,0.3};
      vector<TH1F*> hosb4;
      TF1* bgsb4 = fit::fit_1pbg_SB(dx_fit_range,
				    reject_points,
				    Opoly,
				    fit::GetFitParams(f3),
				    h_dxHCAL_data,
				    hosb4);
      hosb4[0]->Draw(); c4->Update(); 
      // grabbing statbox of the fitted histo
      TPaveStats *stsb4 = (TPaveStats*)hosb4[0]->FindObject("stats"); 
      // hosb4[0]->Draw(); customize_ht(hosb4[0]); customize_dx(hosb4[0]);
      // hosb4[1]->Draw("same ep"); customize_hs(hosb4[1]); customize_dx(hosb4[1]);
      // hosb4[2]->Draw("same"); customize_hbg(hosb4[2]); customize_dx(hosb4[2]);

      // Now, let's fit the bg subtracted signal using MC signals
      vector<TH1F*> ho4;
      TF1 *f4 = fit::fit_2hs_nbg_THI(dx_fit_range,
				     hosb4[1],h_dxHCAL_simu_p,h_dxHCAL_simu_n,
				     ho4);
      // grabbing fit params for future use
      chi2.push_back(f4->GetChisquare()); NDF.push_back(f4->GetNDF());
      R_vals.push_back(f4->GetParameter(1)); Rerr_vals.push_back(f4->GetParError(1));
      B_vals.push_back(bgsb4->GetParameter(0)); Berr_vals.push_back(bgsb4->GetParError(0));
      // converting fit fn to a hostogram
      TH1F *hf4 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f4,hf4);
      ho4[0]->Draw(); c4->Update(); 
      // hosb4[0]->Draw("same"); customize_ht(hosb4[0]); customize_dx(hosb4[0]);
      // hosb4[1]->Draw("same HIST"); customize_hs(hosb4[1]); customize_dx(hosb4[1]);
      // hosb4[2]->Draw("same HIST"); customize_hbg(hosb4[2]); customize_dx(hosb4[2]);

      // grabbing statbox of the fitted histo
      TPaveStats *st4 = (TPaveStats*)ho4[0]->FindObject("stats");
      // getting pads for pull plot
      std::vector<TPad*> p4 = util_pd::GetPadsForPullPlot(c4);
      //
      // preparing the pad for data/MC fit
      //
      p4[0]->cd();
      // drawing all the histograms
      hosb4[0]->Draw(); customize_ht(hosb4[0]); customize_dx(hosb4[0]);
      hosb4[2]->Draw("same"); customize_hbg(hosb4[2]); customize_dx(hosb4[2]);    
      // h_dxHCAL_data->Draw("E"); customize_data(h_dxHCAL_data);
      hosb4[1]->Draw("same E"); customize_data(hosb4[1]); hosb4[1]->SetStats(0);
      hf4->Draw("same HIST"); customize_gfit(hf4);
      ho4[3]->Draw("same HIST"); customize_psig(ho4[3]);
      ho4[4]->Draw("same HIST"); customize_nsig(ho4[4]);
      //ho4[2]->Draw("same HIST"); customize_hbg(ho4[2]);
      // redrawing the stat boxes
      stsb4->SetX1NDC(0.62); stsb4->SetX2NDC(0.9); stsb4->SetY2NDC(0.9);
      stsb4->Draw("same");
      st4->SetX1NDC(0.62); st4->SetX2NDC(0.9); st4->SetY2NDC(0.6);
      st4->Draw("same");
      // drawing a legend
      TLegend *l4=new TLegend(0.10,0.6,0.36,0.9);
      l4->SetTextFont(42);
      l4->AddEntry(hosb4[0],"Data","l");
      l4->AddEntry(hosb4[2],Form("Bg w/ SB fit (poly. ord. %d)",Opoly),"lf");
      //l4->AddEntry(f4,"Fit (MC + poly. bg.)","l");
      l4->AddEntry(hosb4[1],"Data (bg subtracted)","p");
      l4->AddEntry(hf4,"Fit (QE MC + bg.)","lf");
      l4->AddEntry(ho4[3],"p signal (from MC)","lf");
      l4->AddEntry(ho4[4],"n signal (from MC)","lf");
      //l4->AddEntry(ho4[2],"Residual","p");
      if (is_vary_cut) AddCutToLegend(l4,cuts_2[i].c_str());
      //l4->SetFillStyle(0); // makes legend box transparent
      l4->Draw();
      // //
      // // preparing the pad for residual
      // p4[1]->cd();
      // ho4[2]->Draw(); customize_residual(ho4[2]);
      // // drawing a horizontal line at y = 0
      // util_pd::DrawZeroLine(p4[1],dx_fit_range[0],dx_fit_range[1]);

      // Writing out fit parameters
      std::cout << "\n--- Reporting fit params ---\n";
      std::cout << "cut,min,max,R0,R0err,R1,R1err,R2,R2err,R3,R3err,R4,R4err\n";
      std::cout << cuts_2[i] << "," << minval[i] << "," << maxval[i] << ",";
      outdata << cuts_2[i] << "," << minval[i] << "," << maxval[i] << ",";
      for(size_t i=0; i < R_vals.size(); i++){
	std::cout << R_vals[i] << "," << Rerr_vals[i] << ",";
	//outdata << R_vals[i] << "," << Rerr_vals[i] << ",";
	outdata << Form("%.1f,%.1f,%.4f,%.4f,%.4f,%.4f,",chi2[i],NDF[i],R_vals[i],Rerr_vals[i],B_vals[i],Berr_vals[i]);
      }
      outdata << "\n";
      std::cout << "\n------\n";

      // further customization of the data histo 
      h_dxHCAL_data->SetStats(0); //[IMPORTANT!]
      h_dxHCAL_data->SetTitle(Form("dx {%s}",cuts_for_signal_data.c_str()));
      h_dxHCAL_data->GetYaxis()->SetRangeUser(0,h_dxHCAL_data->GetMaximum()*1.1);

      // writing out the canvases
      c0->Update(); c0->Write(); if (!is_vary_cut&&i==0) c0->SaveAs(Form("%s[",outPlot.Data())); 
      c0->SaveAs(Form("%s",outPlot.Data())); 
      c1->Update(); c1->Write(); c1->SaveAs(Form("%s",outPlot.Data())); 
      c2->Update(); c2->Write(); c2->SaveAs(Form("%s",outPlot.Data())); 
      c3->Update(); c3->Write(); c3->SaveAs(Form("%s",outPlot.Data())); 
      c4->Update(); c4->Write(); c4->SaveAs(Form("%s",outPlot.Data())); 
      if (i==iter-1) c4->SaveAs(Form("%s]",outPlot.Data())); 
    } // QE
  } // for, cut vairation

  if (is_vary_cut) {
    // creating infinite loop gif
    cCut->SaveAs(Form("%s++5",outGIF.Data())); 
    // writing additional memory to file
    // Canvas : Plotting cut variation summary
    TCanvas *cGist = util_pd::TC("cGist",1,1); cGist->cd();
    cGist->SetBottomMargin(0.3);
    customize_hsummary(htemp,cuts_2);
    htemp->Draw(); htemp->Write();
    cGist->Write();
  }

  // reporting output files
  std::cout << "------" << std::endl;
  std::cout << " Fit params  : " << outData << std::endl;
  std::cout << " Summary plots  : " << outPlot << std::endl;
  std::cout << " Output ROOT file  : " << outFile << std::endl;
  std::cout << "------" << std::endl;

  sw->Stop();
  std::cout << "CPU time = " << sw->CpuTime() << "s. Real time = " << sw->RealTime() << "s.\n\n";

  return 0;
}
