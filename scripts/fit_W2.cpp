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

void customize_sig(TH1F* h)
{
  h->SetLineColor(kGreen+2);
  h->SetLineStyle(8);
  h->SetLineWidth(3);
  //***
  h->SetFillColor(kGreen+2);
  h->SetFillColorAlpha(kGreen+2,0.3);
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

void gStyleFitCanvas() 
{
  gStyle->SetOptStat("e"); gStyle->SetOptFit(1); 
  gStyle->SetErrorX(0);
}
// ** ^ New

void customize_W2(TH1F* h)
{
  h->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
}

// Path of a simulation analysis output (Tout) for a given file prefix.
std::string SimuFileName(std::string const &prefix, char const *process, std::string const &gen,
			 int conf, int sbsmag, int model)
{
  std::string pre = prefix.empty() ? "" : prefix + "_";
  return Form("simulation/siout/%s%s_ana_%s_sbs%d_sbs%dp_model%d.root",
	      pre.c_str(), process, gen.c_str(), conf, sbsmag, model);
}

// MC W2 histograms from a weighted mix of simulation files (a single file is a mix with weight 1).
// Each file is filled with its per-event "weight", which is already normalized by that sample's own
// number of generated events, so the mix weights are simply the fractions of the mix (e.g. the
// charge fractions of the shield-on and shield-off periods at SBS-11).
// One histogram is returned per entry of extra_cuts (applied on top of base_cut), all filled in one
// pass over each file. W2 is shifted by W2_off before filling.
std::vector<TH1F*> MixedSimuW2Hists(std::vector<std::string> const &files, std::vector<double> const &mix_weights,
				    std::string const &base_cut, std::vector<std::string> const &extra_cuts,
				    std::vector<std::string> const &names, double W2_off, std::vector<double> const &hbin)
{
  std::vector<TH1F*> hsum;
  for (size_t ih=0; ih<names.size(); ih++) {
    hsum.push_back(new TH1F(names[ih].c_str(), "", int(hbin[0]), hbin[1], hbin[2]));
    hsum.back()->Sumw2();
  }
  std::string W2_shifted = "W2+" + std::to_string(W2_off);
  for (size_t ifile=0; ifile<files.size(); ifile++) {
    ROOT::RDataFrame rdf("Tout", files[ifile]);
    auto rdf_cut = rdf.Filter(base_cut).Define("W2_shifted", W2_shifted.c_str());
    // book all histograms first so the file is read only once
    std::vector<ROOT::RDF::RResultPtr<TH1D>> booked;
    for (size_t ih=0; ih<names.size(); ih++)
      booked.push_back(rdf_cut.Filter(extra_cuts[ih])
		       .Histo1D({Form("%s_f%zu",names[ih].c_str(),ifile), "", int(hbin[0]), hbin[1], hbin[2]}, "W2_shifted", "weight"));
    for (size_t ih=0; ih<names.size(); ih++) hsum[ih]->Add(booked[ih].GetPtr(), mix_weights[ifile]);
  }
  return hsum;
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

// void customize_hbg(TH1F* h) 
// {
//   h->SetMarkerColor(kRed);
//   h->SetMarkerStyle(29);
//   h->SetLineColor(kRed);
// }

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

  // MC signal sample: by default the single file given by simu_file_prefix. Optionally a weighted
  // mix of several files (e.g. SBS-11: Al-shield + no-shield MC, weighted by the charge fractions):
  //   "simu_mix_prefixes": ["<prefix 1>", "<prefix 2>"], "simu_mix_weights": [w1, w2]
  // If present, these override simu_file_prefix.
  std::vector<std::string> simu_prefixes; jmgr->GetVectorFromSubKey<std::string>(key,"simu_mix_prefixes",simu_prefixes);
  std::vector<double> simu_mix_weights; jmgr->GetVectorFromSubKey<double>(key,"simu_mix_weights",simu_mix_weights);
  if (simu_prefixes.empty()) { simu_prefixes = {sfprefix}; simu_mix_weights = {1.}; }
  if (simu_prefixes.size() != simu_mix_weights.size()) {
    std::cerr << "[fit_W2] simu_mix_prefixes and simu_mix_weights must have the same length!\n";
    return 1;
  }
  std::vector<std::string> simu_files;
  for (auto const &p : simu_prefixes) simu_files.push_back(SimuFileName(p,process,gen,conf,sbsmag,model));

  // MC background (inelastic) file for canvas 3, relative to simulation/siout/
  std::string inel_file = "simulation/siout/" + jmgr->GetValueFromSubKey_str(key,"inel_file");

  // defining output files
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d.root",filebase.c_str(),conf,sbsmag,model);

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/%s%s_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfprefix.c_str(),process,conf,sbsmag,model,pass));
  ROOT::RDataFrame inel_rdf("Tout",inel_file);
  std::cout << "\nMC signal sample:\n";
  for (size_t i=0; i<simu_files.size(); i++) std::cout << Form("  weight %.3f : %s\n", simu_mix_weights[i], simu_files[i].c_str());
  std::cout << "MC background: " << inel_file << "\n";

  // Applying cuts
  std::string cuts_for_inclusive_W2 = jmgr->GetValueFromSubKey_str(key,"cuts_for_inclusive_W2");
  std::string cuts_for_signal_data = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_data");
  std::string cuts_for_signal_simu = jmgr->GetValueFromSubKey_str(key,"cuts_for_signal_simu");
  std::string cuts_for_bg_data = jmgr->GetValueFromSubKey_str(key,"cuts_for_bg_data");
  std::string cuts_for_bg_simu = jmgr->GetValueFromSubKey_str(key,"cuts_for_bg_simu");
  std::string coinT_cut = jmgr->GetValueFromSubKey_str(key,"coinT_cut");
  int Opoly = jmgr->GetValueFromSubKey<int>(key,"Order_of_poly_bg");
  auto data_rdf_filtered = data_rdf.Filter(cuts_for_signal_data);
  auto inel_rdf_filtered_1 = inel_rdf.Filter(cuts_for_bg_simu);
  auto bg_rdf_filtered = data_rdf.Filter(cuts_for_bg_data);
  // Applying offsets
  std::vector<double> W2Off_range; jmgr->GetVectorFromSubKey<double>(key,"vary_W2Off_ranges",W2Off_range);
  bool is_vary_W2Off = W2Off_range[0];
  // define the offset values
  double W2_off_MC_sig = is_vary_W2Off ? 0 : jmgr->GetValueFromSubKey<double>(key,"W2_peakOff_for_MC_sig");
  double W2_off_MC_bg = is_vary_W2Off ? 0 : jmgr->GetValueFromSubKey<double>(key,"W2_peakOff_for_MC_bg");
  std::string W2_shifted_MC_bg = "W2+" + std::to_string(W2_off_MC_bg);
  auto inel_rdf_filtered = inel_rdf_filtered_1.Define("W2_shifted",W2_shifted_MC_bg.c_str());


  // Creating important histograms
  vector<double> h_W2; jmgr->GetVectorFromSubKey<double>(key,"h_W2",h_W2);
  TH1F *h_W2_data_raw = (TH1F*)data_rdf.Filter(cuts_for_inclusive_W2).Histo1D({"h_W2_data_raw","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2")->Clone();
  TH1F *h_W2_data = (TH1F*)data_rdf_filtered.Histo1D({"h_W2_data","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2")->Clone();
  TH1F *h_W2_data_CT = (TH1F*)data_rdf_filtered.Filter(coinT_cut.c_str()).Histo1D({"h_W2_data_CT","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2")->Clone();
  // MC signal: all, proton only, neutron only (weighted mix if several files are given)
  std::vector<TH1F*> h_simu = MixedSimuW2Hists(simu_files, simu_mix_weights, cuts_for_signal_simu,
					       {"true","mc_fnucl==1","mc_fnucl==0"},
					       {"h_W2_simu","h_W2_simu_p","h_W2_simu_n"}, W2_off_MC_sig, h_W2);
  TH1F *h_W2_simu = h_simu[0];
  TH1F *h_W2_simu_p = h_simu[1];
  TH1F *h_W2_simu_n = h_simu[2];
  TH1F *h_W2_inel = (TH1F*)inel_rdf_filtered.Histo1D({"h_W2_inel","",int(h_W2[0]),h_W2[1],h_W2[2]},"W2_shifted","weight")->Clone();
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
  TLegend *l0=new TLegend(0.10,0.77,0.3,0.9);
  l0->SetTextFont(42);
  l0->AddEntry(ho[0],"Data","l");
  l0->AddEntry(f0,"Fit","l");
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
  // std::vector<double> xOff_range = {0,-0.1};
  // TF1 *f1 = fit::fit_1hs_1pbg_THI_xOffVary(W2_fit_range,
  // 					   h_W2_data_raw,h_W2_simu,Opoly,xOff_range,
  // 					   ho1);
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
  l1->AddEntry(f1,"Fit","l");
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
  // std::vector<double> xOff_range = {0,-0.05};
  // TF1 *f2 = fit::fit_1hs_1pbg_THI_xOffVary(W2_fit_range,
  // 					   h_W2_data_raw,h_W2_data,Opoly,xOff_range,
  // 					   ho2);
  ho2[0]->Draw(); customize_ht(ho2[0]); customize_W2(ho2[0]);
  ho2[1]->Draw("same"); customize_hs(ho2[1]); customize_W2(ho2[1]);
  ho2[2]->Draw("same"); customize_hbg(ho2[2]);
  TLegend *l2=new TLegend(0.10,0.77,0.3,0.9);
  l2->SetTextFont(42);
  l2->AddEntry(ho2[0],"Data","l");
  l2->AddEntry(f2,"Fit","l");
  l2->AddEntry(ho2[1],"Signal (from data)","lep");
  l2->AddEntry(ho2[2],Form("Bg. (%d^{th} order poly.)",Opoly),"lep");
  l2->Draw();
  std::cout << Form("\nTotal # elastics (from fitted data signal): %d\n",
		    int(ho2[1]->Integral(ho2[1]->FindBin(W2_fit_range[0]),ho2[1]->FindBin(W2_fit_range[1]))));

  // Canvas 3 : Fitting w/ signal and background from MC
  TCanvas *c3 = util_pd::TC("c3",1,1); gStyleFitCanvas();
  c3->cd(); gStyle->SetOptFit(1);
  vector<TH1F*> ho3;
  TF1 *f3;
  if (is_vary_W2Off) {
    f3 = fit::fit_1hs_1hbg_THI_xOffVary(W2_fit_range,
					h_W2_data_raw,h_W2_simu,h_W2_inel,W2Off_range,
					ho3);
  } else {
    f3 = fit::fit_1hs_1hbg_THI(W2_fit_range,
			       h_W2_data_raw,h_W2_simu,h_W2_inel,
			       ho3);
  }
  // converting fit fn to a hostogram
  TH1F *hf3 = (TH1F*)h_W2_data_raw->Clone(); util_pd::TF1toTH1F(f3,hf3);
  ho3[0]->Draw(); c3->Update(); 
  // grabbing statbox of the fitted histo
  TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
  
  h_W2_data_raw->Draw("E"); customize_data(h_W2_data_raw);
  hf3->Draw("same HIST"); customize_gfit(hf3); //customize_W2(ho3[0]);
  ho3[3]->Draw("same HIST"); customize_sig(ho3[3]); //customize_W2(ho3[3]);
  ho3[4]->Draw("same HIST"); customize_hbg(ho3[4]); //customize_W2(ho3[4]);
  // redrawing the stat box
  st3->SetX1NDC(0.62); st3->SetX2NDC(0.9); st3->SetY2NDC(0.9);
  st3->Draw("same");
  // drawing the legend
  TLegend *l3=new TLegend(0.10,0.77,0.3,0.9);
  l3->SetTextFont(42);
  l3->AddEntry(h_W2_data_raw,"Data","ep");
  l3->AddEntry(hf3,"Fit","lf");
  l3->AddEntry(ho3[3],"Signal (from MC)","lf");
  l3->AddEntry(ho3[4],"Bg. (from MC)","lf");
  l3->Draw();
  // ---

  // further customization of the data histo 
  h_W2_data_raw->SetStats(0); //[IMPORTANT!]
  h_W2_data_raw->GetXaxis()->SetTitle("W^{2} (GeV^{2})");
  h_W2_data_raw->GetXaxis()->CenterTitle();

  // save all canvases to one PDF: <filebase>_sbs<conf>_sbs<mag>p_model<model>.pdf
  TString pdfFile = outFile; pdfFile.ReplaceAll(".root",".pdf");
  c0->Print(pdfFile + "[");
  for (TCanvas *c : {c0, c1, c2, c3}) c->Print(pdfFile);
  c3->Print(pdfFile + "]");
  std::cout << "Canvases saved to " << pdfFile << "\n";

  return 0;
}
