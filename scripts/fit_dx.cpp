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
#include "TList.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"

#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>

#include "../include/gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

bool temp = 0;

//______________________________________________________________________________
void CalcRnumBinEdges(int rmin, int rmax, bool debug, std::vector<double> &binEdges) {
  /* Calculates explicit bin edges for run histo */
  int nbins = rmax - rmin;
  for (int i = rmin; i <= rmax; ++i) {
    binEdges.push_back((double)i-0.5);
  }
  if (debug) util_pd::PrintVector(binEdges);
}

//______________________________________________________________________________
void gStyleFitCanvas() 
{
  gStyle->SetOptStat("e"); gStyle->SetOptFit(1); 
  gStyle->SetErrorX(0);
}

//______________________________________________________________________________
void customize_dx(TH1F* h)
{
  h->GetXaxis()->SetTitle("x_{HCAL}^{obs} - x_{HCAL}^{exp} (m)");
}

//______________________________________________________________________________
void customize_ht(TH1F* h) 
{
  h->SetLineColor(kBlue);
  h->SetLineWidth(2);
}

//______________________________________________________________________________
void customize_hs(TH1F* h) 
{
  h->SetMarkerColor(kRed);
  h->SetMarkerSize(0.8);
  h->SetMarkerStyle(22);
  h->SetLineColor(kRed);
}

//______________________________________________________________________________
void customize_hcut(TH1F* h)
{
  h->SetLineWidth(2);
  h->SetLineColor(kBlack);
  h->SetStats(0);
  h->GetXaxis()->CenterTitle(true);
}

//______________________________________________________________________________
void customize_hcut_p(TH1F* h)
{
  customize_hcut(h);
  h->SetLineColor(kBlue);
}

//______________________________________________________________________________
void customize_hcut_n(TH1F* h)
{
  customize_hcut(h);
  h->SetLineColor(kGreen+2);
}

//______________________________________________________________________________
void customize_h2fiduCut(TH2F* h2, char const * nORp, char const * DataORMC, double sbs_kick) 
{
  h2->SetTitle(Form("%s Envelope (%s)",nORp,DataORMC));
  h2->GetXaxis()->SetTitle("y_{HCAL}^{exp} (m)");
  std::string ytitle = sbs_kick==0 ? "x_{HCAL}^{exp} (m)" : "x_{HCAL}^{exp} - " + std::to_string(sbs_kick) + " (m)"; 
  h2->GetYaxis()->SetTitle(ytitle.c_str());
}

//______________________________________________________________________________
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

//______________________________________________________________________________
void AddCutToLegend(TLegend *leg, std::string cut) {
  TLegendEntry* legE = leg->AddEntry((TObject*)0,Form("%s",cut.c_str()),"");
  legE->SetTextColor(2);
}

//______________________________________________________________________________
void customize_text(TText *tl) {
  tl->SetTextFont(42);
  tl->SetTextSize(0.04);
  tl->SetTextColor(kRed);
}

//______________________________________________________________________________
void FurtherCustoizeDataHisto(TH1F *h_dxHCAL_data, std::string const &cuts_for_signal_data) {
  h_dxHCAL_data->SetStats(0);
  h_dxHCAL_data->SetTitle(Form("dx {%s}",cuts_for_signal_data.c_str()));
  //h_dxHCAL_data->GetYaxis()->SetMaxDigits(3);
  h_dxHCAL_data->GetYaxis()->SetRangeUser(-5,h_dxHCAL_data->GetMaximum()*1.1);
}

// void AddFiduCutToLegend(TLegend *leg, std::vector<double> hcal_AR, std::vector<double> hcal_SM) {
//   TLegendEntry* legE0 = leg->AddEntry((TObject*)0,"HCAL Boundary","");
//   legE0->SetTextColor(kGreen+2);
//   TLegendEntry* legE0 = leg->AddEntry((TObject*)0,Form("%.1fb(x), %.1fb(y)",),"");
// }

//______________________________________________________________________________
void PlotFiduCut(int pass, std::vector<double> hcal_AR, std::vector<double> hcal_SM) {
  std::vector<double> hcal_area = cut::hcal_active_area_data(0,0,pass); 
  //std::vector<double> hcal_AR = cut::hcal_active_area_data(AR_w[0],AR_w[1],pass); 
  //std::vector<double> hcal_SM = cut::hcal_safety_margin(SM_w[0],SM_w[1],SM_w[2],hcal_AR);
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  util_pd::DrawArea(hcal_SM,4,4,9);
}

//______________________________________________________________________________
double total_fit (double * x, double * par) {
  FitFn *ffn = new FitFn(6);
  return ffn->ffn_gaus(x,&par[0]) + ffn->ffn_poly(x,&par[3]);
}

//______________________________________________________________________________
int fit_dx (const char *configfilename, 
	    bool is_elastic = 1) // 1=>Yes, 0=>QE 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  //TH1::SetDefaultSumw2();

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  char const * key = is_elastic ? "elas" : "qelas";

  // Get the keys of the JSON object
  // std::vector<std::string> subkeys; jmgr->GetSubKeys(key,subkeys);
  // util_pd::PrintVector(subkeys);

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
  auto simu_rdf_filtered = simu_rdf.Filter(cuts_for_signal_simu);
  auto inel_rdf_filtered = inel_rdf.Filter(cuts_for_bg_simu);
  auto bg_data_rdf_filtered = data_rdf.Filter(cuts_for_bg_data);

  // Forming a custom Fiducial Cut
  bool use_custom_fiduCut = jmgr->GetValueFromSubKey<int>(key,"use_custom_fiduCut");
  bool apply_to_bg_data = jmgr->GetValueFromSubKey<int>(key,"apply_to_bg_data");
  bool apply_to_bg_simu = jmgr->GetValueFromSubKey<int>(key,"apply_to_bg_simu");
  std::vector<double> AR_w; jmgr->GetVectorFromSubKey<double>(key,"AR_width_x_y",AR_w);
  double sbs_kick = jmgr->GetValueFromSubKey<double>(key,"sbs_kick");
  std::vector<double> SM_w; jmgr->GetVectorFromSubKey<double>(key,"SM_width_xp_xn_y",SM_w);
  std::string target = is_elastic ? "LH2" : "LD2";  
  std::vector<double> hcal_area = cut::hcal_active_area_data(0,0,pass); 
  std::vector<double> hcal_AR = cut::hcal_active_area_data(AR_w[0],AR_w[1],pass); 
  std::vector<double> hcal_SM = cut::hcal_safety_margin(SM_w[0],SM_w[1],SM_w[2],hcal_AR);

  // Calculating R_MC
  SBSconfig sbsconf(conf,sbsmag);
  double Rp_MC = kine::sigmaBorn_ratio_MC(sbsconf,"np");
  std::cout << "\nRp_MC = " << Rp_MC << "\n\n";

  // Now its time to define new columns
  // 1. to account for dx_p and dx_n shifts
  // 2. to account for xHCAL_exp shift for proton tracks
  // First, let's check if the user wanna vary x offsets as free parameters or not
  std::vector<double> pnXOff_range; jmgr->GetVectorFromSubKey<double>(key,"vary_pn_xOff_ranges",pnXOff_range);
  bool is_vary_pnXOff = pnXOff_range[0]; 
  // define the offset values
  double dx_offset_p = is_vary_pnXOff ? 0 : jmgr->GetValueFromSubKey<double>(key,"dx_offset_MC_for_p");
  std::string dx_shifted_p = "dx+" + std::to_string(dx_offset_p);
  double dx_offset_n = (is_elastic||is_vary_pnXOff) ? 0 : jmgr->GetValueFromSubKey<double>(key,"dx_offset_MC_for_n");
  std::string dx_shifted_n = "dx+" + std::to_string(dx_offset_n);
  std::string xExp_shifted = "xHCAL_exp-" + std::to_string(sbs_kick);
  // define new columns 
  data_rdf_filtered = data_rdf_filtered
    .Define("xExp_shifted",xExp_shifted.c_str());
  simu_rdf_filtered = simu_rdf_filtered
    .Define("dx_shifted_p",dx_shifted_p.c_str())
    .Define("dx_shifted_n",dx_shifted_n.c_str())
    .Define("xExp_shifted",xExp_shifted.c_str());
  inel_rdf_filtered = inel_rdf_filtered
    .Define("dx_shifted_p",dx_shifted_p.c_str())
    .Define("dx_shifted_n",dx_shifted_n.c_str());
  /*
    Since we keep the data histogram untouched for fitting, I think we should 
    do the same for bg histogram from data. Hence, commenting out the following lines.
  */
  // bg_data_rdf_filtered = bg_data_rdf_filtered
  //   .Define("dx_shifted",dx_shifted_p.c_str());

  // defining output files
  std::string filebase = jmgr->GetValueFromSubKey_str(key,"output_filebase");
  char const * confmag = Form("sbs%dsbs%dp",conf,sbsmag);
  char const * outdir = !is_vary_cut ? Form("pdout/fits/%s/",confmag) : Form("pdout/fits/%s/sysstdy/",confmag);
  TString outFile = Form("%s%s_fit_dx_%s_pass%d_%s_sbs%d_sbs%dp_model%d.root",outdir,filebase.c_str(),key,pass,gen.c_str(),conf,sbsmag,model);
  TString outPlot = outFile; outPlot.ReplaceAll(".root",".pdf");
  TString outData = outFile; outData.ReplaceAll(".root",".csv"); ofstream outdata; outdata.open(outData);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // Creating important histograms
  vector<double> h_dx; jmgr->GetVectorFromSubKey<double>(key,"h_dx",h_dx);
  TH1F *h_dxHCAL_data = new TH1F("h_dxHCAL_data","",int(h_dx[0]),h_dx[1],h_dx[2]); h_dxHCAL_data->Sumw2();
  TH1F *h_dxHCAL_data_norm = new TH1F("h_dxHCAL_data_norm","",int(h_dx[0]),h_dx[1],h_dx[2]); // charge normalized & live time corrected 
  TH1F *h_dxHCAL_data_CT = new TH1F("h_dxHCAL_data_CT","",int(h_dx[0]),h_dx[1],h_dx[2]);
  TH1F *h_dxHCAL_simu_p = new TH1F("h_dxHCAL_simu_p","",int(h_dx[0]),h_dx[1],h_dx[2]); h_dxHCAL_simu_p->Sumw2();
  TH1F *h_dxHCAL_simu_n = new TH1F("h_dxHCAL_simu_n","",int(h_dx[0]),h_dx[1],h_dx[2]); h_dxHCAL_simu_n->Sumw2();
  TH1F *h_dxHCAL_simu_p_norm = new TH1F("h_dxHCAL_simu_p_norm","",int(h_dx[0]),h_dx[1],h_dx[2]); // Charge normalized
  TH1F *h_dxHCAL_simu_n_norm = new TH1F("h_dxHCAL_simu_n_norm","",int(h_dx[0]),h_dx[1],h_dx[2]); // Charge normalized
  TH1F *h_dxHCAL_bg_data = new TH1F("h_dxHCAL_bg_data","",int(h_dx[0]),h_dx[1],h_dx[2]); h_dxHCAL_bg_data->Sumw2();
  TH1F *h_dxHCAL_bg_inel_p = new TH1F("h_dxHCAL_bg_inel_p","",int(h_dx[0]),h_dx[1],h_dx[2]); h_dxHCAL_bg_inel_p->Sumw2();
  TH1F *h_dxHCAL_bg_inel_n = new TH1F("h_dxHCAL_bg_inel_n","",int(h_dx[0]),h_dx[1],h_dx[2]); h_dxHCAL_bg_inel_n->Sumw2();
  TH1F *h_dxHCAL_bg_inel = new TH1F("h_dxHCAL_bg_inel","",int(h_dx[0]),h_dx[1],h_dx[2]); h_dxHCAL_bg_inel->Sumw2();
  // kinematic histo "true"
  TH1F *h_vQ2 = new TH1F("h_vQ2","",300,0,16);
  TH1F *h_vetheta= new TH1F("h_vetheta","",300,0.35,1.05);
  // vs Run number
  // ## Explicit x bins for Rnum histos -- Needed to avoid round off error introduced by ROOT's default way of calculating bin edges
  int nruns = -1; 
  std::vector<CodaRun> cruns; util_pd::ReadRunList("../DB",nruns,conf,target,pass,sbsmag,0,cruns);
  int minRnum = cruns[0].runnum-1; int maxRnum = cruns[cruns.size()-1].runnum+2;
  std::vector<double> xbinsRnum; CalcRnumBinEdges(minRnum,maxRnum,0,xbinsRnum);
  int nbinRnum = xbinsRnum.size()-1;
  //
  TH2F *h2_dxHCAL_vs_rnum = new TH2F("h2_dxHCAL_vs_rnum","",nbinRnum,&(xbinsRnum)[0],int(h_dx[0]),h_dx[1],h_dx[2]);
  TH2F *h2_dxHCAL_vs_rnum_norm = new TH2F("h2_dxHCAL_vs_rnum","",nbinRnum,&(xbinsRnum)[0],int(h_dx[0]),h_dx[1],h_dx[2]); // charge normalized & live time corrected

  // Cut variation **************
  // Although not necessary, keeping this loop separate gives more control and clarity
  // forming the cuts
  std::string param_to_vary = jmgr->GetValueFromSubKey_str(key,"param_to_vary");
  vector<double> cut_range; jmgr->GetVectorFromSubKey<double>(key,"cut_iter_min_width",cut_range);
  vector<double> h_cut_param; jmgr->GetVectorFromSubKey<double>(key,"h_cut_param",h_cut_param);
  double min = cut_range[1], width = cut_range[2]; 
  int iter = is_vary_cut ? (int)cut_range[0] : 1;
  // draawing cut histo 
  // ** for data
  std::string cuts_sig_data_modified = use_custom_fiduCut ? cuts_for_signal_data + "&&fiduCut" : cuts_for_signal_data;
  TH1F *hcut = (TH1F*)data_rdf_filtered.Histo1D({"hcut","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]},param_to_vary)->Clone();
  customize_hcut(hcut); hcut->SetTitle(Form("%s {%s}",param_to_vary.c_str(),cuts_sig_data_modified.c_str()));
  TH1F *hcut_p = (TH1F*)data_rdf_filtered.Filter("pCut").Histo1D({"hcut_p","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]},param_to_vary)->Clone();
  customize_hcut_p(hcut_p); hcut_p->SetTitle(Form("%s {pCut&&%s}",param_to_vary.c_str(),cuts_sig_data_modified.c_str()));
  TH1F *hcut_n = new TH1F("hcut_n","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]);
  if (!is_elastic) {
    hcut_n = (TH1F*)data_rdf_filtered.Filter("nCut").Histo1D({"hcut_n","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]},param_to_vary)->Clone();
    customize_hcut_n(hcut_n); hcut_n->SetTitle(Form("%s {nCut&&%s}",param_to_vary.c_str(),cuts_sig_data_modified.c_str()));
  }
  // ** for simu
  std::string cuts_sig_simu_modified = use_custom_fiduCut ? cuts_for_signal_simu + "&&fiduCut" : cuts_for_signal_simu;
  TH1F *hcut_simu = (TH1F*)simu_rdf_filtered.Histo1D({"hcut_simu","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]},param_to_vary,"weight")->Clone();
  customize_hcut(hcut_simu); hcut_simu->SetTitle(Form("%s {%s}",param_to_vary.c_str(),cuts_sig_simu_modified.c_str()));
  TH1F *hcut_p_simu = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"hcut_p_simu","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]},param_to_vary,"weight")->Clone();
  customize_hcut_p(hcut_p_simu); hcut_p_simu->SetTitle(Form("%s {pCut&&%s}",param_to_vary.c_str(),cuts_sig_simu_modified.c_str()));
  TH1F *hcut_n_simu = new TH1F("hcut_n_simu","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]);
  if (!is_elastic) {
    hcut_n_simu = (TH1F*)simu_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"hcut_n_simu","",int(h_cut_param[0]),h_cut_param[1],h_cut_param[2]},param_to_vary,"weight")->Clone();
    customize_hcut_n(hcut_n_simu); hcut_n_simu->SetTitle(Form("%s {nCut&&%s}",param_to_vary.c_str(),cuts_sig_simu_modified.c_str()));
  }
  // choosing cut variation style
  int cut_vary_style = jmgr->GetValueFromSubKey<int>(key,"cut_vary_style");
  // 0 -> scan from low to high with fixed width
  //      NOTE: cut_range[0] = # slices, [1] = low, [2] = slice width, in this case
  // 1 -> scan around a fixed mean with increasing width.
  //      NOTE: cut_range[0] = # slices, [1] = mean, [2] = width, in this case
  // 2 -> scan from low to high with slices of equal statistics.
  //      NOTE: cut_range[0] = # slices, [1] = low, [2] = high, in this case
  // 3 -> increase threshold by a fixed amount
  //      NOTE: cut_range[0] = # slices, [1] = low, [2] = threshold increment, in this case  
  // 4 -> vary fidu cut.
  //      NOTE: Only cut_range[0] matters, sets the range but fidu_vary_* dictates variation
  bool apply_to_data_only = jmgr->GetValueFromSubKey<int>(key,"apply_to_data_only");
  double low = cut_vary_style==1 ? min-width : min;
  double high = cut_vary_style==3 ? h_cut_param[2] : min+width; 
  // fidu cut variation (Cut style 4 -- Very different than the others)
  std::vector<int> fvary_xp; jmgr->GetVectorFromSubKey<int>(key,"fidu_vary_xp",fvary_xp);
  std::vector<int> fvary_xn; jmgr->GetVectorFromSubKey<int>(key,"fidu_vary_xn",fvary_xn);
  std::vector<int> fvary_y; jmgr->GetVectorFromSubKey<int>(key,"fidu_vary_y",fvary_y);
  std::vector<std::vector<double>> hcal_SMs;
  // --
  std::vector<double> minval, maxval;
  std::vector<std::string> cuts, cuts_p, cuts_n, cuts_2;  
  if (is_vary_cut) {
    std::vector<double> xrangeEqStat;
    
    if (cut_vary_style==2) {
      int ndiv = cut_range[0]; //desired # of equi-stat slices
      low = cut_range[1]; high = cut_range[2];
      std::cout << "Varying cut w/ equi-stat slices..\n";
      util_pd::FindEqualStatBins(hcut,low,high,ndiv,1,xrangeEqStat);
      high = xrangeEqStat[1]; //initializing for the first slice
    }
    
    for (int i=0; i<iter; i++) {
      minval.push_back(low); maxval.push_back(high);
      std::string cut, cut_2;
      char low_buff[20]; std::snprintf(low_buff,20,"%.3f",low); std::string low_str(low_buff);
      char high_buff[20]; std::snprintf(high_buff,20,"%.3f",high); std::string high_str(high_buff);
      if (cut_vary_style==0 || cut_vary_style==1 || cut_vary_style==2) {
	cut = param_to_vary+">"+low_str+"&&"+param_to_vary+"<="+high_str; 
	cut_2 = low_str+"<"+param_to_vary+"<="+high_str;
	hcal_SMs.push_back(hcal_SM);
      } 
      else if (cut_vary_style==3) { 
	cut = param_to_vary+">"+low_str; cut_2 = cut; 
	hcal_SMs.push_back(hcal_SM);
      }
      else if (cut_vary_style==4) {
	cut = "1"; cut_2 = "ARCut&&SMCut";
	double SM_xp, SM_xn, SM_y;
	SM_xp = fvary_xp[0] ? SM_w[0]*(1.+0.01*fvary_xp[i+1]) : SM_w[0];
	SM_xn = fvary_xn[0] ? SM_w[1]*(1.+0.01*fvary_xn[i+1]) : SM_w[1];
	SM_y = fvary_y[0] ? SM_w[2]*(1.+0.01*fvary_y[i+1]) : SM_w[2];
	std::vector<double> hcal_SM_i = cut::hcal_safety_margin(SM_xp,SM_xn,SM_y,hcal_AR);
	hcal_SMs.push_back(hcal_SM_i);
      }
      cuts.push_back(cut); cuts_2.push_back(cut_2);
      // defining cuts for MC
      if (apply_to_data_only) cut = "1";
      std::string cut_p = cut + "&&mc_fnucl==1"; cuts_p.push_back(cut_p);
      std::string cut_n = cut + "&&mc_fnucl==0"; cuts_n.push_back(cut_n);
      // update high and low
      if (cut_vary_style==0) { low = high; high += width; }
      else if (cut_vary_style==1) { low -= width; high += width; }
      else if (cut_vary_style==2) { low = xrangeEqStat[i+1]; high = xrangeEqStat[i+2]; }
      else if (cut_vary_style==3) { low += width; high = h_cut_param[2]; }
    }
  } else {
    cuts.push_back("1"); cuts_p.push_back("mc_fnucl==1"); cuts_n.push_back("mc_fnucl==0"); 
    cuts_2.push_back("N/A"); minval.push_back(0); maxval.push_back(0);
    hcal_SMs.push_back(hcal_SM);
  }

  // Defining all the output files
  TString outGIF = outFile; outGIF.ReplaceAll(".root",".gif");
  TString outPNG = outFile; outPNG.ReplaceAll(".root","");

  // Fits
  vector<double> dx_fit_range; jmgr->GetVectorFromSubKey<double>(key,"dx_fit_range",dx_fit_range);
  vector<double> reject_points; jmgr->GetVectorFromSubKey<double>(key,"SB_reject_points",reject_points);

  // summary histo
  TH1F *hgist_cv_2 = new TH1F("hgist_cv_2","",iter,-0.5,iter-0.5);
  // Canvas to plot cut region
  TCanvas *cCut = new TCanvas("cCut","cCut",1400,800);
  if (cut_vary_style==4) cCut->Divide(2,2);
  else cCut->Divide(2,1);

  // ---- 
  for (int i=0; i<iter; i++) {

    // forming the fiducial cut
    auto fiduCut = [&](double x,double y,double xExp,double yExp) {
      return cut::inHCAL_activeA(x,y,hcal_AR) && cut::inHCAL_safety_margin(target,xExp,yExp,sbs_kick,hcal_SMs[i]);
    };
    
    // Filling physics histograms with appropriate cuts -----
    if (!use_custom_fiduCut) { // don't use custom fidu cut
      h_dxHCAL_data = (TH1F*)data_rdf_filtered.Filter(cuts[i]).Histo1D({"h_dxHCAL_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      if (temp) h_dxHCAL_data_norm = (TH1F*)data_rdf_filtered.Filter(cuts[i]).Histo1D({"h_dxHCAL_data_norm","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight_norm")->Clone();
      h_dxHCAL_data_CT = (TH1F*)data_rdf_filtered.Filter(coinT_cut.c_str()).Histo1D({"h_dxHCAL_data_CT","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      h_dxHCAL_simu_p = (TH1F*)simu_rdf_filtered.Filter(cuts_p[i]).Histo1D({"h_dxHCAL_simu_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight")->Clone();
      h_dxHCAL_simu_p_norm = (TH1F*)simu_rdf_filtered.Filter(cuts_p[i]).Histo1D({"h_dxHCAL_simu_p_norm","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight_norm")->Clone();
      if (!is_elastic) {
	h_dxHCAL_simu_n = (TH1F*)simu_rdf_filtered.Filter(cuts_n[i]).Histo1D({"h_dxHCAL_simu_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight")->Clone();
	h_dxHCAL_simu_n_norm = (TH1F*)simu_rdf_filtered.Filter(cuts_n[i]).Histo1D({"h_dxHCAL_simu_n_norm","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight_norm")->Clone();
      }
      // bg histos
      h_dxHCAL_bg_data = (TH1F*)bg_data_rdf_filtered.Histo1D({"h_dxHCAL_bg_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      h_dxHCAL_bg_inel_p = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_bg_inel_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight")->Clone();
      if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_bg_inel_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight")->Clone();
      h_dxHCAL_bg_inel = (TH1F*)h_dxHCAL_bg_inel_p->Clone();
      if (!is_elastic) h_dxHCAL_bg_inel->Add(h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n);
      // h_dxHCAL_bg_inel_p = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_bg_inel_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight")->Clone();
      // if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_bg_inel_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight")->Clone();
      // h_dxHCAL_bg_inel = (TH1F*)h_dxHCAL_bg_inel_p->Clone(); h_dxHCAL_bg_inel->Add(h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n);
      // kinematic histos "true"
      if (!apply_to_data_only||!is_elastic) {
	h_vQ2 =  (TH1F*)simu_rdf_filtered.Filter(cuts[i]).Histo1D({"h_vQ2","",300,0,16},"vQ2","weight")->Clone();
	h_vetheta = (TH1F*)simu_rdf_filtered.Filter(cuts[i]).Histo1D({"h_vetheta","",300,0.35,1.05},"vetheta","weight")->Clone();
      }
      // ** vs Run number histos **
      h2_dxHCAL_vs_rnum = (TH2F*)data_rdf_filtered.Filter(cuts[i]).Histo2D({"h2_dxHCAL_vs_rnum","",nbinRnum,&(xbinsRnum)[0],int(h_dx[0]),h_dx[1],h_dx[2]},"rnum","dx")->Clone();
      if (temp) h2_dxHCAL_vs_rnum_norm = (TH2F*)data_rdf_filtered
		  .Filter(cuts[i]).Histo2D({"h2_dxHCAL_vs_rnum_norm","",nbinRnum,&(xbinsRnum)[0],int(h_dx[0]),h_dx[1],h_dx[2]},"rnum","dx","weight_norm")->Clone();

    } else { // use custom fidu cut
      h_dxHCAL_data = (TH1F*)data_rdf_filtered
	.Filter(cuts[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	.Histo1D({"h_dxHCAL_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      if (temp) h_dxHCAL_data_norm = (TH1F*)data_rdf_filtered
		  .Filter(cuts[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
		  .Histo1D({"h_dxHCAL_data_norm","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight_norm")->Clone();
      h_dxHCAL_data_CT = (TH1F*)data_rdf_filtered.Filter(coinT_cut.c_str())
	.Histo1D({"h_dxHCAL_data_CT","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      h_dxHCAL_simu_p = (TH1F*)simu_rdf_filtered
	.Filter(cuts_p[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	.Histo1D({"h_dxHCAL_simu_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight")->Clone();
      h_dxHCAL_simu_p_norm = (TH1F*)simu_rdf_filtered
	.Filter(cuts_p[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	.Histo1D({"h_dxHCAL_simu_p_norm","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight_norm")->Clone();
      if (!is_elastic) {
	h_dxHCAL_simu_n = (TH1F*)simu_rdf_filtered
	  .Filter(cuts_n[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	  .Histo1D({"h_dxHCAL_simu_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight")->Clone();
	h_dxHCAL_simu_n_norm = (TH1F*)simu_rdf_filtered
	  .Filter(cuts_n[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	  .Histo1D({"h_dxHCAL_simu_n_norm","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight_norm")->Clone();
      }
      // kinematic histos "true"
      if (!apply_to_data_only||!is_elastic) {
	h_vQ2 =  (TH1F*)simu_rdf_filtered.Filter(cuts[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"}).Histo1D({"h_vQ2","",300,0,16},"vQ2","weight")->Clone();
	h_vetheta = (TH1F*)simu_rdf_filtered.Filter(cuts[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"}).Histo1D({"h_vetheta","",300,0.35,1.05},"vetheta","weight")->Clone();
      }
      // ** vs Run number histos **
      h2_dxHCAL_vs_rnum = (TH2F*)data_rdf_filtered
	.Filter(cuts[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	.Histo2D({"h2_dxHCAL_vs_rnum","",nbinRnum,&(xbinsRnum)[0],int(h_dx[0]),h_dx[1],h_dx[2]},"rnum","dx")->Clone();
      if (temp) h2_dxHCAL_vs_rnum_norm = (TH2F*)data_rdf_filtered
		  .Filter(cuts[i]).Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
		  .Histo2D({"h2_dxHCAL_vs_rnum_norm","",nbinRnum,&(xbinsRnum)[0],int(h_dx[0]),h_dx[1],h_dx[2]},"rnum","dx","weight_norm")->Clone();
      
      // bg histos ----- 
      // (data)
      if (apply_to_bg_data) // applying custom fiduCut to data bg
	h_dxHCAL_bg_data = (TH1F*)bg_data_rdf_filtered
	  .Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	  .Histo1D({"h_dxHCAL_bg_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      else
	h_dxHCAL_bg_data = (TH1F*)bg_data_rdf_filtered.Histo1D({"h_dxHCAL_bg_data","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx")->Clone();
      // (simu)
      if (!apply_to_bg_simu) { // applying custom fiduCut to simu bg
	h_dxHCAL_bg_inel_p = (TH1F*)inel_rdf_filtered
	  .Filter("mc_fnucl==1").Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	  .Histo1D({"h_dxHCAL_bg_inel_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight")->Clone();
	if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)inel_rdf_filtered
			   .Filter("mc_fnucl==0").Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
			   .Histo1D({"h_dxHCAL_bg_inel_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight")->Clone();
	// h_dxHCAL_bg_inel_p = (TH1F*)inel_rdf_filtered
	//   .Filter("mc_fnucl==1").Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	//   .Histo1D({"h_dxHCAL_bg_inel_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight")->Clone();
	// if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)inel_rdf_filtered
	// 		   .Filter("mc_fnucl==0").Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	// 		   .Histo1D({"h_dxHCAL_bg_inel_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight")->Clone();
      } 
      else {
	h_dxHCAL_bg_inel_p = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_bg_inel_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_p","weight")->Clone();
	if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_bg_inel_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx_shifted_n","weight")->Clone();
	// h_dxHCAL_bg_inel_p = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==1").Histo1D({"h_dxHCAL_bg_inel_p","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight")->Clone();
	// if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)inel_rdf_filtered.Filter("mc_fnucl==0").Histo1D({"h_dxHCAL_bg_inel_n","",int(h_dx[0]),h_dx[1],h_dx[2]},"dx","weight")->Clone();
      }
    }    
    h_dxHCAL_bg_inel = (TH1F*)h_dxHCAL_bg_inel_p->Clone();
    if (!is_elastic) h_dxHCAL_bg_inel->Add(h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n);
    // ------------
    
    // Now that we have all the important histograms formed, let's write them to the output
    // tree for further analysis.
    if (!is_vary_cut) {
      h_dxHCAL_data->Write();
      if (temp) h_dxHCAL_data_norm->Write(); 
      h_dxHCAL_data_CT->Write();
      h_dxHCAL_simu_p->Write();
      h_dxHCAL_simu_p_norm->Write();
      h_dxHCAL_bg_data->Write();
      h_dxHCAL_bg_inel_p->Write();
      if (!is_elastic) {
	h_dxHCAL_simu_n->Write();
	h_dxHCAL_simu_n_norm->Write();
	h_dxHCAL_bg_inel->Write();
	h_dxHCAL_bg_inel_n->Write();
      }
      // kinematic histos "true"
      if (!apply_to_data_only||!is_elastic) {
	h_vQ2->Write();
	h_vetheta->Write();
      }
      // vs Runnum histos
      h2_dxHCAL_vs_rnum->Write();
      if (temp) h2_dxHCAL_vs_rnum_norm->Write();
    }
    // ------------

    if (is_vary_cut) {
      // drawing cut histos
      cCut->cd();
      gStyle->SetOptStat("e");
      if (cut_vary_style!=4) { // Anything other than fidu cut
	// data
	cCut->cd(1);
	gPad->SetGridx();
	hcut->Draw();
	hcut_p->Draw("same");
	if (!is_elastic) hcut_n->Draw("same");
	util_pd::PlotCutRegion(minval[i],maxval[i]);
	TLegend *lCut = new TLegend(0.6,0.75,0.9,0.9);
	lCut->SetTextFont(42); //lCut->SetFillStyle(0);
	lCut->AddEntry(hcut,Form("%s (Data)",param_to_vary.c_str()),"l");
	lCut->AddEntry(hcut_p,Form("%s w/ pCut",param_to_vary.c_str()),"l");
	if (!is_elastic) lCut->AddEntry(hcut_n,Form("%s w/ nCut",param_to_vary.c_str()),"l");
	AddCutToLegend(lCut,cuts_2[i].c_str());
	lCut->Draw();
	// simu
	cCut->cd(2);
	gPad->SetGridx();
	hcut_simu->Draw();
	hcut_p_simu->Draw("same");
	if (!is_elastic) hcut_n_simu->Draw("same");
	util_pd::PlotCutRegion(minval[i],maxval[i]);
	TLegend *lCut_simu = new TLegend(0.6,0.75,0.9,0.9);
	lCut_simu->SetTextFont(42); //lCut->SetFillStyle(0);
	lCut_simu->AddEntry(hcut,Form("%s (MC)",param_to_vary.c_str()),"l");
	lCut_simu->AddEntry(hcut_p,Form("%s w/ pCut",param_to_vary.c_str()),"l");
	if (!is_elastic) lCut_simu->AddEntry(hcut_n,Form("%s w/ nCut",param_to_vary.c_str()),"l");
	AddCutToLegend(lCut_simu,cuts_2[i].c_str());
	lCut_simu->Draw();
      } else {
	TH2F *h2_fiduCut_n = (TH2F*)data_rdf_filtered
	  .Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	  .Histo2D({"h2_fiduCut_n","",200,-1,1,200,-3,1.5},"yHCAL_exp","xHCAL_exp")->Clone(); 
	TH2F *h2_fiduCut_p = (TH2F*)data_rdf_filtered
	  .Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	  .Histo2D({"h2_fiduCut_p","",200,-1,1,200,-3,1.5},"yHCAL_exp","xExp_shifted")->Clone();
	TH2F *h2_fiduCut_n_MC = (TH2F*)simu_rdf_filtered
	  .Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	  .Histo2D({"h2_fiduCut_n","",200,-1,1,200,-3,1.5},"yHCAL_exp","xHCAL_exp")->Clone(); 
	TH2F *h2_fiduCut_p_MC = (TH2F*)simu_rdf_filtered
	  .Filter(fiduCut,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
	  .Histo2D({"h2_fiduCut_p","",200,-1,1,200,-3,1.5},"yHCAL_exp","xExp_shifted")->Clone();
	//-- histo customization
	customize_h2fiduCut(h2_fiduCut_n,"n","Data",0); customize_h2fiduCut(h2_fiduCut_p,"p","Data",sbs_kick);
	customize_h2fiduCut(h2_fiduCut_n_MC,"n","MC",0); customize_h2fiduCut(h2_fiduCut_p_MC,"p","MC",sbs_kick);
	//---
	cCut->cd(1); //gStyle->SetOptStat(0); 
	h2_fiduCut_n->Draw("colz");
	PlotFiduCut(pass,hcal_AR,hcal_SMs[i]);
	cCut->cd(2); //gStyle->SetOptStat(0); 
	h2_fiduCut_p->Draw("colz");
	PlotFiduCut(pass,hcal_AR,hcal_SMs[i]);
	cCut->cd(3); //gStyle->SetOptStat(0); 
	h2_fiduCut_n_MC->Draw("colz");
	PlotFiduCut(pass,hcal_AR,hcal_SMs[i]);
	cCut->cd(4); //gStyle->SetOptStat(0); 
	h2_fiduCut_p_MC->Draw("colz");
	PlotFiduCut(pass,hcal_AR,hcal_SMs[i]);
      }
      cCut->SetTickx(); cCut->SetTicky(); 
      cCut->Update(); cCut->Write(); 
      if (i==0) {
	cCut->SaveAs(Form("%s+150",outGIF.Data()));
	cCut->SaveAs(Form("%s[",outPlot.Data())); 
      }
      cCut->SaveAs(Form("%s",outPlot.Data())); 
      //cCut->SaveAs(Form("%s_cCut_%d.png",outPNG.Data(),i)); 
      cCut->SaveAs(Form("%s+150",outGIF.Data()));
    }
    // ***

    /*######################################
      ## Fitting elastic dx distributions ##
      ###################################### */
    std::vector<double> R_vals,Rerr_vals,chi2,NDF,B_vals,Berr_vals;
    std::vector<double> pCnt,pCnt_err,nCnt,nCnt_err,bgCnt,bgCnt_err;
    if (is_elastic) {
      // Canvas 0 : Fitting data/MC w/o any background
      TCanvas *c0 = util_pd::TC("c0",1,1); gStyleFitCanvas();
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
      // TCanvas *c1 = util_pd::TC("c1",1,1); gStyleFitCanvas();
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
      TCanvas *c2 = util_pd::TC("c2",1,1); gStyleFitCanvas();
      c2->cd(); gStyle->SetOptFit(1);
      vector<TH1F*> ho2;
      TF1 *f2;
      if (is_vary_pnXOff) {
	f2 = fit::fit_1hs_1pbg_THI_xOffVary(dx_fit_range,
					    h_dxHCAL_data,h_dxHCAL_simu_p,Opoly,pnXOff_range,
					    ho2);
      } else {
	f2 = fit::fit_1hs_1pbg_THI(dx_fit_range,
				   h_dxHCAL_data,h_dxHCAL_simu_p,Opoly,
				   ho2);
      }
      // converting fit fn to a hostogram
      TH1F *hf2 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f2,hf2);
      ho2[0]->Draw(); c2->Update(); 
      // grabbing statbox of the fitted histo
      TPaveStats *st2 = (TPaveStats*)ho2[0]->FindObject("stats");
      ho2[0]->SetBit(TH1::kNoStats); // Sets up the stat box for later modification
      // getting pads for pull plot
      std::vector<TPad*> p2 = util_pd::GetPadsForPullPlot(c2);
      //
      // preparing the pad for data/MC fit
      //
      p2[0]->cd();
      // plotting the histos
      h_dxHCAL_data->Draw("E"); util_pd::customize_data(h_dxHCAL_data);
      hf2->Draw("same HIST"); util_pd::customize_gfit(hf2,1);
      ho2[1]->Draw("same HIST"); util_pd::customize_nsig(ho2[1],1);
      ho2[2]->Draw("same HIST"); util_pd::customize_hbg(ho2[2],1);
      // redrawing the stat box
      st2->SetX1NDC(0.62); st2->SetX2NDC(0.9); st2->SetY2NDC(0.9);
      st2->Draw("same");
      // drawing legend
      TLegend *l2=new TLegend(0.10,0.74,0.33,0.9);
      l2->SetTextFont(42);
      l2->AddEntry(h_dxHCAL_data,"Data","p");
      l2->AddEntry(hf2,"Fit (MC + bg.)","lf");
      l2->AddEntry(ho2[1],"Signal (from MC)","lf");
      l2->AddEntry(ho2[2],Form("Bg. (Poly. %d)",Opoly),"lf");
      l2->Draw();
      // preparing the pad for residual
      p2[1]->cd();
      ho2[3]->Draw(); util_pd::customize_residual(ho2[3]);
      // drawing a horizontal line at y = 0
      util_pd::DrawZeroLine(p2[1],dx_fit_range[0],dx_fit_range[1]);
      // --------

      // Canvas 3 : Sideband fit (distribution from data)
      TCanvas *c3 = util_pd::TC("c3",1,1); gStyleFitCanvas();
      c3->cd(); gStyle->SetOptFit(1);
      //vector<double> reject_points{-1.3,-0.15};
      vector<TH1F*> ho3;
      TF1* bg3 = fit::fit_1pbg_SB(dx_fit_range,
				  reject_points,
				  3,//Opoly,
				  fit::GetFitParams(f2),
				  h_dxHCAL_data,
				  ho3);
      ho3[0]->Draw(); c3->Update(); 
      // grabbing statbox of the fitted histo
      TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
      ho3[0]->SetBit(TH1::kNoStats); // Sets up the stat box for later modification
      // // getting pads for pull plot
      // std::vector<TPad*> p3 = util_pd::GetPadsForPullPlot(c3);
      // //
      // // preparing the pad for data/MC fit
      // //
      // p3[0]->cd();
      // plotting the histos
      ho3[0]->Draw("E"); util_pd::customize_data(ho3[0]);
      ho3[1]->Draw("same HIST"); util_pd::customize_nsig(ho3[1],1);
      ho3[2]->Draw("same HIST"); util_pd::customize_hbg(ho3[2],1);
      // redrawing the stat box
      st3->SetX1NDC(0.62); st3->SetX2NDC(0.9); st3->SetY2NDC(0.9);
      st3->Draw("same");
      // // preparing the pad for residual
      // p3[1]->cd();
      // ho3[3]->Draw(); util_pd::customize_residual(ho3[3]);
      // // drawing a horizontal line at y = 0
      // util_pd::DrawZeroLine(p3[1],dx_fit_range[0],dx_fit_range[1]);
      // drawing legend
      TLegend *l3=new TLegend(0.10,0.74,0.33,0.9);
      l3->SetTextFont(42);
      l3->AddEntry(h_dxHCAL_data,"Data","p");
      l3->AddEntry(ho3[2],Form("Bg. (Side Band)"),"lf");
      l3->AddEntry(ho3[1],"Signal (Data-Bg.)","lf");
      //l3->AddEntry(ho3[3],"Residual","p");
      l3->Draw();

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
      TCanvas *c4 = util_pd::TC("c4",1,1); gStyleFitCanvas();
      c4->cd(); gStyle->SetOptFit(1);
      vector<TH1F*> ho4;
      TF1 *f4 = fit::fit_1hs_1hbg_THI(dx_fit_range,
				      h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_bg_inel,
				      ho4);
      ho4[0]->Draw(); customize_ht(ho4[0]); customize_dx(ho4[0]);
      ho4[1]->Draw("same"); customize_hs(ho4[1]); customize_dx(ho4[1]);
      ho4[2]->Draw("same"); util_pd::customize_hbg(ho4[2],1); customize_dx(ho4[2]);
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

      // further customization of the data histo
      FurtherCustoizeDataHisto(h_dxHCAL_data,cuts_for_signal_data);
      
      // writing out the canvases
      c0->Update(); c0->Write(); if (!is_vary_cut&&i==0) c0->SaveAs(Form("%s[",outPlot.Data())); 
      c0->SaveAs(Form("%s",outPlot.Data())); 
      //c1->Update(); c1->Write(); c1->SaveAs(Form("%s",outPlot.Data())); 
      c2->Update(); c2->Write(); c2->SaveAs(Form("%s",outPlot.Data())); 
      c3->Update(); c3->Write(); c3->SaveAs(Form("%s",outPlot.Data())); 
      c4->Update(); c4->Write(); c4->SaveAs(Form("%s",outPlot.Data())); 
      //c5->Update(); c5->Write(); c5->SaveAs(Form("%s",outPlot.Data())); 
      if (i==iter-1) c4->SaveAs(Form("%s]",outPlot.Data()));      
    } // elastic

    /*#################################
      ## Fitting QE dx distributions ##
      ################################# */
    else {
      // various outputs
      if (i==0) {
	outdata << "cut,min,max,RMCnf,RMCnferr,RpMC,"
		<< "chi20,NDF0,R0,R0err,"
		<< "chi21,NDF1,R1,R1err,B1,B1err,Yp1,Yp1err,Yn1,Yn1err,Ybg1,Ybg1err,"
		<< "chi22,NDF2,R2,R2err,B2,B2err,Yp2,Yp2err,Yn2,Yn2err,Ybg2,Ybg2err,"
		<< "chi23,NDF3,R3,R3err,B3,B3err,Yp3,Yp3err,Yn3,Yn3err,Ybg3,Ybg3err,"
		<< "chi24,NDF4,R4,R4err,B4,B4err,Yp4,Yp4err,Yn4,Yn4err,Ybg4,Ybg4err,"
		<< "chi25,NDF5,R5,R5err,B5,B5err,Yp5,Yp5err,Yn5,Yn5err,Ybg5,Ybg5err,"
		<< "\n";
      }
      
      // calculating MC ratio (before fit) for comparisons
      double count_p_MC, error_p_MC;
      double count_n_MC, error_n_MC;
      count_p_MC = h_dxHCAL_simu_p->IntegralAndError(1,h_dxHCAL_simu_p->GetNbinsX(),error_p_MC);
      count_n_MC = h_dxHCAL_simu_n->IntegralAndError(1,h_dxHCAL_simu_n->GetNbinsX(),error_n_MC);
      double R_MC_nofit = count_n_MC / count_p_MC;
      double R_MC_nofit_err = R_MC_nofit*sqrt(pow(error_p_MC/count_p_MC,2) + pow(error_n_MC/count_n_MC,2));
      std::cout << "\nR_MC_nofit = " << R_MC_nofit << " | R_MC_nofit_err = " << R_MC_nofit_err << "\n\n";
    
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
      pCnt.push_back(0); pCnt_err.push_back(0);
      nCnt.push_back(0); nCnt_err.push_back(0);
      bgCnt.push_back(0); bgCnt_err.push_back(0);
      //ho[0]->Draw("E"); util_pd::customize_data(ho[0]); customize_dx(ho[0]);
      ho[1]->Draw(); customize_hs(ho[1]); ho[1]->SetStats(0);
      h_dxHCAL_data->Draw("E same"); util_pd::customize_data(h_dxHCAL_data);
      //f0->Draw("same");
      TLegend *l0 = new TLegend(0.10,0.73,0.35,0.9);
      l0->SetTextFont(42);
      l0->AddEntry(h_dxHCAL_data,"Data","p");
      //l0->AddEntry(f0,"Fit","l");
      l0->AddEntry(ho[1],"MC (Signal)","p");
      l0->AddEntry((TObject*)0,Form("MC p peak offset: %.3fm",dx_offset_p),"");
      l0->AddEntry((TObject*)0,Form("MC n peak offset: %.3fm",dx_offset_n),"");
      if (is_vary_cut) AddCutToLegend(l0,cuts_2[i].c_str());
      l0->Draw();
      // --- 

      // Canvas 1 : Fitting data/MC w/ background from data
      TCanvas *c1 = util_pd::TC("c1",1,1); gStyleFitCanvas();
      c1->cd();
      // performing the fit
      vector<TH1F*> ho1;
      TF1 *f1;
      if (is_vary_pnXOff) {
	f1 = fit::fit_2hs_1hbg_THI_xOffVary(dx_fit_range,
					    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_data,pnXOff_range,
					    ho1);
      } else {
	f1 = fit::fit_2hs_1hbg_THI(dx_fit_range,
				   h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_data,
				   ho1);
      }
      // grabbing fit params for future use
      chi2.push_back(f1->GetChisquare()); NDF.push_back(f1->GetNDF());
      R_vals.push_back(f1->GetParameter(1)); Rerr_vals.push_back(f1->GetParError(1));
      B_vals.push_back(f1->GetParameter(2)); Berr_vals.push_back(f1->GetParError(2));
      // converting fit fn to a hostogram
      TH1F *hf1 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f1,hf1);
      ho1[0]->Draw(); c1->Update(); 
      // grabbing statbox of the fitted histo
      TPaveStats *st1 = (TPaveStats*)ho1[0]->FindObject("stats");
      ho1[0]->SetBit(TH1::kNoStats); // Sets up the stat box for later modification
      // getting pads for pull plot
      std::vector<TPad*> p1 = util_pd::GetPadsForPullPlot(c1);
      //
      // preparing the pad for data/MC fit
      //
      p1[0]->cd();
      // drawing all the histograms
      //ho1[0]->Draw(); customize_ht(ho1[0]); customize_dx(ho1[0]);
      h_dxHCAL_data->Draw("E"); util_pd::customize_data(h_dxHCAL_data);
      hf1->Draw("same HIST"); util_pd::customize_gfit(hf1,1);
      ho1[4]->Draw("same HIST"); util_pd::customize_psig(ho1[4],1);
      ho1[5]->Draw("same HIST"); util_pd::customize_nsig(ho1[5],1);
      ho1[2]->Draw("same HIST"); util_pd::customize_hbg(ho1[2],1);
      // calculating yields
      std::vector<double> yo1; util_pd::GetYields(f1,ho1[4],ho1[5],ho1[2],yo1);
      pCnt.push_back(yo1[0]); pCnt_err.push_back(yo1[1]);
      nCnt.push_back(yo1[2]); nCnt_err.push_back(yo1[3]);
      bgCnt.push_back(yo1[4]); bgCnt_err.push_back(yo1[5]);
      // redrawing the stat box
      st1->SetX1NDC(0.6); st1->SetX2NDC(0.9); st1->SetY2NDC(0.9);
      // Modifying it to add yield ratio
      double yRatio1 = f1->GetParameter(1)*R_MC_nofit;
      double yRatio1_err = f1->GetParError(1)*R_MC_nofit;
      TText *t1 = st1->AddText(Form("R_{Yield} = %.4f #pm %.4f",yRatio1,yRatio1_err));
      customize_text(t1);
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
      ho1[3]->Draw(); util_pd::customize_residual(ho1[3]);
      // drawing a horizontal line at y = 0
      util_pd::DrawZeroLine(p1[1],dx_fit_range[0],dx_fit_range[1]);
      // ** ----- ***

      // Canvas 2 : Fitting data/MC w/ background from MC
      TCanvas *c2 = util_pd::TC("c2",1,1); gStyleFitCanvas();
      c2->cd();
      // performing the fit
      vector<TH1F*> ho2;
      TF1 *f2;
      if (is_vary_pnXOff) {
	f2 = fit::fit_2hs_2hbg_THI_xOffVary(dx_fit_range,
					    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n,pnXOff_range,
					    ho2);
      } else {
	f2 = fit::fit_2hs_1hbg_THI(dx_fit_range,
				   h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_inel,
				   ho2);
	// f2 = fit::fit_2hs_2hbg_THI(dx_fit_range,
	// 			   h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n,
	// 			   ho2);
      }
      if (is_vary_cut) {
	hgist_cv_2->SetBinContent(i+1,f2->GetParameter(1));
	hgist_cv_2->SetBinError(i+1,f2->GetParError(1));
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
      ho2[0]->SetBit(TH1::kNoStats); // Sets up the stat box for later modification
      // getting pads for pull plot
      std::vector<TPad*> p2 = util_pd::GetPadsForPullPlot(c2);
      //
      // preparing the pad for data/MC fit
      //
      p2[0]->cd();
      // drawing all the histograms
      //ho2[0]->Draw(); customize_ht(ho2[0]); customize_dx(ho2[0]);
      h_dxHCAL_data->Draw("E"); util_pd::customize_data(h_dxHCAL_data);
      hf2->Draw("same HIST"); util_pd::customize_gfit(hf2,1);
      ho2[4]->Draw("same HIST"); util_pd::customize_psig(ho2[4],1);
      ho2[5]->Draw("same HIST"); util_pd::customize_nsig(ho2[5],1);
      ho2[2]->Draw("same HIST"); util_pd::customize_hbg(ho2[2],1);
      // calculating yields
      std::vector<double> yo2; util_pd::GetYields(f2,ho2[4],ho2[5],ho2[2],yo2);
      pCnt.push_back(yo2[0]); pCnt_err.push_back(yo2[1]);
      nCnt.push_back(yo2[2]); nCnt_err.push_back(yo2[3]);
      bgCnt.push_back(yo2[4]); bgCnt_err.push_back(yo2[5]);
      // redrawing the stat box
      st2->SetX1NDC(0.6); st2->SetX2NDC(0.9); st2->SetY2NDC(0.9);
      // Modifying it to add yield ratio
      double yRatio2 = f2->GetParameter(1)*R_MC_nofit;
      double yRatio2_err = f2->GetParError(1)*R_MC_nofit;
      TText *t2 = st2->AddText(Form("R_{Yield} = %.4f #pm %.4f",yRatio2,yRatio2_err));
      customize_text(t2);
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
      ho2[3]->Draw(); util_pd::customize_residual(ho2[3]);
      // drawing a horizontal line at y = 0
      util_pd::DrawZeroLine(p2[1],dx_fit_range[0],dx_fit_range[1]);
      // ** ----- ***
      //c2->SaveAs(Form("%s_c2_%d.png",outPNG.Data(),i)); 
    
      // Canvas 3 : Fitting data/MC w/ polynomial background (2nd order - fixed)
      TCanvas *c3 = util_pd::TC("c3",1,1); gStyleFitCanvas();
      c3->cd(); 
      // performing the fit
      vector<TH1F*> ho3;
      TF1 *f3;
      if (is_vary_pnXOff) {
	f3 = fit::fit_2hs_1pbg_THI_xOffVary(dx_fit_range,
					    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,2,pnXOff_range,
					    ho3);
      } else {
	f3 = fit::fit_2hs_1pbg_THI(dx_fit_range,
				   h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,2,
				   ho3);
      }
      // grabbing fit params for future use
      chi2.push_back(f3->GetChisquare()); NDF.push_back(f3->GetNDF());
      R_vals.push_back(f3->GetParameter(1)); Rerr_vals.push_back(f3->GetParError(1));
      B_vals.push_back(f3->GetParameter(2)); Berr_vals.push_back(f3->GetParError(2));
      // converting fit fn to a hostogram
      TH1F *hf3 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f3,hf3);
      ho3[0]->Draw(); c3->Update(); 
      // grabbing statbox of the fitted histo
      TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
      ho3[0]->SetBit(TH1::kNoStats); // Sets up the stat box for later modification
      // getting pads for pull plot
      std::vector<TPad*> p3 = util_pd::GetPadsForPullPlot(c3);
      //
      // preparing the pad for data/MC fit
      //
      p3[0]->cd();
      // drawing all the histograms
      //ho3[0]->Draw(); customize_ht(ho3[0]); customize_dx(ho3[0]);
      h_dxHCAL_data->Draw("E"); util_pd::customize_data(h_dxHCAL_data);
      hf3->Draw("same HIST"); util_pd::customize_gfit(hf3,1);
      ho3[4]->Draw("same HIST"); util_pd::customize_psig(ho3[4],1);
      ho3[5]->Draw("same HIST"); util_pd::customize_nsig(ho3[5],1);
      ho3[2]->Draw("same HIST"); util_pd::customize_hbg(ho3[2],1);
      // calculating yields
      std::vector<double> yo3; util_pd::GetYields(f3,ho3[4],ho3[5],ho3[2],yo3);
      pCnt.push_back(yo3[0]); pCnt_err.push_back(yo3[1]);
      nCnt.push_back(yo3[2]); nCnt_err.push_back(yo3[3]);
      bgCnt.push_back(yo3[4]); bgCnt_err.push_back(yo3[5]);
      // redrawing the stat box
      st3->SetX1NDC(0.62); st3->SetX2NDC(0.9); st3->SetY2NDC(0.9);
      // Modifying it to add yield ratio
      double yRatio3 = f3->GetParameter(1)*R_MC_nofit;
      double yRatio3_err = f3->GetParError(1)*R_MC_nofit;
      TText *t3 = st3->AddText(Form("R_{Yield} = %.4f #pm %.4f",yRatio3,yRatio3_err));
      customize_text(t3);
      st3->Draw("same");
      // drawing a legend
      TLegend *l3=new TLegend(0.10,0.6,0.36,0.9);
      l3->SetTextFont(42);
      l3->AddEntry(h_dxHCAL_data,"Data","p");
      //l3->AddEntry(f3,"Fit (MC + poly. bg.)","l");
      l3->AddEntry(hf3,"Fit (QE MC + bg.)","lf");
      l3->AddEntry(ho3[4],"p signal (from MC)","lf");
      l3->AddEntry(ho3[5],"n signal (from MC)","lf");
      l3->AddEntry(ho3[2],Form("Bg. (poly. of order %d)",2),"lf");
      l3->AddEntry(ho3[3],"Residual","p");
      if (is_vary_cut) AddCutToLegend(l3,cuts_2[i].c_str());
      //l3->SetFillStyle(0); // makes legend box transparent
      l3->Draw();
      //
      // preparing the pad for residual
      //  
      p3[1]->cd();
      ho3[3]->Draw(); util_pd::customize_residual(ho3[3]);
      // drawing a horizontal line at y = 0
      util_pd::DrawZeroLine(p3[1],dx_fit_range[0],dx_fit_range[1]);
      // *******************
      //c3->SaveAs(Form("%s_c3_%d.png",outPNG.Data(),i)); 

      // Canvas 4 : Fitting data/MC w/ polynomial background
      TCanvas *c4 = util_pd::TC("c4",1,1); gStyleFitCanvas();
      c4->cd(); 
      // performing the fit
      vector<TH1F*> ho4;
      TF1 *f4;
      std::vector<double> gfit_params; jmgr->GetVectorFromSubKey<double>(key,"Par_guess_of_gaus_bg",gfit_params);
      if (is_vary_pnXOff) {
	f4 = fit::fit_2hs_1gbg_THI_xOffVary(dx_fit_range,
					    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,pnXOff_range,gfit_params,
					    ho4);
      } else {
	f4 = fit::fit_2hs_1gbg_THI(dx_fit_range,
				   h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,gfit_params,
				   ho4);
      }
      // grabbing fit params for future use
      chi2.push_back(f4->GetChisquare()); NDF.push_back(f4->GetNDF());
      R_vals.push_back(f4->GetParameter(1)); Rerr_vals.push_back(f4->GetParError(1));
      B_vals.push_back(f4->GetParameter(2)); Berr_vals.push_back(f4->GetParError(2));
      // converting fit fn to a hostogram
      TH1F *hf4 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f4,hf4);
      ho4[0]->Draw(); c4->Update(); 
      // grabbing statbox of the fitted histo
      TPaveStats *st4 = (TPaveStats*)ho4[0]->FindObject("stats");
      ho4[0]->SetBit(TH1::kNoStats); // Sets up the stat box for later modification
      // getting pads for pull plot
      std::vector<TPad*> p4 = util_pd::GetPadsForPullPlot(c4);
      //
      // preparing the pad for data/MC fit
      //
      p4[0]->cd();
      // drawing all the histograms
      //ho4[0]->Draw(); customize_ht(ho4[0]); customize_dx(ho4[0]);
      h_dxHCAL_data->Draw("E"); util_pd::customize_data(h_dxHCAL_data);
      hf4->Draw("same HIST"); util_pd::customize_gfit(hf4,1);
      ho4[4]->Draw("same HIST"); util_pd::customize_psig(ho4[4],1);
      ho4[5]->Draw("same HIST"); util_pd::customize_nsig(ho4[5],1);
      ho4[2]->Draw("same HIST"); util_pd::customize_hbg(ho4[2],1);
      // calculating yields
      std::vector<double> yo4; util_pd::GetYields(f4,ho4[4],ho4[5],ho4[2],yo4);
      pCnt.push_back(yo4[0]); pCnt_err.push_back(yo4[1]);
      nCnt.push_back(yo4[2]); nCnt_err.push_back(yo4[4]);
      bgCnt.push_back(yo4[4]); bgCnt_err.push_back(yo4[5]);
      // redrawing the stat box
      st4->SetX1NDC(0.62); st4->SetX2NDC(0.9); st4->SetY2NDC(0.9);
      // Modifying it to add yield ratio
      double yRatio4 = f4->GetParameter(1)*R_MC_nofit;
      double yRatio4_err = f4->GetParError(1)*R_MC_nofit;
      TText *t4 = st4->AddText(Form("R_{Yield} = %.4f #pm %.4f",yRatio4,yRatio4_err));
      customize_text(t4);
      st4->Draw("same");
      // drawing a legend
      TLegend *l4=new TLegend(0.10,0.6,0.36,0.9);
      l4->SetTextFont(42);
      l4->AddEntry(h_dxHCAL_data,"Data","p");
      //l4->AddEntry(f4,"Fit (MC + poly. bg.)","l");
      l4->AddEntry(hf4,"Fit (QE MC + bg.)","lf");
      l4->AddEntry(ho4[4],"p signal (from MC)","lf");
      l4->AddEntry(ho4[5],"n signal (from MC)","lf");
      l4->AddEntry(ho4[2],"Bg. (Gaussian)","lf");
      l4->AddEntry(ho4[3],"Residual","p");
      if (is_vary_cut) AddCutToLegend(l4,cuts_2[i].c_str());
      //l4->SetFillStyle(0); // makes legend box transparent
      l4->Draw();
      //
      // preparing the pad for residual
      //  
      p4[1]->cd();
      ho4[3]->Draw(); util_pd::customize_residual(ho4[3]);
      // drawing a horizontal line at y = 0
      util_pd::DrawZeroLine(p4[1],dx_fit_range[0],dx_fit_range[1]);
      // *******************
      //c4->SaveAs(Form("%s_c4_%d.png",outPNG.Data(),i)); 

      // Canvas 5 : Fitting data/MC w/ polynomial background
      TCanvas *c5 = util_pd::TC("c5",1,1); gStyleFitCanvas();
      c5->cd(); 
      // performing the fit
      vector<TH1F*> ho5;
      TF1 *f5;
      if (is_vary_pnXOff) {
	f5 = fit::fit_2hs_1pbg_THI_xOffVary(dx_fit_range,
					    h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,Opoly,pnXOff_range,
					    ho5);
      } else {
	f5 = fit::fit_2hs_1pbg_THI(dx_fit_range,
				   h_dxHCAL_data,h_dxHCAL_simu_p,h_dxHCAL_simu_n,Opoly,
				   ho5);
      }
      // grabbing fit params for future use
      chi2.push_back(f5->GetChisquare()); NDF.push_back(f5->GetNDF());
      R_vals.push_back(f5->GetParameter(1)); Rerr_vals.push_back(f5->GetParError(1));
      B_vals.push_back(f5->GetParameter(2)); Berr_vals.push_back(f5->GetParError(2));
      // converting fit fn to a hostogram
      TH1F *hf5 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f5,hf5);
      ho5[0]->Draw(); c5->Update(); 
      // grabbing statbox of the fitted histo
      TPaveStats *st5 = (TPaveStats*)ho5[0]->FindObject("stats");
      ho5[0]->SetBit(TH1::kNoStats); // Sets up the stat box for later modification
      // getting pads for pull plot
      std::vector<TPad*> p5 = util_pd::GetPadsForPullPlot(c5);
      //
      // preparing the pad for data/MC fit
      //
      p5[0]->cd();
      // drawing all the histograms
      //ho5[0]->Draw(); customize_ht(ho5[0]); customize_dx(ho5[0]);
      h_dxHCAL_data->Draw("E"); util_pd::customize_data(h_dxHCAL_data);
      hf5->Draw("same HIST"); util_pd::customize_gfit(hf5,1);
      ho5[4]->Draw("same HIST"); util_pd::customize_psig(ho5[4],1);
      ho5[5]->Draw("same HIST"); util_pd::customize_nsig(ho5[5],1);
      ho5[2]->Draw("same HIST"); util_pd::customize_hbg(ho5[2],1);
      // calculating yields
      std::vector<double> yo5; util_pd::GetYields(f5,ho5[4],ho5[5],ho5[2],yo5);
      pCnt.push_back(yo5[0]); pCnt_err.push_back(yo5[1]);
      nCnt.push_back(yo5[2]); nCnt_err.push_back(yo5[3]);
      bgCnt.push_back(yo5[4]); bgCnt_err.push_back(yo5[5]);
      // redrawing the stat box
      st5->SetX1NDC(0.62); st5->SetX2NDC(0.9); st5->SetY2NDC(0.9);
      // Modifying it to add yield ratio
      double yRatio5 = f5->GetParameter(1)*R_MC_nofit;
      double yRatio5_err = f5->GetParError(1)*R_MC_nofit;
      TText *t5 = st5->AddText(Form("R_{Yield} = %.4f #pm %.4f",yRatio5,yRatio5_err));
      customize_text(t5);
      st5->Draw("same");
      // drawing a legend
      TLegend *l5=new TLegend(0.10,0.6,0.36,0.9);
      l5->SetTextFont(42);
      l5->AddEntry(h_dxHCAL_data,"Data","p");
      //l5->AddEntry(f5,"Fit (MC + poly. bg.)","l");
      l5->AddEntry(hf5,"Fit (QE MC + bg.)","lf");
      l5->AddEntry(ho5[4],"p signal (from MC)","lf");
      l5->AddEntry(ho5[5],"n signal (from MC)","lf");
      l5->AddEntry(ho5[2],Form("Bg. (poly. of order %d)",Opoly),"lf");
      l5->AddEntry(ho5[3],"Residual","p");
      if (is_vary_cut) AddCutToLegend(l5,cuts_2[i].c_str());
      //l5->SetFillStyle(0); // makes legend box transparent
      l5->Draw();
      //
      // preparing the pad for residual
      //  
      p5[1]->cd();
      ho5[3]->Draw(); util_pd::customize_residual(ho5[3]);
      // drawing a horizontal line at y = 0
      util_pd::DrawZeroLine(p5[1],dx_fit_range[0],dx_fit_range[1]);
      // *******************
      //c5->SaveAs(Form("%s_c5_%d.png",outPNG.Data(),i)); 

      // ** --
      // // Canvas 5 : Fitting w/ polynomial background using side band method
      // // Steps: 1. Subrtact bg using sideband method, 2. Perform data/MC fit w/o bg
      // TCanvas *c5 = util_pd::TC("c5",1,1); gStyleFitCanvas();
      // c5->cd(); 
      // // Let's perform the sideband fit first
      // //vector<double> reject_points{-1.2,0.3};
      // vector<TH1F*> hosb5;
      // TF1* bgsb5 = fit::fit_1pbg_SB(dx_fit_range,
      // 				    reject_points,
      // 				    Opoly,
      // 				    fit::GetFitParams(f3),
      // 				    h_dxHCAL_data,
      // 				    hosb5);
      // hosb5[0]->Draw(); c5->Update(); 
      // // grabbing statbox of the fitted histo
      // TPaveStats *stsb5 = (TPaveStats*)hosb5[0]->FindObject("stats"); 
      // // hosb5[0]->Draw(); customize_ht(hosb5[0]); customize_dx(hosb5[0]);
      // // hosb5[1]->Draw("same ep"); customize_hs(hosb5[1]); customize_dx(hosb5[1]);
      // // hosb5[2]->Draw("same"); util_pd::customize_hbg(hosb5[2],1); customize_dx(hosb5[2]);

      // // Now, let's fit the bg subtracted signal using MC signals
      // vector<TH1F*> ho5;
      // TF1 *f5 = fit::fit_2hs_nbg_THI(dx_fit_range,
      // 				     hosb5[1],h_dxHCAL_simu_p,h_dxHCAL_simu_n,
      // 				     ho5);
      // // grabbing fit params for future use
      // chi2.push_back(f5->GetChisquare()); NDF.push_back(f5->GetNDF());
      // R_vals.push_back(f5->GetParameter(1)); Rerr_vals.push_back(f5->GetParError(1));
      // B_vals.push_back(bgsb5->GetParameter(0)); Berr_vals.push_back(bgsb5->GetParError(0));
      // // converting fit fn to a hostogram
      // TH1F *hf5 = (TH1F*)h_dxHCAL_data->Clone(); util_pd::TF1toTH1F(f5,hf5);
      // ho5[0]->Draw(); c5->Update(); 
      // // hosb5[0]->Draw("same"); customize_ht(hosb5[0]); customize_dx(hosb5[0]);
      // // hosb5[1]->Draw("same HIST"); customize_hs(hosb5[1]); customize_dx(hosb5[1]);
      // // hosb5[2]->Draw("same HIST"); util_pd::customize_hbg(hosb5[2],1); customize_dx(hosb5[2]);

      // // grabbing statbox of the fitted histo
      // TPaveStats *st5 = (TPaveStats*)ho5[0]->FindObject("stats");
      // // getting pads for pull plot
      // std::vector<TPad*> p5 = util_pd::GetPadsForPullPlot(c5);
      // //
      // // preparing the pad for data/MC fit
      // //
      // p5[0]->cd();
      // // drawing all the histograms
      // hosb5[0]->Draw(); customize_ht(hosb5[0]); customize_dx(hosb5[0]);
      // hosb5[2]->Draw("same"); util_pd::customize_hbg(hosb5[2],1); customize_dx(hosb5[2]);    
      // // h_dxHCAL_data->Draw("E"); util_pd::customize_data(h_dxHCAL_data);
      // hosb5[1]->Draw("same E"); util_pd::customize_data(hosb5[1]); hosb5[1]->SetStats(0);
      // hf5->Draw("same HIST"); util_pd::customize_gfit(hf5,1);
      // ho5[3]->Draw("same HIST"); util_pd::customize_psig(ho5[3],1);
      // ho5[4]->Draw("same HIST"); util_pd::customize_nsig(ho5[4],1);
      // //ho5[2]->Draw("same HIST"); util_pd::customize_hbg(ho5[2],1);
      // // redrawing the stat boxes
      // stsb5->SetX1NDC(0.62); stsb5->SetX2NDC(0.9); stsb5->SetY2NDC(0.9);
      // stsb5->Draw("same");
      // st5->SetX1NDC(0.62); st5->SetX2NDC(0.9); st5->SetY2NDC(0.6);
      // st5->Draw("same");
      // // drawing a legend
      // TLegend *l5=new TLegend(0.10,0.6,0.36,0.9);
      // l5->SetTextFont(42);
      // l5->AddEntry(hosb5[0],"Data","l");
      // l5->AddEntry(hosb5[2],Form("Bg w/ SB fit (poly. ord. %d)",Opoly),"lf");
      // //l5->AddEntry(f5,"Fit (MC + poly. bg.)","l");
      // l5->AddEntry(hosb5[1],"Data (bg subtracted)","p");
      // l5->AddEntry(hf5,"Fit (QE MC + bg.)","lf");
      // l5->AddEntry(ho5[3],"p signal (from MC)","lf");
      // l5->AddEntry(ho5[4],"n signal (from MC)","lf");
      // //l5->AddEntry(ho5[2],"Residual","p");
      // if (is_vary_cut) AddCutToLegend(l5,cuts_2[i].c_str());
      // //l5->SetFillStyle(0); // makes legend box transparent
      // l5->Draw();
      // // //
      // // // preparing the pad for residual
      // // p5[1]->cd();
      // // ho5[2]->Draw(); util_pd::customize_residual(ho5[2]);
      // // // drawing a horizontal line at y = 0
      // // util_pd::DrawZeroLine(p5[1],dx_fit_range[0],dx_fit_range[1]);
      // ----

      // Writing out fit parameters
      std::cout << "\n--- Reporting fit params ---\n";
      std::cout << "cut,min,max,R0,R0err,R1,R1err,R2,R2err,R3,R3err,R4,R4err\n";
      std::cout << cuts_2[i] << "," << minval[i] << "," << maxval[i] << ",";
      outdata << cuts_2[i] << "," << minval[i] << "," << maxval[i] << "," << R_MC_nofit << "," << R_MC_nofit_err << "," << Rp_MC << ",";
      for(size_t i=0; i < R_vals.size(); i++){
	std::cout << R_vals[i] << "," << Rerr_vals[i] << ",";
	//outdata << Form("%.1f,%.1f,%.4f,%.4f,%.4f,%.4f,",chi2[i],NDF[i],R_vals[i],Rerr_vals[i],B_vals[i],Berr_vals[i]);
	outdata << Form("%.1f,%.0f,",chi2[i],NDF[i]);
	outdata << Form("%.4f,%.4f,",R_vals[i],Rerr_vals[i]);
	if (i>0) {
	  outdata << Form("%.4f,%.4f,",B_vals[i],Berr_vals[i]);
	  outdata << Form("%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,",pCnt[i],pCnt_err[i],nCnt[i],nCnt_err[i],bgCnt[i],bgCnt_err[i]);
	}
      }
      // outdata << Form("%.0f,%.0f,%.0f,%.0f,%.0f,%.0f",yo2[0],yo2[1],yo2[2],yo2[3],yo2[4],yo2[5]);
      outdata << "\n";
      std::cout << "\n------\n";

      // further customization of the data histo
      FurtherCustoizeDataHisto(h_dxHCAL_data,cuts_for_signal_data);

      // writing out the canvases
      c0->Update(); c0->Write(); if (!is_vary_cut&&i==0) c0->SaveAs(Form("%s[",outPlot.Data())); 
      c0->SaveAs(Form("%s",outPlot.Data())); 
      c1->Update(); c1->Write(); c1->SaveAs(Form("%s",outPlot.Data())); 
      c2->Update(); c2->Write(); c2->SaveAs(Form("%s",outPlot.Data())); 
      c3->Update(); c3->Write(); c3->SaveAs(Form("%s",outPlot.Data())); 
      c4->Update(); c4->Write(); c4->SaveAs(Form("%s",outPlot.Data())); 
      c5->Update(); c5->Write(); c5->SaveAs(Form("%s",outPlot.Data())); 
      if (i==iter-1) c5->SaveAs(Form("%s]",outPlot.Data())); 
    } // QE
  } // for, cut vairation

  if (is_vary_cut) {
    // creating infinite loop gif
    cCut->SaveAs(Form("%s++10",outGIF.Data())); 
    // writing additional memory to file
    // Canvas : Plotting cut variation summary
    TCanvas *cGist = util_pd::TC("cGist",1,1); cGist->cd();
    cGist->SetBottomMargin(0.3);
    customize_hsummary(hgist_cv_2,cuts_2);
    hgist_cv_2->Draw(); hgist_cv_2->Write();
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
