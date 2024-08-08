/*
  Script to look at the stability of n/p ratio
*/

#include "gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

void custom_num(TH1F *h) {
  h->SetLineColor(kRed);
  h->SetLineWidth(2);
}

void custom_denom(TH1F *h) {
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
}

void custom_ratio(TH1F *h) {
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(8);
}

int pDE_data(const char *configfilename) {

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  TH1::SetDefaultSumw2();
 
  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  char const * key = "prod";

  // reading in proper data and simu output files
  int conf = jmgr->GetValueFromSubKey<int>(key,"SBS_config");
  int sbsmag = jmgr->GetValueFromSubKey<int>(key,"SBS_magnet_percent");
  int model = jmgr->GetValueFromSubKey<int>(key,"model");
  int pass = jmgr->GetValueFromSubKey<int>(key,"pass");
  std::string dfprefix = jmgr->GetValueFromSubKey_str(key,"data_file_prefix");

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/%s_elas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfprefix.c_str(),conf,sbsmag,model,pass));

  // Applying cuts
  std::string global_cut = jmgr->GetValueFromSubKey_str(key,"global_cut");
  std::string w2_cut = jmgr->GetValueFromSubKey_str(key,"w2_cut");
  std::string earm_cut = global_cut + "&&" + w2_cut;
  std::string harm_cut = jmgr->GetValueFromSubKey_str(key,"harm_cut");
  // std::string proton_cut = jmgr->GetValueFromSubKey_str(key,"proton_cut");
  // std::string neutron_cut = jmgr->GetValueFromSubKey_str(key,"neutron_cut");
  // std::string good_p_ev_cut = global_cut+"&&"+proton_cut;
  // std::string good_n_ev_cut = global_cut+"&&"+neutron_cut;
  std::string eNharm_cut = earm_cut + "&&" + harm_cut;
  auto earm_rdf = data_rdf.Filter(earm_cut.c_str());
  auto eNharm_rdf = data_rdf.Filter(eNharm_cut.c_str());
  
  // fiducial cut
  bool apply_fidu_cut = jmgr->GetValueFromSubKey<int>(key,"apply_fidu_cut");
  std::vector<double> AR_w; jmgr->GetVectorFromSubKey<double>(key,"AR_width_x_y",AR_w);
  double sbs_kick = jmgr->GetValueFromSubKey<double>(key,"sbs_kick");
  std::vector<double> SM_w; jmgr->GetVectorFromSubKey<double>(key,"SM_width_xp_xn_y",SM_w);
  std::vector<double> hcal_AR = cut::hcal_active_area_data(AR_w[0],AR_w[1],pass); 
  std::vector<double> hcal_SM = cut::hcal_safety_margin(SM_w[0],SM_w[1],SM_w[2],hcal_AR);
  // forming the fiducial cut
  auto fiduCut = [&](double x,double y,double xExp,double yExp) {
    return cut::inHCAL_activeA(x,y,hcal_AR) && cut::inHCAL_safety_margin("LH2",xExp,yExp,sbs_kick,hcal_SM);
  };

  // xHCAL_exp
  std::vector<double> h_xexp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_xexp_lim",h_xexp_lim);
  std::vector<double> h_xexp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_xexp_lim",h_xexp_fitR);
  TH1F *h_xexp_earm = (TH1F*)earm_rdf.Histo1D({"h_xexp_earm","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xHCAL_exp")->Clone();
  TH1F *h_xexp_eNharm = (TH1F*)eNharm_rdf.Histo1D({"h_xexp_eNharm","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xHCAL_exp")->Clone();
  custom_denom(h_xexp_earm);
  custom_num(h_xexp_eNharm);
  TH1F *h_xexp_pDE = new TH1F("h_xexp_pDE","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]);
  custom_ratio(h_xexp_pDE);
  h_xexp_pDE->Divide(h_xexp_eNharm,h_xexp_earm);
  //
  TCanvas *cxexp = util_pd::TC("cxexp",1,2);
  cxexp->cd(1);
  h_xexp_earm->Draw();
  h_xexp_eNharm->Draw("same");
  cxexp->cd(2);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  //h_xexp_pDE->SetStats(0);
  h_xexp_pDE->Draw("E");
  h_xexp_pDE->GetYaxis()->SetRangeUser(0,1.2);
  //
  TF1 *fxexp = new TF1("fxexp","pol0",h_xexp_fitR[0],h_xexp_fitR[1]);
  fxexp->SetNpx(2000);
  h_xexp_pDE->Fit("fxexp","R");
  fxexp->Draw("same");
  //--

  // yHCAL_exp
  std::vector<double> h_yexp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_yexp_lim",h_yexp_lim);
  std::vector<double> h_yexp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_yexp_lim",h_yexp_fitR);
  TH1F *h_yexp_earm = (TH1F*)earm_rdf.Histo1D({"h_yexp_earm","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yexp_eNharm = (TH1F*)eNharm_rdf.Histo1D({"h_yexp_eNharm","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]},"yHCAL_exp")->Clone();
  custom_denom(h_yexp_earm);
  custom_num(h_yexp_eNharm);
  TH1F *h_yexp_pDE = new TH1F("h_yexp_pDE","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]);
  custom_ratio(h_yexp_pDE);
  h_yexp_pDE->Divide(h_yexp_eNharm,h_yexp_earm);
  //
  TCanvas *cyexp = util_pd::TC("cyexp",1,2);
  cyexp->cd(1);
  h_yexp_earm->Draw();
  h_yexp_eNharm->Draw("same");
  cyexp->cd(2);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  h_yexp_pDE->Draw("E");
  h_yexp_pDE->GetYaxis()->SetRangeUser(0,1.2);
  //
  TF1 *fyexp = new TF1("fyexp","[0]",h_yexp_fitR[0],h_yexp_fitR[1]);
  fyexp->SetNpx(2000);
  h_yexp_pDE->Fit("fyexp","R");
  fyexp->Draw("same");
  //--

  // visualizing n and p spots
  TCanvas *cdxdy = util_pd::TC("cdxdy",2,2);
  cdxdy->cd(1);
  TH2F *h_dxdy_earm = (TH2F*)earm_rdf.Histo2D({"h_dxdy_earm","",200,-1.5,1.5,200,-3,2},"dy","dx")->Clone();
  h_dxdy_earm->Draw("colz");
  cdxdy->cd(2);
  TH2F *h_dxdy_eNharm = (TH2F*)eNharm_rdf.Histo2D({"h_dxdy_eNharm","",200,-1.5,1.5,200,-3,2},"dy","dx")->Clone();
  h_dxdy_eNharm->Draw("colz");
  cdxdy->cd(3);
  std::vector<double> h_w2_lim; jmgr->GetVectorFromSubKey<double>(key,"h_w2_lim",h_w2_lim);
  TH1F *h_w2_earm = (TH1F*)data_rdf.Filter(global_cut.c_str()).Histo1D({"h_w2_earm","",(int)h_w2_lim[0],h_w2_lim[1],h_w2_lim[2]},"W2")->Clone();
  TH1F *h_w2_eNharm = (TH1F*)data_rdf.Filter(global_cut.c_str()).Filter(harm_cut.c_str()).Histo1D({"h_w2_eNharm","",(int)h_w2_lim[0],h_w2_lim[1],h_w2_lim[2]},"W2")->Clone();
  custom_denom(h_w2_earm);
  h_w2_earm->Draw("HIST");
  custom_num(h_w2_eNharm);
  h_w2_eNharm->Draw("HIST same");
  //--  

  
  return 0;
}
