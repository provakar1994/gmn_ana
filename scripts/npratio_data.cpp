/*
  Script to look at the stability of n/p ratio
*/

#include "gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

int npratio_data(const char *configfilename) {

  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

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
  ROOT::RDataFrame data_rdf("Tout",Form("pdout/%s_qelas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfprefix.c_str(),conf,sbsmag,model,pass));

  // Applying cuts
  std::string global_cut = jmgr->GetValueFromSubKey_str(key,"global_cut");
  std::string proton_cut = jmgr->GetValueFromSubKey_str(key,"proton_cut");
  std::string neutron_cut = jmgr->GetValueFromSubKey_str(key,"neutron_cut");
  std::string good_p_ev_cut = global_cut+"&&"+proton_cut;
  std::string good_n_ev_cut = global_cut+"&&"+neutron_cut;
  auto p_rdf = data_rdf.Filter(good_p_ev_cut.c_str());
  auto n_rdf = data_rdf.Filter(good_n_ev_cut.c_str());
  
  // fiducial cut
  bool apply_fidu_cut = jmgr->GetValueFromSubKey<int>(key,"apply_fidu_cut");
  std::vector<double> AR_w; jmgr->GetVectorFromSubKey<double>(key,"AR_width_x_y",AR_w);
  double sbs_kick = jmgr->GetValueFromSubKey<double>(key,"sbs_kick");
  std::vector<double> SM_w; jmgr->GetVectorFromSubKey<double>(key,"SM_width_xp_xn_y",SM_w);
  std::vector<double> hcal_AR = cut::hcal_active_area_data(AR_w[0],AR_w[1],pass); 
  std::vector<double> hcal_SM = cut::hcal_safety_margin(SM_w[0],SM_w[1],SM_w[2],hcal_AR);
  // forming the fiducial cut
  auto fiduCut = [&](double x,double y,double xExp,double yExp) {
    return cut::inHCAL_activeA(x,y,hcal_AR) && cut::inHCAL_safety_margin("LD2",xExp,yExp,sbs_kick,hcal_SM);
  };

  // xHCAL_exp
  std::vector<double> h_xexp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_xexp_lim",h_xexp_lim);
  std::vector<double> h_xexp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_xexp_lim",h_xexp_fitR);
  TH1F *h_xexp_p = (TH1F*)p_rdf.Histo1D({"h_xexp_p","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xHCAL_exp")->Clone();
  TH1F *h_xexp_n = (TH1F*)n_rdf.Histo1D({"h_xexp_n","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xHCAL_exp")->Clone();
  TH1F *h_xexp_npratio = new TH1F("h_xexp_npratio","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]);
  h_xexp_p->Sumw2();
  h_xexp_n->Sumw2();
  h_xexp_npratio->Divide(h_xexp_n,h_xexp_p);
  //
  TCanvas *cxexp = util_pd::TC("cxexp",2,2);
  cxexp->cd(1);
  h_xexp_p->Draw();
  cxexp->cd(2);
  h_xexp_n->Draw();
  cxexp->cd(3);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  //h_xexp_npratio->SetStats(0);
  h_xexp_npratio->Draw("E");
  h_xexp_npratio->GetYaxis()->SetRangeUser(0,1.0);
  //
  TF1 *fxexp = new TF1("fxexp","[0]",h_xexp_fitR[0],h_xexp_fitR[1]);
  h_xexp_npratio->Fit("fxexp","R");
  fxexp->Draw("same");
  //--

  // yHCAL_exp
  std::vector<double> h_yexp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_yexp_lim",h_yexp_lim);
  std::vector<double> h_yexp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_yexp_lim",h_yexp_fitR);
  TH1F *h_yexp_p = (TH1F*)p_rdf.Histo1D({"h_yexp_p","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yexp_n = (TH1F*)n_rdf.Histo1D({"h_yexp_n","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yexp_npratio = new TH1F("h_yexp_npratio","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]);
  h_yexp_p->Sumw2();
  h_yexp_n->Sumw2();
  h_yexp_npratio->Divide(h_yexp_n,h_yexp_p);
  //
  TCanvas *cyexp = util_pd::TC("cyexp",2,2);
  cyexp->cd(1);
  h_yexp_p->Draw();
  cyexp->cd(2);
  h_yexp_n->Draw();
  cyexp->cd(3);
  gStyle->SetOptFit(1);
  gStyle->SetOptStat(0);
  //h_yexp_npratio->SetStats(0);
  h_yexp_npratio->Draw("E");
  h_yexp_npratio->GetYaxis()->SetRangeUser(0,1.0);
  //
  TF1 *fyexp = new TF1("fyexp","[0]",h_yexp_fitR[0],h_yexp_fitR[1]);
  h_yexp_npratio->Fit("fyexp","R");
  fyexp->Draw("same");
  //--

  // visualizing n and p spots
  TCanvas *cdxdy = util_pd::TC("cdxdy",1,2);
  cdxdy->cd(1);
  TH2F *h_dxdy_n = (TH2F*)n_rdf.Histo2D({"h_dxdy_n","",200,-1.5,1.5,200,-3,2},"dy","dx")->Clone();
  h_dxdy_n->Draw("colz");
  cdxdy->cd(2);
  TH2F *h_dxdy_p = (TH2F*)p_rdf.Histo2D({"h_dxdy_p","",200,-1.5,1.5,200,-3,2},"dy","dx")->Clone();
  h_dxdy_p->Draw("colz");
  //--
  
  return 0;
}
