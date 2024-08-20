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

void custom_statbox_effi(TPaveStats *st) {
  st->SetBit(TH1::kNoStats);
  st->SetX1NDC(0.25); st->SetY1NDC(0.25); st->SetX2NDC(0.85); st->SetY2NDC(0.45);
}

int pDE_map_simu(const char *configfilename) {

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
  std::string gen = jmgr->GetValueFromSubKey_str(key,"generator");
  std::string sfprefix = jmgr->GetValueFromSubKey_str(key,"simu_file_prefix");

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame simu_rdf_raw("Tout",Form("siout/%s_elas_ana_%s_sbs%d_sbs%dp_model%d.root",sfprefix.c_str(),gen.c_str(),conf,sbsmag,model));
  auto simu_rdf = simu_rdf_raw.Define("xExp_shifted","xHCAL_exp-p_def").Define("dx_def","dx+p_def");
 
  // Applying cuts
  std::string global_cut = jmgr->GetValueFromSubKey_str(key,"global_cut");
  std::string w2_cut = jmgr->GetValueFromSubKey_str(key,"w2_cut");
  std::string SM_cut_x = jmgr->GetValueFromSubKey_str(key,"SM_cut_x");
  std::string SM_cut_y = jmgr->GetValueFromSubKey_str(key,"SM_cut_y");
  // 
  std::string earm_cut = global_cut + "&&" + w2_cut; // w/o SM cuts
  std::string harm_cut = jmgr->GetValueFromSubKey_str(key,"harm_cut");
  std::string eNharm_cut = earm_cut + "&&" + harm_cut;
  //
  std::string earm_cut_wSMx = global_cut + "&&" + w2_cut + "&&" + SM_cut_x;
  std::string eNharm_cut_wSMx = earm_cut_wSMx + "&&" + harm_cut;
  std::string earm_cut_wSMy = global_cut + "&&" + w2_cut + "&&" + SM_cut_y;
  std::string eNharm_cut_wSMy = earm_cut_wSMy + "&&" + harm_cut;
  //
  std::string eNantiharm_cut = earm_cut + "&&!(" + harm_cut + ")";
  std::string eNantiharm_cut_wSM = earm_cut + SM_cut_x + SM_cut_y + "&&!(" + harm_cut + ")";
  //auto earm_rdf = simu_rdf.Filter(earm_cut.c_str());
  //auto eNharm_rdf = simu_rdf.Filter(eNharm_cut.c_str());
  
  // // fiducial cut
  // bool apply_fidu_cut = jmgr->GetValueFromSubKey<int>(key,"apply_fidu_cut");
  // std::vector<double> AR_w; jmgr->GetVectorFromSubKey<double>(key,"AR_width_x_y",AR_w);
  // double sbs_kick = jmgr->GetValueFromSubKey<double>(key,"sbs_kick");
  // std::vector<double> SM_w; jmgr->GetVectorFromSubKey<double>(key,"SM_width_xp_xn_y",SM_w);
  // std::vector<double> hcal_AR = cut::hcal_active_area_data(AR_w[0],AR_w[1],2); 
  // std::vector<double> hcal_SM = cut::hcal_safety_margin(SM_w[0],SM_w[1],SM_w[2],hcal_AR);
  // // forming the fiducial cut
  // auto fiduCut = [&](double x,double y,double xExp,double yExp) {
  //   return cut::inHCAL_activeA(x,y,hcal_AR) && cut::inHCAL_safety_margin("LH2",xExp,yExp,sbs_kick,hcal_SM);
  // };

  // output file
  std::string filebase = jmgr->GetValueFromSubKey_str(key,"outfile_prefix");
  filebase = filebase.empty() ? "" : filebase + "_";
  std::string outfilebase = "siout/pDE/" + sfprefix + "_pDE_simu_" + Form("%s_sbs%d_sbs%dp_model%d",gen.c_str(),conf,sbsmag,model);
  TFile *fout = new TFile(Form("%s.root",outfilebase.c_str()), "RECREATE");
  
  // call the canvas customizer
  PlotCustomizer pcust;
 
  // **************
  // visualizing n and p spots
  // **************
  TString xtitledxdy = "#font[32]{#Deltay} (m)";
  TString ytitledxdy = "#font[32]{#Deltax} (m)";
  TCanvas *cdxdy = util_pd::TC("cdxdy",2,2);
  //
  cdxdy->cd(1);
  gPad->SetLogz();
  gStyle->SetPalette(kRainbow);
  TH2F *h_dxdy_earm = (TH2F*)simu_rdf.Filter(earm_cut.c_str()).Filter("eHCAL>0").Histo2D({"h_dxdy_earm","",200,-1.5,1.5,200,-3,2},"dy","dx","weight")->Clone();
  util_pd::SetAxTitles(h_dxdy_earm,ytitledxdy,xtitledxdy);
  h_dxdy_earm->Draw("colz");
  //
  cdxdy->cd(2);
  gPad->SetLogz();
  gStyle->SetPalette(kRainbow);
  TH2F *h_dxdy_eNharm = (TH2F*)simu_rdf.Filter(eNharm_cut.c_str()).Filter("eHCAL>0").Histo2D({"h_dxdy_eNharm","",200,-1.5,1.5,200,-3,2},"dy","dx","weight_effic")->Clone();
  util_pd::SetAxTitles(h_dxdy_eNharm,ytitledxdy,xtitledxdy);
  h_dxdy_eNharm->Draw("colz");
  //
  cdxdy->cd(3);
  std::vector<double> h_w2_lim; jmgr->GetVectorFromSubKey<double>(key,"h_w2_lim",h_w2_lim);
  TH1F *h_w2_earm = (TH1F*)simu_rdf.Filter(global_cut.c_str()).Histo1D({"h_w2_earm","",(int)h_w2_lim[0],h_w2_lim[1],h_w2_lim[2]},"W2","weight")->Clone();
  TH1F *h_w2_eNharm = (TH1F*)simu_rdf.Filter(global_cut.c_str()).Filter(harm_cut.c_str()).Histo1D({"h_w2_eNharm","",(int)h_w2_lim[0],h_w2_lim[1],h_w2_lim[2]},"W2","weight_effic")->Clone();
  util_pd::SetAxTitles(h_w2_earm,"","#font[32]{W^{2}} (GeV^{2})");
  custom_denom(h_w2_earm);
  h_w2_earm->Draw("HIST");
  custom_num(h_w2_eNharm);
  h_w2_eNharm->Draw("HIST same");
  TLegend *lw2 = new TLegend(0.15,0.75,0.55,0.9);
  lw2->AddEntry(h_w2_earm,"All","lp");
  lw2->AddEntry(h_w2_eNharm,"With #font[32]{#Deltax-#Deltay} Cut","lp");
  lw2->SetFillStyle(0);
  lw2->Draw();
  //
  cdxdy->cd(4);
  std::string dx_w_def_cut = earm_cut + "&&eHCAL>0&&abs(dy)<0.3";
  TH1F *h_dx_w_def = (TH1F*)simu_rdf.Filter(dx_w_def_cut.c_str()).Histo1D({"h_dx_w_def","",200,-1,1},"dx_def","weight_effic")->Clone();
  h_dx_w_def->Draw("HIST");
  h_dx_w_def->GetXaxis()->SetTitle("#font[32]{#Deltax} + proton_deflection (m)");
  std::vector<double> hdx_fitR{-0.5,0.5,1.2,1.2};
  TF1 *fdx = fit::fit_1gs_nbg(hdx_fitR,h_dx_w_def);
  double dxM = fdx->GetParameter(1), dxMerr = fdx->GetParError(1);
  double dxS = fdx->GetParameter(2), dxSerr = fdx->GetParError(2);
  pcust.AddTitleText(Form("#mu = %.4f #pm %.4f, #sigma = %.4f",dxM,dxMerr,dxS),0.95);
  fdx->Draw("same");
  //
  pcust.customize_canvas(cdxdy);
  cdxdy->Update();
  cdxdy->Write();
  cdxdy->SaveAs(Form("%s.pdf[",outfilebase.c_str()));
  cdxdy->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  cdxdy->SaveAs(Form("%s_1.png",outfilebase.c_str()));
  //--


  // **************
  // Envelopes
  // **************
  std::vector<double> hcal_area = cut::hcal_active_area_data(0,0,2);
  TH2F *h2_xyexp_all_nodef = (TH2F*)simu_rdf.Filter(earm_cut.c_str()).Histo2D({"h2_xyexp_all_nodef","All (No Deflection)",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xHCAL_exp","weight")->Clone();
  TH2F *h2_xyexp_all = (TH2F*)simu_rdf.Filter(earm_cut.c_str()).Histo2D({"h2_xyexp_all","All",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xExp_shifted","weight")->Clone();
  TH2F *h2_xyexp_pass = (TH2F*)simu_rdf.Filter(eNharm_cut.c_str()).Histo2D({"h2_xyexp_pass","Passed HCAL",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xExp_shifted","weight_effic")->Clone();
  TH2F *h2_xyexp_fail = (TH2F*)simu_rdf.Filter(eNantiharm_cut.c_str()).Histo2D({"h2_xyexp_fail","Failed HCAL",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xExp_shifted","weight_effic")->Clone();
  // TH2F *h2_xyexp_pass = (TH2F*)simu_rdf.Filter(eNharm_cut.c_str()).Histo2D({"h2_xyexp_pass","Passed HCAL",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xExp_shifted","weight")->Clone();
  // TH2F *h2_xyexp_fail = (TH2F*)simu_rdf.Filter(eNantiharm_cut.c_str()).Histo2D({"h2_xyexp_fail","Failed HCAL",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xExp_shifted","weight")->Clone();  
  TCanvas *cenv = util_pd::TC("cenv",2,2);
  cenv->cd(1); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_all_nodef->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->cd(2); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_all->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->cd(3); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_pass->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->cd(4); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_fail->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->Write();
  cenv->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  cenv->SaveAs(Form("%s_2.png",outfilebase.c_str()));
  //--

  /*
    #############################
    ## Efficiency Calculations
    #############################
  */
  PlotCustomizer pcust_wStat{1,1};
  // xHCAL_exp *** -- \\//
  TString xtitlexexp = "#font[32]{x^{exp}_{HCAL}} (m)";
  std::vector<double> h_xexp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_xexp_lim",h_xexp_lim);
  std::vector<double> h_xexp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_xexp_fitR",h_xexp_fitR);
  TH1F *h_xexp_earm = (TH1F*)simu_rdf.Filter(earm_cut_wSMy.c_str()).Histo1D({"h_xexp_earm","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xExp_shifted","weight")->Clone();
  TH1F *h_xexp_eNharm = (TH1F*)simu_rdf.Filter(eNharm_cut_wSMy.c_str()).Histo1D({"h_xexp_eNharm","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xExp_shifted","weight_effic")->Clone();
   // TH1F *h_xexp_eNharm = (TH1F*)simu_rdf.Filter(eNharm_cut_wSMy.c_str()).Histo1D({"h_xexp_eNharm","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xExp_shifted","weight")->Clone();
  util_pd::SetAxTitles(h_xexp_earm,"",xtitlexexp);
  h_xexp_earm->SetStats(0);
  util_pd::SetAxTitles(h_xexp_eNharm,"",xtitlexexp);
  custom_denom(h_xexp_earm);
  custom_num(h_xexp_eNharm);
  TH1F *h_xexp_pDE = new TH1F("h_xexp_pDE","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]);
  util_pd::SetAxTitles(h_xexp_pDE,"Efficiency",xtitlexexp);
  custom_ratio(h_xexp_pDE);
  h_xexp_pDE->Divide(h_xexp_eNharm,h_xexp_earm);
  // Binomial error
  // for( int i=1; i<=h_xexp_pDE->GetNbinsX(); i++ ){
  //   double effi = h_xexp_pDE->GetBinContent(i);
  //   //prevent divide-by-zero errors
  //   double N = std::max(1.0,h_xexp_earm->GetBinContent(i));
  //   h_xexp_pDE->SetBinError(i,sqrt(effi*(1.0-effi)/N));
  // }
  for(int i = 1; i <= h_xexp_pDE->GetNbinsX(); i++) {
    double effi = h_xexp_pDE->GetBinContent(i);

    if (effi>0) {

      // Sum of weights and sum of squared weights
      double sum_w = h_xexp_earm->GetBinContent(i);  // This assumes h_xexp_earm holds the sum of weights
      double sum_w2 = h_xexp_earm->GetBinError(i);   // Assuming the bin errors in h_xexp_earm store sum of squared weights

      // Prevent divide-by-zero errors
      double N = std::max(1.0, sum_w);
      double weighted_error = sqrt(effi * (1.0 - effi) / N) * sqrt(sum_w2 / sum_w);

      h_xexp_pDE->SetBinError(i, weighted_error);
    }
  }    
  // yHCAL_exp *** -- \\//
  TString xtitleyexp = "#font[32]{y^{exp}_{HCAL}} (m)";
  std::vector<double> h_yexp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_yexp_lim",h_yexp_lim);
  std::vector<double> h_yexp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_yexp_fitR",h_yexp_fitR);
  TH1F *h_yexp_earm = (TH1F*)simu_rdf.Filter(earm_cut_wSMx.c_str()).Histo1D({"h_yexp_earm","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]},"yHCAL_exp","weight")->Clone();
  TH1F *h_yexp_eNharm = (TH1F*)simu_rdf.Filter(eNharm_cut_wSMx.c_str()).Histo1D({"h_yexp_eNharm","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]},"yHCAL_exp","weight_effic")->Clone();
  util_pd::SetAxTitles(h_yexp_earm,"",xtitleyexp);
  util_pd::SetAxTitles(h_yexp_eNharm,"",xtitleyexp);
  custom_denom(h_yexp_earm);
  custom_num(h_yexp_eNharm);
  TH1F *h_yexp_pDE = new TH1F("h_yexp_pDE","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]);
  util_pd::SetAxTitles(h_yexp_pDE,"Efficiency",xtitleyexp);
  custom_ratio(h_yexp_pDE);
  h_yexp_pDE->Divide(h_yexp_eNharm,h_yexp_earm);
  // Binomial error
  // for( int i=1; i<=h_yexp_pDE->GetNbinsX(); i++ ){
  //   double effi = h_yexp_pDE->GetBinContent(i);
  //   //prevent divide-by-zero errors
  //   double N = std::max(1.0,h_yexp_earm->GetBinContent(i));
  //   h_yexp_pDE->SetBinError(i,sqrt(effi*(1.0-effi)/N));
  // }
  for(int i = 1; i <= h_yexp_pDE->GetNbinsX(); i++) {
    double effi = h_yexp_pDE->GetBinContent(i);

    if (effi>0) {

      // Sum of weights and sum of squared weights
      double sum_w = h_yexp_earm->GetBinContent(i);  // This assumes h_yexp_earm holds the sum of weights
      double sum_w2 = h_yexp_earm->GetBinError(i);   // Assuming the bin errors in h_yexp_earm store sum of squared weights

      // Prevent divide-by-zero errors
      double N = std::max(1.0, sum_w);
      double weighted_error = sqrt(effi * (1.0 - effi) / N) * sqrt(sum_w2 / sum_w);

      h_yexp_pDE->SetBinError(i, weighted_error);
    }
  }  
  
  // 
  TCanvas *cxyexp = util_pd::TC("cxyexp",2,2);
  cxyexp->cd(1);
  h_xexp_earm->Draw();
  h_xexp_eNharm->Draw("same");
  TLegend *lxexp = new TLegend(0.25,0.25,0.85,0.45);
  lxexp->AddEntry(h_xexp_earm,"All","lp");
  lxexp->AddEntry(h_xexp_eNharm,"With #font[32]{#Deltax-#Deltay} Cut","lp");
  lxexp->SetFillStyle(0);
  lxexp->Draw();
  //
  cxyexp->cd(2);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_xexp_pDE->Draw("E");
  h_xexp_pDE->GetYaxis()->SetRangeUser(0,1.2);
  TF1 *fxexp = new TF1("fxexp","pol0",h_xexp_fitR[0],h_xexp_fitR[1]);
  fxexp->SetNpx(2000);
  h_xexp_pDE->Fit("fxexp","R");
  cxyexp->Update();
  TPaveStats *stxexp = (TPaveStats*)h_xexp_pDE->FindObject("stats");
  custom_statbox_effi(stxexp);
  fxexp->Draw("same");
  //
  cxyexp->cd(3);
  h_yexp_earm->Draw();
  h_yexp_eNharm->Draw("same");
  TLegend *lyexp = new TLegend(0.25,0.25,0.85,0.45);
  lyexp->AddEntry(h_yexp_earm,"All","lp");
  lyexp->AddEntry(h_yexp_eNharm,"With #font[32]{#Deltax-#Deltay} Cut","lp");
  lyexp->SetFillStyle(0);
  lyexp->Draw();  
  //
  cxyexp->cd(4);
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_yexp_pDE->Draw("E");
  h_yexp_pDE->GetYaxis()->SetRangeUser(0,1.2);
  TF1 *fyexp = new TF1("fyexp","[0]",h_yexp_fitR[0],h_yexp_fitR[1]);
  fyexp->SetRange(h_yexp_fitR[0], h_yexp_fitR[1]);
  fyexp->SetNpx(2000);
  h_yexp_pDE->Fit("fyexp","R");
  cxyexp->Update();
  TPaveStats *styexp = (TPaveStats*)h_yexp_pDE->FindObject("stats");
  custom_statbox_effi(styexp);  
  fyexp->Draw("same");
  //
  pcust_wStat.customize_canvas(cxyexp);
  cxyexp->Update();
  cxyexp->Write();
  cxyexp->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  cxyexp->SaveAs(Form("%s_3.png",outfilebase.c_str()));
  //--


  // Efficiency Map
  // TH2F *h2_xyexp_earm = (TH2F*)simu_rdf.Filter(earm_cut.c_str()).Histo2D({"h2_xyexp_earm","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2],int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"yHCAL_exp","xExp_shifted")->Clone();
  // TH2F *h2_xyexp_eNharm = (TH2F*)simu_rdf.Filter(eNharm_cut.c_str()).Histo2D({"h2_xyexp_eNharm","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2],int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"yHCAL_exp","xExp_shifted")->Clone();
  // TH2F *h2_effi_map = new TH2F("h2_effi_map","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2],int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]);
  TH2F *h2_xyexp_earm = (TH2F*)simu_rdf.Filter(earm_cut.c_str()).Histo2D({"h2_xyexp_earm","",65,-1.25,1.25,126,-3.25,1.75},"yHCAL_exp","xExp_shifted","weight")->Clone();
  TH2F *h2_xyexp_eNharm = (TH2F*)simu_rdf.Filter(eNharm_cut.c_str()).Histo2D({"h2_xyexp_eNharm","",65,-1.25,1.25,126,-3.25,1.75},"yHCAL_exp","xExp_shifted","weight_effic")->Clone();
  TH2F *h2_effi_map = new TH2F("h2_effi_map","",65,-1.25,1.25,126,-3.25,1.75);
  h2_effi_map->Divide(h2_xyexp_eNharm,h2_xyexp_earm);
  // Binomial error
  // for( int i=1; i<=h2_effi_map->GetNbinsX(); i++ ){
  //   for( int j=1; j<=h2_effi_map->GetNbinsY(); j++ ){
  //     int bin = h2_effi_map->GetBin(i,j);
  //     double effi = h2_effi_map->GetBinContent(bin);
  //     double N = std::max(1.0,h2_xyexp_earm->GetBinContent(bin));
  //     h2_effi_map->SetBinError(bin,sqrt(effi*(1.0-effi)/N));
  //   }
  // }
  for( int i=1; i<=h2_effi_map->GetNbinsX(); i++ ){
    for( int j=1; j<=h2_effi_map->GetNbinsY(); j++ ){
      int bin = h2_effi_map->GetBin(i,j);
      double effi = h2_effi_map->GetBinContent(bin);

      if (effi>0) {

	// Sum of weights and sum of squared weights
	double sum_w = h2_xyexp_earm->GetBinContent(i);  // This assumes h_yexp_earm holds the sum of weights
	double sum_w2 = h2_xyexp_earm->GetBinError(i);   // Assuming the bin errors in h_yexp_earm store sum of squared weights

	// Prevent divide-by-zero errors
	double N = std::max(1.0, sum_w);
	double weighted_error = sqrt(effi * (1.0 - effi) / N) * sqrt(sum_w2 / sum_w);

	h2_effi_map->SetBinError(i, weighted_error);
      }

    }
  }
  TCanvas *cefmap = util_pd::TC("cefmap",1,1);
  cefmap->cd(); gStyle->SetPalette(kRainbow); gStyle->SetNumberContours(50);
  h2_effi_map->Draw("colz");
  pcust.customize_canvas(cefmap);
  cefmap->Update();
  cefmap->Write();
  cefmap->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  cefmap->SaveAs(Form("%s.pdf]",outfilebase.c_str()));
  cefmap->SaveAs(Form("%s_4.png",outfilebase.c_str()));

  fout->Write();
  
  return 0;
}
