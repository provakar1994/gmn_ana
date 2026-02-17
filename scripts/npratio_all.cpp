/*
  Script to look at the stability of n/p ratio
*/

#include "gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

bool const print_canv_ppar = false; // if true, prints canvas per parameter
bool const print_png = false; // if true, prints png version of all plots

double xNDCoffset = 0.05;
double y1NDCoffset = 0.05;

void drawcutrange (double x1, double x2, bool iseffi=false, bool isthresh=false) {
  double x1NDC = util_pd::GetxNDC(x1)+xNDCoffset;
  double x2NDC = util_pd::GetxNDC(x2)+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  double y2NDC = iseffi ? 0.7 : 0.9;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,y2NDC);
  if (!isthresh) {
    TLine L2;
    L2.SetLineColor(2); L2.SetLineWidth(4); //L2.SetLineStyle(9);
    L2.DrawLineNDC(x2NDC,y1NDC,x2NDC,y2NDC);
  }
}

void custom_extra(TH1F *h) {
  h->GetXaxis()->SetTitleSize(0.075);
  h->GetXaxis()->SetTitleOffset(0.85);
  h->GetXaxis()->SetLabelSize(0.06);

  h->GetYaxis()->SetTitleSize(0.09);
  h->GetYaxis()->SetTitleOffset(0.8);
  h->GetYaxis()->SetLabelSize(0.06);

}

void custom_num(TH1F *h) {
  h->SetLineColor(kGreen+2);
  h->SetLineWidth(2);
  h->SetStats(0);
}

void custom_denom(TH1F *h) {
  h->SetLineColor(kBlue);
  h->SetLineWidth(2);
  h->SetStats(0);
}

void custom_ratio(TH1F *h) {
  h->GetYaxis()->SetRangeUser(0,1.5);
  h->SetLineColor(kBlack);
  h->SetLineWidth(2);
  h->SetMarkerColor(kBlack);
  h->SetMarkerStyle(8);
}

void custom_statbox_effi(TPaveStats *st) {
  st->SetBit(TH1::kNoStats);
  st->SetX1NDC(0.25); st->SetY1NDC(0.7); st->SetX2NDC(0.85); st->SetY2NDC(0.9);
}

// void calc_binomial_error(TH1F *hdenom, TH1F *heffi) {
//   // Binomial error
//   for( int i=1; i<=heffi->GetNbinsX(); i++ ){
//     double effi = heffi->GetBinContent(i);
//     //prevent divide-by-zero errors
//     double N = std::max(1.0,hdenom->GetBinContent(i));
//     heffi->SetBinError(i,sqrt(effi*(1.0-effi)/N));
//   }  
// }

void calc_binomial_error(TH1F *hdenom, TH1F *heffi) {
    // Check that both histograms have the same number of bins
    if (hdenom->GetNbinsX() != heffi->GetNbinsX()) {
        std::cerr << "Error: Histograms must have the same number of bins." << std::endl;
        return;
    }

    // Binomial error calculation
    for (int i = 1; i <= heffi->GetNbinsX(); i++) {
        double effi = heffi->GetBinContent(i);
        double N = hdenom->GetBinContent(i);

        // Ensure there is no divide-by-zero and N is valid
        if (N > 0) {
            double error = sqrt(effi * (1.0 - effi) / N);
            heffi->SetBinError(i, error);
        } else {
            heffi->SetBinError(i, 0); // Set error to zero if N is zero or invalid
        }

        // Optional: Set efficiency to zero if denominator is zero to avoid misleading content
        if (N == 0) {
            heffi->SetBinContent(i, 0);
        }
    }
}

int npratio_all(const char *configfilename) {

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
  ROOT::RDataFrame data_rdf_raw("Tout",Form("pdout/%s_qelas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",dfprefix.c_str(),conf,sbsmag,model,pass));
  auto data_rdf = data_rdf_raw
    .Define("xHCAL_exp_def","xHCAL_exp-p_def")
    .Define("dx_def","dx+p_def")
    .Define("xbb","fpX-0.9*fpTh")
    .Define("ybb","fpY-0.9*fpPh");    

  // Applying cuts
  std::string globalcut = jmgr->GetValueFromSubKey_str(key,"globalcut");
  std::string vzcut = jmgr->GetValueFromSubKey_str(key,"vzcut");
  std::string EovPcut = jmgr->GetValueFromSubKey_str(key,"EovPcut");
  std::string trchi2cut = jmgr->GetValueFromSubKey_str(key,"trchi2cut");
  std::string ePScut = jmgr->GetValueFromSubKey_str(key,"ePScut");
  std::string W2cut = jmgr->GetValueFromSubKey_str(key,"W2cut");
  std::string xbbcut = jmgr->GetValueFromSubKey_str(key,"xbbcut");
  std::string ybbcut = jmgr->GetValueFromSubKey_str(key,"ybbcut");
  std::string eHCALcut = jmgr->GetValueFromSubKey_str(key,"eHCALcut");
  std::string coinTcut = jmgr->GetValueFromSubKey_str(key,"coinTcut");
  std::string coinTvar = coinTcut.find("TOF") != std::string::npos ? "coinT_ADC_TOF_c" : "coinT_ADC_c";
  //
  std::string SMcut_xp = jmgr->GetValueFromSubKey_str(key,"SMcut_xp");
  std::string SMcut_xn = jmgr->GetValueFromSubKey_str(key,"SMcut_xn");
  std::string SMcut_y = jmgr->GetValueFromSubKey_str(key,"SMcut_y");
  //
  std::string pcut = jmgr->GetValueFromSubKey_str(key,"pcut");
  std::string ncut = jmgr->GetValueFromSubKey_str(key,"ncut");
  //
  std::string allcut = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+xbbcut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string allearm = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+xbbcut+"&&"+ybbcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string allearmNOsmx = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+xbbcut+"&&"+ybbcut+"&&"+SMcut_y;  
  std::string allearmNOsm = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+xbbcut+"&&"+ybbcut;
  std::string allharm = eHCALcut+"&&"+coinTcut+"&&("+pcut+"||"+ncut+")";
  std::string alleNharm = allearm + "&&" + allharm;
  std::string eNharm = allearmNOsm + "&&" + allharm;
  std::string eNantiharm = allearmNOsm + "&&!(" + allharm + ")";
  //
  std::string onlyp = pcut+"&&!("+ncut+")";
  std::string onlyn = ncut+"&&!("+pcut+")";
  //  std::cout << eNantiharm << "\n";
  //
  std::string novz = globalcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+xbbcut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string noEovP = globalcut+"&&"+vzcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+xbbcut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string notrchi2ndf = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+ePScut+"&&"+W2cut+"&&"+xbbcut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string noePS = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+W2cut+"&&"+xbbcut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string noW2 = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+xbbcut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut;
  std::string noxbb = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut;
  std::string noybb = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string noeHCAL = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+ybbcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string nocoinT = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+SMcut_xp+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string nosmxp = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xn+"&&"+SMcut_y;
  std::string nosmxn = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_y;
  std::string nosmx = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_y;
  std::string nosmy = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut+"&&"+SMcut_xp+"&&"+SMcut_xn;
  std::string nosm = globalcut+"&&"+vzcut+"&&"+EovPcut+"&&"+trchi2cut+"&&"+ePScut+"&&"+W2cut+"&&"+ybbcut+"&&"+eHCALcut+"&&"+coinTcut;
  //
  std::string good_p_evcut = globalcut+"&&"+pcut;
  std::string good_n_evcut = globalcut+"&&"+ncut;
  auto prdf = data_rdf.Filter(good_p_evcut.c_str());
  auto nrdf = data_rdf.Filter(good_n_evcut.c_str());

  // fiducial cut
  std::vector<double> AR_w; jmgr->GetVectorFromSubKey<double>(key,"AR_width_x_y",AR_w);
  std::vector<double> SM_w; jmgr->GetVectorFromSubKey<double>(key,"SM_width_xp_xn_y",SM_w);
  std::vector<double> hcal_area = cut::hcal_active_area_data(0,0,2); 
  std::vector<double> hcal_AR = cut::hcal_active_area_data(AR_w[0],AR_w[1],2); 
  std::vector<double> hcal_SM = cut::hcal_safety_margin(SM_w[0],SM_w[1],SM_w[2],hcal_AR);  

  // output file
  std::string filebase = jmgr->GetValueFromSubKey_str(key,"outfile_prefix");
  filebase = filebase.empty() ? "" : filebase + "_";
  std::string outfilebase = "pdout/npR/" + dfprefix + "_npR_data_" + Form("sbs%d_sbs%dp_model%d_pass%d",conf,sbsmag,model,pass);
  //std::string outfilebase = "pdout/pDE/" + filebase + dfprefix + "_pDE_data_" + Form("sbs%d_sbsALL_model%d_pass%d",conf,model,pass);
  TFile *fout = new TFile(Form("%s.root",outfilebase.c_str()), "RECREATE");


  // Axes titles
  TString tdx = "#font[32]{#Deltay} (m)";
  TString tdy = "#font[32]{#Deltax} (m)";  
  TString txexp = "#font[32]{x^{exp}_{HCAL}} (m)";
  TString txexp_p = "#font[32]{x^{exp}_{HCAL} - #deltax_{SBS}} (m)";
  TString tyexp = "#font[32]{y^{exp}_{HCAL}} (m)";
  TString tnpr = "#font[32]{R^{QE}}";
  
  // call the canvas customizer
  PlotCustomizer pcust{1,1};
  PlotCustomizer pcustnostat;

#if 1
  
  // **************
  // visualizing n and p spots
  // **************
  TCanvas *cdxdy = util_pd::TC("cdxdy",1,2);
  //
  cdxdy->cd(1);
  gPad->SetLogz();
  gStyle->SetPalette(kRainbow);
  TH2F *h_dxdy_earm = (TH2F*)data_rdf.Filter(allearm.c_str()).Filter("eHCAL>0").Histo2D({"h_dxdy_earm","",200,-1.5,1.5,200,-3,2},"dy","dx")->Clone();
  util_pd::SetAxTitles(h_dxdy_earm,tdx,tdy);
  h_dxdy_earm->Draw("colz");
  //
  cdxdy->cd(2);
  gPad->SetLogz();
  gStyle->SetPalette(kRainbow);
  TH2F *h_dxdy_eNharm = (TH2F*)data_rdf.Filter(alleNharm.c_str()).Filter("eHCAL>0").Histo2D({"h_dxdy_eNharm","",200,-1.5,1.5,200,-3,2},"dy","dx")->Clone();
  util_pd::SetAxTitles(h_dxdy_eNharm,tdx,tdy);
  h_dxdy_eNharm->Draw("colz");
  // //
  // cdxdy->cd(3);
  // std::vector<double> h_w2_lim; jmgr->GetVectorFromSubKey<double>(key,"h_w2_lim",h_w2_lim);
  // TH1F *h_w2_earm = (TH1F*)data_rdf.Filter(noW2.c_str()).Histo1D({"h_w2_earm","",(int)h_w2_lim[0],h_w2_lim[1],h_w2_lim[2]},"W2")->Clone();
  // TH1F *h_w2_eNharm = (TH1F*)data_rdf.Filter(noW2.c_str()).Filter(harm_cut.c_str()).Histo1D({"h_w2_eNharm","",(int)h_w2_lim[0],h_w2_lim[1],h_w2_lim[2]},"W2")->Clone();
  // util_pd::SetAxTitles(h_w2_earm,"","#font[32]{W^{2}} (GeV^{2})");
  // custom_denom(h_w2_earm);
  // h_w2_earm->Draw("HIST");
  // custom_num(h_w2_eNharm);
  // h_w2_eNharm->Draw("HIST same");
  // TLegend *lw2 = new TLegend(0.15,0.75,0.55,0.9);
  // lw2->AddEntry(h_w2_earm,"All","lp");
  // lw2->AddEntry(h_w2_eNharm,"With #font[32]{#Deltax-#Deltay} Cut","lp");
  // lw2->SetFillStyle(0);
  // lw2->Draw();
  // //
  // cdxdy->cd(4);
  // std::string dx_w_def_cut = earm_cut + "&&eHCAL>0&&abs(dy)<0.3";
  // TH1F *h_dx_w_def = (TH1F*)data_rdf.Filter(dx_w_def_cut.c_str()).Histo1D({"h_dx_w_def","",200,-1,1},"dx_def")->Clone();
  // h_dx_w_def->Draw("HIST");
  // h_dx_w_def->GetXaxis()->SetTitle("#font[32]{#Deltax} + proton_deflection (m)");
  // std::vector<double> hdx_fitR{-0.5,0.5,1.2,1.2};
  // TF1 *fdx = fit::fit_1gs_nbg(hdx_fitR,h_dx_w_def);
  // double dxM = fdx->GetParameter(1), dxMerr = fdx->GetParError(1);
  // double dxS = fdx->GetParameter(2), dxSerr = fdx->GetParError(2);
  // pcust.AddTitleText(Form("#mu = %.4f #pm %.4f, #sigma = %.4f",dxM,dxMerr,dxS),0.95);
  // fdx->Draw("same");
  
  pcustnostat.customize_canvas(cdxdy);
  cdxdy->Update();
  cdxdy->Write();
  cdxdy->SaveAs(Form("%s.pdf[",outfilebase.c_str()));
  cdxdy->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cdxdy->SaveAs(Form("%s_dxdy.png",outfilebase.c_str()));
  //--
  
  // W2
  std::vector<double> h_W2_lim; jmgr->GetVectorFromSubKey<double>(key,"h_W2_lim",h_W2_lim);
  std::vector<double> h_W2_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_W2_fitR",h_W2_fitR);
  TH1F *h_W2_p = (TH1F*)prdf.Filter(noW2.c_str()).Histo1D({"h_W2_p","",int(h_W2_lim[0]),h_W2_lim[1],h_W2_lim[2]},"W2")->Clone();
  TH1F *h_W2_n = (TH1F*)nrdf.Filter(noW2.c_str()).Histo1D({"h_W2_n","",int(h_W2_lim[0]),h_W2_lim[1],h_W2_lim[2]},"W2")->Clone();
  TH1F *h_W2_npratio = new TH1F("h_W2_npratio","",int(h_W2_lim[0]),h_W2_lim[1],h_W2_lim[2]);
  h_W2_npratio->Divide(h_W2_n,h_W2_p);
  calc_binomial_error(h_W2_p,h_W2_npratio);
  custom_denom(h_W2_p);
  TString tW2 = "#font[32]{W^{2}} (GeV^{2})";
  util_pd::SetAxTitles(h_W2_p,"",tW2);
  custom_num(h_W2_n);
  util_pd::SetAxTitles(h_W2_n,"",tW2);
  custom_ratio(h_W2_npratio);
  util_pd::SetAxTitles(h_W2_npratio,tnpr,tW2);
  TCanvas *cW2 = util_pd::TC("cW2",1,2);
  cW2->cd(1); //
  gStyle->SetOptStat(0);
  h_W2_p->Draw();
  h_W2_n->Draw("same");
  cW2->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_W2_npratio->Draw();
  TF1 *fW2 = new TF1("fW2","pol0",h_W2_fitR[0],h_W2_fitR[1]);
  fW2->SetNpx(2000);
  h_W2_npratio->Fit("fW2","QR");
  cW2->Update();
  TPaveStats *stW2 = (TPaveStats*)h_W2_npratio->FindObject("stats");
  custom_statbox_effi(stW2);
  fW2->SetLineWidth(4);
  fW2->Draw("same");
  //
  pcust.customize_canvas(cW2);
  cW2->Update();
  cW2->Write();
  cW2->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cW2->SaveAs(Form("%s_W2.png",outfilebase.c_str()));  
  // --

  // dy
  std::vector<double> h_dy_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dy_lim",h_dy_lim);
  std::vector<double> h_dy_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dy_fitR",h_dy_fitR);
  TH1F *h_dy_p = (TH1F*)data_rdf.Filter(allearm.c_str()).Filter("pdx_nS<1").Histo1D({"h_dy_p","",int(h_dy_lim[0]),h_dy_lim[1],h_dy_lim[2]},"dy")->Clone();
  TH1F *h_dy_n = (TH1F*)data_rdf.Filter(allearm.c_str()).Filter("ndx_nS<1").Histo1D({"h_dy_n","",int(h_dy_lim[0]),h_dy_lim[1],h_dy_lim[2]},"dy")->Clone();
  TH1F *h_dy_npratio = new TH1F("h_dy_npratio","",int(h_dy_lim[0]),h_dy_lim[1],h_dy_lim[2]);
  h_dy_npratio->Divide(h_dy_n,h_dy_p);
  calc_binomial_error(h_dy_p,h_dy_npratio);
  custom_denom(h_dy_p);
  util_pd::SetAxTitles(h_dy_p,"",tdy);
  custom_num(h_dy_n);
  util_pd::SetAxTitles(h_dy_n,"",tdy);
  custom_ratio(h_dy_npratio);
  util_pd::SetAxTitles(h_dy_npratio,tnpr,tdy);
  TCanvas *cdy = util_pd::TC("cdy",1,2);
  cdy->cd(1); //
  gStyle->SetOptStat(0);
  h_dy_p->Draw();
  h_dy_n->Draw("same");
  cdy->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_dy_npratio->Draw();
  TF1 *fdy = new TF1("fdy","pol0",h_dy_fitR[0],h_dy_fitR[1]);
  fdy->SetNpx(2000);
  h_dy_npratio->Fit("fdy","QR");
  cdy->Update();
  TPaveStats *stdy = (TPaveStats*)h_dy_npratio->FindObject("stats");
  custom_statbox_effi(stdy);
  fdy->SetLineWidth(4);
  fdy->Draw("same");
  //
  pcust.customize_canvas(cdy);
  cdy->Update();
  cdy->Write();
  cdy->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cdy->SaveAs(Form("%s_dy.png",outfilebase.c_str()));  
  // --  

  // vz
  std::vector<double> h_vz_lim; jmgr->GetVectorFromSubKey<double>(key,"h_vz_lim",h_vz_lim);
  std::vector<double> h_vz_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_vz_fitR",h_vz_fitR);
  TH1F *h_vz_p = (TH1F*)prdf.Filter(novz.c_str()).Histo1D({"h_vz_p","",int(h_vz_lim[0]),h_vz_lim[1],h_vz_lim[2]},"vz")->Clone();
  TH1F *h_vz_n = (TH1F*)nrdf.Filter(novz.c_str()).Histo1D({"h_vz_n","",int(h_vz_lim[0]),h_vz_lim[1],h_vz_lim[2]},"vz")->Clone();
  TH1F *h_vz_npratio = new TH1F("h_vz_npratio","",int(h_vz_lim[0]),h_vz_lim[1],h_vz_lim[2]);
  h_vz_npratio->Divide(h_vz_n,h_vz_p);
  calc_binomial_error(h_vz_p,h_vz_npratio);
  custom_denom(h_vz_p);
  TString tvz = "#font[32]{v_{z}} (m)";
  util_pd::SetAxTitles(h_vz_p,"",tvz);
  custom_num(h_vz_n);
  util_pd::SetAxTitles(h_vz_n,"",tvz);
  custom_ratio(h_vz_npratio);
  util_pd::SetAxTitles(h_vz_npratio,tnpr,tvz);
  TCanvas *cvz = util_pd::TC("cvz",1,2);
  cvz->cd(1); //
  gStyle->SetOptStat(0);
  h_vz_p->Draw();
  h_vz_n->Draw("same");
  cvz->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_vz_npratio->Draw();
  TF1 *fvz = new TF1("fvz","pol0",h_vz_fitR[0],h_vz_fitR[1]);
  fvz->SetNpx(2000);
  h_vz_npratio->Fit("fvz","QR");
  cvz->Update();
  TPaveStats *stvz = (TPaveStats*)h_vz_npratio->FindObject("stats");
  custom_statbox_effi(stvz);
  fvz->SetLineWidth(4);
  fvz->Draw("same");
  //
  pcust.customize_canvas(cvz);
  cvz->Update();
  cvz->Write();
  if(print_canv_ppar) cvz->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_canv_ppar&&print_png) cvz->SaveAs(Form("%s_vz.png",outfilebase.c_str()));  
  // --

  // EovP
  std::vector<double> h_EovP_lim; jmgr->GetVectorFromSubKey<double>(key,"h_EovP_lim",h_EovP_lim);
  std::vector<double> h_EovP_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_EovP_fitR",h_EovP_fitR);
  TH1F *h_EovP_p = (TH1F*)prdf.Filter(noEovP.c_str()).Histo1D({"h_EovP_p","",int(h_EovP_lim[0]),h_EovP_lim[1],h_EovP_lim[2]},"EovP")->Clone();
  TH1F *h_EovP_n = (TH1F*)nrdf.Filter(noEovP.c_str()).Histo1D({"h_EovP_n","",int(h_EovP_lim[0]),h_EovP_lim[1],h_EovP_lim[2]},"EovP")->Clone();
  TH1F *h_EovP_npratio = new TH1F("h_EovP_npratio","",int(h_EovP_lim[0]),h_EovP_lim[1],h_EovP_lim[2]);
  h_EovP_npratio->Divide(h_EovP_n,h_EovP_p);
  calc_binomial_error(h_EovP_p,h_EovP_npratio);
  custom_denom(h_EovP_p);
  TString tEovP = "#font[32]{E_{BBCAL}/p}";
  util_pd::SetAxTitles(h_EovP_p,"",tEovP);
  custom_num(h_EovP_n);
  util_pd::SetAxTitles(h_EovP_n,"",tEovP);
  custom_ratio(h_EovP_npratio);
  util_pd::SetAxTitles(h_EovP_npratio,tnpr,tEovP);
  TCanvas *cEovP = util_pd::TC("cEovP",1,2);
  cEovP->cd(1); //
  gStyle->SetOptStat(0);
  h_EovP_p->Draw();
  h_EovP_n->Draw("same");
  cEovP->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_EovP_npratio->Draw();
  TF1 *fEovP = new TF1("fEovP","pol0",h_EovP_fitR[0],h_EovP_fitR[1]);
  fEovP->SetNpx(2000);
  h_EovP_npratio->Fit("fEovP","QR");
  cEovP->Update();
  TPaveStats *stEovP = (TPaveStats*)h_EovP_npratio->FindObject("stats");
  custom_statbox_effi(stEovP);
  fEovP->SetLineWidth(4);
  fEovP->Draw("same");
  //
  pcust.customize_canvas(cEovP);
  cEovP->Update();
  cEovP->Write();
  if(print_canv_ppar) cEovP->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_canv_ppar&&print_png) cEovP->SaveAs(Form("%s_EovP.png",outfilebase.c_str()));  
  // --

  // trchi2ndf
  std::vector<double> h_trchi2ndf_lim; jmgr->GetVectorFromSubKey<double>(key,"h_trchi2ndf_lim",h_trchi2ndf_lim);
  std::vector<double> h_trchi2ndf_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_trchi2ndf_fitR",h_trchi2ndf_fitR);
  TH1F *h_trchi2ndf_p = (TH1F*)prdf.Filter(notrchi2ndf.c_str()).Histo1D({"h_trchi2ndf_p","",int(h_trchi2ndf_lim[0]),h_trchi2ndf_lim[1],h_trchi2ndf_lim[2]},"trchi2ndf")->Clone();
  TH1F *h_trchi2ndf_n = (TH1F*)nrdf.Filter(notrchi2ndf.c_str()).Histo1D({"h_trchi2ndf_n","",int(h_trchi2ndf_lim[0]),h_trchi2ndf_lim[1],h_trchi2ndf_lim[2]},"trchi2ndf")->Clone();
  TH1F *h_trchi2ndf_npratio = new TH1F("h_trchi2ndf_npratio","",int(h_trchi2ndf_lim[0]),h_trchi2ndf_lim[1],h_trchi2ndf_lim[2]);
  h_trchi2ndf_npratio->Divide(h_trchi2ndf_n,h_trchi2ndf_p);
  calc_binomial_error(h_trchi2ndf_p,h_trchi2ndf_npratio);
  custom_denom(h_trchi2ndf_p);
  TString ttrchi2ndf = "Track #font[32]{#chi^{2}/NDF} (m)";
  util_pd::SetAxTitles(h_trchi2ndf_p,"",ttrchi2ndf);
  custom_num(h_trchi2ndf_n);
  util_pd::SetAxTitles(h_trchi2ndf_n,"",ttrchi2ndf);
  custom_ratio(h_trchi2ndf_npratio);
  util_pd::SetAxTitles(h_trchi2ndf_npratio,tnpr,ttrchi2ndf);
  TCanvas *ctrchi2ndf = util_pd::TC("ctrchi2ndf",1,2);
  ctrchi2ndf->cd(1); //
  gStyle->SetOptStat(0);
  h_trchi2ndf_p->Draw();
  h_trchi2ndf_n->Draw("same");
  ctrchi2ndf->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_trchi2ndf_npratio->Draw();
  TF1 *ftrchi2ndf = new TF1("ftrchi2ndf","pol0",h_trchi2ndf_fitR[0],h_trchi2ndf_fitR[1]);
  ftrchi2ndf->SetNpx(2000);
  h_trchi2ndf_npratio->Fit("ftrchi2ndf","QR");
  ctrchi2ndf->Update();
  TPaveStats *sttrchi2ndf = (TPaveStats*)h_trchi2ndf_npratio->FindObject("stats");
  custom_statbox_effi(sttrchi2ndf);
  ftrchi2ndf->SetLineWidth(4);
  ftrchi2ndf->Draw("same");
  //
  pcust.customize_canvas(ctrchi2ndf);
  ctrchi2ndf->Update();
  ctrchi2ndf->Write();
  if(print_canv_ppar) ctrchi2ndf->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_canv_ppar&&print_png) ctrchi2ndf->SaveAs(Form("%s_trchi2ndf.png",outfilebase.c_str()));  
  // --      

  // ePS
  std::vector<double> h_ePS_lim; jmgr->GetVectorFromSubKey<double>(key,"h_ePS_lim",h_ePS_lim);
  std::vector<double> h_ePS_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_ePS_fitR",h_ePS_fitR);
  TH1F *h_ePS_p = (TH1F*)prdf.Filter(noePS.c_str()).Histo1D({"h_ePS_p","",int(h_ePS_lim[0]),h_ePS_lim[1],h_ePS_lim[2]},"ePS")->Clone();
  TH1F *h_ePS_n = (TH1F*)nrdf.Filter(noePS.c_str()).Histo1D({"h_ePS_n","",int(h_ePS_lim[0]),h_ePS_lim[1],h_ePS_lim[2]},"ePS")->Clone();
  TH1F *h_ePS_npratio = new TH1F("h_ePS_npratio","",int(h_ePS_lim[0]),h_ePS_lim[1],h_ePS_lim[2]);
  h_ePS_npratio->Divide(h_ePS_n,h_ePS_p);
  calc_binomial_error(h_ePS_p,h_ePS_npratio);
  custom_denom(h_ePS_p);
  TString tePS = "#font[32]{E_{PS}} (GeV)";
  util_pd::SetAxTitles(h_ePS_p,"",tePS);
  custom_num(h_ePS_n);
  util_pd::SetAxTitles(h_ePS_n,"",tePS);
  custom_ratio(h_ePS_npratio);
  util_pd::SetAxTitles(h_ePS_npratio,tnpr,tePS);
  TCanvas *cePS = util_pd::TC("cePS",1,2);
  cePS->cd(1); //
  gStyle->SetOptStat(0);
  h_ePS_p->Draw();
  h_ePS_n->Draw("same");
  cePS->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_ePS_npratio->Draw();
  TF1 *fePS = new TF1("fePS","pol0",h_ePS_fitR[0],h_ePS_fitR[1]);
  fePS->SetNpx(2000);
  h_ePS_npratio->Fit("fePS","QR");
  cePS->Update();
  TPaveStats *stePS = (TPaveStats*)h_ePS_npratio->FindObject("stats");
  custom_statbox_effi(stePS);
  fePS->SetLineWidth(4);
  fePS->Draw("same");
  //
  pcust.customize_canvas(cePS);
  cePS->Update();
  cePS->Write();
  if(print_canv_ppar) cePS->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_canv_ppar&&print_png) cePS->SaveAs(Form("%s_ePS.png",outfilebase.c_str()));  
  // --

  // xbb
  std::vector<double> h_xbb_lim; jmgr->GetVectorFromSubKey<double>(key,"h_xbb_lim",h_xbb_lim);
  std::vector<double> h_xbb_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_xbb_fitR",h_xbb_fitR);
  TH1F *h_xbb_p = (TH1F*)prdf.Filter(noxbb.c_str()).Histo1D({"h_xbb_p","",int(h_xbb_lim[0]),h_xbb_lim[1],h_xbb_lim[2]},"xbb")->Clone();
  TH1F *h_xbb_n = (TH1F*)nrdf.Filter(noxbb.c_str()).Histo1D({"h_xbb_n","",int(h_xbb_lim[0]),h_xbb_lim[1],h_xbb_lim[2]},"xbb")->Clone();
  TH1F *h_xbb_npratio = new TH1F("h_xbb_npratio","",int(h_xbb_lim[0]),h_xbb_lim[1],h_xbb_lim[2]);
  h_xbb_npratio->Divide(h_xbb_n,h_xbb_p);
  calc_binomial_error(h_xbb_p,h_xbb_npratio);
  custom_denom(h_xbb_p);
  TString txbb = "#font[32]{x_{BB}} (m)";
  util_pd::SetAxTitles(h_xbb_p,"",txbb);
  custom_num(h_xbb_n);
  util_pd::SetAxTitles(h_xbb_n,"",txbb);
  custom_ratio(h_xbb_npratio);
  util_pd::SetAxTitles(h_xbb_npratio,tnpr,txbb);
  TCanvas *cxbb = util_pd::TC("cxbb",1,2);
  cxbb->cd(1); //
  gStyle->SetOptStat(0);
  h_xbb_p->Draw();
  h_xbb_n->Draw("same");
  cxbb->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_xbb_npratio->Draw();
  TF1 *fxbb = new TF1("fxbb","pol0",h_xbb_fitR[0],h_xbb_fitR[1]);
  fxbb->SetNpx(2000);
  h_xbb_npratio->Fit("fxbb","QR");
  cxbb->Update();
  TPaveStats *stxbb = (TPaveStats*)h_xbb_npratio->FindObject("stats");
  custom_statbox_effi(stxbb);
  fxbb->SetLineWidth(4);
  fxbb->Draw("same");
  //
  pcust.customize_canvas(cxbb);
  cxbb->Update();
  cxbb->Write();
  if(print_canv_ppar) cxbb->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_canv_ppar&&print_png) cxbb->SaveAs(Form("%s_xbb.png",outfilebase.c_str()));  
  // --

  // ybb
  std::vector<double> h_ybb_lim; jmgr->GetVectorFromSubKey<double>(key,"h_ybb_lim",h_ybb_lim);
  std::vector<double> h_ybb_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_ybb_fitR",h_ybb_fitR);
  TH1F *h_ybb_p = (TH1F*)prdf.Filter(noybb.c_str()).Histo1D({"h_ybb_p","",int(h_ybb_lim[0]),h_ybb_lim[1],h_ybb_lim[2]},"ybb")->Clone();
  TH1F *h_ybb_n = (TH1F*)nrdf.Filter(noybb.c_str()).Histo1D({"h_ybb_n","",int(h_ybb_lim[0]),h_ybb_lim[1],h_ybb_lim[2]},"ybb")->Clone();
  TH1F *h_ybb_npratio = new TH1F("h_ybb_npratio","",int(h_ybb_lim[0]),h_ybb_lim[1],h_ybb_lim[2]);
  h_ybb_npratio->Divide(h_ybb_n,h_ybb_p);
  calc_binomial_error(h_ybb_p,h_ybb_npratio);
  custom_denom(h_ybb_p);
  TString tybb = "#font[32]{y_{BB}} (m)";
  util_pd::SetAxTitles(h_ybb_p,"",tybb);
  custom_num(h_ybb_n);
  util_pd::SetAxTitles(h_ybb_n,"",tybb);
  custom_ratio(h_ybb_npratio);
  util_pd::SetAxTitles(h_ybb_npratio,tnpr,tybb);
  TCanvas *cybb = util_pd::TC("cybb",1,2);
  cybb->cd(1); //
  gStyle->SetOptStat(0);
  h_ybb_p->Draw();
  h_ybb_n->Draw("same");
  cybb->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_ybb_npratio->Draw();
  TF1 *fybb = new TF1("fybb","pol0",h_ybb_fitR[0],h_ybb_fitR[1]);
  fybb->SetNpx(2000);
  h_ybb_npratio->Fit("fybb","QR");
  cybb->Update();
  TPaveStats *stybb = (TPaveStats*)h_ybb_npratio->FindObject("stats");
  custom_statbox_effi(stybb);
  fybb->SetLineWidth(4);
  fybb->Draw("same");
  //
  pcust.customize_canvas(cybb);
  cybb->Update();
  cybb->Write();
  if(print_canv_ppar) cybb->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_canv_ppar&&print_png) cybb->SaveAs(Form("%s_ybb.png",outfilebase.c_str()));  
  // --
  
  // eHCAL
  std::vector<double> h_eHCAL_lim; jmgr->GetVectorFromSubKey<double>(key,"h_eHCAL_lim",h_eHCAL_lim);
  std::vector<double> h_eHCAL_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_eHCAL_fitR",h_eHCAL_fitR);
  TH1F *h_eHCAL_p = (TH1F*)prdf.Filter(noeHCAL.c_str()).Histo1D({"h_eHCAL_p","",int(h_eHCAL_lim[0]),h_eHCAL_lim[1],h_eHCAL_lim[2]},"eHCAL")->Clone();
  TH1F *h_eHCAL_n = (TH1F*)nrdf.Filter(noeHCAL.c_str()).Histo1D({"h_eHCAL_n","",int(h_eHCAL_lim[0]),h_eHCAL_lim[1],h_eHCAL_lim[2]},"eHCAL")->Clone();
  TH1F *h_eHCAL_npratio = new TH1F("h_eHCAL_npratio","",int(h_eHCAL_lim[0]),h_eHCAL_lim[1],h_eHCAL_lim[2]);
  h_eHCAL_npratio->Divide(h_eHCAL_n,h_eHCAL_p);
  calc_binomial_error(h_eHCAL_p,h_eHCAL_npratio);
  custom_denom(h_eHCAL_p);
  TString teHCAL = "#font[32]{E_{HCAL}} (GeV)";
  util_pd::SetAxTitles(h_eHCAL_p,"",teHCAL);
  custom_num(h_eHCAL_n);
  util_pd::SetAxTitles(h_eHCAL_n,"",teHCAL);
  custom_ratio(h_eHCAL_npratio);
  util_pd::SetAxTitles(h_eHCAL_npratio,tnpr,teHCAL);
  TCanvas *ceHCAL = util_pd::TC("ceHCAL",1,2);
  ceHCAL->cd(1); //
  gStyle->SetOptStat(0);
  h_eHCAL_p->Draw();
  h_eHCAL_n->Draw("same");
  ceHCAL->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_eHCAL_npratio->Draw();
  TF1 *feHCAL = new TF1("feHCAL","pol0",h_eHCAL_fitR[0],h_eHCAL_fitR[1]);
  feHCAL->SetNpx(2000);
  h_eHCAL_npratio->Fit("feHCAL","QR");
  ceHCAL->Update();
  TPaveStats *steHCAL = (TPaveStats*)h_eHCAL_npratio->FindObject("stats");
  custom_statbox_effi(steHCAL);
  feHCAL->SetLineWidth(4);
  feHCAL->Draw("same");
  //
  pcust.customize_canvas(ceHCAL);
  ceHCAL->Update();
  ceHCAL->Write();
  if(print_canv_ppar) ceHCAL->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_canv_ppar&&print_png) ceHCAL->SaveAs(Form("%s_eHCAL.png",outfilebase.c_str()));  
  // --

  // coinT
  std::vector<double> h_coinT_lim; jmgr->GetVectorFromSubKey<double>(key,"h_coinT_lim",h_coinT_lim);
  std::vector<double> h_coinT_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_coinT_fitR",h_coinT_fitR);
  TH1F *h_coinT_p = (TH1F*)prdf.Filter(nocoinT.c_str()).Histo1D({"h_coinT_p","",int(h_coinT_lim[0]),h_coinT_lim[1],h_coinT_lim[2]},coinTvar)->Clone();
  TH1F *h_coinT_n = (TH1F*)nrdf.Filter(nocoinT.c_str()).Histo1D({"h_coinT_n","",int(h_coinT_lim[0]),h_coinT_lim[1],h_coinT_lim[2]},coinTvar)->Clone();
  TH1F *h_coinT_npratio = new TH1F("h_coinT_npratio","",int(h_coinT_lim[0]),h_coinT_lim[1],h_coinT_lim[2]);
  h_coinT_npratio->Divide(h_coinT_n,h_coinT_p);
  //calc_binomial_error(h_coinT_p,h_coinT_npratio);
  custom_denom(h_coinT_p);
  TString tcoinT = "#font[32]{t^{ADC}_{coin}} (ns)";
  util_pd::SetAxTitles(h_coinT_p,"",tcoinT);
  custom_num(h_coinT_n);
  util_pd::SetAxTitles(h_coinT_n,"",tcoinT);
  custom_ratio(h_coinT_npratio);
  h_coinT_npratio->GetYaxis()->SetRangeUser(0,5);
  util_pd::SetAxTitles(h_coinT_npratio,tnpr,tcoinT);
  TCanvas *ccoinT = util_pd::TC("ccoinT",1,2);
  ccoinT->cd(1); //
  gStyle->SetOptStat(0);
  h_coinT_p->Draw();
  h_coinT_n->Draw("same");
  ccoinT->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_coinT_npratio->GetYaxis()->SetRangeUser(-0.2,5);
  h_coinT_npratio->Draw();
  TF1 *fcoinT = new TF1("fcoinT","pol0",h_coinT_fitR[0],h_coinT_fitR[1]);
  fcoinT->SetNpx(2000);
  h_coinT_npratio->Fit("fcoinT","QR");
  ccoinT->Update();
  TPaveStats *stcoinT = (TPaveStats*)h_coinT_npratio->FindObject("stats");
  custom_statbox_effi(stcoinT);
  fcoinT->SetLineWidth(4);
  fcoinT->Draw("same");
  //
  pcust.customize_canvas(ccoinT);
  ccoinT->Update();
  ccoinT->Write();
  if(print_canv_ppar) ccoinT->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_canv_ppar&&print_png) ccoinT->SaveAs(Form("%s_coinT.png",outfilebase.c_str()));  
  // --
  
  // **********
  // All Cuts
  TCanvas *callcut = util_pd::TC("callcut",3,3);
  callcut->cd(1); //
  h_W2_npratio->Draw();
  fW2->Draw("same");
  callcut->cd(2); //
  h_vz_npratio->Draw();
  fvz->Draw("same");
  callcut->cd(3); //
  h_EovP_npratio->Draw();
  fEovP->Draw("same");
  callcut->cd(4); //
  h_trchi2ndf_npratio->Draw();
  ftrchi2ndf->Draw("same");
  callcut->cd(5); //
  h_ePS_npratio->Draw();
  fePS->Draw("same");
  callcut->cd(6); //
  h_xbb_npratio->Draw();
  fxbb->Draw("same");
  callcut->cd(7); //
  h_ybb_npratio->Draw();
  fybb->Draw("same");
  callcut->cd(8); //
  h_eHCAL_npratio->Draw();
  feHCAL->Draw("same");
  callcut->cd(9); //
  h_coinT_npratio->Draw();
  fcoinT->Draw("same");    
  //
  pcust.customize_canvas(callcut);
  callcut->Update();
  callcut->Write();
  callcut->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) callcut->SaveAs(Form("%s_allcut.png",outfilebase.c_str()));  
  // --
  
  // **************
  // earm cuts for thesis
  TCanvas *cearm = new TCanvas("cearm", "cearm", 1350, 1200);
  cearm->Divide(2,3);
  // TCanvas *cearm = util_pd::TC("cearm",3,4);
  cearm->cd(1); //
  h_vz_p->Draw();
  h_vz_n->Draw("same");
  drawcutrange(h_vz_fitR[0],h_vz_fitR[1]);
  cearm->cd(2); //
  h_vz_npratio->Draw();
  fvz->Draw("same");
  fvz->SetLineColor(kMagenta);
  drawcutrange(h_vz_fitR[0],h_vz_fitR[1],1);
  cearm->cd(3); //
  h_xbb_p->Draw();
  h_xbb_n->Draw("same");
  drawcutrange(h_xbb_fitR[0],h_xbb_fitR[1]);
  cearm->cd(4); //
  h_xbb_npratio->Draw();
  fxbb->Draw("same");
  fxbb->SetLineColor(kMagenta);
  drawcutrange(h_xbb_fitR[0],h_xbb_fitR[1],1);
  cearm->cd(5); //
  h_ybb_p->Draw();
  h_ybb_n->Draw("same");
  drawcutrange(h_ybb_fitR[0],h_ybb_fitR[1]);
  cearm->cd(6); //
  h_ybb_npratio->Draw();
  fybb->Draw("same");
  fybb->SetLineColor(kMagenta);
  drawcutrange(h_ybb_fitR[0],h_ybb_fitR[1],1);  
  //
  pcust.customize_canvas(cearm);
  custom_extra(h_vz_p);
  custom_extra(h_vz_npratio);
  custom_extra(h_xbb_p);
  custom_extra(h_xbb_npratio);
  custom_extra(h_ybb_p);
  custom_extra(h_ybb_npratio);
  custom_extra(h_EovP_p);
  custom_extra(h_EovP_npratio);
  
  
  cearm->Update();
  cearm->Write();
  cearm->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cearm->SaveAs(Form("%s_earm.pdf",outfilebase.c_str()));  
  // --

  // earm1 cuts for thesis
  TCanvas *cearm1 = new TCanvas("cearm1", "cearm1", 1350, 1200);
  cearm1->Divide(2,3);
  // TCanvas *cearm1 = util_pd::TC("cearm1",3,4);
  cearm1->cd(1); //
  h_trchi2ndf_p->Draw();
  h_trchi2ndf_n->Draw("same");
  drawcutrange(h_trchi2ndf_fitR[1],0,0,1);
  cearm1->cd(2); //
  h_trchi2ndf_npratio->Draw();
  ftrchi2ndf->Draw("same");
  ftrchi2ndf->SetLineColor(kMagenta);
  drawcutrange(h_trchi2ndf_fitR[1],0,1,1);
  cearm1->cd(3); //
  h_ePS_p->Draw();
  h_ePS_n->Draw("same");
  drawcutrange(h_ePS_fitR[0],h_ePS_fitR[1],0,1);
  cearm1->cd(4); //
  h_ePS_npratio->Draw();
  fePS->Draw("same");
  fePS->SetLineColor(kMagenta);
  drawcutrange(h_ePS_fitR[0],h_ePS_fitR[1],1,1);
  cearm1->cd(5); //
  h_eHCAL_p->Draw();
  h_eHCAL_n->Draw("same");
  drawcutrange(h_eHCAL_fitR[0],h_eHCAL_fitR[1],0,1);
  cearm1->cd(6); //
  h_eHCAL_npratio->Draw();
  feHCAL->Draw("same");
  feHCAL->SetLineColor(kMagenta);
  drawcutrange(h_eHCAL_fitR[0],h_eHCAL_fitR[1],1,1);
  //
  pcust.customize_canvas(cearm1);
  custom_extra(h_trchi2ndf_p);
  custom_extra(h_trchi2ndf_npratio);
  custom_extra(h_ePS_p);
  custom_extra(h_ePS_npratio);
  
  cearm1->Update();
  cearm1->Write();
  cearm1->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  //cearm1->SaveAs(Form("%s_earm1.pdf",outfilebase.c_str()));  
  // --

  // cearm2 cuts for thesis
  TCanvas *cearm2 = new TCanvas("cearm2", "cearm2", 1350, 800);
  cearm2->Divide(2,2);
  // TCanvas *cearm2 = util_pd::TC("cearm2",3,4);
  cearm2->cd(1); //
  h_coinT_p->Draw();
  h_coinT_n->Draw("same");
  drawcutrange(h_coinT_fitR[0],h_coinT_fitR[1]);
  cearm2->cd(2); //
  h_coinT_npratio->Draw();
  fcoinT->Draw("same");
  fcoinT->SetLineColor(kMagenta);
  drawcutrange(h_coinT_fitR[0],h_coinT_fitR[1]);
  cearm2->cd(3); //
  h_EovP_p->Draw();
  h_EovP_n->Draw("same");
  drawcutrange(h_EovP_fitR[0],h_EovP_fitR[1]);
  cearm2->cd(4); //
  h_EovP_npratio->Draw();
  fEovP->Draw("same");
  fEovP->SetLineColor(kMagenta);
  drawcutrange(h_EovP_fitR[0],h_EovP_fitR[1],1);  
  //
  pcust.customize_canvas(cearm2);
  custom_extra(h_coinT_p);
  custom_extra(h_coinT_npratio);
  custom_extra(h_eHCAL_p);
  custom_extra(h_eHCAL_npratio);   
  
  cearm2->Update();
  cearm2->Write();
  cearm2->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  //cearm2->SaveAs(Form("%s_earm2.pdf",outfilebase.c_str()));  
  // --  
  
  // // cearm2 cuts for thesis
  // TCanvas *cearm2 = new TCanvas("cearm2", "cearm2", 1350, 1200);
  // cearm2->Divide(2,3);
  // // TCanvas *cearm2 = util_pd::TC("cearm2",3,4);
  // cearm2->cd(1); //
  // h_coinT_p->Draw();
  // h_coinT_n->Draw("same");
  // drawcutrange(-5.1,5.1);
  // cearm2->cd(2); //
  // h_coinT_npratio->Draw();
  // // fcoinT->Draw("same");
  // // fcoinT->SetLineColor(kMagenta);
  // drawcutrange(-5.1,5.1);
  // cearm2->cd(3); //
  // h_trchi2ndf_p->Draw();
  // h_trchi2ndf_n->Draw("same");
  // drawcutrange(h_trchi2ndf_fitR[1],0,0,1);
  // cearm2->cd(4); //
  // h_trchi2ndf_npratio->Draw();
  // ftrchi2ndf->Draw("same");
  // ftrchi2ndf->SetLineColor(kMagenta);
  // drawcutrange(h_trchi2ndf_fitR[1],0,1,1);
  // // cearm2->cd(5); //
  // // h_ybb_p->Draw();
  // // h_ybb_n->Draw("same");
  // // drawcutrange(h_ybb_fitR[0],h_ybb_fitR[1]);
  // // cearm2->cd(6); //
  // // h_ybb_npratio->Draw();
  // // fybb->Draw("same");
  // // fybb->SetLineColor(kMagenta);
  // // drawcutrange(h_ybb_fitR[0],h_ybb_fitR[1],1);
  // //
  // pcust.customize_canvas(cearm2);
  // custom_extra(h_coinT_p);
  // custom_extra(h_coinT_npratio);
  // custom_extra(h_trchi2ndf_p);
  // custom_extra(h_trchi2ndf_npratio);
  // // custom_extra(h_ybb_p);
  // // custom_extra(h_ybb_npratio); 
  
  // cearm2->Update();
  // cearm2->Write();
  // cearm2->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  // cearm2->SaveAs(Form("%s_earm2.png",outfilebase.c_str()));  
  // // --
  
  
  // ******************************************************************************************
  //
  // **************
  // HCAL Related 
  // **************
#endif
  
  // effi in xexp
  // p
  std::vector<double> h_xefp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_xefp_lim",h_xefp_lim);
  std::vector<double> h_xefp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_xefp_fitR",h_xefp_fitR);
  TH1F *h_xefp_de = (TH1F*)data_rdf.Filter(allearmNOsmx.c_str()).Histo1D({"h_xefp_de","",int(h_xefp_lim[0]),h_xefp_lim[1],h_xefp_lim[2]},"xHCAL_exp_def")->Clone();
  TH1F *h_xefp_nu = (TH1F*)data_rdf.Filter(allearmNOsmx.c_str()).Filter(onlyp.c_str()).Histo1D({"h_xefp_nu","",int(h_xefp_lim[0]),h_xefp_lim[1],h_xefp_lim[2]},"xHCAL_exp_def")->Clone();
  TH1F *h_xefp_ef = new TH1F("h_xefp_ef","",int(h_xefp_lim[0]),h_xefp_lim[1],h_xefp_lim[2]);
  h_xefp_ef->Divide(h_xefp_nu,h_xefp_de);
  calc_binomial_error(h_xefp_de,h_xefp_ef);
  custom_denom(h_xefp_de);
  util_pd::SetAxTitles(h_xefp_de,"",txexp_p);
  custom_num(h_xefp_nu);
  util_pd::SetAxTitles(h_xefp_nu,"",txexp_p);
  custom_ratio(h_xefp_ef);
  util_pd::SetAxTitles(h_xefp_ef,"Proton Count/No. of Expected Nucleon",txexp_p);
  // n
  std::vector<double> h_xefn_lim; jmgr->GetVectorFromSubKey<double>(key,"h_xefn_lim",h_xefn_lim);
  std::vector<double> h_xefn_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_xefn_fitR",h_xefn_fitR);
  TH1F *h_xefn_de = (TH1F*)data_rdf.Filter(allearmNOsmx.c_str()).Histo1D({"h_xefn_de","",int(h_xefn_lim[0]),h_xefn_lim[1],h_xefn_lim[2]},"xHCAL_exp")->Clone();
  TH1F *h_xefn_nu = (TH1F*)data_rdf.Filter(allearmNOsmx.c_str()).Filter(onlyn.c_str()).Histo1D({"h_xefn_nu","",int(h_xefn_lim[0]),h_xefn_lim[1],h_xefn_lim[2]},"xHCAL_exp")->Clone();
  TH1F *h_xefn_ef = new TH1F("h_xefn_ef","",int(h_xefn_lim[0]),h_xefn_lim[1],h_xefn_lim[2]);
  h_xefn_ef->Divide(h_xefn_nu,h_xefn_de);
  calc_binomial_error(h_xefn_de,h_xefn_ef);
  custom_denom(h_xefn_de);
  util_pd::SetAxTitles(h_xefn_de,"",txexp);
  custom_num(h_xefn_nu);
  util_pd::SetAxTitles(h_xefn_nu,"",txexp);
  custom_ratio(h_xefn_ef);
  util_pd::SetAxTitles(h_xefn_ef,"Neutron Count/No. of Expected Nucleon",txexp);
  //
  TCanvas *cxef = util_pd::TC("cxef",2,2);
  cxef->cd(1); //
  gStyle->SetOptStat(0);
  h_xefp_de->Draw();
  h_xefp_nu->Draw("same");
  cxef->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_xefp_ef->GetYaxis()->SetRangeUser(0,0.6);
  h_xefp_ef->Draw();
  TF1 *fxefp = new TF1("fxefp","pol0",h_xefp_fitR[0],h_xefp_fitR[1]);
  fxefp->SetNpx(2000);
  h_xefp_ef->Fit("fxefp","QR");
  cxef->Update();
  TPaveStats *stxefp = (TPaveStats*)h_xefp_ef->FindObject("stats");
  custom_statbox_effi(stxefp);
  fxefp->SetLineWidth(4);
  fxefp->Draw("same");
  cxef->cd(3); //
  gStyle->SetOptStat(0);
  h_xefn_de->Draw();
  h_xefn_nu->Draw("same");
  cxef->cd(4); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_xefn_ef->GetYaxis()->SetRangeUser(0,0.6);
  h_xefn_ef->Draw();
  TF1 *fxefn = new TF1("fxefn","pol0",h_xefn_fitR[0],h_xefn_fitR[1]);
  fxefn->SetNpx(2000);
  h_xefn_ef->Fit("fxefn","QR");
  cxef->Update();
  TPaveStats *stxefn = (TPaveStats*)h_xefn_ef->FindObject("stats");
  custom_statbox_effi(stxefn);
  fxefn->SetLineWidth(4);
  fxefn->Draw("same");  
  //
  pcust.customize_canvas(cxef);
  cxef->Update();
  cxef->Write();
  cxef->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cxef->SaveAs(Form("%s_xefp.png",outfilebase.c_str()));  
  // --

  // effi in yexp
  // p
  std::vector<double> h_yefp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_yefp_lim",h_yefp_lim);
  std::vector<double> h_yefp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_yefp_fitR",h_yefp_fitR);
  TH1F *h_yefp_de = (TH1F*)data_rdf.Filter(allearmNOsmx.c_str()).Histo1D({"h_yefp_de","",int(h_yefp_lim[0]),h_yefp_lim[1],h_yefp_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yefp_nu = (TH1F*)data_rdf.Filter(allearmNOsmx.c_str()).Filter(onlyp.c_str()).Histo1D({"h_yefp_nu","",int(h_yefp_lim[0]),h_yefp_lim[1],h_yefp_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yefp_ef = new TH1F("h_yefp_ef","",int(h_yefp_lim[0]),h_yefp_lim[1],h_yefp_lim[2]);
  h_yefp_ef->Divide(h_yefp_nu,h_yefp_de);
  calc_binomial_error(h_yefp_de,h_yefp_ef);
  custom_denom(h_yefp_de);
  util_pd::SetAxTitles(h_yefp_de,"",tyexp);
  custom_num(h_yefp_nu);
  util_pd::SetAxTitles(h_yefp_nu,"",tyexp);
  custom_ratio(h_yefp_ef);
  util_pd::SetAxTitles(h_yefp_ef,"Proton Count/No. of Expected Nucleon",tyexp);
  // n
  std::vector<double> h_yefn_lim; jmgr->GetVectorFromSubKey<double>(key,"h_yefn_lim",h_yefn_lim);
  std::vector<double> h_yefn_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_yefn_fitR",h_yefn_fitR);
  TH1F *h_yefn_de = (TH1F*)data_rdf.Filter(allearmNOsmx.c_str()).Histo1D({"h_yefn_de","",int(h_yefn_lim[0]),h_yefn_lim[1],h_yefn_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yefn_nu = (TH1F*)data_rdf.Filter(allearmNOsmx.c_str()).Filter(onlyn.c_str()).Histo1D({"h_yefn_nu","",int(h_yefn_lim[0]),h_yefn_lim[1],h_yefn_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yefn_ef = new TH1F("h_yefn_ef","",int(h_yefn_lim[0]),h_yefn_lim[1],h_yefn_lim[2]);
  h_yefn_ef->Divide(h_yefn_nu,h_yefn_de);
  calc_binomial_error(h_yefn_de,h_yefn_ef);
  custom_denom(h_yefn_de);
  util_pd::SetAxTitles(h_yefn_de,"",tyexp);
  custom_num(h_yefn_nu);
  util_pd::SetAxTitles(h_yefn_nu,"",tyexp);
  custom_ratio(h_yefn_ef);
  util_pd::SetAxTitles(h_yefn_ef,"Neutron Count/No. of Expected Nucleon",tyexp);
  //
  TCanvas *cyef = util_pd::TC("cyef",2,2);
  cyef->cd(1); //
  gStyle->SetOptStat(0);
  h_yefp_de->Draw();
  h_yefp_nu->Draw("same");
  cyef->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_yefp_ef->GetYaxis()->SetRangeUser(0,0.6);
  h_yefp_ef->Draw();
  TF1 *fyefp = new TF1("fyefp","pol0",h_yefp_fitR[0],h_yefp_fitR[1]);
  fyefp->SetNpx(2000);
  h_yefp_ef->Fit("fyefp","QR");
  cyef->Update();
  TPaveStats *styefp = (TPaveStats*)h_yefp_ef->FindObject("stats");
  custom_statbox_effi(styefp);
  fyefp->SetLineWidth(4);
  fyefp->Draw("same");
  cyef->cd(3); //
  gStyle->SetOptStat(0);
  h_yefn_de->Draw();
  h_yefn_nu->Draw("same");
  cyef->cd(4); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_yefn_ef->GetYaxis()->SetRangeUser(0,0.6);
  h_yefn_ef->Draw();
  TF1 *fyefn = new TF1("fyefn","pol0",h_yefn_fitR[0],h_yefn_fitR[1]);
  fyefn->SetNpx(2000);
  h_yefn_ef->Fit("fyefn","QR");
  cyef->Update();
  TPaveStats *styefn = (TPaveStats*)h_yefn_ef->FindObject("stats");
  custom_statbox_effi(styefn);
  fyefn->SetLineWidth(4);
  fyefn->Draw("same");  
  //
  pcust.customize_canvas(cyef);
  cyef->Update();
  cyef->Write();
  cyef->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cyef->SaveAs(Form("%s_yefp.png",outfilebase.c_str()));  
  // --  
  
  // npratio vs. xexp and yexp 
  std::vector<double> h_xexp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_xexp_lim",h_xexp_lim);
  std::vector<double> h_xexp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_xexp_fitR",h_xexp_fitR);
  TH1F *h_xexp_p = (TH1F*)prdf.Filter(nosmx.c_str()).Histo1D({"h_xexp_p","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xHCAL_exp")->Clone();
  TH1F *h_xexp_n = (TH1F*)nrdf.Filter(nosmx.c_str()).Histo1D({"h_xexp_n","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]},"xHCAL_exp")->Clone();
  TH1F *h_xexp_npratio = new TH1F("h_xexp_npratio","",int(h_xexp_lim[0]),h_xexp_lim[1],h_xexp_lim[2]);
  h_xexp_npratio->Divide(h_xexp_n,h_xexp_p);
  calc_binomial_error(h_xexp_p,h_xexp_npratio);
  custom_denom(h_xexp_p);
  util_pd::SetAxTitles(h_xexp_p,"",txexp);
  custom_num(h_xexp_n);
  util_pd::SetAxTitles(h_xexp_n,"",txexp);
  custom_ratio(h_xexp_npratio);
  util_pd::SetAxTitles(h_xexp_npratio,tnpr,txexp);
  // y
  std::vector<double> h_yexp_lim; jmgr->GetVectorFromSubKey<double>(key,"h_yexp_lim",h_yexp_lim);
  std::vector<double> h_yexp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_yexp_fitR",h_yexp_fitR);
  TH1F *h_yexp_p = (TH1F*)prdf.Filter(nosmy.c_str()).Histo1D({"h_yexp_p","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yexp_n = (TH1F*)nrdf.Filter(nosmy.c_str()).Histo1D({"h_yexp_n","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]},"yHCAL_exp")->Clone();
  TH1F *h_yexp_npratio = new TH1F("h_yexp_npratio","",int(h_yexp_lim[0]),h_yexp_lim[1],h_yexp_lim[2]);
  h_yexp_npratio->Divide(h_yexp_n,h_yexp_p);
  calc_binomial_error(h_yexp_p,h_yexp_npratio);
  custom_denom(h_yexp_p);
  util_pd::SetAxTitles(h_yexp_p,"",tyexp);
  custom_num(h_yexp_n);
  util_pd::SetAxTitles(h_yexp_n,"",tyexp);
  custom_ratio(h_yexp_npratio);
  util_pd::SetAxTitles(h_yexp_npratio,tnpr,tyexp);
  //
  TCanvas *cxyexp = util_pd::TC("cxyexp",2,2);
  cxyexp->cd(1); //
  gStyle->SetOptStat(0);
  h_xexp_p->Draw();
  h_xexp_n->Draw("same");
  cxyexp->cd(2); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_xexp_npratio->GetYaxis()->SetRangeUser(0,0.8);
  h_xexp_npratio->Draw();
  TF1 *fxexp = new TF1("fxexp","pol0",h_xexp_fitR[0],h_xexp_fitR[1]);
  fxexp->SetNpx(2000);
  h_xexp_npratio->Fit("fxexp","QR");
  cxyexp->Update();
  TPaveStats *stxexp = (TPaveStats*)h_xexp_npratio->FindObject("stats");
  custom_statbox_effi(stxexp);
  fxexp->SetLineWidth(4);
  fxexp->Draw("same");
  cxyexp->cd(3); //
  gStyle->SetOptStat(0);
  h_yexp_p->Draw();
  h_yexp_n->Draw("same");
  cxyexp->cd(4); //
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(1);
  h_yexp_npratio->Draw();
  TF1 *fyexp = new TF1("fyexp","pol0",h_yexp_fitR[0],h_yexp_fitR[1]);
  fyexp->SetNpx(2000);
  h_yexp_npratio->Fit("fyexp","QR");
  cxyexp->Update();
  TPaveStats *styexp = (TPaveStats*)h_yexp_npratio->FindObject("stats");
  custom_statbox_effi(styexp);
  fyexp->SetLineWidth(4);
  fyexp->Draw("same");  
  //
  pcust.customize_canvas(cxyexp);
  cxyexp->Update();
  cxyexp->Write();
  cxyexp->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cxyexp->SaveAs(Form("%s_xyexp.png",outfilebase.c_str()));  
  // --

  // **************
  // Envelopes
  // **************
  // std::vector<double> hcal_area = cut::hcal_active_area_data(0,0,pass);
  TH2F *h2_xyexp_all_n = (TH2F*)data_rdf.Filter(allearmNOsm.c_str()).Histo2D({"h2_xyexp_all_n","All (Neutron)",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xHCAL_exp")->Clone();
  TH2F *h2_xyexp_all_p = (TH2F*)data_rdf.Filter(allearmNOsm.c_str()).Histo2D({"h2_xyexp_all_p","All (Proton)",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xHCAL_exp_def")->Clone();
  TH2F *h2_xyexp_pass_n = (TH2F*)data_rdf.Filter(eNharm.c_str()).Histo2D({"h2_xyexp_pass","Passed HCAL (Neutron)",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xHCAL_exp")->Clone();
  TH2F *h2_xyexp_pass_p = (TH2F*)data_rdf.Filter(eNharm.c_str()).Histo2D({"h2_xyexp_pass","Passed HCAL (Proton)",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xHCAL_exp_def")->Clone();
  TH2F *h2_xyexp_fail_n = (TH2F*)data_rdf.Filter(eNantiharm.c_str()).Histo2D({"h2_xyexp_pass","Failed HCAL (Neutron)",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xHCAL_exp")->Clone();
  TH2F *h2_xyexp_fail_p = (TH2F*)data_rdf.Filter(eNantiharm.c_str()).Histo2D({"h2_xyexp_fail","Failed HCAL (Proton)",200,-1.25,1.25,200,-3.25,2.5},"yHCAL_exp","xHCAL_exp_def")->Clone();
  TCanvas *cenv = util_pd::TC("cenv",3,2);
  gStyle->SetOptStat(0);
  cenv->cd(1); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_all_n->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->cd(2); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_all_p->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->cd(3); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_pass_n->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->cd(4); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_pass_p->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);  
  cenv->cd(5); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_fail_n->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->cd(6); //
  gPad->SetLogz();
  gPad->SetGridx();
  gPad->SetGridy();
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_fail_p->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  cenv->Write();
  cenv->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cenv->SaveAs(Form("%s_env.png",outfilebase.c_str()));
  //--
  
  // Efficiency Map
  int nBinY = 32; //65; <- pass2
  double hMinY = -1.25;
  double hMaxY = 1.25;
  int nBinX = 63; //126; <- pass2
  double hMinX = -3.25;
  double hMaxX = 1.75;    
  TH2F *h2_xyexp_earm_n = (TH2F*)data_rdf.Filter(allearmNOsm.c_str()).Histo2D({"h2_xyexp_earm_n","",nBinY,hMinY,hMaxY,nBinX,hMinX,hMaxX},"yHCAL_exp","xHCAL_exp")->Clone();
  TH2F *h2_xyexp_eNharm_n = (TH2F*)data_rdf.Filter(eNharm.c_str()).Histo2D({"h2_xyexp_eNharm_n","",nBinY,hMinY,hMaxY,nBinX,hMinX,hMaxX},"yHCAL_exp","xHCAL_exp")->Clone();
  TH2F *h2_effi_map_n = new TH2F("h2_effi_map_n","",nBinY,hMinY,hMaxY,nBinX,hMinX,hMaxX);
  h2_effi_map_n->Divide(h2_xyexp_eNharm_n,h2_xyexp_earm_n);
  h2_effi_map_n->GetXaxis()->SetTitle(tyexp.Data());
  h2_effi_map_n->GetYaxis()->SetTitle(txexp.Data());  
  //
  TH2F *h2_xyexp_earm = (TH2F*)data_rdf.Filter(allearmNOsm.c_str()).Histo2D({"h2_xyexp_earm","",nBinY,hMinY,hMaxY,nBinX,hMinX,hMaxX},"yHCAL_exp","xHCAL_exp_def")->Clone();
  TH2F *h2_xyexp_eNharm = (TH2F*)data_rdf.Filter(eNharm.c_str()).Histo2D({"h2_xyexp_eNharm","",nBinY,hMinY,hMaxY,nBinX,hMinX,hMaxX},"yHCAL_exp","xHCAL_exp_def")->Clone();
  TH2F *h2_effi_map = new TH2F("h2_effi_map","",nBinY,hMinY,hMaxY,nBinX,hMinX,hMaxX);
  h2_effi_map->Divide(h2_xyexp_eNharm,h2_xyexp_earm);
  h2_effi_map->GetXaxis()->SetTitle(tyexp.Data());
  h2_effi_map->GetYaxis()->SetTitle(txexp_p.Data());
  // // Binomial error
  // for( int i=1; i<=h2_effi_map->GetNbinsX(); i++ ){
  //   for( int j=1; j<=h2_effi_map->GetNbinsY(); j++ ){
  //     int bin = h2_effi_map->GetBin(i,j);
  //     double effi = h2_effi_map->GetBinContent(bin);
  //     double N = std::max(1.0,h2_xyexp_earm->GetBinContent(bin));
  //     h2_effi_map->SetBinError(bin,sqrt(effi*(1.0-effi)/N));
  //   }
  // }
  TCanvas *cefmap = util_pd::TC("cefmap",1,2);
  gStyle->SetOptStat(0);
  cefmap->cd(1); //
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_effi_map_n->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  util_pd::DrawArea(hcal_SM,4,4,9);
  pcust.AddTitleText("Neutron Envelope",0.55);
  cefmap->cd(2); //
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_effi_map->Draw("colz");
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  util_pd::DrawArea(hcal_SM,4,4,9);
  pcust.AddTitleText("Proton Envelope",0.52);
  //
  pcust.customize_canvas(cefmap);
  cefmap->Update();
  cefmap->Write();
  cefmap->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  if(print_png) cefmap->SaveAs(Form("%s_efmap.png",outfilebase.c_str()));
  //--

  TCanvas *cenvth = util_pd::TC("cenvth",1,2);
  //cenvth->Divide(2,1);
  gStyle->SetOptStat(0);
  cenvth->cd(1); //
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_all_n->Draw("colz");
  h2_xyexp_all_n->GetXaxis()->SetTitle(tyexp.Data());
  h2_xyexp_all_n->GetYaxis()->SetTitle(txexp.Data());  
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  util_pd::DrawArea(hcal_SM,4,4,9);
  pcust.AddTitleText("Neutron Envelope",0.56);
  cenvth->cd(2); //
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);
  h2_xyexp_all_p->Draw("colz");
  h2_xyexp_all_p->GetXaxis()->SetTitle(tyexp.Data());
  h2_xyexp_all_p->GetYaxis()->SetTitle(txexp_p.Data());  
  util_pd::DrawArea(hcal_area,kGreen+2,2,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  util_pd::DrawArea(hcal_SM,4,4,9);
  pcust.AddTitleText("Proton Envelope",0.53);
  //
  pcust.customize_canvas(cenvth);
  cenvth->Update();
  cenvth->Write();
  cenvth->SaveAs(Form("%s.pdf",outfilebase.c_str()));
  cenvth->SaveAs(Form("%s.pdf]",outfilebase.c_str()));
  //cenvth->SaveAs(Form("%s_envth.pdf",outfilebase.c_str()));
  //--  
  
  fout->Write();

  std::cout << "------" << std::endl;
  std::cout << " Summary plots  : " << Form("%s.pdf",outfilebase.c_str()) << std::endl;
  std::cout << " Output ROOT file  : " << Form("%s.root",outfilebase.c_str()) << std::endl;
  std::cout << "------" << std::endl << std::endl;  
  
  return 0;
}
