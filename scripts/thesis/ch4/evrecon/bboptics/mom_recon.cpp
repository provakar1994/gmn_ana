/*
  Script to generate momentum reconstruction results plot in chapter 4 for thesis
  The following two types of plots are made:
  1. pethbend vs xptg 
  2. dpel 
  The choice of the above two is controlled by a flag named draw_choice. 
*/

#include "TF1.h"

#include "gmn_ana.h"

//______________________________________________________________________________
void GetDpelHistos(TString const & infile, TString const & elascuts, std::vector<TH1F*> &histos) {
  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame rdf("Tout",infile);

  int nbins = 200;
  double hlow = -0.15;
  double hhi = 0.15;

  char const * xtitle = "p/p_{elas} - 1";

  // create and customize inclusive dpel
  TH1F *h_all = (TH1F*)rdf.Histo1D({"h_all","",nbins,hlow,hhi},"dpel")->Clone();
  h_all->GetXaxis()->SetTitle(xtitle);
  h_all->SetLineWidth(2);
  h_all->SetLineColor(kGray+2);
  h_all->SetFillColor(kGray);
  
  // create and customize elastic dpel
  TH1F *h_ep = (TH1F*)rdf.Filter(elascuts.Data()).Histo1D({"h_ep","",nbins,hlow,hhi},"dpel")->Clone();
  h_ep->GetXaxis()->SetTitle(xtitle);
  h_ep->SetLineWidth(2);
  h_ep->SetLineColor(kGray+3);
  h_ep->SetFillColor(kGray+3);

  // filling the output vector
  histos = {h_all,h_ep};
}

//______________________________________________________________________________
void GetpThxpTgHisto(TString const & infile, TString const & elascuts, std::vector<TGraph*> &obj, int nevents, bool isData) {
  // reading ROOT files as df
  //ROOT::EnableImplicitMT();
  ROOT::RDataFrame rdf("Tout",infile);
  char const * pethbend = isData ? "pelas*ethbend" : "pelas*ethbend";
  //char const * pethbend = isData ? "trP_corr*ethbend" : "trP*ethbend";
  auto rdf_m = rdf.Define("pethbend",pethbend);

  int nxbins = 200;
  double hxlow = -0.3;
  double hxhi = 0.3;
  int nybins = 200;
  double hylow = 0;
  double hyhi = 0.45;

  char const * xtitle = "x'_{tg} (rad)";
  char const * ytitle = "p_{elas}#theta_{bend} (rad*GeV/c)";

  // Create a TGraph
  auto x = rdf_m.Filter(elascuts.Data()).Range(0,nevents).Take<double>("tgTh");
  auto y = rdf_m.Filter(elascuts.Data()).Range(0,nevents).Take<double>("pethbend");
  TGraph *h_pth = new TGraph(x->size(), x->data(), y->data());
  
  // // create and customize inclusive dpel
  // // TH2F *h_pth = (TH2F*)rdf_m.Filter(elascuts.Data())
  // //   .Histo2D({"h_pth","",nxbins,hxlow,hxhi,nybins,hylow,hyhi},"tgTh","pethbend")->Clone();
  // TProfile *h_pth = (TProfile*)rdf_m.Filter(elascuts.Data())
  //   .Profile1D({"h_pth","",nxbins,hxlow,hxhi},"tgTh","pethbend")->Clone();

  h_pth->GetXaxis()->SetTitle(xtitle);
  h_pth->GetXaxis()->SetLimits(hxlow,hxhi);
  h_pth->GetYaxis()->SetTitle(ytitle);
  h_pth->GetYaxis()->SetRangeUser(hylow,hyhi);
  h_pth->SetMarkerStyle(8);
  h_pth->SetMarkerSize(0.3);
  //h_pth->SetMarkerColor(kGreen-6);

  // // filling the output vector
  obj = {h_pth};
}

//______________________________________________________________________________
double pfitfn (double *x, double *par) {
  return par[0] * (1. + par[1]*x[0]);
}

//______________________________________________________________________________
void mom_recon(int draw_choice) // 1. pethbend vs xptg, 2. dpel 
{
  
  // **
  // PLEASE UPDATE
  // -----
  // Date modified:
  // 07/11/24 (choice 1), 07/12/24 (choice 2)
  // **

  TString infile; 
  TString elcuts;
  bool is_dpel=0;

  // // call the canvas customizer
  PlotCustomizer pcust;

  if (draw_choice==1) {
    // ********________*********
    // Code to pethbend vs thtg plots from both data and MC
    // ---- ** ----
    // sbs 4
    elcuts = "thpq_p<0.008";

    int nevents = 10000; // controls the number of events to be plotted
    double xMin = -0.3;
    double xMax = 0.3;

    // Data
    infile = "~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs30p_model2_pass2.root";

    std::vector<TGraph*> hs4_d; GetpThxpTgHisto(infile,elcuts,hs4_d,nevents,1);
    TCanvas *c4_d = util_pd::TC("c4_d",1,1);
    c4_d->cd(); gStyle->SetErrorX(0);

    hs4_d[0]->GetXaxis()->SetLimits(-0.3, 0.3);
    hs4_d[0]->Draw("AP");
  
    //c4_d->Update();
    TF1* f4_d = new TF1("f4_d",pfitfn,-0.16,0.15,2);
    f4_d->SetNpx(500);
    f4_d->SetLineWidth(4);
    hs4_d[0]->Fit("f4_d","R");
    auto fpm4_d = fit::GetFitParamANDError(f4_d);
    
    TString t4_d = Form("A = %.3f (%d) rad*GeV/c, B = %.3f (%d) rad^{-1}",
			fpm4_d[0].first,util_pd::GetSigDigit(fpm4_d[0].second,3),
			fpm4_d[1].first,util_pd::GetSigDigit(fpm4_d[1].second,3));
    pcust.AddTitleText(t4_d,0.95);
    pcust.customize_canvas(c4_d);

    c4_d->SaveAs("mom_recon_4_pthbnd_data.pdf");

    // MC
    infile = "~/gmn_ana/scripts/simulation/siout/rcOT0_0p65zoff_0p39sf_elas_ana_simc_sbs4_sbs30p_model2.root";

    std::vector<TGraph*> hs4_s; GetpThxpTgHisto(infile,elcuts,hs4_s,nevents,0);
    TCanvas *c4_s = util_pd::TC("c4_s",1,1);
    c4_s->cd(); gStyle->SetErrorX(0);

    hs4_s[0]->GetXaxis()->SetLimits(-0.3, 0.3);
    hs4_s[0]->Draw("AP");
  
    TF1* f4_s = new TF1("f4_s",pfitfn,-0.16,0.15,2);
    f4_s->SetNpx(500);
    f4_s->SetLineWidth(4);
    hs4_s[0]->Fit("f4_s","R");
    auto fpm4_s = fit::GetFitParamANDError(f4_s);
    
    TString t4_s = Form("A = %.3f (%d) rad*GeV/c, B = %.3f (%d) rad^{-1}",
			fpm4_s[0].first,util_pd::GetSigDigit(fpm4_s[0].second,3),
			fpm4_s[1].first,util_pd::GetSigDigit(fpm4_s[1].second,3));
    pcust.AddTitleText(t4_s,0.95);
    pcust.customize_canvas(c4_s);

    c4_s->SaveAs("mom_recon_4_pthbnd_simu.pdf");
  }
  // *********** \/\/\/\/\/ ************

  else if (draw_choice==2) {
    // ********________*********
    // Code to generate dpel plots  
    // ---- ** ----
    TString fbase = "mom_recon";
    // sbs 4
    infile = "~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs30p_model2_pass2.root";
    elcuts = "pCut";
    std::vector<TH1F*> hs4; GetDpelHistos(infile,elcuts,hs4);
    TCanvas *c4 = util_pd::TC("c4",1,1);
    c4->cd();
  
    std::vector<double> hs4_fitR{-0.02,0.04,2.2,2.2};
    TF1 *f4 = fit::fit_1gs_nbg(hs4_fitR,hs4[1]);
    f4->SetLineWidth(2);

    hs4[0]->Draw();
    hs4[1]->Draw("same");
    f4->Draw("same");
  
    TString t4 = Form("Mean = %.2f, #sigma = %.1f %%",f4->GetParameter(1),f4->GetParameter(2)*100);
    pcust.AddTitleText(t4,0.57);
  
    pcust.customize_canvas(c4);

    TLegend *l4=new TLegend(0.64,0.69,0.95,0.9);
    l4->SetTextFont(62);
    l4->AddEntry(hs4[0],"Inclusive ep","f");
    //l4->AddEntry(f0,"Fit","l");
    l4->AddEntry(hs4[1],"Elastic ep","f");
    l4->Draw();

    c4->SaveAs("mom_recon_4.pdf");
    // -----

    // ---- ** ----  
    // sbs 8
    infile = "~/gmn_ana/scripts/pdout/0p77zoff_elas_ana_data_sbs8_sbs70p_model2_pass2.root";
    elcuts = "thpq_p<0.018";
    std::vector<TH1F*> hs8; GetDpelHistos(infile,elcuts,hs8);
    TCanvas *c8 = util_pd::TC("c8",1,1);
    c8->cd();

    std::vector<double> hs8_fitR{-0.02,0.04,2.2,2.2};
    TF1 *f8 = fit::fit_1gs_nbg(hs8_fitR,hs8[1]);
    f8->SetLineWidth(2);

    hs8[0]->Draw();
    hs8[0]->GetYaxis()->SetRangeUser(0,90000);
    hs8[0]->GetYaxis()->SetMaxDigits(6);
    hs8[1]->Draw("same");
    f8->Draw("same");
  
    TString t8 = Form("Mean = %.2f, #sigma = %.1f %%",f8->GetParameter(1),f8->GetParameter(2)*100);
    pcust.AddTitleText(t8,0.57);
  
    pcust.customize_canvas(c8);

    c8->SaveAs("mom_recon_8.pdf");
    // -----

    // ---- ** ----  
    // sbs 9
    infile = "~/gmn_ana/scripts/pdout/0p7zoff_elas_ana_data_sbs9_sbs70p_model2_pass2.root";
    elcuts = "thpq_p<0.018";
    std::vector<TH1F*> hs9; GetDpelHistos(infile,elcuts,hs9);
    TCanvas *c9 = util_pd::TC("c9",1,1);
    c9->cd();

    std::vector<double> hs9_fitR{-0.02,0.04,2.2,2.2};
    TF1 *f9 = fit::fit_1gs_nbg(hs9_fitR,hs9[1]);
    f9->SetLineWidth(2);

    hs9[0]->Draw();
    // hs9[0]->GetYaxis()->SetRangeUser(0,120000);
    // hs9[0]->GetYaxis()->SetMaxDigits(6);
    hs9[1]->Draw("same");
    f9->Draw("same");
  
    TString t9 = Form("Mean = %.2f, #sigma = %.1f %%",f9->GetParameter(1),f9->GetParameter(2)*100);
    pcust.AddTitleText(t9,0.57);
  
    pcust.customize_canvas(c9);

    c9->SaveAs("mom_recon_9.pdf");
    // -----

    // ---- ** ----
    // sbs 11
    infile = "~/gmn_ana/scripts/pdout/1p27zoff_elas_ana_data_sbs11_sbs100p_model2_pass2.root";
    elcuts = "thpq_p<0.009";
    std::vector<TH1F*> hs11; GetDpelHistos(infile,elcuts,hs11);
    TCanvas *c11 = util_pd::TC("c11",1,1);
    c11->cd(); c11->SetLogy();

    std::vector<double> hs11_fitR{-0.02,0.04,2,2.6};
    TF1 *f11 = fit::fit_1gs_nbg(hs11_fitR,hs11[1]);
    f11->SetLineWidth(2);

    hs11[0]->Draw();

    // Set the y-axis range
    hs11[0]->SetMinimum(1e-1);  // Set minimum to a positive value (e.g., 0.1)
    hs11[0]->SetMaximum(1e4);   // Set maximum to desired value (e.g., 10000)
  
    hs11[1]->Draw("same");
    f11->Draw("same");
  
    TString t11 = Form("Mean = %.2f, #sigma = %.1f %%",f11->GetParameter(1),f11->GetParameter(2)*100);
    pcust.AddTitleText(t11,0.57);
  
    pcust.customize_canvas(c11);

    c11->SaveAs("mom_recon_11.pdf");
  }
  // -----
}
