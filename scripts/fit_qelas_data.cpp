/* 
   ...
   E.g. Config. File: ...
   -----
   P. Datta  Created  01-23-2023 
*/
#include <iostream>

#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"

#include "../include/gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

/* <============ Defining fit fuctions ============> */
//gaussian signal peak
double psignal_fit (double *x, double *par) {
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.));
}
double nsignal_fit (double *x, double *par) {
  return par[0]*exp(-0.5*pow((x[0]-par[1])/par[2],2.));
}
//poly for background
double bg_fit (double *x, double *par) {
  return par[0] + par[1]*x[0] + par[2]*pow(x[0],2) + par[3]*pow(x[0],3) + par[4]*pow(x[0],4);
}
//global fit function
double total_fit (double *x, double *par) {
  return psignal_fit(x,&par[0]) + nsignal_fit(x,&par[3]) + bg_fit(x,&par[6]);
}
/* <====^^^^==== Defining fit fuctions ====^^^^====> */

int fit_qelas_data (const char *configfilename, std::string filebase="pdout/fit_qelas_data")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  std::cout << sbsconf;

  // creating the name of input file
  int model = jmgr->GetValueFromKey<int>("model");
  int pass = jmgr->GetValueFromKey<int>("replay_pass"); 
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::string rfname = Form("%s/qelas_ana_data_sbs%d_sbs%dp_model%d_pass%d*.root", 
			    rootfile_dir.c_str(), sbsconf.GetSBSconf(), sbsconf.GetSBSmag(), model, pass);

  // parsing trees
  TChain *C = new TChain("Tout");
  C->Add(rfname.c_str());
  if (C->GetEntries()==0) {std::cerr << "*!* No ROOT file!" << std::endl; throw;}
  
  // getting useful tree branches
  bool WCut;              C->SetBranchAddress("WCut", &WCut);
  bool pCut;              C->SetBranchAddress("pCut", &pCut);
  bool nCut;              C->SetBranchAddress("nCut", &nCut);
  bool fiduCut;           C->SetBranchAddress("fiduCut", &fiduCut);
  double W;               C->SetBranchAddress("W", &W);
  double W2;              C->SetBranchAddress("W2", &W2);
  double thetapq_p;       C->SetBranchAddress("thetapq_p", &thetapq_p);
  double thetapq_n;       C->SetBranchAddress("thetapq_n", &thetapq_n);
  double xHCAL;           C->SetBranchAddress("xHCAL", &xHCAL);
  double yHCAL;           C->SetBranchAddress("yHCAL", &yHCAL);
  double xHCAL_exp;       C->SetBranchAddress("xHCAL_exp", &xHCAL_exp);
  double yHCAL_exp;       C->SetBranchAddress("yHCAL_exp", &yHCAL_exp);
  double dx;              C->SetBranchAddress("dx", &dx);
  double dy;              C->SetBranchAddress("dy", &dy);

  // defining the oitput root file
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d_pass%d.root", 
			 filebase.c_str(), sbsconf.GetSBSconf(), sbsconf.GetSBSmag(), model, pass);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining interesting histograms
  vector<double> hW2_lim; jmgr->GetVectorFromKey<double>("h_W2_lims", hW2_lim);
  TH1F *h_W2 = new TH1F("h_W2",";W^{2} (GeV^{2})",int(hW2_lim[0]),hW2_lim[1],hW2_lim[2]);
  TH1F *h_W2_nofCut = new TH1F("h_W2_nofCut",";W^{2} (GeV^{2})",int(hW2_lim[0]),hW2_lim[1],hW2_lim[2]);
  TH1F *h_W2_n = new TH1F("h_W2_n",";W^{2} (GeV^{2})",int(hW2_lim[0]),hW2_lim[1],hW2_lim[2]);
  TH1F *h_W2_p = new TH1F("h_W2_p",";W^{2} (GeV^{2})",int(hW2_lim[0]),hW2_lim[1],hW2_lim[2]);
  vector<double> hdx_lim; jmgr->GetVectorFromKey<double>("h_dxHCAL_lims", hdx_lim);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; xHCAL_{obs} - xHCAL_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);  
  TH1F *h_dxHCAL_nofCut = new TH1F("h_dxHCAL_nofCut","; xHCAL_{obs} - xHCAL_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);  
  vector<double> hdy_lim; jmgr->GetVectorFromKey<double>("h_dyHCAL_lims", hdy_lim);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; yHCAL_{obs} - yHCAL_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_dyHCAL_nofCut = new TH1F("h_dyHCAL_nofCut","; yHCAL_{obs} - yHCAL_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH2F *h2_dxdyHCAL = util_pd::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_xyHCAL_p = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_p");
  TH2F *h2_xyHCAL_n = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_n");
  TH2F *h2_xyHCAL_p_nf = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_p_nf");
  TH2F *h2_xyHCAL_n_nf = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_n_nf");

  // reading Cuts
  vector<double> W2cut; jmgr->GetVectorFromKey<double>("W2_cut", W2cut);
  vector<double> dycut; jmgr->GetVectorFromKey<double>("dy_cut", dycut);
  double thetapq_n_min = jmgr->GetValueFromKey<double>("thetapq_n_min_rad");
  double thetapq_p_min = jmgr->GetValueFromKey<double>("thetapq_p_min_rad");

  // defining HCAL cuts
  double sbs_kick = jmgr->GetValueFromKey<double>("sbs_kick");
  vector<double> dx_p; jmgr->GetVectorFromKey<double>("dx_p", dx_p);
  vector<double> dy_p; jmgr->GetVectorFromKey<double>("dy_p", dy_p);
  double Nsigma_cut_dx_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_p");
  double Nsigma_cut_dy_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_p");
  vector<double> dx_n; jmgr->GetVectorFromKey<double>("dx_n", dx_n);
  vector<double> dy_n; jmgr->GetVectorFromKey<double>("dy_n", dy_n);
  double Nsigma_cut_dx_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_n");
  double Nsigma_cut_dy_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_n");
  vector<double> hcal_active_area = cut::hcal_active_area_data(1,1); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p[1], dx_n[1], dy_p[1], hcal_active_area);

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  long nevent = 0, nevents = C->GetEntries(); 
  int treenum = 0, currenttreenum = 0;
  while (C->GetEntry(nevent++)) {
   
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    // filling dx histos
    bool W2Cut = W2 > W2cut[0] && W2 < W2cut[1];
    bool dyCut = dy > dycut[0] && dy < dycut[1];
    // W2 cut
    if (W2Cut) {
      // dy cut
      if (dyCut) {
	// fiducial cut
	if (fiduCut) {
	  h_dxHCAL->Fill(dx);
	} 
	h_dxHCAL_nofCut->Fill(dx);
      }
      // fiducial cut
      if (fiduCut) {
	h_dyHCAL->Fill(dy);
	h2_dxdyHCAL->Fill(dy, dx);
	if (pCut) h2_xyHCAL_p->Fill(yHCAL_exp, xHCAL_exp - sbs_kick);
	if (nCut) h2_xyHCAL_n->Fill(yHCAL_exp, xHCAL_exp);
      }
      h_dyHCAL_nofCut->Fill(dy);
      if (pCut) h2_xyHCAL_p_nf->Fill(yHCAL_exp, xHCAL_exp - sbs_kick);
      if (nCut) h2_xyHCAL_n_nf->Fill(yHCAL_exp, xHCAL_exp);
    }

    // filling W2 histos
    bool thetapq_pCut = thetapq_p < thetapq_p_min;
    bool thetapq_nCut = thetapq_n < thetapq_n_min;
    // fiducial cut
    if (fiduCut) {
      h_W2->Fill(W2);
      // thetapq_p cut
      if (thetapq_pCut) {
	h_W2_p->Fill(W2);
      } 
      // thetapq_n cut
      if (thetapq_nCut) {
	h_W2_n->Fill(W2);
      }
    }
    h_W2_nofCut->Fill(W2);
    
  } // event loop
  std::cout << std::endl << std::endl;

  TCanvas *cdx = new TCanvas("cdx", "cdx", 1200, 600);
  cdx->Divide(2,1);
  
  cdx->cd(1);
  h_dxHCAL->SetTitle(Form("SBS-%d",sbsconf.GetSBSconf()));
  h_dxHCAL->SetLineColor(kBlack); h_dxHCAL->SetLineWidth(2);
  h_dxHCAL->SetLineStyle(9); h_dxHCAL->SetFillColor(18);
  h_dxHCAL->Draw();

  // let's try fitting =========== ******* =========
  double par[11], parf[11];
  vector<double> parGuess; jmgr->GetVectorFromKey<double>("initial_fit_parmas", parGuess);
  TF1* total = new TF1("total", total_fit, hdx_lim[1], hdx_lim[2], 11);
  total->SetNpx(500);
  total->SetLineColor(kRed);
  total->SetParameters(&parGuess[0]);
  h_dxHCAL->Fit("total", "RV+", "ep");
  total->GetParameters(parf);

  // plotting the weighted signal peaks and bg peak
  TF1* psignal = new TF1("psignal", psignal_fit, hdx_lim[1], hdx_lim[2], 3);
  psignal->SetNpx(500);
  psignal->SetParameters(&parf[0]);
  psignal->SetLineColor(kBlue); psignal->SetFillColorAlpha(kBlue, 0.35);
  int pCount = (int)psignal->Integral(hdx_lim[1],hdx_lim[2])/h_dxHCAL->GetBinWidth(1);
  TF1* nsignal = new TF1("nsignal", nsignal_fit, hdx_lim[1], hdx_lim[2], 3);
  nsignal->SetNpx(500);
  nsignal->SetParameters(&parf[3]);
  nsignal->SetLineColor(kGreen);
  int nCount = (int)nsignal->Integral(hdx_lim[1],hdx_lim[2])/h_dxHCAL->GetBinWidth(1);
  TF1* bg = new TF1("bg", bg_fit, hdx_lim[1], hdx_lim[2], 5);
  bg->SetNpx(500);
  bg->SetParameters(&parf[6]);
  bg->SetLineColor(kMagenta);
  double tot_elas = nCount + pCount;

  // draw the legend
  TLegend *legend=new TLegend(0.60,0.66,0.88,0.85);
  legend->SetTextFont(42);
  //legend->SetTextSize(0.02);
  legend->AddEntry(h_dxHCAL,"Data","lf");
  legend->AddEntry(total,"Global Fit","l");
  legend->AddEntry(psignal,"p Signal Fit","l");
  legend->AddEntry(nsignal,"n Signal Fit","l");
  legend->AddEntry(bg,"Background fit","l");

  cdx->cd(2);
  gPad->SetTickx(); gPad->SetTicky(); 
  h_dxHCAL->Draw();
  psignal->Draw("same");
  nsignal->Draw("same");
  bg->Draw("same");
  legend->Draw();

  // ********* === No fiducial cut
  TCanvas *cdx_nf = new TCanvas("cdx_nf", "cdx_nf", 1200, 600);
  cdx_nf->Divide(2,1);
  
  cdx_nf->cd(1);
  h_dxHCAL_nofCut->SetTitle(Form("SBS-%d (No Fiducial Cut)",sbsconf.GetSBSconf()));
  h_dxHCAL_nofCut->SetLineColor(kBlack); h_dxHCAL_nofCut->SetLineWidth(2);
  h_dxHCAL_nofCut->Draw();

  // let's try fitting =========== ******* =========
  TF1* total_nf = new TF1("total_nf", total_fit, hdx_lim[1], hdx_lim[2], 11);
  total_nf->SetNpx(500);
  total_nf->SetLineColor(kRed);
  total_nf->SetParameters(&parGuess[0]);
  h_dxHCAL_nofCut->Fit("total_nf", "RV+", "ep");
  total_nf->GetParameters(parf);

  // plotting the weighted signal peaks and bg peak
  TF1* psignal_nf = new TF1("psignal_nf", psignal_fit, hdx_lim[1], hdx_lim[2], 3);
  psignal_nf->SetNpx(500);
  psignal_nf->SetParameters(&parf[0]);
  psignal_nf->SetLineColor(kBlue);
  int pCount_nf = (int)psignal_nf->Integral(hdx_lim[1],hdx_lim[2])/h_dxHCAL_nofCut->GetBinWidth(1);
  TF1* nsignal_nf = new TF1("nsignal_nf", nsignal_fit, hdx_lim[1], hdx_lim[2], 3);
  nsignal_nf->SetNpx(500);
  nsignal_nf->SetParameters(&parf[3]);
  nsignal_nf->SetLineColor(kGreen);
  int nCount_nf = (int)nsignal_nf->Integral(hdx_lim[1],hdx_lim[2])/h_dxHCAL_nofCut->GetBinWidth(1);
  TF1* bg_nf = new TF1("bg_nf", bg_fit, hdx_lim[1], hdx_lim[2], 5);
  bg_nf->SetNpx(500);
  bg_nf->SetParameters(&parf[6]);
  bg_nf->SetLineColor(kMagenta);
  double tot_elas_nf = nCount_nf + pCount_nf;

  // draw the legend
  TLegend *legend_nf=new TLegend(0.60,0.66,0.88,0.85);
  legend_nf->SetTextFont(42);
  //legend_nf->SetTextSize(0.02);
  legend_nf->AddEntry(h_dxHCAL_nofCut,"Data","l");
  legend_nf->AddEntry(total_nf,"Global Fit","l");
  legend_nf->AddEntry(psignal_nf,"p Signal Fit","l");
  legend_nf->AddEntry(nsignal_nf,"n Signal Fit","l");
  legend_nf->AddEntry(bg_nf,"Background fit","l");

  cdx_nf->cd(2);
  gPad->SetTickx(); gPad->SetTicky(); 
  h_dxHCAL_nofCut->Draw();
  psignal_nf->Draw("same");
  nsignal_nf->Draw("same");
  bg_nf->Draw("same");
  legend_nf->Draw();

  // ********* === Let's plot W2 now
  TCanvas *cW2 = new TCanvas("cW2", "cW2", 1200, 600);
  cW2->Divide(2,1);
  
  cW2->cd(1);
  gPad->SetTickx(); gPad->SetTicky();
  h_W2->SetTitle(Form("SBS-%d",sbsconf.GetSBSconf()));
  h_W2->SetLineColor(kBlack);
  h_W2->Draw();
  h_W2_p->Draw("same"); h_W2_n->Draw("same");

  cW2->cd(2);
  gPad->SetTickx(); gPad->SetTicky();
  h_W2_p->SetTitle(Form("SBS-%d",sbsconf.GetSBSconf()));
  h_W2_p->SetLineColor(kGreen); h_W2_p->SetFillColorAlpha(30, 0.1);
  h_W2_p->Draw();
  h_W2_n->SetLineColor(kRed); h_W2_n->SetFillColorAlpha(46, 0.1);
  h_W2_n->Draw("same");

  // plotting p and n spots at the face of HCAL
  TCanvas *cxyHCAL = new TCanvas("cxyHCAL", "cxyHCAL", 1200, 600);
  cxyHCAL->Divide(2,1);
  
  cxyHCAL->cd(1);
  gPad->SetTickx(); gPad->SetTicky();
  h2_dxdyHCAL->SetTitle(Form("SBS-%d",sbsconf.GetSBSconf()));
  h2_dxdyHCAL->Draw("colz");

  cxyHCAL->cd(2);
  gPad->SetTickx(); gPad->SetTicky();
  h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p[0], dx_p[0], Nsigma_cut_dy_p*dy_p[1], Nsigma_cut_dx_p*dx_p[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n[0], dx_n[0], Nsigma_cut_dy_n*dy_n[1], Nsigma_cut_dx_n*dx_n[1], 0,360,0);


  // ******** === let's plot elastic envelopes at the face of HCAL
  TCanvas *celEnv = new TCanvas("celEnv", "celEnv", 1200, 1000);
  celEnv->Divide(2,2);

  celEnv->cd(1);
  h2_xyHCAL_p->SetTitle(Form("p | SBS-%d",sbsconf.GetSBSconf()));
  h2_xyHCAL_p->Draw("colz");
  util_pd::DrawArea(hcal_active_area,2,4,9);
  util_pd::DrawArea(hcal_safety_margin,4,4,9);

  celEnv->cd(2);
  h2_xyHCAL_n->SetTitle(Form("n | SBS-%d",sbsconf.GetSBSconf()));
  h2_xyHCAL_n->Draw("colz");
  util_pd::DrawArea(hcal_active_area,2,4,9);
  util_pd::DrawArea(hcal_safety_margin,4,4,9);

  celEnv->cd(3);
  h2_xyHCAL_p_nf->SetTitle(Form("p | SBS-%d (No Fiducial Cut)",sbsconf.GetSBSconf()));
  h2_xyHCAL_p_nf->Draw("colz");
  util_pd::DrawArea(hcal_active_area,2,4,9);
  util_pd::DrawArea(hcal_safety_margin,4,4,9);

  celEnv->cd(4);
  h2_xyHCAL_n_nf->SetTitle(Form("n | SBS-%d (No Fiducial Cut)",sbsconf.GetSBSconf()));
  h2_xyHCAL_n_nf->Draw("colz");
  util_pd::DrawArea(hcal_active_area,2,4,9);
  util_pd::DrawArea(hcal_safety_margin,4,4,9);

  // time to summarize the findings
  double nCount_ploss = (nCount_nf-nCount)*100.0 / nCount_nf;
  double pCount_ploss = (pCount_nf-pCount)*100.0 / pCount_nf;
  double total_ploss = (tot_elas_nf-tot_elas)*100.0 / tot_elas_nf;

  TCanvas *cgist = new TCanvas("cgist");
  cgist->cd();

  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText(Form("SBS-%d SBS %dp Field Data",sbsconf.GetSBSconf(),sbsmag));
  pt->AddText(Form(" n count = %d",nCount));
  pt->AddText(Form(" p count = %d",pCount));
  pt->AddText(Form(" n count (no fiducial cut) = %d",nCount_nf));
  pt->AddText(Form(" p count (no fiducial cut) = %d",pCount_nf));
  pt->AddText(Form(" loss of n events = %.1fp",nCount_ploss));
  pt->AddText(Form(" loss of p events = %.1fp",pCount_ploss));
  pt->AddText(Form(" loss of total elastics = %.1fp",total_ploss));
  TText *t1 = pt->GetLineWith("SBS");
  t1->SetTextColor(kBlue);
  pt->Draw();

  cout << endl << "------" << endl;
  cout << " SBS-" << sbsconf.GetSBSconf() << " SBS " << sbsmag << "% data.." << endl;
  cout << " n count = " << nCount << endl;
  cout << " p count = " << pCount << endl;
  cout << " n count (no fiducial cut) = " << nCount_nf << endl;
  cout << " p count (no fiducial cut) = " << pCount_nf << endl;
  cout << " % loss (n events) = " << nCount_ploss << endl;
  cout << " % loss (p events) = " << pCount_ploss << endl;
  cout << " % loss (total elastics) = " << total_ploss << endl;
  cout << " Output file : " << outFile << endl;
  cout << "------" << endl << endl;

  cdx->Write();
  cdx_nf->Write();
  cW2->Write();
  cxyHCAL->Write();
  celEnv->Write();
  cgist->Write();
  fout->Write();
  delete jmgr;
  return 0;
}
