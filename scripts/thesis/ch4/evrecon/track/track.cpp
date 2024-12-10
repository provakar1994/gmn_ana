#include "gmn_ana.h"

// Plot track residual
void PlotTrackResidual(TFile * f, int sbsconf) {
  // call the canvas customizer
  PlotCustomizer pcust;
  
  // grab histos
  TH1F* h_resid_all = (TH1F*)f->Get("hbb_gem_residu_allhits");

  TCanvas *c = util_pd::TC("c",1,1);
  c->cd();

  h_resid_all->GetXaxis()->SetTitle("Track residuals (m)");
  h_resid_all->SetLineWidth(2);
  h_resid_all->SetLineColor(kBlack);
  h_resid_all->SetFillColor(kGreen-5);
  h_resid_all->Draw();

  TF1* f1 = new TF1("f1","gaus",-0.00015,0.00015);
  f1->SetNpx(500);
  f1->SetLineWidth(3);
  h_resid_all->Fit("f1","R");
  auto fpm = fit::GetFitParamANDError(f1);
    
  TString t = Form("#sigma = %.1f (%d) #mum",
		      fpm[2].first*1e6,util_pd::GetSigDigit(fpm[2].second*1e6,1));
  pcust.AddTitleText(t,0.43);
  pcust.customize_canvas(c);  
  c->SaveAs(Form("trackresid_%d.pdf",sbsconf));
}

// plot GEM cluster stuffs
void PlotGEMClusterInfo(TFile *f, int sbsconf) {
  // call the canvas customizer
  PlotCustomizer pcust;

  // grab histos
  TH1F* h2_clmult_ux = (TH1F*)f->Get("hbb_gem_NclustU_layer");
  TH1F* h_clsize_ux = (TH1F*)f->Get("hbb_gem_clustwidthU");

  // cluster multiplicity
  TCanvas *c1 = util_pd::TC("c1",1,1);
  c1->cd();
  gStyle->SetPalette(kRainBow);
  c1->SetLogz();
  util_pd::SetAxTitles(h2_clmult_ux,"No. of 1D clusters","BB GEM layer");
  h2_clmult_ux->GetYaxis()->SetRangeUser(0,80);
  h2_clmult_ux->Draw("colz");
  pcust.customize_canvas(c1);
  c1->SaveAs(Form("gemclmult_%d.pdf",sbsconf));

  // cluster size
  TCanvas *c2 = util_pd::TC("c2",1,1);
  c2->cd();
   h_clsize_ux->SetFillColor(kOrange);
  util_pd::SetAxTitles(h_clsize_ux,"","No. of strips in 1D cluster");
  h_clsize_ux->Draw();
  pcust.customize_canvas(c2);
  c2->SaveAs(Form("gemclsize_%d.pdf",sbsconf));
}

void track() {

  TString infile14 = "~/gmn_ana/scripts/pdout/e1209019_fullreplay_13243_stream0_segALL_ALL.root";  
  TFile *f14 = util_pd::ReadRootFile(infile14);

  PlotTrackResidual(f14,14);

  //PlotGEMClusterInfo(f14,14);

}
