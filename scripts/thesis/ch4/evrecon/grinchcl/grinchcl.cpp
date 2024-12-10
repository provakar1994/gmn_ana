#include "gmn_ana.h"

void grinchcl() {

  // call the canvas customizer
  PlotCustomizer pcust;
  
  TString infile = "~/gmn_ana/scripts/pdout/e1209019_fullreplay_13697_stream0_segALL_ALL.root";

  TFile *f9 = util_pd::ReadRootFile(infile);
  TTree *T9 = (TTree*)f9->Get("T");

  TCanvas *c9 = util_pd::TC("c9",1,1);
  c9->cd(); gStyle->SetPalette(kRainBow);
  
  TH2F *h_grinchcl = new TH2F("h_grinchcl","",200,-0.8,0.8,200,-1.0,1.0);
  T9->Draw("bb.grinch_tdc.clus.x_mean:bb.tr.x[0]+bb.tr.th[0]*0.48>>h_grinchcl","bb.grinch_tdc.nclus>0&&bb.ps.e>0.2&&abs(e.kine.W2-0.88)<0.4");

  // h_grinchcl->GetXaxis()->SetTitle("Best track projection to GRINCH entry window in x (m)");
  // h_grinchcl->GetYaxis()->SetTitle("Best GRINCH cluster mean x (m)");
  h_grinchcl->GetXaxis()->SetTitle("Track projection (m)");
  h_grinchcl->GetYaxis()->SetTitle("GRINCH cluster position (m)");
  h_grinchcl->Draw("colz");

  pcust.customize_canvas(c9);  
  c9->SaveAs("grinchcl_9.pdf");

}
