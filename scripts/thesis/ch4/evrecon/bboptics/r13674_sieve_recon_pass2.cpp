#include "gmn_ana.h"

//______________________________________________________________________________
void r13674_sieve_recon_pass2() {
  TString fbase = "r13674_sieve_recon_pass2";

  TString inFile = Form("%s.root",fbase.Data());
  TString outPlot = inFile; outPlot.ReplaceAll(".root","");
  
  // reading input file
  TFile *f = util_pd::ReadRootFile(inFile);

  // get the histos
  TH1F *h_sieve = (TH1F*)f->Get("h_sieve");
  TH1F *h_vz = (TH1F*)f->Get("h_vz");

  // call the canvas customizer
  PlotCustomizer pcust{0};

  // first canvas
  // h_sieve
  TCanvas *c = util_pd::TC("c",1,1);
  c->cd(); gStyle->SetPalette(kRainBow);

  h_sieve->SetMinimum(-0.5);
  //h_sieve->SetTitle("");
  h_sieve->GetXaxis()->SetTitle("y_{sieve} (m)");
  h_sieve->GetYaxis()->SetTitle("-x_{sieve} (m)");
  h_sieve->Draw("colz");
  
  pcust.customize_canvas(c);
  c->SaveAs(Form("%s_0.png",outPlot.Data()));

  // call the canvas customizer
  PlotCustomizer pcust1;
  
  // h_vz
  TCanvas *c1 = util_pd::TC("c1",1,1);
  c1->cd(); 
  
  //h_vz->SetTitle("");
  h_vz->GetXaxis()->SetTitle("y_{tg} (m)");
  h_vz->SetLineWidth(2);
  h_vz->SetLineColor(kBlack);
  h_vz->SetFillColor(30);
  h_vz->Draw();
  
  std::vector<double> hvz_fitR{-0.02,0.02,1.5,1.5};
  TF1 *fh_vz = fit::fit_1gs_nbg(hvz_fitR,h_vz);
  fh_vz->SetLineWidth(2);

  TString t1 = Form("Data, #sigma = %.1f mm",fh_vz->GetParameter(2)*1000);
  pcust1.AddTitleText(t1,0.53);
  
  pcust1.customize_canvas(c1);
  c1->SaveAs(Form("%s_1.pdf",outPlot.Data()));
  
}

/*
// root commands to regenerate the histogram
cd /w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS9/Optics/rootfiles
ana
TChain c("T")
c.Add("*13674*")
TFile f("~/gmn_ana/scripts/thesis/ch4/r13674_sieve_recon_pass2.root","RECREATE")
f.cd()
// ** h_sieve
c.Draw("-(bb.tr.tg_x+bb.tr.tg_th*1.172):bb.tr.tg_y+bb.tr.tg_ph*1.172>>h_sieve(400,-0.35,0.35,400,-0.35,0.35)","bb.ps.e>0.2&&fEvtHdr.fTrigBits==1&&bb.tr.n==1&&bb.gem.track.nhits>3&&bb.gem.track.ngoodhits[0]>2&&bb.gem.track.chi2ndf[0]<15","colz")
h_sieve->Write()
// ** h_vz
c.Draw("bb.tr.tg_y>>h_vz(150,-.15,.15)","bb.ps.e>0.2&&fEvtHdr.fTrigBits==1&&bb.tr.n==1&&bb.gem.track.nhits>3&&bb.gem.track.ngoodhits[0]>2&&bb.gem.track.chi2ndf[0]<15")
h_vz->Write()
.q
*/
