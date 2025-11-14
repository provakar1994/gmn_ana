#include "gmn_ana.h"

//______________________________________________________________________________
void customize_nblk_h(TH1F* h) {
  h->SetLineWidth(2);
  h->SetLineColor(kBlack);
  // h->SetFillColor(kOrange); // used for thesis
  h->SetFillColor(30); // used for NIM
  h->GetXaxis()->SetTitle("No. of Blocks in Cluster");
}

//______________________________________________________________________________
void customize_mult_h(TH1F* h) {
  h->SetLineWidth(2);
  h->SetLineColor(kBlack);
  // h->SetFillColor(kOrange);
  h->SetFillColor(30);  // used for NIM
  h->GetXaxis()->SetTitle("No. of Clusters in an Event");
}

//______________________________________________________________________________
void bbcalclus() {

  TString infile, infile2; 
  TString elcuts, elcuts2;
  
  // call the canvas customizer
  PlotCustomizer pcust;
  

  // PS
  infile2 = "/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS8/LH2/rootfiles/*13486*";
  elcuts2 = "fEvtHdr.fTrigBits==1&&bb.tr.n>0&&bb.gem.track.nhits>2&&abs(bb.tr.vz[0])<.08&&bb.tr.p[0]>2.9&&bb.ps.e>0.2";
  
  TChain *ChPS = new TChain("T");
  ChPS->Add(infile2);
  
  TCanvas *cPS = util_pd::TC("cPS",1,2);
  gStyle->SetLineScalePS(2.5);
  cPS->cd(1);
  gPad->SetLogy();
  TH1F *hPS_mult = new TH1F("hPS_mult","",8,0.5,8.5);  
  ChPS->Draw("bb.ps.nclus>>hPS_mult",elcuts2.Data());
  customize_mult_h(hPS_mult); hPS_mult->Draw();  
 
  cPS->cd(2);
  TH1F *hPS_nblk = new TH1F("hPS_nblk","",8,0.5,8.5);  
  ChPS->Draw("bb.ps.nblk>>hPS_nblk",elcuts.Data());
  customize_nblk_h(hPS_nblk); hPS_nblk->Draw();

  pcust.customize_canvas(cPS);
  cPS->SaveAs("NIM_bbcalclus_PS.pdf");
  // ----  


  // SH
  TChain *ChSH = new TChain("T");
  ChSH->Add(infile2);
  
  TCanvas *cSH = util_pd::TC("cSH",1,2);
  gStyle->SetLineScalePS(2.5);
  cSH->cd(1);
  gPad->SetLogy();
  TH1F *hSH_mult = new TH1F("hSH_mult","",8,0.5,8.5);  
  ChSH->Draw("bb.sh.nclus>>hSH_mult",elcuts2.Data());
  customize_mult_h(hSH_mult); hSH_mult->Draw();  
 
  cSH->cd(2);
  TH1F *hSH_nblk = new TH1F("hSH_nblk","",8,0.5,8.5);  
  ChSH->Draw("bb.sh.nblk>>hSH_nblk",elcuts.Data());
  customize_nblk_h(hSH_nblk); hSH_nblk->Draw();

  pcust.customize_canvas(cSH);
  cSH->SaveAs("NIM_bbcalclus_SH.pdf");
  // ----  
 


}
