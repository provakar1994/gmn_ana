#include "gmn_ana.h"

//______________________________________________________________________________
void customize_nblk_h(TH1F* h) {
  h->SetLineWidth(2);
  h->SetLineColor(kBlack);
  h->SetFillColor(kOrange);
  h->GetXaxis()->SetTitle("No. of Blocks in Cluster");
}

//______________________________________________________________________________
void customize_mult_h(TH1F* h) {
  h->SetLineWidth(2);
  h->SetLineColor(kBlack);
  h->SetFillColor(kOrange);
  h->GetXaxis()->SetTitle("No. of Clusters in an Event");
}

//______________________________________________________________________________
void hcalclus() {

  TString infile, infile2; 
  TString elcuts, elcuts2;
  
  // call the canvas customizer
  PlotCustomizer pcust;
  
  // // sbs4
  // infile = "~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs30p_model2_pass2.root";
  // elcuts = "idblkHCAL_aclN!=0&&WCut&&abs(dy)<0.3&&nblkHCAL_acl[0]<20&&eHCAL>0";

  // TFile *f4 = util_pd::ReadRootFile(infile);
  // TTree *T4 = (TTree*)f4->Get("Tout");

  // TCanvas *c4 = util_pd::TC("c4",1,1);
  // c4->cd();

  // TH1F *h4_nblk = new TH1F("h4_nblk","",20,0.5,20.5);  
  // T4->Draw("nblkHCAL_acl[idclHCAL_intime]>>h4_nblk",elcuts.Data());
  // customize_nblk_h(h4_nblk); h4_nblk->Draw();

  // pcust.customize_canvas(c4);
  // c4->SaveAs("hcalclus_4.pdf");
  // // ----

  // sbs4
  infile = "~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs30p_model2_pass2.root";
  elcuts = "idblkHCAL_aclN!=0&&WCut&&abs(dy)<0.3&&nblkHCAL_acl[0]<20&&eHCAL>0";

  infile2 = "/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS4/LH2/rootfiles/*11547*";
  elcuts2 = "fEvtHdr.fTrigBits==1&&bb.tr.n>0&&bb.gem.track.nhits>2&&abs(bb.tr.vz[0])<.08&&bb.tr.p[0]>2&&bb.ps.e>0.2&&bb.gem.track.ngoodhits[0]>2&&bb.gem.track.chi2ndf[0]<15&&sbs.hcal.e>0";
  
  TFile *f4 = util_pd::ReadRootFile(infile);
  TTree *T4 = (TTree*)f4->Get("Tout");

  TChain *Ch4 = new TChain("T");
  Ch4->Add(infile2);
  
  TCanvas *c4 = util_pd::TC("c4",1,2);
  c4->cd(1);
  TH1F *h4_mult = new TH1F("h4_mult","",20,0.5,20.5);  
  Ch4->Draw("sbs.hcal.nclus>>h4_mult",elcuts2.Data());
  customize_mult_h(h4_mult); h4_mult->Draw();  
 
  c4->cd(2);
  TH1F *h4_nblk = new TH1F("h4_nblk","",20,0.5,20.5);  
  T4->Draw("nblkHCAL_acl[idclHCAL_intime]>>h4_nblk",elcuts.Data());
  customize_nblk_h(h4_nblk); h4_nblk->Draw();

  pcust.customize_canvas(c4);
  c4->SaveAs("hcalclus_4.pdf");
  // ----  

  // // sbs11
  // infile = "~/gmn_ana/scripts/pdout/1p27zoff_elas_ana_data_sbs11_sbs100p_model2_pass2.root";
  // elcuts = "idblkHCAL_aclN!=0&&WCut&&pCut&&nblkHCAL_acl[0]<20&&abs(coinT_ADC_c)<3.5*1.4";

  // TFile *f11 = util_pd::ReadRootFile(infile);
  // TTree *T11 = (TTree*)f11->Get("Tout");

  // TCanvas *c11 = util_pd::TC("c11",1,1);
  // c11->cd();

  // TH1F *h11_nblk = new TH1F("h11_nblk","",20,0.5,20.5);  
  // T11->Draw("nblkHCAL_acl[idclHCAL_intime]>>h11_nblk",elcuts.Data());
  // customize_nblk_h(h11_nblk); h11_nblk->Draw();

  // pcust.customize_canvas(c11);  
  // c11->SaveAs("hcalclus_11.pdf");
  // // -----


  // sbs11
  infile = "~/gmn_ana/scripts/pdout/1p27zoff_elas_ana_data_sbs11_sbs100p_model2_pass2.root";
  elcuts = "idblkHCAL_aclN!=0&&WCut&&pCut&&nblkHCAL_acl[0]<20&&abs(coinT_ADC_c)<3.5*1.4";

  infile2 = "/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS11/LH2/rootfiles/*12525*";
  elcuts2 = "fEvtHdr.fTrigBits==1&&bb.tr.n>0&&bb.gem.track.nhits>2&&abs(bb.tr.vz[0])<.08&&bb.tr.p[0]>2&&bb.ps.e>0.2&&bb.gem.track.ngoodhits[0]>2&&bb.gem.track.chi2ndf[0]<15&&sbs.hcal.e>0";
  
  TFile *f11 = util_pd::ReadRootFile(infile);
  TTree *T11 = (TTree*)f11->Get("Tout");

  TChain *Ch11 = new TChain("T");
  Ch11->Add(infile2);
  
  TCanvas *c11 = util_pd::TC("c11",1,2);
  c11->cd(1);
  TH1F *h11_mult = new TH1F("h11_mult","",71,59.5,130.5);  
  Ch11->Draw("sbs.hcal.nclus>>h11_mult",elcuts2.Data());
  customize_mult_h(h11_mult); h11_mult->Draw();  

  c11->cd(2);
  TH1F *h11_nblk = new TH1F("h11_nblk","",20,0.5,20.5);  
  T11->Draw("nblkHCAL_acl[idclHCAL_intime]>>h11_nblk",elcuts.Data());
  customize_nblk_h(h11_nblk); h11_nblk->Draw();

  pcust.customize_canvas(c11);  
  c11->SaveAs("hcalclus_11.pdf");
  // -----
  

  
  // sbs14
  infile = "~/gmn_ana/scripts/pdout/1p13zoff_elas_ana_data_sbs14_sbs70p_model2_pass2.root";
  elcuts = "idblkHCAL_aclN!=0&&WCut&&abs(dy)<0.3&&ARCut&&coinTADCCut&&eHCAL>0.2";

  TFile *f14 = util_pd::ReadRootFile(infile);
  TTree *T14 = (TTree*)f14->Get("Tout");

  TCanvas *c14 = util_pd::TC("c14",1,1);
  c14->cd();

  TH1F *h14_dx_intime = new TH1F("h14_dx_intime","",120,-1.2,-0.4);  
  T14->Draw("xHCAL_acl[idclHCAL_intime]-xHCAL_exp>>h14_dx_intime",elcuts.Data());
  TH1F *h14_dx_sthpq = new TH1F("h14_dx_sthpq","",120,-1.2,-0.4);  
  T14->Draw("xHCAL_acl[idclHCAL_sthpq_p]-xHCAL_exp>>h14_dx_sthpq",elcuts.Data());
  TH1F *h14_dx_htote = new TH1F("h14_dx_htote","",120,-1.2,-0.4);  
  T14->Draw("xHCAL_acl[idclHCAL_htote]-xHCAL_exp>>h14_dx_htote",elcuts.Data());

  //customize_nblk_h(h14_nblk);
  h14_dx_intime->Draw();
  h14_dx_intime->SetLineWidth(2);
  h14_dx_intime->SetLineColor(kBlack);
  h14_dx_intime->SetFillColor(kGray);
  h14_dx_intime->GetXaxis()->SetTitle("x_{HCAL}^{obs}-x_{HCAL}^{exp} (m)");
  h14_dx_sthpq->Draw("same");
  h14_dx_sthpq->SetLineWidth(2);
  h14_dx_sthpq->SetLineColor(kBlack);
  h14_dx_sthpq->SetFillColor(kGray+1);
  h14_dx_htote->Draw("same");
  h14_dx_htote->SetLineWidth(2);
  h14_dx_htote->SetLineColor(kBlack);
  h14_dx_htote->SetFillColor(kGray+3);
  
  pcust.customize_canvas(c14);  

  TLegend *l14=new TLegend(0.64,0.69,0.95,0.9);
  l14->SetTextFont(62);
  l14->AddEntry(h14_dx_intime,"In-time","f");
  l14->AddEntry(h14_dx_sthpq,"Smallest #theta_{pq}","f");
  l14->AddEntry(h14_dx_htote,"Highest total energy","f");
  l14->Draw();

  c14->SaveAs("hcalclus_14.pdf");
  // -----  

}
