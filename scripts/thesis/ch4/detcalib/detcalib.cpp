#include "gmn_ana.h"

double xNDCoffset = 0.05;
double y1NDCoffset = 0.05;

void custom_statbox_effi(TPaveStats *st) {
  st->SetBit(TH1::kNoStats);
  st->SetX1NDC(0.64); st->SetY1NDC(0.69); st->SetX2NDC(0.95); st->SetY2NDC(0.9);
}

//____________________________________________________
void PlotBBCALEngResSBS8()
{
  // call the canvas customizer
  PlotCustomizer pcust{1,1};

  TFile *f1 = util_pd::ReadRootFile("/lustre24/expphy/volatile/halla/sbs/pdbforce/bbcalib/hist/pass2/sbs8-sbs70p-set2_prepass2_bbcal_eng_calib_elcut_1.root");
  TH1F *h1= (TH1F*)f1->Get("h_EovP_calib");
  // Reading in ROOT file with cosmic calibration
  TChain *ch = new TChain("T");
  ch->Add("/w/halla-scshelf2102/sbs/pdbforce/hist/e1209019_fullreplay_13486_stream0_ALL_Online_PrePass0.root");  
  
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  h1->SetStats(1);  
  h1->Draw();
  // drawing the second histo
  TH1F *h2 = new TH1F("h2","",200,0.4,1.6);
  ch->Draw("(bb.sh.e+bb.ps.e)/bb.tr.p[0]>>h2","bb.ps.e>0.2","same HIST");
  h2->Scale(1/2.8);
  //  
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","#font[32]{E_{BBCAL}/E'_{e}}");
  //
  h2->SetFillColorAlpha(kRed,0.4);
  h2->SetLineColor(kRed+1); 
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",0.87,1.1);
  h1->Fit(ffn,"R");
  ffn->SetLineColor(kRed);
  ffn->SetLineWidth(4);  
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);  
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  // custom_statbox_effi(s1c);
  // s1c->SetX1NDC(0.62);
  // s1c->SetY1NDC(0.5);
  s1c->SetX1NDC(0.1+xNDCoffset);
  s1c->SetY1NDC(0.5);
  s1c->SetX2NDC(0.48);
  s1c->SetY2NDC(0.9);
  //
  TLegend *l1=new TLegend(0.6,0.75,0.95,0.9);
  l1->SetTextFont(62);
  l1->AddEntry(h2,"Calibration w/ Cosmic Data","f");
  l1->AddEntry(h1,"Calibration w/ #font[32]{l}H_{2} Data","f");
  l1->Draw();  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_bbcaleovp_8.pdf"));        
}

//____________________________________________________
void PlotBBCALEngResSBS9()
{
  // call the canvas customizer
  PlotCustomizer pcust{1,1};

  TFile *f1 = util_pd::ReadRootFile("/lustre24/expphy/volatile/halla/sbs/pdbforce/bbcalib/hist/pass2/sbs9-sbs70p_prepass2_bbcal_eng_calib_elcut.root");
  TH1F *h1= (TH1F*)f1->Get("h_EovP_calib");
  
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  h1->SetStats(1);  
  h1->Draw();
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","#font[32]{E_{BBCAL}/E'_{e}}");
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",0.85,1.1);
  h1->Fit(ffn,"R");
  ffn->SetLineColor(kRed);
  ffn->SetLineWidth(4);  
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);  
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  // custom_statbox_effi(s1c);
  // s1c->SetX1NDC(0.62);
  // s1c->SetY1NDC(0.5);
  s1c->SetX1NDC(0.1+xNDCoffset);
  s1c->SetY1NDC(0.5);
  s1c->SetX2NDC(0.48);
  s1c->SetY2NDC(0.9);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_bbcaleovp_9.pdf"));        
}

//____________________________________________________
void PlotHCALdyRes()
{
  int sbsconf = 14;

  // call the canvas customizer
  PlotCustomizer pcust{1,1};

  TChain *ch = new TChain("Tout");
  ch->Add("/lustre24/expphy/volatile/halla/sbs/pdbforce/gmn_ana/pdout/1p13zoff_elas_ana_data_sbs14_sbs70p_model1_pass2.root"); //1369
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-1.1,1.1);
  ch->Draw("dy>>h1","WCut&&eHCAL>0&&bbfiduCut&&ARCut");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","#font[32]{#Deltay} (m)");
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",-0.08,0.08);
  h1->Fit(ffn,"R");
  ffn->SetLineColor(kRed);
  ffn->SetLineWidth(4);  
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);  
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  custom_statbox_effi(s1c);
  s1c->SetY1NDC(0.6);
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_hcal_eng_posresolution_dy_%d.pdf",sbsconf));        

}

//____________________________________________________
void PlotHCALdxRes()
{
  int sbsconf = 14;

  // call the canvas customizer
  PlotCustomizer pcust{1,1};

  TChain *ch = new TChain("Tout");
  ch->Add("/lustre24/expphy/volatile/halla/sbs/pdbforce/gmn_ana/pdout/1p13zoff_elas_ana_data_sbs14_sbs70p_model1_pass2.root"); //1369
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-2.5,1);
  ch->Draw("dx>>h1","WCut&&eHCAL>0&&abs(dy)<0.3&&bbfiduCut&&ARCut");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","#font[32]{#Deltax} (m)");
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",-0.92,-0.72);
  h1->Fit(ffn,"R");
  ffn->SetLineColor(kRed);
  ffn->SetLineWidth(4);  
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);  
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  custom_statbox_effi(s1c);
  s1c->SetY1NDC(0.6);
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_hcal_eng_posresolution_dx_%d.pdf",sbsconf));        

}

//____________________________________________________
void PlotHCALEngAlign()
{
  int sbsconf = 8;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch = new TChain("Tout");
  ch->Add("~/gmn_ana/scripts/pdout/0p77zoff_elas_ana_data_sbs8_sbs70p_model1_pass2.root");
  //
  TCanvas *cv = new TCanvas("cv","cv",1200,1000); //util_pd::TC("cv",1,1);
  cv->Divide(1,2);
  //gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd(1);
  TH2F *h1 = new TH2F("h1","",24,0,24,200,0,0.25);
  ch->Draw("eHCAL/nu:rblkHCAL>>h1","eHCAL>0&&WCut&&pCut");
  util_pd::SetAxTitles(h1,"#font[32]{f_{cal} = E_{HCAL}/E^{Kin}_{N}}","HCAL Row");
  // profile histo
  TProfile *pfx1 = h1->ProfileX("pfx1", h1->GetYaxis()->FindBin(0.07), h1->GetYaxis()->FindBin(0.12), "S");
  pfx1->Draw("same");
  pfx1->SetMarkerStyle(20);
  pfx1->SetMarkerColor(kRed);
  //
  cv->cd(2);
  TH2F *h2 = new TH2F("h2","",12,0,12,200,0,0.25);
  ch->Draw("eHCAL/nu:cblkHCAL>>h2","eHCAL>0&&WCut&&pCut");
  util_pd::SetAxTitles(h2,"#font[32]{f_{cal} = E_{HCAL}/E^{Kin}_{N}}","HCAL Column");
  // profile histo
  TProfile *pfx2 = h2->ProfileX("pfx2", h2->GetYaxis()->FindBin(0.07), h2->GetYaxis()->FindBin(0.12), "S");
  pfx2->Draw("same");
  pfx2->SetMarkerStyle(20);
  pfx2->SetMarkerColor(kRed);  
  //
  pcust.customize_canvas(cv);  
  //
  h1->GetYaxis()->SetTitleSize(0.07);
  h1->GetYaxis()->SetTitleOffset(0.8);  
  h2->GetYaxis()->SetTitleSize(0.07);
  h2->GetYaxis()->SetTitleOffset(0.8);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_hcal_engalign_%d.pdf",sbsconf));      
}

//____________________________________________________
void PlotHCALtimeRes()
{
  int sbsconf = 9;

  // call the canvas customizer
  PlotCustomizer pcust{1,1};

  TChain *ch = new TChain("T");
  ch->Add("/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS9/LH2/rootfiles/*"); //1369
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-30,30);
  ch->Draw("sbs.hcal.clus_blk.atime-sbs.hcal.atimeblk>>h1","bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.065&&sbs.hcal.clus_blk.id!=sbs.hcal.idblk&&bb.gem.track.nhits>3&&bb.tr.p[0]>1.2&&abs(e.kine.W2-0.88)<0.4");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","Time Difference (ns)");
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",-2.5,2.5);
  h1->Fit(ffn,"R");
  ffn->SetLineColor(kRed);
  ffn->SetLineWidth(4);  
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);  
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  custom_statbox_effi(s1c);
  s1c->SetY1NDC(0.6);
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_hcalatime_resolution_%d.pdf",sbsconf));        

}

//____________________________________________________
void PlotHCALAtimeAlign()
{
  int sbsconf = 7;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch = new TChain("Tout");
  ch->Add("~/gmn_ana/scripts/pdout/0p925zoff_elas_ana_data_sbs7_sbs85p_model1_pass2.root");
  //
  TCanvas *cv = new TCanvas("cv","cv",1200,800); //util_pd::TC("cv",1,1);
  //gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH2F *h1 = new TH2F("h1","",288,1,288,200,-35,35);
  ch->Draw("atimeHCAL:idblkHCAL>>h1","eHCAL>0");
  util_pd::SetAxTitles(h1,"#font[32]{t^{ADC}_{HCAL}} (ns)","HCAL Module ID");
  // profile histo
  TProfile *pfx1 = h1->ProfileX("pfx1", h1->GetYaxis()->FindBin(-3), h1->GetYaxis()->FindBin(3), "S");
  pfx1->Draw("same");
  pfx1->SetMarkerStyle(20);
  pfx1->SetMarkerColor(kRed);
  //
  pcust.customize_canvas(cv);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_hcal_atimealign_%d.pdf",sbsconf));      
}

//____________________________________________________
void PlotOutofTimeSHHits()
{
  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch = new TChain("T");
  ch->Add("/cache/halla/sbs/prod/gmn/pass0/SBS4/LH2/rootfiles/*11436*");
  ch->Add("/cache/halla/sbs/prod/gmn/pass0/SBS4/LH2/rootfiles/*11500*");
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH2F *h1 = new TH2F("h1","",200,0,1,200,-40,40);
  ch->Draw("bb.sh.clus_blk.atime-bb.sh.atimeblk:bb.sh.clus_blk.e/bb.sh.eblk>>h1","bb.sh.clus_blk.id!=bb.sh.idblk","colz");
  util_pd::SetAxTitles(h1,"Time Difference (ns)","Energy Fraction");
  //
  pcust.customize_canvas(cv);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_shatime_issue3_ev1_%d.pdf",sbsconf));       
  //
  TCanvas *cv2 = util_pd::TC("cv2",1,1);
  gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv2->cd();
  TH2F *h2 = new TH2F("h2","",200,-0.1,0.1,200,-40,40);
  ch->Draw("bb.sh.clus_blk.atime-bb.sh.atimeblk:bb.tr.y[0]+bb.z_bcp[0]*bb.tr.ph[0]>>h2","bb.sh.clus_blk.id!=bb.sh.idblk&&bb.sh.clus_blk.row==bb.sh.rowblk&&bb.sh.colblk==3","colz");
  util_pd::SetAxTitles(h2,"Time Difference (ns)","#font[32]{y_{Tr}} Projected at Shower (m)");
  //
  h2->GetXaxis()->SetNdivisions(9); 
  //
  pcust.customize_canvas(cv2);  
  //
  cv2->Update();
  cv2->SaveAs(Form("detcalib_shatime_issue3_ev2_%d.pdf",sbsconf));     

}

//____________________________________________________
void PlotRFTime()
{
  int sbsconf = 8;

  // call the canvas customizer
  PlotCustomizer pcust{1,1};

  TChain *ch = new TChain("T");
  ch->Add("/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS8/LH2/rootfiles/*13486*"); //1369
  ch->Add("/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS8/LH2/rootfiles/*13487*"); //1369
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-30,30);
  ch->Draw("bb.hodotdc.clus.tmean[0]+bb.gem.trigtime-fmod(bb.tdctrig.tdc,4.0)-0.53>>h1","bb.tdctrig.tdcelemID==4&&bb.hodotdc.clus.id[0]==44&&bb.tdctrig.tdc>0");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","#font[32]{t^{TDC}_{TH} + t^{TDC}_{BBtrig} -} fmod(#font[32]{t^{TDC}_{RF}}, 4 ns) (ns)");
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",-1.4,1.2);
  h1->Fit(ffn,"R");
  ffn->SetLineColor(kRed);
  ffn->SetLineWidth(4);  
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);  
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  custom_statbox_effi(s1c);
  s1c->SetY1NDC(0.6);
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_thrfstructure_%d.pdf",sbsconf));        

}

//____________________________________________________
void PlotSHPosRes()
{
  int sbsconf = 8;

  // call the canvas customizer
  PlotCustomizer pcust{1,0,0};  

  // grab the trees
  TFile *f1 = util_pd::ReadRootFile("/lustre24/expphy/volatile/halla/sbs/pdbforce/bbcalib/hist/pass2/sbs8-sbs50p_prepass2_bbcal_eng_calib_elcut.root");
  // TTree *Tlhloq = (TTree*)flhloq->Get("Tout");
  TCanvas *c2 = (TCanvas*)f1->Get("c4");

  // Accessing each sub-pad
  TPad *p1 = (TPad*)c2->GetPad(3); // Sub-pad 1 of Pad 1
  TPad *p2 = (TPad*)c2->GetPad(4); // Sub-pad 2 of Pad 1  

  // Create a new canvas for the contents of subpad1
  TCanvas *cv = util_pd::TC("cv",1,2);
  cv->cd(1);
  // Loop through all objects in subpad1
  TIter next(p1->GetListOfPrimitives());
  TObject *obj;
  bool firstHistogram = true;
  std::vector<TH1*> h1s;         // Vector to store histogram pointers
  while ((obj = next())) {
    TObject *objClone = obj->Clone();  // Clone each object
    
    // Determine the draw options
    if (objClone->InheritsFrom("TH1")) {       // For histograms
      TH1* hist = (TH1*)objClone;
      h1s.push_back(hist);        
      if (firstHistogram) {
	objClone->Draw();                  // Draw the first histogram normally
	firstHistogram = false;
      } else {
	objClone->Draw("same");            // Overlay additional histograms
      }
    } else if (objClone->InheritsFrom("TGraph")) {   // For graphs
      objClone->Draw("same");                // Draw graphs on top
    } else {                                   // For other types of objects
	objClone->Draw();                      // Draw normally without "same"
    }
  }
  util_pd::SetAxTitles(h1s[0],"","#font[32]{x_{SH} - x_{Tr}} (m)");
  pcust.AddTitleText("Dispersive, #font[32]{#sigma =} 1.21 cm",0.68);
  //
  cv->cd(2);  
  // Loop through all objects in subpad2
  TIter next2(p2->GetListOfPrimitives());
  TObject *obj2;
  bool firstHistogram2 = true;
  std::vector<TH1*> h1s2;         // Vector to store histogram pointers
  while ((obj2 = next2())) {
    TObject *objClone2 = obj2->Clone();  // Clone each object
    
    // Determine the draw options
    if (objClone2->InheritsFrom("TH1")) {       // For histograms
      TH1* hist = (TH1*)objClone2;
      h1s2.push_back(hist);        
      if (firstHistogram2) {
	objClone2->Draw();                  // Draw the first histogram normally
	firstHistogram2 = false;
      } else {
	objClone2->Draw("same");            // Overlay additional histograms
      }
    } else if (objClone2->InheritsFrom("TGraph")) {   // For graphs
      objClone2->Draw("same");                // Draw graphs on top
    } else {                                   // For other types of objects
	objClone2->Draw();                      // Draw normally without "same"
    }
  }
  util_pd::SetAxTitles(h1s2[0],"","#font[32]{y_{SH} - y_{Tr}} (m)");
  pcust.AddTitleText("Transverse, #font[32]{#sigma =} 1.24 cm",0.69);
  //
  pcust.customize_canvas(cv);
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_bbcaleng_posres_%d.pdf",sbsconf));
}

//____________________________________________________
void PlotPSAtimeAlign()
{
  int sbsconf = 7;

  // call the canvas customizer
  PlotCustomizer pcust{1,0,0};  

  // grab the trees
  TFile *f1 = util_pd::ReadRootFile("/w/halla-scshelf2102/sbs/pdbforce/hist/sbs7_prepass2_atimeOff.root");
  // TTree *Tlhloq = (TTree*)flhloq->Get("Tout");
  TCanvas *c2 = (TCanvas*)f1->Get("c3");

  // Accessing each sub-pad
  TPad *p1 = (TPad*)c2->GetPad(1); // Sub-pad 1 of Pad 1
  TPad *p2 = (TPad*)c2->GetPad(2); // Sub-pad 2 of Pad 1  

  // Create a new canvas for the contents of subpad1
  TCanvas *cv = new TCanvas("cv","cv",1000,800);
  cv->Divide(1,2);
  cv->cd(1);
  // Loop through all objects in subpad1
  TIter next(p1->GetListOfPrimitives());
  TObject *obj;
  bool firstHistogram = true;
  std::vector<TH1*> h1s;         // Vector to store histogram pointers
  while ((obj = next())) {
    TObject *objClone = obj->Clone();  // Clone each object
    
    // Determine the draw options
    if (objClone->InheritsFrom("TH1")) {       // For histograms
      TH1* hist = (TH1*)objClone;
      h1s.push_back(hist);        
      if (firstHistogram) {
	objClone->Draw();                  // Draw the first histogram normally
	firstHistogram = false;
      } else {
	objClone->Draw("same");            // Overlay additional histograms
      }
    } else if (objClone->InheritsFrom("TGraph")) {   // For graphs
      objClone->Draw("same");                // Draw graphs on top
    } else {                                   // For other types of objects
	objClone->Draw();                      // Draw normally without "same"
    }
  }
  util_pd::SetAxTitles(h1s[0],"#font[32]{t^{TDC}_{TH} - t^{ADC}_{PS}} (ns)","PS Module ID");
  pcust.AddTitleText("Before Correction",0.32);
  //
  cv->cd(2);  
  // Loop through all objects in subpad2
  TIter next2(p2->GetListOfPrimitives());
  TObject *obj2;
  bool firstHistogram2 = true;
  std::vector<TH1*> h1s2;         // Vector to store histogram pointers
  while ((obj2 = next2())) {
    TObject *objClone2 = obj2->Clone();  // Clone each object
    
    // Determine the draw options
    if (objClone2->InheritsFrom("TH1")) {       // For histograms
      TH1* hist = (TH1*)objClone2;
      h1s2.push_back(hist);        
      if (firstHistogram2) {
	objClone2->Draw();                  // Draw the first histogram normally
	firstHistogram2 = false;
      } else {
	objClone2->Draw("same");            // Overlay additional histograms
      }
    } else if (objClone2->InheritsFrom("TGraph")) {   // For graphs
      objClone2->Draw("same");                // Draw graphs on top
    } else {                                   // For other types of objects
	objClone2->Draw();                      // Draw normally without "same"
    }
  }
  util_pd::SetAxTitles(h1s2[0],"#font[32]{t^{TDC}_{TH} - t^{ADC}_{PS}} (ns)","PS Module ID");
  pcust.AddTitleText("After Correction",0.3);
  //
  pcust.customize_canvas(cv);
  //
  h1s[0]->GetYaxis()->SetTitleSize(0.07);
  h1s[0]->GetYaxis()->SetTitleOffset(0.8);
  h1s2[0]->GetYaxis()->SetTitleSize(0.07);
  h1s2[0]->GetYaxis()->SetTitleOffset(0.8);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_psatime_%d.pdf",sbsconf));
}

//____________________________________________________
void PlotSHAtimeAlign()
{
  int sbsconf = 7;

  // call the canvas customizer
  PlotCustomizer pcust{1,0,0};  

  // grab the trees
  TFile *f1 = util_pd::ReadRootFile("/w/halla-scshelf2102/sbs/pdbforce/hist/sbs7_prepass2_atimeOff.root");
  // TTree *Tlhloq = (TTree*)flhloq->Get("Tout");
  TCanvas *c2 = (TCanvas*)f1->Get("c2");

  // Accessing each sub-pad
  TPad *p1 = (TPad*)c2->GetPad(1); // Sub-pad 1 of Pad 1
  TPad *p2 = (TPad*)c2->GetPad(2); // Sub-pad 2 of Pad 1  

  // Create a new canvas for the contents of subpad1
  TCanvas *cv = new TCanvas("cv","cv",1000,800);
  cv->Divide(1,2);
  cv->cd(1);
  // Loop through all objects in subpad1
  TIter next(p1->GetListOfPrimitives());
  TObject *obj;
  bool firstHistogram = true;
  std::vector<TH1*> h1s;         // Vector to store histogram pointers
  while ((obj = next())) {
    TObject *objClone = obj->Clone();  // Clone each object
    
    // Determine the draw options
    if (objClone->InheritsFrom("TH1")) {       // For histograms
      TH1* hist = (TH1*)objClone;
      h1s.push_back(hist);        
      if (firstHistogram) {
	objClone->Draw();                  // Draw the first histogram normally
	firstHistogram = false;
      } else {
	objClone->Draw("same");            // Overlay additional histograms
      }
    } else if (objClone->InheritsFrom("TGraph")) {   // For graphs
      objClone->Draw("same");                // Draw graphs on top
    } else {                                   // For other types of objects
	objClone->Draw();                      // Draw normally without "same"
    }
  }
  util_pd::SetAxTitles(h1s[0],"#font[32]{t^{TDC}_{TH} - t^{ADC}_{SH}} (ns)","SH Module ID");
  pcust.AddTitleText("Before Correction",0.32);
  //
  cv->cd(2);  
  // Loop through all objects in subpad2
  TIter next2(p2->GetListOfPrimitives());
  TObject *obj2;
  bool firstHistogram2 = true;
  std::vector<TH1*> h1s2;         // Vector to store histogram pointers
  while ((obj2 = next2())) {
    TObject *objClone2 = obj2->Clone();  // Clone each object
    
    // Determine the draw options
    if (objClone2->InheritsFrom("TH1")) {       // For histograms
      TH1* hist = (TH1*)objClone2;
      h1s2.push_back(hist);        
      if (firstHistogram2) {
	objClone2->Draw();                  // Draw the first histogram normally
	firstHistogram2 = false;
      } else {
	objClone2->Draw("same");            // Overlay additional histograms
      }
    } else if (objClone2->InheritsFrom("TGraph")) {   // For graphs
      objClone2->Draw("same");                // Draw graphs on top
    } else {                                   // For other types of objects
	objClone2->Draw();                      // Draw normally without "same"
    }
  }
  util_pd::SetAxTitles(h1s2[0],"#font[32]{t^{TDC}_{TH} - t^{ADC}_{SH}} (ns)","SH Module ID");
  pcust.AddTitleText("After Correction",0.3);
  //
  pcust.customize_canvas(cv);
  //
  h1s[0]->GetYaxis()->SetTitleSize(0.07);
  h1s[0]->GetYaxis()->SetTitleOffset(0.8);
  h1s2[0]->GetYaxis()->SetTitleSize(0.07);
  h1s2[0]->GetYaxis()->SetTitleOffset(0.8);
  //  
  cv->Update();
  cv->SaveAs(Form("detcalib_shatime_%d.pdf",sbsconf));
}

//____________________________________________________
void PlotPSAtimeRes()
{
  int sbsconf = 9;

  // call the canvas customizer
  PlotCustomizer pcust{1,1};

  TChain *ch = new TChain("T");
  ch->Add("/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS9/LH2/rootfiles/*"); //1369
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-30,30);
  ch->Draw("bb.ps.clus_blk.atime-bb.ps.atimeblk>>h1","bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.065&&bb.ps.clus_blk.id!=bb.ps.idblk&&bb.gem.track.nhits>3&&bb.tr.p[0]>1.2&&abs(e.kine.W2-0.88)<0.6");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","Time Difference (ns)");
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",-2.5,2.5);
  h1->Fit(ffn,"R");
  ffn->SetLineColor(kRed);
  ffn->SetLineWidth(4);  
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);  
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  custom_statbox_effi(s1c);
  s1c->SetY1NDC(0.6);
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_psatime_resolution_%d.pdf",sbsconf));        

}

//____________________________________________________
void PlotSHAtimeRes()
{
  int sbsconf = 9;

  // call the canvas customizer
  PlotCustomizer pcust{1,1};

  TChain *ch = new TChain("T");
  ch->Add("/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS9/LH2/rootfiles/*"); //1369
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-30,30);
  ch->Draw("bb.sh.clus_blk.atime-bb.sh.atimeblk>>h1","bb.ps.e>0.2&&abs(bb.tr.vz[0])<0.065&&bb.sh.clus_blk.id!=bb.sh.idblk&&bb.gem.track.nhits>3&&bb.tr.p[0]>1.2&&abs(e.kine.W2-0.88)<0.6");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","Time Difference (ns)");
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",-2.5,2.5);
  h1->Fit(ffn,"R");
  ffn->SetLineColor(kRed);
  ffn->SetLineWidth(4);  
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);  
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  custom_statbox_effi(s1c);
  s1c->SetY1NDC(0.6);
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_shatime_resolution_%d.pdf",sbsconf));        

}

//____________________________________________________
void PlotBBCALCalibPlots()
{
  int sbsconf = 7;
  
  // call the canvas customizer
  PlotCustomizer pcust;

  TString rfile = sbsconf==7 ? "sbs7-sbs85p-set2_prepass2_bbcal_eng_calib_elcut.root" : "sbs8-sbs70p-set1_prepass2_bbcal_eng_calib_elcut_1.root";
  
  // grab the trees
  TFile *f1 = util_pd::ReadRootFile(Form("/lustre24/expphy/volatile/halla/sbs/pdbforce/bbcalib/hist/pass2/%s",rfile.Data()));
  // TTree *Tlhloq = (TTree*)flhloq->Get("Tout");
  TCanvas *cel7 = (TCanvas*)f1->Get("c7");
  TPad* p7 = (TPad*)cel7->GetPad(1);

  // Accessing each sub-pad
  TPad *subpad1 = (TPad*)p7->GetPad(1); // Sub-pad 1 of Pad 1
  TPad *subpad2 = (TPad*)p7->GetPad(2); // Sub-pad 2 of Pad 1

  // Create a new canvas for the contents of subpad1
  TCanvas *cv = util_pd::TC("cv",1,2);
  cv->cd(1);
  // Loop through all objects in subpad1
  TIter next(subpad1->GetListOfPrimitives());
  TObject *obj;
  bool firstHistogram = true;
  std::vector<TH1*> histos;         // Vector to store histogram pointers
  std::vector<TLine*> lines;
  while ((obj = next())) {
    TObject *objClone = obj->Clone();  // Clone each object

    if (objClone->InheritsFrom("TLine")) {    // Check if the object is a line
        TLine *line = (TLine*)objClone;
        TLine *lineClone = (TLine*)line->Clone();  // Clone the line

        // Shift line positions
        double x1 = lineClone->GetX1() + xNDCoffset;
        double y1 = lineClone->GetY1() + y1NDCoffset;
        double x2 = lineClone->GetX2() + xNDCoffset;
        double y2 = lineClone->GetY2();
        
        lineClone->SetX1(x1);
        lineClone->SetY1(y1);
        lineClone->SetX2(x2);
        lineClone->SetY2(y2);

        lines.push_back(lineClone);      // Store the shifted line for later use
        lineClone->Draw();               // Draw the shifted line on the new canvas
    }
    
    // Determine the draw options
    if (objClone->InheritsFrom("TH1")) {       // For histograms
      TH1* hist = (TH1*)objClone;
      histos.push_back(hist);        
      if (firstHistogram) {
	objClone->Draw();                  // Draw the first histogram normally
	firstHistogram = false;
      } else {
	objClone->Draw("same");            // Overlay additional histograms
      }
    } else if (objClone->InheritsFrom("TGraph")) {   // For graphs
      objClone->Draw("same");                // Draw graphs on top
    } else {                                   // For other types of objects
      if (!objClone->InheritsFrom("TLine")) 
	objClone->Draw();                      // Draw normally without "same"
    }
  }
  //
  cv->cd(2);
  // Loop through all objects in subpad1
  TIter next2(subpad2->GetListOfPrimitives());
  TObject *obj2;
  bool firstHistogram2 = true;
  while ((obj2 = next2())) {
    TObject *objClone = obj2->Clone();  // Clone each object

    // Determine the draw options
    if (objClone->InheritsFrom("TH1")) {       // For histograms
      if (firstHistogram2) {
	objClone->Draw();                  // Draw the first histogram normally
	firstHistogram2 = false;
      } else {
	objClone->Draw("same");            // Overlay additional histograms
      }
    } else if (objClone->InheritsFrom("TGraph")) {   // For graphs
      objClone->Draw("same");                // Draw graphs on top
    } else {                                   // For other types of objects
      objClone->Draw();                      // Draw normally without "same"
    }
  }
  //
  pcust.customize_canvas(cv);
  //
  cv->cd(1);
  histos[0]->SetLineColor(kGray+2);
  histos[0]->SetFillColor(kGray);
  histos[1]->SetLineColor(kBlack);
  histos[1]->SetFillColor(kGray+3);  
  TLegend *l1=new TLegend(0.36,0.76,0.95,0.9);
  l1->SetTextFont(62);
  l1->AddEntry(histos[0],"All","f");
  l1->AddEntry(histos[1],"With 5#sigma-6#sigma #Deltax-#Deltay Cut","f");
  l1->Draw();  
  //
  cv->Update();       
  cv->SaveAs(Form("detcalib_bbcaleng_elcut_%d.pdf",sbsconf));       
}

//____________________________________________________
void PlotSHAtimeIssue2Solved()
{
  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch = new TChain("Tout");
  ch->Add("~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs50p_model1_pass2.root");
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH2F *h1 = new TH2F("h1","",200,-30,30,200,-30,30);
  ch->Draw("atimeHCAL:atimeSH>>h1","WCut&&thpq_p<0.01&&bbfiduCut&&SMCut&&coinTADCCut");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"#font[32]{t^{ADC}_{HCAL}} (ns)","#font[32]{t^{ADC}_{SH}} (ns)");
  //
  pcust.customize_canvas(cv);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_shatime_issue2_solved_%d.pdf",sbsconf));     

}

//____________________________________________________
void PlotSHAtimeIssue2()
{
  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch = new TChain("Tout");
  ch->Add("~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs50p_model1_pass2.root");
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-20,20);
  ch->Draw("atimeSH-1>>h1","WCut&&thpq_p<0.01&&bbfiduCut&&SMCut&&coinTADCCut");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","#font[32]{t^{ADC}_{SH}} (ns)");
  //
  pcust.customize_canvas(cv);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_shatime_issue2_%d.pdf",sbsconf));     

}

//____________________________________________________
void PlotTHHCALCoinRaw()
{
  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch = new TChain("Tout");
  ch->Add("~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs50p_model1_pass2.root");
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,340,400);
  ch->Draw("cltmeanHODO+bbT_trig-atimeHCAL>>h1","eHCAL>0&&WCut&&pCut&&abs(atimeHCAL-atimeSH)<10");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","#font[32]{t^{TDC}_{TH} + t^{TDC}_{BBtrig} - t^{ADC}_{HCAL}} (ns)");
  //
  pcust.customize_canvas(cv);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_thissue2_%d.pdf",sbsconf));     

}

//____________________________________________________
void PlotTHHCALCoin()
{
  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch = new TChain("Tout");
  ch->Add("~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs50p_model1_pass2.root");
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-18,22);
  ch->Draw("cltmeanHODO-atimeHCAL+1>>h1","eHCAL>0&&WCut&&pCut");
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","#font[32]{t^{TDC}_{TH} - t^{ADC}_{HCAL}} (ns)");
  //
  pcust.customize_canvas(cv);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_thissue1_%d.pdf",sbsconf));     

}

//____________________________________________________
void PlotTHRes2()
{
  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust{1,1};
  
  TChain *ch = new TChain("Tout");
  ch->Add("~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs30p_model1_pass2.root");
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h1 = new TH1F("h1","",200,-20,20);
  ch->Draw("cltmeanHODO-0.87>>h1");
  // statbox manipulation
  gStyle->SetOptStat(1110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)h1->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  h1->SetStats(0);
  //
  h1->SetLineColor(kBlack);
  h1->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(h1,"","TH Cluster Mean Time (ns)");
  //
  pcust.customize_canvas(cv);
  //
  s1c->Draw("same");
  custom_statbox_effi(s1c);  
  //
  cv->Update();
  cv->SaveAs(Form("detcalib_thresolution2_%d.pdf",sbsconf));   
}


//_____________________________ ************ _________________________
void detcalib()
{
  // // Draw TH mean time resolution
  // PlotTHRes2();

  // // Draw TH-HCAL ADC issue 1
  // PlotTHHCALCoin();

  // // Draw TH-HCAL ADC issue 2
  // PlotTHHCALCoinRaw();

  // // Draw SH Atime issue 2
  // PlotSHAtimeIssue2();    

  // // Draw SH Atime issue 2 solved
  // PlotSHAtimeIssue2Solved();    
  
  // // Draw BBCAL Calib plots
  // PlotBBCALCalibPlots();

  // // Draw Intrinsic SH resolution
  // PlotSHAtimeRes(); 
  
  // // Draw Intrinsic PS resolution
  // PlotPSAtimeRes();  

  // // Draw SH atime alignment
  // PlotSHAtimeAlign();

  // // Draw PS atime alignment
  // PlotPSAtimeAlign();

  // Draw SH position resolutions
  PlotSHPosRes();

  // // Plot RF time
  // PlotRFTime();

  // // Plot out of time sh events in a cluster
  // PlotOutofTimeSHHits();

  // // Plot HCAL ADC time alignmnt
  // PlotHCALAtimeAlign();

  // // Plot intrinsic time resolution of HCAL
  // PlotHCALtimeRes();

  // // Plot HCAL Energy alignment
  // PlotHCALEngAlign();

  // // Plot HCAL dx resolution
  // PlotHCALdxRes();

  // // Plot HCAL dy resolution
  // PlotHCALdyRes();

  // Plot BBCAL energy calib resolution
  //PlotBBCALEngResSBS9();
  // PlotBBCALEngResSBS8();
}
