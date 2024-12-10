#include "gmn_ana.h"

double xNDCoffset = 0.05;
double y1NDCoffset = 0.05;

void custom_statbox_effi(TPaveStats *st) {
  st->SetBit(TH1::kNoStats);
  st->SetX1NDC(0.64); st->SetY1NDC(0.69); st->SetX2NDC(0.95); st->SetY2NDC(0.9);
}

void custom_effi_histo(TH1 *h) {
  h->SetLineWidth(0);
  h->SetMarkerStyle(20);
  h->SetLineColor(kWhite);

  // Set bin errors to zero
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    h->SetBinError(i, 0);  // Remove error for each bin
  }  
}


//____________________________________________________
void PlotHCALTrueNDE()
{
  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust;

  TFile *fp = util_pd::ReadRootFile("/w/halla-scshelf2102/sbs/pdbforce/siout/hcal_det_effi_p_gun_sbs4_sbs0p.root");
  TH1F *hp= (TH1F*)fp->Get("h_effi");
  TFile *fn = util_pd::ReadRootFile("/w/halla-scshelf2102/sbs/pdbforce/siout/hcal_det_effi_n_gun_sbs4_sbs0p.root");
  TH1F *hn= (TH1F*)fn->Get("h_effi");

  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  hp->Draw("p");
  util_pd::SetAxTitles(hp,"Efficiency (%)","Nucleon Momentum (GeV/c)");
  custom_effi_histo(hp);
  hn->Draw("same p");
  custom_effi_histo(hn);
  //
  pcust.customize_canvas(cv);
  // plot legend
  TLegend *l1 = new TLegend(0.42,0.32,0.67,0.52);
  l1->SetTextFont(62);
  l1->AddEntry(hp,"Proton","p");
  l1->AddEntry(hn,"Neutron","p");
  //l1->SetBorderSize(0);
  l1->Draw();  
  //
  cv->Update();
  cv->SaveAs(Form("hcalnde_true_1.pdf"));
}

//____________________________________________________
void PlotpNnEngDistforpGun()
{
  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust;

  TFile *fp = util_pd::ReadRootFile("/w/halla-scshelf2102/sbs/pdbforce/siout/hcal_det_effi_p_gun_sbs4_sbs0p.root");
  TH1F *hp= (TH1F*)fp->Get("h_eHCAL_22");
  TFile *fn = util_pd::ReadRootFile("/w/halla-scshelf2102/sbs/pdbforce/siout/hcal_det_effi_n_gun_sbs4_sbs0p.root");
  TH1F *hn= (TH1F*)fn->Get("h_eHCAL_22");

  //
  TCanvas *cv = util_pd::TC("cv",1,2);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd(1);
  hp->Draw();
  //
  hp->SetLineColor(kBlack);
  hp->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(hp,"","#font[32]{E_{HCAL}} (GeV)");
  pcust.AddTitleText("Proton",0.315);  
  // Fitting
  TF1 *ffn = new TF1("ffn","gaus",0.13,0.32);
  hp->Fit(ffn,"R");
  ffn->SetLineColor(kBlue);
  ffn->SetLineWidth(4);
  ffn->Draw("same");
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s1 = (TPaveStats*)hp->FindObject("stats");
  TPaveStats *s1c = (TPaveStats*)s1->Clone("s1c");  // Clone the stats box
  hp->SetStats(0);
  // drawing cut ranges
  double x1NDC = util_pd::GetxNDC(ffn->GetParameter(1)/4.0)+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);    
  //
  cv->cd(2);
  hn->Draw();
  //
  hn->SetLineColor(kBlack);
  hn->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(hn,"","#font[32]{E_{HCAL}} (GeV)");
  pcust.AddTitleText("Neutron",0.34);    
  // Fitting
  TF1 *ffn1 = new TF1("ffn1","gaus",0.13,0.32);
  hn->Fit(ffn1,"R");
  ffn1->SetLineColor(kBlue);
  ffn1->SetLineWidth(4);
  ffn1->Draw("same");
  // statbox manipulation  
  gStyle->SetOptStat("e");
  gStyle->SetOptFit(01110);
  cv->Update();
  TPaveStats *s2 = (TPaveStats*)hn->FindObject("stats");
  TPaveStats *s2c = (TPaveStats*)s2->Clone("s2c");  // Clone the stats box
  hn->SetStats(0);
  // drawing cut ranges
  double x1NDC1 = util_pd::GetxNDC(ffn1->GetParameter(1)/4.0)+xNDCoffset;
  double y1NDC1 = 0.1+y1NDCoffset;
  TLine *L2 = new TLine();
  L2->SetLineColor(2); L2->SetLineWidth(4); //L1.SetLineStyle(9);
  L2->DrawLineNDC(x1NDC1,y1NDC1,x1NDC1,0.9);    
  //
  pcust.customize_canvas(cv);
  //
  cv->cd(1);
  s1c->Draw("same");
  custom_statbox_effi(s1c);
  s1c->SetX1NDC(0.52);
  s1c->SetY1NDC(0.5);
  //
  TLegend *l1=new TLegend(0.55,0.36,0.95,0.495);
  l1->SetTextFont(62);
  l1->AddEntry(ffn,"Fit","l");
  l1->AddEntry(L1,Form("#font[32]{E_{HCAL}} > %0.3f GeV",ffn->GetParameter(1)/4.0),"l");
  l1->Draw();  
  //
  cv->cd(2);
  s2c->Draw("same");
  custom_statbox_effi(s2c);
  s2c->SetX1NDC(0.52);
  s2c->SetY1NDC(0.5);
  //
  TLegend *l2=new TLegend(0.55,0.36,0.95,0.495);
  l2->SetTextFont(62);
  l2->AddEntry(ffn1,"Fit","l");
  l2->AddEntry(L2,Form("#font[32]{E_{HCAL}} > %0.3f GeV",ffn1->GetParameter(1)/4.0),"l");
  l2->Draw();    
  //
  hp->SetLineWidth(1);
  hn->SetLineWidth(1);
  //
  cv->Update();
  cv->SaveAs(Form("hcalnde_true_2.pdf"));    
}


//_____________________________ ************ _________________________
void hcalnde()
{
  
  // // Draw HCAL eng distributions for p and n guns
  // PlotpNnEngDistforpGun();

  // Draw HCAL true NDE
  PlotHCALTrueNDE();
}
