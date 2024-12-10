#include "gmn_ana.h"

double xNDCoffset = 0.05;
double y1NDCoffset = 0.05;

void custom_statbox_effi(TPaveStats *st) {
  st->SetBit(TH1::kNoStats);
  st->SetX1NDC(0.64); st->SetY1NDC(0.69); st->SetX2NDC(0.95); st->SetY2NDC(0.9);
}

//____________________________________________________
void PlotdxInelasticShapeSBS7()
{
  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch7 = new TChain("Tout");
  ch7->Add("~/gmn_ana/scripts/simulation/siout/0p925zoff_0p815sf_inel_qelas_ana_g4sbs_sbs7_sbs85p_model2.root");  
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  //gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h7_1 = new TH1F("h7_1","",200,-4,3);
  ch7->Draw("dx>>h7_1","weight*(W2>0.8&&trP>1.2&&eHCAL>0)","");
  util_pd::SetAxTitles(h7_1,"","#font[32]{#Deltax} (m)");
  //
  TH1F *h7_2 = new TH1F("h7_2","",200,-4,3);
  ch7->Draw("dx>>h7_2","weight*(W2>0.8&&trP>1.2&&eHCAL>0&&mc_fnucl==0)","same");
  //
  TH1F *h7_3 = new TH1F("h7_3","",200,-4,3);
  ch7->Draw("dx>>h7_3","weight*(W2>0.8&&trP>1.2&&eHCAL>0&&mc_fnucl==1)","same");  
  //
  pcust.customize_canvas(cv);
  //
  h7_1->SetLineColor(kBlack);
  h7_1->SetLineWidth(1);
  h7_1->SetMarkerColor(kBlack);
  h7_1->SetMarkerStyle(27);
  h7_2->SetLineColor(kBlue);
  h7_2->SetLineWidth(1);
  h7_2->SetMarkerColor(kBlue);
  h7_2->SetMarkerStyle(27);
  h7_3->SetLineColor(kRed);
  h7_3->SetLineWidth(1);  
  h7_3->SetMarkerColor(kRed);
  h7_3->SetMarkerStyle(27);  
  //
  TLegend *l4=new TLegend(0.72,0.70,0.95,0.9);
  l4->SetTextFont(62);
  l4->AddEntry(h7_1,"All","ep");
  l4->AddEntry(h7_2,"Neutron","ep");
  l4->AddEntry(h7_3,"Proton","ep");
  l4->Draw();  
  //
  cv->Update();
  cv->SaveAs(Form("mc_inelbgshape_7.pdf"));     
  
}

//____________________________________________________
void PlotdxInelasticShapeSBS4()
{
  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *ch4 = new TChain("Tout");
  ch4->Add("~/gmn_ana/scripts/simulation/siout/0p39sf_inel_qelas_ana_g4sbs_sbs4_sbs30p_model2.root");   
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  //gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *h4_1 = new TH1F("h4_1","",200,-4,3);
  ch4->Draw("dx>>h4_1","weight*(W2>0.8&&trP>1.2&&eHCAL>0)","");
  util_pd::SetAxTitles(h4_1,"","#font[32]{#Deltax} (m)");
  //
  TH1F *h4_2 = new TH1F("h4_2","",200,-4,3);
  ch4->Draw("dx>>h4_2","weight*(W2>0.8&&trP>1.2&&eHCAL>0&&mc_fnucl==0)","same");
  //
  TH1F *h4_3 = new TH1F("h4_3","",200,-4,3);
  ch4->Draw("dx>>h4_3","weight*(W2>0.8&&trP>1.2&&eHCAL>0&&mc_fnucl==1)","same"); 
  //
  pcust.customize_canvas(cv);
  //
  h4_1->SetLineColor(kBlack);
  h4_1->SetLineWidth(1);
  h4_1->SetMarkerColor(kBlack);
  h4_1->SetMarkerStyle(27);
  h4_2->SetLineColor(kBlue);
  h4_2->SetLineWidth(1);
  h4_2->SetMarkerColor(kBlue);
  h4_2->SetMarkerStyle(27);
  h4_3->SetLineColor(kRed);
  h4_3->SetLineWidth(1);  
  h4_3->SetMarkerColor(kRed);
  h4_3->SetMarkerStyle(27);
  //
  TLegend *l4=new TLegend(0.72,0.70,0.95,0.9);
  l4->SetTextFont(62);
  l4->AddEntry(h4_1,"All","ep");
  l4->AddEntry(h4_2,"Neutron","ep");
  l4->AddEntry(h4_3,"Proton","ep");
  l4->Draw();  
  //
  cv->Update();
  cv->SaveAs(Form("mc_inelbgshape_4.pdf"));     
  
}

//____________________________________________________
void PlotHCALSFComp()
{
  int sbsconf = 14;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *chd = new TChain("Tout");
  chd->Add("~/gmn_ana/scripts/pdout/bsf_1p13zoff_qelas_ana_data_sbs14_sbs70p_model2_pass2.root");
  TChain *chs = new TChain("Tout");
  chs->Add("~/gmn_ana/scripts/simulation/siout/bsf_1p13zoff_0p674sf_qelas_ana_simc_sbs14_sbs70p_model2.root");  
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  //gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *hd = new TH1F("hd","",200,0,0.5);
  chd->Draw("eHCAL/nu>>hd","eHCAL>0&&WCut&&(thpq_p<0.01||thpq_n<0.01)");
  util_pd::SetAxTitles(hd,"","#font[32]{f_{cal} = E_{HCAL}/E^{Kin}_{N}}");
  //
  TH1F *hs = new TH1F("hs","",200,0,0.5);
  chs->Draw("eHCAL/nu-0.012>>hs","eHCAL>0&&WCut&&(thpq_p<0.01||thpq_n<0.01)");
  hs->Scale(1/4.);
  //
  hd->SetLineColor(kBlack);
  hd->SetFillColor(kBlack);
  hd->Draw();
  hs->SetLineColor(kGray+2);
  hs->SetFillColorAlpha(kGray+2,0.7);
  hs->Draw("same HIST");
  //
  pcust.customize_canvas(cv);
  //
  TLegend *l1=new TLegend(0.6,0.75,0.95,0.9);
  l1->SetTextFont(62);
  l1->AddEntry(hd,"Data","f");
  l1->AddEntry(hs,"MC","f");
  l1->Draw();  
  //
  cv->Update();
  cv->SaveAs(Form("mc_hcalsfcomp_%d.pdf",sbsconf));     
}

//____________________________________________________
void PlotHCALEngComp()
{
  int sbsconf = 14;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *chd = new TChain("Tout");
  chd->Add("~/gmn_ana/scripts/pdout/bsf_1p13zoff_qelas_ana_data_sbs14_sbs70p_model2_pass2.root");
  TChain *chs = new TChain("Tout");
  chs->Add("~/gmn_ana/scripts/simulation/siout/bsf_1p13zoff_0p674sf_qelas_ana_simc_sbs14_sbs70p_model2.root");  
  //
  TCanvas *cv = util_pd::TC("cv",1,1);
  //gStyle->SetPalette(kRainbow);
  gStyle->SetLineScalePS(2.5);
  //
  cv->cd();
  TH1F *hd = new TH1F("hd","",200,0,1.8);
  chd->Draw("eHCAL>>hd","eHCAL>0&&WCut&&(thpq_p<0.01||thpq_n<0.01)");
  util_pd::SetAxTitles(hd,"","#font[32]{E_{HCAL}} (GeV)");
  //
  TH1F *hs = new TH1F("hs","",200,0,1.8);
  chs->Draw("eHCAL-0.0432>>hs","eHCAL>0&&WCut&&(thpq_p<0.01||thpq_n<0.01)");
  hs->Scale(1/4.);
  //
  hd->SetLineColor(kBlack);
  hd->SetFillColor(kBlack);
  hd->Draw();
  hs->SetLineColor(kGray+2);
  hs->SetFillColorAlpha(kGray+2,0.7);
  hs->Draw("same HIST");
  //
  pcust.customize_canvas(cv);
  //
  TLegend *l1=new TLegend(0.6,0.75,0.95,0.9);
  l1->SetTextFont(62);
  l1->AddEntry(hd,"Data","f");
  l1->AddEntry(hs,"MC","f");
  l1->Draw();  
  //
  cv->Update();
  cv->SaveAs(Form("mc_ehcalcomp_%d.pdf",sbsconf));     
}

//_____________________________ ************ _________________________
void datamc()
{
  // // Plot eHCAL comp for data/MC
  // PlotHCALEngComp();

  // // Plot eHCAL comp for data/MC
  // PlotHCALSFComp();

  // Plot inelastic bg shapes from MC
  PlotdxInelasticShapeSBS4();  
  //PlotdxInelasticShapeSBS7();  
}
