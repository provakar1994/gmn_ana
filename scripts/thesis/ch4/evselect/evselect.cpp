#include "gmn_ana.h"

double xNDCoffset = 0.05;
double y1NDCoffset = 0.05;

//____________________________________________________
void PlotthpqCut()
{
  int sbsconf = 14;

  TString xtitle = "#font[32]{#theta_{pq}} (rad)";
  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *cthpqc = new TChain("Tout");
  cthpqc->Add("~/gmn_ana/scripts/pdout/bsf_1p13zoff_qelas_ana_data_sbs14_sbs70p_model2_pass2.root");
  //
  TCanvas *cvthpqc = util_pd::TC("cvthpqc",2,2);
  gStyle->SetLineScalePS(2.5);
  gStyle->SetPalette(kRainbow);
  //
  cvthpqc->cd(1);
  TH1F *hthpqn = new TH1F("hthpqn","",200,0,0.4);
  cthpqc->Draw("thpq_n>>hthpqn","WCut&&bbfiduCut&&fiduCut");
  hthpqn->SetLineColor(kBlack);
  hthpqn->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(hthpqn,"",xtitle);
  pcust.AddTitleText("Neutron Hypothesis",0.65);
  // drawing cut ranges
  double x1NDC = util_pd::GetxNDC(0.02)+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);    
  //
  cvthpqc->cd(2);
  TH1F *hthpqp = new TH1F("hthpqp","",200,0,0.4);
  cthpqc->Draw("thpq_p>>hthpqp","WCut&&bbfiduCut&&fiduCut");
  hthpqp->SetLineColor(kBlack);
  hthpqp->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(hthpqp,"",xtitle);
  pcust.AddTitleText("Proton Hypothesis",0.63);  
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);    
  //
  cvthpqc->cd(3);
  TH2F *hnothpq = new TH2F("hnothpq","",200,-1,1,200,-2,1);
  cthpqc->Draw("dx:dy>>hnothpq","WCut&&bbfiduCut&&fiduCut&&eHCAL>0","colz");
  util_pd::SetAxTitles(hnothpq,"#font[32]{#Deltax} (m)","#font[32]{#Deltay} (m)");
  pcust.AddTitleText("No Cut on #font[32]{#theta_{pq}}",0.48);  
  //
  cvthpqc->cd(4);
  TH2F *hwthpq = new TH2F("hwthpq","",200,-1,1,200,-2,1);
  cthpqc->Draw("dx:dy>>hwthpq","WCut&&bbfiduCut&&fiduCut&&eHCAL>0&&(thpq_n<0.02||thpq_p<0.02)","colz");
  util_pd::SetAxTitles(hwthpq,"#font[32]{#Deltax} (m)","#font[32]{#Deltay} (m)");
  pcust.AddTitleText("With #font[32]{#theta_{pq}} Cut",0.48);  
  //  
  pcust.customize_canvas(cvthpqc);
  //
  TLegend *lthpq=new TLegend(0.58,0.83,0.95,0.9);
  lthpq->SetTextFont(62);
  lthpq->AddEntry(L1,Form("#font[32]{#theta_{pq}} < 0.02 rad"),"l");
  cvthpqc->cd(1); lthpq->Draw();
  cvthpqc->cd(2); lthpq->Draw();
  //
  cvthpqc->Update();
  cvthpqc->SaveAs(Form("evselect_thpqc_%d.pdf",sbsconf));    
}

//____________________________________________________
void PlotOpticsValidityCut()
{
  int nbinsx = 200;
  int nbinsy = 200;
  double hminx1 = -0.5;
  double hmaxx1 = 0.55;
  double hminy1 = -1.8;
  double hmaxy1 = 5;
  double hminx2 = -0.125;
  double hmaxx2 = 0.125;
  double hminy2 = -1.8;
  double hmaxy2 = 5;

  TString xtitle1 = "#font[32]{x_{BB}} (m)";
  TString xtitle2 = "#font[32]{y_{BB}} (m)";  
  TString ytitle = "#font[32]{W^{2}} (GeV^{2})";

  int sbsconf = 14;
    
  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *copval = new TChain("Tout");
  copval->Add("~/gmn_ana/scripts/pdout/1p13zoff_elas_ana_data_sbs14_sbs70p_model1_pass2.root");  
  //
  TString cuts = "thpq_p<0.03";
  //
  TCanvas *cvopval = util_pd::TC("cvopval",1,2);
  gStyle->SetLineScalePS(2.5);
  gStyle->SetPalette(kRainbow);
  cvopval->cd(1);
  //
  TH2F *hw2xbb = new TH2F("hw2xbb","",nbinsx,hminx1,hmaxx1,nbinsy,hminy1,hmaxy1);
  copval->Draw("W2:fpX-0.9*fpTh>>hw2xbb",cuts.Data(),"");
  util_pd::SetAxTitles(hw2xbb,ytitle,xtitle1);
  // drawing cut ranges
  double x1NDC = util_pd::GetxNDC(-0.25)+xNDCoffset;
  double x2NDC = util_pd::GetxNDC(0.3)+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);
  TLine L2;
  L2.SetLineColor(2); L2.SetLineWidth(4); //L2.SetLineStyle(9);
  L2.DrawLineNDC(x2NDC,y1NDC,x2NDC,0.9);  
  //
  cvopval->cd(2);
  //
  TH2F *hw2ybb = new TH2F("hw2ybb","",nbinsx,hminx2,hmaxx2,nbinsy,hminy2,hmaxy2);
  copval->Draw("W2:fpY-0.9*fpPh>>hw2ybb",cuts.Data(),"");
  util_pd::SetAxTitles(hw2ybb,ytitle,xtitle2);
  // drawing cut ranges
  x1NDC = util_pd::GetxNDC(-0.08)+xNDCoffset;
  x2NDC = util_pd::GetxNDC(0.1)+xNDCoffset;
  y1NDC = 0.1+y1NDCoffset;
  TLine *L3 = new TLine();
  L3->SetLineColor(2); L3->SetLineWidth(4); //L3.SetLineStyle(9);
  L3->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);
  TLine L4;
  L4.SetLineColor(2); L4.SetLineWidth(4); //L4.SetLineStyle(9);
  L4.DrawLineNDC(x2NDC,y1NDC,x2NDC,0.9);  
  //  
  pcust.customize_canvas(cvopval);
  //
  //
  cvopval->Update();
  cvopval->SaveAs(Form("evselect_opval_%d.pdf",sbsconf));  
}

//____________________________________________________
void PlotdyCutEffect()
{
  int nbin = 200;
  double hmin = -3.5;
  double hmax = 2;

  int sbsconf = 14;

  // call the canvas customizer
  PlotCustomizer pcust;
  
  TChain *cdyceff = new TChain("Tout");
  cdyceff->Add("~/gmn_ana/scripts/pdout/bsf_1p13zoff_qelas_ana_data_sbs14_sbs70p_model2_pass2.root");  
  //
  TString cuts1 = "W>0.5&&W<1.15&&eHCAL>0.1&&bbfiduCut&&abs(coinT_ADC_c)<1.34*3.5&&fiduCut";
  TString cuts2 = "W>0.5&&W<1.15&&abs(dy)<0.3&&eHCAL>0.1&&bbfiduCut&&abs(coinT_ADC_c)<1.34*3.5&&fiduCut";
  //
  TCanvas *cvdyceff = util_pd::TC("cvdyceff",1,1);
  gStyle->SetLineScalePS(2.5);
  //gPad->SetLogy();
  cvdyceff->cd();
  //
  TH1F *h1 = new TH1F("h1","",nbin,hmin,hmax);
  cdyceff->Draw("dx>>h1",cuts1.Data(),"");
  h1->SetLineColor(kGray+2);
  h1->SetFillColor(kGray);  
  util_pd::SetAxTitles(h1,"","#font[32]{#Deltax} (m)");
  //
  TH1F *hdyceff = new TH1F("hdyceff","",nbin,hmin,hmax);
  cdyceff->Draw("dx>>hdyceff",cuts2.Data(),"same");  
  hdyceff->SetLineColor(kBlack);
  hdyceff->SetFillColor(kGray+3);  
  //
  pcust.customize_canvas(cvdyceff);
  //
  TLegend *ldyceff=new TLegend(0.62,0.76,0.95,0.9);
  ldyceff->SetTextFont(62);
  ldyceff->AddEntry(h1,"No Cut on #font[32]{#Deltay}","f");
  ldyceff->AddEntry(hdyceff,"#left|#font[32]{#Deltay}#right| #leq 0.3 m","f");
  ldyceff->Draw();
  //
  cvdyceff->Update();
  cvdyceff->SaveAs(Form("evselect_dyceff_%d.pdf",sbsconf));    
}

//____________________________________________________
void PlotdxdyCutEffect()
{
  int nbin = 200;
  double hmin = 0;
  double hmax = 2;

  int sbsconf = 14;

  // call the canvas customizer
  PlotCustomizer pcust;
  
  TChain *cdxdyceff = new TChain("Tout");
  cdxdyceff->Add("~/gmn_ana/scripts/pdout/1p13zoff_elas_ana_data_sbs14_sbs70p_model2_pass2.root");  
  //
  TString cuts1 = "bbfiduCut&&eHCAL>0.1&&abs(coinT_ADC_c)<1.3*3.5";
  TString cuts2 = "bbfiduCut&&eHCAL>0.1&&abs(coinT_ADC_c)<1.3*3.5&&pCut";
  //
  TCanvas *cvdxdyceff = util_pd::TC("cvdxdyceff",1,1);
  gStyle->SetLineScalePS(2.5);
  //gPad->SetLogy();
  cvdxdyceff->cd();
  //
  TH1F *h1 = new TH1F("h1","",nbin,hmin,hmax);
  cdxdyceff->Draw("W2>>h1",cuts1.Data(),"");
  h1->SetLineColor(kGray+2);
  h1->SetFillColor(kGray);  
  util_pd::SetAxTitles(h1,"","#font[32]{W^{2}} (GeV^{2})");
  //
  TH1F *hdxdyceff = new TH1F("hdxdyceff","",nbin,hmin,hmax);
  cdxdyceff->Draw("W2>>hdxdyceff",cuts2.Data(),"same");  
  hdxdyceff->SetLineColor(kBlack);
  hdxdyceff->SetFillColor(kGray+3);  
  //
  pcust.customize_canvas(cvdxdyceff);
  //
  TLegend *ldxdyceff=new TLegend(0.1+xNDCoffset,0.76,0.65,0.9);
  ldxdyceff->SetTextFont(62);
  ldxdyceff->AddEntry(h1,"No Cut on #font[32]{#Deltax-#Deltay} Correlation","f");
  ldxdyceff->AddEntry(hdxdyceff,"With 2#sigma #font[32]{#Deltax-#Deltay} Correlation Cut","f");
  ldxdyceff->Draw();
  //
  cvdxdyceff->Update();
  cvdxdyceff->SaveAs(Form("evselect_dxdyceff_%d.pdf",sbsconf));    
}

//____________________________________________________
void PlotdxwW2Cut()
{
  int nbin = 200;
  double hmin = -3.5;
  double hmax = 2.5;

  int sbsconf = 7;

  // call the canvas customizer
  PlotCustomizer pcust;

  TChain *cdxW2 = new TChain("Tout");
  cdxW2->Add("~/gmn_ana/scripts/pdout/bsf_0p925zoff_qelas_ana_data_sbs7_sbs85p_model2_pass2.root");  
  //
  TString cuts_noW2 = "abs(dy+0.05)<0.3&&eHCAL>0.1&&bbfiduCut&&abs(vz)<0.075&&abs(EovP-1)<0.2&&abs(coinT_ADC_c)<1.3*3.5&&fiduCut";
  TString cuts_wW2 = "abs(W-0.8)<0.35&&abs(dy+0.05)<0.3&&eHCAL>0.1&&bbfiduCut&&abs(vz)<0.075&&abs(EovP-1)<0.2&&abs(coinT_ADC_c)<1.3*3.5&&fiduCut";
  //
  TCanvas *cvdxW2 = util_pd::TC("cvdxW2",1,1);
  gStyle->SetLineScalePS(2.5);
  //gPad->SetLogy();
  cvdxW2->cd();
  //
  TH1F *hdxnoW2 = new TH1F("hdxnoW2","",nbin,hmin,hmax);
  cdxW2->Draw("dx>>hdxnoW2",cuts_noW2.Data(),"");
  hdxnoW2->SetLineColor(kGray+2);
  hdxnoW2->SetFillColor(kGray);  
  util_pd::SetAxTitles(hdxnoW2,"","#font[32]{#Deltax} (m)");
  //
  TH1F *hdxW2 = new TH1F("hdxW2","",nbin,hmin,hmax);
  cdxW2->Draw("dx>>hdxW2",cuts_wW2.Data(),"same");  
  hdxW2->SetLineColor(kBlack);
  hdxW2->SetFillColor(kGray+3);  
  //
  pcust.customize_canvas(cvdxW2);
  //
  TLegend *ldxW2=new TLegend(0.62,0.76,0.95,0.9);
  ldxW2->SetTextFont(62);
  ldxW2->AddEntry(hdxnoW2,"No Cut on #font[32]{W^{2}}","f");
  ldxW2->AddEntry(hdxW2,"#left|#font[32]{W^{2}} - 0.88#right| #leq 0.35 GeV^{2}","f");
  ldxW2->Draw();
  //
  cvdxW2->Update();
  cvdxW2->SaveAs(Form("evselect_dxW2_%d.pdf",sbsconf));    
}
   
//____________________________________________________
void PlotdyAll()
{
  int nbin = 200;
  double hmin = -1.5;
  double hmax = 1.5;

  int sbsconf = 14;

  // call the canvas customizer
  PlotCustomizer pcust;  

  TChain *clh70 = new TChain("Tout");
  clh70->Add("~/gmn_ana/scripts/pdout/1p13zoff_elas_ana_data_sbs14_sbs70p_model1_pass2.root");
  TChain *cld70 = new TChain("Tout");
  cld70->Add("~/gmn_ana/scripts/pdout/bsf_1p13zoff_qelas_ana_data_sbs14_sbs70p_model2_pass2.root");
  //
  TString cuts = "WCut&&eHCAL>0.05&&abs(coinT_ADC_c)<3.5*1.3&&ARCut";
  //
  TCanvas *cvdyall = util_pd::TC("cvdyall",1,1);
  gStyle->SetLineScalePS(2.5);
  //gPad->SetLogy();
  cvdyall->cd();  
  //
  TH1F *hlh70 = new TH1F("hlh70","",nbin,hmin,hmax);
  clh70->Draw("dy>>hlh70",cuts.Data(),"");
  //hlh70->Scale(1/5.48);
  //hlh70->Draw("HIST");
  hlh70->SetFillColor(kGreen+2);
  hlh70->SetLineColor(kGreen+3);
  util_pd::SetAxTitles(hlh70,"","#font[32]{#Deltay} (m)");
  TH1F *hld70 = new TH1F("hld70","",nbin,hmin,hmax);
  cld70->Draw("dy>>hld70",cuts.Data(),"HIST same");
  hld70->Scale(1/2.0);
  //hld70->Draw("HIST same");
  hld70->SetFillColorAlpha(kRed,0.55);
  hld70->SetLineColor(kRed+1); 
  //
  pcust.customize_canvas(cvdyall);
  //
  TLegend *ldyall=new TLegend(0.62,0.76,0.95,0.9);
  ldyall->SetTextFont(62);
  ldyall->AddEntry(hlh70,"#font[32]{l}H_{2}, 50% SBS Field","f");
  ldyall->AddEntry(hld70,"#font[32]{l}D_{2}, 50% SBS Field","f");
  ldyall->Draw();
  //
  cvdyall->Update();
  cvdyall->SaveAs(Form("evselect_dyall_%d.pdf",sbsconf));      
}

//____________________________________________________
void PlotdxAll()
{
  int nbin = 200;
  double hmin = -2.5;
  double hmax = 2;

  int sbsconf = 4;

  // call the canvas customizer
  PlotCustomizer pcust;  

  TChain *clh0 = new TChain("Tout");
  clh0->Add("~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs0p_model1_pass2.root");
  TChain *clh50 = new TChain("Tout");
  clh50->Add("~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs50p_model1_pass2.root");
  TChain *cld50 = new TChain("Tout");
  cld50->Add("~/gmn_ana/scripts/pdout/0p65zoff_qelas_ana_data_sbs4_sbs50p_model2_pass2.root");  

  //
  TCanvas *cvdxall = util_pd::TC("cvdxall",1,1);
  gStyle->SetLineScalePS(2.5);
  //gPad->SetLogy();
  cvdxall->cd();  
  //
  TH1F *hlh0 = new TH1F("hlh0","",nbin,hmin,hmax);
  clh0->Draw("dx>>hlh0","WCut&&eHCAL>0&&abs(dy)<0.3&&abs(coinT_ADC_c)<2*1.3&&ARCut","");
  hlh0->Scale(1/5.48);
  hlh0->Draw("HIST");
  hlh0->SetFillColor(kGreen+2);
  hlh0->SetLineColor(kGreen+3);
  util_pd::SetAxTitles(hlh0,"","#font[32]{#Deltax} (m)");
  TH1F *hlh50 = new TH1F("hlh50","",nbin,hmin,hmax);
  clh50->Draw("dx>>hlh50","WCut&&eHCAL>0&&abs(dy)<0.3&&abs(coinT_ADC_c)<2*1.3&&ARCut","same");
  hlh50->SetFillColor(kBlue);
  hlh50->SetLineColor(kBlue+1);  
  TH1F *hld50 = new TH1F("hld50","",nbin,hmin,hmax);
  cld50->Draw("dx>>hld50","WCut&&eHCAL>0&&abs(dy)<0.3&&abs(coinT_ADC_c)<2*1.3&&ARCut","same");
  hld50->SetFillColorAlpha(kRed,0.55);
  hld50->SetLineColor(kRed+1); 
  //
  pcust.customize_canvas(cvdxall);
  //
  TLegend *ldxall=new TLegend(0.62,0.69,0.95,0.9);
  ldxall->SetTextFont(62);
  ldxall->AddEntry(hlh0,"#font[32]{l}H_{2}, SBS OFF","f");
  ldxall->AddEntry(hlh50,"#font[32]{l}H_{2}, 50% SBS Field","f");
  ldxall->AddEntry(hld50,"#font[32]{l}D_{2}, 50% SBS Field","f");
  ldxall->Draw();
  //
  cvdxall->Update();
  cvdxall->SaveAs(Form("evselect_dxall_%d.pdf",sbsconf));    
}

//____________________________________________________
void PlotCoinTime(int sbsconf, TString const & infile, TString const & cuts, double coincutrange)
{
  int nbin = 200;
  double hmin = -30;
  double hmax = 30;

  // call the canvas customizer
  PlotCustomizer pcust;
  
  TChain *ccoin = new TChain("T");
  ccoin->Add(infile.Data());
  //
  TCanvas *cvcoin = util_pd::TC("cvcoin",1,1);
  gStyle->SetLineScalePS(2.5);
  //gPad->SetLogy();
  cvcoin->cd();
  //
  TH1F *hcoin = new TH1F("hcoin","hcoin",nbin,hmin,hmax);
  ccoin->Draw("sbs.hcal.atimeblk-bb.sh.atimeblk-0.77>>hcoin",cuts.Data());
  //
  hcoin->SetLineColor(kBlack);
  hcoin->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(hcoin,"","#font[32]{t_{coin}} (ns)");
  hcoin->Draw();
  // drawing cut ranges
  double x1NDC = util_pd::GetxNDC(-coincutrange)+xNDCoffset;
  double x2NDC = util_pd::GetxNDC(coincutrange)+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);
  TLine L2;
  L2.SetLineColor(2); L2.SetLineWidth(4); //L2.SetLineStyle(9);
  L2.DrawLineNDC(x2NDC,y1NDC,x2NDC,0.9);  
  //
  pcust.customize_canvas(cvcoin);
  //
  TLegend *lcoin=new TLegend(0.64,0.83,0.95,0.9);
  lcoin->SetTextFont(62);
  lcoin->AddEntry(L1,Form("#left|#font[32]{t_{coin}}#right| #leq %.0f ns",coincutrange),"l");
  lcoin->Draw();
  
  cvcoin->Update();
  cvcoin->SaveAs(Form("evselect_coin_%d.pdf",sbsconf));  
}

//____________________________________________________
void PlotHCALClEng(int sbsconf, TString const & infile, TString const & cuts, double thresh)
{
  int nbin = 200;
  double hmin = 0;
  double hmax = 1.2;

  // call the canvas customizer
  PlotCustomizer pcust;
  
  TChain *chcl = new TChain("Tout");
  chcl->Add(infile.Data());
  //
  TCanvas *cvhcl = util_pd::TC("cvhcl",1,1);
  gStyle->SetLineScalePS(2.5);
  cvhcl->cd();
  //  
  TH1F *hhcl = new TH1F("hhcl","hhcl",nbin,hmin,hmax);
  chcl->Draw("eHCAL>>hhcl",cuts.Data());
  //
  hhcl->SetLineColor(kBlack);
  hhcl->SetFillColor(kGreen-8);
  util_pd::SetAxTitles(hhcl,"","#font[32]{E_{HCAL}} (GeV)");
  hhcl->Draw();
  //
  // drawing cut ranges
  double x1NDC = util_pd::GetxNDC(thresh)+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);  
  //
  pcust.customize_canvas(cvhcl);
  //
  TLegend *lhcl=new TLegend(0.64,0.83,0.95,0.9);
  lhcl->SetTextFont(62);
  // lhcl->AddEntry(h_mip,"MIP peak","f");
  // lhcl->AddEntry(fmip,"Fit","l");
  lhcl->AddEntry(L1,Form("#font[32]{E_{HCAL}} > %.2f GeV",thresh),"l");
  lhcl->Draw();
  //
  cvhcl->Update();
  cvhcl->SaveAs(Form("evselect_hcalcl_%d.pdf",sbsconf));
}

//____________________________________________________
void PlotFiduCut(int sbsconf, TString const &infile_ld, double sbs_kick, std::vector<double> const &SM_w, TString const &cut)
{
  int nbinsx = 200;
  int nbinsy = 200;
  double hminx = -1.5;
  double hmaxx = 1.5;
  double hminy = -3;
  double hmaxy = 2.5;

  TString xtitle = "#font[32]{y^{exp}_{HCAL}} (m)";
  TString ytitlen = "#font[32]{x^{exp}_{HCAL}} (m)";
  //TString ytitlep = Form("#font[32]{x^{exp}_{HCAL}} - %.2f (m)",sbs_kick);
  TString ytitlep = "#font[32]{x^{exp}_{HCAL} - #deltax_{SBS}} (m)";
  
  // call the canvas customizer
  PlotCustomizer pcust{0};

  // HCAL Area, AR, and SM
  std::vector<double> hcal_AR0 = cut::hcal_active_area_data(0,0,2); 
  std::vector<double> hcal_AR = cut::hcal_active_area_data(1,1,2); //1 block from all sides 
  std::vector<double> hcal_SM0 = cut::hcal_safety_margin(0,0,0,hcal_AR); // effectively no safety margin
  std::vector<double> hcal_SM = cut::hcal_safety_margin(SM_w[0],SM_w[1],SM_w[2],hcal_AR);

  // reading ROOT files as df
  ROOT::EnableImplicitMT();
  ROOT::RDataFrame rdf("Tout",infile_ld);
  auto rdf_filtered = rdf.Filter(cut.Data());
  std::string xHCAL_exp_shifted = "xHCAL_exp-p_def"; //"xHCAL_exp-" + std::to_string(sbs_kick);
  rdf_filtered = rdf_filtered.Define("xHCAL_exp_shifted",xHCAL_exp_shifted.c_str());

  // No fiducial cut (entire region)
  // create and plot the stuff related to the first canvas
  TCanvas *cfid1 = util_pd::TC("cfid1",1,2);
  gStyle->SetLineScalePS(2.5);
  gStyle->SetPalette(kRainbow);
  cfid1->cd(1);
  TH2F *h2_n = (TH2F*)rdf_filtered.Histo2D({"h2_n","",nbinsx,hminx,hmaxx,nbinsy,hminy,hmaxy,},"yHCAL_exp","xHCAL_exp")->Clone();
  h2_n->Draw("colz"); 
  util_pd::SetAxTitles(h2_n,ytitlen,xtitle);
  util_pd::DrawArea(hcal_AR0,kGreen+2,4,1);
  pcust.AddTitleText("Neutron Envelope",0.56);
  cfid1->cd(2);
  TH2F *h2_p = (TH2F*)rdf_filtered.Histo2D({"h2_p","",nbinsx,hminx,hmaxx,nbinsy,hminy,hmaxy,},"yHCAL_exp","xHCAL_exp_shifted")->Clone();
  h2_p->Draw("colz");
  util_pd::SetAxTitles(h2_p,ytitlep,xtitle);
  util_pd::DrawArea(hcal_AR0,kGreen+2,4,1);
  //
  pcust.AddTitleText("Proton Envelope",0.53);
  pcust.customize_canvas(cfid1);
  //pcust.AddText(Form("Events %d",h2_n->GetEntries()),0.1+xNDCoffset,0.81,0.53,0.9);
  // --

  // Fiducial cut with no safety margin
  // forming the fiducial cut w/ no SM
  auto fiduCut0 = [&](double x,double y,double xExp,double yExp) {
    return cut::inHCAL_activeA(x,y,hcal_AR) && cut::inHCAL_safety_margin("LD2",xExp,yExp,sbs_kick,hcal_SM0);
  };
  // create and plot the stuff related to the second canvas
  TCanvas *cfid2 = util_pd::TC("cfid2",1,2);
  gStyle->SetLineScalePS(2.5);
  //gStyle->SetPalette(kRainbow);
  cfid2->cd(1);
  TH2F *h2_n_wAR = (TH2F*)rdf_filtered.Filter(fiduCut0,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
    .Histo2D({"h2_n_wAR","",nbinsx,hminx,hmaxx,nbinsy,hminy,hmaxy,},"yHCAL_exp","xHCAL_exp")->Clone();
  h2_n_wAR->Draw("colz");
  util_pd::SetAxTitles(h2_n_wAR,ytitlen,xtitle);
  util_pd::DrawArea(hcal_AR0,kGreen+2,4,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  pcust.AddTitleText("Neutron Envelope",0.56);
  cfid2->cd(2);
  TH2F *h2_p_wAR = (TH2F*)rdf_filtered.Filter(fiduCut0,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
    .Histo2D({"h2_p_wAR","",nbinsx,hminx,hmaxx,nbinsy,hminy,hmaxy,},"yHCAL_exp","xHCAL_exp_shifted")->Clone();
  h2_p_wAR->Draw("colz");
  util_pd::SetAxTitles(h2_p_wAR,ytitlep,xtitle);
  util_pd::DrawArea(hcal_AR0,kGreen+2,4,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  pcust.AddTitleText("Proton Envelope",0.53);
  pcust.customize_canvas(cfid2);

  // Fiducial cut with no safety margin
  // forming the fiducial cut w/ no SM
  auto fiduCut1 = [&](double x,double y,double xExp,double yExp) {
    return cut::inHCAL_activeA(x,y,hcal_AR) && cut::inHCAL_safety_margin("LD2",xExp,yExp,sbs_kick,hcal_SM);
  };
  // create and plot the stuff related to the second canvas
  TCanvas *cfid3 = util_pd::TC("cfid3",1,2);
  gStyle->SetLineScalePS(2.5);  
  //gStyle->SetPalette(kRainbow);
  cfid3->cd(1);
  TH2F *h2_n_wARSM = (TH2F*)rdf_filtered.Filter(fiduCut1,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
    .Histo2D({"h2_n_wARSM","",nbinsx,hminx,hmaxx,nbinsy,hminy,hmaxy,},"yHCAL_exp","xHCAL_exp")->Clone();
  h2_n_wARSM->Draw("colz");
  util_pd::SetAxTitles(h2_n_wARSM,ytitlen,xtitle);
  util_pd::DrawArea(hcal_AR0,kGreen+2,4,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  util_pd::DrawArea(hcal_SM,4,4,9);
  pcust.AddTitleText("Neutron Envelope",0.56);
  cfid3->cd(2);
  TH2F *h2_p_wARSM = (TH2F*)rdf_filtered.Filter(fiduCut1,{"xHCAL","yHCAL","xHCAL_exp","yHCAL_exp"})
    .Histo2D({"h2_p_wARSM","",nbinsx,hminx,hmaxx,nbinsy,hminy,hmaxy,},"yHCAL_exp","xHCAL_exp_shifted")->Clone();
  h2_p_wARSM->Draw("colz");
  util_pd::SetAxTitles(h2_p_wARSM,ytitlep,xtitle);
  util_pd::DrawArea(hcal_AR0,kGreen+2,4,1);
  util_pd::DrawArea(hcal_AR,2,4,9);
  util_pd::DrawArea(hcal_SM,4,4,9);
  pcust.AddTitleText("Proton Envelope",0.53);
  pcust.customize_canvas(cfid3);

  // saving the canvases
  cfid1->Update();  cfid2->Update();  cfid3->Update();
  cfid1->SaveAs(Form("evselect_fidu1_%d.pdf",sbsconf));
  cfid2->SaveAs(Form("evselect_fidu2_%d.pdf",sbsconf));
  cfid3->SaveAs(Form("evselect_fidu3_%d.pdf",sbsconf));
}

//____________________________________________________
void PlotW2Cut(int sbsconf_loq, int sbsconf_hiq, TString const &infile_lh_loq, TString const &infile_lh_hiq, TString const &infile_ld, TString const &cut)
{
  int nbins = 160;
  double hmin = 0;
  double hmax = 2;
  double scale_lh_lonhiq = 11;
  double scale_lh_ld = 1.5;
  double height_lh_lonhiq = 25000;
  double height_lh_ld = 32000;

  // call the canvas customizer
  PlotCustomizer pcust;

  // grab the trees
  TFile *flhloq = util_pd::ReadRootFile(infile_lh_loq);
  TTree *Tlhloq = (TTree*)flhloq->Get("Tout");
  //
  TFile *flhhiq = util_pd::ReadRootFile(infile_lh_hiq);
  TTree *Tlhhiq = (TTree*)flhhiq->Get("Tout");
  //
  TFile *fld = util_pd::ReadRootFile(infile_ld);
  TTree *Tld = (TTree*)fld->Get("Tout");

  // create and plot the stuff related to the first canvas
  TCanvas *cw21 = util_pd::TC("cw21",1,1);
  gStyle->SetLineScalePS(2.5);
  cw21->cd();
  // histo
  TH1F *h_lhloq = new TH1F("h_lhloq","",nbins,hmin,hmax);
  Tlhloq->Draw("W2>>h_lhloq",cut.Data());
  TH1F *h_lhhiq = new TH1F("h_lhhiq","",nbins,hmin,hmax);
  Tlhhiq->Draw("W2>>h_lhhiq",cut.Data());
  h_lhhiq->Scale(scale_lh_lonhiq);

  // plotting histos
  h_lhloq->SetLineColor(kBlack);
  h_lhloq->Draw();
  h_lhhiq->SetLineColor(kBlack);
  h_lhhiq->Draw("HIST same");
  h_lhloq->Draw("same");

  // drawing cut ranges
  double x1NDC = util_pd::GetxNDC(pow(constant::Mp,2))+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(2); L1->SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);  

  // customize
  h_lhloq->GetYaxis()->SetRangeUser(0.,height_lh_lonhiq);
  h_lhloq->SetFillColor(kGray+3);
  h_lhhiq->SetFillColor(kGray);
  pcust.customize_canvas(cw21);
  // further customization
  util_pd::SetAxTitles(h_lhloq,"","#font[32]{W^{2}} (GeV^{2})");

  //
  TLegend *lw21=new TLegend(0.1+xNDCoffset,0.69,0.86,0.9);
  lw21->SetTextFont(62);
  lw21->AddEntry(h_lhhiq,"Inclusive #font[32]{H(e,e'p), Q^{2} = 7.4} (GeV/c)^{2}","f");
  lw21->AddEntry(h_lhloq,"Inclusive #font[32]{H(e,e'p), Q^{2} = 3} (GeV/c)^{2}","f");
  lw21->AddEntry(L1,Form("#font[32]{W^{2} = M_{p}^{2}#approx} %.2f GeV^{2}",pow(constant::Mp,2)),"l");
  lw21->Draw();  
    
  // create and plot the stuff related to the first canvas
  TCanvas *cw22 = util_pd::TC("cw22",1,1);
  gStyle->SetLineScalePS(2.5);
  cw22->cd();
  // histo
  TH1F *h_ld = new TH1F("h_ld","",nbins,hmin,hmax);
  Tld->Draw("W2>>h_ld",cut.Data());
  h_ld->Scale(scale_lh_ld);

  // plotting histos
  h_lhloq->SetLineColor(kBlack);
  h_lhloq->Draw();
  h_ld->SetLineColor(kBlack);
  h_ld->Draw("same HIST");
  h_lhloq->Draw("same");
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);  
  
  // customize
  h_lhloq->GetYaxis()->SetRangeUser(0.,height_lh_ld);
  h_lhloq->SetFillColor(kGray+3);
  h_ld->SetFillColor(kGray);
  pcust.customize_canvas(cw22);
  // further customization
  util_pd::SetAxTitles(h_lhloq,"","#font[32]{W^{2}} (GeV^{2})");

  //
  TLegend *lw22=new TLegend(0.1+xNDCoffset,0.69,0.56,0.9);
  lw22->SetTextFont(62);
  lw22->AddEntry(h_ld,"Inclusive #font[32]{D(e,e'N)}","f");
  lw22->AddEntry(h_lhloq,"Inclusive #font[32]{H(e,e'p)}","f");
  lw22->AddEntry(L1,Form("#font[32]{W^{2} = M_{p}^{2}#approx} %.2f GeV^{2}",pow(constant::Mp,2)),"l");
  lw22->Draw();

  // saving the canvases
  cw21->Update(); cw22->Update();
  cw21->SaveAs(Form("evselect_w2cut1_%d_%d.pdf",sbsconf_loq,sbsconf_hiq));
  cw22->SaveAs(Form("evselect_w2cut2_%d_%d.pdf",sbsconf_loq,sbsconf_hiq));  
}

//____________________________________________________
void PlotGRINCHCut(int sbsconf, TString const &infile, TString const &cut)
{
  TString wgrinchcl = cut + "&&(bb.grinch_tdc.clus.trackindex==0&&bb.grinch_tdc.clus.size>2)";
  TString wogrinchcl = cut + "&&(bb.grinch_tdc.clus.trackindex!=0||bb.grinch_tdc.clus.size==1)";
  int nbins = 160;
  double hmin = 0;
  double hmax = 1.4;

  // call the canvas customizer
  PlotCustomizer pcust;

  TCanvas *cgrinch = util_pd::TC("cgrinch",1,1);
  gStyle->SetLineScalePS(2.5);  
  cgrinch->cd();
  // grab cryo target related tree
  TFile *f = util_pd::ReadRootFile(infile);
  TTree *T = (TTree*)f->Get("T");

  // histo
  TH1F *h_pse = new TH1F("h_pse","",nbins,hmin,hmax);
  T->Draw("bb.ps.e>>h_pse",cut.Data());
  TH1F *h_wgrinch = new TH1F("h_wgrinch","",nbins,hmin,hmax);
  T->Draw("bb.ps.e>>h_wgrinch",wgrinchcl.Data());
  TH1F *h_wogrinch = new TH1F("h_wogrinch","",nbins,hmin,hmax);
  T->Draw("bb.ps.e>>h_wogrinch",wogrinchcl.Data());

  // plotting histos
  h_pse->SetLineColor(kBlack);
  h_pse->Draw();
  h_wgrinch->Draw("same");
  h_wogrinch->Draw("same");

  // customize
  pcust.customize_canvas(cgrinch);
  // further customization
  util_pd::SetAxTitles(h_pse,"","#font[32]{E_{PS}} (GeV)");
  h_wgrinch->SetLineColor(kGreen+2);
  h_wgrinch->SetFillColorAlpha(kGreen+2,0.7);
  h_wogrinch->SetLineColor(kRed);
  h_wogrinch->SetFillColorAlpha(kRed,0.7);

  //
  TLegend *lgrinch=new TLegend(0.58,0.69,0.95,0.9);
  lgrinch->SetTextFont(62);
  lgrinch->AddEntry(h_pse,"All","l");
  lgrinch->AddEntry(h_wgrinch,"With Good GRINCH Cluster","f");
  lgrinch->AddEntry(h_wogrinch,"W/O Good GRINCH Cluster","f");
  lgrinch->Draw();

  cgrinch->SaveAs(Form("evselect_grinch_%d.pdf",sbsconf));
}

//____________________________________________________
void PlotPSEngCut(int sbsconf, TString const &infile, TString const &cut, double const thresh)
{
  // h_pse limits
  // NOTE: make sure (hmax-hmin)/nbins is a multiple of 0.0125
  //       otherwise, the superposition of h_mip will not work
  int nbins = 224;
  double hmin = 0;
  double hmax = 1.4;
  
  // call the canvas customizer
  PlotCustomizer pcust;

  TCanvas *cpse = util_pd::TC("cpse",1,1);
  gStyle->SetLineScalePS(2.5);
  cpse->cd();
  // grab cryo target related tree
  TFile *f = util_pd::ReadRootFile(infile);
  TTree *T = (TTree*)f->Get("T");

  // histo
  TH1F *h_pse = new TH1F("h_pse","",nbins,hmin,hmax);
  T->Draw("bb.ps.e>>h_pse",cut.Data());
  double binw = h_pse->GetBinWidth(1);
  TH1F *h_mip = new TH1F("h_mip","",0.125/binw,0.025,0.15);
  T->Draw("bb.ps.e>>h_mip",cut.Data());

  // Fit
  std::vector<double> hmip_fitR{0.025,0.15,1.7,1.4};
  TF1 *fmip = fit::fit_1gs_nbg(hmip_fitR,h_mip);
  auto fpm = fit::GetFitParamANDError(fmip);
  fmip->SetLineColor(kBlue);
  fmip->SetLineWidth(4);
  
  // plotting histos
  h_pse->SetLineColor(kBlack);
  h_pse->Draw();
  h_mip->Draw("same");
  h_pse->Draw("same");
  fmip->Draw("same");
  
  // drawing cut ranges
  double x1NDC = util_pd::GetxNDC(thresh)+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);
  
  // customize
  TString tmip = Form("Mean = %.1f (%d) MeV, #sigma = %.1f (%d) MeV",
		      fpm[1].first*1e3,util_pd::GetSigDigit(fpm[1].second*1e3,1),
		      fpm[2].first*1e3,util_pd::GetSigDigit(fpm[2].second*1e3,1));
  pcust.AddTitleText(tmip,0.85);
  pcust.customize_canvas(cpse);
  // further customization
  util_pd::SetAxTitles(h_pse,"","#font[32]{E_{PS}} (GeV)");
  h_mip->SetLineWidth(1);
  h_mip->SetLineColor(kGreen+2);
  h_mip->SetFillColorAlpha(kGreen+2,0.7);

  //
  TLegend *lpse=new TLegend(0.64,0.69,0.95,0.9);
  lpse->SetTextFont(62);
  lpse->AddEntry(h_mip,"MIP peak","f");
  lpse->AddEntry(fmip,"Fit","l");
  lpse->AddEntry(L1,Form("#font[32]{E_{PS}} > %.2f GeV",thresh),"l");
  lpse->Draw();
  
  cpse->SaveAs(Form("evselect_pse_%d.pdf",sbsconf));
}

//____________________________________________________
void PlotOpticsValCut(int sbsconf, TString const &infile, TString const &elcut, std::vector<double> const &cutRange) {
  // call the canvas customizer
  PlotCustomizer pcust;

  // dispersive
  TCanvas *copvX = util_pd::TC("copvX",1,1);
  copvX->cd(); gStyle->SetPalette(kRainBow);
  // grab cryo target related tree
  TFile *fX = util_pd::ReadRootFile(infile);
  TTree *TX = (TTree*)fX->Get("Tout");
  // histo
  TH2F *h_opvX = new TH2F("h_opvX","",200,-0.15,0.15,200,-0.5,0.7);
  TX->Draw("fpX-0.9*fpTh:vz>>h_opvX",elcut.Data(),"colz");
  h_opvX->Draw("colz");
  // Add title text
  TString topvX = "Dispersive Direction";
  pcust.AddTitleText(topvX,0.54);
  // drawing cut ranges
  double x1NDC = 0.1+xNDCoffset;
  double x2NDC = 0.9-0.01;
  double y1C = cutRange[0] - cutRange[1];
  double y2C = cutRange[0] + cutRange[1];
  double y1NDC = util_pd::GetyNDC(y1C)+y1NDCoffset-0.005; // Ad-hoc offset of 0.005 is applied!
  double y2NDC = util_pd::GetyNDC(y2C)+0.005;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x2NDC,y1NDC);
  TLine L2;
  L2.SetLineColor(2); L2.SetLineWidth(4); //L2.SetLineStyle(9);
  L2.DrawLineNDC(x1NDC,y2NDC,x2NDC,y2NDC);
  // customize
  pcust.customize_canvas(copvX);
  // further customization
  util_pd::SetAxTitles(h_opvX,"#font[32]{x_{BB}} (m)","#font[32]{v_{z}} (m)");
  // plot legend
  TLegend *lopvX = new TLegend(0.1+xNDCoffset,0.81,0.53,0.9);
  lopvX->SetTextFont(62);
  lopvX->AddEntry(L1,Form("|#font[32]{x_{BB}}-%.3f| #leq %.3f",cutRange[0],cutRange[1]),"l");
  //lopvX->SetBorderSize(0);
  lopvX->Draw();
  // **** ----

  // non-dispersive
  TCanvas *copvY = util_pd::TC("copvY",1,1);
  copvY->cd(); gStyle->SetPalette(kRainBow);
  // grab cryo target related tree
  TFile *fY = util_pd::ReadRootFile(infile);
  TTree *TY = (TTree*)fY->Get("Tout");
  // histo
  TH2F *h_opvY = new TH2F("h_opvY","",200,-0.15,0.15,200,-0.15,0.15);
  TY->Draw("fpY-0.9*fpPh:vz>>h_opvY",elcut.Data(),"colz");
  h_opvY->Draw("colz");
  // Add title text
  TString topvY = "Non-Dispersive Direction";
  pcust.AddTitleText(topvY,0.62);
  // customize
  pcust.customize_canvas(copvY);
  // further customization
  util_pd::SetAxTitles(h_opvY,"#font[32]{y_{BB}} (m)","#font[32]{v_{z}} (m)");
  // **** ----

  // saving the canvases
  copvX->Update(); copvY->Update();
  copvX->SaveAs(Form("evselect_opvX_%d.pdf",sbsconf));
  copvY->SaveAs(Form("evselect_opvY_%d.pdf",sbsconf));
}

//____________________________________________________
void PlotvzCut(int sbsconf, TString &infile_ld2, TString const &elcut, double vzcutrange) {
  // call the canvas customizer
  PlotCustomizer pcust;

  // Cryo stuff
  TCanvas *cvzC = util_pd::TC("cvzC",1,1);
  gStyle->SetLineScalePS(2.5);
  cvzC->cd();
  // grab cryo target related tree
  TFile *fC = util_pd::ReadRootFile(infile_ld2);
  TTree *TC = (TTree*)fC->Get("Tout");
  // histo
  TH1F *h_vzC = new TH1F("h_vzC","",200,-0.15,0.15);
  TC->Draw("vz>>h_vzC",elcut.Data());
  h_vzC->Draw();

  // Dummy stuff
  TString infile_dummy = infile_ld2.ReplaceAll("qelas","dummy_qelas");
  // //
  // TCanvas *cvzD = util_pd::TC("cvzD",1,1);
  // cvzD->cd();  
  // grab dummy target related tree
  TFile *fD = util_pd::ReadRootFile(infile_dummy);
  TTree *TD = (TTree*)fD->Get("Tout");
  // histo
  TH1F *h_vzD = new TH1F("h_vzD","",200,-0.15,0.15);
  TD->Draw("vz>>h_vzD",elcut.Data());

  // drawing histos
  h_vzC->SetLineColor(kBlack);
  h_vzC->Draw();
  h_vzD->Draw("same");
  // drawing cut ranges
  double x1NDC = util_pd::GetxNDC(-vzcutrange)+xNDCoffset;
  double x2NDC = util_pd::GetxNDC(vzcutrange)+xNDCoffset;
  double y1NDC = 0.1+y1NDCoffset;
  TLine *L1 = new TLine();
  L1->SetLineColor(2); L1->SetLineWidth(4); //L1.SetLineStyle(9);
  L1->DrawLineNDC(x1NDC,y1NDC,x1NDC,0.9);
  TLine L2;
  L2.SetLineColor(2); L2.SetLineWidth(4); //L2.SetLineStyle(9);
  L2.DrawLineNDC(x2NDC,y1NDC,x2NDC,0.9);

  // customize
  pcust.customize_canvas(cvzC);
  // further customization
  util_pd::SetAxTitles(h_vzC,"","#font[32]{v_{z}} (m)");
  h_vzD->SetLineWidth(1);
  h_vzD->SetLineColor(kGreen+2);
  h_vzD->SetFillColorAlpha(kGreen+2,0.7);

  // plot legend
  TLegend *lvz = new TLegend(0.42,0.27,0.67,0.52);
  lvz->SetTextFont(62);
  lvz->AddEntry(h_vzC,"#font[32]{l}D_{2}","l");
  lvz->AddEntry(h_vzD,"Dummy","f");
  lvz->AddEntry(L1,Form("#left|#font[32]{v_{z}}#right| #leq %.3f m",vzcutrange),"l");
  lvz->SetBorderSize(0);
  lvz->Draw();

  cvzC->SaveAs(Form("evselect_vz_%d.pdf",sbsconf));
}

//_____________________________ ************ _________________________
void evselect()
{
#if 0
#endif
  // Drawing vz cuts
  int sbsconf = 4;
  TString infile_ld2 = "~/gmn_ana/scripts/pdout/0p65zoff_qelas_ana_data_sbs4_sbs50p_model2_pass2.root";
  TString elcut = "W2>0.25&&W2<1.2&&abs(dy)<0.3&&bbfiduCut&&eHCAL>0&&abs(coinT_ADC_c)<1.38*3.5&&fiduCut";
  double vzcutrange = 0.075; 
  PlotvzCut(sbsconf,infile_ld2,elcut,vzcutrange);

  // // Drawing optics validity cut
  // int sbsconf1 = 14;
  // TString infile1 = "~/gmn_ana/scripts/pdout/widevz_1p13zoff_qelas_ana_data_sbs14_sbs70p_model2_pass2.root";
  // TString elcut1 = "WCut&&(pCut||nCut)";
  // std::vector<double> cutRange{0.035,0.335};
  // PlotOpticsValCut(sbsconf1,infile1,elcut1,cutRange);

  // // Drawing ps energy cut
  // int sbsconf2 = 9;
  // TString infile2 = "~/gmn_ana/scripts/pdout/e1209019_fullreplay_13697_stream0_segALL_ALL.root";
  // TString cut2 = "fEvtHdr.fTrigBits==1&&bb.tr.n>0&&abs(bb.tr.vz[0])<.075&&bb.gem.track.nhits>3&&bb.tr.p[0]>1.2";
  // double psethresh = 0.2;
  // PlotPSEngCut(sbsconf2,infile2,cut2,psethresh);

  // // Drawing grinch cluster cut
  // int sbsconf3 = 9;
  // TString infile3 = "~/gmn_ana/scripts/pdout/e1209019_fullreplay_13697_stream0_segALL_ALL.root";
  // TString cut3 = "fEvtHdr.fTrigBits==1&&bb.tr.n>0&&abs(bb.tr.vz[0])<.075&&bb.gem.track.nhits>3&&bb.tr.p[0]>1.2";
  // PlotGRINCHCut(sbsconf3,infile3,cut3);

  // // Drawing W2 cut
  // int sbsconf_loq = 4;
  // int sbsconf_hiq = 14;
  // TString infile_lh_loq = "~/gmn_ana/scripts/pdout/0p65zoff_elas_ana_data_sbs4_sbs0p_model2_pass2.root";
  // TString infile_lh_hiq = "~/gmn_ana/scripts/pdout/1p13zoff_elas_ana_data_sbs14_sbs70p_model2_pass2.root";
  // TString infile_ld = "~/gmn_ana/scripts/pdout/0p65zoff_qelas_ana_data_sbs4_sbs0p_model2_pass2.root";
  // TString cut4 = "bbfiduCut";
  // PlotW2Cut(sbsconf_loq,sbsconf_hiq,infile_lh_loq,infile_lh_hiq,infile_ld,cut4);

  // // Drawing Fiducial cut
  // int sbsconf5 = 4;
  // TString infile_ld5 = "~/gmn_ana/scripts/pdout/bsf_0p65zoff_qelas_ana_data_sbs4_sbs30p_model2_pass2.root";
  // double sbs_kick = 0.68; //m
  // std::vector<double> SM_w{0.19,0.19,0.24};
  // TString cut5 = "WCut&&bbfiduCut&&coinTADCCut";
  // PlotFiduCut(sbsconf5,infile_ld5,sbs_kick,SM_w,cut5);

  // // Drawing HCAL cluster energy
  // int sbsconf6 = 14;
  // TString infile6 = "~/gmn_ana/scripts/pdout/bsf_1p13zoff_qelas_ana_data_sbs14_sbs70p_model2_pass2.root";
  // TString cuts6 = "eHCAL>0&&WCut&&(pCut||nCut)&&coinTADCCut&&fiduCut";
  // double threshold6 = 0.05;
  // PlotHCALClEng(sbsconf6,infile6,cuts6,threshold6);

  // // Drawing HCAL-SH coin time
  // int sbsconf7 = 9;
  // TString infile7 = "~/gmn_ana/scripts/pdout/e1209019_fullreplay_13697_stream0_segALL_ALL.root";
  // TString cuts7 = "fEvtHdr.fTrigBits==1&&bb.tr.n>0&&abs(bb.tr.vz[0])<.08&&bb.gem.track.nhits>3&&bb.ps.e>0.2&&bb.tr.p[0]>1.2&&sbs.hcal.e>0.05&&abs(e.kine.W2-0.88)<1";
  // double coincutrange = 5;
  // PlotCoinTime(sbsconf7,infile7,cuts7,coincutrange);

  // // Draw dx w/ various target and field
  // PlotdxAll();

  // // Draw dy w/ various target and field
  // PlotdyAll();

  // // Draw dx w/ and w/ W2 cut 
  // PlotdxwW2Cut();

  // // Draw dx-dy cut effect
  // PlotdxdyCutEffect();

  // // Draw dy cut effect
  // PlotdyCutEffect();

  // // Draw optics validity cut
  // PlotOpticsValidityCut();

  // // Draw thepq cuts
  // PlotthpqCut();
}
