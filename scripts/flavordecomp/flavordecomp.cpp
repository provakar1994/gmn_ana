#include <iostream>

#include "gmn_ana.h"


double const ku = 1.67;
double const kd = -2.03;

//_______________________________________
double ErrPropAplusB(double const AErr, double const BErr) {
  return sqrt(AErr*AErr + BErr*BErr);
}

//_______________________________________
double ErrPropAovB(double const A, double const AErr, double const B, double const BErr) {
  double R = A/B;
  return fabs(R) * sqrt( pow(AErr/A,2.) + pow(BErr/B,2.) );
}

//_______________________________________
std::vector<double> CalcF1(double Q2, double const GE, double const GEErr, double const GM, double const GMErr, std::string const & ntype) {

  double tau = ntype.compare("p")==0 ? kine::tau(Q2,"p") : kine::tau(Q2,"n");

  double numer = tau*GM + GE;
  double numerErr = ErrPropAplusB(tau*GMErr,GEErr);
  double denom = 1. + tau;
  double denomErr = 0.; //ErrPropAplusB(1.,tau); 
  
  return {numer/denom, ErrPropAovB(numer,numerErr,denom,denomErr)};
}

//_______________________________________
std::vector<double> CalcF2(double Q2, double const GE, double const GEErr, double const GM, double const GMErr, std::string const & ntype) {

  double tau = ntype.compare("p")==0 ? kine::tau(Q2,"p") : kine::tau(Q2,"n");

  double numer = GM - GE;
  double numerErr = ErrPropAplusB(GMErr,GEErr);
  double denom = 1. + tau;
  double denomErr = 0.; //ErrPropAplusB(1.,tau);
  
  return {numer/denom, ErrPropAovB(numer,numerErr,denom,denomErr)};
}

//_______________________________________
std::vector<double> CalcF1Quark(double const F1p, double const F1pErr,  double const F1n, double const F1nErr, std::string const & qtype) {

  double F1u = 2.*F1p + F1n;
  double F1uErr = ErrPropAplusB(2.*F1pErr,F1nErr);
  double F1d = 2.*F1n + F1p;
  double F1dErr = ErrPropAplusB(2.*F1nErr,F1pErr);

  std::vector<double> result;
  if (qtype.compare("u")==0) result = {F1u,F1uErr};
  else result = {F1d,F1dErr}; 
    
  return result; 
}

//_______________________________________
std::vector<double> CalcF2Quark(double const F2p, double const F2pErr, double const F2n, double const F2nErr, std::string const & qtype) {

  double F2u = 2.*F2p + F2n;
  double F2uErr = ErrPropAplusB(2.*F2pErr,F2nErr);  
  double F2d = 2.*F2n + F2p;
  double F2dErr = ErrPropAplusB(2.*F2nErr,F2pErr);
  
  std::vector<double> result;
  if (qtype.compare("u")==0) result = {F2u,F2uErr};
  else result = {F2d,F2dErr}; 
    
  return result;   
}

//_______________________________________
void GetFQuark(double Q2, double GMn, double GMnErr, int verbose, std::vector<double> &FQ) {

  // EMFF fits
  Ye2017 yefit;
  // Kelly2004 kellyfit;
  // Seamus20XX seamusfit;
  // Christy2022 christyfit;
  // Arrington2007 arfit;
  // Galster1971 galfit;

  // True
  std::vector<double> temp = yefit.GetFFwErr(G_t::kGEp,Q2);
  double GEp_ye = temp[0], GEp_yeErr = temp[1];
  temp = yefit.GetFFwErr(G_t::kGMp,Q2);
  double GMp_ye = temp[0], GMp_yeErr = temp[1];
  temp = yefit.GetFFwErr(G_t::kGEn,Q2);
  double GEn_ye = temp[0], GEn_yeErr = temp[1];
  temp = yefit.GetFFwErr(G_t::kGMn,Q2);
  double GMn_ye = GMn>0 ? temp[0] : GMn;
  double GMn_yeErr = GMn>0 ? temp[1] : GMnErr;  

  // double GEp_ye = yefit.GetFF(G_t::kGEp,Q2);
  // double GMp_ye = yefit.GetFF(G_t::kGMp,Q2);
  // double GEn_ye = yefit.GetFF(G_t::kGEn,Q2);
  // double GMn_ye = GMn>0 ? yefit.GetFF(G_t::kGMn,Q2) : GMn;

  temp = CalcF1(Q2,GEp_ye,GEp_yeErr,GMp_ye,GMp_yeErr,"p");
  double F1p = temp[0], F1pErr = temp[1];
  temp = CalcF2(Q2,GEp_ye,GEp_yeErr,GMp_ye,GMp_yeErr,"p");
  double F2p = temp[0], F2pErr = temp[1];
  temp = CalcF1(Q2,GEn_ye,GEn_yeErr,GMn_ye,GMn_yeErr,"n");
  double F1n = temp[0], F1nErr = temp[1];
  temp = CalcF2(Q2,GEn_ye,GEn_yeErr,GMn_ye,GMn_yeErr,"n");
  double F2n = temp[0], F2nErr = temp[1];

  temp = CalcF1Quark(F1p,F1pErr,F1n,F1nErr,"u");
  double F1u = temp[0], F1uErr = temp[1];
  temp = CalcF1Quark(F1p,F1pErr,F1n,F1nErr,"d");  
  double F1d = temp[0], F1dErr = temp[1];
  temp = CalcF2Quark(F2p,F2pErr,F2n,F2nErr,"u");
  double F2u = temp[0], F2uErr = temp[1];
  temp = CalcF2Quark(F2p,F2pErr,F2n,F2nErr,"d");
  double F2d = temp[0], F2dErr = temp[1];

  if (verbose>0) {
    std::cout << Form("GEp: %f, GMp: %f, GEn: %f, GMn: %f\n",GEp_ye,GMp_ye,GEn_ye,GMn_ye);
    std::cout << Form("F1p: %f, F2p: %f, F1n: %f, F2n: %f\n",F1p,F2p,F1n,F2n);
    std::cout << Form("F1u: %f, F2u: %f, F1d: %f, F2d: %f\n",F1u,F2u,F1d,F2d);
    std::cout << Form("Q4*F1u: %f, Q4*F2u/ku: %f, Q4*F1d*2.5: %f, Q4*F2d*0.75/kd: %f\n",Q2*Q2*F1u,Q2*Q2*F2u/ku,Q2*Q2*F1d*2.5,Q2*Q2*F2d*0.75/kd);
  }

  FQ = {F1u,F1uErr,F2u,F2uErr,F1d,F1dErr,F2d,F2dErr};
}

//_______________________________________
double Get_Q4F1u_cont(double *x, double *par) {
  /* Return a continuous curve for Q4*F1u using Ye 2017 fit */
  double Q2 = x[0];
  std::vector<double> FQ;

  GetFQuark(Q2,100,0,0,FQ);

  return Q2*Q2*FQ[0];
}

//_______________________________________
double Get_Q4F1d2p5_cont(double *x, double *par) {
  /* Return a continuous curve for Q4*F1d*2.5 using Ye 2017 fit */
  double Q2 = x[0];
  std::vector<double> FQ;

  GetFQuark(Q2,100,0,0,FQ);

  return Q2*Q2*FQ[4]*2.5;
}

//_______________________________________
double Get_Q4F2uovku_cont(double *x, double *par) {
  /* Return a continuous curve for Q4*F2u/ku using Ye 2017 fit */
  double Q2 = x[0];
  std::vector<double> FQ;

  GetFQuark(Q2,100,0,0,FQ);

  return Q2*Q2*FQ[2]/ku;
}

//_______________________________________
double Get_Q4F2d0p75ovkd_cont(double *x, double *par) {
  /* Return a continuous curve for Q4*F2d*0.75/kd using Ye 2017 fit */
  double Q2 = x[0];
  std::vector<double> FQ;

  GetFQuark(Q2,100,0,0,FQ);

  return Q2*Q2*FQ[6]*0.75/kd;
}

//_______________________________________
void customize_hframe(TH1 *h) {
  h->SetStats(0);
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
  // h->GetXaxis()->SetLabelSize(0.04);
  // h->GetXaxis()->SetTitleSize(0.05);
  // // h->GetXaxis()->SetTitleOffset(0.9);
  // h->GetYaxis()->SetLabelSize(0.04);
  // h->GetYaxis()->SetTitleSize(0.05);  
  h->GetYaxis()->SetNdivisions(7);
      
  h->Draw();
  h->GetXaxis()->CenterTitle();
  h->GetYaxis()->CenterTitle();
}

//_______________________________________
void customize_gfit(TGraphErrors *g, std::string qtype) {
  g->SetLineWidth(2);
  g->SetLineStyle(1);
  g->SetFillStyle( 1001 );
  if (qtype.compare("u")==0) {
    g->SetLineColor(1);
    g->SetFillColorAlpha(1,0.25);
  } else {
    g->SetLineColor(2);    
    g->SetFillColorAlpha(2,0.25);
  }
}

//_______________________________________
void customize_sbsgmn(TGraphErrors *g, std::string qtype) {
  g->SetMarkerStyle(20);
  if (qtype.compare("u")==0) {
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
  } else {
    g->SetMarkerColor(kRed);
    g->SetLineColor(kRed);
  }
}

//_______________________________________
int flavordecomp() {

  double Q2min_ye = 1.e-5;
  double Q2max_ye = 15.0;
  int npoints_ye = 2000;
  double Q2step_ye = (Q2max_ye-Q2min_ye)/double(npoints_ye);

  double Q2ye[npoints_ye+1];
  double Q2yeErr[npoints_ye+1];  
  
  double F1u[npoints_ye+1];
  double F1uErr[npoints_ye+1];  
  double F2u[npoints_ye+1];
  double F2uErr[npoints_ye+1];  
  double F1d[npoints_ye+1];
  double F1dErr[npoints_ye+1];  
  double F2d[npoints_ye+1];
  double F2dErr[npoints_ye+1];

  // Reading Tyler's GMn table
  LookUpTableReader reader;
  std::string filename = "GMn_lookup.csv";
  reader.readCSV(filename);

  double F1u_tyler[npoints_ye+1];
  double F1uErr_tyler[npoints_ye+1];  
  double F2u_tyler[npoints_ye+1];
  double F2uErr_tyler[npoints_ye+1];  
  double F1d_tyler[npoints_ye+1];
  double F1dErr_tyler[npoints_ye+1];  
  double F2d_tyler[npoints_ye+1];
  double F2dErr_tyler[npoints_ye+1];        

  for( int i=0; i<=npoints_ye; i++ ){
    double Q2i = Q2min_ye + i*Q2step_ye;
    Q2ye[i] = Q2i; Q2yeErr[i] = 0.;
    
    std::vector<double> FQ;
    GetFQuark(Q2i,100,0,0,FQ);

    F1u[i] = Q2i*Q2i*FQ[0]; F1uErr[i] = Q2i*Q2i*FQ[1]; 
    F2u[i] = Q2i*Q2i*FQ[2]/ku; F2uErr[i] = Q2i*Q2i*FQ[3]/ku;
    F1d[i] = Q2i*Q2i*FQ[4]*2.5; F1dErr[i] = Q2i*Q2i*FQ[5]*2.5;     
    F2d[i] = Q2i*Q2i*FQ[6]*0.75/kd; F2dErr[i] = Q2i*Q2i*FQ[7]*0.75/kd;

    double GMn_tyler = reader.GetClosestValueByKey(Q2i,0)*EMFFFits::GetGDip(Q2i)*constant::mun;
    double GMn_tyler_err = reader.GetClosestValueByKey(Q2i,1)*EMFFFits::GetGDip(Q2i)*constant::mun;     
    GetFQuark(Q2i,GMn_tyler,GMn_tyler_err,0,FQ);    
    
    F1u_tyler[i] = Q2i*Q2i*FQ[0]; F1uErr_tyler[i] = Q2i*Q2i*FQ[1]; 
    F2u_tyler[i] = Q2i*Q2i*FQ[2]/ku; F2uErr_tyler[i] = Q2i*Q2i*FQ[3]/ku;
    F1d_tyler[i] = Q2i*Q2i*FQ[4]*2.5; F1dErr_tyler[i] = Q2i*Q2i*FQ[5]*2.5;     
    F2d_tyler[i] = Q2i*Q2i*FQ[6]*0.75/kd; F2dErr_tyler[i] = Q2i*Q2i*FQ[7]*0.75/kd;     
    
    // F1u[i] = FQ[0]; F1uErr[i] = FQ[1]; 
    // F2u[i] = FQ[2]; F2uErr[i] = FQ[3];
    // F1d[i] = FQ[4]; F1dErr[i] = FQ[5];     
    // F2d[i] = FQ[6]; F2dErr[i] = FQ[7]; 
  }  

  // Using only Ye 2018
  TGraphErrors *F1uTG = new TGraphErrors( npoints_ye+1, Q2ye, F1u, Q2yeErr, F1uErr);
  customize_gfit(F1uTG,"u");
  TGraphErrors *F1dTG = new TGraphErrors( npoints_ye+1, Q2ye, F1d, Q2yeErr, F1dErr);
  customize_gfit(F1dTG,"d");
  TGraphErrors *F2uTG = new TGraphErrors( npoints_ye+1, Q2ye, F2u, Q2yeErr, F2uErr);
  customize_gfit(F2uTG,"u");
  TGraphErrors *F2dTG = new TGraphErrors( npoints_ye+1, Q2ye, F2d, Q2yeErr, F2dErr);
  customize_gfit(F2dTG,"d");
  
  // Using Tyler GMn + Ye 2018
  TGraphErrors *F1uTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, F1u_tyler, Q2yeErr, F1uErr_tyler);
  customize_gfit(F1uTG_tyler,"u");
  TGraphErrors *F1dTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, F1d_tyler, Q2yeErr, F1dErr_tyler);
  customize_gfit(F1dTG_tyler,"d");
  TGraphErrors *F2uTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, F2u_tyler, Q2yeErr, F2uErr_tyler);
  customize_gfit(F2uTG_tyler,"u");
  TGraphErrors *F2dTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, F2d_tyler, Q2yeErr, F2dErr_tyler);
  customize_gfit(F2dTG_tyler,"d");  

  // Using GMn data
  //std::vector<double> Q2v{3.0,4.5,7.4,9.9,13.5};
  std::vector<double> Q2v{2.9890,4.4880,7.4640,9.8340,13.4650};
  std::vector<double> SBSGMnovMuGD{0.9774,0.9763,0.9071,0.8473,0.7582};
  std::vector<double> SBSGMnErrovMuGD{0.0145,0.0164,0.0174,0.0245,0.0226};

  size_t npoints_sbsgmn = Q2v.size();
  double Q2_sbsgmn[npoints_sbsgmn];
  double Q2Err_sbsgmn[npoints_sbsgmn];    
  double F1u_sbsgmn[npoints_sbsgmn];
  double F1uErr_sbsgmn[npoints_sbsgmn];  
  double F2u_sbsgmn[npoints_sbsgmn];
  double F2uErr_sbsgmn[npoints_sbsgmn];  
  double F1d_sbsgmn[npoints_sbsgmn];
  double F1dErr_sbsgmn[npoints_sbsgmn];  
  double F2d_sbsgmn[npoints_sbsgmn];
  double F2dErr_sbsgmn[npoints_sbsgmn];   
  
  for (size_t i=0; i<npoints_sbsgmn; i++) {
    Q2_sbsgmn[i] = Q2v[i];
    Q2Err_sbsgmn[i] = 0.;
    
    std::vector<double> FQ;
    GetFQuark(Q2v[i],SBSGMnovMuGD[i]*EMFFFits::GetGDip(Q2v[i])*constant::mun,SBSGMnErrovMuGD[i]*EMFFFits::GetGDip(Q2v[i])*constant::mun,0,FQ);

    F1u_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[0]; F1uErr_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[1]; 
    F2u_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[2]/ku; F2uErr_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[3]/ku;
    F1d_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[4]*2.5; F1dErr_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[5]*2.5;     
    F2d_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[6]*0.75/kd; F2dErr_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[7]*0.75/kd;     
  }

  TGraphErrors *F1uTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F1u_sbsgmn, Q2Err_sbsgmn, F1uErr_sbsgmn);
  customize_sbsgmn(F1uTG_sbsgmn,"u");
  TGraphErrors *F1dTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F1d_sbsgmn, Q2Err_sbsgmn, F1dErr_sbsgmn);
  customize_sbsgmn(F1dTG_sbsgmn,"d");
  TGraphErrors *F2uTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F2u_sbsgmn, Q2Err_sbsgmn, F2uErr_sbsgmn);
  customize_sbsgmn(F2uTG_sbsgmn,"u");
  TGraphErrors *F2dTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F2d_sbsgmn, Q2Err_sbsgmn, F2dErr_sbsgmn);  
  customize_sbsgmn(F2dTG_sbsgmn,"d");
  
  //**************
  // Plot F1 with SBS-GMn + Ye 2018 and Ye 2018
  TCanvas *c1 = util_pd::TC("c1",1,1);
  c1->cd(); c1->SetGridy();

  TH2D *hframe_f1 = new TH2D("hframe_f1",";Q^{2} (GeV/c)^{2};Q^{4}F_{1}^{q}",500,0,15,500,-0.6,1.6);
  customize_hframe(hframe_f1);
  
  F1uTG->Draw("C3 SAME");
  F1dTG->Draw("C3 SAME");
  F1uTG_sbsgmn->Draw("P SAME");
  F1dTG_sbsgmn->Draw("P SAME");  

  TLegend *l1 = new TLegend(0.1,0.1,0.49,0.3);
  l1->SetTextFont(42);
  l1->AddEntry(F1uTG_sbsgmn,"u quark (SBSGMn+Ye 2018)","ep");
  l1->AddEntry(F1uTG,"u quark (Ye 2018)","lf");
  l1->AddEntry(F1dTG_sbsgmn,"d quark x 2.5 (SBSGMn+Ye 2018)","ep");
  l1->AddEntry(F1dTG,"d quark x 2.5 (Ye 2018)","lf");
  l1->Draw();
  
  //**************
  // Plot F2 with SBS-GMn + Ye 2018 and Ye 2018  
  TCanvas *c2 = util_pd::TC("c2",1,1);
  c2->cd(); c2->SetGridy();

  TH2D *hframe_f2 = new TH2D("hframe_f2",";Q^{2} (GeV/c)^{2};#kappa_{q}^{-1}Q^{4}F_{2}^{q}",500,0,15,500,-0.02,0.3);
  customize_hframe(hframe_f2);
  
  F2uTG->Draw("C3 SAME");
  F2dTG->Draw("C3 SAME");    
  F2uTG_sbsgmn->Draw("P SAME");
  F2dTG_sbsgmn->Draw("P SAME");  

  TLegend *l2 = new TLegend(0.14,0.1,0.54,0.3);
  l2->SetTextFont(42);
  l2->AddEntry(F2uTG_sbsgmn,"u quark (SBSGMn+Ye 2018)","ep");
  l2->AddEntry(F2uTG,"u quark (Ye 2018)","lf");
  l2->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (SBSGMn+Ye 2018)","ep");
  l2->AddEntry(F2dTG,"d quark x 0.75 (Ye 2018)","lf");
  l2->Draw();

  //**************
  // Plot F1 with SBS-GMn + Ye 2018 and Tyler + Ye 2018
  TCanvas *c3 = util_pd::TC("c3",1,1);
  c3->cd(); c3->SetGridy();

  TH2D *hframe_f3 = new TH2D("hframe_f3",";Q^{2} (GeV/c)^{2};Q^{4}F_{1}^{q}",500,0,15,500,-0.6,1.6);
  customize_hframe(hframe_f3);
  
  F1uTG_tyler->Draw("C3 SAME");
  F1dTG_tyler->Draw("C3 SAME");
  F1uTG_sbsgmn->Draw("P SAME");
  F1dTG_sbsgmn->Draw("P SAME");  

  TLegend *l3 = new TLegend(0.1,0.1,0.49,0.3);
  l3->SetTextFont(42);
  l3->AddEntry(F1uTG_sbsgmn,"u quark (SBSGMn+Ye 2018)","ep");
  l3->AddEntry(F1uTG_tyler,"u quark (Tyler+Ye 2018)","lf");
  l3->AddEntry(F1dTG_sbsgmn,"d quark x 2.5 (SBSGMn+Ye 2018)","ep");
  l3->AddEntry(F1dTG_tyler,"d quark x 2.5 (Tyler+Ye 2018)","lf");
  l3->Draw();
  
  //**************
  // Plot F2 with SBS-GMn + Ye 2018 and Tyler + Ye 2018  
  TCanvas *c4 = util_pd::TC("c4",1,1);
  c4->cd(); c4->SetGridy();

  TH2D *hframe_f4 = new TH2D("hframe_f4",";Q^{2} (GeV/c)^{2};#kappa_{q}^{-1}Q^{4}F_{2}^{q}",500,0,15,500,-0.02,0.3);
  customize_hframe(hframe_f4);
  
  F2uTG_tyler->Draw("C3 SAME");
  F2dTG_tyler->Draw("C3 SAME");    
  F2uTG_sbsgmn->Draw("P SAME");
  F2dTG_sbsgmn->Draw("P SAME");  

  TLegend *l4 = new TLegend(0.14,0.1,0.54,0.3);
  l4->SetTextFont(42);
  l4->AddEntry(F2uTG_sbsgmn,"u quark (SBSGMn+Ye 2018)","ep");
  l4->AddEntry(F2uTG_tyler,"u quark (Tyler+Ye 2018)","lf");
  l4->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (SBSGMn+Ye 2018)","ep");
  l4->AddEntry(F2dTG_tyler,"d quark x 0.75 (Tyler+Ye 2018)","lf");
  l4->Draw();  
  
  // TF1 *f1u = new TF1("f1u",Get_Q4F1u_cont,0.1,15,0);
  // f1u->SetLineColor(kBlack);  
  // f1u->Draw();
  // TF1 *f1d = new TF1("f1d",Get_Q4F1d2p5_cont,0.1,15,0);
  // f1d->Draw("same");

  // //
  // TCanvas *c2 = util_pd::TC("c2",1,1);
  // c2->cd();

  // TF1 *f2u = new TF1("f2u",Get_Q4F2uovku_cont,0.1,15,0);
  // f2u->SetLineColor(kBlack);
  // f2u->Draw();
  // TF1 *f2d = new TF1("f2d",Get_Q4F2d0p75ovkd_cont,0.1,15,0);
  // f2d->Draw("same");    


  ////***************
  // TString outData = "flavordecomp.csv";
  // ofstream outdata; outdata.open(outData);
  // outdata << "Q2,Q4*F1u,Q4*F2u/ku,Q4*F1d*2.5,Q4*F2d*0.75/kd\n";
  // for (int i=0; i<7; i++) {    
  //   std::vector<double> FQ(4);
  //   GetFQuark(Q2,100,0,FQ);
  //   outdata << Q2 << "," << Q2*Q2*FQ[0] << "," << Q2*Q2*FQ[1]/ku << "," << Q2*Q2*FQ[2]*2.5 << "," << Q2*Q2*FQ[3]*0.75/kd << "\n";

  //   Q2 += 0.5;
  // }

  // TString outData = "flavordecomp_1.csv";
  // ofstream outdata; outdata.open(outData);

  // std::vector<double> Q2v{2.9890,4.4880,7.4640,9.8340,13.4650};
  
  // outdata << "Q2,Q4*F1u,Q4*F2u/ku,Q4*F1d*2.5,Q4*F2d*0.75/kd\n";
  // for (double item : Q2v) {
  //   Q2 = item;
  //   std::vector<double> FQ(4);
  //   GetFQuark(Q2,100,0,FQ);
  //   outdata << Q2 << "," << Q2*Q2*FQ[0] << "," << Q2*Q2*FQ[1]/ku << "," << Q2*Q2*FQ[2]*2.5 << "," << Q2*Q2*FQ[3]*0.75/kd << "\n";
  // }

  // TString outData = "flavordecomp_2.csv";
  // ofstream outdata; outdata.open(outData);

  // //std::vector<double> Q2v{3.0,4.5,7.4,9.9,13.5};
  // std::vector<double> Q2v{2.9890,4.4880,7.4640,9.8340,13.4650};
  // std::vector<double> SBSGMnovMuGD{0.9774,0.9763,0.9071,0.8473,0.7582};
  
  // outdata << "Q2,Q4*F1u,Q4*F2u/ku,Q4*F1d*2.5,Q4*F2d*0.75/kd\n";
  // for (size_t i = 0; i < Q2v.size(); ++i) {
  //   Q2 = Q2v[i];
  //   double GMn = SBSGMnovMuGD[i]*EMFFFits::GetGDip(Q2)*constant::mun;
  //   std::vector<double> FQ(4);
  //   GetFQuark(Q2,GMn,0,FQ);
  //   outdata << Q2 << "," << Q2*Q2*FQ[0] << "," << Q2*Q2*FQ[1]/ku << "," << Q2*Q2*FQ[2]*2.5 << "," << Q2*Q2*FQ[3]*0.75/kd << "\n";
  // }
  
  return 0;
}
