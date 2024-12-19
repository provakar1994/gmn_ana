#include <iostream>

#include "gmn_ana.h"


double const ku = 1.67;
double const kd = -2.03;

double ErrPropAplusB(double const AErr, double const BErr) {
  return sqrt(AErr*AErr + BErr*BErr);
}

double ErrPropAovB(double const A, double const AErr, double const B, double const BErr) {
  double R = A/B;
  return fabs(R) * sqrt( pow(AErr/A,2.) + pow(BErr/B,2.) );
}

std::vector<double> CalcF1(double Q2, double const GE, double const GEErr, double const GM, double const GMErr, std::string const & ntype) {

  double tau = ntype.compare("p")==0 ? kine::tau(Q2,"p") : kine::tau(Q2,"n");

  double numer = tau*GM + GE;
  double numerErr = ErrPropAplusB(tau*GMErr,GEErr);
  double denom = 1. + tau;
  double denomErr = 0.; //ErrPropAplusB(1.,tau); 
  
  return {numer/denom, ErrPropAovB(numer,numerErr,denom,denomErr)};
}

std::vector<double> CalcF2(double Q2, double const GE, double const GEErr, double const GM, double const GMErr, std::string const & ntype) {

  double tau = ntype.compare("p")==0 ? kine::tau(Q2,"p") : kine::tau(Q2,"n");

  double numer = GM - GE;
  double numerErr = ErrPropAplusB(GMErr,GEErr);
  double denom = 1. + tau;
  double denomErr = 0.; //ErrPropAplusB(1.,tau);
  
  return {numer/denom, ErrPropAovB(numer,numerErr,denom,denomErr)};
}

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



double Get_Q4F1u_cont(double *x, double *par) {
  /* Return a continuous curve for Q4*F1u using Ye 2017 fit */
  double Q2 = x[0];
  std::vector<double> FQ;

  GetFQuark(Q2,100,0,0,FQ);

  return Q2*Q2*FQ[0];
}

double Get_Q4F1d2p5_cont(double *x, double *par) {
  /* Return a continuous curve for Q4*F1d*2.5 using Ye 2017 fit */
  double Q2 = x[0];
  std::vector<double> FQ;

  GetFQuark(Q2,100,0,0,FQ);

  return Q2*Q2*FQ[4]*2.5;
}

double Get_Q4F2uovku_cont(double *x, double *par) {
  /* Return a continuous curve for Q4*F2u/ku using Ye 2017 fit */
  double Q2 = x[0];
  std::vector<double> FQ;

  GetFQuark(Q2,100,0,0,FQ);

  return Q2*Q2*FQ[2]/ku;
}

double Get_Q4F2d0p75ovkd_cont(double *x, double *par) {
  /* Return a continuous curve for Q4*F2d*0.75/kd using Ye 2017 fit */
  double Q2 = x[0];
  std::vector<double> FQ;

  GetFQuark(Q2,100,0,0,FQ);

  return Q2*Q2*FQ[6]*0.75/kd;
}

int flavordecomp(double Q2) {

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

  for( int i=0; i<=npoints_ye; i++ ){
    double Q2i = Q2min_ye + i*Q2step_ye;
    Q2ye[i] = Q2i; Q2yeErr[i] = 0.;
    
    std::vector<double> FQ;

    GetFQuark(Q2i,100,0,0,FQ);

    F1u[i] = Q2i*Q2i*FQ[0]; F1uErr[i] = Q2i*Q2i*FQ[1]; 
    F2u[i] = Q2i*Q2i*FQ[2]/ku; F2uErr[i] = Q2i*Q2i*FQ[3]/ku;
    F1d[i] = Q2i*Q2i*FQ[4]*2.5; F1dErr[i] = Q2i*Q2i*FQ[5]*2.5;     
    F2d[i] = Q2i*Q2i*FQ[6]*0.75/kd; F2dErr[i] = Q2i*Q2i*FQ[7]*0.75/kd; 
    
    // F1u[i] = FQ[0]; F1uErr[i] = FQ[1]; 
    // F2u[i] = FQ[2]; F2uErr[i] = FQ[3];
    // F1d[i] = FQ[4]; F1dErr[i] = FQ[5];     
    // F2d[i] = FQ[6]; F2dErr[i] = FQ[7]; 

  }  

  // std::cout << F1u[1] << " " << F1uErr[1] << "\n";
  
  TGraphErrors *F1uTG = new TGraphErrors( npoints_ye+1, Q2ye, F1u, Q2yeErr, F1uErr);
  TGraphErrors *F1dTG = new TGraphErrors( npoints_ye+1, Q2ye, F1d, Q2yeErr, F1dErr);
  TGraphErrors *F2uTG = new TGraphErrors( npoints_ye+1, Q2ye, F2u, Q2yeErr, F2uErr);
  TGraphErrors *F2dTG = new TGraphErrors( npoints_ye+1, Q2ye, F2d, Q2yeErr, F2dErr);

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
  TGraphErrors *F1dTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F1d_sbsgmn, Q2Err_sbsgmn, F1dErr_sbsgmn);
  TGraphErrors *F2uTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F2u_sbsgmn, Q2Err_sbsgmn, F2uErr_sbsgmn);
  TGraphErrors *F2dTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F2d_sbsgmn, Q2Err_sbsgmn, F2dErr_sbsgmn);  
  
  //**************
  TCanvas *c1 = util_pd::TC("c1",1,1);
  c1->cd();

  TH2D *hframe_f1 = new TH2D("hframe_f1",";Q^{2} (GeV/c)^{2};G_{M}^{n}/(#mu_{n}G_{D})",500,0,15,500,-0.6,1.6);

 // p->SetTopMargin(0.10);
  hframe_f1->SetStats(0);
  hframe_f1->GetXaxis()->CenterTitle();
  hframe_f1->GetYaxis()->CenterTitle();
  hframe_f1->GetXaxis()->SetLabelSize(0.04);
  hframe_f1->GetXaxis()->SetTitleSize(0.05);
  // hframe_f1->GetXaxis()->SetTitleOffset(0.9);
  hframe_f1->GetYaxis()->SetLabelSize(0.04);
  hframe_f1->GetYaxis()->SetTitleSize(0.05);  
  hframe_f1->GetYaxis()->SetNdivisions(7);
      
  hframe_f1->Draw();
  hframe_f1->GetXaxis()->CenterTitle();
  hframe_f1->GetYaxis()->CenterTitle();
  
  F1uTG->SetFillStyle( 1001 );
  F1uTG->SetFillColorAlpha( 1, 0.25 );
  F1uTG->SetLineWidth(2);
  F1uTG->SetLineColor(1);
  F1uTG->SetLineStyle(1);

  F1dTG->SetFillStyle( 1001 );
  F1dTG->SetFillColorAlpha( 2, 0.25 );
  F1dTG->SetLineWidth(2);
  F1dTG->SetLineColor(2);
  F1dTG->SetLineStyle(1);

  F1uTG_sbsgmn->SetMarkerColor(kBlack);
  F1uTG_sbsgmn->SetMarkerStyle(20);

  F1dTG_sbsgmn->SetMarkerColor(kRed);
  F1dTG_sbsgmn->SetMarkerStyle(20);
  F1dTG_sbsgmn->SetLineColor(kRed);
  
  F1uTG->Draw("C3 SAME");
  F1dTG->Draw("C3 SAME");

  F1uTG_sbsgmn->Draw("P SAME");
  F1dTG_sbsgmn->Draw("P SAME");  

  //**************
  TCanvas *c2 = util_pd::TC("c2",1,1);
  c2->cd();

  TH2D *hframe_f2 = new TH2D("hframe_f2",";Q^{2} (GeV/c)^{2};G_{M}^{n}/(#mu_{n}G_{D})",500,0,15,500,-0.02,0.3);

 // p->SetTopMargin(0.10);
  hframe_f2->SetStats(0);
  hframe_f2->GetXaxis()->CenterTitle();
  hframe_f2->GetYaxis()->CenterTitle();
  hframe_f2->GetXaxis()->SetLabelSize(0.04);
  hframe_f2->GetXaxis()->SetTitleSize(0.05);
  // hframe_f2->GetXaxis()->SetTitleOffset(0.9);
  hframe_f2->GetYaxis()->SetLabelSize(0.04);
  hframe_f2->GetYaxis()->SetTitleSize(0.05);  
  hframe_f2->GetYaxis()->SetNdivisions(7);
      
  hframe_f2->Draw();
  hframe_f2->GetXaxis()->CenterTitle();
  hframe_f2->GetYaxis()->CenterTitle();
  
  F2uTG->SetFillStyle( 1001 );
  F2uTG->SetFillColorAlpha( 1, 0.25 );
  F2uTG->SetLineWidth(2);
  F2uTG->SetLineColor(1);
  F2uTG->SetLineStyle(1);

  F2dTG->SetFillStyle( 1001 );
  F2dTG->SetFillColorAlpha( 2, 0.25 );
  F2dTG->SetLineWidth(2);
  F2dTG->SetLineColor(2);
  F2dTG->SetLineStyle(1);

  F2uTG_sbsgmn->SetMarkerColor(kBlack);
  F2uTG_sbsgmn->SetMarkerStyle(20);

  F2dTG_sbsgmn->SetMarkerColor(kRed);
  F2dTG_sbsgmn->SetMarkerStyle(20);
  F2dTG_sbsgmn->SetLineColor(kRed);  
  
  F2uTG->Draw("C3 SAME");
  F2dTG->Draw("C3 SAME");    

  F2uTG_sbsgmn->Draw("P SAME");
  F2dTG_sbsgmn->Draw("P SAME");  
  
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
