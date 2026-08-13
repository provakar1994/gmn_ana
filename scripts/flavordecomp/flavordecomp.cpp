#include <iostream>
#include "TSpline.h"

#include "gmn_ana.h"


double const ku = 1.67;
double const kd = -2.03;

//_______________________________________
double ErrPropAplusB(double const AErr, double const BErr) {
  return sqrt(AErr*AErr + BErr*BErr);
}

//_______________________________________
double ErrPropAplusBplusC(double const AErr, double const BErr, double const CErr) {
  return sqrt(AErr*AErr + BErr*BErr + CErr*CErr);
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
  double numerErr = ErrPropAplusB(tau*abs(GMErr),abs(GEErr));
  double denom = 1. + tau;
  double denomErr = 0.; //ErrPropAplusB(1.,tau); 
  
  return {numer/denom, ErrPropAovB(numer,numerErr,denom,denomErr)};
}

//_______________________________________
std::vector<double> CalcF2(double Q2, double const GE, double const GEErr, double const GM, double const GMErr, std::string const & ntype) {

  double tau = ntype.compare("p")==0 ? kine::tau(Q2,"p") : kine::tau(Q2,"n");

  double numer = GM - GE;
  double numerErr = ErrPropAplusB(abs(GMErr),abs(GEErr));
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
  //std::cout << Form("F2nErr,F2pErr,F2dErr = %f,%f,%f \n",F2nErr,F2pErr,F2dErr);
  
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
  //std::cout << F2p << "," << F2pErr << "," << F2n << "," << F2nErr << "\n";
  double F2d = temp[0], F2dErr = temp[1];

  if (verbose==1) {
    double GD = EMFFFits::GetGDip(Q2);
    std::cout << Form("\nQ2 = %f GeV2 \n",Q2);
    std::cout << Form("GEp: %f #pm %f \n",0.2*GEp_ye/GD,0.2*GEp_yeErr/GD);
    std::cout << Form("GMp: %f #pm %f \n",0.8*GMp_ye/GD,0.8*GMp_yeErr/GD);
    std::cout << Form("GEn: %f #pm %f \n",0.4*GEn_ye/GD,0.4*GEn_yeErr/GD);
    std::cout << Form("GMn: %f #pm %f\n",1.6*GMn_ye/GD,1.6*GMn_yeErr/GD);          
    //
    // std::cout << Form("GEp: %f #pm %f \n",GEp_ye,GEp_yeErr);
    // std::cout << Form("GMp: %f #pm %f \n",GMp_ye,GMp_yeErr);
    // std::cout << Form("GEn: %f #pm %f \n",GEn_ye,GEn_yeErr);
    // std::cout << Form("GMn: %f #pm %f\n",GMn_ye,GMn_yeErr);          
    //std::cout << Form("GMn: %f #pm %f\n",GMn_ye/(EMFFFits::GetGDip(Q2)*constant::mun),fabs(GMn_yeErr/(EMFFFits::GetGDip(Q2)*constant::mun)));      
  }
  if (verbose>1) {
    std::cout << Form("GEp: %f, GMp: %f, GEn: %f, GMn: %f\n",GEp_ye,GMp_ye,GEn_ye,GMn_ye);
    std::cout << Form("F1p: %f, F2p: %f, F1n: %f, F2n: %f\n",F1p,F2p,F1n,F2n);
    std::cout << Form("F1u: %f, F2u: %f, F1d: %f, F2d: %f\n",F1u,F2u,F1d,F2d);
    std::cout << Form("Q4*F1u: %f, Q4*F2u/ku: %f, Q4*F1d*2.5: %f, Q4*F2d*0.75/kd: %f\n",Q2*Q2*F1u,Q2*Q2*F2u/ku,Q2*Q2*F1d*2.5,Q2*Q2*F2d*0.75/kd);
  }

  FQ = {F1u,F1uErr,F2u,F2uErr,F1d,F1dErr,F2d,F2dErr};
}

//_______________________________________
std::vector<double> CalcGMQuark(double const GMp, double const GMpErr,  double const GMn, double const GMnErr, std::string const & qtype) {

  double GMu = 2.*GMp + GMn;
  double GMuErr = ErrPropAplusB(2.*GMpErr,GMnErr);
  double GMd = 2.*GMn + GMp;
  double GMdErr = ErrPropAplusB(2.*GMnErr,GMpErr);

  std::vector<double> result;
  if (qtype.compare("u")==0) result = {GMu,GMuErr};
  else result = {GMd,GMdErr}; 
    
  return result; 
}

//_______________________________________
std::vector<double> CalcGEQuark(double const GEp, double const GEpErr, double const GEn, double const GEnErr, std::string const & qtype) {

  double GEu = 2.*GEp + GEn;
  double GEuErr = ErrPropAplusB(2.*GEpErr,GEnErr);  
  double GEd = 2.*GEn + GEp;
  double GEdErr = ErrPropAplusB(2.*GEnErr,GEpErr);
  
  std::vector<double> result;
  if (qtype.compare("u")==0) result = {GEu,GEuErr};
  else result = {GEd,GEdErr}; 
    
  return result;   
}

//_______________________________________
void GetGQuark(double Q2, double GMn, double GMnErr, int verbose, std::vector<double> &GQ) {

  // EMFF fits
  Ye2017 yefit;
  // Kelly2004 kellyfit;
  // Seamus20XX seamusfit;
  // Christy2022 christyfit;
  // Arrington2007 arfit;
  // Galster1971 galfit;

  // True
  std::vector<double> temp = yefit.GetFFwErr(G_t::kGEp,Q2);
  double GEp = temp[0], GEpErr = temp[1];
  temp = yefit.GetFFwErr(G_t::kGMp,Q2);
  double GMp = temp[0], GMpErr = temp[1];
  temp = yefit.GetFFwErr(G_t::kGEn,Q2);
  double GEn = temp[0], GEnErr = temp[1];
  temp = yefit.GetFFwErr(G_t::kGMn,Q2);
  double GMn_new = GMn>0 ? temp[0] : GMn;
  double GMnErr_new = GMn>0 ? temp[1] : GMnErr;  

  temp = CalcGMQuark(GMp,GMpErr,GMn_new,GMnErr_new,"u");
  double GMu = temp[0], GMuErr = temp[1];
  temp = CalcGMQuark(GMp,GMpErr,GMn_new,GMnErr_new,"d");  
  double GMd = temp[0], GMdErr = temp[1];
  temp = CalcGEQuark(GEp,GEpErr,GEn,GEnErr,"u");
  double GEu = temp[0], GEuErr = temp[1];
  temp = CalcGEQuark(GEp,GEpErr,GEn,GEnErr,"d");
  double GEd = temp[0], GEdErr = temp[1];

  if (verbose>0) {
    //std::cout << Form("GEp: %f, GMp: %f, GEn: %f, GMn: %f\n",GEp,GMp,GEn,GMn);
    std::cout << Form("GMp: %f, GEp: %f, GMn: %f, GEn: %f\n",GMp,GEp,GMn,GEn);
    std::cout << Form("GMu: %f, GEu: %f, GMd: %f, GEd: %f\n",GMu,GEu,GMd,GEd);
    // std::cout << Form("Q4*F1u: %f, Q4*F2u/ku: %f, Q4*F1d*2.5: %f, Q4*F2d*0.75/kd: %f\n",Q2*Q2*F1u,Q2*Q2*F2u/ku,Q2*Q2*F1d*2.5,Q2*Q2*F2d*0.75/kd);
  }

  GQ = {GMu,GMuErr,GEu,GEuErr,GMd,GMdErr,GEd,GEdErr};
}

//_______________________________________
void GetGQFromFQ(double Q2, double F1u, double F1uErr, double F2u, double F2uErr, double F1d, double F1dErr, double F2d, double F2dErr, int verbose, std::vector<double> &GQ) {
  double tau_p = kine::tau(Q2,"p");
  double tau_n = kine::tau(Q2,"n");

  double GEu = F1u - (1./3.)*(4.*tau_p-tau_n)*F2u + (2./3.)*(tau_p-tau_n)*F2d;
  double GEuErr = ErrPropAplusBplusC(F1uErr,(1./3.)*(4.*tau_p-tau_n)*F2uErr,(2./3.)*(tau_p-tau_n)*F2dErr);

  double GMu = F1u + F2u;
  double GMuErr = ErrPropAplusB(F1uErr,F2uErr);

  double GEd = F1d - (1./3.)*(4.*tau_n-tau_p)*F2d + (2./3.)*(tau_n-tau_p)*F2u;
  double GEdErr = ErrPropAplusBplusC(F1dErr,(1./3.)*(4.*tau_n-tau_p)*F2dErr,(2./3.)*(tau_n-tau_p)*F2uErr);
  
  double GMd = F1d + F2d;
  double GMdErr = ErrPropAplusB(F1dErr,F2dErr);

  GQ = {GMu,GMuErr,GEu,GEuErr,GMd,GMdErr,GEd,GEdErr};
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
void customize_cates(TGraphErrors *g, std::string qtype) {
  g->SetMarkerStyle(24);
  if (qtype.compare("u")==0) {
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
  } else {
    g->SetMarkerColor(kRed);
    g->SetLineColor(kRed);
  }
}

//------------------------------------------------------------------------------
// MakeSmoothedGraphErrors:
//   • builds a cubic spline through (x,y) 
//   • builds two splines through (x,y+ey) and (x,y–ey)
//   • returns a new TGraphErrors with nSteps points whose
//     y = spline(x)  and  ey = ½[spline_up(x)–spline_dn(x)]
//------------------------------------------------------------------------------
TGraphErrors* MakeSmoothedGraphErrors(const TGraphErrors* graw,
                                      Int_t nSteps = 200)
{
  // 1) read raw points
  const Int_t n0 = graw->GetN();
  if (n0 < 2 || nSteps < 2) {
    ::Error("MakeSmoothedGraphErrors","need at least 2 points and nSteps>=2");
    return nullptr;
  }
  std::vector<Double_t> x0(n0), y0(n0), yup(n0), ydn(n0);
  for (Int_t i = 0; i < n0; ++i) {
    Double_t ex, ey;
    graw->GetPoint(i, x0[i], y0[i]);
    ey = graw->GetErrorY(i);
    yup[i] = y0[i] + ey;
    ydn[i] = y0[i] - ey;
  }

  // 2) build splines
  TSpline3* splC = new TSpline3("splC", &x0[0], &y0[0],  n0);
  TSpline3* splU = new TSpline3("splU", &x0[0], &yup[0], n0);
  TSpline3* splD = new TSpline3("splD", &x0[0], &ydn[0], n0);

  // 3) prepare arrays for the new graph
  std::vector<Double_t> xs(nSteps), ys(nSteps), exs(nSteps, 0.0), eys(nSteps);
  Double_t xmin = x0.front(), xmax = x0.back();
  Double_t dx   = (xmax - xmin) / (nSteps - 1);

  for (Int_t i = 0; i < nSteps; ++i) {
    Double_t xx = xmin + i * dx;
    xs[i]  = xx;
    Double_t yc = splC->Eval(xx);
    Double_t yu = splU->Eval(xx);
    Double_t yd = splD->Eval(xx);
    ys[i]  = yc;
    eys[i] = 0.5*(yu - yd);              // half‐width of the band
  }

  // 4) build and return the smoothed TGraphErrors
  TGraphErrors* gsm = new TGraphErrors(nSteps,
                                       &xs[0], &ys[0],
                                       &exs[0], &eys[0]);
  return gsm;
}

//_______________________________________
int flavordecomp() {

  // Initializing parameters and arrays to store Ye 2018 fit values
  double Q2min_ye = 1.e-5;
  double Q2max_ye = 15.0;
  int npoints_ye = 2000;
  double Q2step_ye = (Q2max_ye-Q2min_ye)/double(npoints_ye);
  //
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
  //
  double GMu[npoints_ye+1];
  double GMuOVGD[npoints_ye+1];
  double GMuOVmupGD[npoints_ye+1];  
  double GMuErr[npoints_ye+1];
  double GMuOVGDErr[npoints_ye+1];
  double GMuOVmupGDErr[npoints_ye+1];  
  double GEu[npoints_ye+1];
  double GEuOVGMu[npoints_ye+1];
  double mupGEuOVGMu[npoints_ye+1];  
  double GEuErr[npoints_ye+1];
  double GEuOVGMuErr[npoints_ye+1];
  double mupGEuOVGMuErr[npoints_ye+1];    
  double GMd[npoints_ye+1];
  double GMdOVGD[npoints_ye+1];
  double GMdOVmunGD[npoints_ye+1];  
  double GMdErr[npoints_ye+1];
  double GMdOVGDErr[npoints_ye+1];
  double GMdOVmunGDErr[npoints_ye+1];  
  double GEd[npoints_ye+1];
  double GEdOVGMd[npoints_ye+1];
  double munGEdOVGMd[npoints_ye+1];  
  double GEdErr[npoints_ye+1];
  double GEdOVGMdErr[npoints_ye+1];
  double munGEdOVGMdErr[npoints_ye+1];  
  
  // Reading Hague's GMn table
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
  //
  double GMu_tyler[npoints_ye+1];
  double GMuOVGD_tyler[npoints_ye+1];
  double GMuOVmupGD_tyler[npoints_ye+1];  
  double GMuErr_tyler[npoints_ye+1];  
  double GMuOVGDErr_tyler[npoints_ye+1];
  double GMuOVmupGDErr_tyler[npoints_ye+1];  
  double GEu_tyler[npoints_ye+1];
  double GEuOVGMu_tyler[npoints_ye+1];
  double mupGEuOVGMu_tyler[npoints_ye+1];   
  double GEuErr_tyler[npoints_ye+1];  
  double GEuOVGMuErr_tyler[npoints_ye+1];
  double mupGEuOVGMuErr_tyler[npoints_ye+1];   
  double GMd_tyler[npoints_ye+1];
  double GMdOVGD_tyler[npoints_ye+1];
  double GMdOVmunGD_tyler[npoints_ye+1];  
  double GMdErr_tyler[npoints_ye+1];  
  double GMdOVGDErr_tyler[npoints_ye+1];
  double GMdOVmunGDErr_tyler[npoints_ye+1];  
  double GEd_tyler[npoints_ye+1];
  double GEdOVGMd_tyler[npoints_ye+1];
  double munGEdOVGMd_tyler[npoints_ye+1];  
  double GEdErr_tyler[npoints_ye+1];
  double GEdOVGMdErr_tyler[npoints_ye+1];
  double munGEdOVGMdErr_tyler[npoints_ye+1];

  // sanity checks
  // 1. reproducing GMn from GMu and GMd
  double GMnSC[npoints_ye+1];
  double GMnSCErr[npoints_ye+1];
  double GMnOVmunGDSC[npoints_ye+1];
  double GMnOVmunGDSCErr[npoints_ye+1];
  

  for( int i=0; i<=npoints_ye; i++ ){
    double Q2i = Q2min_ye + i*Q2step_ye;
    Q2ye[i] = Q2i; Q2yeErr[i] = 0.;
    
    std::vector<double> FQ;
    GetFQuark(Q2i,100,0,0,FQ);

    F1u[i] = Q2i*Q2i*FQ[0]; F1uErr[i] = Q2i*Q2i*FQ[1]; 
    F2u[i] = Q2i*Q2i*FQ[2]/ku; F2uErr[i] = Q2i*Q2i*FQ[3]/ku;
    F1d[i] = Q2i*Q2i*FQ[4]*2.5; F1dErr[i] = Q2i*Q2i*FQ[5]*2.5;     
    F2d[i] = Q2i*Q2i*FQ[6]*0.75/kd; F2dErr[i] = Q2i*Q2i*FQ[7]*0.75/abs(kd);

    std::vector<double> GQ;
    GetGQuark(Q2i,100,0,0,GQ);
    
    GMu[i] = GQ[0]; GMuErr[i] = GQ[1]; 
    GMuOVGD[i] = GQ[0]/EMFFFits::GetGDip(Q2i); GMuOVGDErr[i] = GQ[1]/EMFFFits::GetGDip(Q2i);
    GMuOVmupGD[i] = GQ[0]/(EMFFFits::GetGDip(Q2i)*constant::mup); GMuOVmupGDErr[i] = GQ[1]/(EMFFFits::GetGDip(Q2i)*constant::mup);    
    GEu[i] = GQ[2]; GEuErr[i] = GQ[3];
    GEuOVGMu[i] = GQ[2]/GQ[0]; GEuOVGMuErr[i] = ErrPropAovB(GQ[2],GQ[3],GQ[0],GQ[1]);    
    mupGEuOVGMu[i] = constant::mup*GQ[2]/GQ[0]; mupGEuOVGMuErr[i] = ErrPropAovB(constant::mup*GQ[2],constant::mup*GQ[3],GQ[0],GQ[1]);    
    GMd[i] = GQ[4]; GMdErr[i] = GQ[5];     
    GMdOVGD[i] = GQ[4]/EMFFFits::GetGDip(Q2i); GMdOVGDErr[i] = GQ[5]/EMFFFits::GetGDip(Q2i);
    GMdOVmunGD[i] = GQ[4]/(EMFFFits::GetGDip(Q2i)*constant::mun); GMdOVmunGDErr[i] = GQ[5]/(EMFFFits::GetGDip(Q2i)*constant::mun);    
    GEd[i] = GQ[6]; GEdErr[i] = GQ[7];
    GEdOVGMd[i] = GQ[6]/GQ[4]; GEdOVGMdErr[i] = ErrPropAovB(GQ[6],GQ[7],GQ[4],GQ[5]);    
    munGEdOVGMd[i] = constant::mun*GQ[6]/GQ[4]; munGEdOVGMdErr[i] = ErrPropAovB(constant::mun*GQ[6],constant::mun*GQ[7],GQ[4],GQ[5]);    

    double GMn_tyler = reader.GetClosestValueByKey(Q2i,0)*EMFFFits::GetGDip(Q2i)*constant::mun;
    double GMn_tyler_err = reader.GetClosestValueByKey(Q2i,1)*EMFFFits::GetGDip(Q2i)*constant::mun;     
    GetFQuark(Q2i,GMn_tyler,GMn_tyler_err,0,FQ);
    F1u_tyler[i] = Q2i*Q2i*FQ[0]; F1uErr_tyler[i] = Q2i*Q2i*FQ[1]; 
    F2u_tyler[i] = Q2i*Q2i*FQ[2]/ku; F2uErr_tyler[i] = Q2i*Q2i*FQ[3]/ku;
    F1d_tyler[i] = Q2i*Q2i*FQ[4]*2.5; F1dErr_tyler[i] = Q2i*Q2i*FQ[5]*2.5;     
    F2d_tyler[i] = Q2i*Q2i*FQ[6]*0.75/kd; F2dErr_tyler[i] = Q2i*Q2i*FQ[7]*0.75/abs(kd);     

    GetGQuark(Q2i,GMn_tyler,GMn_tyler_err,0,GQ);        
    GMu_tyler[i] = GQ[0]; GMuErr_tyler[i] = GQ[1]; 
    GMuOVGD_tyler[i] = GQ[0]/EMFFFits::GetGDip(Q2i); GMuOVGDErr_tyler[i] = GQ[1]/EMFFFits::GetGDip(Q2i);
    GMuOVmupGD_tyler[i] = GQ[0]/(EMFFFits::GetGDip(Q2i)*constant::mup); GMuOVmupGDErr_tyler[i] = GQ[1]/(EMFFFits::GetGDip(Q2i)*constant::mup);    
    GEu_tyler[i] = GQ[2]; GEuErr_tyler[i] = GQ[3];
    GEuOVGMu_tyler[i] = GQ[2]/GQ[0]; GEuOVGMuErr_tyler[i] = ErrPropAovB(GQ[2],GQ[3],GQ[0],GQ[1]);    
    mupGEuOVGMu_tyler[i] = constant::mup*GQ[2]/GQ[0]; mupGEuOVGMuErr_tyler[i] = ErrPropAovB(constant::mup*GQ[2],constant::mup*GQ[3],GQ[0],GQ[1]);    
    GMd_tyler[i] = GQ[4]; GMdErr_tyler[i] = GQ[5];     
    GMdOVGD_tyler[i] = GQ[4]/EMFFFits::GetGDip(Q2i); GMdOVGDErr_tyler[i] = GQ[5]/EMFFFits::GetGDip(Q2i);
    GMdOVmunGD_tyler[i] = GQ[4]/(EMFFFits::GetGDip(Q2i)*constant::mun); GMdOVmunGDErr_tyler[i] = GQ[5]/(EMFFFits::GetGDip(Q2i)*constant::mun);    
    GEd_tyler[i] = GQ[6]; GEdErr_tyler[i] = GQ[7];     
    GEdOVGMd_tyler[i] = GQ[6]/GQ[4]; GEdOVGMdErr_tyler[i] = ErrPropAovB(GQ[6],GQ[7],GQ[4],GQ[5]);    
    munGEdOVGMd_tyler[i] = constant::mun*GQ[6]/GQ[4]; munGEdOVGMdErr_tyler[i] = ErrPropAovB(constant::mun*GQ[6],constant::mun*GQ[7],GQ[4],GQ[5]);    

    // sanity checks
    // 1. reproducing GMn from GMu and GMd
    GMnSC[i] = (2./3.)*GMd[i] - (1./3.)*GMu[i];
    GMnSCErr[i] = ErrPropAplusB((2./3.)*GMdErr[i], (1./3.)*GMuErr[i]); 
    GMnOVmunGDSC[i] = GMnSC[i]/(EMFFFits::GetGDip(Q2i)*constant::mun);
    GMnOVmunGDSCErr[i] = GMnSCErr[i]/(EMFFFits::GetGDip(Q2i)*constant::mun); 

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
  //
  TGraphErrors *GMuTG = new TGraphErrors( npoints_ye+1, Q2ye, GMu, Q2yeErr, GMuErr);
  customize_gfit(GMuTG,"u");
  TGraphErrors *GMuOVGDTG = new TGraphErrors( npoints_ye+1, Q2ye, GMuOVGD, Q2yeErr, GMuOVGDErr);
  customize_gfit(GMuOVGDTG,"u");
  TGraphErrors *GMuOVmupGDTG = new TGraphErrors( npoints_ye+1, Q2ye, GMuOVmupGD, Q2yeErr, GMuOVmupGDErr);
  customize_gfit(GMuOVmupGDTG,"u");
  TGraphErrors *GMdTG = new TGraphErrors( npoints_ye+1, Q2ye, GMd, Q2yeErr, GMdErr);
  customize_gfit(GMdTG,"d");
  TGraphErrors *GMdOVGDTG = new TGraphErrors( npoints_ye+1, Q2ye, GMdOVGD, Q2yeErr, GMdOVGDErr);
  customize_gfit(GMdOVGDTG,"d");
  TGraphErrors *GMdOVmunGDTG = new TGraphErrors( npoints_ye+1, Q2ye, GMdOVmunGD, Q2yeErr, GMdOVmunGDErr);
  customize_gfit(GMdOVmunGDTG,"d");  
  TGraphErrors *GEuTG = new TGraphErrors( npoints_ye+1, Q2ye, GEu, Q2yeErr, GEuErr);
  customize_gfit(GEuTG,"u");
  TGraphErrors *GEuOVGMuTG = new TGraphErrors( npoints_ye+1, Q2ye, GEuOVGMu, Q2yeErr, GEuOVGMuErr);
  customize_gfit(GEuOVGMuTG,"u");
  TGraphErrors *mupGEuOVGMuTG = new TGraphErrors( npoints_ye+1, Q2ye, mupGEuOVGMu, Q2yeErr, mupGEuOVGMuErr);
  customize_gfit(mupGEuOVGMuTG,"u");
  TGraphErrors *GEdTG = new TGraphErrors( npoints_ye+1, Q2ye, GEd, Q2yeErr, GEdErr);
  customize_gfit(GEdTG,"d");
  TGraphErrors *GEdOVGMdTG = new TGraphErrors( npoints_ye+1, Q2ye, GEdOVGMd, Q2yeErr, GEdOVGMdErr);
  customize_gfit(GEdOVGMdTG,"d");
  TGraphErrors *munGEdOVGMdTG = new TGraphErrors( npoints_ye+1, Q2ye, munGEdOVGMd, Q2yeErr, munGEdOVGMdErr);
  customize_gfit(munGEdOVGMdTG,"d");
  
  // Using Hague GMn + Ye 2018
  TGraphErrors *F1uTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, F1u_tyler, Q2yeErr, F1uErr_tyler);
  //TGraphErrors *F1uTG_tyler = MakeSmoothedGraphErrors(F1uTG_tyler_coarse, 5000);
  customize_gfit(F1uTG_tyler,"u");
  TGraphErrors *F1dTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, F1d_tyler, Q2yeErr, F1dErr_tyler);
  //TGraphErrors *F1dTG_tyler = MakeSmoothedGraphErrors(F1dTG_tyler_coarse, 5000);
  customize_gfit(F1dTG_tyler,"d");
  TGraphErrors *F2uTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, F2u_tyler, Q2yeErr, F2uErr_tyler);
  customize_gfit(F2uTG_tyler,"u");
  TGraphErrors *F2dTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, F2d_tyler, Q2yeErr, F2dErr_tyler);
  customize_gfit(F2dTG_tyler,"d");  
  //
  TGraphErrors *GMuTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GMu_tyler, Q2yeErr, GMuErr_tyler);
  customize_gfit(GMuTG_tyler,"u");
  TGraphErrors *GMuOVGDTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GMuOVGD_tyler, Q2yeErr, GMuOVGDErr_tyler);
  customize_gfit(GMuOVGDTG_tyler,"u");
  TGraphErrors *GMuOVmupGDTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GMuOVmupGD_tyler, Q2yeErr, GMuOVmupGDErr_tyler);
  customize_gfit(GMuOVmupGDTG_tyler,"u");  
  TGraphErrors *GMdTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GMd_tyler, Q2yeErr, GMdErr_tyler);
  customize_gfit(GMdTG_tyler,"d");
  TGraphErrors *GMdOVGDTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GMdOVGD_tyler, Q2yeErr, GMdOVGDErr_tyler);
  customize_gfit(GMdOVGDTG_tyler,"d");
  TGraphErrors *GMdOVmunGDTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GMdOVmunGD_tyler, Q2yeErr, GMdOVmunGDErr_tyler);
  customize_gfit(GMdOVmunGDTG_tyler,"d");  
  TGraphErrors *GEuTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GEu_tyler, Q2yeErr, GEuErr_tyler);
  customize_gfit(GEuTG_tyler,"u");
  TGraphErrors *GEuOVGMuTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GEuOVGMu_tyler, Q2yeErr, GEuOVGMuErr_tyler);
  customize_gfit(GEuOVGMuTG_tyler,"u");
  TGraphErrors *mupGEuOVGMuTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, mupGEuOVGMu_tyler, Q2yeErr, mupGEuOVGMuErr_tyler);
  customize_gfit(mupGEuOVGMuTG_tyler,"u");
  TGraphErrors *GEdTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GEd_tyler, Q2yeErr, GEdErr_tyler);
  customize_gfit(GEdTG_tyler,"d");
  TGraphErrors *GEdOVGMdTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, GEdOVGMd_tyler, Q2yeErr, GEdOVGMdErr_tyler);
  customize_gfit(GEdOVGMdTG_tyler,"d");
  TGraphErrors *munGEdOVGMdTG_tyler = new TGraphErrors( npoints_ye+1, Q2ye, munGEdOVGMd_tyler, Q2yeErr, munGEdOVGMdErr_tyler);
  customize_gfit(munGEdOVGMdTG_tyler,"d");  

  // sanity checks
  // 1.
  TGraphErrors *GMnSCTG = new TGraphErrors( npoints_ye+1, Q2ye, GMnSC, Q2yeErr, GMnSCErr);
  customize_gfit(GMnSCTG,"u");
  TGraphErrors *GMnOVmunGDSCTG = new TGraphErrors( npoints_ye+1, Q2ye, GMnOVmunGDSC, Q2yeErr, GMnOVmunGDSCErr);
  customize_gfit(GMnOVmunGDSCTG,"u");
  
  // Using GMn data
  //std::vector<double> Q2v{3.0,4.5,7.4,9.9,13.5};
  std::vector<double> Q2v{2.9890,4.4880,7.4640,9.8340,13.4650}; 
  // std::vector<double> SBSGMnovMuGD{0.9774,0.9763,0.9071,0.8473,0.7582}; // thesis
  // std::vector<double> SBSGMnErrovMuGD{0.0145,0.0164,0.0174,0.0245,0.0226}; // thesis
  // std::vector<double> SBSGMnovMuGD{0.9696,0.9533,0.8870,0.8185,0.7314}; // MCp2
  // std::vector<double> SBSGMnErrovMuGD{0.0180,0.0217,0.0240,0.0278,0.0277}; // MCp2
  std::vector<double> SBSGMnovMuGD{0.9674,0.9563,0.8829,0.8257,0.7347}; // pass3
  std::vector<double> SBSGMnErrovMuGD{0.0173,0.0210,0.0227,0.0286,0.0282}; // pass3

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
  //
  double GMu_sbsgmn[npoints_sbsgmn];
  double GMuOVGD_sbsgmn[npoints_sbsgmn];
  double GMuOVmupGD_sbsgmn[npoints_sbsgmn];  
  double GMuErr_sbsgmn[npoints_sbsgmn];
  double GMuOVGDErr_sbsgmn[npoints_sbsgmn];
  double GMuOVmupGDErr_sbsgmn[npoints_sbsgmn];  
  double GEu_sbsgmn[npoints_sbsgmn];
  double GEuOVGMu_sbsgmn[npoints_sbsgmn];
  double mupGEuOVGMu_sbsgmn[npoints_sbsgmn];  
  double GEuErr_sbsgmn[npoints_sbsgmn];
  double GEuOVGMuErr_sbsgmn[npoints_sbsgmn];
  double mupGEuOVGMuErr_sbsgmn[npoints_sbsgmn];  
  double GMd_sbsgmn[npoints_sbsgmn];
  double GMdOVGD_sbsgmn[npoints_sbsgmn];
  double GMdOVmunGD_sbsgmn[npoints_sbsgmn];  
  double GMdErr_sbsgmn[npoints_sbsgmn];
  double GMdOVGDErr_sbsgmn[npoints_sbsgmn];
  double GMdOVmunGDErr_sbsgmn[npoints_sbsgmn];    
  double GEd_sbsgmn[npoints_sbsgmn];
  double GEdOVGMd_sbsgmn[npoints_sbsgmn];
  double munGEdOVGMd_sbsgmn[npoints_sbsgmn];   
  double GEdErr_sbsgmn[npoints_sbsgmn];
  double GEdOVGMdErr_sbsgmn[npoints_sbsgmn];
  double munGEdOVGMdErr_sbsgmn[npoints_sbsgmn];    
  
  for (size_t i=0; i<npoints_sbsgmn; i++) {
    Q2_sbsgmn[i] = Q2v[i];
    Q2Err_sbsgmn[i] = 0.;
    
    std::vector<double> FQ;
    GetFQuark(Q2v[i],SBSGMnovMuGD[i]*EMFFFits::GetGDip(Q2v[i])*constant::mun,fabs(SBSGMnErrovMuGD[i]*EMFFFits::GetGDip(Q2v[i])*constant::mun),0,FQ);

    F1u_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[0]; F1uErr_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[1]; 
    F2u_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[2]/ku; F2uErr_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[3]/ku;
    F1d_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[4]*2.5; F1dErr_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[5]*2.5;     
    F2d_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[6]*0.75/kd; F2dErr_sbsgmn[i] = Q2v[i]*Q2v[i]*FQ[7]*0.75/abs(kd);

    std::vector<double> GQ;
    GetGQuark(Q2v[i],SBSGMnovMuGD[i]*EMFFFits::GetGDip(Q2v[i])*constant::mun,fabs(SBSGMnErrovMuGD[i]*EMFFFits::GetGDip(Q2v[i])*constant::mun),0,GQ);
    
    GMu_sbsgmn[i] = GQ[0]; GMuErr_sbsgmn[i] = GQ[1];
    GMuOVGD_sbsgmn[i] = GQ[0]/EMFFFits::GetGDip(Q2v[i]); GMuOVGDErr_sbsgmn[i] = GQ[1]/EMFFFits::GetGDip(Q2v[i]);
    GMuOVmupGD_sbsgmn[i] = GQ[0]/(EMFFFits::GetGDip(Q2v[i])*constant::mup); GMuOVmupGDErr_sbsgmn[i] = GQ[1]/(EMFFFits::GetGDip(Q2v[i])*constant::mup);        
    GEu_sbsgmn[i] = GQ[2]; GEuErr_sbsgmn[i] = GQ[3];
    GEuOVGMu_sbsgmn[i] = GQ[2]/GQ[0]; GEuOVGMuErr_sbsgmn[i] = ErrPropAovB(GQ[2],GQ[3],GQ[0],GQ[1]);    
    mupGEuOVGMu_sbsgmn[i] = constant::mup*GQ[2]/GQ[0]; mupGEuOVGMuErr_sbsgmn[i] = ErrPropAovB(constant::mup*GQ[2],constant::mup*GQ[3],GQ[0],GQ[1]);        
    GMd_sbsgmn[i] = GQ[4]; GMdErr_sbsgmn[i] = GQ[5];
    GMdOVGD_sbsgmn[i] = GQ[4]/EMFFFits::GetGDip(Q2v[i]); GMdOVGDErr_sbsgmn[i] = GQ[5]/EMFFFits::GetGDip(Q2v[i]);
    GMdOVmunGD_sbsgmn[i] = GQ[4]/(EMFFFits::GetGDip(Q2v[i])*constant::mun); GMdOVmunGDErr_sbsgmn[i] = GQ[5]/(EMFFFits::GetGDip(Q2v[i])*constant::mun);
    GEd_sbsgmn[i] = GQ[6]; GEdErr_sbsgmn[i] = GQ[7];
    GEdOVGMd_sbsgmn[i] = GQ[6]/GQ[4]; GEdOVGMdErr_sbsgmn[i] = ErrPropAovB(GQ[6],GQ[7],GQ[4],GQ[5]);    
    munGEdOVGMd_sbsgmn[i] = constant::mun*GQ[6]/GQ[4]; munGEdOVGMdErr_sbsgmn[i] = ErrPropAovB(constant::mun*GQ[6],constant::mun*GQ[7],GQ[4],GQ[5]);

    // std::cout << Q2_sbsgmn[i] << ",";
    // std::cout << F1u_sbsgmn[i] << "," << F1uErr_sbsgmn[i] << ",";
    // std::cout << F1d_sbsgmn[i] << "," << F1dErr_sbsgmn[i] << ",";
    // std::cout << F2u_sbsgmn[i] << "," << F2uErr_sbsgmn[i] << ",";
    // std::cout << F2d_sbsgmn[i] << "," << F2dErr_sbsgmn[i] << "\n";
  }

  TGraphErrors *F1uTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F1u_sbsgmn, Q2Err_sbsgmn, F1uErr_sbsgmn);
  customize_sbsgmn(F1uTG_sbsgmn,"u");
  TGraphErrors *F1dTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F1d_sbsgmn, Q2Err_sbsgmn, F1dErr_sbsgmn);
  customize_sbsgmn(F1dTG_sbsgmn,"d");
  TGraphErrors *F2uTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F2u_sbsgmn, Q2Err_sbsgmn, F2uErr_sbsgmn);
  customize_sbsgmn(F2uTG_sbsgmn,"u");
  TGraphErrors *F2dTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, F2d_sbsgmn, Q2Err_sbsgmn, F2dErr_sbsgmn);  
  customize_sbsgmn(F2dTG_sbsgmn,"d");

  TGraphErrors *GMuTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, GMu_sbsgmn, Q2Err_sbsgmn, GMuErr_sbsgmn);
  customize_sbsgmn(GMuTG_sbsgmn,"u");
  TGraphErrors *GMuOVGDTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, GMuOVGD_sbsgmn, Q2Err_sbsgmn, GMuOVGDErr_sbsgmn);
  customize_sbsgmn(GMuOVGDTG_sbsgmn,"u");
  TGraphErrors *GMuOVmupGDTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, GMuOVmupGD_sbsgmn, Q2Err_sbsgmn, GMuOVmupGDErr_sbsgmn);
  customize_sbsgmn(GMuOVmupGDTG_sbsgmn,"u");  
  TGraphErrors *GMdTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, GMd_sbsgmn, Q2Err_sbsgmn, GMdErr_sbsgmn);
  customize_sbsgmn(GMdTG_sbsgmn,"d");
  TGraphErrors *GMdOVGDTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, GMdOVGD_sbsgmn, Q2Err_sbsgmn, GMdOVGDErr_sbsgmn);
  customize_sbsgmn(GMdOVGDTG_sbsgmn,"d");
  TGraphErrors *GMdOVmunGDTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, GMdOVmunGD_sbsgmn, Q2Err_sbsgmn, GMdOVmunGDErr_sbsgmn);
  customize_sbsgmn(GMdOVmunGDTG_sbsgmn,"d");  
  TGraphErrors *GEuTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, GEu_sbsgmn, Q2Err_sbsgmn, GEuErr_sbsgmn);
  customize_sbsgmn(GEuTG_sbsgmn,"u");
  TGraphErrors *GEuOVGMuTG_sbsgmn = new TGraphErrors( npoints_sbsgmn,  Q2_sbsgmn, GEuOVGMu_sbsgmn, Q2Err_sbsgmn, GEuOVGMuErr_sbsgmn);
  customize_sbsgmn(GEuOVGMuTG_sbsgmn,"u");
  TGraphErrors *mupGEuOVGMuTG_sbsgmn = new TGraphErrors( npoints_sbsgmn,  Q2_sbsgmn, mupGEuOVGMu_sbsgmn, Q2Err_sbsgmn, mupGEuOVGMuErr_sbsgmn);
  customize_sbsgmn(mupGEuOVGMuTG_sbsgmn,"u");  
  TGraphErrors *GEdTG_sbsgmn = new TGraphErrors( npoints_sbsgmn, Q2_sbsgmn, GEd_sbsgmn, Q2Err_sbsgmn, GEdErr_sbsgmn);  
  customize_sbsgmn(GEdTG_sbsgmn,"d");
  TGraphErrors *GEdOVGMdTG_sbsgmn = new TGraphErrors( npoints_sbsgmn,  Q2_sbsgmn, GEdOVGMd_sbsgmn, Q2Err_sbsgmn, GEdOVGMdErr_sbsgmn);
  customize_sbsgmn(GEdOVGMdTG_sbsgmn,"d");
  TGraphErrors *munGEdOVGMdTG_sbsgmn = new TGraphErrors( npoints_sbsgmn,  Q2_sbsgmn, munGEdOVGMd_sbsgmn, Q2Err_sbsgmn, munGEdOVGMdErr_sbsgmn);
  customize_sbsgmn(munGEdOVGMdTG_sbsgmn,"d");


  // ++++++++++++++++++++++++++++++++++++ ==============
  // *************************************
  // Cates et al
  // ++++++++++++++++++++++++++++++++++++ ==============
  // *************************************
  LookUpTableReader reader_cates;
  std::string filename_cates = "f1f2_cates.csv";
  reader_cates.readCSV(filename_cates);  
  int npoints_cates = reader_cates.getRowCount();

  double Q2_cates[npoints_cates+1];
  double Q2Err_cates[npoints_cates+1];  
  double F1u_cates[npoints_cates+1];
  double F1uErr_cates[npoints_cates+1];  
  double F2u_cates[npoints_cates+1];
  double F2uErr_cates[npoints_cates+1];  
  double F1d_cates[npoints_cates+1];
  double F1dErr_cates[npoints_cates+1];  
  double F2d_cates[npoints_cates+1];
  double F2dErr_cates[npoints_cates+1];
  //
  double GMu_cates[npoints_cates+1];
  double GMuOVGD_cates[npoints_cates+1];
  double GMuOVmupGD_cates[npoints_cates+1];  
  double GMuErr_cates[npoints_cates+1];  
  double GMuOVGDErr_cates[npoints_cates+1];
  double GMuOVmupGDErr_cates[npoints_cates+1];  
  double GEu_cates[npoints_cates+1];
  double GEuOVGMu_cates[npoints_cates+1];
  double mupGEuOVGMu_cates[npoints_cates+1];   
  double GEuErr_cates[npoints_cates+1];  
  double GEuOVGMuErr_cates[npoints_cates+1];
  double mupGEuOVGMuErr_cates[npoints_cates+1];   
  double GMd_cates[npoints_cates+1];
  double GMdOVGD_cates[npoints_cates+1];
  double GMdOVmunGD_cates[npoints_cates+1];  
  double GMdErr_cates[npoints_cates+1];  
  double GMdOVGDErr_cates[npoints_cates+1];
  double GMdOVmunGDErr_cates[npoints_cates+1];  
  double GEd_cates[npoints_cates+1];
  double GEdOVGMd_cates[npoints_cates+1];
  double munGEdOVGMd_cates[npoints_cates+1];  
  double GEdErr_cates[npoints_cates+1];
  double GEdOVGMdErr_cates[npoints_cates+1];
  double munGEdOVGMdErr_cates[npoints_cates+1];  
  
  for (int i = 0; i < npoints_cates; ++i) {
    double Q2 = reader_cates.GetValueByRowAndColumnName(i, "Q2");
    Q2_cates[i] = Q2;
    Q2Err_cates[i] = 0.;
    double Q4 = Q2_cates[i]*Q2_cates[i];
    //
    double F1u = reader_cates.GetValueByRowAndColumnName(i, "F1u");
    F1u_cates[i] = F1u*Q4;
    double F1uErr = reader_cates.GetValueByRowAndColumnName(i, "F1u_err");
    F1uErr_cates[i] = F1uErr*Q4;
    double F1d = reader_cates.GetValueByRowAndColumnName(i, "F1d");
    F1d_cates[i] = F1d*Q4*2.5;
    double F1dErr = reader_cates.GetValueByRowAndColumnName(i, "F1d_err");
    F1dErr_cates[i] = F1dErr*Q4*2.5;
    double F2u = reader_cates.GetValueByRowAndColumnName(i, "F2u");
    F2u_cates[i] = F2u*Q4/ku;
    double F2uErr = reader_cates.GetValueByRowAndColumnName(i, "F2u_err");
    F2uErr_cates[i] = F2uErr*Q4/ku;
    double F2d = reader_cates.GetValueByRowAndColumnName(i, "F2d");
    F2d_cates[i] = F2d*Q4*0.75/kd;
    double F2dErr = reader_cates.GetValueByRowAndColumnName(i, "F2d_err");
    F2dErr_cates[i] = F2dErr*Q4*0.75/abs(kd);        

    std::vector<double> GQ;
    GetGQFromFQ(Q2,F1u,F1uErr,F2u,F2uErr,F1d,F1dErr,F2d,F2dErr,0,GQ);


    GMu_cates[i] = GQ[0]; GMuErr_cates[i] = GQ[1];
    GMuOVGD_cates[i] = GQ[0]/EMFFFits::GetGDip(Q2); GMuOVGDErr_cates[i] = GQ[1]/EMFFFits::GetGDip(Q2);
    GMuOVmupGD_cates[i] = GQ[0]/(EMFFFits::GetGDip(Q2)*constant::mup); GMuOVmupGDErr_cates[i] = GQ[1]/(EMFFFits::GetGDip(Q2)*constant::mup);        
    GEu_cates[i] = GQ[2]; GEuErr_cates[i] = GQ[3];
    GEuOVGMu_cates[i] = GQ[2]/GQ[0]; GEuOVGMuErr_cates[i] = ErrPropAovB(GQ[2],GQ[3],GQ[0],GQ[1]);    
    mupGEuOVGMu_cates[i] = constant::mup*GQ[2]/GQ[0]; mupGEuOVGMuErr_cates[i] = ErrPropAovB(constant::mup*GQ[2],constant::mup*GQ[3],GQ[0],GQ[1]);        
    GMd_cates[i] = GQ[4]; GMdErr_cates[i] = GQ[5];
    GMdOVGD_cates[i] = GQ[4]/EMFFFits::GetGDip(Q2); GMdOVGDErr_cates[i] = GQ[5]/EMFFFits::GetGDip(Q2);
    GMdOVmunGD_cates[i] = GQ[4]/(EMFFFits::GetGDip(Q2)*constant::mun); GMdOVmunGDErr_cates[i] = GQ[5]/(EMFFFits::GetGDip(Q2)*constant::mun);
    GEd_cates[i] = GQ[6]; GEdErr_cates[i] = GQ[7];
    GEdOVGMd_cates[i] = GQ[6]/GQ[4]; GEdOVGMdErr_cates[i] = ErrPropAovB(GQ[6],GQ[7],GQ[4],GQ[5]);    
    munGEdOVGMd_cates[i] = constant::mun*GQ[6]/GQ[4]; munGEdOVGMdErr_cates[i] = ErrPropAovB(constant::mun*GQ[6],constant::mun*GQ[7],GQ[4],GQ[5]);    
  }

  TGraphErrors *F1uTG_cates = new TGraphErrors( npoints_cates, Q2_cates, F1u_cates, Q2Err_cates, F1uErr_cates);
  customize_cates(F1uTG_cates,"u");
  TGraphErrors *F1dTG_cates = new TGraphErrors( npoints_cates, Q2_cates, F1d_cates, Q2Err_cates, F1dErr_cates);
  customize_cates(F1dTG_cates,"d");
  TGraphErrors *F2uTG_cates = new TGraphErrors( npoints_cates, Q2_cates, F2u_cates, Q2Err_cates, F2uErr_cates);
  customize_cates(F2uTG_cates,"u");
  TGraphErrors *F2dTG_cates = new TGraphErrors( npoints_cates, Q2_cates, F2d_cates, Q2Err_cates, F2dErr_cates);  
  customize_cates(F2dTG_cates,"d");
  //
  TGraphErrors *GMuTG_cates = new TGraphErrors( npoints_cates, Q2_cates, GMu_cates, Q2Err_cates, GMuErr_cates);
  customize_cates(GMuTG_cates,"u");
  TGraphErrors *GMuOVGDTG_cates = new TGraphErrors( npoints_cates, Q2_cates, GMuOVGD_cates, Q2Err_cates, GMuOVGDErr_cates);
  customize_cates(GMuOVGDTG_cates,"u");
  TGraphErrors *GMuOVmupGDTG_cates = new TGraphErrors( npoints_cates, Q2_cates, GMuOVmupGD_cates, Q2Err_cates, GMuOVmupGDErr_cates);
  customize_cates(GMuOVmupGDTG_cates,"u");  
  TGraphErrors *GMdTG_cates = new TGraphErrors( npoints_cates, Q2_cates, GMd_cates, Q2Err_cates, GMdErr_cates);
  customize_cates(GMdTG_cates,"d");
  TGraphErrors *GMdOVGDTG_cates = new TGraphErrors( npoints_cates, Q2_cates, GMdOVGD_cates, Q2Err_cates, GMdOVGDErr_cates);
  customize_cates(GMdOVGDTG_cates,"d");
  TGraphErrors *GMdOVmunGDTG_cates = new TGraphErrors( npoints_cates, Q2_cates, GMdOVmunGD_cates, Q2Err_cates, GMdOVmunGDErr_cates);
  customize_cates(GMdOVmunGDTG_cates,"d");  
  TGraphErrors *GEuTG_cates = new TGraphErrors( npoints_cates, Q2_cates, GEu_cates, Q2Err_cates, GEuErr_cates);
  customize_cates(GEuTG_cates,"u");
  TGraphErrors *GEuOVGMuTG_cates = new TGraphErrors( npoints_cates,  Q2_cates, GEuOVGMu_cates, Q2Err_cates, GEuOVGMuErr_cates);
  customize_cates(GEuOVGMuTG_cates,"u");
  TGraphErrors *mupGEuOVGMuTG_cates = new TGraphErrors( npoints_cates,  Q2_cates, mupGEuOVGMu_cates, Q2Err_cates, mupGEuOVGMuErr_cates);
  customize_cates(mupGEuOVGMuTG_cates,"u");  
  TGraphErrors *GEdTG_cates = new TGraphErrors( npoints_cates, Q2_cates, GEd_cates, Q2Err_cates, GEdErr_cates);  
  customize_cates(GEdTG_cates,"d");
  TGraphErrors *GEdOVGMdTG_cates = new TGraphErrors( npoints_cates,  Q2_cates, GEdOVGMd_cates, Q2Err_cates, GEdOVGMdErr_cates);
  customize_cates(GEdOVGMdTG_cates,"d");
  TGraphErrors *munGEdOVGMdTG_cates = new TGraphErrors( npoints_cates,  Q2_cates, munGEdOVGMd_cates, Q2Err_cates, munGEdOVGMdErr_cates);
  customize_cates(munGEdOVGMdTG_cates,"d");  
  
  //**************

  // call the canvas customizer
  PlotCustomizer pcust{0};  
  
  // Plot F1 with SBS-GMn + Ye 2018 and Ye 2018
  TCanvas *c1 = util_pd::TC("c1",1,1);
  c1->cd(); c1->SetGridy();

  TH2D *hframe_f1 = new TH2D("hframe_f1",";Q^{2} (GeV/c)^{2};Q^{4}F_{1}^{q}",500,0,15,500,-0.6,1.6);
  customize_hframe(hframe_f1);
  
  F1uTG->Draw("C3 SAME");
  F1dTG->Draw("C3 SAME");
  F1uTG_cates->Draw("P SAME");
  F1dTG_cates->Draw("P SAME");  
  F1uTG_sbsgmn->Draw("P SAME");
  F1dTG_sbsgmn->Draw("P SAME");  
  
  TLegend *l1 = new TLegend(0.1,0.1,0.49,0.3);
  l1->SetTextFont(42);
  l1->AddEntry(F1uTG_cates,"u quark (Cates 2011)","ep");
  l1->AddEntry(F1uTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l1->AddEntry(F1uTG,"u quark (Ye2018)","lf");
  l1->AddEntry(F1dTG_cates,"d quark x 2.5 (Cates 2011","ep");
  l1->AddEntry(F1dTG_sbsgmn,"d quark x 2.5 (SBSGMn+Ye2018)","ep");
  l1->AddEntry(F1dTG,"d quark x 2.5 (Ye2018)","lf");
  l1->Draw();

  c1->SaveAs("c1_pass3.pdf");
  
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
  l2->AddEntry(F2uTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l2->AddEntry(F2uTG,"u quark (Ye2018)","lf");
  l2->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (SBSGMn+Ye2018)","ep");
  l2->AddEntry(F2dTG,"d quark x 0.75 (Ye2018)","lf");
  l2->Draw();

  c2->SaveAs("c2_pass3.pdf");
  
  //**************
  // Plot F1 with SBS-GMn + Ye 2018 and Hague + Ye 2018
  // TCanvas *c3 = util_pd::TC("c3",1,1);
  // c3->cd(); c3->SetGridy();

  TCanvas *c3 = new TCanvas("c3","c3",1500,600);
  c3->Divide(2,1);
  c3->cd(1); c3->SetGridy();

  
  TH2D *hframe_f3 = new TH2D("hframe_f3",";Q^{2} (GeV/c)^{2};Q^{4}F_{1}^{q}",500,0,15,500,-0.6,1.6);
  customize_hframe(hframe_f3);
  
  F1uTG_tyler->Draw("C3 SAME");
  F1dTG_tyler->Draw("C3 SAME");
  F1uTG_cates->Draw("P SAME");
  F1dTG_cates->Draw("P SAME");      
  F1uTG_sbsgmn->Draw("P SAME");
  F1dTG_sbsgmn->Draw("P SAME");  

  //TLegend *l3 = new TLegend(0.1,0.1,0.49,0.3);
  TLegend *l3 = new TLegend(0.15,0.15,0.54,0.35);  
  l3->SetTextFont(42);
  l3->AddEntry(F1uTG_sbsgmn,"u quark (Cates 2011)","ep");  
  l3->AddEntry(F1uTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l3->AddEntry(F1uTG_tyler,"u quark (Hague-GMn+Ye2018)","lf");
  l3->AddEntry(F1dTG_sbsgmn,"d quark x 2.5 (Cates 2011)","ep");  
  l3->AddEntry(F1dTG_sbsgmn,"d quark x 2.5 (SBSGMn+Ye2018)","ep");
  l3->AddEntry(F1dTG_tyler,"d quark x 2.5 (Hague-GMn+Ye2018)","lf");
  l3->Draw();

  c3->cd(2);
  TH2D *hframe_f3_2 = new TH2D("hframe_f3_2",";Q^{2} (GeV/c)^{2};Q^{4}F_{1}^{q}",500,0,5,500,-0.6,1.6);
  customize_hframe(hframe_f3_2);
  
  F1uTG_tyler->Draw("C3 SAME");
  F1dTG_tyler->Draw("C3 SAME");
  F1uTG_cates->Draw("P SAME");
  F1dTG_cates->Draw("P SAME");      
  F1uTG_sbsgmn->Draw("P SAME");
  F1dTG_sbsgmn->Draw("P SAME");  

  TLegend *l3_2 = new TLegend(0.15,0.15,0.54,0.35);
  l3_2->SetTextFont(42);
  l3_2->AddEntry(F1uTG_sbsgmn,"u quark (Cates 2011)","ep");  
  l3_2->AddEntry(F1uTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l3_2->AddEntry(F1uTG_tyler,"u quark (Hague-GMn+Ye2018)","lf");
  l3_2->AddEntry(F1dTG_sbsgmn,"d quark x 2.5 (Cates 2011)","ep");  
  l3_2->AddEntry(F1dTG_sbsgmn,"d quark x 2.5 (SBSGMn+Ye2018)","ep");
  l3_2->AddEntry(F1dTG_tyler,"d quark x 2.5 (Hague-GMn+Ye2018)","lf");
  l3_2->Draw();  

  pcust.customize_canvas(c3);  
  c3->SaveAs("c3_pass3.pdf");
  
  //**************
  // Plot F2 with SBS-GMn + Ye 2018 and Hague + Ye 2018  
  // TCanvas *c4 = util_pd::TC("c4",1,1);
  // c4->cd(); c4->SetGridy();

  TCanvas *c4 = new TCanvas("c4","c4",1500,600);
  c4->Divide(2,1);
  c4->cd(1); c3->SetGridy();
  

  TH2D *hframe_f4 = new TH2D("hframe_f4",";Q^{2} (GeV/c)^{2};#kappa_{q}^{-1}Q^{4}F_{2}^{q}",500,0,15,500,-0.02,0.3);
  customize_hframe(hframe_f4);
  
  F2uTG_tyler->Draw("C3 SAME");
  F2dTG_tyler->Draw("C3 SAME");
  F2uTG_cates->Draw("P SAME");
  F2dTG_cates->Draw("P SAME");    
  F2uTG_sbsgmn->Draw("P SAME");
  F2dTG_sbsgmn->Draw("P SAME");  

  //TLegend *l4 = new TLegend(0.14,0.1,0.54,0.3);
  TLegend *l4 = new TLegend(0.19,0.15,0.59,0.35);  
  l4->SetTextFont(42);
  l4->AddEntry(F2uTG_sbsgmn,"u quark (Cates 2011)","ep");  
  l4->AddEntry(F2uTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l4->AddEntry(F2uTG_tyler,"u quark (Hague-GMn+Ye2018)","lf");
  l4->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (Cates 2011)","ep");  
  l4->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (SBSGMn+Ye2018)","ep");
  l4->AddEntry(F2dTG_tyler,"d quark x 0.75 (Hague-GMn+Ye2018)","lf");
  l4->Draw();

  c4->cd(2);
  TH2D *hframe_f4_2 = new TH2D("hframe_f4_2",";Q^{2} (GeV/c)^{2};#kappa_{q}^{-1}Q^{4}F_{2}^{q}",500,0,5,500,-0.02,0.3);
  customize_hframe(hframe_f4_2);
  
  F2uTG_tyler->Draw("C3 SAME");
  F2dTG_tyler->Draw("C3 SAME");
  F2uTG_cates->Draw("P SAME");
  F2dTG_cates->Draw("P SAME");    
  F2uTG_sbsgmn->Draw("P SAME");
  F2dTG_sbsgmn->Draw("P SAME");  

  TLegend *l4_2 = new TLegend(0.24,0.15,0.64,0.35);
  l4_2->SetTextFont(42);
  l4_2->AddEntry(F2uTG_sbsgmn,"u quark (Cates 2011)","ep");  
  l4_2->AddEntry(F2uTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l4_2->AddEntry(F2uTG_tyler,"u quark (Hague-GMn+Ye2018)","lf");
  l4_2->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (Cates 2011)","ep");  
  l4_2->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (SBSGMn+Ye2018)","ep");
  l4_2->AddEntry(F2dTG_tyler,"d quark x 0.75 (Hague-GMn+Ye2018)","lf");
  l4_2->Draw();  

  pcust.customize_canvas(c4);    
  c4->SaveAs("c4_pass3.pdf");  

  // -------- ###############
  // Plotting GE and GM
  // --------
  TCanvas *c5 = util_pd::TC("c5",1,2);
  //**************

  // Plot GMu with SBS-GMn + Ye 2018 and Ye 2018
  c5->cd(1); gPad->SetGridy();
  TH2D *hframe_f5_1 = new TH2D("hframe_f5_1",";Q^{2} (GeV/c)^{2};G_{M}^{u}/G_{D}",500,0,15,500,3,5);
  customize_hframe(hframe_f5_1);
  GMuOVGDTG->Draw("C3 SAME");
  GMuOVGDTG_sbsgmn->Draw("P SAME");

  //
  TLegend *l5_1 = new TLegend(0.3,0.75,0.95,0.9); //0.15,0.15,0.80,0.3);
  l5_1->SetTextFont(42);
  l5_1->AddEntry(GMuOVmupGDTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l5_1->AddEntry(GMuOVmupGDTG,"u quark (Ye2018)","lf");
  l5_1->Draw();

  // Plot GMd with SBS-GMn + Ye 2018 and Ye 2018
  c5->cd(2); gPad->SetGridy();
  TH2D *hframe_f5_2 = new TH2D("hframe_f5_2",";Q^{2} (GeV/c)^{2};G_{M}^{d}/G_{D}",500,0,15,500,-1.5,1.5);
  customize_hframe(hframe_f5_2);  
  GMdOVGDTG->Draw("C3 SAME");
  GMdOVGDTG_sbsgmn->Draw("P SAME");  
  //
  TLegend *l5_2 = new TLegend(0.3,0.75,0.95,0.9);
  l5_2->SetTextFont(42);
  l5_2->AddEntry(GMdOVmunGDTG_sbsgmn,"d quark (SBSGMn+Ye2018)","ep");
  l5_2->AddEntry(GMdOVmunGDTG,"d quark (Ye2018)","lf");
  l5_2->Draw();

  // // Plot mupGEu/GMu with SBS-GMn + Ye2018 and Ye2018
  // c5->cd(3); gPad->SetGridy();
  // TH2D *hframe_f5_3 = new TH2D("hframe_f5_3",";Q^{2} (GeV/c)^{2};G_{E}^{u}/G_{M}^{u}",500,0,15,500,-0.1,0.8);
  // customize_hframe(hframe_f5_3);
  // GEuOVGMuTG->Draw("C3 SAME");
  // GEuOVGMuTG_sbsgmn->Draw("P SAME");  
  // //
  // TLegend *l5_3 = new TLegend(0.3,0.75,0.95,0.9);
  // l5_3->SetTextFont(42);
  // l5_3->AddEntry(mupGEuOVGMuTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  // l5_3->AddEntry(mupGEuOVGMuTG,"u quark (Ye2018)","lf");
  // l5_3->Draw();

  // // Plot munGEd/GMd with SBS-GMn + Ye2018 and Ye2018
  // c5->cd(4); gPad->SetGridy();
  // TH2D *hframe_f5_4 = new TH2D("hframe_f5_4",";Q^{2} (GeV/c)^{2};G_{E}^{d}/G_{M}^{d}",500,0,15,500,-10,3);
  // customize_hframe(hframe_f5_4);
  // GEdOVGMdTG->Draw("C3 SAME");
  // GEdOVGMdTG_sbsgmn->Draw("P SAME");

  // //
  // TLegend *l5_4 = new TLegend(0.3,0.75,0.95,0.9);
  // l5_4->SetTextFont(42);
  // l5_4->AddEntry(munGEdOVGMdTG_sbsgmn,"d quark (SBSGMn+Ye2018)","ep");
  // l5_4->AddEntry(munGEdOVGMdTG,"d quark (Ye2018)","lf");
  // l5_4->Draw();  

  pcust.customize_canvas(c5);
  c5->SaveAs("c5_pass3.pdf");

  
  TCanvas *c6 = util_pd::TC("c6",2,2);
  //**************

  // Plot GMu with SBS-GMn + Ye2018 and Ye2018
  c6->cd(1); gPad->SetGridy();
  TH2D *hframe_f6_1 = new TH2D("hframe_f6_1",";Q^{2} (GeV/c)^{2};G_{M}^{u}/(#mu_{p}G_{D})",500,0,15,500,1,1.7);
  customize_hframe(hframe_f6_1);
  GMuOVmupGDTG->Draw("C3 SAME");
  GMuOVmupGDTG_sbsgmn->Draw("P SAME");

  //
  TLegend *l6_1 = new TLegend(0.3,0.75,0.95,0.9); //0.15,0.15,0.80,0.3);
  l6_1->SetTextFont(42);
  l6_1->AddEntry(GMuOVmupGDTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l6_1->AddEntry(GMuOVmupGDTG,"u quark (Ye2018)","lf");
  l6_1->Draw();

  // Plot GMd with SBS-GMn + Ye2018 and Ye2018
  c6->cd(2); gPad->SetGridy();
  TH2D *hframe_f6_2 = new TH2D("hframe_f6_2",";Q^{2} (GeV/c)^{2};G_{M}^{d}/(#mu_{n}G_{D})",500,0,15,500,0,0.6);
  customize_hframe(hframe_f6_2);  
  GMdOVmunGDTG->Draw("C3 SAME");
  GMdOVmunGDTG_sbsgmn->Draw("P SAME");
  //
  TLegend *l6_2 = new TLegend(0.3,0.75,0.95,0.9);
  l6_2->SetTextFont(42);
  l6_2->AddEntry(GMdOVmunGDTG_sbsgmn,"d quark (SBSGMn+Ye2018)","ep");
  l6_2->AddEntry(GMdOVmunGDTG,"d quark (Ye2018)","lf");
  l6_2->Draw();

  // Plot mupGEu/GMu with SBS-GMn + Ye2018 and Ye2018
  c6->cd(3); gPad->SetGridy();
  TH2D *hframe_f6_3 = new TH2D("hframe_f6_3",";Q^{2} (GeV/c)^{2};#mu_{p}G_{E}^{u}/G_{M}^{u}",500,0,15,500,0,1.6);
  customize_hframe(hframe_f6_3);
  mupGEuOVGMuTG->Draw("C3 SAME");
  mupGEuOVGMuTG_sbsgmn->Draw("P SAME");
  //
  TLegend *l6_3 = new TLegend(0.3,0.75,0.95,0.9);
  l6_3->SetTextFont(42);
  l6_3->AddEntry(mupGEuOVGMuTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l6_3->AddEntry(mupGEuOVGMuTG,"u quark (Ye2018)","lf");
  l6_3->Draw();

  // Plot munGEd/GMd with SBS-GMn + Ye2018 and Ye2018
  c6->cd(4); gPad->SetGridy();
  TH2D *hframe_f6_4 = new TH2D("hframe_f6_4",";Q^{2} (GeV/c)^{2};#mu_{n}G_{E}^{d}/G_{M}^{d}",500,0,15,500,0,15);
  customize_hframe(hframe_f6_4);
  munGEdOVGMdTG->Draw("C3 SAME");
  munGEdOVGMdTG_sbsgmn->Draw("P SAME");

  //
  TLegend *l6_4 = new TLegend(0.3,0.75,0.95,0.9);
  l6_4->SetTextFont(42);
  l6_4->AddEntry(munGEdOVGMdTG_sbsgmn,"d quark (SBSGMn+Ye2018)","ep");
  l6_4->AddEntry(munGEdOVGMdTG,"d quark (Ye2018)","lf");
  l6_4->Draw();  

  pcust.customize_canvas(c6);
  c6->SaveAs("c6_pass3.pdf");


  // sanity checks
  // 1. let's reproduce GMn from GMu and GMd
  // Plot F2 with SBS-GMn + Ye2018 and Hague + Ye2018  
  TCanvas *csc1 = util_pd::TC("csc1",1,1);
  csc1->cd(); csc1->SetGridy();

  TH2D *hframe_fsc1 = new TH2D("hframe_fsc1",";Q^{2} (GeV/c)^{2};G_{M}^{n}/(#muG_{D})",500,0,15,500,0,1.2);
  customize_hframe(hframe_fsc1);

  GMnOVmunGDSCTG->Draw("C3 SAME");  

  // TLegend *lsc1 = new TLegend(0.14,0.1,0.54,0.3);
  // lsc1->SetTextFont(42);
  // lsc1->AddEntry(F2uTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  // lsc1->AddEntry(F2uTG_tyler,"u quark (Hague-GMn+Ye2018)","lf");
  // lsc1->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (SBSGMn+Ye2018)","ep");
  // lsc1->AddEntry(F2dTG_tyler,"d quark x 0.75 (Hague-GMn+Ye2018)","lf");
  // lsc1->Draw();


  // -------- ###############
  // Plotting GM - Hague+Ye
  // --------
  TCanvas *c7 = util_pd::TC("c7",1,1);
  //**************

  // Plot GMu with SBS-GMn + Ye 2018 and Ye 2018
  c7->cd(1); gPad->SetGridy();
  TH2D *hframe_f7 = new TH2D("hframe_f7",";Q^{2} (GeV/c)^{2};G_{M}^{u}/G_{D}",500,0,15,500,3,5);
  customize_hframe(hframe_f7);
  GMuOVGDTG_tyler->Draw("C3 SAME");
  GMuOVGDTG_cates->Draw("P SAME");  
  GMuOVGDTG_sbsgmn->Draw("P SAME");

  //
  TLegend *l7 = new TLegend(0.3,0.75,0.95,0.9); //0.15,0.15,0.80,0.3);
  l7->SetTextFont(42);
  l7->AddEntry(GMuOVmupGDTG_cates,"u quark (Cates 2011)","ep");  
  l7->AddEntry(GMuOVmupGDTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l7->AddEntry(GMuOVmupGDTG_tyler,"u quark (Hague-GMn+Ye2018)","lf");
  l7->Draw();

  pcust.customize_canvas(c7);
  c7->SaveAs("c7_pass3.pdf");
  
  // -------- ###############
  // Plotting GE - Hague+Ye
  // --------
  TCanvas *c8 = util_pd::TC("c8",1,1);
  //**************

  // Plot GMu with SBS-GMn + Ye 2018 and Ye 2018
  c8->cd(1); gPad->SetGridy();
  TH2D *hframe_f8 = new TH2D("hframe_f8",";Q^{2} (GeV/c)^{2};G_{M}^{d}/G_{D}",500,0,15,500,-1.5,1.5);
  customize_hframe(hframe_f8);
  GMdOVGDTG_tyler->Draw("C3 SAME");
  GMdOVGDTG_cates->Draw("P SAME");  
  GMdOVGDTG_sbsgmn->Draw("P SAME");

  //
  TLegend *l8 = new TLegend(0.3,0.75,0.95,0.9); //0.15,0.15,0.80,0.3);
  l8->SetTextFont(42);
  l8->AddEntry(GMdOVGDTG_cates,"d quark (Cates 2011)","ep");  
  l8->AddEntry(GMdOVGDTG_sbsgmn,"d quark (SBSGMn+Ye2018)","ep");
  l8->AddEntry(GMdOVGDTG_tyler,"d quark (Hague-GMn+Ye2018)","lf");
  l8->Draw();

  pcust.customize_canvas(c8);
  c8->SaveAs("c8_pass3.pdf");
  
  // -------- ###############
  // Plotting GE and GM - Hague+Ye
  // --------
  TCanvas *c9 = util_pd::TC("c9",1,1);
  //**************

  // Plot GMu with SBS-GMn + Ye 2018 and Ye 2018
  c9->cd(1); gPad->SetGridy();
  TH2D *hframe_f9 = new TH2D("hframe_f9",";Q^{2} (GeV/c)^{2};G_{E}^{u}/G_{M}^{u}",500,0,15,500,-0.1,0.8);
  customize_hframe(hframe_f9);
  GEuOVGMuTG_tyler->Draw("C3 SAME");
  GEuOVGMuTG_cates->Draw("P SAME");  
  GEuOVGMuTG_sbsgmn->Draw("P SAME");

  //
  TLegend *l9 = new TLegend(0.3,0.75,0.95,0.9); //0.15,0.15,0.80,0.3);
  l9->SetTextFont(42);
  l9->AddEntry(GEuOVGMuTG_cates,"u quark (Cates 2011)","ep");  
  l9->AddEntry(GEuOVGMuTG_sbsgmn,"u quark (SBSGMn+Ye2018)","ep");
  l9->AddEntry(GEuOVGMuTG_tyler,"u quark (Hague-GMn+Ye2018)","lf");
  l9->Draw();

  pcust.customize_canvas(c9);
  c9->SaveAs("c9_pass3.pdf");

  // --------
  TCanvas *c10 = util_pd::TC("c10",1,1);
  //**************  

  // Plot GMu with SBS-GMn + Ye 2018 and Ye 2018
  c10->cd(1); gPad->SetGridy();
  TH2D *hframe_f10 = new TH2D("hframe_f10",";Q^{2} (GeV/c)^{2};G_{E}^{d}/G_{M}^{d}",500,0,15,500,-10,3);
  customize_hframe(hframe_f10);
  GEdOVGMdTG_tyler->Draw("C3 SAME");
  GEdOVGMdTG_cates->Draw("P SAME");  
  GEdOVGMdTG_sbsgmn->Draw("P SAME");

  //
  TLegend *l10 = new TLegend(0.3,0.75,0.95,0.9); //0.15,0.15,0.80,0.3);
  l10->SetTextFont(42);
  l10->AddEntry(GEdOVGMdTG_cates,"d quark (Cates 2011)","ep");  
  l10->AddEntry(GEdOVGMdTG_sbsgmn,"d quark (SBSGMn+Ye2018)","ep");
  l10->AddEntry(GEdOVGMdTG_tyler,"d quark (Hague-GMn+Ye2018)","lf");
  l10->Draw();

  pcust.customize_canvas(c10);
  c10->SaveAs("c10_pass3.pdf");  

  // // Cates F1, F2 + SBS F1, F2
  // //**************
  // // Plot F1 with SBS-GMn + Ye 2018 and Ye 2018
  // TCanvas *ccatesf1 = util_pd::TC("ccatesf1",1,1);
  // ccatesf1->cd(); ccatesf1->SetGridy();

  // TH2D *hframe_fcatesf1 = new TH2D("hframe_fcatesf1",";Q^{2} (GeV/c)^{2};Q^{4}F_{1}^{q}",500,0,15,500,-0.6,1.6);
  // customize_hframe(hframe_fcatesf1);
  
  // F1uTG_cates->Draw("P SAME");
  // F1dTG_cates->Draw("P SAME");
  // F1uTG_sbsgmn->Draw("P SAME");
  // F1dTG_sbsgmn->Draw("P SAME");  

  // TLegend *lcatesf1 = new TLegend(0.1,0.1,0.49,0.3);
  // lcatesf1->SetTextFont(42);
  // lcatesf1->AddEntry(F1uTG_sbsgmn,"u quark (SBSGMn)","ep");
  // lcatesf1->AddEntry(F1uTG_cates,"u quark (Cates 2011)","ep");
  // lcatesf1->AddEntry(F1dTG_sbsgmn,"d quark x 2.5 (SBSGMn)","ep");
  // lcatesf1->AddEntry(F1dTG_cates,"d quark x 2.5 (Cates 2011)","ep");
  // lcatesf1->Draw();

  // ccatesf1->SaveAs("ccatesf1_pass3_pass3.pdf");
  
  // //**************
  // // Plot F2 with SBS-GMn + Ye 2018 and Ye 2018  
  // TCanvas *ccatesf2 = util_pd::TC("ccatesf2",1,1);
  // ccatesf2->cd(); ccatesf2->SetGridy();

  // TH2D *hframe_fcatesf2 = new TH2D("hframe_fcatesf2",";Q^{2} (GeV/c)^{2};#kappa_{q}^{-1}Q^{4}F_{2}^{q}",500,0,15,500,-0.02,0.3);
  // customize_hframe(hframe_fcatesf2);
  
  // F2uTG_cates->Draw("P SAME");
  // F2dTG_cates->Draw("P SAME");    
  // F2uTG_sbsgmn->Draw("P SAME");
  // F2dTG_sbsgmn->Draw("P SAME");  

  // TLegend *lcatesf2 = new TLegend(0.14,0.1,0.54,0.3);
  // lcatesf2->SetTextFont(42);
  // lcatesf2->AddEntry(F2uTG_sbsgmn,"u quark (SBSGMn)","ep");
  // lcatesf2->AddEntry(F2uTG_cates,"u quark (Cates 2011)","ep");
  // lcatesf2->AddEntry(F2dTG_sbsgmn,"d quark x 0.75 (SBSGMn)","ep");
  // lcatesf2->AddEntry(F2dTG_cates,"d quark x 0.75 (Cates 2011)","ep");
  // lcatesf2->Draw();

  // ccatesf2->SaveAs("ccatesf2_pass3.pdf");  

  
  return 0;
}
