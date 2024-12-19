#ifndef FIT_FNS_H
#define FIT_FNS_H

#include "TH1F.h"

class FitFn {
 private:
  int fpoly{0};          // order of polynomial
  double frp1{0};        // 1st reject point (lower limit)
  double frp2{0};        // 2nd reject point (upper limit)
  TH1F *fhs1{nullptr};   // Signal histo 1
  TH1F *fhs2{nullptr};   // Signal histo 2
  TH1F *fhbg{nullptr};   // Background histo
  TH1F *fhbg1{nullptr};  // Background histo 1
  TH1F *fhbg2{nullptr};  // Background histo 2
  TH1F *fhbg3{nullptr};  // Background histo 3

 public:
  FitFn() {}
  FitFn(TH1F *hs1): fhs1(hs1) {}
  FitFn(int poly): fpoly(poly) {}
  FitFn(int poly, double rp1, double rp2): fpoly(poly),frp1(rp1),frp2(rp2) {}
  FitFn(TH1F *hs1, int poly): fhs1(hs1),fpoly(poly) {}
  FitFn(TH1F *hs1, TH1F *hs2): fhs1(hs1),fhs2(hs2) {}
  FitFn(TH1F *hs1, TH1F *hs2, TH1F *hbg): fhs1(hs1),fhs2(hs2),fhbg(hbg) {}
  FitFn(TH1F *hs1, TH1F *hs2, TH1F *hbg1, TH1F *hbg2): fhs1(hs1),fhs2(hs2),fhbg1(hbg1),fhbg2(hbg2) {}
  FitFn(TH1F *hs1, TH1F *hs2, TH1F *hbg1, TH1F *hbg2, TH1F *hbg3): fhs1(hs1),fhs2(hs2),fhbg1(hbg1),fhbg2(hbg2),fhbg3(hbg3) {}
  FitFn(TH1F *hs1, TH1F *hs2, int poly): fhs1(hs1),fhs2(hs2),fpoly(poly) {}

  // returns Gaussian fit function
  double ffn_gaus (double *x, double *par) const {
    return par[0]*std::exp(-0.5*std::pow((x[0]-par[1])/par[2],2.));
  }

  // returns polynomial fit function of order fpoly
  double ffn_poly (double *x, double *par) const {
    double poly = 0; for (int i=0; i<fpoly+1; i++) poly += par[i]*pow(x[0],i);
    return poly;
  }

  // double ffn_poly (double *x, double *par) const {
  //   double poly = 0;
  //   // Ensure the highest-order term has a negative coefficient to enforce concavity (concave down)
  //   poly += -fabs(par[fpoly]) * pow(x[0], fpoly);
  //   for (int i = 0; i < fpoly; i++) {
  //       poly += par[i] * pow(x[0], i);
  //   }
  //   return poly;
  // }

  // returns polynomial fit fn of order fpoly considering the reject points (Sideband method)
  double ffn_1pbg_sb_2rp(double *x, double *par) const {
    if (x[0]>frp1 && x[0]<frp2){
      TF1::RejectPoint();
      return 0;
    }
    return ffn_poly(x,par);
  }

  // fits using just 1 signal histo (1 param)
  double ffn_1hs_nbg (double *x, double *par) const {
    double Norm = par[0];
    return Norm*fhs1->Interpolate(x[0]);
  }

  // fits using 1 signal histo and 1 poly bg (1+fpoly+1 params)
  double ffn_1hs_1pbg (double *x, double *par) const {
    double Norm = par[0];
    return Norm*fhs1->Interpolate(x[0]) + ffn_poly(x,&par[1]);
  }

  // fits using 1 signal histo (vary x offset) and 1 poly bg (2+fpoly+1 params)
  double ffn_1hs_1pbg_xOffVary (double *x, double *par) const {
    double Norm = par[0];
    double hsxOff = par[1];
    //
    double modx = x[0]-hsxOff;
    return Norm*fhs1->Interpolate(modx) + ffn_poly(x,&par[1]);
  }

  // fits using just 2 signal histos (2 params)
  double ffn_2hs_nbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0]));
  }

  // fits using just 2 signal histos (vary x offsets for both) (2+2 params)
  double ffn_2hs_nbg_xOffVary (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double hs1xOff = par[2];
    double hs2xOff = par[3];
    //
    double modx1 = x[0]-hs1xOff;
    double modx2 = x[0]-hs2xOff;
    return Norm*(fhs1->Interpolate(modx1)+R*fhs2->Interpolate(modx2));
    // return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(modx1));
  }

  // fits using 2 signal histos and 1 bg histo (3 params)
  double ffn_2hs_1hbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double B = par[2];
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0]))+B*fhbg->Interpolate(x[0]);
  }

  // fits using 2 signal histos (vary x offsets for both) and 1 bg histo (5 params)
  double ffn_2hs_1hbg_xOffVary (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double B = par[2];
    double hs1xOff = par[3];
    double hs2xOff = par[4];
    //
    double modx1 = x[0]-hs1xOff;
    double modx2 = x[0]-hs2xOff;
    return Norm*(fhs1->Interpolate(modx1)+R*fhs2->Interpolate(modx2)+B*fhbg->Interpolate(x[0]));
  }

  // fits using 2 signal histos and 2 bg histo (3 params)
  double ffn_2hs_2hbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double B = par[2];
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0]))+B*(fhbg1->Interpolate(x[0])+fhbg2->Interpolate(x[0]));
  }

  // fits using 2 signal histos (vary x offsets for both) and 2 bg histo (5 params)
  double ffn_2hs_2hbg_xOffVary (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double B = par[2];
    double hs1xOff = par[3];
    double hs2xOff = par[4];
    //
    double modx1 = x[0]-hs1xOff;
    double modx2 = x[0]-hs2xOff;
    return Norm*(fhs1->Interpolate(modx1)+R*fhs2->Interpolate(modx2))+B*(fhbg1->Interpolate(modx1)+fhbg2->Interpolate(modx2));
  }

  // fits using 2 signal histos and 3 bg histos. B2 belongs to the third hbg (4 params)
  double ffn_2hs_3hbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double B1 = par[2];
    double B2 = par[3];    
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0]))+B1*(fhbg1->Interpolate(x[0])+fhbg2->Interpolate(x[0]))+B2*fhbg3->Interpolate(x[0]);
  }

  // fits using 2 signal histos and 3 bg histos (6 params)
  // x offsets are varied for both signal histos and the first and second hbg
  // B2 belongs to the third hbg
  double ffn_2hs_3hbg_xOffVary (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double B1 = par[2];
    double hs1xOff = par[3];
    double hs2xOff = par[4];
    double B2 = par[5];    
    //
    double modx1 = x[0]-hs1xOff;
    double modx2 = x[0]-hs2xOff;
    return Norm*(fhs1->Interpolate(modx1)+R*fhs2->Interpolate(modx2))+B1*(fhbg1->Interpolate(modx1)+fhbg2->Interpolate(modx2))+B2*fhbg3->Interpolate(x[0]);
  }  

  // fits using 2 signal histos and 1 poly bg (2+fpoly+1 params)
  double ffn_2hs_1pbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0])) + ffn_poly(x,&par[2]);
  }

  // fits using 2 signal (vary x offsets for both) histos and 1 poly bg (4+fpoly+1 params)
  double ffn_2hs_1pbg_xOffVary (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double hs1xOff = par[2];
    double hs2xOff = par[3];
    //
    double modx1 = x[0]-hs1xOff;
    double modx2 = x[0]-hs2xOff;
    return Norm*(fhs1->Interpolate(modx1)+R*fhs2->Interpolate(modx2)) + ffn_poly(x,&par[4]);
  }

  // fits using 2 signal histos and 1 Gaussian bg (2+3 params)
  double ffn_2hs_1gbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0])) + ffn_gaus(x,&par[2]);
  }

  // fits using 2 signal (vary x offsets for both) histos and 1 poly bg (4+fpoly+1 params)
  double ffn_2hs_1gbg_xOffVary (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double hs1xOff = par[2];
    double hs2xOff = par[3];
    //
    double modx1 = x[0]-hs1xOff;
    double modx2 = x[0]-hs2xOff;
    return Norm*(fhs1->Interpolate(modx1)+R*fhs2->Interpolate(modx2)) + ffn_gaus(x,&par[4]);
  }
};

#endif
