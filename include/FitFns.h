#ifndef FIT_FNS_H
#define FIT_FNS_H

class FitFn {
 private:
  int fpoly=0;       // order of polynomial
  TH1D *fhs1=NULL;   // Signal histo 1
  TH1D *fhs2=NULL;   // Signal histo 2
  TH1D *fhbg=NULL;   // Background histo

 public:
  FitFn(TH1D *hs1): fhs1(hs1) {}
  FitFn(TH1D *hs1, int poly): fhs1(hs1),fpoly(poly) {}
  FitFn(TH1D *hs1, TH1D *hs2): fhs1(hs1),fhs2(hs2) {}
  FitFn(TH1D *hs1, TH1D *hs2, TH1D *hbg): fhs1(hs1),fhs2(hs2),fhbg(hbg) {}
  FitFn(TH1D *hs1, TH1D *hs2, int poly): fhs1(hs1),fhs2(hs2),fpoly(poly) {}

  // fits using just 1 signal histo (1 param)
  double ffn_1hs_nbg (double *x, double *par) const {
    double Norm = par[0];
    return Norm*fhs1->Interpolate(x[0]);
  }

  // fits using 1 signal histo and 1 poly bg (1+poly params)
  double ffn_1hs_1pbg (double *x, double *par) const {
    double Norm = par[0];
    double bg = 0; for (int i=1; i<fpoly+2; i++) bg += par[i]*pow(x[0],i-1);
    return Norm*fhs1->Interpolate(x[0]) + bg;
  }

  // fits using just 2 signal histos (2 params)
  double ffn_2hs_nbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0]));
  }

  // fits using 2 signal histos and 1 bg histo (3 params)
  double ffn_2hs_1hbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double B = par[2];
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0])+B*fhbg->Interpolate(x[0]));
  }

  // fits using 2 signal histos and 1 poly bg (2+poly params)
  double ffn_2hs_1pbg (double *x, double *par) const {
    double Norm = par[0];
    double R = par[1];
    double bg = 0; for (int i=2; i<fpoly+3; i++) bg += par[i]*pow(x[0],i-2);
    return Norm*(fhs1->Interpolate(x[0])+R*fhs2->Interpolate(x[0])) + bg;
  }
};

#endif
