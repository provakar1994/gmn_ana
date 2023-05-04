#include "../include/Fit.h"

namespace fit {

  // TF1* fit_1hs_nbg_THI (vector<double> fit_range,
  // 		    TH1D* ht,           // total histo to fit 
  // 		    TH1D* hs,           // signal histo for
  // 		    vector<TH1D*> &ho) // Output: ht,hs (fitted w/ proper scaling)
  // /* TH Interpolation fit using 1 signal histo & no background */
  // {
    
  // }

  TF1* fit_2hs_nbg_THI (vector<double> fit_range,
			TH1D* ht,             // total histo to fit 
			TH1D* hs1,            // 1st signal histo for fit
			TH1D* hs2,            // 2nd signal histo for fit
			vector<TH1D*> &ho)    // Output: ht,hs1,hs2,hst (fitted w/ proper scaling)
  /* TH Interpolation fit using 2 signal histo & no background */
  {
    const int npars = 2;
    vector<double> setpars{1,1};

    TH1D *ht_cp = (TH1D*)ht->Clone(); 
    TH1D *hs1_cp1 = (TH1D*)hs1->Clone(); 
    TH1D *hs2_cp1 = (TH1D*)hs2->Clone();

    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_nbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");

    ht_cp->Fit(f1); 
    vector<double> pars;
    for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1D *hs1_cp2 = (TH1D*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1D *hs2_cp2 = (TH1D*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1D *hst = (TH1D*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    ho = {ht_cp,hst,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_2hs_1hbg_THI (vector<double> fit_range,
			 TH1D* ht,             // total histo to fit 
			 TH1D* hs1,            // 1st signal histo for fit
			 TH1D* hs2,            // 2nd signal histo for fit
			 TH1D* hbg,            // bg histo for fit
			 vector<TH1D*> &ho)    // Output: ht,hs1,hs2,hst (fitted w/ proper scaling)
  /* TH Interpolation fit using 2 signal histo & 1 bg histo */
  {
    const int npars = 3;
    vector<double> setpars{1,1,0};

    TH1D *ht_cp = (TH1D*)ht->Clone(); 
    TH1D *hs1_cp1 = (TH1D*)hs1->Clone(); 
    TH1D *hs2_cp1 = (TH1D*)hs2->Clone();
    TH1D *hbg_cp1 = (TH1D*)hbg->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1,hbg_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_1hbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    f1->SetParName(2,"B");

    ht_cp->Fit(f1);
    vector<double> pars;
    for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1D *hs1_cp2 = (TH1D*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1D *hs2_cp2 = (TH1D*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1D *hst = (TH1D*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1D *hbg_cp2 = (TH1D*)hs1_cp1->Clone(); hbg_cp2->Add(ht_cp,hst,1,-1); 
    ho = {ht_cp,hst,hbg_cp2,hs1_cp2,hs2_cp2};
    
    return f1;
  }

}
