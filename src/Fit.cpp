#include "../include/Fit.h"

namespace fit {

  void set_poly_par_names (TF1* f1,
			   int const Spar,
			   int const Opoly)
  /* Sets the name of polynomial parameters */
  {
    for (int i=Spar; i<Opoly+Spar+1; i++) {
      std::string s = "p" + std::to_string(i-Spar);
      f1->SetParName(i,s.c_str());
    }
  }

  void set_gaus_par_names (TF1* f1,
			   int const Spar)
  /* Sets the name of Gaussian fit parameters */
  {
    f1->SetParName(Spar,"Const");
    f1->SetParName(Spar+1,"Mean");
    f1->SetParName(Spar+2,"Sigma");
  }

  std::vector<double> GetFitParams(TF1 * const f1)
  /* Returns a vector filled with fit pramater values from f1 */
  {
    int npars = f1->GetNumberFreeParameters();
    std::vector<double> pars;
    for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}
    return pars;
  }

  std::vector<double> GetFitParamErrors(TF1 * const f1)
  /* Returns a vector filled with fit pramater errors from f1 */
  {
    int npars = f1->GetNumberFreeParameters();
    std::vector<double> parerrs;
    for (int i=0;i<npars;i++) {parerrs.push_back(f1->GetParError(i));}
    return parerrs;
  }

  TF1* fit_1pbg_SB (std::vector<double> const & fit_range,
		    std::vector<double> const & reject_points,
		    int Opoly,          // Order of poly to fit bg
		    std::vector<double> const & initial_guesses,
		    TH1F* ht,           // total histo to fit
		    std::vector<TH1F*> &ho)    // Output: ht,hs,hbg
  /* Side band (SB) fit using 1 polynomial and 2 reject points (Opoly+1 pars) */
  /* NOTE: This fit method just fits the bg */
  {
    const int npars = Opoly+1;

    TH1F *ht_cp = (TH1F*)ht->Clone(); 

    FitFn *ffn = new FitFn(Opoly,reject_points[0],reject_points[1]);
    TF1* f1 = new TF1("f1",ffn,&FitFn::ffn_1pbg_sb_2rp,fit_range[0],fit_range[1],Opoly+1);
    f1->SetNpx(1000);
    f1->SetParameters(&initial_guesses[0]);

    ht_cp->Fit(f1,"R");

    // drawing the bg
    TF1* bg = new TF1("bg",ffn,&FitFn::ffn_poly,fit_range[0],fit_range[1],Opoly+1);
    bg->SetNpx(500);
    // bg->SetParameters(&pars[0]);
    bg->SetParameters(&GetFitParams(f1)[0]);
    bg->SetLineColor(kGreen+2);

    // let's extract the signal histo
    TH1F *hbg = (TH1F*)ht_cp->Clone(); util_pd::TF1toTH1F(bg,hbg);
    TH1F *hs = (TH1F*)ht_cp->Clone(); hs->Add(ht_cp,hbg,1,-1);
    ho = {ht_cp,hs,hbg};

    return bg;
  }

  TF1* fit_1hs_nbg_THI (std::vector<double> const & fit_range,
			TH1F* ht,           // total histo to fit 
			TH1F* hs,           // signal histo for
			std::vector<TH1F*> &ho)  // Output: ht,hs (fitted w/ proper scaling)
  /* TH Interpolation fit using 1 signal histo & no background (1 par) */
  {
    const int npars = 1;
    std::vector<double> setpars{1};

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs_cp1 = (TH1F*)hs->Clone(); 

    FitFn *ffn = new FitFn(hs_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_1hs_nbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");

    ht_cp->Fit(f1,"R"); 
    // std::vector<double> pars;
    // for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1F *hs_cp2 = (TH1F*)hs_cp1->Clone(); hs_cp2->Scale(GetFitParams(f1)[0]); //hs_cp2->Scale(pars[0]);
    ho = {ht_cp,hs_cp2};
    
    return f1;    
  }

  TF1* fit_1hs_1pbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs,             // 1st signal histo for fit
			 int Opoly,            // Order of poly to fit bg
			 std::vector<TH1F*> &ho)    // Output: ht,hs,hbg
  /* TH Interpolation fit using 1 signal histo & 1 poly bg (1+Opoly+1 pars) */
  {
    const int npars = 1+Opoly+1;
    std::vector<double> setpars{1}; for (int i=1;i<npars;i++) setpars.push_back(0);

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs_cp1 = (TH1F*)hs->Clone(); 
 
    FitFn *ffn = new FitFn(hs_cp1,Opoly);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_1hs_1pbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    set_poly_par_names(f1,1,Opoly);

    ht_cp->Fit(f1,"R");
    // std::vector<double> pars;
    // for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1F *hs_cp2 = (TH1F*)hs_cp1->Clone(); hs_cp2->Scale(GetFitParams(f1)[0]); //hs_cp2->Scale(pars[0]);
    TH1F *hbg = (TH1F*)ht_cp->Clone(); hbg->Add(ht_cp,hs_cp2,1,-1); 
    ho = {ht_cp,hs_cp2,hbg};
    
    return f1;
  }

  TF1* fit_1hs_1hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs,             // signal histo for fit
			 TH1F* hbg,            // bg histo for fit
			 std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,N*hs (N=par[0],B=par[1])
  /* TH Interpolation fit using 1 signal histo & 1 bg histo (2 pars) */
  {
    const int npars = 2;
    std::vector<double> setpars{1,0};

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs_cp1 = (TH1F*)hs->Clone(); 
    TH1F *hbg_cp1 = (TH1F*)hbg->Clone();
 
    FitFn *ffn = new FitFn(hs_cp1,hbg_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_nbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"B");

    ht_cp->Fit(f1,"R");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1F *hs_cp2 = (TH1F*)hs_cp1->Clone(); hs_cp2->Scale(pars[0]);
    TH1F *hst = (TH1F*)hs_cp2->Clone();
    TH1F *hbg_cp2 = (TH1F*)hbg_cp1->Clone(); hbg_cp2->Scale(pars[0]*pars[1]); 
    ho = {ht_cp,hst,hbg_cp2,hs_cp2};
    
    return f1;
  }

  TF1* fit_2hs_nbg_THI (std::vector<double> const & fit_range,
			TH1F* ht,             // total histo to fit 
			TH1F* hs1,            // 1st signal histo for fit
			TH1F* hs2,            // 2nd signal histo for fit
			std::vector<TH1F*> &ho)    // Output: ht,hs,hres,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histo & no background (2 pars) */
  {
    const int npars = 2;
    std::vector<double> setpars{1,1};

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();

    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_nbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");

    ht_cp->Fit(f1,"RWL"); 
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hres = (TH1F*)ht_cp->Clone(); hres->Add(ht_cp,hst,1,-1); 
    ho = {ht_cp,hst,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_2hs_1hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 TH1F* hbg,            // bg histo for fit
			 std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,hres,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histo & 1 bg histo (3 pars) */
  {
    const int npars = 3;
    std::vector<double> setpars{1,1,0};

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();
    TH1F *hbg_cp1 = (TH1F*)hbg->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1,hbg_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_1hbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    f1->SetParName(2,"B");

    ht_cp->Fit(f1,"RWL");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hbg_sc = (TH1F*)hbg_cp1->Clone(); hbg_sc->Scale(pars[0]*pars[2]); 
    TH1F *hsANDbg = (TH1F*)hs1_cp2->Clone(); hsANDbg->Add(hst,hbg_sc);
    TH1F *hres = (TH1F*)hs2_cp1->Clone(); hres->Add(ht_cp,hsANDbg,1,-1); 
    ho = {ht_cp,hst,hbg_sc,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  
  void ShiftHistogramX(TH1F* hist, double shiftAmount) {
    /* Shifts the histogram x by a constant amount
       WARNING: Doesn't work very well if the shift amount is
                much smaller than the bin width!
    */
    // Get the number of bins in the histogram
    int nBins = hist->GetNbinsX();
    
    // Create a temporary histogram to store the shifted content
    TH1F* tempHist = (TH1F*)hist->Clone("tempHist");
    
    // Loop over the bins and shift the bin contents and centers
    for (int i = 1; i <= nBins; ++i) {
      double binCenter = hist->GetBinCenter(i);
      double shiftedBinCenter = binCenter + shiftAmount;
        
      // Calculate the shifted bin number
      int shiftedBin = hist->GetXaxis()->FindBin(shiftedBinCenter);
        
      // Check if the shifted bin is within the histogram range
      if (shiftedBin >= 1 && shiftedBin <= nBins) {
  	tempHist->SetBinContent(shiftedBin, hist->GetBinContent(i));
  	tempHist->SetBinError(shiftedBin, hist->GetBinError(i));
      }
    }
    
    // Copy the shifted content back to the original histogram
    hist->Reset("ICESM");  // Reset the original histogram
    hist->Add(tempHist);  // Add the shifted content to the original histogram
    
    delete tempHist;  // Delete the temporary histogram
  }

  TF1* fit_2hs_1hbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  TH1F* hbg,            // bg histo for fit
				  std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				  std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,hres,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histo & 1 bg histo (3+2 pars) 
     This time x position of the signal histos are also free parameters
  */
  {
    const int npars = 5;
    std::vector<double> setpars{1,1,0,xOff_range[1],xOff_range[3]};

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();
    TH1F *hbg_cp1 = (TH1F*)hbg->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1,hbg_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_1hbg_xOffVary,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    f1->SetParName(2,"B");
    f1->SetParName(3,"pOff");
    f1->SetParName(4,"nOff");

    // setting limits for p and n peak offset varitations
    f1->SetParLimits(3,xOff_range[1],xOff_range[2]);
    f1->SetParLimits(4,xOff_range[3],xOff_range[4]);

    ht_cp->Fit(f1,"SQER");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    // ShiftHistogramX(hs1_cp1,-pars[3]);
    // ShiftHistogramX(hs2_cp1,-pars[4]);

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]); 
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hbg_sc = (TH1F*)hbg_cp1->Clone(); hbg_sc->Scale(pars[0]*pars[2]); 
    TH1F *hsANDbg = (TH1F*)hs1_cp2->Clone(); hsANDbg->Add(hst,hbg_sc);
    TH1F *hres = (TH1F*)hs2_cp1->Clone(); hres->Add(ht_cp,hsANDbg,1,-1); 
    ho = {ht_cp,hst,hbg_sc,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_2hs_2hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 TH1F* hbg1,           // 1st bg histo for fit
			 TH1F* hbg2,           // 2nd bg histo for fit
			 std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,hres,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histo & 2 bg histo (3 pars) */
  {
    const int npars = 3;
    std::vector<double> setpars{1,1,0};

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();
    TH1F *hbg1_cp1 = (TH1F*)hbg1->Clone();
    TH1F *hbg2_cp1 = (TH1F*)hbg2->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1,hbg1_cp1,hbg2_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_2hbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    f1->SetParName(2,"B");

    ht_cp->Fit(f1,"RWL");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hbg1_cp2 = (TH1F*)hbg1_cp1->Clone(); hbg1_cp2->Scale(pars[0]*pars[2]); 
    TH1F *hbg2_cp2 = (TH1F*)hbg2_cp1->Clone(); hbg2_cp2->Scale(pars[0]*pars[2]); 
    TH1F *hbg_sc = (TH1F*)hbg1_cp2->Clone(); hbg_sc->Add(hbg1_cp2,hbg2_cp2); 
    TH1F *hsANDbg = (TH1F*)hs1_cp2->Clone(); hsANDbg->Add(hst,hbg_sc);
    TH1F *hres = (TH1F*)hs2_cp1->Clone(); hres->Add(ht_cp,hsANDbg,1,-1); 
    ho = {ht_cp,hst,hbg_sc,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_2hs_2hbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  TH1F* hbg1,           // 1st bg histo for fit
				  TH1F* hbg2,           // 2nd bg histo for fit
				  std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				  std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,hres,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histo & 2 bg histo (3+2 pars)
     This time x position of the signal histos are also free parameters
  */
  {
    const int npars = 5;
    std::vector<double> setpars{1,1,0,xOff_range[1],xOff_range[3]};

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();
    TH1F *hbg1_cp1 = (TH1F*)hbg1->Clone();
    TH1F *hbg2_cp1 = (TH1F*)hbg2->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1,hbg1_cp1,hbg2_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_2hbg_xOffVary,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    f1->SetParName(2,"B");
    f1->SetParName(3,"pOff");
    f1->SetParName(4,"nOff");

    // setting limits for p and n peak offset varitations
    f1->SetParLimits(3,xOff_range[1],xOff_range[2]);
    f1->SetParLimits(4,xOff_range[3],xOff_range[4]);

    ht_cp->Fit(f1,"R");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hbg1_cp2 = (TH1F*)hbg1_cp1->Clone(); hbg1_cp2->Scale(pars[0]*pars[2]); 
    TH1F *hbg2_cp2 = (TH1F*)hbg2_cp1->Clone(); hbg2_cp2->Scale(pars[0]*pars[2]); 
    TH1F *hbg_sc = (TH1F*)hbg1_cp2->Clone(); hbg_sc->Add(hbg1_cp2,hbg2_cp2); 
    TH1F *hsANDbg = (TH1F*)hs1_cp2->Clone(); hsANDbg->Add(hst,hbg_sc);
    TH1F *hres = (TH1F*)hs2_cp1->Clone(); hres->Add(ht_cp,hsANDbg,1,-1); 
    ho = {ht_cp,hst,hbg_sc,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_2hs_1pbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 int Opoly,            // Order of poly to fit bg
			 std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histos & 1 poly bg (2+Opoly+1 pars) */
  {
    const int npars = 2+Opoly+1;
    std::vector<double> setpars{1,1}; for (int i=2;i<npars;i++) setpars.push_back(0);

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1,Opoly);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_1pbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    set_poly_par_names(f1,2,Opoly);

    ht_cp->Fit(f1,"R");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    // drawing the polynomial background function
    FitFn *ffn2 = new FitFn(Opoly);
    TF1* bgf = new TF1("bgf",ffn2,&FitFn::ffn_poly,fit_range[0],fit_range[1],Opoly+1);
    bgf->SetNpx(500);
    bgf->SetParameters(&fit::GetFitParams(f1)[2]);
    //bgf->SetLineColor(46);

    // now get a bg histo from bgf
    TH1F *hbg_sc = (TH1F*)hs1_cp1->Clone(); util_pd::TF1toTH1F(bgf,hbg_sc);
    hbg_sc->SetMarkerColor(46); hbg_sc->SetLineColor(46); //hbg_sc->SetLineWidth(2);

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hsANDbg = (TH1F*)hs1_cp2->Clone(); hsANDbg->Add(hst,hbg_sc);
    TH1F *hres = (TH1F*)hs2_cp1->Clone(); hres->Add(ht_cp,hsANDbg,1,-1); 
    ho = {ht_cp,hst,hbg_sc,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_2hs_1pbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  int Opoly,            // Order of poly to fit bg
				  std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				  std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histos & 1 poly bg (2+2+Opoly+1 pars)
     This time x position of the signal histos are also free parameters
  */
  {
    const int npars = 4+Opoly+1;
    std::vector<double> setpars{1,1,xOff_range[1],xOff_range[3]}; 
    for (int i=4;i<npars;i++) setpars.push_back(0);

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1,Opoly);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_1pbg_xOffVary,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    f1->SetParName(2,"pOff");
    f1->SetParName(3,"nOff");
    set_poly_par_names(f1,4,Opoly);

    // setting limits for p and n peak offset varitations
    f1->SetParLimits(2,xOff_range[1],xOff_range[2]);
    f1->SetParLimits(3,xOff_range[3],xOff_range[4]);

    ht_cp->Fit(f1,"R");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    // drawing the polynomial background function
    FitFn *ffn2 = new FitFn(Opoly);
    TF1* bgf = new TF1("bgf",ffn2,&FitFn::ffn_poly,fit_range[0],fit_range[1],Opoly+1);
    bgf->SetNpx(500);
    bgf->SetParameters(&fit::GetFitParams(f1)[4]);
    //bgf->SetLineColor(46);

    // now get a bg histo from bgf
    TH1F *hbg_sc = (TH1F*)hs1_cp1->Clone(); util_pd::TF1toTH1F(bgf,hbg_sc);
    hbg_sc->SetMarkerColor(46); hbg_sc->SetLineColor(46); //hbg_sc->SetLineWidth(2);

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hsANDbg = (TH1F*)hs1_cp2->Clone(); hsANDbg->Add(hst,hbg_sc);
    TH1F *hres = (TH1F*)hs2_cp1->Clone(); hres->Add(ht_cp,hsANDbg,1,-1); 
    ho = {ht_cp,hst,hbg_sc,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_2hs_1gbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 std::vector<double> const & gfit_params, // Gaussian fit params
			 std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histos & 1 poly bg (2+Opoly+1 pars) */
  {
    const int npars = 5;
    std::vector<double> setpars{1,1}; for (int i=2;i<npars;i++) setpars.push_back(gfit_params[i-2]);
    //std::vector<double> setpars{1,1,10,-0.36,0.4}; 

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_1gbg,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    set_gaus_par_names(f1,2);

    ht_cp->Fit(f1,"R");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    // drawing the polynomial background function
    FitFn *ffn2 = new FitFn();
    TF1* bgf = new TF1("bgf",ffn2,&FitFn::ffn_gaus,fit_range[0],fit_range[1],3);
    bgf->SetNpx(500);
    bgf->SetParameters(&fit::GetFitParams(f1)[2]);
    //bgf->SetLineColor(46);

    // now get a bg histo from bgf
    TH1F *hbg_sc = (TH1F*)hs1_cp1->Clone(); util_pd::TF1toTH1F(bgf,hbg_sc);
    hbg_sc->SetMarkerColor(46); hbg_sc->SetLineColor(46); //hbg_sc->SetLineWidth(2);

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hsANDbg = (TH1F*)hs1_cp2->Clone(); hsANDbg->Add(hst,hbg_sc);
    TH1F *hres = (TH1F*)hs2_cp1->Clone(); hres->Add(ht_cp,hsANDbg,1,-1); 
    ho = {ht_cp,hst,hbg_sc,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_2hs_1gbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				  std::vector<double> const & gfit_params, // Gaussian fit params
				  std::vector<TH1F*> &ho)    // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  /* TH Interpolation fit using 2 signal histos & 1 poly bg (2+2+Opoly+1 pars)
     This time x position of the signal histos are also free parameters
  */
  {
    const int npars = 4+3;
    std::vector<double> setpars{1,1,xOff_range[1],xOff_range[3]};
    for (int i=4;i<npars;i++) setpars.push_back(gfit_params[i-2]);

    TH1F *ht_cp = (TH1F*)ht->Clone(); 
    TH1F *hs1_cp1 = (TH1F*)hs1->Clone(); 
    TH1F *hs2_cp1 = (TH1F*)hs2->Clone();
 
    FitFn *ffn = new FitFn(hs1_cp1,hs2_cp1);
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_2hs_1gbg_xOffVary,fit_range[0],fit_range[1],npars);
    f1->SetNpx(2000);
    f1->SetParameters(&setpars[0]);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"R");
    f1->SetParName(2,"pOff");
    f1->SetParName(3,"nOff");
    set_gaus_par_names(f1,4);

    // setting limits for p and n peak offset varitations
    f1->SetParLimits(2,xOff_range[1],xOff_range[2]);
    f1->SetParLimits(3,xOff_range[3],xOff_range[4]);

    ht_cp->Fit(f1,"R");
    std::vector<double> pars = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars.push_back(f1->GetParameter(i));}

    // drawing the polynomial background function
    FitFn *ffn2 = new FitFn();
    TF1* bgf = new TF1("bgf",ffn2,&FitFn::ffn_gaus,fit_range[0],fit_range[1],3);
    bgf->SetNpx(500);
    bgf->SetParameters(&fit::GetFitParams(f1)[4]);
    //bgf->SetLineColor(46);

    // now get a bg histo from bgf
    TH1F *hbg_sc = (TH1F*)hs1_cp1->Clone(); util_pd::TF1toTH1F(bgf,hbg_sc);
    hbg_sc->SetMarkerColor(46); hbg_sc->SetLineColor(46); //hbg_sc->SetLineWidth(2);

    TH1F *hs1_cp2 = (TH1F*)hs1_cp1->Clone(); hs1_cp2->Scale(pars[0]);
    TH1F *hs2_cp2 = (TH1F*)hs2_cp1->Clone(); hs2_cp2->Scale(pars[0]*pars[1]);
    TH1F *hst = (TH1F*)hs1_cp2->Clone(); hst->Add(hs1_cp2,hs2_cp2);
    TH1F *hsANDbg = (TH1F*)hs1_cp2->Clone(); hsANDbg->Add(hst,hbg_sc);
    TH1F *hres = (TH1F*)hs2_cp1->Clone(); hres->Add(ht_cp,hsANDbg,1,-1); 
    ho = {ht_cp,hst,hbg_sc,hres,hs1_cp2,hs2_cp2};
    
    return f1;
  }

  TF1* fit_1gs_nbg (std::vector<double> const & fit_range, // [0]=>xmin,[1]=>xmax (for 1st fit)
                                                           // [2]=>nSLow,[3]=>nSHi (for 2nd fit)
  		    TH1F const * ht)                       // input histogram
  /* Fitting signal peak using a Gaussian (3 pars). Fits twice for optimization. */
  {
    int const npars = 3;
    
    // let's not edit the original histogram
    TH1F *ht_cp = (TH1F*)ht->Clone();

    // define fit function
    FitFn *ffn = new FitFn();
    TF1 *f1 = new TF1("f1",ffn,&FitFn::ffn_gaus,fit_range[0],fit_range[1],npars);
    f1->SetNpx(1000);
    f1->SetParName(0,"Norm");
    f1->SetParName(1,"Mean"); 
    f1->SetParName(2,"Sigma"); 
    f1->SetRange(fit_range[0],fit_range[1]);
    
    // first fit
    ht_cp->GetXaxis()->SetRangeUser(fit_range[0],fit_range[1]);
    double norm = ht_cp->GetMaximum();
    double mean = ht_cp->GetMean();
    double sigma = ht_cp->GetStdDev();
    
    f1->SetRange(fit_range[0],fit_range[1]);
    f1->SetParameter(0,norm);
    f1->SetParameter(1,mean);
    f1->SetParameter(2,sigma);

    ht_cp->Fit(f1,"R");
    std::vector<double> pars1 = GetFitParams(f1);
    //for (int i=0;i<npars;i++) {pars1.push_back(f1->GetParameter(i));}

    // Second fit with tailored range
    double llim = pars1[1] - fit_range[2]*pars1[2];
    double hlim = pars1[1] + fit_range[3]*pars1[2];
    f1->SetParameters(&pars1[0]);
    f1->SetRange(llim,hlim);
    ht_cp->Fit(f1,"R");

    // re-adjust histogram range
    ht_cp->GetXaxis()->SetRangeUser(ht->GetXaxis()->GetXmin(),
				    ht->GetXaxis()->GetXmax());

    return f1;
  }


} // namespace fit
