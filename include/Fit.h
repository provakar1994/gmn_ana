#ifndef FIT_H
#define FIT_H

#include <vector>

#include "TF1.h"

#include "FitFns.h"

/* IMPORTANT: Be sure to make a canvas available before calling a fit function. */

namespace fit {

  /* // TH Interpolation fit using 1 signal histo & no background */
  /* TF1* fit_1hs_nbg_THI (vector<double> fit_range, */
  /* 			TH1D* ht,           // total histo to fit  */
  /* 			TH1D* hs,           // signal histo for fit */
  /* 			vector<TH1D*> &ho); // Output: ht,hs (fitted w/ proper scaling) */

  // TH Interpolation fit using 1 signal histo & no background
  TF1* fit_2hs_nbg_THI (vector<double> fit_range,
			TH1D* ht,             // total histo to fit 
			TH1D* hs1,            // 1st signal histo for fit
			TH1D* hs2,            // 2nd signal histo for fit
			vector<TH1D*> &ho);   // Output: ht,hs (fitted w/ proper scaling)

  // TH Interpolation fit using 1 signal histo & 1 bg histo
  TF1* fit_2hs_1hbg_THI (vector<double> fit_range,
			 TH1D* ht,             // total histo to fit 
			 TH1D* hs1,            // 1st signal histo for fit
			 TH1D* hs2,            // 2nd signal histo for fit
			 TH1D* hbg,            // bg histo for fit
			 vector<TH1D*> &ho);   // Output: ht,hs (fitted w/ proper scaling)

}


#endif
