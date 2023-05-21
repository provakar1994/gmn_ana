#ifndef FIT_H
#define FIT_H

#include <vector>

#include "TF1.h"

#include "FitFns.h"

/* IMPORTANT: Be sure to make a canvas available before calling a fit function. */

namespace fit {

  // TH Interpolation fit using 1 signal histo & no background (1 par)
  TF1* fit_1hs_nbg_THI (std::vector<double> fit_range,
  			TH1D* ht,           // total histo to fit
  			TH1D* hs,           // signal histo for fit
  			std::vector<TH1D*> &ho); // Output: ht,hs (fitted w/ proper scaling)

  // TH Interpolation fit using 2 signal histos & no background (2 pars)
  TF1* fit_2hs_nbg_THI (std::vector<double> fit_range,
			TH1D* ht,             // total histo to fit 
			TH1D* hs1,            // 1st signal histo for fit
			TH1D* hs2,            // 2nd signal histo for fit
			std::vector<TH1D*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 1 bg histo (3 pars)
  TF1* fit_2hs_1hbg_THI (std::vector<double> fit_range,
			 TH1D* ht,             // total histo to fit 
			 TH1D* hs1,            // 1st signal histo for fit
			 TH1D* hs2,            // 2nd signal histo for fit
			 TH1D* hbg,            // bg histo for fit
			 std::vector<TH1D*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 1 signal histo & 1 poly bg (1+Opoly+1 pars)
  TF1* fit_2hs_nbg_THI (std::vector<double> fit_range,
			TH1D* ht,             // total histo to fit 
			TH1D* hs,             // signal histo for fit
			int Opoly,            // Order of poly to fit bg
			std::vector<TH1D*> &ho);   // Output: ht,hs,hbg

  // TH Interpolation fit using 2 signal histos & 1 poly bg (2+Opoly+1 pars)
  TF1* fit_2hs_1pbg_THI (std::vector<double> fit_range,
			 TH1D* ht,             // total histo to fit 
			 TH1D* hs1,            // 1st signal histo for fit
			 TH1D* hs2,            // 2nd signal histo for fit
			 int Opoly,            // Order of poly to fit bg
			 std::vector<TH1D*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])
}


#endif
