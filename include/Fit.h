#ifndef FIT_H
#define FIT_H

#include <vector>

#include "TF1.h"

#include "FitFns.h"

/* IMPORTANT: Be sure to make a canvas available before calling a fit function. */

namespace fit {

  // Sets the name of polynomial parameters
  void set_poly_par_names (TF1* f1,
			   int const Spar,  // starting parameter
			   int const Opoly);// Order of polynomial fit

  // Returns a vector filled with fit pramater values from f1
  std::vector<double> GetFitParams(TF1 * const f1);

  // Returns a vector filled with fit pramater errors from f1
  std::vector<double> GetFitParamErrors(TF1 * const f1);

  // Side band (SB) fit using 1 polynomial and 2 reject points (Opoly+1 pars)
  TF1* fit_1pbg_SB (std::vector<double> const & fit_range,
		    std::vector<double> const & reject_points,
		    int Opoly,          // Order of poly to fit bg
		    std::vector<double> const & initial_guesses,
		    TH1F* ht);          // total histo to fit

  // TH Interpolation fit using 1 signal histo & no background (1 par)
  TF1* fit_1hs_nbg_THI (std::vector<double> const & fit_range,
  			TH1F* ht,           // total histo to fit
  			TH1F* hs,           // signal histo for fit
  			std::vector<TH1F*> &ho); // Output: ht,hs (fitted w/ proper scaling)

  // TH Interpolation fit using 1 signal histo & 1 poly bg (1+Opoly+1 pars)
  TF1* fit_1hs_1pbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs,             // 1st signal histo for fit
			 int Opoly,            // Order of poly to fit bg
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg

  // TH Interpolation fit using 1 signal histo & 1 bg histo (2 pars)
  TF1* fit_1hs_1hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs,             // signal histo for fit
			 TH1F* hbg,            // bg histo for fit
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs (N=par[0],B=par[1])

  // TH Interpolation fit using 2 signal histos & no background (2 pars)
  TF1* fit_2hs_nbg_THI (std::vector<double> const & fit_range,
			TH1F* ht,             // total histo to fit 
			TH1F* hs1,            // 1st signal histo for fit
			TH1F* hs2,            // 2nd signal histo for fit
			std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 1 bg histo (3 pars)
  TF1* fit_2hs_1hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 TH1F* hbg,            // bg histo for fit
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 1 signal histo & 1 poly bg (1+Opoly+1 pars)
  TF1* fit_2hs_nbg_THI (std::vector<double> const & fit_range,
			TH1F* ht,             // total histo to fit 
			TH1F* hs,             // signal histo for fit
			int Opoly,            // Order of poly to fit bg
			std::vector<TH1F*> &ho);   // Output: ht,hs,hbg

  // TH Interpolation fit using 2 signal histos & 1 poly bg (2+Opoly+1 pars)
  TF1* fit_2hs_1pbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 int Opoly,            // Order of poly to fit bg
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // Fitting signal peak using a Gaussian (3 pars)
  TF1* fit_1gs_nbg (std::vector<double> const & fit_range, // [0]=>xmin,[1]=>xmax (for 1st fit)
                                                           // [2]=>nSLow,[3]=>nSHi (for 2nd fit)
		    TH1F const * ht);                      // input histogram
}


#endif
