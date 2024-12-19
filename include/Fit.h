#ifndef FIT_H
#define FIT_H

#include <vector>

#include "TF1.h"

#include "FitFns.h"
#include "Utilities.h"

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

  // Returns a vector filled with pair of fit pramater values and errors from f1
  std::vector<std::pair<double,double>> GetFitParamANDError(TF1 * const f1);

  // Side band (SB) fit using 1 polynomial and 2 reject points (Opoly+1 pars
  // NOTE: This fit method just fits the bg
  TF1* fit_1pbg_SB (std::vector<double> const & fit_range,
		    std::vector<double> const & reject_points,
		    int Opoly,          // Order of poly to fit bg
		    std::vector<double> const & initial_guesses,
		    TH1F* ht,          // total histo to fit
		    std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs 

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

  // TH Interpolation fit using 1 signal histo & 1 poly bg (1+Opoly+1 pars)
  TF1* fit_1hs_1pbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs,             // 1st signal histo for fit
				  int Opoly,            // Order of poly to fit bg
				  std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				  std::vector<TH1F*> &ho);   // Output: ht,hs,hbg

  // TH Interpolation fit using 1 signal histo & 1 bg histo (2 pars)
  TF1* fit_1hs_1hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs,             // signal histo for fit
			 TH1F* hbg,            // bg histo for fit
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs (N=par[0],B=par[1])

  // TH Interpolation fit using 1 signal histo (vary x offsets for both) & 1 bg histo (2+2 pars)
  // This time x position of the signal histo is also free parameter
  TF1* fit_1hs_1hbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs,             // signal histo for fit
				  TH1F* hbg,            // bg histo for fit
				  std::vector<double> const & xOff_range,  // ranges to vary hs & hbg offsets
				  std::vector<TH1F*> &ho);    // Output: ht,hs,hres,N*hs1,N*R*hs2 (N=par[0],R=par[1])
  
  // TH Interpolation fit using 2 signal histos & no background (2 pars)
  TF1* fit_2hs_nbg_THI (std::vector<double> const & fit_range,
			TH1F* ht,             // total histo to fit 
			TH1F* hs1,            // 1st signal histo for fit
			TH1F* hs2,            // 2nd signal histo for fit
			std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos (vary x offsets for both) & no background (2+2 pars)
  // This time x position of the signal histo is also free parameter
  TF1* fit_2hs_nbg_THI_xOffVary (std::vector<double> const & fit_range,
				 TH1F* ht,             // total histo to fit 
				 TH1F* hs1,            // 1st signal histo for fit
				 TH1F* hs2,            // 2nd signal histo for fit
				 std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 1 bg histo (3 pars)
  TF1* fit_2hs_1hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 TH1F* hbg,            // bg histo for fit
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])s

  // TH Interpolation fit using 2 signal histos & 1 bg histo & vary hs1 and hs2 x positions (5 pars)
  TF1* fit_2hs_1hbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  TH1F* hbg,            // bg histo for fit
				  std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				  std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 2 bg histos (3 pars)
  TF1* fit_2hs_2hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 TH1F* hbg1,           // 1st bg histo for fit
			 TH1F* hbg2,           // 2nd bg histo for fit
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 2 bg histos (5 pars)
  // x offsets are varied for both signal and background histos
  TF1* fit_2hs_2hbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  TH1F* hbg1,           // 1st bg histo for fit
				  TH1F* hbg2,           // 2nd bg histo for fit
				  std::vector<double> const & xOff_range,  // ranges to vary hs1/hbg1 & hs2/hbg2 offsets
				  std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 3 bg histos (4 pars)  
  TF1* fit_2hs_3hbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 TH1F* hbg1,           // 1st bg histo for fit
			 TH1F* hbg2,           // 2nd bg histo for fit
			 TH1F* hbg3,           // 3rd bg histo for fit (for param B2)
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 3 bg histos (5 pars)
  // x offsets are varied for both signal histos and the first and second hbg
  // B2 belongs to the third hbg  
  TF1* fit_2hs_3hbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  TH1F* hbg1,           // 1st bg histo for fit
				  TH1F* hbg2,           // 2nd bg histo for fit
				  TH1F* hbg3,           // 3rd bg histo for fit (for param B2)				  
				  std::vector<double> const & xOff_range,  // ranges to vary hs1/hbg1 & hs2/hbg2 offsets
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
			 std::vector<double> &setpars,
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 1 poly bg (2+2+Opoly+1 pars)
  // This time x positions of the signal histos are also free parameters
  TF1* fit_2hs_1pbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  int Opoly,            // Order of poly to fit bg
				  std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				  std::vector<double> &setpars,
				  std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 1 Gaussian bg (2+3 pars)
  TF1* fit_2hs_1gbg_THI (std::vector<double> const & fit_range,
			 TH1F* ht,             // total histo to fit 
			 TH1F* hs1,            // 1st signal histo for fit
			 TH1F* hs2,            // 2nd signal histo for fit
			 std::vector<double> const & gfit_params, // Gaussian fit params
			 std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // TH Interpolation fit using 2 signal histos & 1 poly bg (2+2+3 pars)
  // This time x positions of the signal histos are also free parameters
  TF1* fit_2hs_1gbg_THI_xOffVary (std::vector<double> const & fit_range,
				  TH1F* ht,             // total histo to fit 
				  TH1F* hs1,            // 1st signal histo for fit
				  TH1F* hs2,            // 2nd signal histo for fit
				  std::vector<double> const & xOff_range,  // ranges to vary hs1 & hs2 offsets
				  std::vector<double> const & gfit_params, // Gaussian fit params
				  std::vector<TH1F*> &ho);   // Output: ht,hs,hbg,N*hs1,N*R*hs2 (N=par[0],R=par[1])

  // Fitting signal peak using a Gaussian (3 pars)
  TF1* fit_1gs_nbg (std::vector<double> const & fit_range, // [0]=>xmin,[1]=>xmax (for 1st fit)
                                                           // [2]=>nSLow,[3]=>nSHi (for 2nd fit)
		    TH1F const * ht);                      // input histogram
}


#endif
