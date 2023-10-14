/*
  Class to calculate nucleon EMFFs using the parametrization of choice.
  ----
  NOTE:
  1. Fit functions return FF values as, GEp/GD, GMp/(GD*mu_p), GEn/GD, 
     GMn/(GD*mu_n). But GetFF() and GetFFwErr() fns return pure FF values.
     StripGDandMu() fn helps us achieve that.
  2. In g4sbs and SIMC we have used Seamus's fit for GEn parametrization
     and Kelly's fit for all the remaining EMFF parametrizations. 
  ----
  To-do:
  1. Implement StripGDandMu for error values
*/
#ifndef EMFF_FITS_H
#define EMFF_FITS_H

#include <iostream>
#include <vector>
#include <cmath>

#include "Constants.h"
#include "KinematicVar.h"

enum class G_t {kGEp = 1, kGMp = 2, kGEn = 3, kGMn = 4};

class EMFFFits {
 public:
  EMFFFits(){};
  // Fns to get pure EMFF values
  virtual double GetFF(G_t kG, double Q2) = 0;
  virtual std::vector<double> GetFFwErr(G_t kG, double const Q2) = 0;
  // Utility fns
  static double GetGDip(double const Q2) {return pow(1./(1. + Q2/0.71), 2);}
  static void StripGDandMu(G_t const kG, double const Q2, double & value);
  static G_t StrToG_t(const std::string& str);
};

//##############################
class Kelly2004 : public EMFFFits {
  // J. J. Kelly: PHYSICAL REVIEW C 70, 068202 (2004)
 public:
  Kelly2004(){};
  double GetFF(G_t const kG, double const Q2) override {
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);
    int err = KellyFit(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    return GNGD_Fit[0];
  }
  std::vector<double> GetFFwErr(G_t kG, double const Q2) override {
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error */
    std::cout << "WARNING: Error calculation is yet to be added for Kelly Fit (2004)\n";
    std::vector<double> result{Kelly2004::GetFF(kG,Q2),-1000};
    return result;
  }
 private:
  int KellyFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
};

//##############################
class Seamus20XX : public EMFFFits {
  // Seamus Riordan: Couldn't find the corresponding paper
 public:
  Seamus20XX(){};
  double GetFF(G_t const kG, double const Q2) override {
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);
    int err = SeamusFit(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    return GNGD_Fit[0];
  }
  std::vector<double> GetFFwErr(G_t kG, double const Q2) override {
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error */
    std::cout << "WARNING: Error calculation is yet to be added for Seamus Fit (2004)\n";
    std::vector<double> result{Seamus20XX::GetFF(kG,Q2),-1000};
    return result;
  }
 private:
  int SeamusFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
};

//##############################
class Ye2017 : public EMFFFits {
  // Z. Ye et al: Physics Letters B 777 (2018) 8–15
 public:
  Ye2017(){};
  double GetFF(G_t const kG, double const Q2) override {
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);		
    int err = YeFit(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    return GNGD_Fit[0];
  }
  std::vector<double> GetFFwErr(G_t kG, double const Q2) override {
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error */
    std::vector<double> result;
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);
    int err = YeFit(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    result = {GNGD_Fit[0],GNGD_Err[0]};
    return result;
  }
 private:
  int YeFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
};


#endif
