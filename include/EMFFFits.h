/*
  Class to calculate nucleon EMFFs using the parametrization of choice.
  ----
  NOTE:
  1. Fit functions return FF values as, GEp/GD, GMp/(GD*mu_p), GEn/GD, 
     GMn/(GD*mu_n). But GetFF() and GetFFwErr() fns return pure FF values.
     StripGDandMu() fn helps us achieve that.
  2. In g4sbs and SIMC we have used Seamus's fit for GEn parametrization
     and Kelly's fit for all the remaining EMFF parametrizations.
  3. All fits were checked against their original papers (and Ye2017 against the
     official code and lookup tables in DB/FF_fits/Ye2018/) on 2026-10-28.
*/
#ifndef EMFF_FITS_H
#define EMFF_FITS_H

#include <iostream>
#include <vector>
#include <cmath>

#include "Constants.h"
#include "KinematicVar.h"
#include "LookUpTableReader.h"

enum class G_t {kGEp = 1, kGMp = 2, kGEn = 3, kGMn = 4};

class EMFFFits {
 public:
  EMFFFits(){};
  // Fns to get pure EMFF values
  virtual double GetFF(G_t kG, double Q2) = 0;
  virtual std::vector<double> GetFFwErr(G_t kG, double const Q2) = 0;
  // Fns to get effective FF values (ie w/o TPE corrections)
  virtual double GetFFNoTPE(G_t kG, double Q2) = 0;
  virtual std::vector<double> GetFFwErrNoTPE(G_t kG, double const Q2) = 0;
  // Utility fns
  static double GetGDip(double const Q2) {return pow(1./(1. + Q2/0.71), 2);}
  static void StripGDandMu(G_t const kG, double const Q2, double & value);
  static G_t StrToG_t(const std::string& str);
  // Full path of a file under the DB directory (e.g. "FF_fits/Ye2018/x.csv"). The DB
  // directory is set at build time and can be overridden w/ env. variable GMNANA_DB_DIR.
  static std::string GetDBPath(std::string const& relPath);
};

//##############################
class Galster1971 : public EMFFFits {
  // S. Galster: Nuclear Physics B32 (1971) 221-237
 public:
  Galster1971(){};
  double GetFF(G_t const kG, double const Q2) override {
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);
    int err = GalsterFit(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    return GNGD_Fit[0];
  }
  std::vector<double> GetFFwErr(G_t kG, double const Q2) override {
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error */
    std::cout << "WARNING: Error calculation is yet to be added for Galster Fit (1971)\n";
    std::vector<double> result{Galster1971::GetFF(kG,Q2),-1000};
    return result;
  }
  double GetFFNoTPE(G_t kG, double Q2) override {
    std::cout << "WARNING! No implementation provided!\n";
    return -1000;
  }
  std::vector<double> GetFFwErrNoTPE(G_t kG, double const Q2) override {
    std::cout << "WARNING! No implementation provided!\n";
    std::vector<double> result{-1000,-1000};
    return result;
  }
 private:
  int GalsterFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
};

//##############################
class Kelly2004 : public EMFFFits {
  // J. J. Kelly: PHYSICAL REVIEW C 70, 068202 (2004)
  // NOTE: TPE "mixed": GEp from polarization data (~Born), GMp from cross sections w/o TPE
  // corrections (effective). See scripts/gmn_extraction/README.md (Caveats on the fits).
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
  double GetFFNoTPE(G_t kG, double Q2) override {
    std::cout << "WARNING! No implementation provided!\n";
    return -1000;
  }
  std::vector<double> GetFFwErrNoTPE(G_t kG, double const Q2) override {
    std::cout << "WARNING! No implementation provided!\n";
    std::vector<double> result{-1000,-1000};
    return result;
  }
 private:
  int KellyFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
};

//##############################
class Seamus20XX : public EMFFFits {
  // GEn fit used in g4sbs & SIMC (presumably by Seamus Riordan). Unpublished, and not the
  // fit of S. Riordan et al., PRL 105, 262302 (2010). Use it only to reproduce the MC.
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
    std::cout << "WARNING: Error calculation is yet to be added for Seamus Fit\n";
    std::vector<double> result{Seamus20XX::GetFF(kG,Q2),-1000};
    return result;
  }
  double GetFFNoTPE(G_t kG, double Q2) override {
    std::cout << "WARNING! No implementation provided!\n";
    return -1000;
  }
  std::vector<double> GetFFwErrNoTPE(G_t kG, double const Q2) override {
    std::cout << "WARNING! No implementation provided!\n";
    std::vector<double> result{-1000,-1000};
    return result;
  }
 private:
  int SeamusFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
};

//##############################
class SIMCbigbite : public EMFFFits {
  // EMFFs exactly as implemented in the SIMC generator used for the GMn MC:
  // https://github.com/MarkKJones/simc_gfortran/blob/bigbite/physics_proton.f
  // (fofa_best_fit_p: Kelly GEp & GMp; fofa_best_fit_n: Seamus GEn & Kelly GMn).
  // It differs from Kelly2004 only by SIMC's rounded magnetic moments, GMp(0) = 2.79 and
  // GMn(0) = -1.913. Use it to compute the MC Born CS ratio (kine::sigmaBorn_ratio_MC),
  // so that data/MC corrections use the same FFs as the generator.
 public:
  SIMCbigbite(){};
  static constexpr double kGMp0 = 2.79;   // GMp(Q2=0) used in SIMC
  static constexpr double kGMn0 = -1.913; // GMn(Q2=0) used in SIMC
  double GetFF(G_t const kG, double const Q2) override;
  std::vector<double> GetFFwErr(G_t kG, double const Q2) override {
    std::cout << "WARNING: No errors available for the SIMC EMFFs\n";
    std::vector<double> result{SIMCbigbite::GetFF(kG,Q2),-1000};
    return result;
  }
  double GetFFNoTPE(G_t kG, double Q2) override {
    std::cout << "WARNING! No implementation provided!\n";
    return -1000;
  }
  std::vector<double> GetFFwErrNoTPE(G_t kG, double const Q2) override {
    std::cout << "WARNING! No implementation provided!\n";
    std::vector<double> result{-1000,-1000};
    return result;
  }
 private:
  Kelly2004 kellyfit;
  Seamus20XX seamusfit;
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
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error
       The value is from the z-expansion parametrization (agrees w/ Ye's lookup tables to
       <~1e-3 of the error). The error is the full fit error dG from Ye's lookup tables
       (log-interpolated in Q2), not the parametrized one (see GetFFwParErr). */
    double FF_err;
    int err = YeTableError(kG, Q2, FF_err);
    std::vector<double> result = {Ye2017::GetFF(kG,Q2), FF_err};
    return result;
  }
  std::vector<double> GetFFwParErr(G_t kG, double const Q2) {
    /* Same as GetFFwErr, but w/ Ye's parametrization of the error (dG_Par in the lookup
       tables), a smooth fit to dG. It differs from dG by up to ~15% for 0.01 < Q2 < 20 GeV2. */
    std::vector<double> result;
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);
    int err = YeFit(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    // Ye's error parametrization returns dG/(G(0)*GD) w/ G(0) = 1 (GEp, GEn), mu_p (GMp), mu_n (GMn),
    // so StripGDandMu() converts it to dG. Verified against dG_Par in Ye's lookup tables.
    StripGDandMu(kG,Q2,GNGD_Err[0]);
    result = {GNGD_Fit[0],fabs(GNGD_Err[0])};
    return result;
  }
  // No TPE. !! PLACEHOLDER: currently the "noExtraTPE" fit, which is still
  // SIFF-Blunden TPE corrected (see DB/FF_fits/Ye2018/README.md)
  double GetFFNoTPE(G_t const kG, double const Q2) override {
    double FF, FF_err;
    const int kID = static_cast<int>(kG);		
    int err = YeFitNoTPE(kID, Q2, FF, FF_err);
    return FF;
  }
  std::vector<double> GetFFwErrNoTPE(G_t kG, double const Q2) override {
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error */
    std::vector<double> result;
    double FF, FF_err;
    const int kID = static_cast<int>(kG);
    int err = YeFitNoTPE(kID, Q2, FF, FF_err);
    result = {FF,FF_err};
    return result;
  }
 private:
  int YeFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
  int YeTableError(G_t const kG, double const kQ2, double &FF_err);
  int YeFitNoTPE(const int kID, const double kQ2, double &FF, double &FF_err);
};

//##############################
class Christy2022 : public EMFFFits {
  // M. E. Christy et al: PHYSICAL REVIEW LETTERS 128, 102002 (2022)
  // The global fit (Eq. 3) is to cross sections only, w/o hard TPE corrections, i.e. it
  // gives effective (TPE-uncorrected) proton FFs. Hence only GetFFNoTPE() is provided.
  // NOTE: it uses an updated radiative correction prescription, different from SIMC's.
 public:
  Christy2022(){};
  double GetFF(G_t const kG, double const Q2) override {
    std::cout << "WARNING: Christy fit (2022) is TPE uncorrected, use GetFFNoTPE()\n";
    return -1000;
  }
  std::vector<double> GetFFwErr(G_t kG, double const Q2) override {
    std::cout << "WARNING: Christy fit (2022) is TPE uncorrected, use GetFFwErrNoTPE()\n";
    std::vector<double> result{-1000,-1000};
    return result;
  }
  double GetFFNoTPE(G_t const kG, double const Q2) override {
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);
    int err = ChristyFit(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    return GNGD_Fit[0];
  }
  std::vector<double> GetFFwErrNoTPE(G_t kG, double const Q2) override {
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error */
    std::cout << "WARNING: Error calculation is yet to be added for Christy Fit (2022)\n";
    std::vector<double> result{Christy2022::GetFFNoTPE(kG,Q2),-1000};
    return result;
  }
 private:
  int ChristyFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
};

//##############################
class Arrington2007 : public EMFFFits {
  // J. Arrington, W. Melnitchouk, J. A. Tjon: PHYSICAL REVIEW C 76, 035205 (2007)
 public:
  Arrington2007(){};
  double GetFF(G_t const kG, double const Q2) override {
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);
    int err = ArringtonFit(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    return GNGD_Fit[0];
  }
  std::vector<double> GetFFwErr(G_t kG, double const Q2) override {
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error */
    std::cout << "WARNING: Error calculation is yet to be added for Arrington Fit (2007)\n";
    std::vector<double> result{Arrington2007::GetFF(kG,Q2),-1000};
    return result;
  }
  double GetFFNoTPE(G_t const kG, double const Q2) override {
    double GNGD_Fit[1], GNGD_Err[1];
    const int kID = static_cast<int>(kG);
    int err = ArringtonFitNoTPE(kID, Q2, GNGD_Fit, GNGD_Err);
    StripGDandMu(kG,Q2,GNGD_Fit[0]);
    return GNGD_Fit[0];
  }
  std::vector<double> GetFFwErrNoTPE(G_t kG, double const Q2) override {
    /* Output: GetFFnErr[0] = Fit value, GetFFnErr[1] = Fit error */
    std::cout << "WARNING: Error calculation is yet to be added for Arrington Fit (2007)\n";
    std::vector<double> result{Arrington2007::GetFFNoTPE(kG,Q2),-1000};
    return result;
  }
 private:
  int ArringtonFit(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
  int ArringtonFitNoTPE(const int kID, const double kQ2, double *GNGD_Fit, double *GNGD_Err);
};

#endif
