#include <iostream>

#include "EMFFFits.h"

int extractGMn(double const etheta, double const Q2, double const R_fit, double const R_fit_err, double const R_fit_err_sys) {

  // EMFF fits
  Ye2017 yefit;
  Kelly2004 kellyfit;
  Seamus20XX seamusfit;
  Christy2022 christyfit;
  Arrington2007 arfit;
  Galster1971 galfit;

  // EMFF extraction using same parametrization used in SIMC generator
  double GEp_kelly = kellyfit.GetFF(G_t::kGEp,Q2);
  double GMp_kelly = kellyfit.GetFF(G_t::kGMp,Q2);
  double GEn_seamus = seamusfit.GetFF(G_t::kGEn,Q2);
  double GMn_kelly = kellyfit.GetFF(G_t::kGMn,Q2);

  // EMFF extraction using Ye fit (2017), to get GMn from data/MC fit
  // -- p FFs
  // Effective
  double GEp_christy = christyfit.GetFF(G_t::kGEp,Q2);
  double GMp_christy = christyfit.GetFF(G_t::kGMp,Q2);
  // True
  double GEp_ye = yefit.GetFF(G_t::kGEp,Q2);
  double GMp_ye = yefit.GetFF(G_t::kGMp,Q2);
  double GEp_ar = arfit.GetFF(G_t::kGEp,Q2);
  double GMp_ar = arfit.GetFF(G_t::kGMp,Q2);
  // -- n FFs
  double GEn_ye = yefit.GetFF(G_t::kGEn,Q2);
  double GEn_kelly = kellyfit.GetFF(G_t::kGEn,Q2);
  double GEn_gal = galfit.GetFF(G_t::kGEn,Q2);

  // Getting sigma Born ratio (n/p) for the given Q2 and etheta 
  // (basically, mimicing the calculation done in MC generator)
  double sigmaBorn_Ratio_MC = kine::sigmaBorn_ratio(etheta,Q2,GEp_kelly,GMp_kelly,GEn_seamus,GMn_kelly);
  //std::cout << "Born ratio MC " << sigmaBorn_Ratio_MC << "\n";

  // Assuming that the discrepancy in n/p ratio we get by fitting MC to data
  // is totally attributable to Born cross section ratio. Hece we calculate the 
  // corrected Born cross section ratio by multiplying the MC Born cross sectio
  // ratio with the data/MC fit parameter "R".
  double sigmaBorn_Ratio_corr = sigmaBorn_Ratio_MC * R_fit;
  double sigmaBorn_Ratio_corr_error = abs(sigmaBorn_Ratio_MC) * R_fit_err;
  double sigmaBorn_Ratio_corr_error_sys = abs(sigmaBorn_Ratio_MC) * R_fit_err_sys;
  double sigmaBorn_Ratio_corr_error_tot = sqrt(pow(sigmaBorn_Ratio_corr_error,2) + pow(sigmaBorn_Ratio_corr_error_sys,2));
  

  // -- No TPE Check - beta phase [04/26/2024]
  double GMn_data_christy_ye = kine::ExtractGMn(etheta,Q2,GEp_christy,GMp_christy,GEn_ye,sigmaBorn_Ratio_corr); 
  double GMn_ov_munGD_data_christy_ye = GMn_data_christy_ye/(EMFFFits::GetGDip(Q2)*constant::mun);
  // calculating error
  double GMn_err_christy_ye = kine::CalcGMnError(etheta,Q2,GEp_christy,GMp_christy,GEn_ye,sigmaBorn_Ratio_corr,sigmaBorn_Ratio_corr_error); 
  double GMn_ov_munGD_err_christy_ye = GMn_err_christy_ye/abs(EMFFFits::GetGDip(Q2)*constant::mun);

  // Extract GMn from data using the corrected born CS ratio and the remaining EMFF
  // values got from some parametrization. Here we have used Z. Ye et al's parametrization (2017).
  // -- Primary
  double GMn_data_ye_ye = kine::ExtractGMn(etheta,Q2,GEp_ye,GMp_ye,GEn_ye,sigmaBorn_Ratio_corr); 
  double GMn_ov_munGD_data_ye_ye = GMn_data_ye_ye/(EMFFFits::GetGDip(Q2)*constant::mun);
  // calculating error
  double GMn_err_ye_ye = kine::CalcGMnError(etheta,Q2,GEp_ye,GMp_ye,GEn_ye,sigmaBorn_Ratio_corr,sigmaBorn_Ratio_corr_error); 
  double GMn_ov_munGD_err_ye_ye = GMn_err_ye_ye/abs(EMFFFits::GetGDip(Q2)*constant::mun);
  // -- Vary GEn w.r.t. Primary --
  // 1
  double GMn_data_ye_kelly = kine::ExtractGMn(etheta,Q2,GEp_ye,GMp_ye,GEn_kelly,sigmaBorn_Ratio_corr); 
  double GMn_ov_munGD_data_ye_kelly = GMn_data_ye_kelly/(EMFFFits::GetGDip(Q2)*constant::mun);
  // calculating error
  double GMn_err_ye_kelly = kine::CalcGMnError(etheta,Q2,GEp_ye,GMp_ye,GEn_kelly,sigmaBorn_Ratio_corr,sigmaBorn_Ratio_corr_error); 
  double GMn_ov_munGD_err_ye_kelly = GMn_err_ye_kelly/abs(EMFFFits::GetGDip(Q2)*constant::mun); 
  // 2
  double GMn_data_ye_gal = kine::ExtractGMn(etheta,Q2,GEp_ye,GMp_ye,GEn_gal,sigmaBorn_Ratio_corr); 
  double GMn_ov_munGD_data_ye_gal = GMn_data_ye_gal/(EMFFFits::GetGDip(Q2)*constant::mun);
  // calculating error
  double GMn_err_ye_gal = kine::CalcGMnError(etheta,Q2,GEp_ye,GMp_ye,GEn_gal,sigmaBorn_Ratio_corr,sigmaBorn_Ratio_corr_error); 
  double GMn_ov_munGD_err_ye_gal = GMn_err_ye_gal/abs(EMFFFits::GetGDip(Q2)*constant::mun); 
  // -- Vary sigma_p w.r.t. Primary --
  // 1
  double GMn_data_kelly_ye = kine::ExtractGMn(etheta,Q2,GEp_kelly,GMp_kelly,GEn_ye,sigmaBorn_Ratio_corr); 
  double GMn_ov_munGD_data_kelly_ye = GMn_data_kelly_ye/(EMFFFits::GetGDip(Q2)*constant::mun);
  // calculating error
  double GMn_err_kelly_ye = kine::CalcGMnError(etheta,Q2,GEp_kelly,GMp_kelly,GEn_ye,sigmaBorn_Ratio_corr,sigmaBorn_Ratio_corr_error); 
  double GMn_ov_munGD_err_kelly_ye = GMn_err_kelly_ye/abs(EMFFFits::GetGDip(Q2)*constant::mun);
  // 2
  double GMn_data_ar_ye = kine::ExtractGMn(etheta,Q2,GEp_ar,GMp_ar,GEn_ye,sigmaBorn_Ratio_corr); 
  double GMn_ov_munGD_data_ar_ye = GMn_data_ar_ye/(EMFFFits::GetGDip(Q2)*constant::mun);
  // calculating error
  double GMn_err_ar_ye = kine::CalcGMnError(etheta,Q2,GEp_ar,GMp_ar,GEn_ye,sigmaBorn_Ratio_corr,sigmaBorn_Ratio_corr_error); 
  double GMn_ov_munGD_err_ar_ye = GMn_err_ar_ye/abs(EMFFFits::GetGDip(Q2)*constant::mun);
  // --

  std::cout << "\n----------\n";
  std::cout << Form("Q2: %.1f, n/p ratio: %.4f +/- %.4f +/- %.4f, Tot err: %.4f\n\n",Q2,sigmaBorn_Ratio_corr,sigmaBorn_Ratio_corr_error,sigmaBorn_Ratio_corr_error_sys,sigmaBorn_Ratio_corr_error_tot);
  std::cout << "No TPE Corr (Beta): Christy, Ye \n";
  std::cout << Form("etheta: %.1f, Q2: %.1f, GMn: %.4f, GMn/(mun*GD): %.4f +/- %.4f\n",etheta,Q2,GMn_data_christy_ye,GMn_ov_munGD_data_christy_ye,GMn_ov_munGD_err_christy_ye);
  std::cout << "----- || ----- \n";
  std::cout << "Primary: Ye, Ye \n";
  std::cout << Form("etheta: %.1f, Q2: %.1f, GMn: %.4f, GMn/(mun*GD): %.4f +/- %.4f\n",etheta,Q2,GMn_data_ye_ye,GMn_ov_munGD_data_ye_ye,GMn_ov_munGD_err_ye_ye);
  std::cout << "Vary GEn w.r.t. Primary: \n";
  std::cout << "1. Ye, Kelly \n";
  std::cout << Form("etheta: %.1f, Q2: %.1f, GMn: %.4f, GMn/(mun*GD): %.4f +/- %.4f\n",etheta,Q2,GMn_data_ye_kelly,GMn_ov_munGD_data_ye_kelly,GMn_ov_munGD_err_ye_kelly);
  std::cout << "2. Ye, Galster \n";
  std::cout << Form("etheta: %.1f, Q2: %.1f, GMn: %.4f, GMn/(mun*GD): %.4f +/- %.4f\n",etheta,Q2,GMn_data_ye_gal,GMn_ov_munGD_data_ye_gal,GMn_ov_munGD_err_ye_gal);
  std::cout << "Vary sigma_p w.r.t. Primary: \n";
  std::cout << "1. Kelly, Ye \n";
  std::cout << Form("etheta: %.1f, Q2: %.1f, GMn: %.4f, GMn/(mun*GD): %.4f +/- %.4f\n",etheta,Q2,GMn_data_kelly_ye,GMn_ov_munGD_data_kelly_ye,GMn_ov_munGD_err_kelly_ye);
  std::cout << "2. Arrington, Ye \n";
  std::cout << Form("etheta: %.1f, Q2: %.1f, GMn: %.4f, GMn/(mun*GD): %.4f +/- %.4f\n",etheta,Q2,GMn_data_ar_ye,GMn_ov_munGD_data_ar_ye,GMn_ov_munGD_err_ar_ye);
  std::cout << "----------\n";

  return 0;
}
