#include <iostream>

#include "EMFFFits.h"

int extractGMn(double const etheta, double const Q2, double const R_fit) {

  // EMFF fits
  Ye2017 yefit;
  Kelly2004 kellyfit;
  Seamus20XX seamusfit;
  Christy2022 christyfit;

  // EMFF extraction using same parametrization used in SIMC generator
  double GEp_kelly = kellyfit.GetFF(G_t::kGEp,Q2);
  double GMp_kelly = kellyfit.GetFF(G_t::kGMp,Q2);
  double GEn_seamus = seamusfit.GetFF(G_t::kGEn,Q2);
  double GMn_kelly = kellyfit.GetFF(G_t::kGMn,Q2);

  // EMFF extraction using Ye fit (2017), to get GMn from data/MC fit
  double GEp_christy = christyfit.GetFF(G_t::kGEp,Q2);
  double GMp_christy = christyfit.GetFF(G_t::kGMp,Q2);
  double GEn_ye = yefit.GetFF(G_t::kGEn,Q2);

  // Getting sigma Born ratio (n/p) for the given Q2 and etheta 
  // (basically, mimicing the calculation done in MC generator)
  double sigmaBorn_Ratio_MC = kine::sigmaBorn_ratio(etheta,Q2,GEp_kelly,GMp_kelly,GEn_seamus,GMn_kelly);
  //std::cout << "Born ratio MC " << sigmaBorn_Ratio_MC << "\n";

  // Assuming that the discrepancy in n/p ratio we get by fitting MC to data
  // is totally attributable to Born cross section ratio. Hece we calculate the 
  // corrected Born cross section ratio by multiplying the MC Born cross sectio
  // ratio with the data/MC fit parameter "R".
  double sigmaBorn_Ratio_corr = sigmaBorn_Ratio_MC * R_fit;

  // Extract GMn from data using the corrected born CS ratio and the remaining EMFF
  // values got from some parametrization. Here we have used Z. Ye et al's parametrization (2017).
  double GMn_data = kine::ExtractGMn(etheta,Q2,GEp_christy,GMp_christy,GEn_ye,sigmaBorn_Ratio_corr); 
  double GMn_ov_munGD_data = GMn_data/(EMFFFits::GetGDip(Q2)*constant::mun);

  std::cout << "\n----------\n";
  std::cout << Form("etheta: %.1f, Q2: %.1f, GMn: %.4f, GMn/(mun*GD): %.4f\n",etheta,Q2,GMn_data,GMn_ov_munGD_data);
  std::cout << "----------\n";

  return 0;
}
