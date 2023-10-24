#include <iostream>

//#include "../include/gmn_ana.h"
#include "EMFFFits.h"

int get_ff(double const Q2, std::string G) {

  Ye2017 yefit;
  Kelly2004 kellyfit;
  Seamus20XX seamusfit;
  Christy2022 christyfit;
  Galster1971 galsterfit;

  G_t kG = EMFFFits::StrToG_t(G);

  // std::cout << yefit.GetFFwErr(stringToEnum(G),Q2)[0] << std::endl;
  // std::cout << yefit.GetFFwErr(stringToEnum(G),Q2)[1] << std::endl;

  // std::cout << kellyfit.GetFFwErr(stringToEnum(G),Q2)[0] << std::endl;
  // std::cout << kellyfit.GetFFwErr(stringToEnum(G),Q2)[1] << std::endl;

  // std::cout << seamusfit.GetFFwErr(stringToEnum(G),Q2)[0] << std::endl;
  // std::cout << seamusfit.GetFFwErr(stringToEnum(G),Q2)[1] << std::endl;

  // std::cout << "Ye2017: " << G << ": " << yefit.GetFF(kG,Q2) << std::endl;
  // std::cout << "Kelly2004: " << G << ": " << kellyfit.GetFF(kG,Q2) << std::endl;
  // std::cout << "Christy2022: " << G << ": " << christyfit.GetFF(kG,Q2) << std::endl;

  // std::cout << EMFFFits::GetGDip(Q2) << std::endl;
  // std::cout << kellyfit.GetGDip(Q2) << std::endl;
  // std::cout << seamusfit.GetGDip(Q2) << std::endl;

  TString outFile;
  outFile = Form("%s_parametrization.csv",G.c_str());
  ofstream outData;
  outData.open(outFile);
  for (int i=0; i<141; i++) {
    double Q2 = 1; Q2 += 0.1*i;
    if (kG == G_t::kGEp) {
      //std::cout << Q2 << "," << kellyfit.GetFF(kG,Q2) << "," << yefit.GetFF(kG,Q2) << "," << christyfit.GetFF(kG,Q2) << "\n";
      if (i==0) outData << "Q2" << "," << "Kelly2004" << "," << "Ye2017" << "," << "Christy2022" << "\n";
      outData << Q2 << "," << kellyfit.GetFF(kG,Q2) << "," << yefit.GetFF(kG,Q2) << "," << christyfit.GetFF(kG,Q2) << "\n";
    } else if (kG == G_t::kGMp) {
      //std::cout << Q2 << "," << kellyfit.GetFF(kG,Q2) << "," << yefit.GetFF(kG,Q2) << "," << christyfit.GetFF(kG,Q2) << "\n";
      if (i==0) outData << "Q2" << "," << "Kelly2004" << "," << "Ye2017" << "," << "Christy2022" << "\n";
      outData << Q2 << "," << kellyfit.GetFF(kG,Q2) << "," << yefit.GetFF(kG,Q2) << "," << christyfit.GetFF(kG,Q2) << "\n";
    } else if (kG == G_t::kGEn) {
      //std::cout << Q2 << "," << kellyfit.GetFF(kG,Q2) << "," << yefit.GetFF(kG,Q2) << "," << seamusfit.GetFF(kG,Q2) << "," << galster.GetFF(kG,Q2) << "\n";
      if (i==0) outData << "Q2" << "," << "Kelly2004" << "," << "Ye2017" << "," << "Seamus20XX" << "," << "Galster1971" << "\n";
      outData << Q2 << "," << kellyfit.GetFF(kG,Q2) << "," << yefit.GetFF(kG,Q2) << "," << seamusfit.GetFF(kG,Q2) << "," << galsterfit.GetFF(kG,Q2) << "\n";
    } else {
      //std::cout << Q2 << "," << kellyfit.GetFF(kG,Q2) << "," << yefit.GetFF(kG,Q2) << "\n";
      if (i==0) outData << "Q2" << "," << "Kelly2004" << "," << "Ye2017" << "\n";
      outData << Q2 << "," << kellyfit.GetFF(kG,Q2) << "," << yefit.GetFF(kG,Q2) << "\n";
    }
  }

  return 0;
}
