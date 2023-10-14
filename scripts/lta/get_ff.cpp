#include <iostream>

//#include "../include/gmn_ana.h"
#include "EMFFFits.h"

int get_ff(double const Q2, std::string G) {

  Ye2017 yefit;
  Kelly2004 kellyfit;
  Seamus20XX seamusfit;

  // std::cout << yefit.GetFFwErr(stringToEnum(G),Q2)[0] << std::endl;
  // std::cout << yefit.GetFFwErr(stringToEnum(G),Q2)[1] << std::endl;

  // std::cout << kellyfit.GetFFwErr(stringToEnum(G),Q2)[0] << std::endl;
  // std::cout << kellyfit.GetFFwErr(stringToEnum(G),Q2)[1] << std::endl;

  // std::cout << seamusfit.GetFFwErr(stringToEnum(G),Q2)[0] << std::endl;
  // std::cout << seamusfit.GetFFwErr(stringToEnum(G),Q2)[1] << std::endl;

  std::cout << "Ye2017: " << G << ": " << yefit.GetFF(EMFFFits::StrToG_t(G),Q2) << std::endl;
  std::cout << "Kelly2004: " << G << ": " << kellyfit.GetFF(EMFFFits::StrToG_t(G),Q2) << std::endl;

  // std::cout << EMFFFits::GetGDip(Q2) << std::endl;
  // std::cout << kellyfit.GetGDip(Q2) << std::endl;
  // std::cout << seamusfit.GetGDip(Q2) << std::endl;

  return 0;
}
