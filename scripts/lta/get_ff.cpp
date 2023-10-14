#include <iostream>

//#include "../include/gmn_ana.h"
#include "../include/EMFFFits.h"

// Function to convert string to enum class
G_t stringToEnum(const std::string& str) {
    if (str == "GEp") {
        return G_t::kGEp;
    } else if (str == "GMp") {
        return G_t::kGMp;
    } else if (str == "GEn") {
        return G_t::kGEn;
    } else if (str == "GMn") {
        return G_t::kGMn;
    } else {
        throw std::invalid_argument("Invalid string for conversion to enum");
    }
}

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

  std::cout << "Ye2017: " << G << ": " << yefit.GetFF(stringToEnum(G),Q2) << std::endl;
  std::cout << "Kelly2004: " << G << ": " << kellyfit.GetFF(stringToEnum(G),Q2) << std::endl;
  //std::cout << yefit.GetFFwErr(stringToEnum(G),Q2)[0] << std::endl;

  // std::cout << yefit.GetGDip(Q2) << std::endl;
  // std::cout << kellyfit.GetGDip(Q2) << std::endl;
  // std::cout << seamusfit.GetGDip(Q2) << std::endl;

  return 0;
}
