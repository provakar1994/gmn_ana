#ifndef PD_UTIL_H
#define PD_UTIL_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"

#include "../include/CodaRun.h"
#include "../include/ExpConstants.h"

namespace util_pd {

  /* #################################################
     ##                HCAL Histograms              ##  
     ################################################# */
  TH2F *TH2FHCALface_rc(std::string name);      // returns TH2F for HCAL face (row,col)
  TH2F *TH2FHCALface_xy_data(std::string name); // returns TH2F for HCAL face (x,y) [Data]
  TH2F *TH2FHCALface_xy_simu(std::string name); // returns TH2F for HCAL face (x,y) [Simu]
  TH2F *TH2FdxdyHCAL(std::string name);         // returns TH2F for dxdyHCAL

  // draws rectangular cut regions
  void DrawArea(vector<double> dimensions,      // a vector with extreme points
		int lcolor,  // Default = 2 
		int lwidth,  // Default = 4
		int lstyle); // Default = 9


  /* #################################################
     ##              Kinematic Histograms           ##  
     ################################################# */
  TH1F *TH1FhW(std::string name);   // returns W histogram
  TH1F *TH1FhQ2(std::string name,   // returns Q2 histogram
		int conf);   // SBS config

  /* ###########################################################
     ##   Function to standard read CSV files with run info   ##  
     ########################################################## */
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   vector<CodaRun> &crun);    // Output: Vector of CodaRun structs

  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   vector<CodaRun> &crun);    // Output: Vector of CodaRun structs

  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   int bbmag,                 // BB magnet current (in %)
		   vector<CodaRun> &crun);    // Output: Vector of CodaRun structs
}

#endif
