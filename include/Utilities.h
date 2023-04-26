#ifndef PD_UTIL_H
#define PD_UTIL_H

#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <dirent.h>

#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"

#include "../include/CodaRun.h"
#include "../include/SimuJob.h"
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
		   int verbose,               // verbosity
		   vector<CodaRun> &crun);    // Output: Vector of CodaRun structs

  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   int verbose,               // verbosity
		   vector<CodaRun> &crun);    // Output: Vector of CodaRun structs

  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   int bbmag,                 // BB magnet current (in %)
		   int verbose,               // verbosity
		   vector<CodaRun> &crun);    // Output: Vector of CodaRun structs

  /* ####################################################
     ## Functions to read ROOT files by sorted segment ##  
     #################################################### */
  int LoadROOTTree(std::string path,          // ROOT file directory path
		   std::vector<CodaRun> crun, // CodaRun objects with run info
		   bool sort,                 // Sort by segments before parsing?
		   int verbose,               // verbosity
		   TChain* &C);               // Output: TChain with data

  int LoadROOTTree(std::string path,          // ROOT file directory path
		   CodaRun crun,              // CodaRun object with run info
		   bool sort,                 // Sort by segments before parsing?
		   int verbose,               // verbosity
		   TChain* &C);               // Output: TChain with data

  /* ##################################################
     ##   Function to read MC replay summary files   ##  
     ################################################## */
  /* Reads simulation job specifics from summary files and loads the values to SimuJob objects.
     This function has the standard naming conventions of output simulation and summary files 
     hard coded. Please make sure the files to analyze have names compatible with the standard
     naming convention for successful execution.
     ---------
     Standard naming conventions:
     1. Filebase: sbs<sbsconf>_sbs<sbsmag> for g4sbs gen., sbs<sbsconf>_sbs<sbsmag>_simc for SIMC gen.
     2. MC job summary file: <filebase>_summary.csv
     3. Digitized ROOT file: <filebase>_<process>_job_<jobid>.root
     4. Replayed digitized ROOT file: replayed_<filebase>_<process>_job_<jobid>.root */
  void ReadSimuJobSummary(std::string logfile_dir,   // Dir. path containing MC summary files
			  int &njobs,                // # jobs to analyze per process
			  bool issimcgen,            // True=>SIMC genrated events
			  int sbsconf,               // SBS configuration
			  int sbsmag,                // SBS magnet current (in %)
			  std::string target,        // target type
			  int verbose,               // verbosity
			  vector<SimuJob> &sjobs);   // Output: Vector of SimuJob objects

  void LoadSimuROOTTree(std::vector<SimuJob> sjobs,  // SimuJob objects with run info
			int verbose,                 // verbosity
			TChain* &C);                 // Output: TChain with data

  /* ###############################
     ## General Purpose Functions ##  
     ############################### */
  /* Acknowledgement: SplitString is based on a function written by David Flay.*/ 
  int SplitString(const char delim,                // delimiter
		  const std::string myStr,         // input string
		  std::vector<std::string> &out);  // output sub-strings

  double GetTotCharge(std::vector<CodaRun> cruns); // Calc. tot. charge from CodaRun objects
  double GetTotCharge(std::vector<SimuJob> sjobs); // Calc. tot. charge from SimuJob objects

  double GetTotNtries(std::vector<SimuJob> sjobs); // Calc. tot. # tries from SimuJob objects

  /* Calc. # tries & tot. Ch. from SimuJob objects */
  void GetTotNtriesnCh(std::vector<SimuJob> sjobs, // Input: List of SimuJob objects
		       vector<double> &data);      // Output: data[0]=>Tot. Ntries, data[1]=>Tot. Charge

}
#endif
