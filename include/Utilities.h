#ifndef PD_UTIL_H
#define PD_UTIL_H

#include <vector>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <dirent.h>
#include <algorithm>

#include "TF1.h"
#include "TBox.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLine.h"
#include "TChain.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TString.h"

#include "CodaRun.h"
#include "SimuJob.h"
#include "Constants.h"
#include "ExpConstants.h"

namespace util_pd {

  /* #########################################
     ##                General              ##  
     ######################################### */
  // returns TCanvas object with optimized size
  TCanvas *TC(std::string hname,  // name of the canvas
 	      int rdiv,           // # divisions in row
	      int cdiv);          // # divisions in column 

  // splits a given canvas into two pads suitable for pull plots
  std::vector<TPad*> GetPadsForPullPlot(TCanvas *c1);

  // draws a horizontal line at y=0
  void DrawZeroLine(TPad *p1, double xmin, double xmax);

  std::string getDate(); // returns today's date

  // creates a TH1F object from a given TF1 object
  void TF1toTH1F(TF1* const func, TH1F* &hist);

  // Returns normalized coordinate for a given x value
  double GetxNDC(double x);

  // Returns y value for a given y coordinate in NDC
  double UnfoldyNDC(double yNDC);

  // Plots cut region for given x range
  void PlotCutRegion(double xmin, double xmax);

  /* #################################################
     ##                HCAL Histograms              ##  
     ################################################# */
  TH2F *TH2FHCALface_rc(std::string hname);      // returns TH2F for HCAL face (row,col)
  TH2F *TH2FHCALface_xy_data(std::string hname, double sbs_kick, int rpass); // returns TH2F for HCAL face (x,y) [Data]
  TH2F *TH2FHCALface_xy_simu(std::string hname, double sbs_kick); // returns TH2F for HCAL face (x,y) [Simu]
  TH2F *TH2FdxdyHCAL(std::string hname);         // returns TH2F for dxdyHCAL

  // draws rectangular cut regions
  void DrawArea(std::vector<double> dimensions,      // a vector with extreme points
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
		   std::vector<CodaRun> &crun);    // Output: Vector of CodaRun structs

  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   int verbose,               // verbosity
		   std::vector<CodaRun> &crun);    // Output: Vector of CodaRun structs

  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   int bbmag,                 // BB magnet current (in %)
		   int verbose,               // verbosity
		   std::vector<CodaRun> &crun);    // Output: Vector of CodaRun structs

  /* ####################################################
     ## Functions to read ROOT files by sorted segment ##  
     #################################################### */
  int LoadROOTTree(std::string path,          // ROOT file directory path
		   std::vector<CodaRun> crun, // CodaRun objects with run info
		   bool sort,                 // Sort by segments before parsing?
		   int verbose,               // verbosity
		   TChain* &C);               // Output: TChain with data

  int LoadROOTTree(std::string path,          // ROOT file directory path
		   CodaRun &crun,             // CodaRun object with run info
		   bool sort,                 // Sort by segments before parsing?
		   int verbose,               // verbosity
		   bool is_read_log,          // want to read corresponding log files?
		   TChain* &C);               // Output: TChain with data

  /* ##################################################
     ##   Function to read MC replay summary files   ##  
     ################################################## */
  /* Reads simulation job specifics from summary files and loads the values to SimuJob objects.
     This function has the standard naming conventions of output simulation and summary files 
     hard coded. Please make sure the files to analyze have names compatible with the standard
     naming convention for successful execution.
     NOTE: The presence of a summary file is must for this process to work. In case of SIMC 
     D(ee'N) process, the ROOT file directory must have two summary files presents. One for 
     the deep process and one for the deen process. Each file must have at least one entry.
     ---------
     Standard naming conventions:
     1. Filebase: sbs<sbsconf>_sbs<sbsmag>_<generator>
     2. MC job summary file: <filebase>_summary.csv
     3. Digitized ROOT file: <filebase>_<process>_job_<jobid>.root
     4. Replayed digitized ROOT file: replayed_<filebase>_<process>_job_<jobid>.root */
  void ReadSimuJobSummary(std::string logfile_dir,   // Dir. path containing MC summary files
			  std::string prefix,        // prefix to standard filebase (Special case)
			  int sbsconf,               // SBS configuration
			  int sbsmag,                // SBS magnet current (in %)
			  std::string generator,     // simc / g4sbs
			  std::string process,       // reaction process being simulated
			  int &njobs,                // # jobs to analyze per process
			  int verbose,               // verbosity
			  std::vector<SimuJob> &sjobs);   // Output: Vector of SimuJob objects

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
		       std::vector<long double> &data); // Output: data[0]=>Tot. Ntries, data[1]=>Tot. Charge

  /* Calc. # tries & tot. Ch. from SimuJob objects */
  void GetTotNtriesnCh(std::vector<SimuJob> sjobs, // Input: List of SimuJob objects
		       std::string process,        // Input: Reaction process
		       std::vector<long double> &data); // Output: data[0]=>Tot. Ntries, data[1]=>Tot. Charge

  /* Calculates mean energy loss in the target before and after scattering  */
  void GetElossInTgt(std::string const target,  // Target type
		     int const rnum,            // Run number
		     int const sbsconf,         // SBS configuration
		     double const vz,           // m, vertex z co-ordinate
		     double const etheta,       // rad, scattering angle
		     int const verbose,         // verbosity
		     std::vector<double> &data);// Output: data[0]=>before scattering, data[1]=>after scattering,

  /* Luminosity calculation for g4sbs data */
  double Luminosity(double ibeam, std::string targetType);               // ibeam shoud be in A

  /* #########################################
     ## General Purpose Templated Functions ##  
     ######################################### */
  /* Function to sort array indices by their values in descending order */
  template <typename T>
    std::vector<size_t> SortIndices(const T* array, size_t size) {
    // Initialize original indices
    std::vector<size_t> indices(size);
    for (size_t i = 0; i < size; ++i) {
      indices[i] = i;
    }
    // Sort indices based on array values
    std::sort(indices.begin(), indices.end(),
              [array](size_t i1, size_t i2) { return array[i1] > array[i2]; });
    return indices;
  }

}
#endif
