#include "../include/Utilities.h"

namespace util_pd {

  /* #################################################
     ##                HCAL Related                 ##  
     ################################################# */
  //_____________________________________
  TH2F *TH2FHCALface_rc(std::string name) {
    // returns TH2F for HCAL face (row,col)
    /* NOTE: HCAL block id (ibblk) starts from 1 and goes up to 288 but both
       HCAL row (rowblk) and column starts from 0 and goes up to 23 and 
       11, respectively. Extremely annoying! */
    TH2F *h = new TH2F(name.c_str(), ";HCAL columns;HCAL rows",
		       expconst::hcalcol, 0, expconst::hcalcol,
		       expconst::hcalrow, 0, expconst::hcalrow);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_data(std::string name) {
    // returns TH2F for HCAL face (x,y) [Data]
    double y_min = expconst::yHCAL_r_DB - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB + expconst::hcalblk_h/2.;
    TH2F *h = new TH2F(name.c_str(), ";yHCAL_{exp} (m);xHCAL_{exp} (m)",
		       expconst::hcalcol, y_min, y_max,
		       expconst::hcalrow, x_min, x_max);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_simu(std::string name) {
    // returns TH2F for HCAL face (x,y) [Simu]
    double y_min = expconst::yHCAL_r_DB_MC - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB_MC + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB_MC - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB_MC + expconst::hcalblk_h/2.;
    TH2F *h = new TH2F(name.c_str(), ";yHCAL_{exp} (m);xHCAL_{exp} (m)",
		       expconst::hcalcol, y_min, y_max,
		       expconst::hcalrow, x_min, x_max);
    return h;
  }
  //_____________________________________
  TH2F *TH2FdxdyHCAL(std::string name) {
    // returns TH2F for dxdyHCAL
    TH2F *h = new TH2F(name.c_str(), "; yHCAL_{obs} - yHCAL_{exp} (m); xHCAL_{obs} - xHCAL_{exp} (m)",
		       250, -1.25, 1.25, 250, -3.5, 2);
    return h;
  }
  //_____________________________________
  void DrawArea(vector<double> dimensions, int lcolor=2, int lwidth=4, int lstyle=9) {
    /* Draws four lines to represent a rectangular cut area */
    double top = dimensions[0];                 // -X axis
    double bottom = dimensions[1];              // +X axis
    double right = dimensions[2];               // -Y axis
    double left = dimensions[3];                // +Y axis
    TLine line;
    line.SetLineColor(lcolor); 
    line.SetLineWidth(lwidth); 
    line.SetLineStyle(lstyle);
    line.DrawLine(right, bottom, left, bottom); // bottom margin
    line.DrawLine(right, top, left, top);       // top margin
    line.DrawLine(right, top, right, bottom);   // right margin
    line.DrawLine(left, top, left, bottom);     // left margin
  }


  /* #################################################
     ##              Kinematic Histograms           ##  
     ################################################# */
  TH1F *TH1FhW(std::string name) {
    // returns W histogram
    TH1F *h = new TH1F(name.c_str(), "W Distribution (GeV)", 250,0,2);
    return h;
  }
  TH1F *TH1FhQ2(std::string name,       // Name of histogram
		int conf) {             // SBS config
    // returns Q2 histogram
    int nbin=0; double hmin=-100, hmax=-100;
    if (conf==4) { nbin=100; hmin=1.; hmax=4.; } 
    else if (conf==14) { nbin=100; hmin=5.; hmax=10.; }
    else if (conf==7) { nbin=120; hmin=6.; hmax=12.; }
    else if (conf==11) { nbin=200; hmin=8.; hmax=18.; }
    else cerr << "[Utilities::TH1FhQ2] Enter valid SBS config!!" << endl;
    TH1F *h = new TH1F(name.c_str(), "Q^{2} Distribution (GeV^{2})", 
		       nbin, hmin, hmax);
    return h;
  }

  /* #################################################
     ##   Function to read CSV file with run info   ##  
     ################################################# */
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   vector<CodaRun> &crun)     // Output: Vector of CodaRun structs
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/good_runList_GMn_nTPE_"; 
    std::string mid = "_pass_";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    ifstream run_data; run_data.open(run_spreadsheet);
    string readline;
    if(run_data.is_open()){
      std::cout << std::endl 
		<< "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;
      string skip_header; getline(run_data, skip_header);  // skipping column header
      crun.clear();
      while(getline(run_data,readline)){                   // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){       // reading each element of a line
	  string temptoken=token;
	  temp.push_back(temptoken);
	}
	// add relevant info to CodaRun objects
	if (stoi(temp[0]) == sbsconf) {
	  if (crun.size() >= nruns) break;
	  CodaRun temp_cr;
	  temp_cr.SetDataRunSheet(temp);
	  crun.push_back(temp_cr);

	}

	temp.clear();
      }
      // let's update nruns with total no. of runs to analyze
      nruns = crun.size();
    }else{
      std::cerr << " **!**[Utilities::ReadRunList] Error - No file named: " << run_spreadsheet  << std::endl;
      throw;
    }
    run_data.close();
  }
  //_____________________________________
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   vector<CodaRun> &crun)     // Output: Vector of CodaRun structs
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/good_runList_GMn_nTPE_"; 
    std::string mid = "_pass_";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // convert magnet field values from % to A
    sbsmag *= 21;    // 100% SBS magnet current = 2100 A

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    ifstream run_data; run_data.open(run_spreadsheet);
    string readline;
    if(run_data.is_open()){
      std::cout << std::endl 
		<< "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;
      string skip_header; getline(run_data, skip_header); // skipping column header
      while(getline(run_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;
	  temp.push_back(temptoken);
	}
	// add relevant info to CodaRun objects
	if (stoi(temp[0]) == sbsconf && 
	    stoi(temp[3]) == sbsmag) {
	  if (crun.size() >= nruns) break;
	  CodaRun temp_cr;
	  temp_cr.SetDataRunSheet(temp);
	  crun.push_back(temp_cr);
	}

	temp.clear();
      }
      // let's update nruns with total no. of runs to analyze
      nruns = crun.size();
    }else{
      std::cerr << " **!**[Utilities::ReadRunList] Error - No file named: " << run_spreadsheet  << std::endl;
      throw;
    }
    run_data.close();
  }
  //_____________________________________
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   int bbmag,                 // BB magnet current (in %)
		   vector<CodaRun> &crun)     // Output: Vector of CodaRun structs
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/good_runList_GMn_nTPE_"; 
    std::string mid = "_pass_";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // convert magnet field values from % to A
    sbsmag *= 21;    // 100% SBS magnet current = 2100 A
    bbmag *= 7.5;    // 100% BB magnet current = 750 A

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    ifstream run_data; run_data.open(run_spreadsheet);
    string readline;
    if(run_data.is_open()){
      std::cout << std::endl 
		<< "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;     
      string skip_header; getline(run_data, skip_header); // skipping column header
      while(getline(run_data,readline)){                  // reading each line
	istringstream tokenStream(readline);
	string token;
	char delimiter = ',';
	vector<string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  string temptoken=token;
	  temp.push_back(temptoken);
	}
	// add relevant info to CodaRun objects
	if (stoi(temp[0]) == sbsconf 
	    && stoi(temp[3]) == sbsmag 
	    && stoi(temp[4]) == bbmag) {
	  if (crun.size() >= nruns) break;
	  CodaRun temp_cr;
	  temp_cr.SetDataRunSheet(temp);
	  crun.push_back(temp_cr);
	}

	temp.clear();
      }
      // let's update nruns with total no. of runs to analyze
      nruns = crun.size();
    }else{
      std::cerr << " **!**[Utilities::ReadRunList] Error - No file named: " << run_spreadsheet  << std::endl;
      throw;
    }
    run_data.close();
  }

}
