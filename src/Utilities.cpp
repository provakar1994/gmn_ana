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
		   int verbose,               // verbosity
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
      std::cout << "Reading run info from: "<< run_spreadsheet 
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
      if (verbose == 1) {
	std::cout << "First run info:" << std::endl << crun[0];
	std::cout << "Last run info:" << std::endl << crun[nruns-1];
      }
    }else
      throw std::runtime_error("[util_pd::ReadRunList] Run spreadsheet doesn't exist");
    run_data.close();
  }
  //_____________________________________
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int sbsmag,                // SBS magnet current (in %)
		   int verbose,               // verbosity
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
      std::cout << "Reading run info from: "<< run_spreadsheet 
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
      if (verbose == 1) {
	std::cout << "First run info:" << std::endl << crun[0];
	std::cout << "Last run info:" << std::endl << crun[nruns-1];
      }
    }else
      throw std::runtime_error("[util_pd::ReadRunList] Run spreadsheet doesn't exist");
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
		   int verbose,               // verbosity
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
      std::cout << "Reading run info from: "<< run_spreadsheet 
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
      if (verbose == 1) {
	std::cout << "First run info:" << std::endl << crun[0];
	std::cout << "Last run info:" << std::endl << crun[nruns-1];
      }
    }else
      throw std::runtime_error("[util_pd::ReadRunList] Run spreadsheet doesn't exist");
    run_data.close();
  }

  /* ####################################################
     ## Functions to read ROOT files by sorted segment ##  
     #################################################### */
  //______________________________________________________________________________
  void ListDirectory(const char *path,std::vector<std::string> &list){
    struct dirent *entry;
    DIR *dir = opendir(path);
    if(dir==NULL){
      return;
    }
    std::string myStr;
    while ((entry = readdir(dir)) != NULL) {
      // std::cout << " mystr= " << myStr << std::endl;
      myStr = entry->d_name;
      list.push_back(myStr);
    }
    closedir(dir);
  }
  //______________________________________________________________________________
  int SplitString(const char delim, const std::string myStr, std::vector<std::string> &out){
    // split a string by a delimiter
    std::stringstream ss(myStr);
    std::vector<std::string> result;
    while( ss.good() ){
      std::string substr;
      std::getline(ss, substr, delim);
      out.push_back(substr);
    }
    return 0;
  }
  //______________________________________________________________________________
  bool sortbyval(const pair<int, int> &a, const pair<int, int> &b) {
    return (a.first < b.first);
  } 
  //______________________________________________________________________________
  int GetROOTFileMetaData(const char *rfDirPath, int run,
			  std::vector<int> &data, 
			  std::vector<pair<int, int>> &segB_segE,
			  int verbose){
    /* Acknowledgement: Based on a function written by David Flay.
    Determines the beginning and end segment number for a CODA run
    input: 
    - rfDirPath: path to where the ROOT file(s) are located 
    - run: CODA run number 
    output: vector containing:  
    - stream number of the EVIO file associated with the run  
    - number of files associated with the run
    - list of sengment numbers of the run (sorted) */

    int rc=-1; // assume fail 

    // first get list of files in the directory 
    std::vector<std::string> fileList;
    ListDirectory(rfDirPath,fileList); 

    if (verbose > 1) {
      std::cout << "----" << std::endl;
      std::cout << " fileList.size() " << fileList.size() << std::endl;
      for (int i = 0; i < fileList.size(); i++) {
	std::cout << " fileList[" << i << "] = " << fileList[i] << std::endl;
      }
      std::cout << "----" << std::endl;
    }

    // skip first two entries since we get . and .. at top of list
    for(int i=0;i<2;i++) fileList.erase(fileList.begin()); 

    // identify the right index, and find number of files 
    int j=-1, fileCnt=0, theRun=0, stream=0;
    const int NF = fileList.size();
    std::vector<std::string> o1;
    for(int i=0; i<NF; i++){ 

      // lets implement a filter 
      if (fileList[i].find("e1209019_fullreplay") == 0) {

	// split each entry based on a sufficient delimiter 
	SplitString('_', fileList[i], o1);
	// 3rd entry (index 2) is the one that has the run number 
	theRun = std::atoi(o1[2].c_str());
	if(theRun==run){
	  fileCnt++; 

	  // split each entry based on a sufficient delimiter 
	  std::vector<std::string> o2, o3; 
	  // determine the stream (index 3) 
	  std::string theStr = o1[3]; 
	  SplitString('m', theStr, o2); 
	  stream = std::atoi(o2[1].c_str());
	  o2.clear();
	  // determine the segment numbers (index 4 and 5) 
	  theStr = o1[4]; 
	  SplitString('g', theStr, o2); 
	  int bseg = std::atoi(o2[1].c_str());
	  int eseg = std::atoi(o1[5].c_str());
	  segB_segE.push_back(make_pair(bseg,eseg));
	  o2.clear();
	}
	o1.clear();
      } 
    }

    // lets sort the segments by begging segment number 
    // this is necessary for proper beam charge calculation
    sort(segB_segE.begin(), segB_segE.end(), sortbyval);  

    if (fileCnt==0) {
      std::cout << std::endl << "--!!--" << std::endl  
		<< "WARNING! [util_pd::GetROOTFileMetaData]: ROOT file directory is empty!" 
		<< std::endl << "--!!--" << std::endl;
      return -1;
    }

    data.push_back(stream); 
    data.push_back(fileCnt);
 
    return 0;
  } 
  //______________________________________________________________________________
  int LoadROOTTree(std::string path,
		   std::vector<CodaRun> crun,    // multiple CODA runs
		   bool sort,
		   int verbose,
		   TChain* &C) 
  {
    
    if (!crun.empty()) {
      const int nruns = crun.size();
      std::cout << "Parsing ROOT files from " << nruns << " runs.." << std::endl;
      if (!sort) {
	for (int i=0; i<nruns; i++) {
	  std::string rfname = Form("%s/*%d*",path.c_str(),crun[i].runnum);
	  if (verbose > 1) std::cout << rfname << std::endl;
	  C->Add(rfname.c_str());
	}
      } else {
	std::cout << "Sorting by segments.." << std::endl;
	int aRun, aNumFiles, aStream, rc;
	std::vector<int> md;
	std::vector<pair<int, int>> segB_segE;
	// Looping through unique runs
	for (int irun=0; irun<nruns; irun++) {
	  aRun = crun[irun].runnum;
	  rc = GetROOTFileMetaData(path.c_str(),aRun,md,segB_segE,0);
	  if (rc==0) { // non-zero number of segments
	    aStream   = md[0]; //stream no.
	    aNumFiles = md[1]; //total # segments
	    if (verbose > 0) {
	      std::cout << "----" << std::endl;
	      std::cout << Form(" Run %d, Total # of segments %d",aRun,aNumFiles) << std::endl;
	      std::cout << Form(" Beg seg %d-%d, End seg %d-%d",segB_segE[0].first,segB_segE[0].second,
				segB_segE[aNumFiles-1].first,segB_segE[aNumFiles-1].second) << std::endl;
	      std::cout << "----" << std::endl;
	    }
	    // Looping through segments and adding to tree
	    for (int iseg=0; iseg<aNumFiles; iseg++) {
	      std::string rfname = Form("%s/e1209019_fullreplay_%d_stream%d_seg%d_%d.root",path.c_str(),
					aRun,aStream,segB_segE[iseg].first,segB_segE[iseg].second);
	      if (verbose > 1) std::cout << rfname << std::endl;
	      C->Add(rfname.c_str());
	    }	
	  } else {
	    std::cout << std::endl << "--!!--" << std::endl  
		      << "WARNING! [util_pd::LoadROOTTree]: No ROOT file exists for run " << crun[irun].runnum << "!" 
		      << std::endl << "--!!--" << std::endl;
	    return -1;
	  }
	  // getting ready for next run
	  md.clear();
	  segB_segE.clear();
	}
      }
      if (C->GetEntries()==0) 
	throw std::runtime_error("[util_pd::LoadROOTTree] Empty ROOT files Or, they don't exist!");
    }else 
      throw std::runtime_error("[util_pd::LoadROOTTree] CODA run list is empty!");

    return 0;
  }
  //______________________________________________________________________________
  int LoadROOTTree(std::string path,
		   CodaRun crun,           // single CODA run
		   bool sort,
		   int verbose,
		   TChain* &C) 
  {
    
    if (crun.runnum != 0) {
       std::cout << "Parsing ROOT files from run " << crun.runnum << std::endl;
      if (!sort) {
	std::string rfname = Form("%s/*%d*",path.c_str(),crun.runnum);
	if (verbose > 1) std::cout << rfname << std::endl;
	C->Add(rfname.c_str());
      } else {
	std::cout << "Sorting by segments.." << std::endl;
	int aRun, aNumFiles, aStream, rc;
	std::vector<int> md;
	std::vector<pair<int, int>> segB_segE;
	// Looping through unique runs
	aRun = crun.runnum;
	rc = GetROOTFileMetaData(path.c_str(),aRun,md,segB_segE,0);
	if (rc==0) { // non-zero number of segments
	  aStream   = md[0]; //stream no.
	  aNumFiles = md[1]; //total # segments
	  if (verbose > 0) {
	    std::cout << "----" << std::endl;
	    std::cout << Form(" Run %d, Total # of segments %d",aRun,aNumFiles) << std::endl;
	    std::cout << Form(" Beg seg %d-%d, End seg %d-%d",segB_segE[0].first,segB_segE[0].second,
			      segB_segE[aNumFiles-1].first,segB_segE[aNumFiles-1].second) << std::endl;
	    std::cout << "----" << std::endl;
	  }
	  // Looping through segments and adding to tree
	  for (int iseg=0; iseg<aNumFiles; iseg++) {
	    std::string rfname = Form("%s/e1209019_fullreplay_%d_stream%d_seg%d_%d.root",path.c_str(),
				      aRun,aStream,segB_segE[iseg].first,segB_segE[iseg].second);
	    if (verbose > 1) std::cout << rfname << std::endl;
	    C->Add(rfname.c_str());
	  }
	} else {
	  std::cout << std::endl << "--!!--" << std::endl  
		    << "WARNING! [util_pd::LoadROOTTree]: No ROOT file exists for run " << crun.runnum << "!" 
		    << std::endl << "--!!--" << std::endl;
	  return -1;
	}
	// getting ready for next run
	md.clear();
	segB_segE.clear();
      }
    }else {
      throw std::runtime_error("[util_pd::LoadROOTTree] CodaRun object is empty!");
    }
    
    return 0;
  }
  /* ###############################
     ## General Purpose Functions ##  
     ############################### */
  //______________________________________________________________________________
  double GetTotCharge(std::vector<CodaRun> crun) 
  /* Calculates total charge from a vector of CodaRuns */
  {
    double totcharge=0.;
    for (auto & run : crun)
      totcharge += run.charge;

    return totcharge;
  }
}
