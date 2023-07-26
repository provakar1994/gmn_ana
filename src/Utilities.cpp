#include "../include/Utilities.h"

namespace util_pd {

  /* #########################################
     ##                General              ##  
     ######################################### */
  TCanvas *TC(std::string name,   // name of the canvas
 	      int rdiv,           // # divisions in row
	      int cdiv)           // # divisions in column 
  /* returns a customized canvas*/
  {
    int w = 1000; int h = 800;
    if (rdiv<1 || cdiv<1) {rdiv = 1; cdiv = 1;}
    // if (rdiv>1 && cdiv==1) {h = 500*rdiv; w = 500;}
    // if (rdiv==1 && cdiv>1) {h = 500; w = 500*cdiv;}
    // if (rdiv>1 && cdiv>1) {h = 1400; w = 1400;}
    if (rdiv>1 || cdiv>1) {h = 450*rdiv; w = 450*cdiv;}
    if (rdiv>4 || cdiv>4) {h = 200*rdiv; w = 200*cdiv;}
    TCanvas *c = new TCanvas(name.c_str(),name.c_str(),w,h);
    c->Divide(cdiv,rdiv);
    return c;
  }

  std::string getDate()
  /* returns today's date */
  {
    time_t now = time(0);
    tm ltm = *localtime(&now);
    std::string yyyy = std::to_string(1900 + ltm.tm_year);
    std::string mm = std::to_string(1 + ltm.tm_mon);
    std::string dd = std::to_string(ltm.tm_mday);
    return mm + '/' + dd + '/' + yyyy;
  }

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
  TH2F *TH2FHCALface_xy_data(std::string name, double sbs_kick) {
    // returns TH2F for HCAL face (x,y) [Data]
    double y_min = expconst::yHCAL_r_DB - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB + expconst::hcalblk_h/2.;
    std::string ylabel = sbs_kick==0 ? "x_{HCAL}^{exp} (m)" : "x_{HCAL}^{exp}-" + std::to_string(sbs_kick) + " (m)";
    TH2F *h = new TH2F(name.c_str(),Form(";y_{HCAL}^{exp} (m);%s",ylabel.c_str()),
		       expconst::hcalcol, y_min, y_max,
		       expconst::hcalrow, x_min, x_max);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_simu(std::string name, double sbs_kick) {
    // returns TH2F for HCAL face (x,y) [Simu]
    double y_min = expconst::yHCAL_r_DB_MC - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB_MC + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB_MC - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB_MC + expconst::hcalblk_h/2.;
    std::string ylabel = sbs_kick==0 ? "x_{HCAL}^{exp} (m)" : "x_{HCAL}^{exp}-" + std::to_string(sbs_kick) + " (m)";
    TH2F *h = new TH2F(name.c_str(),Form(";y_{HCAL}^{exp} (m);%s",ylabel.c_str()),
		       expconst::hcalcol, y_min, y_max,
		       expconst::hcalrow, x_min, x_max);
    return h;
  }
  //_____________________________________
  TH2F *TH2FdxdyHCAL(std::string name) {
    // returns TH2F for dxdyHCAL
    TH2F *h = new TH2F(name.c_str(), ";y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);x_{HCAL}^{obs} - x_{HCAL}^{exp} (m)",
		       250, -1.25, 1.25, 250, -3.5, 2);
    return h;
  }
  //_____________________________________
  void DrawArea(std::vector<double> dimensions, 
		int lcolor=2, int lwidth=4, int lstyle=9) {
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
    else std::cerr << "[Utilities::TH1FhQ2] Enter valid SBS config!!" << std::endl;
    TH1F *h = new TH1F(name.c_str(), "Q^{2} Distribution (GeV^{2})", 
		       nbin, hmin, hmax);
    return h;
  }

  /* ##################################################
     ##   Functions to read CSV file with run info   ##  
     ################################################## */
  void ReadRunList(std::string runsheet_dir,  // Dir. path containing CSV files with run info
		   int &nruns,                // No. of runs to analyze
		   int sbsconf,               // SBS configuration
		   std::string target,        // target type
		   int replay_pass,           // replay pass
		   int verbose,               // verbosity
		   std::vector<CodaRun> &crun)     // Output: Vector of CodaRun structs
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/good_runList_GMn_nTPE_"; 
    std::string mid = "_pass_";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    std::ifstream run_data; run_data.open(run_spreadsheet);
    std::string readline;
    if(run_data.is_open()){
      std::cout << "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;
      std::string skip_header; getline(run_data, skip_header);  // skipping column header
      crun.clear();
      while(getline(run_data,readline)){                   // reading each line
	std::istringstream tokenStream(readline);
	std::string token;
	char delimiter = ',';
	std::vector<std::string> temp;
	while(getline(tokenStream,token,delimiter)){       // reading each element of a line
	  std::string temptoken=token;
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
		   std::vector<CodaRun> &crun)     // Output: Vector of CodaRun objects
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/good_runList_GMn_nTPE_"; 
    std::string mid = "_pass_";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // convert magnet field values from % to A
    sbsmag *= int(expconst::sbsmaxcurr/100);    // 100% SBS magnet current = 2100 A

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    std::ifstream run_data; run_data.open(run_spreadsheet);
    std::string readline;
    if(run_data.is_open()){
      std::cout << "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;
      std::string skip_header; getline(run_data, skip_header); // skipping column header
      while(getline(run_data,readline)){                  // reading each line
	std::istringstream tokenStream(readline);
	std::string token;
	char delimiter = ',';
	std::vector<std::string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  std::string temptoken=token;
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
		   std::vector<CodaRun> &crun)     // Output: Vector of CodaRun structs
  {
    // Define the name of the relevant run spreadsheet
    std::string fst = "/good_runList_GMn_nTPE_"; 
    std::string mid = "_pass_";
    std::string lst = ".csv";
    if (replay_pass < 2) replay_pass = 1; // single spreadsheet exists for pass 0 & 1
    std::string run_spreadsheet = runsheet_dir + fst + target + mid + std::to_string(replay_pass) + lst;

    // convert magnet field values from % to A
    sbsmag *= int(expconst::sbsmaxcurr/100);  // 100% SBS magnet current = 2100 A
    bbmag *= int(expconst::bbmaxcurr/100);    // 100% BB magnet current = 750 A

    // Reading the spreadsheet
    if (nruns < 0) nruns = 1e6;         // replay all runs if nruns < 0
    std::ifstream run_data; run_data.open(run_spreadsheet);
    std::string readline;
    if(run_data.is_open()){
      std::cout << "Reading run info from: "<< run_spreadsheet 
		<< std::endl << std::endl;     
      std::string skip_header; getline(run_data, skip_header); // skipping column header
      while(getline(run_data,readline)){                  // reading each line
	std::istringstream tokenStream(readline);
	std::string token;
	char delimiter = ',';
	std::vector<std::string> temp;
	while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	  std::string temptoken=token;
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
     ## Functions to load ROOT files by sorted segment ##  
     #################################################### */
  //______________________________________________________________________________
  void ListDirectory(const char *path,
		     std::vector<std::string> &list){
    /* Acknowledgement: Based on a function written by David Flay. */
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
  int SplitString(const char delim, 
		  const std::string myStr, 
		  std::vector<std::string> &out){
    /* Acknowledgement: Based on a function written by David Flay. */
    // split a string by a delimiter
    std::stringstream ss(myStr);
    while( ss.good() ){
      std::string substr;
      std::getline(ss, substr, delim);
      if (!substr.empty()) out.push_back(substr);
    }
    if (out.empty()) std::cerr << "WARNING! No substrings found!\n";
    return 0;
  }
  //______________________________________________________________________________
  bool sortbyval(const std::pair<int, int> &a, const std::pair<int, int> &b) {
    return (a.first < b.first);
  } 
  //______________________________________________________________________________
  int GetROOTFileMetaData(const char *rfDirPath, int run,
			  std::vector<int> &data, 
			  std::vector<std::pair<int, int>> &segB_segE,
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
	  segB_segE.push_back(std::make_pair(bseg,eseg));
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
	std::vector<std::pair<int, int>> segB_segE;
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
	std::vector<std::pair<int, int>> segB_segE;
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

  /* ##################################################
     ##   Function to read MC replay summary files   ##  
     ################################################## */
  //______________________________________________________________________________
  void ReadSimuJobSummary(std::string logfile_dir,   // Dir. path containing MC summary files
			  std::string prefix,        // prefix to standard filebase (Special case handling)
			  int sbsconf,               // SBS configuration
			  int sbsmag,                // SBS magnet current (in %)
			  std::string generator,     // simc / g4sbs
			  std::string target,        // target type
			  int &njobs,                // # jobs to analyze per process
			  int verbose,               // verbosity
			  std::vector<SimuJob> &sjobs)    // Output: Vector of SimuJob objects
  /* Reads simulation job specifics from summary files and loads the values to SimuJob objects.
     This function has the standard naming conventions of output simulation and summary files 
     hard coded. Please make sure the files to analyze have names compatible with the standard
     naming convention for successful execution.
     NOTE: The presence of a summary file is must for this process to work. In case of SIMC 
     D(ee'N) process, the ROOT file directory must have two summary files presents. One for 
     the deep process and one for the deen process. Each file must have at least one entry.
     ---------
     Standard naming conventions:
     1. Filebase: prefix_sbs<sbsconf>_sbs<sbsmag>p_<generator>
     2. MC job summary file: <filebase>_summary.csv
     3. Digitized ROOT file: <filebase>_<process>_job_<jobid>.root
     4. Replayed digitized ROOT file: replayed_<filebase>_<process>_job_<jobid>.root
  */
  {
    TString filebase = Form("sbs%d_sbs%dp_%s",sbsconf,sbsmag,generator.c_str());
    if (!prefix.empty()) filebase = (TString)prefix + "_" + filebase;

    // Define the name of the summary file and corresponding process based on generator and target
    std::vector<TString> simu_logfile, process;
    if (target.compare("LH2") == 0) {
      TString temp = Form("%s_heep_summary.csv",filebase.Data());
      simu_logfile.push_back(temp); process.push_back("heep");
    } 
    else if (target.compare("LD2") == 0) {
      TString temp = "";
      if (generator.compare("simc") == 0) {
	temp = Form("%s_deep_summary.csv",filebase.Data());
	simu_logfile.push_back(temp); process.push_back("deep");
	temp = Form("%s_deen_summary.csv",filebase.Data());
	simu_logfile.push_back(temp); process.push_back("deen");
      } else {
	temp = Form("%s_deeN_summary.csv",filebase.Data());
	simu_logfile.push_back(temp); process.push_back("deeN");
      }
    }
    else
      throw std::invalid_argument("[util_pd::ReadSimuJobSummary] Given target type is invalid!");

    // Reading the summary spreadsheet
    if (njobs < 0) njobs = 1e4;         // replay all runs if njobs < 0
    // looping through list of summary files
    for (int ifile=0; ifile<simu_logfile.size(); ifile++) {
      int inisize = sjobs.size();
      TString logfile_temp = Form("%s/%s",logfile_dir.c_str(),simu_logfile[ifile].Data());
      std::ifstream simu_log; simu_log.open(logfile_temp);
      std::string readline;
      if(simu_log.is_open()){
	std::cout << "Reading summary file: " << logfile_temp << std::endl;
	std::string skip_header; getline(simu_log,skip_header); // skipping column header
	while(getline(simu_log,readline)){                  // reading each line
	  std::istringstream tokenStream(readline);
	  std::string token;
	  char delimiter = ',';
	  std::vector<std::string> temp, data;
	  while(getline(tokenStream,token,delimiter)){      // reading each element of a line
	    std::string temptoken=token;
	    temp.push_back(temptoken);
	  }
	  // add relevant info to SimuJob objects
	  if (sjobs.size() >= njobs) {njobs += njobs; break;}
	  SimuJob temp_sj;
	  int jobid = stoi(temp[0]);
	  TString sfname = Form("%s/%s_%s_job_%d.root",logfile_dir.c_str(),
				filebase.Data(),process[ifile].Data(),jobid); 
	  TString rfname = Form("%s/replayed_%s_%s_job_%d.root",logfile_dir.c_str(),
				filebase.Data(),process[ifile].Data(),jobid); 
	  // handling the sicrepancy in summary file by generator
	  data.push_back(sfname.Data());
	  data.push_back(rfname.Data());
	  data.push_back(generator);
	  // if (generator.compare("simc") == 0) {
	  //   data.push_back(temp[1]); //ngenreq
	  //   data.push_back(temp[2]); //nthrown
	  //   data.push_back(temp[5]); //genvol(MeV*sr2)
	  //   data.push_back(temp[4]); //lumi(ub^-1)
	  //   data.push_back(temp[3]); //charge(mC)
	  // }else	if (generator.compare("g4sbs") == 0) {
	  //   data.push_back(temp[1]); //ngenreq
	  //   data.push_back(temp[2]); //nthrown
	  //   data.push_back(temp[5]); //genvol(sr)
	  //   data.push_back(temp[6]); //lumi(Hz/cm2)
	  //   data.push_back(temp[3]); //ebeam(GeV) 
	  //   data.push_back(temp[4]); //ibeam(A)
	  // }
	  data.push_back(temp[1]); //ngenreq
	  data.push_back(temp[2]); //nthrown
	  data.push_back(temp[3]); //genvol(MeV*sr2)
	  data.push_back(temp[4]); //lumi(ub^-1)
	  data.push_back(temp[5]); //ebeam(GeV)
	  if (generator.compare("simc") == 0) data.push_back(temp[6]);  //charge(mC)
	  if (generator.compare("g4sbs") == 0) data.push_back(temp[6]); //ibeam(muA)

	  temp_sj.SetDataSimuJob(data);
	  sjobs.push_back(temp_sj);

	  temp.clear();
	  data.clear();
	}
	if (verbose > 0) {
	  std::cout << "First job info:" << std::endl << sjobs[inisize];
	  std::cout << "Last job info:" << std::endl << sjobs[sjobs.size()-1];
	}
      }else
	throw std::runtime_error("[util_pd::ReadSimuJobSummary] Summary file doesn't exist");
      simu_log.close();
    }
    std::cout << std::endl;
    // let's update njobs with total no. of runs to analyze
    njobs = sjobs.size();
  }
  //______________________________________________________________________________
  void LoadSimuROOTTree(std::vector<SimuJob> sjobs,  // SimuJob objects with run info
			int verbose,                 // verbosity
			TChain* &C)                  // Output: TChain with data
  {
    if (!sjobs.empty()) {
      const int njobs = sjobs.size();
      std::cout << "Parsing ROOT files from " << njobs << " jobs.." << std::endl;
      for (int ijob=0; ijob<sjobs.size(); ijob++) {
	if (verbose > 1) std::cout << sjobs[ijob].rfname << std::endl;
	C->Add(sjobs[ijob].rfname.c_str());
      }
      if (C->GetEntries()==0) 
	throw std::runtime_error("[util_pd::LoadSimuROOTTree] Empty ROOT files Or, they don't exist!");
    }else {
      throw std::runtime_error("[util_pd::LoadSimuROOTTree] Simu job list is empty!");
    }
  }

  /* ###############################
     ## General Purpose Functions ##  
     ############################### */
  //______________________________________________________________________________
  double GetTotCharge(std::vector<CodaRun> cruns) 
  /* Calculates total charge from a vector of CodaRuns */
  {
    double totcharge=0.;
    for (auto & run : cruns)
      totcharge += run.charge;
    return totcharge;
  }
  //______________________________________________________________________________
  double GetTotCharge(std::vector<SimuJob> sjobs) 
  /* Calculates total charge from a vector of SimuJobs */
  {
    double totcharge=0.;
    for (auto & job : sjobs)
      totcharge += job.charge;
    return totcharge;
  }
  //______________________________________________________________________________
  double GetTotNtries(std::vector<SimuJob> sjobs) 
  /* Calculates total # tries from a vector of SimuJobs */
  {
    double totntries=0.;
    for (auto & job : sjobs)
      totntries += job.ntried;
    return totntries;
  }
  //______________________________________________________________________________
  void GetTotNtriesnCh(std::vector<SimuJob> sjobs, // Input: List of SimuJob objects
		       std::vector<double> &data)  // Output: data[0]=>Tot. Ch., data[1]=>Tot. ntries
  /* Calc. # tries & tot. Ch. from SimuJob objects */
  {
    double totntries=0.,totcharge=0.;
    for (auto & job : sjobs) {
      totntries += job.ntried;
      totcharge += job.charge;
    }
    data = {totntries,totcharge};
  }
  //______________________________________________________________________________
  void GetElossInTgt(std::string const target, // Target type
		    int const rnum,            // Run number
		    int const sbsconf,         // SBS configuration
		    double const vz,           // m, vertex z co-ordinate
		    double const etheta,       // rad, scattering angle
		    int const verbose,         // verbosity
		    std::vector<double> &data) // Output: data[0]=>before scattering, data[1]=>after scattering,
  /* Calculates energy loss in the target before and after scattering */
  {
    double tgt_rho, tgt_dEdx_bs, tgt_dEdx_as, tgt_uwinthick, tgt_celldiam, tgt_cellthick; 
    double Al_dEdx_bs = expconst::GetdEdxCollAl(sbsconf,1), Al_dEdx_as = expconst::GetdEdxCollAl(sbsconf,0);
    double PE_dEdx_bs = expconst::GetdEdxCollPE(sbsconf,1), PE_dEdx_as = expconst::GetdEdxCollPE(sbsconf,0); 
    if (target.compare("LH2")==0) {
      tgt_rho = expconst::lh2_TgtRho;
      //tgt_dEdx = expconst::lh2_dEdx;
      tgt_dEdx_bs = expconst::GetdEdxCollH(sbsconf,1);
      tgt_dEdx_as = expconst::GetdEdxCollH(sbsconf,0);
      tgt_uwinthick = expconst::lh2_uWinThick;
      tgt_celldiam = expconst::lh2_CellDiam;
      tgt_cellthick = expconst::lh2_CellThick;
    }else if (target.compare("LD2")==0) {
      tgt_rho = expconst::ld2_TgtRho;
      //tgt_dEdx = expconst::ld2_dEdx;
      tgt_dEdx_bs = expconst::GetdEdxCollH(sbsconf,1);
      tgt_dEdx_as = expconst::GetdEdxCollH(sbsconf,0);
      tgt_uwinthick = expconst::ld2_uWinThick;
      tgt_celldiam = expconst::ld2_CellDiam;
      tgt_cellthick = expconst::ld2_CellThick;
    }else
      throw std::invalid_argument("[util_pd::GetElossInTgt] Given target type is invalid!");
      
    // Scattering chamber shielding 
    /* NOTE: 
       1. Additional shielding beside the scattering chamber was installed during SBS11 to reduce 
          background rates in the front tracker (I think).
       2. Initially, Bogdan just put a 10mm plastic shielding leaning on the scattering chamber!
          Bogdan did that on 12/04/2021 at 18:20. https://logbooks.jlab.org/entry/3956556
	  First run after the installation: 12556
	  - Andrew suggested to assume the material to be acrylic or polyethylene. I am going with
	    polyethylene (PE). 
       3. Then on 12/08/2021, Andrew (tech) replaced the plastic shielding with a permanent 1/8th
          inch thick Al shield. https://logbooks.jlab.org/entry/3958901 
	  First run after the installation: 12675
    */
    bool PE_shield_in = (rnum<12556 || rnum>=12675) ? false : true;
    bool Al_shield_in = rnum<12675 ? false : true;

    double Eloss_bf = (vz+0.08)*100.0*tgt_rho*tgt_dEdx_bs + tgt_uwinthick*expconst::Al_Rho*Al_dEdx_bs;
    double Eloss_af = tgt_celldiam/2.0/sin(etheta)*tgt_rho*tgt_dEdx_as + tgt_cellthick/sin(etheta)*expconst::Al_Rho*Al_dEdx_as;
    if (PE_shield_in) Eloss_af += expconst::PE_ShieldThick*expconst::PE_Rho*PE_dEdx_as;
    if (Al_shield_in) Eloss_af += expconst::Al_ShieldThick*expconst::Al_Rho*Al_dEdx_as;
    data = {Eloss_bf,Eloss_af};

    if (verbose>1) {
      std::cout << Form("\netheta = %.1f deg, effective vz = %.3f m \n",etheta*TMath::RadToDeg(),vz+0.08);
      std::cout << Form("Plastic shield in place = %d, Al shield in place = %d \n",PE_shield_in,Al_shield_in);
      std::cout << Form("Mean energy loss (MeV) %.1f (before), %.1f (after) \n",Eloss_bf*1E3,Eloss_af*1E3);
    }
  }
}
