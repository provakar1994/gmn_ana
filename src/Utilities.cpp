#include "../include/Utilities.h"

#include <filesystem>

namespace util_pd {

  /* #########################################
     ##                General              ##  
     ######################################### */
  //_____________________________________
  TFile *ReadRootFile(char const * filename)
  /* Open and return a ROOT file */
  {
    // Open the ROOT file
    TFile *file = TFile::Open(filename, "READ");
    // Check if the file is open and the histogram exists
    if (!file || file->IsZombie()){ //|| !file->GetListOfKeys()->Contains(histname)) {
      std::cerr << "Error: Failed to open the file or histogram not found!" << std::endl;
      throw;
    }
    return file;
  }
    
  //_____________________________________
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

  // //_____________________________________
  // void SetAxTitles(TH1* obj, TString const & xtitle, TString const & ytitle) {
  //   obj->GetXaxis()->SetTitle(xtitle);
  //   obj->GetYaxis()->SetTitle(ytitle);
  // }
  // //_____________________________________
  // void SetAxTitles(TH2* obj, TString const & xtitle, TString const & ytitle) {
  //   obj->GetXaxis()->SetTitle(xtitle);
  //   obj->GetYaxis()->SetTitle(ytitle);
  // }  
  // //_____________________________________
  // void SetAxTitles(TGraph* obj, TString const & xtitle, TString const & ytitle) {
  //   obj->GetXaxis()->SetTitle(xtitle);
  //   obj->GetYaxis()->SetTitle(ytitle);
  // }
  
  //_____________________________________
  int GetSigDigit(double value, int position)
  /* returns the desired significant digit */
  {
    // Shift the decimal point to the right by (position) places
    double scaledValue = value * std::pow(10, position);

    // Take the integer part of the scaled value
    int integerPart = static_cast<int>(scaledValue);

    // Extract the digit in the desired position
    int significantDigit = integerPart % 10;

    // Return the significant digit
    return significantDigit;
  }
  
  //_____________________________________
  std::vector<TPad*> GetPadsForPullPlot(TCanvas *c1) 
  /* splits a given canvas into two pads suitable for pull plots */
  // {
  //   c1->cd();
  //   TPad *p1 = new TPad("p1","p1",0,0.25,1,1);    //for fit
  //   TPad *p2 = new TPad("p2","p2",0,0.02,1,0.25); //for residual
  //   p1->SetBottomMargin(0.00001); p1->SetBorderMode(0);
  //   p1->SetTickx(); p1->SetTicky();
  //   p1->SetGridx(); p1->Draw();
  //   p2->SetTopMargin(0.00001); p2->SetBottomMargin(0.2);
  //   p2->SetBorderMode(0);
  //   p2->SetTickx(); p2->SetTicky(); 
  //   p2->SetGridx(); p2->Draw();
  //   std::vector<TPad*> pads{p1,p2};
  //   return pads;
  // }
{
    c1->cd();
    TPad *p1 = new TPad("p1", "p1", 0, 0.3, 1, 1);    // Top pad for the fit
    TPad *p2 = new TPad("p2", "p2", 0, 0.05, 1, 0.3); // Bottom pad for the residual/pull plot

    // Configure p1 (Top pad)
    p1->SetBottomMargin(0.019);  // Adjust this margin to provide space for the lower pad
    p1->SetRightMargin(0.05);
    p1->SetBorderMode(0);
    p1->SetTickx(); p1->SetTicky();
    p1->SetGridx(); p1->Draw();

    // Configure p2 (Bottom pad)
    p2->SetTopMargin(0.015);      // Keep this small to give space to p1
    p2->SetBottomMargin(0.3);    // Increase to provide enough space for the x-axis title
    p2->SetRightMargin(0.05);
    p2->SetBorderMode(0);
    p2->SetTickx(); p2->SetTicky();
    p2->SetGridx(); p2->Draw();

    // Return the pads
    std::vector<TPad*> pads{p1, p2};
    return pads;
}    

  //_____________________________________
  void DrawZeroLine(TPad *p1, double xmin, double xmax)
  /* Draws a horizontal line at y=0 */
  {
    p1->cd();
    TLine* zeroLine = new TLine(xmin,0,xmax,0);
    zeroLine->SetLineColor(kGray+2);
    zeroLine->SetLineWidth(2);
    zeroLine->SetLineStyle(9);
    zeroLine->Draw("same"); 
  }

  //_____________________________________
  double GetxNDC(double x)
  /* Returns normalized coordinate for a given x value */
  {
    gPad->Update();
    return (x - gPad->GetX1())/(gPad->GetX2()-gPad->GetX1());
  }

  //_____________________________________
  double GetyNDC(double y)
  /* Returns normalized coordinate for a given y value */
  {
    gPad->Update();
    return (y - gPad->GetY1())/(gPad->GetY2()-gPad->GetY1());
  }  

  //_____________________________________
  double UnfoldyNDC(double yNDC) 
  /* Returns y value for a given y coordinate in NDC */
  {
    gPad->Update();
    return yNDC*(gPad->GetY2()-gPad->GetY1()) + gPad->GetY1();
  }

  //_____________________________________
  void PlotCutRegion(double xmin, double xmax)
  /* Plots cut region for given x range */
  {
    TBox *cutRegion = new TBox(xmin,UnfoldyNDC(0.1),xmax,UnfoldyNDC(0.9));
    cutRegion->SetFillColorAlpha(kRed, 0.3);
    cutRegion->Draw();
  }
  
  //______________________________________________________________________________
  void PlotFiduCut(int pass, std::vector<double> hcal_AR, std::vector<double> hcal_SM) {
    std::vector<double> hcal_boundary = cut::hcal_active_area_data(0,0,pass); 
    util_pd::DrawArea(hcal_boundary,kGreen+2,2,1);
    util_pd::DrawArea(hcal_AR,2,4,9);
    util_pd::DrawArea(hcal_SM,4,4,9);
  }
  
  //_____________________________________
  void TF1toTH1F(TF1* const func, // TF1 object to mimic
		 TH1F* &hist)     // TH1F object to modify
  /* creates a TH1F object from a given TF1 object */
  {
    for (int i = 1; i <= hist->GetNbinsX(); ++i) {
      double binCenter = hist->GetBinCenter(i);
      double functionValue = func->Eval(binCenter);
      hist->SetBinContent(i, functionValue);
    }
  }

  //_____________________________________
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

  //_____________________________________  
  double GetXSpreadOfTH1(const TH1* hist, double minCounts, int verbose) {
    /* returns the spread of TH1 in X based on a threshold of counts */
    if (!hist || hist->GetNbinsX() == 0)
      return 0.0;

    double xmin = -1, xmax = -1;
    int nbins = hist->GetNbinsX();

    for (int i = 1; i <= nbins; ++i) {
      double content = hist->GetBinContent(i);
      if (content >= minCounts) {
	double xcenter = hist->GetBinCenter(i);
	if (xmin < 0) xmin = xcenter;  // first qualifying bin
	xmax = xcenter;               // always update to latest
      }
    }

    if (verbose>0) {
      std::cout << "xmin: " << xmin << ", xmax: " << xmax << "\n";
    }
    
    if (xmin < 0 || xmax < 0) return 0.0; // no bins met the threshold
    return xmax - xmin;
  }  

  /* #################################################
     ##                HCAL Related                 ##  
     ################################################# */
  //_____________________________________
  TH2F *TH2FHCALface_rc(std::string hname) {
    // returns TH2F for HCAL face (row,col)
    /* NOTE: HCAL block id (ibblk) starts from 1 and goes up to 288 but both
       HCAL row (rowblk) and column starts from 0 and goes up to 23 and 
       11, respectively. Extremely annoying! */
    // TH2F *h = new TH2F(hname.c_str(), ";HCAL columns;HCAL rows",
    // 		       expconst::hcalcol, 0, expconst::hcalcol,
    // 		       expconst::hcalrow, 0, expconst::hcalrow);
    TH2F *h = new TH2F(hname.c_str(), ";HCAL columns;HCAL rows",
		       200, 0, expconst::hcalcol,
		       200, 0, expconst::hcalrow);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_data(std::string hname, double sbs_kick, int rpass) {    
    // returns TH2F for HCAL face (x,y) [Data]
    // sbs_kick = -99 assumes p_deflection has been used
    // block positions from DB
    double xHCAL_t_DB = rpass<2 ? expconst::xHCAL_t_DB_p1 : expconst::xHCAL_t_DB;   //m, center of top row blocks (from DB)
    double xHCAL_b_DB = rpass<2 ? expconst::xHCAL_b_DB_p1 : expconst::xHCAL_b_DB;   //m, center of bottom row blocks (from DB)
    double yHCAL_r_DB = rpass<2 ? expconst::yHCAL_r_DB_p1 : expconst::yHCAL_r_DB;   //m, center of right most blocks (from DB)
    double yHCAL_l_DB = rpass<2 ? expconst::yHCAL_l_DB_p1 : expconst::yHCAL_l_DB;   //m, center of left most blocks (from DB)
    // cut region
    double y_min = yHCAL_r_DB - expconst::hcalblk_w/2.;
    double y_max = yHCAL_l_DB + expconst::hcalblk_w/2.;
    double x_min = xHCAL_t_DB - expconst::hcalblk_h/2.;
    double x_max = xHCAL_b_DB + expconst::hcalblk_h/2.;
    std::string ylabel = sbs_kick==0 ? "x_{HCAL}^{exp} (m)" : "x_{HCAL}^{exp}-" + std::to_string(sbs_kick) + " (m)";
    if (sbs_kick==-99) ylabel = "x_{HCAL}^{exp} - #deltax_{SBS} (m)"; 
    // TH2F *h = new TH2F(hname.c_str(),Form(";y_{HCAL}^{exp} (m);%s",ylabel.c_str()),
    // 		       expconst::hcalcol, y_min, y_max,
    // 		       expconst::hcalrow, x_min, x_max);
    TH2F *h = new TH2F(hname.c_str(),Form(";y_{HCAL}^{exp} (m);%s",ylabel.c_str()),
		       200, -1.25, 1.25,
		       200, -3.25, 2.25);
    return h;
  }
  //_____________________________________
  TH2F *TH2FHCALface_xy_simu(std::string hname, double sbs_kick) {
    // returns TH2F for HCAL face (x,y) [Simu]
    double y_min = expconst::yHCAL_r_DB_MC - expconst::hcalblk_w/2.;
    double y_max = expconst::yHCAL_l_DB_MC + expconst::hcalblk_w/2.;
    double x_min = expconst::xHCAL_t_DB_MC - expconst::hcalblk_h/2.;
    double x_max = expconst::xHCAL_b_DB_MC + expconst::hcalblk_h/2.;
    std::string ylabel = sbs_kick==0 ? "x_{HCAL}^{exp} (m)" : "x_{HCAL}^{exp}-" + std::to_string(sbs_kick) + " (m)";
    if (sbs_kick==-99) ylabel = "x_{HCAL}^{exp} - #deltax_{SBS} (m)"; 
    // TH2F *h = new TH2F(hname.c_str(),Form(";y_{HCAL}^{exp} (m);%s",ylabel.c_str()),
    // 		       expconst::hcalcol, y_min, y_max,
    // 		       expconst::hcalrow, x_min, x_max);
    TH2F *h = new TH2F(hname.c_str(),Form(";y_{HCAL}^{exp} (m);%s",ylabel.c_str()),
		       200, -1.25, 1.25,
		       200, -3.25, 2.25);
    return h;
  }
  //_____________________________________
  TH2F *TH2FdxdyHCAL(std::string hname) {
    // returns TH2F for dxdyHCAL
    TH2F *h = new TH2F(hname.c_str(), ";y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);x_{HCAL}^{obs} - x_{HCAL}^{exp} (m)",
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
  TH1F *TH1FhW(std::string hname) {
    // returns W histogram
    TH1F *h = new TH1F(hname.c_str(), "W Distribution (GeV)", 250,0,2);
    return h;
  }
  TH1F *TH1FhQ2(std::string hname,      // Name of histogram
		int conf) {             // SBS config
    // returns Q2 histogram
    int nbin=0; double hmin=-100, hmax=-100;
    if (conf==4) { nbin=100; hmin=1.; hmax=4.; } 
    else if (conf==14) { nbin=100; hmin=5.; hmax=10.; }
    else if (conf==7) { nbin=120; hmin=6.; hmax=12.; }
    else if (conf==11) { nbin=200; hmin=8.; hmax=18.; }
    else std::cerr << "[Utilities::TH1FhQ2] Enter valid SBS config!!" << std::endl;
    TH1F *h = new TH1F(hname.c_str(), "Q^{2} Distribution (GeV^{2})", 
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
  int ReadLog(const std::string& filename) {
    /* Bare bone function to read Pood generated log files to extract
       accepted BBCal singles trigger. In future, it should be made generalized
       so that it can be used to read other log file entries as well.
    */
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
      std::cerr << "Error: Unable to open file " << filename << std::endl;
      return -999;
    }

    // Read each line of the file
    std::string myStr;
    while (std::getline(file, line)) {
      if (line.find("BLOCK: Physics") != std::string::npos) {
	// Start reading Physics block
	while (std::getline(file, line)) {
	  if (line.find("BBCalSingles") != std::string::npos) {
	    myStr = line; // Return BBCalSingles entry

	    // Tokenizing line to extract passed and called
	    std::vector<std::string> tokens;
	    util_pd::SplitString(' ',myStr,tokens);
	  
	    return stoi(tokens[3]); // 3=> Passed
	  }
	}
      }
    }

    std::cerr << "Error: BBCalSingles entry not found in file " << filename << std::endl;
    return -999;
  }
  //______________________________________________________________________________
  int LoadROOTTree(std::string path,
		   CodaRun &crun,           // single CODA run
		   bool sort,
		   int verbose,
		   bool is_read_log,
		   TChain* &C) 
  /* if (is_read_log) - Reads log file and adds passed BBCALsingles to crun  */
  {
    
    std::string rfpath = path + "/rootfiles";
    std::string lfpath = path + "/logs";

    if (crun.runnum != 0) {
      std::cout << "Parsing ROOT files from run " << crun.runnum << std::endl;
      if (!sort) {
	std::string rfname = Form("/*%d*",rfpath.c_str(),crun.runnum);
	std::string lfname = Form("/*%d*",lfpath.c_str(),crun.runnum);
	if (verbose > 1) std::cout << rfname << std::endl;
	C->Add(rfname.c_str());
      } else {
	std::cout << "Sorting by segments.." << std::endl;
	int aRun, aNumFiles, aStream, rc;
	std::vector<int> md;
	std::vector<std::pair<int, int>> segB_segE;
	// Looping through unique runs
	aRun = crun.runnum;
	rc = GetROOTFileMetaData(rfpath.c_str(),aRun,md,segB_segE,0);
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
	  double BBCalSingles_passed = 0;
	  for (int iseg=0; iseg<aNumFiles; iseg++) {
	    std::string rfname = Form("%s/e1209019_fullreplay_%d_stream%d_seg%d_%d.root",rfpath.c_str(),
				      aRun,aStream,segB_segE[iseg].first,segB_segE[iseg].second);
	    if (verbose > 1) std::cout << rfname << std::endl;
	    C->Add(rfname.c_str());

	    if (is_read_log) {
	      // reading corresponding log files
	      std::string lfname = Form("%s/e1209019_fullreplay_%d_stream%d_seg%d_%d.log",lfpath.c_str(),
					aRun,aStream,segB_segE[iseg].first,segB_segE[iseg].second);
	      int temp = ReadLog(lfname); //system(Form("python3 daq_livetime.py %s",lfname));
	      if (temp<0) std::cerr << Form("[WARNING] %s doesn't exist!! \n",lfname);
	      if (verbose > 1) {
		std::cout << lfname << std::endl;
		std::cout << temp << std::endl;
	      }
	      BBCalSingles_passed += temp;
	    }
	  }
	  // adding total BBCalSingles to crun
	  crun.BBCalSinglesPassed = BBCalSingles_passed;
	  if (verbose > 0) 
	    std::cout << Form("Total passed BBCal singles for run %d: %d\n",crun.runnum,crun.BBCalSinglesPassed);
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
  // //______________________________________________________________________________
  // int LoadROOTTree(std::string path,
  // 		   CodaRun crun,           // single CODA run
  // 		   bool sort,
  // 		   int verbose,
  // 		   TChain* &C) 
  // {
    
  //   if (crun.runnum != 0) {
  //      std::cout << "Parsing ROOT files from run " << crun.runnum << std::endl;
  //     if (!sort) {
  // 	std::string rfname = Form("%s/*%d*",path.c_str(),crun.runnum);
  // 	if (verbose > 1) std::cout << rfname << std::endl;
  // 	C->Add(rfname.c_str());
  //     } else {
  // 	std::cout << "Sorting by segments.." << std::endl;
  // 	int aRun, aNumFiles, aStream, rc;
  // 	std::vector<int> md;
  // 	std::vector<std::pair<int, int>> segB_segE;
  // 	// Looping through unique runs
  // 	aRun = crun.runnum;
  // 	rc = GetROOTFileMetaData(path.c_str(),aRun,md,segB_segE,0);
  // 	if (rc==0) { // non-zero number of segments
  // 	  aStream   = md[0]; //stream no.
  // 	  aNumFiles = md[1]; //total # segments
  // 	  if (verbose > 0) {
  // 	    std::cout << "----" << std::endl;
  // 	    std::cout << Form(" Run %d, Total # of segments %d",aRun,aNumFiles) << std::endl;
  // 	    std::cout << Form(" Beg seg %d-%d, End seg %d-%d",segB_segE[0].first,segB_segE[0].second,
  // 			      segB_segE[aNumFiles-1].first,segB_segE[aNumFiles-1].second) << std::endl;
  // 	    std::cout << "----" << std::endl;
  // 	  }
  // 	  // Looping through segments and adding to tree
  // 	  for (int iseg=0; iseg<aNumFiles; iseg++) {
  // 	    std::string rfname = Form("%s/e1209019_fullreplay_%d_stream%d_seg%d_%d.root",path.c_str(),
  // 				      aRun,aStream,segB_segE[iseg].first,segB_segE[iseg].second);
  // 	    if (verbose > 1) std::cout << rfname << std::endl;
  // 	    C->Add(rfname.c_str());
  // 	  }
  // 	} else {
  // 	  std::cout << std::endl << "--!!--" << std::endl  
  // 		    << "WARNING! [util_pd::LoadROOTTree]: No ROOT file exists for run " << crun.runnum << "!" 
  // 		    << std::endl << "--!!--" << std::endl;
  // 	  return -1;
  // 	}
  // 	// getting ready for next run
  // 	md.clear();
  // 	segB_segE.clear();
  //     }
  //   }else {
  //     throw std::runtime_error("[util_pd::LoadROOTTree] CodaRun object is empty!");
  //   }
    
  //   return 0;
  // }
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


  /* ##################################################
     ##   Function to read MC replay summary files   ##  
     ################################################## */
  //______________________________________________________________________________
  void ReadSimuJobSummary(std::string logfile_dir,   // Dir. path containing MC summary files
			  std::vector<std::string> const &prefixes, // prefix to standard filebase, Convention: [prefix_deep,prefix_deen] 
			  int sbsconf,               // SBS configuration
			  int sbsmag,                // SBS magnet current (in %)
			  std::string generator,     // simc / g4sbs
			  std::string process,       // reaction process being simulated
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
    std::string prefix, prefix_deen;
    bool prefixesGTone = false; 
    if (!prefixes.empty()) {
      if (prefixes.size()>1 && generator.compare("simc")==0 && process.compare("heep")!=0) { //in case deep and deen has different prefixces
	prefix = prefixes[0];
	prefix_deen = prefixes[1];
	prefixesGTone = true;
      } else {
	prefix = prefixes[0];
      }
    }
    
    TString filebase_std = Form("sbs%d_sbs%dp_%s",sbsconf,sbsmag,generator.c_str());
    TString filebase, filebase_deen;
    if (!prefix.empty()) filebase = (TString)prefix + "_" + filebase_std;
    if (prefixesGTone && !prefix_deen.empty())
      filebase_deen = (TString)prefix_deen + "_" + filebase_std;
    
    // Define the name of the summary file
    std::vector<TString> simu_logfile, processes;
    if (generator.compare("simc")==0 && process.compare("deeN")==0) {
      TString temp = Form("simcout/%s_deep_summary.csv",filebase.Data());
      simu_logfile.push_back(temp); processes.push_back("deep");
      temp = !prefixesGTone ? Form("simcout/%s_deen_summary.csv",filebase.Data())
	: Form("simcout/%s_deen_summary.csv",filebase_deen.Data());
      simu_logfile.push_back(temp); processes.push_back("deen");
    } else if (generator.compare("simc")==0 && process.compare("deep")==0) {
      TString temp = Form("simcout/%s_deep_summary.csv",filebase.Data());
      simu_logfile.push_back(temp); processes.push_back("deep");
    } else if (generator.compare("simc")==0 && process.compare("deen")==0) {
      TString temp = !prefixesGTone ? Form("simcout/%s_deen_summary.csv",filebase.Data())
	: Form("simcout/%s_deen_summary.csv",filebase_deen.Data());
      simu_logfile.push_back(temp); processes.push_back("deen");
    } else if (generator.compare("simc")==0 && (process.compare("heep")==0 || process.compare("neen")==0)) {
      TString temp = Form("simcout/%s_%s_summary.csv",filebase.Data(),process.c_str());
      simu_logfile.push_back(temp); processes.push_back(process);
    } else {
      TString temp = Form("%s_%s_summary.csv",filebase.Data(),process.c_str());
      simu_logfile.push_back(temp); processes.push_back(process);
    }

    // Reading the summary spreadsheet
    if (njobs < 0) njobs = 1e4;         // replay all runs if njobs < 0
    // looping through list of summary files
    for (int ifile=0; ifile<simu_logfile.size(); ifile++) {
      int inisize = sjobs.size();
      TString logfile_temp = Form("%s/%s",logfile_dir.c_str(),simu_logfile[ifile].Data());
      std::ifstream simu_log; simu_log.open(logfile_temp);
      std::string readline;
      //bool is_simc = false; // true if at least one SIMC deeN process summary file is present.      
      if(simu_log.is_open()){
	//if (generator.compare("simc")==0 && process.compare("deeN")==0) is_simc = true;
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
	  TString fbase = prefixesGTone&&processes[ifile].CompareTo("deen")==0 ? filebase_deen : filebase;
	  TString sfname = Form("%s/%s_%s_job_%d.root",logfile_dir.c_str(),
				fbase.Data(),processes[ifile].Data(),jobid); 
	  TString rfname = Form("%s/replayed_%s_%s_job_%d.root",logfile_dir.c_str(),
				fbase.Data(),processes[ifile].Data(),jobid);
	  
	  // handling the sicrepancy in summary file by generator
	  data.push_back(sfname.Data());
	  data.push_back(rfname.Data());
	  data.push_back(generator);
	  data.push_back(processes[ifile].Data());
	  data.push_back(temp[1]); //ngenreq
	  data.push_back(temp[2]); //nthrown
	  data.push_back(temp[3]); //genvol
	  data.push_back(temp[4]); //lumi
	  data.push_back(temp[5]); //ebeam(GeV)

	  // let's check if the replayed ROOT file exists or not
	  std::string is_rfile = "0";
	  if (std::filesystem::exists(rfname.Data())) is_rfile = "1";
	  // TFile* file_temp = TFile::Open(rfname.Data(), "READ");
	  // if (file_temp && !file_temp->IsZombie()) is_rfile = "1";
	  data.push_back(is_rfile);

	  if (generator.compare("g4sbs") == 0) data.push_back(temp[6]); //ibeam(muA)
	  if (generator.compare("simc") == 0) {
	    data.push_back(temp[6]);  //charge(mC)
	    if (temp.size()>8) {
	      data.push_back(temp[8]); //using_RS
	      data.push_back(temp[9]); //max_wt_RS
	    }
	  }
	  
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
	throw std::runtime_error(Form("[util_pd::ReadSimuJobSummary] Summary file, %s, doesn't exist",logfile_temp.Data()));
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
      if (job.is_rfile) totcharge += job.charge;
    return totcharge;
  }
  //______________________________________________________________________________
  double GetTotNtries(std::vector<SimuJob> sjobs) 
  /* Calculates total # tries from a vector of SimuJobs */
  {
    double totntries=0.;
    for (auto & job : sjobs)
      if (job.is_rfile) totntries += job.ntried;
    return totntries;
  }
  //______________________________________________________________________________
  void GetTotNtriesnCh(std::vector<SimuJob> sjobs, // Input: List of SimuJob objects
		       std::vector<long double> &data)  // Output: data[0]=>Tot. Ch., data[1]=>Tot. ntries
  /* Calc. # tries & tot. Ch. from SimuJob objects */
  {
    long double totntries=0.,totcharge=0.;
    for (auto & job : sjobs) {
      if (job.is_rfile) {
	totntries += job.ntried;
	totcharge += job.charge;
      }
    }
    data = {totntries,totcharge};
  }
  //______________________________________________________________________________
  void GetTotNtriesnCh(std::vector<SimuJob> sjobs, // Input: List of SimuJob objects
		       std::string process,        // Input: Reaction process
		       std::vector<long double> &data)  // Output: data[0]=>Tot. Ch., data[1]=>Tot. ntries
  /* Calc. # tries & tot. Ch. from SimuJob objects */
  {
    long double totntries=0.,totcharge=0.;
    for (auto & job : sjobs) {
      if (process.compare(job.process)==0) {
	if (job.is_rfile) {
	  totntries += job.ntried;
	  totcharge += job.charge;
	}
      }
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
    // for efficiency and brevity
    bool is_lh2 = 1 ? target.compare("LH2")==0 : 0;
    bool is_ld2 = 1 ? target.compare("LD2")==0 : 0;
    bool is_dummy = 1 ? target.compare("Dummy")==0 : 0;
    
    double tgt_rho, tgt_dEdx_bs, tgt_dEdx_as, tgt_uwinthick, tgt_celldiam, tgt_cellthick;
    double tgt_ufoilthick, tgt_dfoilthick;
    double Al_dEdx_bs = expconst::GetdEdxCollAl(sbsconf,1), Al_dEdx_as = expconst::GetdEdxCollAl(sbsconf,0);
    double PE_dEdx_bs = expconst::GetdEdxCollPE(sbsconf,1), PE_dEdx_as = expconst::GetdEdxCollPE(sbsconf,0); 
    if (is_lh2) {
      tgt_rho = expconst::lh2_TgtRho;
      tgt_dEdx_bs = expconst::GetdEdxCollH(sbsconf,1);
      tgt_dEdx_as = expconst::GetdEdxCollH(sbsconf,0);
      tgt_uwinthick = expconst::lh2_uWinThick;
      tgt_celldiam = expconst::lh2_CellDiam;
      tgt_cellthick = expconst::lh2_CellThick;
    }else if (is_ld2) {
      tgt_rho = expconst::ld2_TgtRho;
      tgt_dEdx_bs = expconst::GetdEdxCollH(sbsconf,1) * 0.5; // Z/A =1 for H2 & 0.5 for D2 (See expconst::GetdEdxCollH)
      tgt_dEdx_as = expconst::GetdEdxCollH(sbsconf,0) * 0.5; // Z/A =1 for H2 & 0.5 for D2 (See expconst::GetdEdxCollH)
      tgt_uwinthick = expconst::ld2_uWinThick;
      tgt_celldiam = expconst::ld2_CellDiam;
      tgt_cellthick = expconst::ld2_CellThick;
    }else if (is_dummy) {
      tgt_rho = expconst::dummy_TgtRho;
      tgt_dEdx_bs = Al_dEdx_bs;
      tgt_dEdx_as = Al_dEdx_bs;
      tgt_ufoilthick = expconst::dummy_uFoilThick;
      tgt_dfoilthick = expconst::dummy_dFoilThick;
    }else
      throw std::invalid_argument("[util_pd::GetElossInTgt] Given target type is invalid!");

    // handling unphysical values ***
    // Scaling vz properly
    double vz_scaled = vz*1E2 + expconst::tgtlen*0.5; //cm
    if (vz_scaled<0) vz_scaled = 0.; // upstream of the tg 
    else if (vz_scaled>expconst::tgtlen) vz_scaled = expconst::tgtlen; //dwnstream of the tg
      
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

    double Eloss_bs = vz_scaled*tgt_rho*tgt_dEdx_bs + tgt_uwinthick*expconst::Al_Rho*Al_dEdx_bs;
    double Eloss_as = tgt_celldiam/2.0/sin(etheta)*tgt_rho*tgt_dEdx_as + tgt_cellthick/sin(etheta)*expconst::Al_Rho*Al_dEdx_as;
    if (is_dummy) { // slightly diff calculations for dummy
      if (vz_scaled<0.5*expconst::tgtlen) { // scattering just from the upstream foil
	Eloss_bs = 0.0; // nothing upstream
	Eloss_as = tgt_ufoilthick*tgt_rho*tgt_dEdx_as;
      } else { // scattering just from the downstream foil
	Eloss_bs = tgt_ufoilthick*tgt_rho*tgt_dEdx_bs; // ufoil is sitting upstream
	Eloss_as = tgt_dfoilthick*tgt_rho*tgt_dEdx_as;
      }
    }
    // add contribution from shield
    if (PE_shield_in) Eloss_as += expconst::PE_ShieldThick*expconst::PE_Rho*PE_dEdx_as;
    if (Al_shield_in) Eloss_as += expconst::Al_ShieldThick*expconst::Al_Rho*Al_dEdx_as;
    data = {Eloss_bs,Eloss_as};

    if (verbose>1) {
      std::cout << Form("\netheta = %.1f deg, vz = %.3f m, vz_scaled = %.3f m \n",etheta*TMath::RadToDeg(),vz,vz_scaled*1E-2);
      std::cout << Form("Plastic shield in place = %d, Al shield in place = %d \n",PE_shield_in,Al_shield_in);
      std::cout << Form("Mean energy loss (MeV) %.1f (before), %.1f (after) \n",Eloss_bs*1E3,Eloss_as*1E3);
    }
  }
  //______________________________________________________________________________
  double Luminosity(double ibeam,           // beam current in A
		    std::string targetType) // Valid options - "LH2" and "LD2"
  /* Calculates luminosity for g4sbs data */
  {
    double lumi = 0.;
    if (targetType.compare("LH2") == 0)
      lumi = ((ibeam/constant::qe)*expconst::tgtlen*expconst::lh2_TgtRho*(constant::N_A/constant::H2_Amass));
    else if (targetType.compare("LD2") == 0)
      lumi = ((ibeam/constant::qe)*expconst::tgtlen*expconst::ld2_TgtRho*(constant::N_A/constant::D2_Amass));
    else
      std::cerr << "[KinematicVar::Luminosity] Enter a valid target type! **!**" << std::endl;
    return lumi;
  }

  /* ############################################
     ##   Functions to manipulate histograms   ##  
     ############################################ */
  // Function to return xrange of equal statistics
  void FindEqualStatBins(TH1F *h,
			 double xlow, double xhi, int ndiv,
			 int verbose,
			 std::vector<double> &xrange) {

    int binlow = h->FindBin(xlow);
    int binhi = h->FindBin(xhi);

    double totstat = h->Integral(binlow, binhi);
    double statpercut = totstat / double(ndiv);

    xrange.push_back(xlow);
    int loweredge = binlow;
    int slices = 0;
    double accumulatedStats = 0.0;

    for (int ibin = binlow + 1; ibin <= binhi; ++ibin) {
      accumulatedStats += h->GetBinContent(ibin);

      if (accumulatedStats >= statpercut && slices < ndiv - 1) {
	double diff1 = std::abs(accumulatedStats - statpercut);
	double diff2 = std::abs((accumulatedStats - h->GetBinContent(ibin)) - statpercut);
	if (diff2 < diff1) {
	  accumulatedStats -= h->GetBinContent(ibin);
	  ibin--;
	}
	xrange.push_back(h->GetBinCenter(ibin));
	slices++;
	accumulatedStats = 0.0;
      }
    }

    xrange.push_back(xhi);

    // Final adjustment to ensure equal statistics in the last bin
    while (xrange.size() - 1 < ndiv) {
      double lastBinCenter = (xrange.back() + xhi) / 2.0;
      xrange.insert(xrange.end() - 1, lastBinCenter);
    }
    util_pd::PrintVector(xrange);

    if (verbose==1) {
      for (size_t i = 0; i < xrange.size() - 1; ++i) {
        std::cout << "Range " << i + 1 << ": " << xrange[i] << " to " << xrange[i + 1]
                  << " with integral " << h->Integral(h->FindBin(xrange[i]), h->FindBin(xrange[i + 1])) << std::endl;
      }
    }
  }

  /* #####################################
     ##   Function toCalculate Yields   ##  
     ##################################### */
    void GetYields(TH1F* hfit_p, TH1F* hfit_bg, std::vector<double> &output) {
    // TODO: Implement better error calculation
    // determining bin width and ranges
    double binW = hfit_p->GetBinWidth(1);
    double xMin = hfit_p->GetXaxis()->GetXmin();
    double xMax = hfit_p->GetXaxis()->GetXmax();
    // // calculating the integral under total fit curve
    // double totCount = gfit->Integral(xMin,xMax) / binW;
    // calculating p counts
    double pCount_err;
    double pCount = hfit_p->IntegralAndError(1,hfit_p->GetNbinsX(),pCount_err);
    pCount_err = sqrt(pCount);
    // calculating bg counts
    double bgCount_err;
    double bgCount = hfit_bg->IntegralAndError(1,hfit_bg->GetNbinsX(),bgCount_err);
    bgCount_err = sqrt(bgCount);
    // summary
    std::cout << "\n---- Various counts ----\n";
    //std::cout << "Total Count: " << totCount << "\n";
    std::cout << "p Count    : " << pCount << " +/- " << pCount_err << "\n";
    std::cout << "bg Count   : " << bgCount << " +/- " << bgCount_err << "\n";
    std::cout << "------------- \n";
    // filling output vector
    output = {pCount,pCount_err,bgCount,bgCount_err};
  }
  //______________________________________________________________________________
  void GetYields(TF1* gfit, TH1F* hfit_p, TH1F* hfit_n, TH1F* hfit_bg, std::vector<double> &output) {
    // TODO: Implement better error calculation
    // determining bin width and ranges
    double binW = hfit_p->GetBinWidth(1);
    double xMin = hfit_p->GetXaxis()->GetXmin();
    double xMax = hfit_p->GetXaxis()->GetXmax();
    // calculating the integral under total fit curve
    double totCount = gfit->Integral(xMin,xMax) / binW;
    // calculating p counts
    double pCount_err;
    double pCount = hfit_p->IntegralAndError(1,hfit_p->GetNbinsX(),pCount_err);
    pCount_err = sqrt(pCount);
    // calculating n counts
    double nCount_err;
    double nCount = hfit_n->IntegralAndError(1,hfit_n->GetNbinsX(),nCount_err);
    nCount_err = sqrt(nCount);
    // calculating total yield and error
    double totSigCount = pCount + nCount;
    double totSigCount_err = sqrt(totSigCount);
    // calculating bg counts
    double bgCount_err;
    double bgCount = hfit_bg->IntegralAndError(1,hfit_bg->GetNbinsX(),bgCount_err);
    bgCount_err = sqrt(bgCount);
    // summary
    std::cout << "\n---- Various counts ----\n";
    std::cout << "Total Count: " << totCount << "\n";
    std::cout << "p Count    : " << pCount << " +/- " << pCount_err << "\n";
    std::cout << "n Count    : " << nCount << " +/- " << nCount_err << "\n";
    std::cout << "bg Count   : " << bgCount << " +/- " << bgCount_err << "\n";
    std::cout << "------------- \n";
    // filling output vector
    output = {pCount,pCount_err,nCount,nCount_err,bgCount,bgCount_err,totSigCount,totSigCount_err};
  }
  //______________________________________________________________________________
  void GetYields(TH1F* gfit, TH1F* hfit_p, TH1F* hfit_n, TH1F* hfit_bg, std::vector<double> &output) {
    /* Calculates normalized yields */
    // TODO: Implement better error calculation
    // calculating the integral under total fit curve
    double totCount_err;
    double totCount = gfit->IntegralAndError(1,gfit->GetNbinsX(),totCount_err);
    // calculating p counts
    double pCount_err;
    double pCount = hfit_p->IntegralAndError(1,hfit_p->GetNbinsX(),pCount_err);
    pCount_err = sqrt(pCount);
    // calculating n counts
    double nCount_err;
    double nCount = hfit_n->IntegralAndError(1,hfit_n->GetNbinsX(),nCount_err);
    nCount_err = sqrt(nCount);
    // calculating total yield and error
    double totSigCount = pCount + nCount;
    double totSigCount_err = sqrt(totSigCount);
    // calculating bg counts
    double bgCount_err;
    double bgCount = hfit_bg->IntegralAndError(1,hfit_bg->GetNbinsX(),bgCount_err);
    bgCount_err = sqrt(bgCount);
    // summary
    std::cout << "\n---- Various counts ----\n";
    std::cout << "Total Count: " << totCount << "\n";
    std::cout << "p Count    : " << pCount << " +/- " << pCount_err << "\n";
    std::cout << "n Count    : " << nCount << " +/- " << nCount_err << "\n";
    std::cout << "bg Count   : " << bgCount << " +/- " << bgCount_err << "\n";
    std::cout << "------------- \n";
    // filling output vector
    output = {pCount,pCount_err,nCount,nCount_err,bgCount,bgCount_err,totSigCount,totSigCount_err};
  }

  /* ######################################################
     ##   Functions to customize fit histos and canvas   ##  
     ###################################################### */
  void customize_data(TH1F* h)
  {
    h->SetMarkerStyle(21);
    h->SetMarkerSize(0.8);
    h->SetMarkerColor(kBlack);
    h->SetLineColor(kBlack);
    h->GetXaxis()->SetLabelOffset(1.5);
    h->GetYaxis()->SetLabelSize(0.04);
  }
  //______________________________________________________________________________
  void customize_psig(TH1F* h, bool isTransp)
  {
    h->SetLineColor(kBlue);
    h->SetLineStyle(4);
    h->SetLineWidth(3);
    //***
    h->SetFillColor(kBlue);
    if (isTransp) h->SetFillColorAlpha(kBlue,0.3);
  }
  //______________________________________________________________________________
  void customize_nsig(TH1F* h, bool isTransp)
  {
    h->SetLineColor(kGreen+2);
    h->SetLineStyle(8);
    h->SetLineWidth(3);
    //***
    h->SetFillColor(kGreen+2);
    if (isTransp) h->SetFillColorAlpha(kGreen+2,0.3);
  }
  //______________________________________________________________________________
  void customize_hbg(TH1F* h, bool isTransp) 
  {
    h->SetMarkerColor(6);
    h->SetMarkerStyle(29);
    h->SetLineColor(6);
    h->SetLineStyle(9);
    h->SetLineWidth(3);
    //***
    h->SetFillColor(6);
    if (isTransp) h->SetFillColorAlpha(6,0.3);
  }
  //______________________________________________________________________________
  void customize_gfit(TH1F* h, bool isTransp)
  {
    h->SetLineColor(kRed);
    h->SetLineStyle(7);
    h->SetLineWidth(3);
    //***
    h->SetFillColor(kRed);
    if (isTransp) h->SetFillColorAlpha(kRed,0.2);
  }
  //______________________________________________________________________________
  void customize_residual(TH1F* h)
  {
    h->GetXaxis()->SetLabelOffset(0.03);
    h->GetXaxis()->SetLabelSize(0.12);
    h->GetXaxis()->SetTitle("#font[32]{#Deltax} (m)");
    h->GetXaxis()->SetTitleOffset(1);
    h->GetXaxis()->SetTitleSize(0.15);
    h->GetXaxis()->CenterTitle();
    h->GetYaxis()->SetLabelOffset(0.005);
    h->GetYaxis()->SetLabelSize(0.11);
    h->SetMarkerStyle(22);
    h->SetMarkerColor(46);
    h->SetLineColor(46);
    h->SetStats(0);
  }
}
