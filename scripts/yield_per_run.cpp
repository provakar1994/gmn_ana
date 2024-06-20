/* This macro will be used to extract yields per run 
   It will be configures by con_fit_dx.json. Will look
   for prodection keys ie qeals_prod or elas_prod. 
   -------
   P. Datta Created 06-19-2024
*/

#include <algorithm>
#include <unordered_map>

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"

#include "gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

void gStyleFitCanvas() 
{
  gStyle->SetOptStat("e"); gStyle->SetOptFit(1); 
  gStyle->SetErrorX(0);
}

//______________________________________________________________________________
TFile * ReadRootFile(char const * filename) {
  // Open the ROOT file
  TFile *file = TFile::Open(filename, "READ");
    
  // Check if the file is open and the histogram exists
  if (!file || file->IsZombie()){ //|| !file->GetListOfKeys()->Contains(histname)) {
    std::cerr << "Error: Failed to open the file or histogram not found!" << std::endl;
    throw;
  }

  return file;

}

//______________________________________________________________________________
void FitGraphWithConstant(TGraph* graph) {
  // Define a constant function
  TF1* constantFit = new TF1("constantFit", "[0]", graph->GetX()[0], graph->GetX()[graph->GetN()-1]);
  // Set initial parameter value (optional)
  constantFit->SetParameter(0, 0.5);  // Initial guess for constant value
  // Fit the graph with the constant function
  graph->Fit("constantFit", "Q");  // "Q" for quiet mode, i.e., no drawing
  // Print fit results
  double constantValue = constantFit->GetParameter(0);
  double constantError = constantFit->GetParError(0);
  std::cout << "Fit Result:" << std::endl;
  std::cout << "Constant Value: " << constantValue << " +/- " << constantError << std::endl;
  // Optionally, draw the graph and fit function on a canvas
  graph->Draw("AP");  // Draw graph with axes and points
  constantFit->Draw("SAME");  // Draw fit function on top of graph
}

//______________________________________________________________________________
void DrawTGraphWithErrors(const std::vector<double>& xValues, const std::vector<double>& yValues, const std::vector<double>& yErrors,
			  std::string gtitle, std::string gxtitle, std::string gytitle, bool doFit=1,
			  bool yrange=0, double ylow=0, double yhi=0) {
  // Check if the sizes of the arrays match
  if (xValues.size() != yValues.size() || xValues.size() != yErrors.size()) {
    std::cerr << "Error: Array sizes do not match!" << std::endl;
    return;
  }
  // Create a TGraphErrors object
  TGraphErrors *graph = new TGraphErrors(xValues.size(), &xValues[0], &yValues[0], 0, &yErrors[0]);
  // Set graph attributes
  graph->SetTitle(gtitle.c_str());
  graph->GetXaxis()->SetTitle(gxtitle.c_str());
  graph->GetYaxis()->SetTitle(gytitle.c_str());
  // Center the axis titles
  graph->GetXaxis()->CenterTitle(true);
  graph->GetYaxis()->CenterTitle(true);
  if (yrange) {
    graph->GetYaxis()->SetRangeUser(ylow,yhi);
  }
  // Create a canvas and draw the graph
  graph->Draw("AP");  // A: Draw axis, P: Draw points
  // Optionally, add a legend or additional customization
  graph->SetMarkerStyle(20);
  graph->GetYaxis()->SetMaxDigits(3);
  // Fit graph
  if (doFit) FitGraphWithConstant(graph);
}

//______________________________________________________________________________
TCanvas *FitYields( std::string const &canvname, std::string const &bgshape,
		    std::vector<double> const &rnums, std::vector<double> const &stats,
		    std::vector<double> const &R_fit, std::vector<double> const &R_fit_err,
		    std::vector<double> const &CNpY, std::vector<double> const &CNpYerr,
		    std::vector<double> const &CNnY, std::vector<double> const &CNnYerr,
		    std::vector<double> const &CNtotY, std::vector<double> const &CNtotYerr) {

  TCanvas *cfit = new TCanvas(canvname.c_str(),bgshape.c_str(), 1200, 1000);
  cfit->Divide(3,2);
  gStyleFitCanvas(); cfit->SetGridy();  cfit->SetGridx();
  // 
  std::vector<double> statsErr;
  for (int i=0; i<rnums.size(); i++) {
    statsErr.push_back(0);
  }
  //
  cfit->cd(1); // Data statistics
  DrawTGraphWithErrors(rnums,stats,statsErr,"# Entries in the Data Histogram","Run Number","# of Entries",0,1,0,70E3);
  cfit->cd(2); // R_fit
  DrawTGraphWithErrors(rnums,R_fit,R_fit_err,Form("R_fit vs Runs | %s",bgshape.c_str()),"Run Number","R_fit");
  cfit->cd(3); // Charge-normalized p yield
  DrawTGraphWithErrors(rnums,CNpY,CNpYerr,Form("Normalized p Yield vs Runs | %s",bgshape.c_str()),"Run Number","Normalized p Yield (1/C)",1,1,0,1E7);
  cfit->cd(4); // Charge-normalized n yield
  DrawTGraphWithErrors(rnums,CNnY,CNnYerr,Form("Normalized n Yield vs Runs | %s",bgshape.c_str()),"Run Number","Normalized n Yield (1/C)",1,1,0,1E7);
  cfit->cd(5); // Charge-normalized n yield
  DrawTGraphWithErrors(rnums,CNtotY,CNtotYerr,Form("Normalized n+p Yield vs Runs | %s",bgshape.c_str()),"Run Number","Normalized n+p Yield (1/C)",1,1,0,1E7);
  return cfit;
}

//______________________________________________________________________________
TCanvas *FitYields( std::string const &canvname, std::string const &bgshape,
		    std::vector<double> const &rnums, std::vector<double> const &stats,
		    std::vector<double> const &CNpY, std::vector<double> const &CNpYerr) {

  TCanvas *cfit = new TCanvas(canvname.c_str(),bgshape.c_str(), 1000, 800);
  cfit->Divide(1,2);
  gStyleFitCanvas(); cfit->SetGridy();  cfit->SetGridx();
  // 
  std::vector<double> statsErr;
  for (int i=0; i<rnums.size(); i++) {
    statsErr.push_back(0);
  }
  //
  cfit->cd(1); // Data statistics
  DrawTGraphWithErrors(rnums,stats,statsErr,"# Entries in the Data Histogram","Run Number","# of Entries",0);
  cfit->cd(2); // Charge-normalized p yield
  DrawTGraphWithErrors(rnums,CNpY,CNpYerr,Form("Normalized p Yield vs Runs | %s",bgshape.c_str()),"Run Number","Normalized p Yield (1/C)",1,1,0,3E7);
  return cfit;
}

//______________________________________________________________________________
void yield_per_run (const char *configfilename, 
		    bool is_elastic = 1) // 1=>Yes, 0=>QE
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  char const * key = is_elastic ? "elas" : "qelas";

  // reading parent config file ---------------------------------------
  std::string config_p = jmgr->GetValueFromSubKey_str(key,"parent_config"); 
  std::string key_p = jmgr->GetValueFromSubKey_str(key,"parent_config_key");
  JSONManager *jmgr_p = new JSONManager(config_p.c_str());  

  // reading in proper data and simu output files based on parent config
  int conf = jmgr_p->GetValueFromSubKey<int>(key_p.c_str(),"SBS_config");
  int sbsmag = jmgr_p->GetValueFromSubKey<int>(key_p.c_str(),"SBS_magnet_percent");
  int model = jmgr_p->GetValueFromSubKey<int>(key_p.c_str(),"model");
  int pass = jmgr_p->GetValueFromSubKey<int>(key_p.c_str(),"pass");
  std::string file_prefix = jmgr_p->GetValueFromSubKey_str(key_p.c_str(),"output_filebase");
  std::string target = is_elastic ? "LH2" : "LD2";

  // reading analysis params from parent config
  vector<double> dx_fit_range; jmgr_p->GetVectorFromSubKey<double>(key_p.c_str(),"dx_fit_range",dx_fit_range);

  // reading stuff from the local config
  int nruns = -1;
  bool has_BC_cut = jmgr->GetValueFromSubKey<int>(key,"has_beamcurr_cut");
  std::vector<int> runs_to_exclude; jmgr->GetVectorFromSubKey<int>(key,"runs_to_exclude",runs_to_exclude);
  //= {11449,11451,11452}; // = {11436,11616};

  // input and output filename format
  char const * in_script = "fit_dx";
  char const * out_script = "yield_per_run";
  std::string rootdir_path = Form("pdout/fits/sbs%dsbs%dp",conf,sbsmag);
  file_prefix = has_BC_cut ? "BCcut_" + file_prefix : file_prefix;
  TString inFile = Form("%s/%s_%s_%s_pass%d_simc_sbs%d_sbs%dp_model%d.root"
			,rootdir_path.c_str(),file_prefix.c_str(),in_script,key,pass,conf,sbsmag,model);
  TString outFile = Form("%s/%s_%s_%s_pass%d_simc_sbs%d_sbs%dp_model%d.root"
			 ,rootdir_path.c_str(),file_prefix.c_str(),out_script,key,pass,conf,sbsmag,model);

  //reading input file
  TFile * file = ReadRootFile(inFile);
  
  // Reading run list
  std::vector<CodaRun> cruns; util_pd::ReadRunList("../DB",nruns,conf,target,pass,sbsmag,0,cruns);
  // implementing hashtable w/ run number as key, for efficiency
  std::unordered_map<int,CodaRun> mrun;
  for (auto & crun: cruns) mrun[crun.runnum] = crun;
  
  // Defining histos to fit
  TH1F * h_dxHCAL_simu_p = (TH1F*)file->Get("h_dxHCAL_simu_p");
  TH1F * h_dxHCAL_simu_n = (TH1F*)file->Get("h_dxHCAL_simu_n");
  TH1F * h_dxHCAL_bg_data = (TH1F*)file->Get("h_dxHCAL_bg_data");
  TH1F * h_dxHCAL_bg_inel_p = (TH1F*)file->Get("h_dxHCAL_bg_inel_p");
  TH1F * h_dxHCAL_bg_inel_n;
  if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)file->Get("h_dxHCAL_bg_inel_n");
  TH1F * h_dxHCAL_bg_inel = (TH1F*)file->Get("h_dxHCAL_bg_inel");
  // vs Rnum histo
  TH2F * h_dxHCAL_vs_rnum = (TH2F*)file->Get("h2_dxHCAL_vs_rnum");

  // Define arrays to store fit results
  int nBinsX = h_dxHCAL_vs_rnum->GetNbinsX();
  std::vector<double> rnums, stats;
  std::vector<double> R2, Rerr2, CNpY2, CNpYerr2, CNnY2, CNnYerr2, CNtotY2, CNtotYerr2; 
  std::vector<double> R3, Rerr3, CNpY3, CNpYerr3, CNnY3, CNnYerr3, CNtotY3, CNtotYerr3;

  // creating output file
  TFile *fout = new TFile(outFile, "RECREATE");
  
  int totc = 0; int gc = 0;
  // Loop over each x-bin
  for (int i = 1; i <= nBinsX; ++i) {
    totc++;

    // skipping empty bins
    if (h_dxHCAL_vs_rnum->GetBinContent(i)!=0) {
      gc++; // counts the # of runs being analyzed

      // extracting the run number for this bin
      int rnum = (int)h_dxHCAL_vs_rnum->GetXaxis()->GetBinCenter(i);
      // check if want to include this run in the analysis or not
      if (!runs_to_exclude.empty()) {
	if (std::find(runs_to_exclude.begin(), runs_to_exclude.end(), rnum) != runs_to_exclude.end())
	  continue;
      }
      rnums.push_back(rnum);

      // getting the charge and DAQ livetime for the run
      double charge = mrun[rnum].charge;
      double daqlt = mrun[rnum].daqlvtm;

      // Project slice along y-axis
      TH1F *sliceY = (TH1F*)h_dxHCAL_vs_rnum->ProjectionY(Form("hfit_run%d",rnum), i, i);
      sliceY->SetTitle(Form("%d",rnum));
      stats.push_back(sliceY->GetEntries());


      /*##############################################
	## Fitting elastic dx distributions per run ##
	############################################## */
      if (is_elastic) {
	TCanvas *cA = util_pd::TC("cA",1,2); gStyleFitCanvas();
	// --------
	// Fitting data/MC w/ polynomial ---
	// --------
	cA->cd(1);
	vector<TH1F*> ho2;
	TF1 *f2 = fit::fit_1hs_1pbg_THI(dx_fit_range,
					sliceY,h_dxHCAL_simu_p,2,
					ho2);
	// recording fit params
	std::cout << Form("Pol2 Bg: Run #: %d, R: %f",rnum,f2->GetParameter(1)) << "\n";
	R2.push_back(f2->GetParameter(1)); Rerr2.push_back(f2->GetParError(1));
	// drawing fit histos
	ho2[0]->Draw("E"); util_pd::customize_data(ho2[0]);
	ho2[1]->Draw("same HIST"); util_pd::customize_nsig(ho2[1],0);
	ho2[2]->Draw("same HIST"); util_pd::customize_hbg(ho2[2],0);
	// calculating yields
	std::vector<double> yo2; util_pd::GetYields(f2,ho2[1],ho2[2],yo2);
	CNpY2.push_back(yo2[0]/charge/daqlt); CNpYerr2.push_back(yo2[1]/charge/daqlt);
	// calculating charge normalized yields
	std::cout << "---- \n Run: " << rnum << "\n";
	std::cout << "Charge: " << charge << " DAQ LT: " << daqlt << "\n";
	std::cout << "Charge norm. p Yield: " << yo2[0]/charge/daqlt << " +/- " << yo2[1]/charge/daqlt << "\n---\n";
	// grad the stat box of the fitted histo
	cA->Update();
	TPaveStats *st2 = (TPaveStats*)ho2[0]->FindObject("stats");
	st2->SetOptStat(1); st2->SetOptFit(1111); 
	//
	ho2[0]->Write();
	// ----
	cA->Write();
      }

      /*#########################################
	## Fitting QE dx distributions per run ##
	######################################### */      
      else {      
	TCanvas *cA = util_pd::TC("cA",1,2); gStyleFitCanvas();
	// --------
	// Fitting data/MC w/ background from MC ---
	// --------
	cA->cd(1);
	vector<TH1F*> ho2;
	TF1 *f2 = fit::fit_2hs_2hbg_THI(dx_fit_range,
					sliceY,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n,
					ho2);
	// recording fit params
	std::cout << Form("MC Bg: Run #: %d, R: %f",rnum,f2->GetParameter(1)) << "\n";
	R2.push_back(f2->GetParameter(1)); Rerr2.push_back(f2->GetParError(1));
	// drawing fit histos
	ho2[0]->SetTitle(Form("MC Bg | Run: %d",rnum));
	ho2[0]->GetXaxis()->SetTitle("#Deltax (m)");
	ho2[0]->Draw("E"); util_pd::customize_data(ho2[0]);
	ho2[4]->Draw("same HIST"); util_pd::customize_psig(ho2[4],0);
	ho2[5]->Draw("same HIST"); util_pd::customize_nsig(ho2[5],0);
	ho2[2]->Draw("same HIST"); util_pd::customize_hbg(ho2[2],0);
	// calculating yields
	std::vector<double> yo2; util_pd::GetYields(f2,ho2[4],ho2[5],ho2[2],yo2);
	CNpY2.push_back(yo2[0]/charge/daqlt); CNpYerr2.push_back(yo2[1]/charge/daqlt);
	CNnY2.push_back(yo2[2]/charge/daqlt); CNnYerr2.push_back(yo2[3]/charge/daqlt);
	CNtotY2.push_back(yo2[6]/charge/daqlt); CNtotYerr2.push_back(yo2[7]/charge/daqlt);
	// grad the stat box of the fitted histo
	cA->Update();
	TPaveStats *st2 = (TPaveStats*)ho2[0]->FindObject("stats");
	st2->SetOptStat(1); st2->SetOptFit(1111); 
	//
	ho2[0]->Write();
	// --------
	// Fitting data/MC w/ polynomial background (2nd order - fixed) ---
	// --------
	cA->cd(2);
	vector<TH1F*> ho3;
	TF1 *f3 = fit::fit_2hs_1pbg_THI(dx_fit_range,
					sliceY,h_dxHCAL_simu_p,h_dxHCAL_simu_n,2,
					ho3);
	// recording fit params
	std::cout << Form("Pol2 Bg: Run #: %d, R: %f",rnum,f3->GetParameter(1)) << "\n";
	R3.push_back(f3->GetParameter(1)); Rerr3.push_back(f3->GetParError(1));
	// drawing fit histos
	ho3[0]->SetTitle(Form("Pol2 Bg | Run: %d",rnum));
	ho3[0]->GetXaxis()->SetTitle("#Deltax (m)");
	ho3[0]->Draw("E"); util_pd::customize_data(ho3[0]);
	ho3[4]->Draw("same HIST"); util_pd::customize_psig(ho3[4],0);
	ho3[5]->Draw("same HIST"); util_pd::customize_nsig(ho3[5],0);
	ho3[2]->Draw("same HIST"); util_pd::customize_hbg(ho3[2],0);
	// calculating yields
	std::vector<double> yo3; util_pd::GetYields(f3,ho3[4],ho3[5],ho3[2],yo3);
	CNpY3.push_back(yo3[0]/charge/daqlt); CNpYerr3.push_back(yo3[1]/charge/daqlt);
	CNnY3.push_back(yo3[2]/charge/daqlt); CNnYerr3.push_back(yo3[3]/charge/daqlt);
	CNtotY3.push_back(yo3[6]/charge/daqlt); CNtotYerr3.push_back(yo3[7]/charge/daqlt);
	// grad the stat box of the fitted histo
	cA->Update();
	TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
	st3->SetOptStat(1); st3->SetOptFit(1111); 
	//
	ho3[0]->Write();
	// ----
	cA->Write();
      }
    }
  } // end of loop over bins ie runs

  // Time to plot and fit yields vs run
  if (is_elastic) {
    TCanvas *cfit2 = FitYields("cfit2","Pol2 Bg",rnums,stats,CNpY2,CNpYerr2);
    cfit2->Update(); cfit2->Write();
  }
  else {
    TCanvas *cfit2 = FitYields("cfit2","MC Bg",rnums,stats,R2,Rerr2,CNpY2,CNpYerr2,CNnY2,CNnYerr2,CNtotY2,CNtotYerr2);
    cfit2->Update(); cfit2->Write();
    TCanvas *cfit3 = FitYields("cfit3","Pol2 Bg",rnums,stats,R3,Rerr3,CNpY3,CNpYerr3,CNnY3,CNnYerr3,CNtotY3,CNtotYerr3);
    cfit3->Update(); cfit3->Write();
  }
  std::cout << totc << " " << gc << "\n";

  // reporting output files
  std::cout << "------" << std::endl;
  //std::cout << " Fit params  : " << outData << std::endl;
  //std::cout << " Summary plots  : " << outPlot << std::endl;
  std::cout << " Output ROOT file  : " << outFile << std::endl;
  std::cout << "------" << std::endl;
}



// calculating charge normalized yields
// std::cout << "---- \n Run: " << rnum << "\n";
// std::cout << "Charge: " << charge << " DAQ LT: " << daqlt << "\n";
// std::cout << "Charge norm. p Yield: " << yo2[0]/charge/daqlt << "\n";
// std::cout << "Charge norm. n Yield: " << yo2[2]/charge/daqlt << "\n----\n";
