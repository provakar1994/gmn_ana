/* This macro will be used to extract yields per run 
   It will be configures by con_fit_dx.json. Will look
   for prodection keys ie qeals_prod or elas_prod. 
   -------
   P. Datta Created 06-19-2024
*/

/*
  NOTES:
  -----
  * Map between canvas number and the backgorund model:
  ** For elastics:
  *** 2 => 2nd order polynomial
  *** 3 => Side band fit
  ** For QE:
  *** 2 => 2nd order polynomial
  *** 3 => Side band if SBS field is ZERO
  *** 3 => Inelastic bg. from MC if SBS field is Non-ZERO
  *** 4 => Bg. shape from data
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

bool temp = 0;

//______________________________________________________________________________
void gStyleFitCanvas() 
{
  gStyle->SetOptStat("e"); gStyle->SetOptFit(1); 
  gStyle->SetErrorX(0);
}

//______________________________________________________________________________
void FurtherCustomizeDxHisto(TH1F* h, std::string const &bgshape, int rnum)
{
  h->SetTitle(Form("%s | Run: %d",bgshape.c_str(),rnum));
  h->GetXaxis()->SetTitle("#Deltax (m)");
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
TF1 * FitGraphWithConstant(TGraph* graph) {
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
  //std::vector<std::double> fitpar = {constantFit->GetParameter(0),constantFit->GetParError(0)};
  return constantFit;
}

//______________________________________________________________________________
TGraph * DrawTGraphWithErrors(const std::vector<double>& xValues,
			      const std::vector<double>& yValues, const std::vector<double>& yErrors,
			      std::string gtitle, std::string gxtitle, std::string gytitle,
			      bool yrange=0, double ylow=0, double yhi=0) {
  // Check if the sizes of the arrays match
  if (xValues.size() != yValues.size() || xValues.size() != yErrors.size()) {
    std::cerr << "Error: Array sizes do not match!" << std::endl;
    throw;
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
  return graph;
}

//______________________________________________________________________________
TF1 * DrawTGraphAndFit(const std::vector<double>& xValues, const std::vector<double>& yValues, const std::vector<double>& yErrors,
		      std::string gtitle, std::string gxtitle, std::string gytitle,
		      bool yrange=0, double ylow=0, double yhi=0) {
  TGraph * graph = DrawTGraphWithErrors(xValues,yValues,yErrors,gtitle,gxtitle,gytitle,yrange,ylow,yhi);
  
  TF1 * fit = FitGraphWithConstant(graph);
  return fit;
}

//______________________________________________________________________________
void DrawSummaryCanv(std::string const &config, std::string const &bgshape,
		     double const Yexp, double const Ydata, double const YdataErr) {
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText(Form(" %s Runs",config.c_str()));
  pt->AddText(Form(" Bg. Shape : %s",bgshape.c_str()));
  pt->AddText(Form(" Total Yield from Data : %.0f #pm %.0f",Ydata,YdataErr));
  pt->AddText(Form(" Total Yield Expected : %.0f",Yexp));
  pt->AddText(Form(" Efficiency : %.1f%%",(Ydata*100)/Yexp));
  TText *t1 = pt->GetLineWith("Runs"); t1->SetTextColor(kRed+2);
  TText *t2 = pt->GetLineWith("Efficiency"); t2->SetTextColor(kBlue);
  pt->Draw();
}

//______________________________________________________________________________
TCanvas *FitYields( std::string const &canvname, std::string const &config,
		    std::string const &bgshape, double const exp_yield,
		    double const h_yield_range_hi, double const h_stats_range_hi,
		    std::vector<double> const &rnums, std::vector<double> const &stats,
		    std::vector<double> const &R_fit, std::vector<double> const &R_fit_err,
		    std::vector<double> const &CNpY, std::vector<double> const &CNpYerr,
		    std::vector<double> const &CNnY, std::vector<double> const &CNnYerr,
		    std::vector<double> const &CNtotY, std::vector<double> const &CNtotYerr) {
  /* Meant for LD2 data */
  TCanvas *cfit = new TCanvas(canvname.c_str(),bgshape.c_str(), 1200, 900);
  cfit->Divide(3,2);
  gStyleFitCanvas(); cfit->SetGridy();  cfit->SetGridx();
  // 
  std::vector<double> statsErr;
  for (int i=0; i<rnums.size(); i++) {
    statsErr.push_back(0);
  }
  //
  cfit->cd(1); // Data statistics
  TGraph *g1 = DrawTGraphWithErrors(rnums,stats,statsErr,"# Entries in the Data Histogram",
				    "Run Number","# of Entries",1,0,h_stats_range_hi);
  cfit->cd(2); // R_fit
  TF1 *f1 = DrawTGraphAndFit(rnums,R_fit,R_fit_err,Form("R_fit vs Runs | %s",bgshape.c_str()),"Run Number","R_fit");
  cfit->cd(3); // Charge-normalized p yield
  TF1 *f2 = DrawTGraphAndFit(rnums,CNpY,CNpYerr,Form("Normalized p Yield vs Runs | %s",bgshape.c_str()),"Run Number",
			     "Normalized p Yield (1/C)",1,0,h_yield_range_hi);
  cfit->cd(4); // Charge-normalized n yield
  TF1 *f3 = DrawTGraphAndFit(rnums,CNnY,CNnYerr,Form("Normalized n Yield vs Runs | %s",bgshape.c_str()),"Run Number",
			     "Normalized n Yield (1/C)",1,0,h_yield_range_hi);
  cfit->cd(5); // Charge-normalized n+p yield
  TF1 *f4 = DrawTGraphAndFit(rnums,CNtotY,CNtotYerr,Form("Normalized n+p Yield vs Runs | %s",bgshape.c_str()),"Run Number",
			     "Normalized n+p Yield (1/C)",1,0,h_yield_range_hi);
  cfit->cd(6);
  double Ydata = f4->GetParameter(0); double YdataErr = f4->GetParError(0);
  DrawSummaryCanv(config,bgshape,exp_yield,Ydata,YdataErr);
  return cfit;
}

//______________________________________________________________________________
TCanvas *FitYields( std::string const &canvname, std::string const &config,
		    std::string const &bgshape, double const exp_yield,
		    double const h_yield_range_hi, double const h_stats_range_hi,
		    std::vector<double> const &rnums, std::vector<double> const &stats,
		    std::vector<double> const &CNpY, std::vector<double> const &CNpYerr) {
  /* Meant for LH2 data */
  TCanvas *cfit = new TCanvas(canvname.c_str(),bgshape.c_str(), 1200, 900);
  cfit->Divide(3,2);
  gStyleFitCanvas(); cfit->SetGridy();  cfit->SetGridx();
  // 
  std::vector<double> statsErr;
  for (int i=0; i<rnums.size(); i++) {
    statsErr.push_back(0);
  }
  //
  cfit->cd(1); // Data statistics
  TGraph * g1 = DrawTGraphWithErrors(rnums,stats,statsErr,"# Entries in the Data Histogram","Run Number","# of Entries");
  cfit->cd(2); // Charge-normalized p yield
  TF1 * f1 = DrawTGraphAndFit(rnums,CNpY,CNpYerr,Form("Normalized p Yield vs Runs | %s",bgshape.c_str()),"Run Number",
			      "Normalized p Yield (1/C)",1,0,h_yield_range_hi);
  cfit->cd(3); // summary
  double Ydata = f1->GetParameter(0); double YdataErr = f1->GetParError(0);
  DrawSummaryCanv(config,bgshape,exp_yield,Ydata,YdataErr);
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
  std::string config = "SBS" + to_string(conf) + "-" + to_string(sbsmag) + "p " + target;

  // reading analysis params from parent config
  vector<double> dx_fit_range; jmgr_p->GetVectorFromSubKey<double>(key_p.c_str(),"dx_fit_range",dx_fit_range);
  vector<double> reject_points; jmgr_p->GetVectorFromSubKey<double>(key_p.c_str(),"SB_reject_points",reject_points);

  double RpMC = 0.4;

  // reading stuff from the local config
  int nruns = -1;
  bool has_BC_cut = jmgr->GetValueFromSubKey<int>(key,"has_beamcurr_cut");
  std::vector<int> runs_to_exclude; jmgr->GetVectorFromSubKey<int>(key,"runs_to_exclude",runs_to_exclude);
  double hR_yield = jmgr->GetValueFromSubKey<double>(key,"h_yield_range_hi");
  double hR_stats = jmgr->GetValueFromSubKey<double>(key,"h_stats_range_hi");

  // ^ 1st entry of the vector is treated as bool

  // input and output filename format
  char const * in_script = "fit_dx";
  char const * out_script = "yield_per_run";
  std::string rootdir_path = Form("pdout/fits/sbs%dsbs%dp",conf,sbsmag);
  file_prefix = has_BC_cut ? "BCcut_" + file_prefix : file_prefix;
  TString inFile = Form("%s/%s_%s_%s_pass%d_simc_sbs%d_sbs%dp_model%d.root"
			,rootdir_path.c_str(),file_prefix.c_str(),in_script,key,pass,conf,sbsmag,model);
  TString outFile = Form("%s/%s_%s_%s_pass%d_simc_sbs%d_sbs%dp_model%d.root"
			 ,rootdir_path.c_str(),file_prefix.c_str(),out_script,key,pass,conf,sbsmag,model);
  TString outPlot = outFile; outPlot.ReplaceAll(".root",".pdf");
  TString outData = outFile; outData.ReplaceAll(".root",".csv");
  ofstream outdata; outdata.open(outData);
  
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
  // charge normalized simu histos
  TH1F * h_dxHCAL_simu_p_norm = (TH1F*)file->Get("h_dxHCAL_simu_p_norm");
  TH1F * h_dxHCAL_simu_n_norm = (TH1F*)file->Get("h_dxHCAL_simu_n_norm");
  TH1F * h_dxHCAL_simu_norm = (TH1F*)h_dxHCAL_simu_p_norm->Clone();
  if (!is_elastic) h_dxHCAL_simu_norm->Add(h_dxHCAL_simu_p_norm,h_dxHCAL_simu_n_norm);
  // bg histos
  TH1F * h_dxHCAL_bg_data = (TH1F*)file->Get("h_dxHCAL_bg_data");
  TH1F * h_dxHCAL_bg_inel_p = (TH1F*)file->Get("h_dxHCAL_bg_inel_p");
  TH1F * h_dxHCAL_bg_inel_n;
  if (!is_elastic) h_dxHCAL_bg_inel_n = (TH1F*)file->Get("h_dxHCAL_bg_inel_n");
  TH1F * h_dxHCAL_bg_inel = (TH1F*)file->Get("h_dxHCAL_bg_inel");
  // vs Rnum histo
  TH2F * h_dxHCAL_vs_rnum = (TH2F*)file->Get("h2_dxHCAL_vs_rnum");
  TH2F * h_dxHCAL_vs_rnum_norm; // charge norm & live time corr data
  if (temp) h_dxHCAL_vs_rnum_norm = (TH2F*)file->Get("h2_dxHCAL_vs_rnum_norm");
  
  // ----
  // Calculating expected yield from MC
  double exp_yield = h_dxHCAL_simu_norm->Integral();
  // ----
  
  // Define arrays to store fit results
  int nBinsX = h_dxHCAL_vs_rnum->GetNbinsX();
  std::vector<double> rnums, stats;
  std::vector<double> R2, Rerr2, CNpY2, CNpYerr2, CNnY2, CNnYerr2, CNtotY2, CNtotYerr2; 
  std::vector<double> R3, Rerr3, CNpY3, CNpYerr3, CNnY3, CNnYerr3, CNtotY3, CNtotYerr3;
  std::vector<double> R4, Rerr4, CNpY4, CNpYerr4, CNnY4, CNnYerr4, CNtotY4, CNtotYerr4;

  // creating output file
  TFile *fout = new TFile(outFile, "RECREATE");

  bool fRun = true;
  //bool is_heading = true;
  int totc = 0; int gc = 0;
  // ******** -------- *********
  // * Loop over each x-bin
  // * Each x-bin correspond to a run and the corresponding y-bin has data dx.
  // * There may be empty bins as the good run numbers aren't continuous
  //
  for (int i = 1; i <= nBinsX; ++i) {
    totc++;

    // skipping empty bins
    if (h_dxHCAL_vs_rnum->GetBinContent(i)!=0) {
      gc++; // counts the # of runs being analyzed
      
      // extracting the run number for this bin
      int rnum = (int)h_dxHCAL_vs_rnum->GetXaxis()->GetBinCenter(i);
      // check if want to include this run in the analysis or not
      if (!runs_to_exclude.empty()) {
	if (runs_to_exclude[0]!=0) { // first entry of the vector is treated a bool
	  if (std::find(runs_to_exclude.begin(), runs_to_exclude.end(), rnum) != runs_to_exclude.end())
	    continue;
	}
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
	// Initializing output csv file for elastics
	if (fRun) {
	  outdata << "runnum,"
		  << "Yp2,Yp2err,"
		  << "Yp3,Yp3err,"
		  << "\n";
	}
	outdata << rnum << ",";
	// --
	// Defining the dx fit canvas
	TCanvas *cA = new TCanvas("cA","cA",1200,1000);
	cA->Divide(2,1); gStyleFitCanvas();
	// --------
	// Fitting data/MC w/ polynomial ---
	// --------
	cA->cd(1);
	vector<TH1F*> ho2;
	TF1 *f2 = fit::fit_1hs_1pbg_THI(dx_fit_range,
					sliceY,h_dxHCAL_simu_p,2,
					ho2);
	// recording fit params
	std::cout << Form("Pol2 Bg: Run #: %d",rnum) << "\n";
	// drawing fit histos
	ho2[0]->Draw("E"); util_pd::customize_data(ho2[0]);
	FurtherCustomizeDxHisto(ho2[0],"Pol2 Bg",rnum);
	ho2[1]->Draw("same HIST"); util_pd::customize_nsig(ho2[1],0);
	ho2[2]->Draw("same HIST"); util_pd::customize_hbg(ho2[2],0);
	// calculating yields
	std::vector<double> yo2; util_pd::GetYields(ho2[1],ho2[2],yo2);
	outdata << Form("%.0f,%.0f,",yo2[0],yo2[1]);
	CNpY2.push_back(yo2[0]/charge/daqlt); CNpYerr2.push_back(yo2[1]/charge/daqlt);
	// grad the stat box of the fitted histo
	cA->Update();
	TPaveStats *st2 = (TPaveStats*)ho2[0]->FindObject("stats");
	st2->SetOptStat(1); st2->SetOptFit(1111); 
	//
	ho2[0]->Write();
	// ----
	// --------
	// Fitting data/MC w/ polynomial ---
	// --------
	cA->cd(2);
	vector<TH1F*> ho3;
	TF1* bg3 = fit::fit_1pbg_SB(dx_fit_range,
				    reject_points,
				    3,//Opoly,
				    fit::GetFitParams(f2),
				    sliceY,
				    ho3);
	// recording fit params
	std::cout << Form("Side Band: Run #: %d",rnum) << "\n";
	// drawing fit histos
	ho3[0]->Draw("E"); util_pd::customize_data(ho3[0]);
	FurtherCustomizeDxHisto(ho3[0],"Side Band",rnum);
	ho3[1]->Draw("same HIST"); util_pd::customize_nsig(ho3[1],0);
	ho3[2]->Draw("same HIST"); util_pd::customize_hbg(ho3[2],0);
	// calculating yields
	std::vector<double> yo3; util_pd::GetYields(ho3[1],ho3[2],yo3);
	outdata << Form("%.0f,%.0f,",yo3[0],yo3[1]);
	CNpY3.push_back(yo3[0]/charge/daqlt); CNpYerr3.push_back(yo3[1]/charge/daqlt);
	// calculating charge normalized yields
	std::cout << "---- \n Run: " << rnum << "\n";
	std::cout << "Charge: " << charge << " DAQ LT: " << daqlt << "\n";
	std::cout << "Charge norm. p Yield: " << yo3[0]/charge/daqlt << " +/- " << yo3[1]/charge/daqlt << "\n---\n";
	// grad the stat box of the fitted histo
	cA->Update();
	TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
	st3->SetOptStat(1); st3->SetOptFit(1111); 
	//
	ho3[0]->Write();
	// ----
	outdata << "\n";
	// ----
	cA->Write();
	if (fRun) { cA->SaveAs(Form("%s[",outPlot.Data())); fRun = false; }
	cA->SaveAs(Form("%s",outPlot.Data()));
      }

      /*#########################################
	## Fitting QE dx distributions per run ##
	######################################### */      
      else {      
	// Initializing output csv file for QE
	if (fRun) {
	  outdata << "runnum,RpMC,"
		  << "R2,R2err,Yp2,Yp2err,Yn2,Yn2err,Ytot2,Ytot2err,"
		  << "R3,R3err,Yp3,Yp3err,Yn3,Yn3err,Ytot3,Ytot3err,";
	  if (sbsmag!=0) {
	    outdata << "R4,R4err,Yp4,Yp4err,Yn4,Yn4err,Ytot4,Ytot4err,"
		    << "R5,R5err,Yp5,Yp5err,Yn5,Yn5err,Ytot4,Ytot4err,";
	  }
	  outdata << "\n";
	}
	outdata << rnum << "," << RpMC << ",";
	// --
	// Defining the dx fit canvas
	TCanvas *cA = new TCanvas("cA","cA",1200,1200);
	cA->Divide(2,2); gStyleFitCanvas();
	// --------
	// Fitting data/MC w/ polynomial background (2nd order - fixed) ---
	// --------
	cA->cd(1);
	vector<TH1F*> ho2;
	TF1 *f2 = fit::fit_2hs_1pbg_THI(dx_fit_range,
					sliceY,h_dxHCAL_simu_p,h_dxHCAL_simu_n,3,
					ho2);
	// recording fit params
	std::cout << Form("Pol2 Bg: Run #: %d, R: %f",rnum,f2->GetParameter(1)) << "\n";
	double r2 = f2->GetParameter(1); double r2err = f2->GetParError(1);
	R2.push_back(r2); Rerr2.push_back(r2err);
	// drawing fit histos
	ho2[0]->Draw("E"); util_pd::customize_data(ho2[0]);
	FurtherCustomizeDxHisto(ho2[0],"Pol2 Bg",rnum);
	ho2[4]->Draw("same HIST"); util_pd::customize_psig(ho2[4],0);
	ho2[5]->Draw("same HIST"); util_pd::customize_nsig(ho2[5],0);
	ho2[2]->Draw("same HIST"); util_pd::customize_hbg(ho2[2],0);
	// calculating yields
	std::vector<double> yo2; util_pd::GetYields(f2,ho2[4],ho2[5],ho2[2],yo2);
	outdata << Form("%.4f,%.4f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,",
			r2,r2err,yo2[0],yo2[1],yo2[2],yo2[3],yo2[6],yo2[7]);
	CNpY2.push_back(yo2[0]/charge/daqlt); CNpYerr2.push_back(yo2[1]/charge/daqlt);
	CNnY2.push_back(yo2[2]/charge/daqlt); CNnYerr2.push_back(yo2[3]/charge/daqlt);
	CNtotY2.push_back(yo2[6]/charge/daqlt); CNtotYerr2.push_back(yo2[7]/charge/daqlt);
	// grad the stat box of the fitted histo
	cA->Update();
	TPaveStats *st2 = (TPaveStats*)ho2[0]->FindObject("stats");
	st2->SetOptStat(1); st2->SetOptFit(1111); 
	//
	ho2[0]->Write();
	// ***** -- \\//
	// If 0 field data, then try simple side band fit instead of MC bg
	if (sbsmag==0) {
	  // --------
	  // Fitting data w/ bg estimate from side band fit ---
	  // --------
	  cA->cd(2);
	  vector<TH1F*> ho3;
	  TF1* bg3 = fit::fit_1pbg_SB(dx_fit_range,
				      reject_points,
				      3,//Opoly,
				      fit::GetFitParams(f2),
				      sliceY,
				      ho3);
	  // recording fit params
	  std::cout << Form("Side Band: Run #: %d",rnum) << "\n";
	  R3.push_back(0); Rerr3.push_back(0);
	  // drawing fit histos
	  ho3[0]->Draw("E"); util_pd::customize_data(ho3[0]);
	  FurtherCustomizeDxHisto(ho3[0],"Side Band",rnum);
	  ho3[1]->Draw("same HIST"); util_pd::customize_nsig(ho3[1],0);
	  ho3[2]->Draw("same HIST"); util_pd::customize_hbg(ho3[2],0);
	  // calculating yields
	  std::vector<double> yo3; util_pd::GetYields(ho3[1],ho3[2],yo3);
	  outdata << Form("-99,-99,-99,-99,-99,-99,%.0f,%.0f,",yo3[0],yo3[1]);
	  CNpY3.push_back(0); CNpYerr3.push_back(0);
	  CNnY3.push_back(0); CNnYerr3.push_back(0);
	  CNtotY3.push_back(yo3[0]/charge/daqlt); CNtotYerr3.push_back(yo3[1]/charge/daqlt);
	  // grad the stat box of the fitted histo
	  cA->Update();
	  TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
	  st3->SetOptStat(1); st3->SetOptFit(1111); 
	  //
	  ho3[0]->Write();	  
	}
	else {
	  // --------
	  // Fitting data/MC w/ background from MC ---
	  // --------
	  cA->cd(2);
	  vector<TH1F*> ho3;
	  TF1 *f3 = fit::fit_2hs_2hbg_THI(dx_fit_range,
					  sliceY,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n,
					  ho3);
	  // recording fit params
	  std::cout << Form("MC Bg: Run #: %d, R: %f",rnum,f3->GetParameter(1)) << "\n";
	  double r3 = f3->GetParameter(1); double r3err = f3->GetParError(1);
	  R3.push_back(r3); Rerr3.push_back(r3err);
	  // drawing fit histos
	  ho3[0]->Draw("E"); util_pd::customize_data(ho3[0]);
	  FurtherCustomizeDxHisto(ho3[0],"MC Bg",rnum);
	  ho3[4]->Draw("same HIST"); util_pd::customize_psig(ho3[4],0);
	  ho3[5]->Draw("same HIST"); util_pd::customize_nsig(ho3[5],0);
	  ho3[2]->Draw("same HIST"); util_pd::customize_hbg(ho3[2],0);
	  // calculating yields
	  std::vector<double> yo3; util_pd::GetYields(f3,ho3[4],ho3[5],ho3[2],yo3);
	  outdata << Form("%.4f,%.4f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,",
			  r3,r3err,yo3[0],yo3[1],yo3[2],yo3[3],yo3[6],yo3[7]);
	  CNpY3.push_back(yo3[0]/charge/daqlt); CNpYerr3.push_back(yo3[1]/charge/daqlt);
	  CNnY3.push_back(yo3[2]/charge/daqlt); CNnYerr3.push_back(yo3[3]/charge/daqlt);
	  CNtotY3.push_back(yo3[6]/charge/daqlt); CNtotYerr3.push_back(yo3[7]/charge/daqlt);
	  // grad the stat box of the fitted histo
	  cA->Update();
	  TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
	  st3->SetOptStat(1); st3->SetOptFit(1111); 
	  //
	  ho3[0]->Write();
	  // --------
	  // Fitting data/MC w/ background from Data ---
	  // --------
	  cA->cd(3);
	  vector<TH1F*> ho4;
	  TF1 *f4 = fit::fit_2hs_1hbg_THI(dx_fit_range,
					  sliceY,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_data,
					  ho4);
	  // recording fit params
	  std::cout << Form("Data Bg: Run #: %d, R: %f",rnum,f4->GetParameter(1)) << "\n";
	  double r4 = f4->GetParameter(1); double r4err = f4->GetParError(1);
	  R4.push_back(r4); Rerr4.push_back(r4err);
	  // drawing fit histos
	  ho4[0]->Draw("E"); util_pd::customize_data(ho4[0]);
	  FurtherCustomizeDxHisto(ho4[0],"Data Bg",rnum);
	  ho4[4]->Draw("same HIST"); util_pd::customize_psig(ho4[4],0);
	  ho4[5]->Draw("same HIST"); util_pd::customize_nsig(ho4[5],0);
	  ho4[2]->Draw("same HIST"); util_pd::customize_hbg(ho4[2],0);
	  // calculating yields
	  std::vector<double> yo4; util_pd::GetYields(f4,ho4[4],ho4[5],ho4[2],yo4);
	  outdata << Form("%.4f,%.4f,%.0f,%.0f,%.0f,%.0f,%.0f,%.0f,",
			r4,r4err,yo4[0],yo4[1],yo4[2],yo4[3],yo4[6],yo4[7]);
	  CNpY4.push_back(yo4[0]/charge/daqlt); CNpYerr4.push_back(yo4[1]/charge/daqlt);
	  CNnY4.push_back(yo4[2]/charge/daqlt); CNnYerr4.push_back(yo4[3]/charge/daqlt);
	  CNtotY4.push_back(yo4[6]/charge/daqlt); CNtotYerr4.push_back(yo4[7]/charge/daqlt);
	  // grad the stat box of the fitted histo
	  cA->Update();
	  TPaveStats *st4 = (TPaveStats*)ho4[0]->FindObject("stats");
	  st4->SetOptStat(1); st4->SetOptFit(1111); 
	  //
	  ho4[0]->Write();
	}
	outdata << "\n";
	// ----
	cA->Write();
	if (fRun) { cA->SaveAs(Form("%s[",outPlot.Data())); fRun = false; }
	cA->SaveAs(Form("%s",outPlot.Data()));
      } // if statement
      fRun = false; //set it to false immediately after the 1st good run is found
    }
  } // end of loop over bins ie runs

  // Time to plot and fit yields vs run
  if (is_elastic) {
    TCanvas *cfit2 = FitYields("cfit2",config,"Pol2 Bg",exp_yield,hR_yield,hR_stats,rnums,stats,CNpY2,CNpYerr2);
    cfit2->Update(); cfit2->Write(); cfit2->SaveAs(Form("%s",outPlot.Data()));
    TCanvas *cfit3 = FitYields("cfit3",config,"Side Band",exp_yield,hR_yield,hR_stats,rnums,stats,CNpY3,CNpYerr3);
    cfit3->Update(); cfit3->Write(); cfit3->SaveAs(Form("%s",outPlot.Data()));
    cfit3->SaveAs(Form("%s]",outPlot.Data()));
  }
  else {
    TCanvas *cfit2 = FitYields("cfit2",config,"Pol2 Bg",exp_yield,hR_yield,hR_stats,rnums,stats,R2,Rerr2,
			       CNpY2,CNpYerr2,CNnY2,CNnYerr2,CNtotY2,CNtotYerr2);
    cfit2->Update(); cfit2->Write(); cfit2->SaveAs(Form("%s",outPlot.Data()));
    if (sbsmag==0) {
      TCanvas *cfit3 = FitYields("cfit3",config,"Side Band",exp_yield,hR_yield,hR_stats,rnums,stats,R3,Rerr3,
				 CNpY3,CNpYerr3,CNnY3,CNnYerr3,CNtotY3,CNtotYerr3);
      cfit3->Update(); cfit3->Write(); cfit3->SaveAs(Form("%s",outPlot.Data()));
      cfit3->SaveAs(Form("%s]",outPlot.Data()));
    }
    else {
      TCanvas *cfit3 = FitYields("cfit3",config,"MC Bg",exp_yield,hR_yield,hR_stats,rnums,stats,R3,Rerr3,
				 CNpY3,CNpYerr3,CNnY3,CNnYerr3,CNtotY3,CNtotYerr3);
      cfit3->Update(); cfit3->Write(); cfit3->SaveAs(Form("%s",outPlot.Data()));
      TCanvas *cfit4 = FitYields("cfit4",config,"Data Bg",exp_yield,hR_yield,hR_stats,rnums,stats,R4,Rerr4,
				 CNpY4,CNpYerr4,CNnY4,CNnYerr4,CNtotY4,CNtotYerr4);
      cfit4->Update(); cfit4->Write(); cfit4->SaveAs(Form("%s",outPlot.Data()));
      cfit4->SaveAs(Form("%s]",outPlot.Data()));
    }
  }
  std::cout << totc << " " << gc << "\n";
  std::cout << "Expected normalized yield: " << exp_yield << "\n";

  // reporting output files
  std::cout << "------" << std::endl;
  std::cout << " Fit params  : " << outData << std::endl;
  std::cout << " Summary plots  : " << outPlot << std::endl;
  std::cout << " Output ROOT file  : " << outFile << std::endl;
  std::cout << "------" << std::endl;
}



// calculating charge normalized yields
// std::cout << "---- \n Run: " << rnum << "\n";
// std::cout << "Charge: " << charge << " DAQ LT: " << daqlt << "\n";
// std::cout << "Charge norm. p Yield: " << yo2[0]/charge/daqlt << "\n";
// std::cout << "Charge norm. n Yield: " << yo2[2]/charge/daqlt << "\n----\n";
