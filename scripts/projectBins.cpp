/* This macro will be used to extract yields per run 
   It will be configures by con_fit_dx.json. Will look
   for prodection keys ie qeals_prod or elas_prod. 
   -------
   P. Datta Created 06-19-2024
*/

#include <unordered_map>

#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"
#include "TLatex.h"

#include "gmn_ana.h"

void gStyleFitCanvas() 
{
  gStyle->SetOptStat("e"); gStyle->SetOptFit(1); 
  gStyle->SetErrorX(0);
}

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

TCanvas *FitYields( std::string const &canvname, std::string const &bgshape,
		    std::vector<double> const &rnums, std::vector<double> const &stats,
		    std::vector<double> const &R_fit, std::vector<double> const &R_fit_err,
		    std::vector<double> const &CNpY, std::vector<double> const &CNpYerr,
		    std::vector<double> const &CNnY, std::vector<double> const &CNnYerr) {

  TCanvas *cfit = new TCanvas(canvname.c_str(),bgshape.c_str(), 1200, 1000);
  cfit->Divide(2,2);
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
  return cfit;
}


void projectBins() {
    
  TString filename = "/lustre19/expphy/volatile/halla/sbs/pdbforce/gmn_ana/pdout/fits/sbs4sbs30p/0p65zoff_0p39sf_fit_dx_qelas_pass2_simc_sbs4_sbs30p_model2.root";
  TFile * file = ReadRootFile(filename);

  int nruns = -1;
  int conf = 4;
  std::string target = "LD2";
  int pass = 2;
  int sbsmag = 30;
  int verbosefn = 0;
  
  // Reading run list
  std::vector<CodaRun> cruns; util_pd::ReadRunList("../DB",nruns,conf,target,pass,sbsmag,verbosefn,cruns);
  // implementing hashtable w/ run number as key, for efficiency
  std::unordered_map<int,CodaRun> mrun;
  for (auto & crun: cruns) mrun[crun.runnum] = crun;
  
  // Defining histos to fit
  TH1F * h_dxHCAL_simu_p = (TH1F*)file->Get("h_dxHCAL_simu_p");
  TH1F * h_dxHCAL_simu_n = (TH1F*)file->Get("h_dxHCAL_simu_n");
  TH1F * h_dxHCAL_bg_data = (TH1F*)file->Get("h_dxHCAL_bg_data");
  TH1F * h_dxHCAL_bg_inel_p = (TH1F*)file->Get("h_dxHCAL_bg_inel_p");
  TH1F * h_dxHCAL_bg_inel_n = (TH1F*)file->Get("h_dxHCAL_bg_inel_n");
  TH1F * h_dxHCAL_bg_inel = (TH1F*)file->Get("h_dxHCAL_bg_inel");
  // vs Rnum histo
  TH2F * h_dxHCAL_vs_rnum = (TH2F*)file->Get("h2_dxHCAL_vs_rnum");

  // Define arrays to store fit results
  int nBinsX = h_dxHCAL_vs_rnum->GetNbinsX();
  std::vector<double> rnums, stats;
  std::vector<double> R2, Rerr2, CNpY2, CNpYerr2, CNnY2, CNnYerr2; 
  std::vector<double> R3, Rerr3, CNpY3, CNpYerr3, CNnY3, CNnYerr3;

  std::vector<double> dx_fit_range{-2,1};

  int totc = 0; int gc = 0;

  TFile *fout = new TFile("projectBins.root", "RECREATE");

  // Loop over each x-bin
  for (int i = 1; i <= nBinsX; ++i) {
    totc++;

    if (h_dxHCAL_vs_rnum->GetBinContent(i)!=0) {
      gc++;

      int rnum = (int)h_dxHCAL_vs_rnum->GetXaxis()->GetBinCenter(i);
      rnums.push_back(rnum);

      // getting the charge and DAQ livetime for the run
      double charge = mrun[rnum].charge;
      double daqlt = mrun[rnum].daqlvtm;

      // Project slice along y-axis
      TH1F *sliceY = (TH1F*)h_dxHCAL_vs_rnum->ProjectionY(Form("hfit_run%d",rnum), i, i);
      sliceY->SetTitle(Form("%d",rnum));
      stats.push_back(sliceY->GetEntries());

      TCanvas *cA = util_pd::TC("cA",1,2); gStyleFitCanvas();
      
      // Fitting data/MC w/ background from MC ---
      cA->cd(1);
      vector<TH1F*> ho2;
      TF1 *f2 = fit::fit_2hs_2hbg_THI(dx_fit_range,
				      sliceY,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n,
				      ho2);
      // recording fit params
      std::cout << Form("MC Bg: Run #: %d, R: %f",rnum,f2->GetParameter(1)) << "\n";
      R2.push_back(f2->GetParameter(1)); Rerr2.push_back(f2->GetParError(1));
      // converting fit fn to a hostogram
      TH1F *hf2 = (TH1F*)sliceY->Clone(); util_pd::TF1toTH1F(f2,hf2);
      ho2[0]->SetTitle(Form("MC Bg | Run: %d",rnum));
      ho2[0]->GetXaxis()->SetTitle("#Deltax (m)");
      ho2[0]->Draw("E"); util_pd::customize_data(ho2[0]);
      ho2[4]->Draw("same HIST"); util_pd::customize_psig(ho2[4],0);
      ho2[5]->Draw("same HIST"); util_pd::customize_nsig(ho2[5],0);
      ho2[2]->Draw("same HIST"); util_pd::customize_hbg(ho2[2],0);
      // calculating yields
      std::vector<double> yo2; util_pd::GetYields(f2,ho2[4],ho2[5],ho2[2],yo2);
      // calculating charge normalized yields
      // std::cout << "---- \n Run: " << rnum << "\n";
      // std::cout << "Charge: " << charge << " DAQ LT: " << daqlt << "\n";
      // std::cout << "Charge norm. p Yield: " << yo2[0]/charge/daqlt << "\n";
      // std::cout << "Charge norm. n Yield: " << yo2[2]/charge/daqlt << "\n----\n";
      CNpY2.push_back(yo2[0]/charge/daqlt); CNpYerr2.push_back(yo2[1]/charge/daqlt);
      CNnY2.push_back(yo2[2]/charge/daqlt); CNnYerr2.push_back(yo2[3]/charge/daqlt);
      // grad the stat box of the fitted histo
      cA->Update();
      TPaveStats *st2 = (TPaveStats*)ho2[0]->FindObject("stats");
      st2->SetOptStat(1);
      st2->SetOptFit(1111); 
      //
      ho2[0]->Write();
      // ----

      // Fitting data/MC w/ polynomial background (2nd order - fixed) ---
      cA->cd(2);
      vector<TH1F*> ho3;
      TF1 *f3 = fit::fit_2hs_1pbg_THI(dx_fit_range,
				      sliceY,h_dxHCAL_simu_p,h_dxHCAL_simu_n,2,
				      ho3);
      // recording fit params
      std::cout << Form("Pol2 Bg: Run #: %d, R: %f",rnum,f3->GetParameter(1)) << "\n";
      R3.push_back(f3->GetParameter(1)); Rerr3.push_back(f3->GetParError(1));
      // converting fit fn to a hostogram
      TH1F *hf3 = (TH1F*)sliceY->Clone(); util_pd::TF1toTH1F(f3,hf3);
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
      // grad the stat box of the fitted histo
      cA->Update();
      TPaveStats *st3 = (TPaveStats*)ho3[0]->FindObject("stats");
      st3->SetOptStat(1);
      st3->SetOptFit(1111); 
      //
      ho3[0]->Write();
      // ----

      cA->Write();
    }

  }


  TCanvas *cfit2 = FitYields("cfit2","MC Bg",rnums,stats,R2,Rerr2,CNpY2,CNpYerr2,CNnY2,CNnYerr2);
  cfit2->Write();

  TCanvas *cfit3 = FitYields("cfit3","Pol2 Bg",rnums,stats,R3,Rerr3,CNpY3,CNpYerr3,CNnY3,CNnYerr3);
  cfit3->Write();
  
  std::cout << totc << " " << gc << "\n"; 
}


