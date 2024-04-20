#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TString.h"

#include "gmn_ana.h"

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

   // Create stat box with fit parameters
    TPaveText* statBox = new TPaveText(0.6, 0.8, 0.9, 0.9, "NDC");
    statBox->AddText(Form("C: %.4f #pm %.4f", constantValue, constantError));
    // statBox->AddText(Form("Error: %.2f", constantError));
    statBox->SetFillColor(0);  // Set background color to white
    statBox->SetBorderSize(1);  // Set border size
    statBox->Draw();
    
    std::cout << "Fit Result:" << std::endl;
    std::cout << "Constant Value: " << constantValue << " +/- " << constantError << std::endl;

    // Optionally, draw the graph and fit function on a canvas
    TCanvas* canvas = new TCanvas("canvas", "Fit Canvas", 800, 600);
    graph->Draw("AP");  // Draw graph with axes and points
    constantFit->Draw("SAME");  // Draw fit function on top of graph
    statBox->Draw("SAME");
    canvas->Update();
}

void DrawTGraphWithErrors(const std::vector<double>& xValues, const std::vector<double>& yValues, const std::vector<double>& yErrors) {
    // Check if the sizes of the arrays match
    if (xValues.size() != yValues.size() || xValues.size() != yErrors.size()) {
        std::cerr << "Error: Array sizes do not match!" << std::endl;
        return;
    }

    // Create a TGraphErrors object
    TGraphErrors *graph = new TGraphErrors(xValues.size(), &xValues[0], &yValues[0], 0, &yErrors[0]);

    // Set graph attributes
    graph->SetTitle("Stability of R_fit Across Runs");
    graph->GetXaxis()->SetTitle("Run Number");
    graph->GetYaxis()->SetTitle("R_fit");

    // Create a canvas and draw the graph
    TCanvas *canvas = new TCanvas("canvas", "TGraphErrors Example", 800, 600);
    graph->Draw("AP");  // A: Draw axis, P: Draw points

    // Optionally, add a legend or additional customization
    graph->SetMarkerStyle(20);

    // Delete the graph and canvas when done to free memory
    // delete graph;
    // delete canvas;

    // Fit graph
    FitGraphWithConstant(graph);
}

void projectBins() {

  TString filename = "/lustre19/expphy/volatile/halla/sbs/pdbforce/gmn_ana/pdout/fits/sbs14sbs70p/0p674sf_fit_dx_qelas_pass2_simc_sbs14_sbs70p_model2.root";

  TFile * file = ReadRootFile(filename);

  // Defining histos to fit
  TH1F * h_dxHCAL_simu_p = (TH1F*)file->Get("h_dxHCAL_simu_p");
  TH1F * h_dxHCAL_simu_n = (TH1F*)file->Get("h_dxHCAL_simu_n");
  TH1F * h_dxHCAL_bg_data = (TH1F*)file->Get("h_dxHCAL_bg_data");
  TH1F * h_dxHCAL_bg_inel_p = (TH1F*)file->Get("h_dxHCAL_bg_inel_p");
  TH1F * h_dxHCAL_bg_inel_n = (TH1F*)file->Get("h_dxHCAL_bg_inel_n");
  TH1F * h_dxHCAL_bg_inel = (TH1F*)file->Get("h_dxHCAL_bg_inel");
  // vs Rnum histo
  TH2F * h_dxHCAL_vs_rnum = (TH2F*)file->Get("h2_dxHCAL_vs_rnum");

  // Create a canvas to draw the histograms
  TCanvas *canvas = new TCanvas("canvas", "Projection per Bin Example", 800, 600);

  // Define arrays to store fit results
  int nBinsX = h_dxHCAL_vs_rnum->GetNbinsX();
  std::vector<double> rnums;
  std::vector<double> R, Rerr;

  std::vector<double> dx_fit_range{-1.5,0.5};

  int totc = 0; int gc = 0;

  // Loop over each x-bin
  for (int i = 1; i <= nBinsX; ++i) {
    totc++;

    if (h_dxHCAL_vs_rnum->GetBinContent(i)!=0) {
      gc++;

      int rnum = (int)h_dxHCAL_vs_rnum->GetXaxis()->GetBinCenter(i);
      rnums.push_back(rnum);

      // Project slice along y-axis
      TH1F *sliceY = (TH1F*)h_dxHCAL_vs_rnum->ProjectionY(Form("hfit_run%d",rnum), i, i);
      sliceY->SetTitle(Form("%d",rnum));

      // fitting
      vector<TH1F*> ho2;
      TF1 *f2 = fit::fit_2hs_2hbg_THI(dx_fit_range,
				      sliceY,h_dxHCAL_simu_p,h_dxHCAL_simu_n,h_dxHCAL_bg_inel_p,h_dxHCAL_bg_inel_n,
				      ho2);
      // recording fit params
      std::cout << Form("Run #: %d, R: %f",rnum,f2->GetParameter(1)) << "\n";
      R.push_back(f2->GetParameter(1));
      Rerr.push_back(f2->GetParError(1));
      // converting fit fn to a hostogram
      TH1F *hf2 = (TH1F*)sliceY->Clone(); util_pd::TF1toTH1F(f2,hf2);
      ho2[0]->Draw(); 

      // sliceY->Draw("E");
      canvas->Update();
    }
  }

  DrawTGraphWithErrors(rnums,R,Rerr);
  std::cout << totc << " " << gc << "\n"; 

}
