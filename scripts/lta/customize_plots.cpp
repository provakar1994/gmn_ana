/*
  This script will contain code to customize plots/canvases for these and possibly
  for publication.
  --------
  P. Datta <pdbforce@jlab.org> 07/01/24
*/

#include <iostream>

#include "gmn_ana.h"

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
Int_t GetNumberOfPads(TCanvas const *canvas) {
  /* Returns # pads in a canvas */
  Int_t numPads = 0;
  TList *primitives = canvas->GetListOfPrimitives();
  TIter next(primitives);
  TObject *obj;
  while ((obj = next())) {
    if (obj->InheritsFrom(TPad::Class())) {
      numPads++;
    }
  }
  std::cout << "# Pads found: " << numPads << "\n";
  return numPads;
}

//______________________________________________________________________________
void GrabAllPads(TCanvas const *c, std::vector<TPad*> & pads) {
  /* Grabs all the pads in a canvas */
  int Npads = GetNumberOfPads(c);
  for (int ipad = 1; ipad <= Npads; ++ipad) {
    TPad *pad = (TPad*)c->GetPad(ipad);
    if (pad) pads.push_back(pad);
  }
}

//______________________________________________________________________________
void customize_title(TPaveText *title) {
  /* Customizes the canvas title */
  title->SetTextSize(0.05);
  title->SetTextFont(62);
}

//______________________________________________________________________________
void customize_axes(TAxis *ax) {
  /* customizes the axes */
  // title
  ax->CenterTitle();
  ax->SetTitleFont(62);
  ax->SetTitleSize(0.05);
  ax->SetTitleOffset(0.8); // gets overwritten in customize_yaxis
  // label
  ax->SetLabelFont(62);
  // ticks
  ax->SetTickLength(0.04);
}

//______________________________________________________________________________
void customize_yaxis(TAxis *ax) {
  /* further customization for the Y axis */
  customize_axes(ax);
  ax->SetMaxDigits(3);
  ax->SetTitleOffset(1);
}

//______________________________________________________________________________
void customize_frame(TFrame *frame) {
  /* customizes the frame */
  frame->SetLineWidth(2);
}

//______________________________________________________________________________
void customize_pad(TPad *p) {
  /* Customizes the following primitives in a pad: 
     TFrame, TPaveText, TH1, TH2, TGraph
   */
  // Add grid lines and xy ticks
  p->SetGridx(); p->SetGridy();
  p->SetTickx(); p->SetTicky();

  // Grabbing list of primitives
  TList *primitives = p->GetListOfPrimitives();
  // iterate through the primitives to do object specific customizations
  TIter next(primitives);
  TObject *obj;
  while ((obj = next())) {
    if (obj->InheritsFrom("TFrame")) {
      TFrame *frame = dynamic_cast<TFrame*>(obj);
      if (frame) {
	customize_frame(frame);
      }
    } else if (obj->InheritsFrom("TPaveText")) {
      TPaveText *title = dynamic_cast<TPaveText*>(obj);
      if (title) {
	customize_title(title);
      }
    } else if (obj->InheritsFrom("TH1")) {
      TH1 *hist = dynamic_cast<TH1*>(obj);
      if (hist) {
	customize_axes(hist->GetXaxis());
	customize_yaxis(hist->GetYaxis());
      }
    } else if (obj->InheritsFrom("TH2")) {
      TH2 *hist = dynamic_cast<TH2*>(obj);
      if (hist) {
	customize_axes(hist->GetXaxis());
	customize_yaxis(hist->GetYaxis());
      }
    } else if (obj->InheritsFrom("TGraph") || obj->InheritsFrom("TGraphWithError")) {
      TGraph *graph = dynamic_cast<TGraph*>(obj);
      if (graph) {
	customize_axes(graph->GetXaxis());
	customize_yaxis(graph->GetYaxis());
      }
    }
  }
}

//______________________________________________________________________________
void customize_canvas(TCanvas *c) {
  /* Customizes all pads within a given canvas */
  // grab all the pads
  std::vector<TPad*> pads;
  GrabAllPads(c,pads);
  if (pads.size()>0) {
    // if multiple pads available then
    // customize them one-by-one
    for (auto & pad : pads) {
      customize_pad(pad);
    }
  } else {
    // otherwise, grab the main canvas
    TPad *pad = (TPad*)c->GetPad(0);
    customize_pad(pad);
  }
}

//______________________________________________________________________________
int customize_plots() {

  //TString inFile = "../pdout/fits/sbs4sbs0p/0p65zoff_yield_per_run_qelas_pass2_simc_sbs4_sbs0p_model2.root";
  TString inFile = "../pdout/fits/sbs4sbs0p/0p65zoff_fit_dx_qelas_pass2_simc_sbs4_sbs0p_model2.root";
  
  // reading input file
  TFile * file = ReadRootFile(inFile);

  // read a canvas
  TCanvas *c = (TCanvas*)file->Get("c1");
  customize_canvas(c);
  c->Draw();

  return 0;
}
