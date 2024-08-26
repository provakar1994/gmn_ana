/*
  Class to customize plots for thesis.
  ---------
  P. Datta CREATED 07/09/2024
*/

#ifndef PLOT_CUSTOMIZER_H
#define PLOT_CUSTOMIZER_H

#include <iostream>
#include <vector>
#include <TCanvas.h>
#include <TPad.h>
#include <TPaveText.h>
#include <TPaveStats.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TList.h>
#include <TObject.h>
#include <TAxis.h>
#include <TFrame.h>
#include <TStyle.h>
#include <TPaletteAxis.h> 

class PlotCustomizer {
 public:
  PlotCustomizer(bool gridON = true, bool statON = false);

  Int_t GetNumberOfPads(TCanvas const *canvas) const;
  void GrabAllPads(TCanvas const *c, std::vector<TPad*> &pads) const;
  void customize_TH1(TH1 *h) const;
  void customize_TH2(TH2 *h) const;
  void customize_TGraph(TGraph *g) const;
  void customize_title(TPaveText *title) const;
  void customize_stats(TPaveText *stats) const;
  void customize_axes(TAxis *ax) const;
  void customize_yaxis(TAxis *ax) const;
  void customize_frame(TFrame *frame) const;
  void customize_margin(TPad *p) const;
  void customize_palette(TPaletteAxis *pal) const;
  void customize_pad(TPad *p) const;
  void customize_canvas(TCanvas *c) const;
  void customize_pvtext(TPaveText *pt) const;

  void AddText(TString const &text,
	       double x1NDC, double y1NDC,
	       double x2NDC, double y2NDC) const;
  void AddTitleText(TString const &text, double x2NDC) const;

 private:
  bool fgridON; // Attribute to control grid lines
  bool fstatON; // Attribute to control statbox lines
  //
  int kFont = 62;
  double kLabelSize = 0.04;
};

#endif // PLOT_CUSTOMIZER_H
