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

class PlotCustomizer {
 public:
  PlotCustomizer(bool gridON = true);

  Int_t GetNumberOfPads(TCanvas const *canvas) const;
  void GrabAllPads(TCanvas const *c, std::vector<TPad*> &pads) const;
  void customize_title(TPaveText *title) const;
  void customize_stats(TPaveText *stats) const;
  void customize_axes(TAxis *ax) const;
  void customize_yaxis(TAxis *ax) const;
  void customize_frame(TFrame *frame) const;
  void customize_margin(TPad *p) const;
  void customize_pad(TPad *p) const;
  void customize_canvas(TCanvas *c) const;

 private:
  bool fgridON; // Attribute to control grid lines
};

#endif // PLOT_CUSTOMIZER_H
