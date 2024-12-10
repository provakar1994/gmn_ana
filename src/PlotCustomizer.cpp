#include "PlotCustomizer.h"

PlotCustomizer::PlotCustomizer(bool gridON, bool statON, bool titleON) : fgridON(gridON), fstatON(statON), ftitleON(titleON) {}
//______________________________________________________________________________
void PlotCustomizer::customize_pvtext(TPaveText *pt) const {
  pt->SetFillColor(0); // Transparent fill
  pt->SetTextAlign(12);
  pt->SetTextFont(kFont);     // Helvetica font
  pt->SetTextSize(0.05);   // Text size
  pt->SetShadowColor(kGray+1);   // Shadow color
  //pt->SetLineWidth(2);     // Border width
}
//______________________________________________________________________________
void PlotCustomizer::AddText(TString const &text,
			     double x1NDC, double y1NDC,
			     double x2NDC, double y2NDC) const {
  /* Adds a TPaveText in the place of the title */
  TPaveText *pt = new TPaveText(x1NDC, y1NDC, x2NDC, y2NDC, "NDC");
  pt->AddText(text);
  customize_pvtext(pt);
  pt->Draw();
}
//______________________________________________________________________________
void PlotCustomizer::AddTitleText(TString const &text, double x2NDC) const {
  /* Adds a TPaveText in the place of the title */
  AddText(text,0.15,0.92,x2NDC,0.98);
}
//______________________________________________________________________________
Int_t PlotCustomizer::GetNumberOfPads(TCanvas const *canvas) const {
  /* Returns # pads in a canvas */
  if (canvas) {
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
  } else {
    std::cout << "ERROR!! Canvas doesn't exist!!\n";
    return -1;
  }
}
//______________________________________________________________________________
void PlotCustomizer::GrabAllPads(TCanvas const *c, std::vector<TPad*> &pads) const {
  /* Grabs all the pads in a canvas */
  int Npads = GetNumberOfPads(c);
  for (int ipad = 1; ipad <= Npads; ++ipad) {
    TPad *pad = (TPad*)c->GetPad(ipad);
    if (pad) pads.push_back(pad);
  }
}
//______________________________________________________________________________
void PlotCustomizer::customize_TH1(TH1 *h) const {
  /* Customizes 1D histos */
  if (!ftitleON) {h->SetTitle(""); h->SetTitleFont(kFont);}
  if (!fstatON) h->SetStats(0);
  h->SetLineWidth(2);
  //h->SetLineColor(kBlack);
  customize_axes(h->GetXaxis());
  customize_axes(h->GetYaxis());
}
//______________________________________________________________________________
void PlotCustomizer::customize_TH2(TH2 *h) const {
  /* Customizes 1D histos */
  if (!ftitleON) {h->SetTitle(""); h->SetTitleFont(kFont);}
  if (!fstatON) h->SetStats(0);
  customize_axes(h->GetXaxis());
  customize_axes(h->GetYaxis());
}
//______________________________________________________________________________
void PlotCustomizer::customize_TGraph(TGraph *g) const {
  /* Customizes 1D histos */
  if (!ftitleON) g->SetTitle("");
  if (!fstatON) g->SetStats(0);
  customize_axes(g->GetXaxis());
  customize_axes(g->GetYaxis());
}
//______________________________________________________________________________
void PlotCustomizer::customize_title(TPaveText *title) const {
  /* Customizes the canvas title */
  title->SetTextSize(0.05);
  title->SetTextFont(kFont);
  //title->Clear();
}
//______________________________________________________________________________
void PlotCustomizer::customize_stats(TPaveText *stats) const {
  /* Customizes the canvas stats */
  stats->Clear();
}
//______________________________________________________________________________
void PlotCustomizer::customize_axes(TAxis *ax) const {
  /* customizes the axes */
  ax->CenterTitle();
  ax->SetTitleFont(kFont);
  ax->SetTitleSize(0.05);
  ax->SetTitleOffset(1.2);
  ax->SetLabelFont(kFont);
  ax->SetLabelSize(kLabelSize);
  ax->SetTickLength(0.04);
}
//______________________________________________________________________________
void PlotCustomizer::customize_yaxis(TAxis *ax) const {
  /* further customization for the Y axis */
  customize_axes(ax);
  ax->SetMaxDigits(3);
  ax->SetTitleOffset(1);
}
//______________________________________________________________________________
void PlotCustomizer::customize_frame(TFrame *frame) const {
  /* customizes the frame */
  frame->SetLineWidth(2);
}
//______________________________________________________________________________
void PlotCustomizer::customize_palette(TPaletteAxis *pal) const {
  /* customizes color palette */
  pal->SetLabelFont(kFont);
  pal->SetLabelSize(kLabelSize);
  pal->SetNdivisions(6);
  pal->SetLineWidth(2);
  // set maxdigits to 3  
  pal->SetMaxDigits(3);
}
//______________________________________________________________________________
void PlotCustomizer::customize_margin(TPad *p) const {
  p->SetLeftMargin(0.15);
  p->SetRightMargin(0.05);
  p->SetBottomMargin(0.15);
  p->SetTopMargin(0.10);

  // Print the list of primitives
  TList *primitives = p->GetListOfPrimitives();
  TIter next(primitives);
  TObject *obj;
  //std::cout << "Primitives in the pad:" << std::endl;
  while ((obj = next())) {
    if (obj->InheritsFrom("TH2")) {
      TH2 *hist = dynamic_cast<TH2*>(obj);
      if (hist) {
	TPaletteAxis *palette = (TPaletteAxis*)hist->GetListOfFunctions()->FindObject("palette");
	if (palette) {
	  std::cout << "TPaletteAxis found, adjusting right margin\n";
	  double paletteWidth = palette->GetX2NDC() - palette->GetX1NDC();
	  p->SetRightMargin(0.062 + paletteWidth);
	  // Adjust the X1NDC of palette accordingly
	  palette->SetX1NDC(1.0-p->GetRightMargin()+0.003);
	  palette->SetX2NDC(palette->GetX1NDC()+paletteWidth);

	  // Adjust the Y1NDC and Y2NDC of the palette to match the pad
	  palette->SetY1NDC(p->GetBottomMargin());
	  palette->SetY2NDC(1.0 - p->GetTopMargin());

	  // further customization
	  customize_palette(palette);
	}
      }
    }
  }
}
//______________________________________________________________________________
void PlotCustomizer::customize_pad(TPad *p) const {
  /* Customizes the following primitives in a pad: 
     TFrame, TPaveText, TH1, TH2, TGraph
  */
  p->cd();
  p->Update();
  //if (!fstatON) gStyle->SetOptStat(0);

  if (fgridON) {
    p->SetGridx();
    p->SetGridy();
  }
    
  p->SetTickx();
  p->SetTicky();
  customize_margin(p);

  TList *primitives = p->GetListOfPrimitives();
  TIter next(primitives);
  TObject *obj;
  while ((obj = next())) {    
    //std::cout << "  " << obj->ClassName() << " - " << obj->GetName() << "\n";
    if (obj->InheritsFrom("TFrame")) {
      TFrame *frame = dynamic_cast<TFrame*>(obj);
      if (frame) customize_frame(frame);
      p->Update();
    } // else if (obj->InheritsFrom("TPaveText")) {
    //   TPaveText *title = dynamic_cast<TPaveText*>(obj);
    //   if (title) customize_title(title);
    // }
    else if (obj->InheritsFrom("TPaveStats")) {
      TPaveStats *stats = dynamic_cast<TPaveStats*>(obj);
      if (stats) {
	customize_stats(stats);
	p->Update();
      }
    } else if (obj->InheritsFrom("TH1")) {
      TH1 *hist = dynamic_cast<TH1*>(obj);
      if (hist) {
	customize_TH1(hist);
	p->Update();
      }
    } else if (obj->InheritsFrom("TH2")) {
      TH2 *hist = dynamic_cast<TH2*>(obj);
      if (hist) {
	customize_TH2(hist);
	p->Update();
      }
    } else if (obj->InheritsFrom("TGraph") || obj->InheritsFrom("TGraphWithError")) {
      TGraph *graph = dynamic_cast<TGraph*>(obj);
      if (graph) {
	customize_TGraph(graph);
	p->Update();
      }
    }
  }
  
  p->Update();
}
//______________________________________________________________________________
void PlotCustomizer::customize_canvas(TCanvas *c) const {
  /* Customizes all pads within a given canvas */
  if (c) {
    std::vector<TPad*> pads;
    GrabAllPads(c, pads);
    if (!pads.empty()) {
      for (auto &pad : pads) {
	customize_pad(pad);
	if (fgridON) pad->RedrawAxis("g");
      }
    } else {
      TPad *pad = (TPad*)c->GetPad(0);
      if (pad) customize_pad(pad);
      pad->Update();
      if (fgridON) {
	pad->SetGridx();
	pad->SetGridy();
      }      
    }
    if (fgridON) c->RedrawAxis("g");
    
  } else {
    std::cout << "ERROR!! Canvas doesn't exist!!\n";
  }
}
