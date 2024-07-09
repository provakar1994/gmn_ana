#include "PlotCustomizer.h"

PlotCustomizer::PlotCustomizer(bool gridON) : fgridON(gridON) {}

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
void PlotCustomizer::customize_title(TPaveText *title) const {
    /* Customizes the canvas title */
    title->SetTextSize(0.05);
    title->SetTextFont(62);
    title->Clear();
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
    ax->SetTitleFont(62);
    ax->SetTitleSize(0.05);
    ax->SetTitleOffset(1.2);
    ax->SetLabelFont(62);
    ax->SetLabelSize(0.04);
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
void PlotCustomizer::customize_margin(TPad *p) const {
    p->SetLeftMargin(0.15);
    p->SetRightMargin(0.05);
    p->SetBottomMargin(0.15);
    p->SetTopMargin(0.10);
}
//______________________________________________________________________________
void PlotCustomizer::customize_pad(TPad *p) const {
    /* Customizes the following primitives in a pad: 
       TFrame, TPaveText, TH1, TH2, TGraph
    */
    p->cd();
    gStyle->SetOptStat(0);

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
        if (obj->InheritsFrom("TFrame")) {
            TFrame *frame = dynamic_cast<TFrame*>(obj);
            if (frame) customize_frame(frame);
        } else if (obj->InheritsFrom("TPaveText")) {
            TPaveText *title = dynamic_cast<TPaveText*>(obj);
            if (title) customize_title(title);
        } else if (obj->InheritsFrom("TPaveStats")) {
            TPaveStats *stats = dynamic_cast<TPaveStats*>(obj);
            if (stats) customize_stats(stats);
        } else if (obj->InheritsFrom("TH1")) {
            TH1 *hist = dynamic_cast<TH1*>(obj);
            if (hist) {
                hist->SetStats(0);
                customize_axes(hist->GetXaxis());
                customize_axes(hist->GetYaxis());
            }
        } else if (obj->InheritsFrom("TH2")) {
            TH2 *hist = dynamic_cast<TH2*>(obj);
            if (hist) {
                hist->SetStats(0);
                customize_axes(hist->GetXaxis());
                customize_axes(hist->GetYaxis());
            }
        } else if (obj->InheritsFrom("TGraph") || obj->InheritsFrom("TGraphWithError")) {
            TGraph *graph = dynamic_cast<TGraph*>(obj);
            if (graph) {
                graph->SetStats(0);
                customize_axes(graph->GetXaxis());
                customize_axes(graph->GetYaxis());
            }
        }
    }
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
            }
        } else {
            TPad *pad = (TPad*)c->GetPad(0);
            if (pad) customize_pad(pad);
        }
    } else {
        std::cout << "ERROR!! Canvas doesn't exist!!\n";
    }
}
