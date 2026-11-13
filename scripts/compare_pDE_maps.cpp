/*
  compare_pDE_maps.cpp
  --------------------
  Compare two or more HCAL proton detection efficiency (pDE) maps made by pDE_map.cpp.

  Why: the map h2_effi_map is used to weight all MC (weight_effic = weight * map(x,y) / const,
  see scripts/simulation/qelas_ana_simu.cpp). Maps built from different data subsets (e.g.
  different SBS field settings, or different cuts) should agree within statistics if the
  subset choice does not matter. This macro quantifies that.

  Inputs: ROOT files written by pDE_map.cpp. Objects used:
    h2_effi_map      : efficiency map, x-axis = yHCAL_exp, y-axis = xHCAL_exp - p_def (m)
    h2_xyexp_earm    : "all" counts (e-arm elastic cuts), same binning as h2_effi_map
    h2_xyexp_eNharm  : "pass" counts (+ HCAL hit cut), same binning as h2_effi_map
    h_xexp_pDE, h_yexp_pDE : 1D efficiencies vs x (with y cut) and vs y (with x cut)
    h2_xyexp_all       : "all" events vs (y, x - p_def) = proton hit position, fine binning
    h2_xyexp_all_nodef : "all" events vs (y, x)         = neutron hit position, fine binning

  Outputs (to <outBase>.pdf and <outBase>.csv):
    page 1   : the efficiency maps
    per pair : pull map (A-B)/sigma, pull distribution (mean/RMS), chi2/ndf of the difference,
               ratio map A/B
    1D page  : efficiency vs x and vs y overlays with pol0 fits, and ratios to the first map
    impact   : acceptance-weighted average of each map at the proton and at the neutron
               hit positions (proxy for the mean MC weight_effic), see impact_on_mc() below
    CSV      : long format "quantity,label,value,error"

  Usage (from scripts/):
    root -l -b -q 'compare_pDE_maps.cpp("a.root,b.root,c.root","A,B,C","pdout/pDE/out", "0:2", 0)'
  Arguments:
    fileList    : comma-separated map files
    labelList   : comma-separated labels (same order)
    outBase     : output path without extension
    nestedPairs : optional, comma-separated "i:j" pairs meaning the event sample of map i
                  is a SUBSET of the sample of map j (e.g. 70p-only vs all fields). Their
                  statistical errors are then correlated and the error of the difference is
                  sigma^2 = |sigma_i^2 - sigma_j^2| instead of sigma_i^2 + sigma_j^2.
    weightIndex : which file provides the acceptance (h2_xyexp_all/_nodef) for the impact
                  estimate (default 0 = first file)
*/

#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TObjArray.h>
#include <TObjString.h>
#include <TLatex.h>
#include <TLine.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <cmath>

// ---------------------------------------------------------------------------
// Named constants
// ---------------------------------------------------------------------------
// pol0 fit ranges used by pDE_map.cpp (conf_pDE_map_pass3.json: h_xexp_fitR, h_yexp_fitR).
// Also used as the "central region" of the 2D map (x = xHCAL_exp - p_def, y = yHCAL_exp).
static const double kFitXmin = -2.2, kFitXmax = 1.0;
static const double kFitYmin = -0.5, kFitYmax = 0.5;
// A 2D bin enters the pull/chi2 only if both maps have at least this many "all" counts;
// below that the binomial (Gaussian) error approximation is unreliable.
static const double kMinAllCounts = 50.;
// Pull map colour range (in sigma).
static const double kPullRange = 5.;
// Ratio map colour range (A/B).
static const double kRatioMin = 0.9, kRatioMax = 1.1;
// Content given to unused pull-map bins (below the colour range -> not drawn).
static const double kEmptyBin = -999.;

// One loaded map with everything we need from its file.
struct PDEMap {
  TString label;
  TH2 *effi = nullptr;      // h2_effi_map
  TH2 *nAll = nullptr;      // h2_xyexp_earm   (may be null for old files)
  TH2 *nPass = nullptr;     // h2_xyexp_eNharm (may be null for old files)
  TH1 *effX = nullptr;      // h_xexp_pDE
  TH1 *effY = nullptr;      // h_yexp_pDE
  TH2 *accP = nullptr;      // h2_xyexp_all       (proton positions)
  TH2 *accN = nullptr;      // h2_xyexp_all_nodef (neutron positions)
};

// Split a comma-separated string into a vector of trimmed strings.
std::vector<TString> split_csv(const char *list) {
  std::vector<TString> out;
  TObjArray *tok = TString(list).Tokenize(",");
  for (int i = 0; i < tok->GetEntries(); i++) {
    TString s = ((TObjString*)tok->At(i))->GetString();
    s = s.Strip(TString::kBoth);
    if (s.Length() > 0) out.push_back(s);
  }
  delete tok;
  return out;
}

// Get a histogram from a file and detach it from the file (returns null if missing).
template <class T> T *get_hist(TFile *f, const char *name, const TString &tag) {
  T *h = dynamic_cast<T*>(f->Get(name));
  if (!h) return nullptr;
  h = (T*)h->Clone(Form("%s_%s", name, tag.Data()));
  h->SetDirectory(nullptr);
  h->GetListOfFunctions()->Delete();  // drop fits/stat boxes stored by pDE_map.cpp
  h->SetStats(0);
  return h;
}

// Statistical variance of the efficiency in one 2D bin.
// Uses the pass/all counts when available: var = e(1-e)/N with the Laplace estimate
// e = (pass+1)/(all+2), so that bins with e = 0 or 1 still get a non-zero error.
// Falls back to the error stored in the map (binomial, from pDE_map.cpp) otherwise.
double bin_variance(const PDEMap &m, int bin) {
  if (m.nAll && m.nPass) {
    double nAll = m.nAll->GetBinContent(bin);
    if (nAll <= 0) return -1.;
    double e = (m.nPass->GetBinContent(bin) + 1.) / (nAll + 2.);
    return e * (1. - e) / nAll;
  }
  double err = m.effi->GetBinError(bin);
  return err * err;
}

// Number of "all" counts in a bin (huge if unknown, so the count threshold is not applied).
double bin_all_counts(const PDEMap &m, int bin) {
  return m.nAll ? m.nAll->GetBinContent(bin) : 1e30;
}

// True if the bin centre (y = HCAL horizontal, x = HCAL vertical) lies in the fit region.
bool in_fit_region(double ycen, double xcen) {
  return xcen > kFitXmin && xcen < kFitXmax && ycen > kFitYmin && ycen < kFitYmax;
}

// pol0 fit of a 1D efficiency in [xmin,xmax]; returns value and error.
std::pair<double,double> fit_pol0(TH1 *h, double xmin, double xmax) {
  if (!h) return {NAN, NAN};
  TF1 f(Form("fpol0_%s", h->GetName()), "pol0", xmin, xmax);
  h->Fit(&f, "QR0N");
  return {f.GetParameter(0), f.GetParError(0)};
}

// Overall efficiency pass/all (optionally only inside the fit region) with binomial error.
std::pair<double,double> overall_eff(const PDEMap &m, bool fitRegionOnly) {
  if (!m.nAll || !m.nPass) return {NAN, NAN};
  double nPass = 0., nAll = 0.;
  for (int ix = 1; ix <= m.nAll->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= m.nAll->GetNbinsY(); iy++) {
      double ycen = m.nAll->GetXaxis()->GetBinCenter(ix);
      double xcen = m.nAll->GetYaxis()->GetBinCenter(iy);
      if (fitRegionOnly && !in_fit_region(ycen, xcen)) continue;
      nAll  += m.nAll->GetBinContent(ix, iy);
      nPass += m.nPass->GetBinContent(ix, iy);
    }
  }
  double e = nAll > 0 ? nPass / nAll : NAN;
  return {e, nAll > 0 ? std::sqrt(e * (1. - e) / nAll) : NAN};
}

// True if every map has at least kMinAllCounts "all" events in the bin containing (y, x),
// i.e. all maps are measured there (a map with no data at a position returns 0 in the MC).
bool covered_by_all(const std::vector<PDEMap> &maps, double ycen, double xcen) {
  for (const PDEMap &m : maps) {
    if (!m.nAll) continue;  // unknown counts: cannot judge, accept
    if (m.nAll->GetBinContent(m.nAll->FindBin(ycen, xcen)) < kMinAllCounts) return false;
  }
  return true;
}

// Acceptance-weighted average of a map: sum_w acc(bin) * map(bin centre) / sum_w acc.
// acc is a fine-binned (y, x) distribution of hit positions; the map is looked up at each
// acc bin centre exactly as the MC does (FindBin + GetBinContent, no interpolation).
// If 'common' is given, only positions measured by all of those maps are used, so that maps
// built from different field settings (different HCAL coverage) are compared on equal footing.
double weighted_map_average(TH2 *map, TH2 *acc, const std::vector<PDEMap> *common = nullptr) {
  double sumW = 0., sumWE = 0.;
  for (int ix = 1; ix <= acc->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= acc->GetNbinsY(); iy++) {
      double w = acc->GetBinContent(ix, iy);
      if (w <= 0) continue;
      double ycen = acc->GetXaxis()->GetBinCenter(ix);
      double xcen = acc->GetYaxis()->GetBinCenter(iy);
      if (common && !covered_by_all(*common, ycen, xcen)) continue;
      sumW  += w;
      sumWE += w * map->GetBinContent(map->FindBin(ycen, xcen));
    }
  }
  return sumW > 0 ? sumWE / sumW : NAN;
}

// Fraction of the acceptance weight lying at positions where the map is measured
// (>= kMinAllCounts "all" events). Positions not covered get efficiency 0 in the MC.
double coverage_fraction(const PDEMap &m, TH2 *acc) {
  if (!m.nAll) return NAN;
  double sumW = 0., sumWin = 0.;
  for (int ix = 1; ix <= acc->GetNbinsX(); ix++) {
    for (int iy = 1; iy <= acc->GetNbinsY(); iy++) {
      double w = acc->GetBinContent(ix, iy);
      if (w <= 0) continue;
      double ycen = acc->GetXaxis()->GetBinCenter(ix);
      double xcen = acc->GetYaxis()->GetBinCenter(iy);
      sumW += w;
      if (m.nAll->GetBinContent(m.nAll->FindBin(ycen, xcen)) >= kMinAllCounts) sumWin += w;
    }
  }
  return sumW > 0 ? sumWin / sumW : NAN;
}

int compare_pDE_maps(const char *fileList, const char *labelList, const char *outBase,
                     const char *nestedPairs = "", int weightIndex = 0) {

  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gStyle->SetPalette(kRainbow);
  gStyle->SetNumberContours(50);

  // ---- read inputs -------------------------------------------------------
  std::vector<TString> files = split_csv(fileList);
  std::vector<TString> labels = split_csv(labelList);
  if (files.size() < 2 || files.size() != labels.size()) {
    std::cerr << "Need >= 2 files and one label per file." << std::endl;
    return 1;
  }
  std::vector<PDEMap> maps;
  for (size_t i = 0; i < files.size(); i++) {
    TFile *f = TFile::Open(files[i]);
    if (!f || f->IsZombie()) { std::cerr << "Cannot open " << files[i] << std::endl; return 1; }
    TString tag = Form("m%zu", i);
    PDEMap m;
    m.label = labels[i];
    m.effi  = get_hist<TH2>(f, "h2_effi_map", tag);
    m.nAll  = get_hist<TH2>(f, "h2_xyexp_earm", tag);
    m.nPass = get_hist<TH2>(f, "h2_xyexp_eNharm", tag);
    m.effX  = get_hist<TH1>(f, "h_xexp_pDE", tag);
    m.effY  = get_hist<TH1>(f, "h_yexp_pDE", tag);
    m.accP  = get_hist<TH2>(f, "h2_xyexp_all", tag);
    m.accN  = get_hist<TH2>(f, "h2_xyexp_all_nodef", tag);
    f->Close();
    if (!m.effi) { std::cerr << "No h2_effi_map in " << files[i] << std::endl; return 1; }
    if (!m.nAll || !m.nPass)
      std::cout << "Warning: no pass/all counts in " << files[i] << ", using stored map errors\n";
    maps.push_back(m);
  }

  // Which pairs are nested (sample i is a subset of sample j)?
  std::set<std::pair<int,int>> nested;
  for (const TString &p : split_csv(nestedPairs)) {
    int i = -1, j = -1;
    if (sscanf(p.Data(), "%d:%d", &i, &j) == 2) { nested.insert({i, j}); nested.insert({j, i}); }
  }

  TString pdf = TString(outBase) + ".pdf";
  std::ofstream csv(TString(outBase) + ".csv");
  csv << "quantity,label,value,error\n";
  TLatex tex; tex.SetNDC(); tex.SetTextSize(0.035);

  // ---- page 1: the maps --------------------------------------------------
  int nMaps = maps.size();
  TCanvas *cMaps = new TCanvas("cMaps", "maps", 600 * nMaps, 600);
  cMaps->Divide(nMaps, 1);
  for (int i = 0; i < nMaps; i++) {
    cMaps->cd(i + 1);
    gPad->SetRightMargin(0.15);
    maps[i].effi->SetTitle(maps[i].label + ";y_{HCAL}^{exp} (m);x_{HCAL}^{exp} - #deltax_{SBS} (m)");
    maps[i].effi->GetZaxis()->SetRangeUser(0.5, 1.0);
    maps[i].effi->Draw("colz");
  }
  cMaps->Print(pdf + "[");
  cMaps->Print(pdf);

  // ---- per pair: pull map, pull distribution, chi2, ratio map ------------
  for (int a = 0; a < nMaps; a++) {
    for (int b = a + 1; b < nMaps; b++) {
      const PDEMap &A = maps[a], &B = maps[b];
      bool isNested = nested.count({a, b}) > 0;
      TString pairLabel = A.label + "_vs_" + B.label;

      // Unused bins get a value below the colour range so that "colz" leaves them blank.
      TH2 *hPull  = (TH2*)A.effi->Clone("hPull_" + pairLabel);   hPull->Reset();
      for (int k = 0; k < hPull->GetNcells(); k++) hPull->SetBinContent(k, kEmptyBin);
      TH2 *hRatio = (TH2*)A.effi->Clone("hRatio_" + pairLabel);  hRatio->Reset();
      TH1F *hPullDist   = new TH1F("hPullDist_" + pairLabel, "", 50, -kPullRange, kPullRange);
      TH1F *hPullDistFR = new TH1F("hPullDistFR_" + pairLabel, "", 50, -kPullRange, kPullRange);
      double chi2 = 0., chi2FR = 0.; int ndf = 0, ndfFR = 0;

      for (int ix = 1; ix <= A.effi->GetNbinsX(); ix++) {
        for (int iy = 1; iy <= A.effi->GetNbinsY(); iy++) {
          int bin = A.effi->GetBin(ix, iy);
          if (bin_all_counts(A, bin) < kMinAllCounts || bin_all_counts(B, bin) < kMinAllCounts) continue;
          double vA = bin_variance(A, bin), vB = bin_variance(B, bin);
          if (vA < 0 || vB < 0) continue;
          // independent samples add variances; for nested samples the shared events cancel
          double var = isNested ? std::fabs(vA - vB) : vA + vB;
          if (var <= 0) continue;
          double eA = A.effi->GetBinContent(bin), eB = B.effi->GetBinContent(bin);
          double pull = (eA - eB) / std::sqrt(var);
          hPull->SetBinContent(bin, pull);
          if (eB > 0) hRatio->SetBinContent(bin, eA / eB);
          hPullDist->Fill(pull);  chi2 += pull * pull;  ndf++;
          double ycen = A.effi->GetXaxis()->GetBinCenter(ix);
          double xcen = A.effi->GetYaxis()->GetBinCenter(iy);
          if (in_fit_region(ycen, xcen)) { hPullDistFR->Fill(pull); chi2FR += pull * pull; ndfFR++; }
        }
      }

      TCanvas *cP = new TCanvas("cP_" + pairLabel, "", 1200, 1000);
      cP->Divide(2, 2);
      cP->cd(1); gPad->SetRightMargin(0.15);
      hPull->SetTitle("Pull (" + A.label + " - " + B.label + ")/#sigma" + (isNested ? " [nested]" : "")
                      + ";y_{HCAL}^{exp} (m);x_{HCAL}^{exp} - #deltax_{SBS} (m)");
      hPull->GetZaxis()->SetRangeUser(-kPullRange, kPullRange);
      hPull->Draw("colz");
      cP->cd(2); gPad->SetRightMargin(0.15);
      hRatio->SetTitle("Ratio " + A.label + " / " + B.label + ";y_{HCAL}^{exp} (m);x_{HCAL}^{exp} - #deltax_{SBS} (m)");
      hRatio->GetZaxis()->SetRangeUser(kRatioMin, kRatioMax);
      hRatio->Draw("colz");
      cP->cd(3);
      hPullDist->SetTitle(Form("All bins with N_{all} #geq %.0f;pull;bins", kMinAllCounts));
      hPullDist->SetLineWidth(2);
      hPullDist->Draw("hist");
      tex.DrawLatex(0.15, 0.85, Form("mean = %.2f, RMS = %.2f", hPullDist->GetMean(), hPullDist->GetRMS()));
      tex.DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2, ndf, ndf ? chi2 / ndf : NAN));
      cP->cd(4);
      hPullDistFR->SetTitle(Form("Fit region x#in[%.1f,%.1f], y#in[%.1f,%.1f];pull;bins",
                                 kFitXmin, kFitXmax, kFitYmin, kFitYmax));
      hPullDistFR->SetLineWidth(2);
      hPullDistFR->Draw("hist");
      tex.DrawLatex(0.15, 0.85, Form("mean = %.2f, RMS = %.2f", hPullDistFR->GetMean(), hPullDistFR->GetRMS()));
      tex.DrawLatex(0.15, 0.80, Form("#chi^{2}/ndf = %.1f/%d = %.2f", chi2FR, ndfFR, ndfFR ? chi2FR / ndfFR : NAN));
      cP->Print(pdf);

      csv << "pull_mean_all," << pairLabel << "," << hPullDist->GetMean() << "," << hPullDist->GetMeanError() << "\n";
      csv << "pull_rms_all," << pairLabel << "," << hPullDist->GetRMS() << ",\n";
      csv << "chi2ndf_all," << pairLabel << "," << (ndf ? chi2 / ndf : NAN) << "," << ndf << "\n";
      csv << "pull_mean_fitreg," << pairLabel << "," << hPullDistFR->GetMean() << "," << hPullDistFR->GetMeanError() << "\n";
      csv << "pull_rms_fitreg," << pairLabel << "," << hPullDistFR->GetRMS() << ",\n";
      csv << "chi2ndf_fitreg," << pairLabel << "," << (ndfFR ? chi2FR / ndfFR : NAN) << "," << ndfFR << "\n";
      std::cout << Form("%-28s %s pull mean/RMS all = %6.2f/%5.2f chi2/ndf = %6.1f/%d ; fit region = %6.2f/%5.2f chi2/ndf = %6.1f/%d",
                        pairLabel.Data(), isNested ? "[nested]" : "        ",
                        hPullDist->GetMean(), hPullDist->GetRMS(), chi2, ndf,
                        hPullDistFR->GetMean(), hPullDistFR->GetRMS(), chi2FR, ndfFR) << std::endl;
    }
  }

  // ---- 1D efficiency overlays + pol0 averages -----------------------------
  const int colors[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7};
  const int nColors = sizeof(colors) / sizeof(int);
  TCanvas *c1D = new TCanvas("c1D", "", 1200, 1000);
  c1D->Divide(2, 2);
  TLegend *leg = new TLegend(0.15, 0.15, 0.6, 0.4);
  leg->SetFillStyle(0);
  for (int i = 0; i < nMaps; i++) {
    PDEMap &m = maps[i];
    int col = colors[i % nColors];
    auto px = fit_pol0(m.effX, kFitXmin, kFitXmax);
    auto py = fit_pol0(m.effY, kFitYmin, kFitYmax);
    auto eAll = overall_eff(m, false);
    auto eFR  = overall_eff(m, true);
    csv << "pol0_x," << m.label << "," << px.first << "," << px.second << "\n";
    csv << "pol0_y," << m.label << "," << py.first << "," << py.second << "\n";
    csv << "pass_over_all," << m.label << "," << eAll.first << "," << eAll.second << "\n";
    csv << "pass_over_all_fitreg," << m.label << "," << eFR.first << "," << eFR.second << "\n";
    std::cout << Form("%-20s pol0 x = %.4f +- %.4f, pol0 y = %.4f +- %.4f, pass/all = %.4f +- %.4f (fit region %.4f +- %.4f), N_all = %.0f",
                      m.label.Data(), px.first, px.second, py.first, py.second,
                      eAll.first, eAll.second, eFR.first, eFR.second,
                      m.nAll ? m.nAll->Integral() : NAN) << std::endl;

    // overlays
    for (TH1 *h : {m.effX, m.effY}) {
      if (!h) continue;
      h->SetLineColor(col); h->SetMarkerColor(col); h->SetMarkerStyle(20); h->SetMarkerSize(0.6);
    }
    const char *opt = (i == 0) ? "E" : "E same";
    if (m.effX) {
      c1D->cd(1); m.effX->SetTitle(";x_{HCAL}^{exp} - #deltax_{SBS} (m);Efficiency");
      m.effX->GetYaxis()->SetRangeUser(0.5, 1.05); m.effX->Draw(opt);
    }
    if (m.effY) {
      c1D->cd(2); m.effY->SetTitle(";y_{HCAL}^{exp} (m);Efficiency");
      m.effY->GetYaxis()->SetRangeUser(0.5, 1.05); m.effY->Draw(opt);
    }
    leg->AddEntry(m.effX ? m.effX : m.effY, Form("%s: x %.4f, y %.4f", m.label.Data(), px.first, py.first), "lp");

    // ratio to the first map (errors treated as independent: visual guide only)
    if (i > 0) {
      for (int k = 0; k < 2; k++) {
        TH1 *h = (k == 0) ? m.effX : m.effY;
        TH1 *h0 = (k == 0) ? maps[0].effX : maps[0].effY;
        if (!h || !h0) continue;
        TH1 *r = (TH1*)h->Clone(Form("%s_ratio", h->GetName()));
        r->Divide(h0);
        r->SetStats(0);
        r->SetTitle(Form(";%s;ratio to %s", k == 0 ? "x_{HCAL}^{exp} - #deltax_{SBS} (m)" : "y_{HCAL}^{exp} (m)",
                         maps[0].label.Data()));
        r->GetYaxis()->SetRangeUser(kRatioMin, kRatioMax);
        c1D->cd(3 + k);
        r->Draw(i == 1 ? "E" : "E same");
        TLine *one = new TLine(r->GetXaxis()->GetXmin(), 1., r->GetXaxis()->GetXmax(), 1.);
        one->SetLineStyle(2); one->Draw();
      }
    }
  }
  c1D->cd(1); leg->Draw();
  c1D->Print(pdf);

  // ---- impact on MC weight_effic ------------------------------------------
  // In the MC, weight_effic = weight * map(y, x_hit) / const, where x_hit is the expected
  // proton position (x - p_def) for protons and the undeflected position x for neutrons.
  // As a proxy for the LD2 MC acceptance we use the LH2 "all" distributions of one input:
  //   proton positions  : h2_xyexp_all       (y, x - p_def)
  //   neutron positions : h2_xyexp_all_nodef (y, x)
  // <map>_p and <map>_n are the mean MC weights for p and n; their ratio enters R = n/p.
  const PDEMap &W = maps[weightIndex];
  if (W.accP && W.accN) {
    TCanvas *cImp = new TCanvas("cImp", "", 1200, 600);
    cImp->Divide(2, 1);
    cImp->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.15);
    W.accP->SetTitle("Proton positions (" + W.label + ");y_{HCAL}^{exp} (m);x_{HCAL}^{exp} - #deltax_{SBS} (m)");
    W.accP->Draw("colz");
    cImp->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.15);
    W.accN->SetTitle("Neutron positions (" + W.label + ");y_{HCAL}^{exp} (m);x_{HCAL}^{exp} (m)");
    W.accN->Draw("colz");
    cImp->Print(pdf);

    TCanvas *cTxt = new TCanvas("cTxt", "", 1200, 800);
    TLatex t2; t2.SetNDC(); t2.SetTextSize(0.024); t2.SetTextFont(82);
    double yPos = 0.94;
    t2.DrawLatex(0.02, yPos, Form("Acceptance-weighted map averages (weights: %s LH2 'all' events)", W.label.Data()));
    yPos -= 0.04;
    t2.DrawLatex(0.02, yPos, "Averages use only positions measured by ALL maps (N_{all} >= threshold);");
    yPos -= 0.04;
    t2.DrawLatex(0.02, yPos, "cov_p/cov_n = fraction of p/n weight where that map alone is measured.");
    yPos -= 0.06;
    t2.DrawLatex(0.02, yPos, "map            cov_p  cov_n   <m>_p    <m>_n   <m>_n/<m>_p  rel.chg vs first p, n, n/p (%)");
    double p0 = weighted_map_average(maps[0].effi, W.accP, &maps);
    double n0 = weighted_map_average(maps[0].effi, W.accN, &maps);
    for (int i = 0; i < nMaps; i++) {
      double covP = coverage_fraction(maps[i], W.accP);
      double covN = coverage_fraction(maps[i], W.accN);
      double avgP = weighted_map_average(maps[i].effi, W.accP, &maps);
      double avgN = weighted_map_average(maps[i].effi, W.accN, &maps);
      double relP = 100. * (avgP / p0 - 1.), relN = 100. * (avgN / n0 - 1.);
      double relNP = 100. * ((avgN / avgP) / (n0 / p0) - 1.);
      csv << "coverage_proton," << maps[i].label << "," << covP << ",\n";
      csv << "coverage_neutron," << maps[i].label << "," << covN << ",\n";
      csv << "wavg_proton_common," << maps[i].label << "," << avgP << ",\n";
      csv << "wavg_neutron_common," << maps[i].label << "," << avgN << ",\n";
      csv << "wavg_n_over_p_common," << maps[i].label << "," << avgN / avgP << ",\n";
      csv << "relchg_pct_vs_" << maps[0].label << "_proton," << maps[i].label << "," << relP << ",\n";
      csv << "relchg_pct_vs_" << maps[0].label << "_neutron," << maps[i].label << "," << relN << ",\n";
      csv << "relchg_pct_vs_" << maps[0].label << "_n_over_p," << maps[i].label << "," << relNP << ",\n";
      TString line = Form("%-14s %.3f  %.3f   %.4f   %.4f   %.4f       %+.3f, %+.3f, %+.3f",
                          maps[i].label.Data(), covP, covN, avgP, avgN, avgN / avgP, relP, relN, relNP);
      yPos -= 0.04;
      t2.DrawLatex(0.02, yPos, line);
      std::cout << "impact: " << line << std::endl;
    }
    cTxt->Print(pdf);
  } else {
    std::cout << "No h2_xyexp_all/_nodef in weight file: impact estimate skipped." << std::endl;
  }

  cMaps->Print(pdf + "]");
  csv.close();
  std::cout << "\n Output: " << pdf << "\n         " << outBase << ".csv" << std::endl;
  return 0;
}
