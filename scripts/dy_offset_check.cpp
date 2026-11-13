/*
  dy_offset_check.cpp : check the MC horizontal HCAL offset (hcal_hoffset) with the data/MC dy peak
                        positions after the analysis cuts, separately for LH2 and LD2.

  Physics / method
  ----------------
  dy = yHCAL - yHCAL_exp is the horizontal distance between the measured HCAL cluster and the
  position predicted from the electron (q-vector) kinematics. The analysis cuts on |dy| (e.g.
  |dy| < 0.3 m in fit_dx), so the data and MC dy peaks must sit at the same place, otherwise the
  cut removes a different fraction of the events in data and MC.

  How hcal_hoffset enters the MC dy (simulation/elas_ana_simu.cpp and qelas_ana_simu.cpp):
    HCAL_origin = (d_HCAL + zoff) z^ + voff x^ + hoff y^              (elas:348, qelas:368)
    yHCAL_exp   = (HCAL_intersect - HCAL_origin) . y^                (src/KinematicVar.cpp:122)
    dy          = yHCAL - yHCAL_exp ; stored dy = dy + dy_offset     (elas:573/578, qelas:647/652)
  The HCAL plane is not moved by hoff (y^ is in the plane), so HCAL_intersect does not change and
  yHCAL_exp(hoff) = yHCAL_exp(0) - hoff, i.e. MC dy(hoff) = dy(0) + hoff exactly; the measured
  yHCAL does not depend on it. Hence, with
    Delta = <dy>_data - <dy>_MC   (MC made with hoff_file),
  the MC needs  hoff_new = hoff_file + Delta.  (The data macros use their own hcal_hoffset in the
  same way; it is a different number and is only printed for information.)
  The current hoff (and dy_offset) are read from the summary canvas stored in each MC file, and
  cross-checked with the MC analysis config.

  Event selection
    LH2 : the strict elastic cuts of w2_offset_check (eHCAL > min, bbfiduCut, |vz| < cut, data
          coincidence-time window, run cut; same files, from configs/conf_w2_offset_check.json),
          plus a W2 window. The thpq_p cut of w2_offset_check is NOT used in the nominal selection:
          thpq_p is the angle between the predicted and the observed nucleon direction, so it
          contains dy and would truncate the dy distribution around 0. It is replaced by a dx
          window on the proton spot; the thpq_p selection is shown as a reference variation.
    LD2 : the production cuts of the "qelas" block of configs/<setting>/conf_fit_dx.json, built
          exactly as fit_dx.cpp does (data = cuts_for_signal_simu && extra_data_signal_cut,
          "&&nblkHCAL>1" when apply_nblkHCAL>1_cut = 1, the custom fiducial cut when
          use_custom_fiduCut = 1; same data/MC files and MC weight). The |dy| cut term itself is
          removed (it would bias the peak position measurement); the peak is measured within fit
          windows instead. The proton and neutron spots are selected with dx windows of
          +-kSpotNsigma sigma around each spot (data and MC each around their own spot centre;
          MC p/n separated with mc_fnucl), "all" is the union of both windows.

  Peak position difference Delta, as in w2_offset_check:
    (a) core Gaussian : iterative Gaussian fit to the dy peak core of data and MC, Delta = mu_d - mu_MC;
    (b) template scan : data = N * MC(dy - delta) + pol2 (background, e.g. inelastic events in LD2),
                        solved by linear least squares for each delta; best delta from a parabola
                        through chi2(delta); "stretched" version also scales the MC width.
    best = mean of (a) and the stretched (b);
    systematic = largest shift over cut variations (+) spread over fit windows (+) spread between
                 methods, added in quadrature. Variations: physics cuts (refilled), spot window
                 width and a central W2 half (from the same histograms). Reference variations
                 (thpq_p for LH2, no ARCut for LD2) are shown but not included.
  dy is correlated with W2 (both come from the electron kinematics), so Delta is also measured in
  W2 slices of equal data population (core-Gaussian method); a straight-line fit gives dDelta/dW2.
  Verdict on |Delta|: "ok" below kOkThreshold, otherwise "FLAG" if also > kFlagNsigma total sigma,
  "check" if not.

  Caveat: the LD2 data files (qelas_ana_data) used here were produced before the fix of the ARCut
  bug (ARCut evaluated at the wrong HCAL cluster, 0.3-7% of the events). The "no_ARCut" reference
  variation shows how sensitive the dy peak is to ARCut.

  Inputs : configs/conf_dy_offset_check.json (+ the configs it points to).
  Outputs (in "out_dir", default simulation/siout/dy_offset_check):
    <out_file_base>.pdf  - notes page, a few pages per setting and target, summary pages
    <out_file_base>.root - the dy histograms and W2-slice graphs
    <out_file_base>.csv  - one row per setting, target and spot
  Usage (from the scripts/ directory):
    root -l -b -q rootlogon.C 'dy_offset_check.cpp("configs/conf_dy_offset_check.json")'
    root -l -b -q rootlogon.C 'dy_offset_check.cpp("configs/conf_dy_offset_check.json","sbs9_70p_LD2")'
  The output file names of a single-block run get the block name appended.
  -------
  P. Datta, 2026-11-12 (modelled on w2_offset_check.cpp)
*/

#include <cmath>
#include <cstdio>
#include <cstring>
#include <cctype>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TPaveText.h"
#include "TText.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "ROOT/RDataFrame.hxx"

#include "../include/gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

// ---------------------------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------------------------
// Histogram axes: dy (m), dx (m), W2 (GeV2). dy bins of 1 cm are used for all fits (the MC
// template is interpolated on them); dx only selects the spots; W2 only for slices.
const int kNbinsDy = 200; const double kDyMin = -1.0, kDyMax = 1.0;
const int kNbinsDx = 140; const double kDxMin = -2.5, kDxMax = 1.0;
const int kNbinsW2 = 80;  const double kW2Min = -0.2, kW2Max = 1.8;

// LH2 data coincidence-time peak, found before the main event loop (same recipe as w2_offset_check).
const int    kNbinsCoinT = 80; const double kCoinTMin = -20, kCoinTMax = 20;   // ns
const double kCoinTPeakSearch = 3.0;                 // ns, half-width around the maximum for the mean
const double kElasW2Lo = 0.6, kElasW2Hi = 1.3;       // GeV2, elastic events used to find the peak
const Long64_t kNentriesForCoinT = 3000000;          // first pass on these entries only ...
const int    kMinCoinTEvents = 2000;                 // ... unless fewer events pass (run-cut files)

// Core Gaussian of a dy peak: maximum searched within |dy| < kDyCoreSearch, fit +- kCoreNsigma sigma.
const double kDyCoreSearch = 0.4;      // m
const double kDyCoreSigmaStart = 0.08; // m
const double kCoreNsigma = 1.2;
const int    kCoreIterations = 6;

// dx spots. Spot centre and width from a core Gaussian of the dx distribution (|dy| < kDySpotFind),
// data spot searched within +- kDxSpotSearch of the expected position.
const double kSpotNsigma = 1.5;                         // nominal dx window half-width (sigma)
const std::vector<double> kSpotNsigmaVar = {1.0, 2.0};  // systematic variations of it
const double kDySpotFind = 0.3;                         // m
const double kDxSpotSearch = 0.3;                       // m
const double kDxCoreSigmaStart = 0.15;                  // m

// Template scan (as in w2_offset_check), delta in m.
const double kScanStep = 0.0005, kScanHalfWidth = 0.04;
const double kScanParabolaDchi2 = 4;
const std::vector<std::pair<double,double>> kFitWindows = {{3,3}, {2,2}, {4,4}, {3,2}, {2,3}};   // data sigma units, first nominal
const double kFitWindowEdge = 0.9;      // m, fit windows are kept within |dy| < this (histogram edge 1.0)
const int    kNominalPolOrder = 2;
const double kStretchMin = 0.70, kStretchMax = 1.40, kStretchStep = 0.02;

// Physics cut variations.
const double kVzTight = 0.05;           // m, LD2 tighter vertex cut
const std::vector<double> kCoinTWindowVariations = {1.5, 6.0};   // ns, LH2 data only
// Central W2 half (data quantiles) used as a variation.
const double kW2CoreQlo = 0.25, kW2CoreQhi = 0.75;

// W2 slices: equal data population; slices with fewer data events are skipped.
const int    kNW2Slices = 5;
const double kMinSliceEvents = 200;

// Verdict on |Delta| (m): "ok" below kOkThreshold (3 mm, a small fraction of the dy cut width);
// above it, "FLAG" if also more than kFlagNsigma total standard deviations, otherwise "check".
const double kOkThreshold = 0.003;
const double kFlagNsigma = 2.0;

// Page layout
const int kPageWidth = 1500, kPageHeight = 650;   // pixels
const double kMinDisplayDchi2 = 25;

// Unique-name counter: ROOT replaces objects of the same name in the current directory.
int gUid = 0;
std::string Uid(const std::string &base) { return base + "_" + std::to_string(gUid++); }

// ---------------------------------------------------------------------------------------------
// Small containers
// ---------------------------------------------------------------------------------------------
struct Peak { double mu = 0, emu = -1, sigma = 0; };

struct ScanResult {
  double delta = 0, err = -1;        // best shift and its error (err < 0: minimum not found)
  double chi2 = 0; int ndf = 0;
  double stretch = 1;                // MC width scale (1 for the plain template)
  double lo = 0, hi = 0; int npol = 0;
  std::vector<double> par;           // N, polynomial coefficients
  std::vector<double> scanDelta, scanChi2;
};

// One selection: data cut, MC cut, and whether it enters the systematic error.
struct CutSet { std::string name, dataCut, mcCut; bool inSys; };

// dx bin ranges of a spot selection (one or two ranges, see BinRanges).
typedef std::vector<std::pair<int,int>> Ranges;

// dx windows of one spot for data and MC: nominal, narrow and wide.
struct SpotWindows { std::string name; Ranges dNom, mNom, dNarrow, mNarrow, dWide, mWide; };

struct SpotResult {
  std::string spot;
  Peak data{}, mc{};
  double nData = 0;
  double offCore = 0, errCore = 0;
  ScanResult tpl, str;
  double best = 0, stat = 0, sysCut = 0, sysFit = 0, sysMethod = 0, sys = 0;
  std::string worstCut, verdict;
  std::vector<std::string> varName; std::vector<double> varOffset, varError; std::vector<bool> varInSys;
  std::vector<double> w2Mean, w2Delta, w2Err, w2MuData, w2MuMC, w2ErrData, w2ErrMC;
  double slope = 0, slopeErr = 0, constChi2 = 0; int constNdf = 0;
  TH1D *hData = nullptr, *hMC = nullptr;
};

// Everything needed to process one block of the config (one setting and target).
struct CaseInput {
  std::string name, setting, target;
  bool isLD2 = false;
  std::string dataFile; std::vector<std::string> simuFiles; std::vector<double> simuWeights;
  std::string mcWeight;
  std::vector<CutSet> cuts;          // [0] = nominal
  std::vector<std::string> removedTerms;
  double dxOffP = 0, dxOffN = 0;     // fit_dx MC->data dx shifts (LD2), to locate the data spots
  // custom fiducial cut of fit_dx (use_custom_fiduCut = 1)
  bool customFidu = false; std::vector<double> fiduAR, fiduSM; double sbsKick = 0;
  // hcal_hoffset bookkeeping (m)
  double hoffFile = 0, dyOffFile = 0, hoffConfig = 0, hoffData = 0;
  bool hoffFromFile = false;
  bool ok = true;
};

struct CaseResult {
  CaseInput in;
  std::vector<SpotResult> spots;
  double dataCenter[2] = {0, 0}, mcCenter[2] = {0, 0}, mcSigma[2] = {0, 0};   // [p, n]
  int nSpots = 1;
};

// ---------------------------------------------------------------------------------------------
// Config / cut-string helpers
// ---------------------------------------------------------------------------------------------
bool HasSubKey(JSONManager *j, const std::string &key, const std::string &sub) {
  std::vector<std::string> subs; j->GetSubKeys(key, subs);
  return std::find(subs.begin(), subs.end(), sub) != subs.end();
}

// Split a cut string at the top-level "&&" (outside parentheses).
std::vector<std::string> SplitAnd(const std::string &cut) {
  std::vector<std::string> terms; std::string cur; int depth = 0;
  for (size_t i = 0; i < cut.size(); i++) {
    char c = cut[i];
    if (c == '(') depth++;
    if (c == ')') depth--;
    if (depth == 0 && c == '&' && i + 1 < cut.size() && cut[i + 1] == '&') { terms.push_back(cur); cur.clear(); i++; continue; }
    cur += c;
  }
  if (!cur.empty()) terms.push_back(cur);
  return terms;
}

// True if "term" uses the variable "var" as a whole identifier (so "dy" does not match "dy_nS").
bool UsesVariable(const std::string &term, const std::string &var) {
  auto isId = [](char c) { return std::isalnum((unsigned char)c) || c == '_'; };
  for (size_t p = term.find(var); p != std::string::npos; p = term.find(var, p + 1)) {
    bool before = (p == 0) || !isId(term[p - 1]);
    bool after = (p + var.size() >= term.size()) || !isId(term[p + var.size()]);
    if (before && after) return true;
  }
  return false;
}

// Remove all top-level terms that use "var"; the removed terms are appended to "removed".
std::string RemoveTerms(const std::string &cut, const std::string &var, std::vector<std::string> &removed) {
  std::string out;
  for (auto &t : SplitAnd(cut)) {
    if (UsesVariable(t, var)) { removed.push_back(t); continue; }
    out += (out.empty() ? "" : "&&") + t;
  }
  return out.empty() ? "1" : out;
}

// hcal_hoffset and dy_offset printed in the summary canvas of an analysis output file.
bool ReadSummaryOffsets(const std::string &file, double &hoff, double &dyOff) {
  TFile f(file.c_str());
  if (f.IsZombie()) return false;
  TCanvas *c = (TCanvas*)f.Get("cSummary");
  if (!c) return false;
  bool found = false;
  TIter next(c->GetListOfPrimitives());
  while (TObject *o = next()) {
    if (!o->InheritsFrom("TPaveText")) continue;
    TIter nl(((TPaveText*)o)->GetListOfLines());
    while (TObject *l = nl()) {
      const char *s = ((TText*)l)->GetTitle();
      double v, h, z;
      if (const char *p = strstr(s, "HCAL offsets:"))
        if (sscanf(p, "HCAL offsets: v = %lf, h = %lf, z = %lf", &v, &h, &z) == 3) { hoff = h; found = true; }
      if (const char *p = strstr(s, "dy offset ="))
        sscanf(p, "dy offset = %lf", &dyOff);
    }
  }
  return found;
}

// Mean of the e-HCAL coincidence-time peak of LH2 elastic events (data).
double ReadCoinTPeak(const std::string &file, const std::string &sel) {
  TFile f(file.c_str());
  TTree *t = (TTree*)f.Get("Tout");
  if (!t) { printf("ERROR: no Tout in %s\n", file.c_str()); return 0; }
  TH1D h("hcoinT_tmp", "", kNbinsCoinT, kCoinTMin, kCoinTMax);
  Long64_t n = t->Draw("coinT_ADC_c>>hcoinT_tmp", sel.c_str(), "goff", kNentriesForCoinT);
  if (n < kMinCoinTEvents) t->Draw("coinT_ADC_c>>hcoinT_tmp", sel.c_str(), "goff");   // whole file
  double cmax = h.GetBinCenter(h.GetMaximumBin());
  h.GetXaxis()->SetRangeUser(cmax - kCoinTPeakSearch, cmax + kCoinTPeakSearch);
  return h.GetMean();
}

// ---------------------------------------------------------------------------------------------
// Histogram helpers
// ---------------------------------------------------------------------------------------------
// Iterative Gaussian fit to the peak core (+-kCoreNsigma sigma around the current mean); the
// maximum is searched in [lo, hi].
Peak CorePeak(TH1D *h, double lo, double hi, double sigmaStart) {
  Peak p;
  if (h->Integral() <= 0) return p;
  h->GetXaxis()->SetRangeUser(lo, hi);
  double mu = h->GetBinCenter(h->GetMaximumBin()), sigma = sigmaStart;
  h->GetXaxis()->SetRange(0, 0);
  TF1 g("gcore", "gaus");
  for (int i = 0; i < kCoreIterations; i++) {
    g.SetRange(mu - kCoreNsigma * sigma, mu + kCoreNsigma * sigma);
    g.SetParameters(h->GetBinContent(h->FindBin(mu)), mu, sigma);
    h->Fit(&g, "QR0");
    const double s = std::abs(g.GetParameter(2));
    if (!(s > 0 && s < (hi - lo))) break;   // keep the last sensible values if a fit runs away
    mu = g.GetParameter(1); sigma = s;
  }
  p.mu = mu; p.emu = g.GetParError(1); p.sigma = sigma;
  return p;
}
Peak DyPeak(TH1D *h) { return CorePeak(h, -kDyCoreSearch, kDyCoreSearch, kDyCoreSigmaStart); }

// Bin ranges of one or more dx windows [lo, hi]: bins whose centre lies inside; overlapping
// windows are merged so that no bin is counted twice.
Ranges BinRanges(const TAxis *ax, std::vector<std::pair<double,double>> win) {
  std::sort(win.begin(), win.end());
  Ranges r;
  for (auto &w : win) {
    int b1 = std::max(1, ax->FindFixBin(w.first)), b2 = std::min(ax->GetNbins(), ax->FindFixBin(w.second));
    if (ax->GetBinCenter(b1) < w.first) b1++;
    if (ax->GetBinCenter(b2) > w.second) b2--;
    if (b2 < b1) continue;
    if (!r.empty() && b1 <= r.back().second + 1) r.back().second = std::max(r.back().second, b2);
    else r.push_back({b1, b2});
  }
  return r;
}

// dy projection of the (dy, dx, W2) histogram over the dx ranges and W2 bins [z1, z2].
TH1D *ProjectDy(const TH3D *h3, const Ranges &rx, int z1, int z2, const std::string &name) {
  TH1D *out = nullptr;
  for (auto &r : rx) {
    TH1D *p = h3->ProjectionX(Uid(name).c_str(), r.first, r.second, z1, z2, "e");
    p->SetDirectory(0);
    if (!out) { out = p; out->SetName(name.c_str()); } else { out->Add(p); delete p; }
  }
  if (!out) { out = h3->ProjectionX(Uid(name).c_str(), 1, 0, z1, z2, "e"); out->Reset(); out->SetDirectory(0); }
  return out;
}
// dy projection of a (dy, dx) histogram over the dx ranges.
TH1D *ProjectDy(const TH2D *h2, const Ranges &rx, const std::string &name) {
  TH1D *out = nullptr;
  for (auto &r : rx) {
    TH1D *p = h2->ProjectionX(Uid(name).c_str(), r.first, r.second, "e");
    p->SetDirectory(0);
    if (!out) { out = p; out->SetName(name.c_str()); } else { out->Add(p); delete p; }
  }
  if (!out) { out = h2->ProjectionX(Uid(name).c_str(), 1, 0, "e"); out->Reset(); out->SetDirectory(0); }
  return out;
}
// W2 projection over the dx ranges (all dy).
TH1D *ProjectW2(const TH3D *h3, const Ranges &rx, const std::string &name) {
  TH1D *out = nullptr;
  for (auto &r : rx) {
    TH1D *p = h3->ProjectionZ(Uid(name).c_str(), 1, h3->GetNbinsX(), r.first, r.second, "e");
    p->SetDirectory(0);
    if (!out) out = p; else { out->Add(p); delete p; }
  }
  return out;
}
// dx projection with |dy| < kDySpotFind (all W2), used to locate the spots.
TH1D *ProjectDx(const TH3D *h3, const std::string &name) {
  const int b1 = h3->GetXaxis()->FindFixBin(-kDySpotFind + 1e-6), b2 = h3->GetXaxis()->FindFixBin(kDySpotFind - 1e-6);
  TH1D *p = h3->ProjectionY(Uid(name).c_str(), b1, b2, 1, h3->GetNbinsZ(), "e");
  p->SetDirectory(0);
  return p;
}
// W2 bin range [first, last] holding the data quantiles [qlo, qhi] of hW2.
std::pair<int,int> QuantileBins(TH1D *hW2, double qlo, double qhi) {
  double p[2] = {qlo, qhi}, q[2];
  hW2->GetQuantiles(2, q, p);
  return {hW2->FindFixBin(q[0]), hW2->FindFixBin(q[1])};
}

// ---------------------------------------------------------------------------------------------
// Template scan (same algorithm as w2_offset_check.cpp)
// ---------------------------------------------------------------------------------------------
// MC template value at dy = x for shift delta and width scale k (stretched about the MC peak muMC).
double TemplateAt(const TH1D *mc, double x, double delta, double k, double muMC) {
  return const_cast<TH1D*>(mc)->Interpolate(muMC + (x - delta - muMC) / k) / k;
}

// chi2 of data = N*template(delta,k) + sum_j c_j x^j in [lo,hi], solved by weighted linear least squares.
double LinearChi2(const TH1D *data, const TH1D *mc, double delta, double k, double muMC,
                  int npol, double lo, double hi, std::vector<double> &par, int &ndf) {
  const int np = npol + 2;
  TMatrixD A(np, np); TVectorD B(np);
  std::vector<std::vector<double>> rows; std::vector<double> ys, ws;
  for (int i = data->GetXaxis()->FindFixBin(lo + 1e-9); i <= data->GetXaxis()->FindFixBin(hi - 1e-9); i++) {
    double x = data->GetBinCenter(i), y = data->GetBinContent(i), e = data->GetBinError(i);
    if (e <= 0) e = 1;                          // empty bin: Poisson-like error of 1
    std::vector<double> f(np);
    f[0] = TemplateAt(mc, x, delta, k, muMC);
    for (int j = 0; j <= npol; j++) f[1 + j] = std::pow(x, j);
    double w = 1 / (e * e);
    for (int a = 0; a < np; a++) { B(a) += w * f[a] * y; for (int b = 0; b < np; b++) A(a, b) += w * f[a] * f[b]; }
    rows.push_back(f); ys.push_back(y); ws.push_back(w);
  }
  // A singular system (e.g. a template that vanishes in the window) gets a huge chi2 so that the
  // scan ignores it; ROOT's error message is silenced for this check only.
  const int oldLevel = gErrorIgnoreLevel; gErrorIgnoreLevel = kFatal;
  double det = 0; A.Invert(&det);
  gErrorIgnoreLevel = oldLevel;
  if (det == 0 || rows.empty()) { par.assign(np, 0); ndf = 0; return 1e30; }
  TVectorD p = A * B;
  par.assign(p.GetMatrixArray(), p.GetMatrixArray() + np);
  double chi2 = 0;
  for (size_t r = 0; r < rows.size(); r++) {
    double model = 0; for (int a = 0; a < np; a++) model += p(a) * rows[r][a];
    chi2 += ws[r] * (ys[r] - model) * (ys[r] - model);
  }
  ndf = (int)rows.size() - np - 1;   // -1 for delta
  return chi2;
}

// Scan delta for fixed width scale k; best delta from a parabola through the chi2 minimum.
ScanResult ScanDelta(const TH1D *data, const TH1D *mc, double muMC, double k, int npol,
                     double lo, double hi, double center) {
  ScanResult r; r.stretch = k; r.lo = lo; r.hi = hi; r.npol = npol;
  int ndf = 0;
  for (double d = center - kScanHalfWidth; d <= center + kScanHalfWidth + 1e-9; d += kScanStep) {
    r.scanDelta.push_back(d);
    r.scanChi2.push_back(LinearChi2(data, mc, d, k, muMC, npol, lo, hi, r.par, ndf));
  }
  const int n = r.scanDelta.size();
  const int im = std::min_element(r.scanChi2.begin(), r.scanChi2.end()) - r.scanChi2.begin();
  TMatrixD M(3, 3); TVectorD V(3);
  for (int j = 0; j < n; j++) {
    if (r.scanChi2[j] > r.scanChi2[im] + kScanParabolaDchi2 && std::abs(j - im) > 3) continue;
    double x = r.scanDelta[j] - r.scanDelta[im], y = r.scanChi2[j];
    double xp[3] = {x * x, x, 1};
    for (int a = 0; a < 3; a++) { V(a) += xp[a] * y; for (int b = 0; b < 3; b++) M(a, b) += xp[a] * xp[b]; }
  }
  M.Invert(); TVectorD q = M * V;
  const bool atEdge = (im == 0 || im == n - 1);
  r.delta = (q(0) > 0) ? r.scanDelta[im] - q(1) / (2 * q(0)) : r.scanDelta[im];
  r.err   = (q(0) > 0 && !atEdge) ? 1 / std::sqrt(q(0)) : -1;
  r.chi2  = LinearChi2(data, mc, r.delta, k, muMC, npol, lo, hi, r.par, r.ndf);
  if (r.err > 0 && r.ndf > 0 && r.chi2 > r.ndf) r.err *= std::sqrt(r.chi2 / r.ndf);   // inflate for shape mismatch
  return r;
}

// Stretched template: also scan the MC width scale k; keep the k with the lowest chi2.
ScanResult ScanDeltaAndStretch(const TH1D *data, const TH1D *mc, double muMC, int npol,
                               double lo, double hi, double center) {
  ScanResult best; best.chi2 = 1e30;
  for (double k = kStretchMin; k <= kStretchMax + 1e-9; k += kStretchStep) {
    ScanResult r = ScanDelta(data, mc, muMC, k, npol, lo, hi, center);
    if (r.err > 0 && r.chi2 < best.chi2) best = r;
  }
  return best;
}

// Fit window [mu - a*sigma, mu + b*sigma] around the data peak, kept inside the histogram.
std::pair<double,double> FitWindow(const Peak &d, size_t i) {
  return {std::max(-kFitWindowEdge, d.mu - kFitWindows[i].first * d.sigma),
          std::min( kFitWindowEdge, d.mu + kFitWindows[i].second * d.sigma)};
}

// MC template shifted by delta, normalized to the data in the peak core (for display).
TH1D *ShiftedMC(const TH1D *data, const TH1D *mc, double delta, const Peak &d, int color) {
  TH1D *h = (TH1D*)data->Clone(Uid("shiftedMC").c_str()); h->Reset(); h->SetDirectory(0);
  for (int i = 1; i <= h->GetNbinsX(); i++) h->SetBinContent(i, TemplateAt(mc, h->GetBinCenter(i), delta, 1, 0));
  const double lo = d.mu - 1.5 * d.sigma, hi = d.mu + 1.5 * d.sigma;
  const double norm = h->Integral(h->FindBin(lo), h->FindBin(hi));
  if (norm > 0) h->Scale(data->Integral(data->GetXaxis()->FindFixBin(lo), data->GetXaxis()->FindFixBin(hi)) / norm);
  h->SetLineColor(color); h->SetLineWidth(2);
  return h;
}

// ---------------------------------------------------------------------------------------------
// Analysis of one spot
// ---------------------------------------------------------------------------------------------
SpotResult AnalyzeSpot(const std::string &tag, const SpotWindows &sw, const TH3D *d3, const TH3D *m3,
                       const std::vector<CutSet> &cuts, const std::vector<TH2D*> &d2, const std::vector<TH2D*> &m2) {
  SpotResult R; R.spot = sw.name;
  const int nz = d3->GetNbinsZ();
  // ---- nominal dy distributions
  R.hData = ProjectDy(d3, sw.dNom, 1, nz, tag + "_" + sw.name + "_dy_data");
  R.hMC   = ProjectDy(m3, sw.mNom, 1, nz, tag + "_" + sw.name + "_dy_mc");
  R.nData = R.hData->Integral();
  R.data = DyPeak(R.hData); R.mc = DyPeak(R.hMC);
  R.offCore = R.data.mu - R.mc.mu; R.errCore = std::hypot(R.data.emu, R.mc.emu);

  // ---- template scans over all fit windows
  double strMin = 1e9, strMax = -1e9;
  for (size_t iw = 0; iw < kFitWindows.size(); iw++) {
    auto [lo, hi] = FitWindow(R.data, iw);
    ScanResult plain = ScanDelta(R.hData, R.hMC, R.mc.mu, 1, kNominalPolOrder, lo, hi, R.offCore);
    ScanResult str   = ScanDeltaAndStretch(R.hData, R.hMC, R.mc.mu, kNominalPolOrder, lo, hi, R.offCore);
    if (iw == 0) { R.tpl = plain; R.str = str; }
    if (str.err > 0) { strMin = std::min(strMin, str.delta); strMax = std::max(strMax, str.delta); }
  }

  // ---- variations (core method)
  auto AddVar = [&](const std::string &name, TH1D *a, TH1D *b, bool inSys) {
    Peak pa = DyPeak(a), pb = DyPeak(b);
    R.varName.push_back(name); R.varOffset.push_back(pa.mu - pb.mu);
    R.varError.push_back(std::hypot(pa.emu, pb.emu)); R.varInSys.push_back(inSys);
    delete a; delete b;
  };
  for (size_t i = 1; i < cuts.size(); i++)
    AddVar(cuts[i].name, ProjectDy(d2[i], sw.dNom, "v_d"), ProjectDy(m2[i], sw.mNom, "v_m"), cuts[i].inSys);
  AddVar(Form("dx win %.1f#sigma", kSpotNsigmaVar[0]), ProjectDy(d3, sw.dNarrow, 1, nz, "v_d"), ProjectDy(m3, sw.mNarrow, 1, nz, "v_m"), true);
  AddVar(Form("dx win %.1f#sigma", kSpotNsigmaVar[1]), ProjectDy(d3, sw.dWide, 1, nz, "v_d"), ProjectDy(m3, sw.mWide, 1, nz, "v_m"), true);
  TH1D *hW2 = ProjectW2(d3, sw.dNom, "w2_data");
  auto [zc1, zc2] = QuantileBins(hW2, kW2CoreQlo, kW2CoreQhi);
  AddVar("central W2 half", ProjectDy(d3, sw.dNom, zc1, zc2, "v_d"), ProjectDy(m3, sw.mNom, zc1, zc2, "v_m"), true);
  for (size_t i = 0; i < R.varName.size(); i++)
    if (R.varInSys[i] && std::abs(R.varOffset[i] - R.offCore) > R.sysCut) { R.sysCut = std::abs(R.varOffset[i] - R.offCore); R.worstCut = R.varName[i]; }

  // ---- combine
  R.best = (R.str.err > 0) ? 0.5 * (R.offCore + R.str.delta) : R.offCore;
  R.stat = std::max(R.errCore, R.str.err);
  R.sysFit = (strMax > strMin) ? 0.5 * (strMax - strMin) : 0;
  R.sysMethod = 0.5 * (std::max({R.offCore, R.tpl.delta, R.str.delta}) - std::min({R.offCore, R.tpl.delta, R.str.delta}));
  R.sys = std::sqrt(R.sysCut * R.sysCut + R.sysFit * R.sysFit + R.sysMethod * R.sysMethod);
  const double total = std::hypot(R.stat, R.sys);
  R.verdict = "ok";
  if (std::abs(R.best) > kOkThreshold) R.verdict = (std::abs(R.best) > kFlagNsigma * total) ? "FLAG" : "check";

  // ---- W2 slices of equal data population (core method)
  std::vector<int> edges = {1};
  for (int k = 1; k < kNW2Slices; k++) {
    double p = double(k) / kNW2Slices, q; hW2->GetQuantiles(1, &q, &p);
    int b = hW2->FindFixBin(q);
    if (b > edges.back()) edges.push_back(b);
  }
  edges.push_back(nz + 1);
  for (size_t k = 0; k + 1 < edges.size(); k++) {
    const int z1 = edges[k], z2 = edges[k + 1] - 1;
    TH1D *a = ProjectDy(d3, sw.dNom, z1, z2, "s_d"), *b = ProjectDy(m3, sw.mNom, z1, z2, "s_m");
    if (a->Integral() >= kMinSliceEvents) {
      Peak pa = DyPeak(a), pb = DyPeak(b);
      hW2->GetXaxis()->SetRange(z1, z2);
      R.w2Mean.push_back(hW2->GetMean()); hW2->GetXaxis()->SetRange(0, 0);
      R.w2MuData.push_back(pa.mu); R.w2ErrData.push_back(pa.emu);
      R.w2MuMC.push_back(pb.mu);   R.w2ErrMC.push_back(pb.emu);
      R.w2Delta.push_back(pa.mu - pb.mu); R.w2Err.push_back(std::hypot(pa.emu, pb.emu));
    }
    delete a; delete b;
  }
  if (R.w2Delta.size() >= 3) {
    TGraphErrors g(R.w2Delta.size());
    for (size_t k = 0; k < R.w2Delta.size(); k++) { g.SetPoint(k, R.w2Mean[k], R.w2Delta[k]); g.SetPointError(k, 0, R.w2Err[k]); }
    TF1 f0("f0", "pol0"), f1("f1", "pol1");
    g.Fit(&f0, "Q0"); g.Fit(&f1, "Q0");
    R.constChi2 = f0.GetChisquare(); R.constNdf = f0.GetNDF();
    R.slope = f1.GetParameter(1); R.slopeErr = f1.GetParError(1);
  }
  delete hW2;
  return R;
}

// ---------------------------------------------------------------------------------------------
// Plots (appended to the open multi-page PDF)
// ---------------------------------------------------------------------------------------------
const int kSpotColor[3] = {kRed + 1, kBlue + 1, kBlack};

// Two pages per spot: dy overlay, template fit + chi2 scans.
void DrawSpot(const CaseResult &C, const SpotResult &R, const std::string &pdf) {
  const std::string title = C.in.name + " spot " + R.spot;
  const double tot = std::hypot(R.stat, R.sys);
  const double recHoff = (C.in.hoffFromFile ? C.in.hoffFile : C.in.hoffConfig) + R.best;
  // Page 1: data vs MC as produced (current hoffset) and shifted by Delta (linear and log)
  TCanvas c1("c1", "", kPageWidth, kPageHeight); c1.Divide(2, 1);
  for (int pad = 1; pad <= 2; pad++) {
    c1.cd(pad);
    TH1D *d = (TH1D*)R.hData->Clone(Uid("d").c_str()); d->SetDirectory(0);
    d->SetMarkerStyle(20); d->SetMarkerSize(0.6); d->SetLineColor(kBlack);
    d->SetTitle(Form("%s%s;dy (m);events", title.c_str(), pad == 2 ? " (log)" : ""));
    if (pad == 2) { gPad->SetLogy(); d->SetMinimum(0.5); }
    d->Draw("E");
    TH1D *h0 = ShiftedMC(d, R.hMC, 0, R.data, kBlue); h0->Draw("hist same");
    TH1D *hb = ShiftedMC(d, R.hMC, R.best, R.data, kRed); hb->Draw("hist same");
    TLegend *l = new TLegend(0.55, 0.66, 0.89, 0.89);
    l->AddEntry(d, Form("data: #mu %+.4f, #sigma %.3f", R.data.mu, R.data.sigma), "lep");
    l->AddEntry(h0, Form("MC (hoff %+.3f): #mu %+.4f, #sigma %.3f", C.in.hoffFile, R.mc.mu, R.mc.sigma), "l");
    l->AddEntry(hb, Form("MC shifted by #Delta = %+.4f #pm %.4f", R.best, tot), "l");
    l->AddEntry((TObject*)0, Form("recommended hoff %+.4f (%s)", recHoff, R.verdict.c_str()), "");
    l->Draw();
  }
  c1.Print(pdf.c_str());

  // Page 2: nominal template fit (left) and chi2(delta) of the plain and stretched scans (right)
  TCanvas c2("c2", "", kPageWidth, kPageHeight); c2.Divide(2, 1);
  c2.cd(1);
  TH1D *d = (TH1D*)R.hData->Clone(Uid("dfit").c_str()); d->SetDirectory(0);
  d->SetMarkerStyle(20); d->SetMarkerSize(0.6); d->SetLineColor(kBlack);
  d->SetTitle(Form("%s: N#upointMC(dy-#delta)+pol%d, #delta = %+.4f #pm %.4f, #chi^{2}/ndf = %.1f/%d;dy (m);events",
                   title.c_str(), R.tpl.npol, R.tpl.delta, R.tpl.err, R.tpl.chi2, R.tpl.ndf));
  d->Draw("E");
  TH1D *fitTot = (TH1D*)d->Clone(Uid("ftot").c_str()), *fitBg = (TH1D*)d->Clone(Uid("fbg").c_str());
  fitTot->Reset(); fitBg->Reset();
  if (!R.tpl.par.empty())
    for (int i = d->FindBin(R.tpl.lo + 1e-9); i <= d->FindBin(R.tpl.hi - 1e-9); i++) {
      double x = d->GetBinCenter(i), bg = 0;
      for (int j = 0; j <= R.tpl.npol; j++) bg += R.tpl.par[1 + j] * std::pow(x, j);
      fitBg->SetBinContent(i, bg);
      fitTot->SetBinContent(i, bg + R.tpl.par[0] * TemplateAt(R.hMC, x, R.tpl.delta, 1, 0));
    }
  fitTot->SetLineColor(kRed); fitTot->SetLineWidth(2); fitTot->Draw("hist same");
  fitBg->SetLineColor(kMagenta); fitBg->SetLineStyle(2); fitBg->Draw("hist same");
  c2.cd(2);
  auto Curve = [](const ScanResult &s, int color) {
    TGraph *g = new TGraph();
    if (s.scanChi2.empty()) return g;
    double cmin = *std::min_element(s.scanChi2.begin(), s.scanChi2.end());
    for (size_t j = 0; j < s.scanDelta.size(); j++) g->SetPoint(j, s.scanDelta[j], s.scanChi2[j] - cmin);
    g->SetLineColor(color); g->SetLineWidth(2); return g;
  };
  TGraph *gt = Curve(R.tpl, kBlack), *gs = Curve(R.str, kRed);
  double chi2Max = kMinDisplayDchi2;
  if (gt->GetN() > 0 && 0 > R.tpl.scanDelta.front() && 0 < R.tpl.scanDelta.back())
    chi2Max = std::max(chi2Max, 1.2 * std::max(gt->Eval(0), gs->GetN() ? gs->Eval(0) : 0.));
  gt->SetTitle(Form("%s: #chi^{2} scan;#delta = data - MC (m);#chi^{2} - #chi^{2}_{min}", title.c_str()));
  gt->SetMinimum(0); gt->SetMaximum(chi2Max); gt->Draw("AL"); if (gs->GetN()) gs->Draw("L");
  TLine *l0 = new TLine(0, 0, 0, chi2Max); l0->SetLineColor(kBlue); l0->SetLineStyle(2);
  if (gt->GetN() > 0 && 0 > R.tpl.scanDelta.front() && 0 < R.tpl.scanDelta.back()) l0->Draw();
  TLegend *l2 = new TLegend(0.50, 0.72, 0.89, 0.89);
  l2->AddEntry(gt, Form("plain template: %+.4f", R.tpl.delta), "l");
  l2->AddEntry(gs, Form("stretched (k=%.2f): %+.4f", R.str.stretch, R.str.delta), "l");
  l2->AddEntry(l0, "current MC (#delta = 0)", "l");
  l2->Draw();
  c2.Print(pdf.c_str());
}

// One page per case: dx spots and windows (left), Delta per variation for every spot (right).
void DrawVariations(const CaseResult &C, TH1D *dxData, TH1D *dxMC, const std::string &pdf) {
  TCanvas c("cv", "", kPageWidth, kPageHeight); c.Divide(2, 1);
  c.cd(1);
  dxData->SetTitle(Form("%s: dx with |dy| < %.2f m, spot windows #pm%.1f#sigma (data solid, MC dashed);dx (m);events",
                        C.in.name.c_str(), kDySpotFind, kSpotNsigma));
  dxData->SetMarkerStyle(20); dxData->SetMarkerSize(0.5); dxData->SetLineColor(kBlack); dxData->Draw("E");
  if (dxMC->Integral() > 0) dxMC->Scale(dxData->Integral() / dxMC->Integral());
  dxMC->SetLineColor(kGreen + 2); dxMC->SetLineWidth(2); dxMC->Draw("hist same");
  const double ymax = dxData->GetMaximum() * 1.05;
  for (int s = 0; s < (C.in.isLD2 ? 2 : 1); s++)
    for (int sign = -1; sign <= 1; sign += 2) {
      TLine *ld = new TLine(C.dataCenter[s] + sign * kSpotNsigma * C.mcSigma[s], 0, C.dataCenter[s] + sign * kSpotNsigma * C.mcSigma[s], ymax);
      TLine *lm = new TLine(C.mcCenter[s] + sign * kSpotNsigma * C.mcSigma[s], 0, C.mcCenter[s] + sign * kSpotNsigma * C.mcSigma[s], ymax);
      ld->SetLineColor(kSpotColor[s]); lm->SetLineColor(kSpotColor[s]); lm->SetLineStyle(2);
      ld->Draw(); lm->Draw();
    }
  c.cd(2);
  const int nv = C.spots[0].varName.size() + 1;   // + nominal
  TH1D *frame = new TH1D(Uid("frame").c_str(), Form("%s: #Delta = <dy>_{data} - <dy>_{MC} per variation (core Gaussian; * = not in sys.);;#Delta (m)", C.in.name.c_str()), nv, 0, nv);
  frame->SetDirectory(0);
  frame->GetXaxis()->SetBinLabel(1, "nominal");
  for (int i = 1; i < nv; i++)
    frame->GetXaxis()->SetBinLabel(i + 1, (C.spots[0].varName[i - 1] + (C.spots[0].varInSys[i - 1] ? "" : " *")).c_str());
  double ymin = 1e9, ymx = -1e9;
  for (auto &R : C.spots) {
    ymin = std::min(ymin, R.offCore); ymx = std::max(ymx, R.offCore);
    for (double v : R.varOffset) { ymin = std::min(ymin, v); ymx = std::max(ymx, v); }
  }
  frame->SetMinimum(std::min(ymin, -kOkThreshold) - 0.01); frame->SetMaximum(std::max(ymx, kOkThreshold) + 0.01);
  frame->Draw("AXIS");
  TLegend *l = new TLegend(0.12, 0.75, 0.45, 0.89);
  for (size_t s = 0; s < C.spots.size(); s++) {
    const SpotResult &R = C.spots[s];
    const double dxs = (double(s) - 0.5 * (C.spots.size() - 1)) * 0.15;
    TGraphErrors *g = new TGraphErrors();
    g->SetPoint(0, 0.5 + dxs, R.offCore); g->SetPointError(0, 0, R.errCore);
    for (size_t i = 0; i < R.varName.size(); i++) { g->SetPoint(i + 1, i + 1.5 + dxs, R.varOffset[i]); g->SetPointError(i + 1, 0, R.varError[i]); }
    const int col = C.in.isLD2 ? kSpotColor[s] : kBlack;
    g->SetMarkerStyle(20 + s); g->SetMarkerColor(col); g->SetLineColor(col); g->Draw("P");
    TLine *lb = new TLine(0, R.best, nv, R.best); lb->SetLineColor(col); lb->SetLineStyle(2); lb->Draw();
    l->AddEntry(g, Form("%s: best %+.4f #pm %.4f #pm %.4f", R.spot.c_str(), R.best, R.stat, R.sys), "p");
  }
  TLine *z = new TLine(0, 0, nv, 0); z->Draw();
  l->Draw();
  c.Print(pdf.c_str());
}

// One page per case: core dy peaks vs W2 (left) and Delta vs W2 (right).
void DrawW2Slices(const CaseResult &C, const std::string &pdf) {
  TCanvas c("cw", "", kPageWidth, kPageHeight); c.Divide(2, 1);
  const SpotResult &A = C.spots.back();   // "all" (LD2) or "p" (LH2)
  c.cd(1);
  TGraphErrors *gd = new TGraphErrors(), *gm = new TGraphErrors();
  for (size_t k = 0; k < A.w2Mean.size(); k++) {
    gd->SetPoint(k, A.w2Mean[k], A.w2MuData[k]); gd->SetPointError(k, 0, A.w2ErrData[k]);
    gm->SetPoint(k, A.w2Mean[k], A.w2MuMC[k]);   gm->SetPointError(k, 0, A.w2ErrMC[k]);
  }
  gd->SetTitle(Form("%s spot %s: dy core peak vs W^{2};W^{2} (GeV^{2});<dy> (m)", C.in.name.c_str(), A.spot.c_str()));
  gd->SetMarkerStyle(20); gm->SetMarkerStyle(24); gm->SetMarkerColor(kBlue); gm->SetLineColor(kBlue);
  TGraph range; for (size_t k = 0; k < A.w2Mean.size(); k++) { range.SetPoint(2 * k, A.w2Mean[k], A.w2MuData[k]); range.SetPoint(2 * k + 1, A.w2Mean[k], A.w2MuMC[k]); }
  if (gd->GetN() > 0) {
    gd->SetMinimum(TMath::MinElement(range.GetN(), range.GetY()) - 0.02); gd->SetMaximum(TMath::MaxElement(range.GetN(), range.GetY()) + 0.02);
    gd->Draw("AP"); gm->Draw("P");
  }
  TLegend *l1 = new TLegend(0.60, 0.78, 0.89, 0.89); l1->AddEntry(gd, "data", "p"); l1->AddEntry(gm, "MC (current hoff)", "p"); l1->Draw();
  c.cd(2);
  TH1D *frame = new TH1D(Uid("w2frame").c_str(), Form("%s: #Delta vs W^{2} (core Gaussian);W^{2} (GeV^{2});#Delta = <dy>_{data} - <dy>_{MC} (m)", C.in.name.c_str()), 10, kW2Min, kW2Max);
  frame->SetDirectory(0);
  double lo = 1e9, hi = -1e9, ymin = 1e9, ymax = -1e9;
  for (auto &R : C.spots) for (size_t k = 0; k < R.w2Mean.size(); k++) {
    lo = std::min(lo, R.w2Mean[k]); hi = std::max(hi, R.w2Mean[k]);
    ymin = std::min(ymin, R.w2Delta[k] - R.w2Err[k]); ymax = std::max(ymax, R.w2Delta[k] + R.w2Err[k]);
  }
  if (lo < hi) frame->GetXaxis()->SetRangeUser(lo - 0.1, hi + 0.1);
  frame->SetMinimum(std::min(ymin, -kOkThreshold) - 0.005); frame->SetMaximum(std::max(ymax, kOkThreshold) + 0.005);
  frame->Draw("AXIS");
  TLegend *l2 = new TLegend(0.12, 0.75, 0.60, 0.89);
  for (size_t s = 0; s < C.spots.size(); s++) {
    const SpotResult &R = C.spots[s];
    const int col = C.in.isLD2 ? kSpotColor[s] : kBlack;
    TGraphErrors *g = new TGraphErrors();
    for (size_t k = 0; k < R.w2Mean.size(); k++) { g->SetPoint(k, R.w2Mean[k] + 0.005 * s, R.w2Delta[k]); g->SetPointError(k, 0, R.w2Err[k]); }
    g->SetMarkerStyle(20 + s); g->SetMarkerColor(col); g->SetLineColor(col); g->Draw("P");
    if (R.w2Mean.size() >= 3) {
      TF1 *f = new TF1(Uid("lin").c_str(), "pol1", lo - 0.05, hi + 0.05);
      g->Fit(f, "Q0"); f->SetLineColor(col); f->SetLineStyle(2); f->Draw("same");
    }
    l2->AddEntry(g, Form("%s: slope %+.4f #pm %.4f m/GeV^{2}, const #chi^{2}/ndf %.1f/%d", R.spot.c_str(), R.slope, R.slopeErr, R.constChi2, R.constNdf), "p");
  }
  TLine *z = new TLine(lo - 0.1, 0, hi + 0.1, 0); z->Draw();
  l2->Draw();
  c.Print(pdf.c_str());
}

// Text page(s).
void DrawText(const std::vector<std::string> &lines, const std::string &pdf, const char *title) {
  const size_t kLinesPerPage = 30;
  for (size_t start = 0; start < lines.size(); start += kLinesPerPage) {
    TCanvas c("ct", "", kPageWidth, kPageHeight);
    TPaveText *pt = new TPaveText(0.02, 0.02, 0.98, 0.98, "NDC");
    pt->SetTextAlign(12); pt->SetTextFont(42); pt->SetTextSize(0.024); pt->SetFillColor(0); pt->SetBorderSize(0);
    pt->AddText(title);
    for (size_t i = start; i < std::min(lines.size(), start + kLinesPerPage); i++) pt->AddText(lines[i].c_str());
    pt->Draw();
    c.Print(pdf.c_str());
  }
}

// Summary pages: Delta of every (setting, target, spot); recommended vs current hoffset per setting.
void DrawSummary(const std::vector<CaseResult> &all, const std::string &pdf) {
  std::vector<std::string> label; std::vector<double> val, err;
  for (auto &C : all) for (auto &R : C.spots) {
    label.push_back(C.in.setting + " " + C.in.target + " " + R.spot);
    val.push_back(R.best); err.push_back(std::hypot(R.stat, R.sys));
  }
  const int n = label.size();
  if (n == 0) return;
  TCanvas c("csum", "", kPageWidth, kPageHeight);
  c.SetBottomMargin(0.25);
  TH1D *frame = new TH1D("summary_frame", "#Delta = <dy>_{data} - <dy>_{MC} (MC with its current hcal_hoffset);;#Delta (m)", n, 0, n);
  frame->SetDirectory(0);
  for (int i = 0; i < n; i++) frame->GetXaxis()->SetBinLabel(i + 1, label[i].c_str());
  frame->GetXaxis()->LabelsOption("v");
  double ymin = *std::min_element(val.begin(), val.end()), ymax = *std::max_element(val.begin(), val.end());
  frame->SetMinimum(std::min(ymin, -0.01) - 0.01); frame->SetMaximum(std::max(ymax, 0.01) + 0.01); frame->Draw("AXIS");
  TGraphErrors *g = new TGraphErrors();
  for (int i = 0; i < n; i++) { g->SetPoint(i, i + 0.5, val[i]); g->SetPointError(i, 0, err[i]); }
  g->SetMarkerStyle(20); g->Draw("P");
  TLine *l0 = new TLine(0, 0, n, 0); l0->Draw();
  TLine *lp = new TLine(0,  kOkThreshold, n,  kOkThreshold); lp->SetLineStyle(2); lp->Draw();
  TLine *lm = new TLine(0, -kOkThreshold, n, -kOkThreshold); lm->SetLineStyle(2); lm->Draw();
  c.Print(pdf.c_str());

  // recommended hoffset per setting: LH2 (p) and LD2 (all), with the current MC values
  std::vector<std::string> sets;
  for (auto &C : all) if (std::find(sets.begin(), sets.end(), C.in.setting) == sets.end()) sets.push_back(C.in.setting);
  const int ns = sets.size();
  TCanvas c2("csum2", "", kPageWidth, kPageHeight);
  TH1D *f2 = new TH1D("summary_frame2", "MC hcal_hoffset: recommended (filled, #pm stat#oplussys) and current (open);;hcal_hoffset (m)", ns, 0, ns);
  f2->SetDirectory(0);
  for (int i = 0; i < ns; i++) f2->GetXaxis()->SetBinLabel(i + 1, sets[i].c_str());
  TGraphErrors *gr[2] = {new TGraphErrors(), new TGraphErrors()}; TGraph *gc[2] = {new TGraph(), new TGraph()};
  double lo = 1e9, hi = -1e9;
  for (auto &C : all) {
    const int i = std::find(sets.begin(), sets.end(), C.in.setting) - sets.begin();
    const int t = C.in.isLD2 ? 1 : 0;
    const SpotResult &R = C.spots.back();
    const double cur = C.in.hoffFromFile ? C.in.hoffFile : C.in.hoffConfig;
    const double rec = cur + R.best;
    gr[t]->SetPoint(gr[t]->GetN(), i + 0.4 + 0.2 * t, rec); gr[t]->SetPointError(gr[t]->GetN() - 1, 0, std::hypot(R.stat, R.sys));
    gc[t]->SetPoint(gc[t]->GetN(), i + 0.4 + 0.2 * t, cur);
    lo = std::min({lo, rec, cur}); hi = std::max({hi, rec, cur});
  }
  f2->SetMinimum(lo - 0.02); f2->SetMaximum(hi + 0.02); f2->Draw("AXIS");
  const int col[2] = {kRed + 1, kBlue + 1};
  TLegend *l = new TLegend(0.12, 0.75, 0.50, 0.89);
  for (int t = 0; t < 2; t++) {
    if (gr[t]->GetN() == 0) continue;   // no block of this target in the run
    gr[t]->SetMarkerStyle(20 + t); gr[t]->SetMarkerColor(col[t]); gr[t]->SetLineColor(col[t]); gr[t]->Draw("P");
    gc[t]->SetMarkerStyle(24 + t); gc[t]->SetMarkerColor(col[t]); gc[t]->Draw("P");
    l->AddEntry(gr[t], t ? "LD2 (p+n spots), recommended" : "LH2, recommended", "p");
    l->AddEntry(gc[t], t ? "LD2 MC, current" : "LH2 MC, current", "p");
  }
  l->Draw();
  c2.Print(pdf.c_str());
}

// ---------------------------------------------------------------------------------------------
// Inputs of one config block
// ---------------------------------------------------------------------------------------------
// LH2: files and strict elastic cuts from the w2_offset_check config, plus the W2 window.
CaseInput SetupLH2(JSONManager *j, JSONManager *w2, const std::string &key,
                   const std::string &dataDir, const std::string &simuDir) {
  CaseInput in; in.name = key; in.target = "LH2"; in.isLD2 = false;
  in.setting = j->GetValueFromSubKey_str(key, "setting");
  const std::string wkey = j->GetValueFromSubKey_str(key, "w2_key");
  in.mcWeight = j->GetValueFromSubKey_str("common", "lh2_mc_weight");
  in.dataFile = dataDir + "/" + w2->GetValueFromSubKey_str(wkey, "data_file");
  std::vector<std::string> files; w2->GetVectorFromSubKey<std::string>(wkey, "simu_files", files);
  for (auto &f : files) in.simuFiles.push_back(simuDir + "/" + f);
  w2->GetVectorFromSubKey<double>(wkey, "simu_weights", in.simuWeights);
  const std::string runCut = w2->GetValueFromSubKey_str(wkey, "data_run_cut");
  const double eHCALmin = w2->GetValueFromSubKey<double>(wkey, "eHCAL_min");
  const double thpqCut = w2->GetValueFromSubKey<double>("common", "thpq_p_cut");
  const double vzCut = w2->GetValueFromSubKey<double>("common", "vz_cut");
  const double coinTWindow = w2->GetValueFromSubKey<double>("common", "coinT_window");
  std::vector<double> w2win; j->GetVectorFromSubKey<double>(key, "W2_window", w2win);

  // coincidence-time peak of elastic events (data), as in w2_offset_check
  // (TTreeFormula does not know "true", the run cut of most settings)
  const std::string runSel = (runCut == "true") ? "1" : runCut;
  const std::string elasSel = Form("(%s)&&W2>%g&&W2<%g&&thpq_p<%g&&eHCAL>%g&&bbfiduCut", runSel.c_str(), kElasW2Lo, kElasW2Hi, thpqCut, eHCALmin);
  const double coinTPeak = ReadCoinTPeak(in.dataFile, elasSel);
  printf("  %s: coinT peak %.2f ns\n", key.c_str(), coinTPeak);
  auto Coin = [&](double w) { return std::string(Form("abs(coinT_ADC_c-(%g))<%g", coinTPeak, w)); };

  auto Base = [&](double eScale, bool useVz) {
    std::string c = Form("eHCAL>%g&&bbfiduCut&&W2>%g&&W2<%g", eScale * eHCALmin, w2win[0], w2win[1]);
    if (useVz) c += Form("&&abs(vz)<%g", vzCut);
    return c;
  };
  auto Data = [&](const std::string &base, double w) { return "(" + runCut + ")&&" + base + "&&" + Coin(w); };
  in.cuts.push_back({"nominal", Data(Base(1, true), coinTWindow), Base(1, true), true});
  in.cuts.push_back({"thpq_p (w2 check)", Data(Base(1, true) + Form("&&thpq_p<%g", thpqCut), coinTWindow), Base(1, true) + Form("&&thpq_p<%g", thpqCut), false});
  in.cuts.push_back({"no vz", Data(Base(1, false), coinTWindow), Base(1, false), true});
  in.cuts.push_back({"eHCAL x2", Data(Base(2, true), coinTWindow), Base(2, true), true});
  in.cuts.push_back({"eHCAL x0.5", Data(Base(0.5, true), coinTWindow), Base(0.5, true), true});
  in.cuts.push_back({"ARCut", Data(Base(1, true) + "&&ARCut", coinTWindow), Base(1, true) + "&&ARCut", true});
  for (double w : kCoinTWindowVariations)
    in.cuts.push_back({Form("coinT #pm%.1f ns", w), Data(Base(1, true), w), Base(1, true), true});

  // current hoffset: config value of every MC file's block; all must agree
  JSONManager mcj(j->GetValueFromSubKey_str(key, "mc_config").c_str());
  std::vector<std::string> blocks; j->GetVectorFromSubKey<std::string>(key, "mc_config_blocks", blocks);
  for (size_t b = 0; b < blocks.size(); b++) {
    double h = mcj.GetValueFromSubKey<double>(blocks[b], "hcal_hoffset");
    if (b == 0) in.hoffConfig = h;
    else if (std::abs(h - in.hoffConfig) > 1e-6) printf("WARNING %s: MC config blocks have different hcal_hoffset\n", key.c_str());
  }
  return in;
}

// LD2: production cuts, files and weights exactly as fit_dx.cpp builds them ("qelas" block).
CaseInput SetupLD2(JSONManager *j, const std::string &key, const std::string &dataDir, const std::string &simuDir) {
  CaseInput in; in.name = key; in.target = "LD2"; in.isLD2 = true;
  in.setting = j->GetValueFromSubKey_str(key, "setting");
  JSONManager fd(j->GetValueFromSubKey_str(key, "fit_dx_config").c_str());
  const char *fk = "qelas";   // fit_dx reads the "qelas" block for LD2 (fit_dx.cpp: key = is_elastic ? "elas" : "qelas")
  const int conf = fd.GetValueFromSubKey<int>(fk, "SBS_config"), sbsmag = fd.GetValueFromSubKey<int>(fk, "SBS_magnet_percent");
  const int model = fd.GetValueFromSubKey<int>(fk, "model"), pass = fd.GetValueFromSubKey<int>(fk, "pass");
  const std::string gen = fd.GetValueFromSubKey_str(fk, "generator");
  std::string dfprefix = fd.GetValueFromSubKey_str(fk, "data_file_prefix");
  std::string sfprefix = fd.GetValueFromSubKey_str(fk, "simu_file_prefix");
  const std::string RCtype = fd.GetValueFromSubKey_str(fk, "RC_prescription");
  sfprefix = Form(sfprefix.c_str(), RCtype.c_str());
  dfprefix = dfprefix.empty() ? "" : dfprefix + "_";
  sfprefix = sfprefix.empty() ? "" : sfprefix + "_";
  in.dataFile = dataDir + "/" + Form("%s%s_ana_data_sbs%d_sbs%dp_model%d_pass%d.root", dfprefix.c_str(), fk, conf, sbsmag, model, pass);
  in.simuFiles = {simuDir + "/" + Form("%s%s_ana_%s_sbs%d_sbs%dp_model%d.root", sfprefix.c_str(), fk, gen.c_str(), conf, sbsmag, model)};
  in.simuWeights = {1.0};
  // MC weight branch (fit_dx: weightn)
  const bool useEffi = fd.GetValueFromSubKey<int>(fk, "use_effi_corr"), useNorm = fd.GetValueFromSubKey<int>(fk, "use_norm_weight");
  in.mcWeight = useEffi ? (useNorm ? "weight_effic_norm" : "weight_effic") : (useNorm ? "weight_norm" : "weight");
  // MC -> data dx shifts used by fit_dx (dx_shifted = dx + offset)
  in.dxOffP = fd.GetValueFromSubKey<double>(fk, "dx_offset_MC_for_p");
  in.dxOffN = fd.GetValueFromSubKey<double>(fk, "dx_offset_MC_for_n");
  // custom fiducial cut (fit_dx: use_custom_fiduCut)
  in.customFidu = fd.GetValueFromSubKey<int>(fk, "use_custom_fiduCut");
  if (in.customFidu) {
    std::vector<double> ARw, SMw; fd.GetVectorFromSubKey<double>(fk, "AR_width_x_y", ARw); fd.GetVectorFromSubKey<double>(fk, "SM_width_xp_xn_y", SMw);
    in.sbsKick = fd.GetValueFromSubKey<double>(fk, "sbs_kick");
    in.fiduAR = cut::hcal_active_area_data(ARw[0], ARw[1], pass);
    in.fiduSM = cut::hcal_safety_margin(SMw[0], SMw[1], SMw[2], in.fiduAR);
  }
  // cuts: data = simu && extra (&& nblkHCAL>1) ; the dy term is removed
  const std::string simu = RemoveTerms(fd.GetValueFromSubKey_str(fk, "cuts_for_signal_simu"), "dy", in.removedTerms);
  const std::string extra = fd.GetValueFromSubKey_str(fk, "extra_data_signal_cut");
  const bool nblk = fd.GetValueFromSubKey<int>(fk, "apply_nblkHCAL>1_cut");
  const std::string fidu = in.customFidu ? "&&fidu_custom" : "";
  const std::string nb = nblk ? "&&nblkHCAL>1" : "";
  const std::string mcBase = simu + fidu, dataBase = simu + "&&" + extra + fidu;
  in.cuts.push_back({"nominal", dataBase + nb, mcBase + nb, true});
  in.cuts.push_back({Form("|vz|<%.2f", kVzTight), dataBase + nb + Form("&&abs(vz)<%g", kVzTight), mcBase + nb + Form("&&abs(vz)<%g", kVzTight), true});
  const double eHard = j->GetValueFromSubKey<double>(key, "eHCAL_hard");
  in.cuts.push_back({Form("eHCAL>%.3g", eHard), dataBase + nb + Form("&&eHCAL>%g", eHard), mcBase + nb + Form("&&eHCAL>%g", eHard), true});
  if (nblk) in.cuts.push_back({"no nblkHCAL>1", dataBase, mcBase, true});
  std::vector<std::string> dummy;
  in.cuts.push_back({"no ARCut", RemoveTerms(dataBase + nb, "ARCut", dummy), RemoveTerms(mcBase + nb, "ARCut", dummy), false});

  JSONManager mcj(j->GetValueFromSubKey_str(key, "mc_config").c_str());
  in.hoffConfig = mcj.GetValueFromSubKey<double>(j->GetValueFromSubKey_str(key, "mc_config_block"), "hcal_hoffset");
  return in;
}

// ---------------------------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------------------------
int dy_offset_check(const char *configfilename = "configs/conf_dy_offset_check.json", std::string only = "") {
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);

  // ---- config
  JSONManager *jmgr = new JSONManager(configfilename);
  const std::string dataDir = jmgr->GetValueFromSubKey_str("common", "data_dir");
  const std::string simuDir = jmgr->GetValueFromSubKey_str("common", "simu_dir");
  const std::string outDir  = jmgr->GetValueFromSubKey_str("common", "out_dir");
  std::string outBase       = jmgr->GetValueFromSubKey_str("common", "out_file_base");
  JSONManager *w2cfg = new JSONManager(jmgr->GetValueFromSubKey_str("common", "w2_offset_config").c_str());
  std::vector<std::string> settings; jmgr->GetVectorFromKey_str("settings", settings);
  if (!only.empty()) { settings = {only}; outBase += "_" + only; }

  // ---- inputs of every block (single-threaded file reads first, as in w2_offset_check)
  std::vector<CaseInput> inputs;
  for (auto &key : settings) {
    printf("==== setting up %s\n", key.c_str());
    const std::string target = jmgr->GetValueFromSubKey_str(key, "target");
    CaseInput in = (target == "LD2") ? SetupLD2(jmgr, key, dataDir, simuDir) : SetupLH2(jmgr, w2cfg, key, dataDir, simuDir);
    if (in.simuFiles.size() != in.simuWeights.size()) { printf("ERROR %s: simu files and weights differ in length\n", key.c_str()); continue; }
    // hoffset (and dy_offset) the MC files were made with, from their summary canvas
    for (size_t f = 0; f < in.simuFiles.size(); f++) {
      double h = 0, dyo = 0;
      bool found = ReadSummaryOffsets(in.simuFiles[f], h, dyo);
      if (f == 0) { in.hoffFromFile = found; in.hoffFile = h; in.dyOffFile = dyo; }
      else if (found != in.hoffFromFile || std::abs(h - in.hoffFile) > 1e-6) printf("WARNING %s: MC files have different hcal_hoffset\n", key.c_str());
    }
    if (!in.hoffFromFile) { in.hoffFile = in.hoffConfig; printf("WARNING %s: hoffset not found in the MC file, assuming the config value %.4f\n", key.c_str(), in.hoffConfig); }
    else if (std::abs(in.hoffFile - in.hoffConfig) > 1e-6) printf("WARNING %s: MC file hoffset %.4f differs from the config value %.4f\n", key.c_str(), in.hoffFile, in.hoffConfig);
    if (std::abs(in.dyOffFile) > 1e-9) printf("WARNING %s: MC dy_offset = %.4f is included in the stored dy\n", key.c_str(), in.dyOffFile);
    double dy0 = 0; ReadSummaryOffsets(in.dataFile, in.hoffData, dy0);
    printf("  data %s\n  MC   %s (weight %s)\n  hoffset: MC file %.4f, MC config %.4f, data %.4f\n",
           in.dataFile.c_str(), in.simuFiles[0].c_str(), in.mcWeight.c_str(), in.hoffFile, in.hoffConfig, in.hoffData);
    for (auto &t : in.removedTerms) printf("  removed cut term: %s\n", t.c_str());
    inputs.push_back(in);
  }
  ROOT::EnableImplicitMT();

  mkdir(outDir.c_str(), 0755);
  const std::string pdf = outDir + "/" + outBase + ".pdf";
  TFile fout((outDir + "/" + outBase + ".root").c_str(), "RECREATE");
  FILE *csv = fopen((outDir + "/" + outBase + ".csv").c_str(), "w");
  fprintf(csv, "setting,target,spot,hoffset_mc_file,hoffset_mc_config,hoffset_data,dy_data,dy_data_err,dy_mc,dy_mc_err,sigma_data,sigma_mc,"
               "core,core_err,template,template_err,stretched,stretched_err,stretch_k,delta,stat,sys,sys_cut,sys_fit,sys_method,worst_cut,"
               "ref_variation_shift,recommended_hoffset,w2_slope,w2_slope_err,w2_const_chi2,w2_const_ndf,n_data,verdict\n");

  TCanvas cOpen("cOpen", "", kPageWidth, kPageHeight); cOpen.Print((pdf + "[").c_str());
  DrawText({"Data/MC dy peak alignment and the MC hcal_hoffset (configs/conf_dy_offset_check.json).",
            "dy = yHCAL - yHCAL_exp. In the MC, HCAL_origin includes hoff * y_HCAL, so yHCAL_exp(hoff) = yHCAL_exp(0) - hoff and MC dy(hoff) = dy(0) + hoff exactly",
            "   (simulation/elas_ana_simu.cpp:348,573; qelas_ana_simu.cpp:368,647; src/KinematicVar.cpp:122). Stored dy = dy + dy_offset (0 in all files used).",
            "#Delta = <dy>_{data} - <dy>_{MC}  =>  recommended MC hcal_hoffset = hoff(file) + #Delta. hoff(file) is read from the MC file's summary canvas.",
            "LH2: strict elastic cuts of w2_offset_check (eHCAL, bbfiduCut, |vz|, coinT, run cut) + W2 window, with a dx window on the proton spot",
            "   instead of thpq_p (thpq_p contains dy and truncates it around 0; shown as a reference variation, not in the systematic error).",
            "LD2: production cuts of the fit_dx qelas block (same files, weight, nblkHCAL>1, custom fiducial cut where used), WITHOUT the |dy| cut term.",
            Form("   Spots: dx windows #pm%.1f#sigma around each spot (data and MC each around their own centre); 'all' = union of the p and n windows.", kSpotNsigma),
            "Methods: core Gaussian and template shift scan (MC(dy-#delta) + pol2 background); best = mean of core and stretched template.",
            "Systematics: physics-cut variations, dx window width, central W2 half, fit windows, method spread. Reference variations (*) not included.",
            Form("Verdict: ok if |#Delta| < %.3f m; otherwise FLAG if |#Delta| > %.0f#sigma_{tot}, else check.", kOkThreshold, kFlagNsigma),
            "W2 dependence: #Delta in W2 slices of equal data population (core Gaussian) with a straight-line fit.",
            "MC hoffset changes also move yHCAL_exp, so fiducial cuts on yHCAL_exp (in the fit_dx cut strings) select slightly different MC events: second order.",
            "CAVEAT: the LD2 data files were made before the ARCut fix in qelas_ana_data (ARCut evaluated at the wrong HCAL cluster, 0.3-7% of events);",
            "   the 'no ARCut' reference variation shows the sensitivity of the dy peak to ARCut."}, pdf, "dy_offset_check: notes");

  std::vector<CaseResult> all;
  std::vector<std::string> tableLines;
  for (auto &in : inputs) {
    printf("==== %s\n", in.name.c_str());
    CaseResult C; C.in = in; C.nSpots = in.isLD2 ? 2 : 1;
    const std::string &key = in.name;
    const ROOT::RDF::TH3DModel m3model(Form("%s_h3", key.c_str()), "", kNbinsDy, kDyMin, kDyMax, kNbinsDx, kDxMin, kDxMax, kNbinsW2, kW2Min, kW2Max);
    const ROOT::RDF::TH2DModel m2model(Form("%s_h2", key.c_str()), "", kNbinsDy, kDyMin, kDyMax, kNbinsDx, kDxMin, kDxMax);
    const ROOT::RDF::TH1DModel m1model(Form("%s_dx", key.c_str()), "", kNbinsDx, kDxMin, kDxMax);

    // the custom fiducial cut of fit_dx, as a column (only defined when used)
    auto AddFidu = [&](ROOT::RDF::RNode node) -> ROOT::RDF::RNode {
      if (!in.customFidu) return node;
      std::vector<double> AR = in.fiduAR, SM = in.fiduSM; double kick = in.sbsKick;
      return node.Define("fidu_custom", [AR, SM, kick](double x, double y, double xe, double ye) {
          return cut::inHCAL_activeA(x, y, AR) && cut::inHCAL_safety_margin("LD2", xe, ye, kick, SM); },
        {"xHCAL", "yHCAL", "xHCAL_exp", "yHCAL_exp"});
    };

    // ---- book all histograms (data and every MC file), then run the event loops together
    std::vector<ROOT::RDF::RResultHandle> handles;
    ROOT::RDataFrame dfData("Tout", in.dataFile);
    ROOT::RDF::RNode dData = AddFidu(dfData);
    auto d3 = dData.Filter(in.cuts[0].dataCut).Histo3D(m3model, "dy", "dx", "W2");
    handles.push_back(d3);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> d2r(in.cuts.size());
    for (size_t i = 1; i < in.cuts.size(); i++) { d2r[i] = dData.Filter(in.cuts[i].dataCut).Histo2D(m2model, "dy", "dx"); handles.push_back(d2r[i]); }

    std::vector<std::unique_ptr<ROOT::RDataFrame>> dfMC;
    std::vector<ROOT::RDF::RResultPtr<TH3D>> m3r;
    std::vector<std::vector<ROOT::RDF::RResultPtr<TH2D>>> m2r;
    std::vector<ROOT::RDF::RResultPtr<TH1D>> mdxp, mdxn;
    for (size_t f = 0; f < in.simuFiles.size(); f++) {
      dfMC.emplace_back(new ROOT::RDataFrame("Tout", in.simuFiles[f]));
      ROOT::RDF::RNode dMC = AddFidu(*dfMC.back()).Define("w_dy", Form("%s*%.10f", in.mcWeight.c_str(), in.simuWeights[f]));
      auto nomMC = dMC.Filter(in.cuts[0].mcCut);
      m3r.push_back(nomMC.Histo3D(m3model, "dy", "dx", "W2", "w_dy")); handles.push_back(m3r.back());
      if (in.isLD2) {   // MC proton / neutron dx spots separately (mc_fnucl: 1 = p, 0 = n)
        mdxp.push_back(nomMC.Filter("mc_fnucl==1").Histo1D(m1model, "dx", "w_dy")); handles.push_back(mdxp.back());
        mdxn.push_back(nomMC.Filter("mc_fnucl==0").Histo1D(m1model, "dx", "w_dy")); handles.push_back(mdxn.back());
      }
      m2r.emplace_back(in.cuts.size());
      for (size_t i = 1; i < in.cuts.size(); i++) { m2r[f][i] = dMC.Filter(in.cuts[i].mcCut).Histo2D(m2model, "dy", "dx", "w_dy"); handles.push_back(m2r[f][i]); }
    }
    ROOT::RDF::RunGraphs(handles);

    // ---- sum the MC files
    TH3D *h3d = (TH3D*)d3->Clone(Uid(key + "_data3").c_str()); h3d->SetDirectory(0);
    TH3D *h3m = (TH3D*)m3r[0]->Clone(Uid(key + "_mc3").c_str()); h3m->SetDirectory(0);
    for (size_t f = 1; f < m3r.size(); f++) h3m->Add(m3r[f].GetPtr());
    std::vector<TH2D*> h2d(in.cuts.size(), nullptr), h2m(in.cuts.size(), nullptr);
    for (size_t i = 1; i < in.cuts.size(); i++) {
      h2d[i] = (TH2D*)d2r[i]->Clone(Uid(key + "_data2").c_str()); h2d[i]->SetDirectory(0);
      h2m[i] = (TH2D*)m2r[0][i]->Clone(Uid(key + "_mc2").c_str()); h2m[i]->SetDirectory(0);
      for (size_t f = 1; f < m2r.size(); f++) h2m[i]->Add(m2r[f][i].GetPtr());
    }
    printf("  events after nominal cuts: data %.0f, MC (weighted) %.4g\n", h3d->Integral(), h3m->Integral());

    // ---- dx spots: centre and width from core Gaussians (MC first, then data around the expectation)
    TH1D *dxData = ProjectDx(h3d, key + "_dx_data"), *dxMC = ProjectDx(h3m, key + "_dx_mc");
    if (in.isLD2) {
      TH1D *hp = (TH1D*)mdxp[0]->Clone(Uid("mdxp").c_str()), *hn = (TH1D*)mdxn[0]->Clone(Uid("mdxn").c_str());
      hp->SetDirectory(0); hn->SetDirectory(0);
      for (size_t f = 1; f < mdxp.size(); f++) { hp->Add(mdxp[f].GetPtr()); hn->Add(mdxn[f].GetPtr()); }
      Peak pp = CorePeak(hp, kDxMin, kDxMax, kDxCoreSigmaStart), pn = CorePeak(hn, kDxMin, kDxMax, kDxCoreSigmaStart);
      C.mcCenter[0] = pp.mu; C.mcSigma[0] = pp.sigma; C.mcCenter[1] = pn.mu; C.mcSigma[1] = pn.sigma;
      const double off[2] = {in.dxOffP, in.dxOffN};
      for (int s = 0; s < 2; s++) {
        const double expect = C.mcCenter[s] + off[s];   // fit_dx compares data with MC dx + offset
        C.dataCenter[s] = CorePeak(dxData, expect - kDxSpotSearch, expect + kDxSpotSearch, C.mcSigma[s]).mu;
      }
      delete hp; delete hn;
    } else {
      Peak pm = CorePeak(dxMC, kDxMin, kDxMax, kDxCoreSigmaStart);
      C.mcCenter[0] = pm.mu; C.mcSigma[0] = pm.sigma;
      C.dataCenter[0] = CorePeak(dxData, pm.mu - kDxSpotSearch, pm.mu + kDxSpotSearch, pm.sigma).mu;
    }
    for (int s = 0; s < C.nSpots; s++)
      printf("  spot %s: dx centre data %+.3f, MC %+.3f, MC sigma %.3f m\n", s ? "n" : "p", C.dataCenter[s], C.mcCenter[s], C.mcSigma[s]);

    // ---- spot windows (data and MC each around their own centre, MC width)
    const TAxis *ax = h3d->GetYaxis();
    auto Win = [&](int s, bool isData, double nsig) {
      const double c = isData ? C.dataCenter[s] : C.mcCenter[s];
      return std::make_pair(c - nsig * C.mcSigma[s], c + nsig * C.mcSigma[s]);
    };
    auto MakeSpot = [&](const std::string &name, std::vector<int> spotIdx) {
      SpotWindows w; w.name = name;
      auto R = [&](bool isData, double nsig) {
        std::vector<std::pair<double,double>> v; for (int s : spotIdx) v.push_back(Win(s, isData, nsig));
        return BinRanges(ax, v);
      };
      w.dNom = R(true, kSpotNsigma); w.mNom = R(false, kSpotNsigma);
      w.dNarrow = R(true, kSpotNsigmaVar[0]); w.mNarrow = R(false, kSpotNsigmaVar[0]);
      w.dWide = R(true, kSpotNsigmaVar[1]); w.mWide = R(false, kSpotNsigmaVar[1]);
      return w;
    };
    std::vector<SpotWindows> spots;
    if (in.isLD2) spots = {MakeSpot("p", {0}), MakeSpot("n", {1}), MakeSpot("all", {0, 1})};
    else          spots = {MakeSpot("p", {0})};

    // ---- analysis per spot
    const double cur = in.hoffFromFile ? in.hoffFile : in.hoffConfig;
    for (auto &sw : spots) {
      SpotResult R = AnalyzeSpot(key, sw, h3d, h3m, in.cuts, h2d, h2m);
      double refShift = 0;   // shift of the reference (not in sys.) variation, for the CSV
      for (size_t i = 0; i < R.varName.size(); i++) if (!R.varInSys[i]) refShift = R.varOffset[i] - R.offCore;
      printf("  RESULT %s %s: dy data %+.4f, MC %+.4f | core %+.4f | template %+.4f | stretched %+.4f (k %.2f) | Delta %+.4f +- %.4f (stat) +- %.4f (sys) [%s] | hoff %.4f -> %.4f | W2 slope %+.4f +- %.4f -> %s\n",
             in.target.c_str(), R.spot.c_str(), R.data.mu, R.mc.mu, R.offCore, R.tpl.delta, R.str.delta, R.str.stretch,
             R.best, R.stat, R.sys, R.worstCut.c_str(), cur, cur + R.best, R.slope, R.slopeErr, R.verdict.c_str());
      fprintf(csv, "%s,%s,%s,%.4f,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.4f,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%s,%.5f,%.4f,%.5f,%.5f,%.1f,%d,%.0f,%s\n",
              in.setting.c_str(), in.target.c_str(), R.spot.c_str(), in.hoffFile, in.hoffConfig, in.hoffData,
              R.data.mu, R.data.emu, R.mc.mu, R.mc.emu, R.data.sigma, R.mc.sigma,
              R.offCore, R.errCore, R.tpl.delta, R.tpl.err, R.str.delta, R.str.err, R.str.stretch,
              R.best, R.stat, R.sys, R.sysCut, R.sysFit, R.sysMethod, R.worstCut.c_str(), refShift,
              cur + R.best, R.slope, R.slopeErr, R.constChi2, R.constNdf, R.nData, R.verdict.c_str());
      tableLines.push_back(Form("%-10s %s %-3s  hoff %+.3f  <dy> data %+.4f  MC %+.4f  #Delta = %+.4f #pm %.4f #pm %.4f  =>  hoff %+.4f  [%s]  dDelta/dW2 = %+.4f #pm %.4f",
                                in.setting.c_str(), in.target.c_str(), R.spot.c_str(), cur, R.data.mu, R.mc.mu, R.best, R.stat, R.sys,
                                cur + R.best, R.verdict.c_str(), R.slope, R.slopeErr));
      // ROOT output: nominal dy histograms and the W2-slice graph
      fout.cd();
      R.hData->Write(); R.hMC->Write();
      TGraphErrors gw(R.w2Mean.size());
      for (size_t k = 0; k < R.w2Mean.size(); k++) { gw.SetPoint(k, R.w2Mean[k], R.w2Delta[k]); gw.SetPointError(k, 0, R.w2Err[k]); }
      gw.Write(Form("%s_%s_delta_vs_W2", key.c_str(), R.spot.c_str()));
      C.spots.push_back(R);
    }

    // ---- plots
    for (auto &R : C.spots) DrawSpot(C, R, pdf);
    DrawVariations(C, dxData, dxMC, pdf);
    DrawW2Slices(C, pdf);
    delete h3d; delete h3m; delete dxData; delete dxMC;
    for (size_t i = 1; i < in.cuts.size(); i++) { delete h2d[i]; delete h2m[i]; }
    all.push_back(C);
  }
  DrawSummary(all, pdf);
  DrawText(tableLines, pdf, "Summary: #Delta = <dy>_{data} - <dy>_{MC} (stat, sys) and recommended MC hcal_hoffset = hoff + #Delta (m)");
  cOpen.Print((pdf + "]").c_str());
  fclose(csv);
  fout.Close();
  printf("Output: %s/%s.{pdf,root,csv}\n", outDir.c_str(), outBase.c_str());
  return 0;
}
