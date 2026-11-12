/*
  w2_offset_check.cpp : check the MC W2_offset of every kinematic setting with LH2 elastics.

  Physics / method
  ----------------
  The MC analysis (simulation/elas_ana_simu.cpp) writes W2 = W2_recon + W2_offset, where W2_offset
  is a constant taken from the config ("baked" into the file). The offset is meant to absorb the
  small data/MC difference in the reconstructed electron kinematics, so that the elastic W2 peak
  of the MC sits where the data peak is.

  To check it we
    1. remove the baked offset from the MC (W2_0 = W2 - baked). The offset enters nowhere else in
       the MC analysis, so this is exactly the analysis with W2_offset = 0;
    2. select strict LH2 elastics with identical cuts on data and MC:
         thpq_p < cut   (proton direction predicted from the electron matches the HCAL cluster)
         eHCAL > min    (HCAL cluster energy)
         bbfiduCut      (BigBite fiducial region)
         |vz| < cut     (target vertex)
       plus, for data only, a window around the e-HCAL coincidence-time peak;
    3. measure the offset needed for the MC, i.e. the shift that moves the zero-offset MC peak
       onto the data peak, with two independent methods:
         (a) core Gaussian : iterative Gaussian fit to the peak core of data and MC separately,
                             offset = mu_data - mu_MC;
         (b) template scan : data = N * MC_0(W2 - delta) + polynomial, scanned in delta. For a
                             fixed delta the model is linear in N and the polynomial coefficients,
                             so it is solved exactly by weighted least squares (no MINUIT). The best
                             delta comes from a parabola through chi2(delta); its error from
                             delta-chi2 = 1. The "stretched" version also scales the MC width by a
                             factor k about its peak, because the MC peak is wider than the data at
                             some settings and a width mismatch would otherwise pull the shift.
    4. best value = mean of (a) and the stretched template (b);
       systematic   = largest shift from cut variations (+) spread over fit windows (+) spread
                      between methods, added in quadrature.

  The result is compared with the baked offset: "ok" if they differ by less than kFlagThreshold,
  otherwise "FLAG" when the difference is also significant (> kFlagNsigma sigma) and "check" when not.

  Inputs : configs/conf_w2_offset_check.json (one block per setting; see "settings" list there).
  Outputs (in "out_dir", default simulation/siout/w2_offset_check):
    <out_file_base>.pdf  - all plots, a few pages per setting plus a summary page
    <out_file_base>.root - the W2 histograms used
    <out_file_base>.csv  - one row of results per setting

  Usage (from the scripts/ directory):
    root -l -b -q 'w2_offset_check.cpp("configs/conf_w2_offset_check.json")'
    root -l -b -q 'w2_offset_check.cpp("configs/conf_w2_offset_check.json","sbs8_70p")'   // one setting
  The output file names of a single-setting run get the setting name appended.
  -------
  P. Datta, 2026-11-12 (method developed in the 2026-11-12 W2_offset study)
*/

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>
#include <algorithm>
#include <sys/stat.h>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMatrixD.h"
#include "TVectorD.h"
#include "ROOT/RDataFrame.hxx"

#include "../dflay/src/JSONManager.cxx"

// ---------------------------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------------------------
// W2 histograms: fine 0.005 GeV2 bins (MC template is interpolated on these); data are rebinned
// by kDataRebin for the template fits.
const int    kNbinsW2 = 600;   const double kW2Min = -0.5, kW2Max = 2.5;   // GeV2
const int    kDataRebin = 2;                                             // -> 0.01 GeV2 data bins
// Coincidence time histogram (data) used to find the e-HCAL coincidence peak.
const int    kNbinsCoinT = 80; const double kCoinTMin = -20, kCoinTMax = 20; // ns
const double kCoinTPeakSearch = 3.0;   // ns, half-width around the maximum used to compute the peak mean
// W2 window of elastic events used to locate the coinT peak (GeV2).
const double kElasW2Lo = 0.6, kElasW2Hi = 1.3;
// Number of entries used to read back the baked offset (it is constant within a file).
const Long64_t kNentriesForBaked = 100000;
const double   kWminForBaked = 0.1;   // W > 0.1 GeV: avoid W = 0 entries (W2_recon < 0)

// Core Gaussian: peak searched in [kCoreSearchLo, kCoreSearchHi], fit range +- kCoreNsigma sigma.
const double kCoreSearchLo = 0.3, kCoreSearchHi = 1.6;   // GeV2
const double kCoreNsigma = 1.2;
const double kCoreSigmaStart = 0.08;                     // GeV2, starting width
const int    kCoreIterations = 6;

// Template scan: delta scanned +- kScanHalfWidth around the core result in steps of kScanStep.
const double kScanStep = 0.0005, kScanHalfWidth = 0.06;  // GeV2
const double kScanParabolaDchi2 = 4;   // points within chi2_min + this are used for the parabola
const double kFitWindowMinW2 = 0.2;    // GeV2, lowest allowed edge of a fit window
// Fit windows in units of the data core sigma: [mu - lo*sigma, mu + hi*sigma]. The first is nominal.
const std::vector<std::pair<double,double>> kFitWindows = {{3,3}, {2,2}, {4,4}, {3,2}, {2.5,4}};
const int    kNominalPolOrder = 2;     // background polynomial order of the nominal template fit
// Width-scale scan of the stretched template.
const double kStretchMin = 0.70, kStretchMax = 1.15, kStretchStep = 0.01;

// Systematic variations (data and MC get the same physics cut; coinT windows are data only).
struct CutVariation { std::string name; double thpq_scale; bool use_vz; double eHCAL_scale; };
const std::vector<CutVariation> kCutVariations = {
  {"nominal",   1.0, true,  1.0},
  {"thpq_x2",   2.0, true,  1.0},   // looser proton direction match
  {"thpq_x0.5", 0.5, true,  1.0},   // tighter
  {"no_vz",     1.0, false, 1.0},   // no target vertex cut
  {"eHCAL_x2",  1.0, true,  2.0},   // harder HCAL energy cut
  {"eHCAL_x0.5",1.0, true,  0.5}};  // softer
const std::vector<double> kCoinTWindowVariations = {1.5, 6.0};   // ns

// Verdict on |best - baked|: "ok" below kFlagThreshold (GeV2, about twice the typical total
// uncertainty of the well-measured settings); above it, "FLAG" if also more than
// kFlagNsigma total standard deviations, otherwise "check" (large but not significant).
const double kFlagThreshold = 0.005;
const double kFlagNsigma = 2.0;

// All PDF pages use one canvas size, so every page is laid out the same way.
const int kPageWidth = 1500, kPageHeight = 650;   // pixels
// chi2 scan display: y range is at least this, enlarged so that the baked offset is visible.
const double kMinDisplayDchi2 = 25;

// ---------------------------------------------------------------------------------------------
// Small containers
// ---------------------------------------------------------------------------------------------
struct Peak { double mu, emu, sigma; };

struct ScanResult {
  double delta = 0, err = -1;        // best shift and its error (err < 0: minimum not found)
  double chi2 = 0; int ndf = 0;
  double stretch = 1;                // MC width scale (1 for the plain template)
  double lo = 0, hi = 0; int npol = 0;
  std::vector<double> par;           // N, polynomial coefficients
  std::vector<double> scanDelta, scanChi2;
};

struct SettingResult {
  std::string name;
  double baked = 0;
  Peak data{}, mc{};
  double offCore = 0, errCore = 0;
  ScanResult tpl, str;
  double best = 0, stat = 0, sysCut = 0, sysFit = 0, sysMethod = 0, sys = 0;
  std::string worstCut;
  std::vector<std::string> varName; std::vector<double> varOffset, varError;
};

// ---------------------------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------------------------

// Offset baked into an MC file: W2 = W*W + offset for every event (W is written without it).
double ReadBakedOffset(const std::string &file) {
  TFile f(file.c_str());
  TTree *t = (TTree*)f.Get("Tout");
  if (!t) { printf("ERROR: no Tout in %s\n", file.c_str()); return 0; }
  Long64_t n = t->Draw("W2-W*W", Form("W>%g", kWminForBaked), "goff", kNentriesForBaked);
  double sum = 0; for (Long64_t i = 0; i < n; i++) sum += t->GetV1()[i];
  return n > 0 ? sum / n : 0;
}

// Selection string shared by data and MC for one cut variation.
std::string ElasticCut(const CutVariation &v, double thpqCut, double vzCut, double eHCALmin) {
  std::string cut = Form("thpq_p<%g&&eHCAL>%g&&bbfiduCut", v.thpq_scale * thpqCut, v.eHCAL_scale * eHCALmin);
  if (v.use_vz) cut += Form("&&abs(vz)<%g", vzCut);
  return cut;
}

// Iterative Gaussian fit to the peak core (+-kCoreNsigma sigma around the current mean).
Peak CorePeak(TH1D *h) {
  h->GetXaxis()->SetRangeUser(kCoreSearchLo, kCoreSearchHi);
  double mu = h->GetBinCenter(h->GetMaximumBin()), sigma = kCoreSigmaStart;
  h->GetXaxis()->SetRange(0, 0);
  TF1 g("gcore", "gaus");
  for (int i = 0; i < kCoreIterations; i++) {
    g.SetRange(mu - kCoreNsigma * sigma, mu + kCoreNsigma * sigma);
    g.SetParameters(h->GetBinContent(h->FindBin(mu)), mu, sigma);
    h->Fit(&g, "QR0");
    mu = g.GetParameter(1); sigma = std::abs(g.GetParameter(2));
  }
  return {mu, g.GetParError(1), sigma};
}

// Data W2 histogram with |coinT - coinT_peak| < window, from the 2D (W2, coinT) histogram.
TH1D *ProjectCoinTWindow(TH2D *h2, double peak, double window, const std::string &name) {
  int b1 = h2->GetYaxis()->FindBin(peak - window + 1e-6), b2 = h2->GetYaxis()->FindBin(peak + window - 1e-6);
  TH1D *h = h2->ProjectionX(name.c_str(), b1, b2, "e");
  h->SetDirectory(0);
  return h;
}

// MC template value at W2 = x for shift delta and width scale k (stretched about the MC peak muMC).
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
  A.Invert(); TVectorD p = A * B;
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
  // least-squares parabola y = a x^2 + b x + c (x measured from the minimum point)
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
  // inflate the error when the shapes do not agree well (chi2/ndf > 1)
  if (r.err > 0 && r.ndf > 0 && r.chi2 > r.ndf) r.err *= std::sqrt(r.chi2 / r.ndf);
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

// Fit window [mu - a*sigma, mu + b*sigma] around the data peak.
std::pair<double,double> FitWindow(const Peak &d, size_t i) {
  return {std::max(kFitWindowMinW2, d.mu - kFitWindows[i].first * d.sigma), d.mu + kFitWindows[i].second * d.sigma};
}

// MC template shifted by delta, normalized to the data in the peak core (for display).
TH1D *ShiftedMC(const TH1D *dataBinning, const TH1D *mc, double delta, const Peak &d, const char *name, int color) {
  TH1D *h = (TH1D*)dataBinning->Clone(name); h->Reset(); h->SetDirectory(0);
  for (int i = 1; i <= h->GetNbinsX(); i++) h->SetBinContent(i, TemplateAt(mc, h->GetBinCenter(i), delta, 1, 0));
  const double lo = d.mu - 1.5 * d.sigma, hi = d.mu + 1.5 * d.sigma;
  const double norm = h->Integral(h->FindBin(lo), h->FindBin(hi));
  if (norm > 0) h->Scale(dataBinning->Integral(dataBinning->GetXaxis()->FindFixBin(lo), dataBinning->GetXaxis()->FindFixBin(hi)) / norm);
  h->SetLineColor(color); h->SetLineWidth(2);
  return h;
}

// ---------------------------------------------------------------------------------------------
// Plots for one setting (appended to the open multi-page PDF)
// ---------------------------------------------------------------------------------------------
void DrawSetting(const SettingResult &R, TH1D *data, const TH1D *mc, const std::string &pdf) {
  const double tot = std::hypot(R.stat, R.sys);
  // Page 1: data vs zero-offset MC shifted by the baked and by the best offset (linear and log)
  TCanvas c1("c1", "", kPageWidth, kPageHeight); c1.Divide(2, 1);
  for (int pad = 1; pad <= 2; pad++) {
    c1.cd(pad);
    TH1D *d = (TH1D*)data->Clone(Form("%s_d%d", R.name.c_str(), pad)); d->SetDirectory(0);
    d->SetMarkerStyle(20); d->SetMarkerSize(0.6); d->SetLineColor(kBlack);
    d->SetTitle(Form("%s LH2 elastics%s;W^{2} (GeV^{2});events", R.name.c_str(), pad == 2 ? " (log)" : ""));
    if (pad == 1) d->GetXaxis()->SetRangeUser(R.data.mu - 5 * R.data.sigma, R.data.mu + 5 * R.data.sigma);
    else { d->GetXaxis()->SetRangeUser(-0.2, 2.0); gPad->SetLogy(); d->SetMinimum(0.5); }
    d->Draw("E");
    TH1D *hb = ShiftedMC(d, mc, R.baked, R.data, Form("%s_mcb%d", R.name.c_str(), pad), kBlue); hb->Draw("hist same");
    TH1D *hr = ShiftedMC(d, mc, R.best, R.data, Form("%s_mcr%d", R.name.c_str(), pad), kRed);   hr->Draw("hist same");
    TLegend *l = new TLegend(0.55, 0.70, 0.89, 0.89);
    l->AddEntry(d, Form("data (#sigma_{core} %.3f)", R.data.sigma), "lep");
    l->AddEntry(hb, Form("MC, baked offset %+.4f", R.baked), "l");
    l->AddEntry(hr, Form("MC, best %+.4f #pm %.4f", R.best, tot), "l");
    l->AddEntry((TObject*)0, Form("MC #sigma_{core} %.3f", R.mc.sigma), "");
    l->Draw();
  }
  c1.Print(pdf.c_str());

  // Page 2: nominal template fit (left) and chi2(delta) of the plain and stretched scans (right)
  TCanvas c2("c2", "", kPageWidth, kPageHeight); c2.Divide(2, 1);
  c2.cd(1);
  TH1D *d = (TH1D*)data->Clone(Form("%s_dfit", R.name.c_str())); d->SetDirectory(0);
  d->SetMarkerStyle(20); d->SetMarkerSize(0.6); d->SetLineColor(kBlack);
  d->GetXaxis()->SetRangeUser(R.data.mu - 5 * R.data.sigma, R.data.mu + 5 * R.data.sigma);
  d->SetTitle(Form("%s: N#upointMC_{0}(W^{2}-#delta)+pol%d, #delta = %+.4f #pm %.4f, #chi^{2}/ndf = %.1f/%d;W^{2} (GeV^{2});events",
                   R.name.c_str(), R.tpl.npol, R.tpl.delta, R.tpl.err, R.tpl.chi2, R.tpl.ndf));
  d->Draw("E");
  TH1D *fitTot = (TH1D*)d->Clone(Form("%s_ftot", R.name.c_str())), *fitBg = (TH1D*)d->Clone(Form("%s_fbg", R.name.c_str()));
  fitTot->Reset(); fitBg->Reset();
  for (int i = d->FindBin(R.tpl.lo + 1e-9); i <= d->FindBin(R.tpl.hi - 1e-9); i++) {
    double x = d->GetBinCenter(i), bg = 0;
    for (int j = 0; j <= R.tpl.npol; j++) bg += R.tpl.par[1 + j] * std::pow(x, j);
    fitBg->SetBinContent(i, bg);
    fitTot->SetBinContent(i, bg + R.tpl.par[0] * TemplateAt(mc, x, R.tpl.delta, 1, 0));
  }
  fitTot->SetLineColor(kRed); fitTot->SetLineWidth(2); fitTot->Draw("hist same");
  fitBg->SetLineColor(kMagenta); fitBg->SetLineStyle(2); fitBg->Draw("hist same");
  c2.cd(2);
  auto Curve = [](const ScanResult &s, int color) {
    TGraph *g = new TGraph();
    double cmin = *std::min_element(s.scanChi2.begin(), s.scanChi2.end());
    for (size_t j = 0; j < s.scanDelta.size(); j++) g->SetPoint(j, s.scanDelta[j], s.scanChi2[j] - cmin);
    g->SetLineColor(color); g->SetLineWidth(2); return g;
  };
  TGraph *gt = Curve(R.tpl, kBlack), *gs = Curve(R.str, kRed);
  // y range: show the curves up to the baked offset (if it is inside the scan), at least kMinDisplayDchi2
  double chi2Max = kMinDisplayDchi2;
  if (R.baked > R.tpl.scanDelta.front() && R.baked < R.tpl.scanDelta.back())
    chi2Max = std::max(chi2Max, 1.2 * std::max(gt->Eval(R.baked), gs->Eval(R.baked)));
  gt->SetTitle(Form("%s: #chi^{2} scan;#delta (GeV^{2});#chi^{2} - #chi^{2}_{min}", R.name.c_str()));
  gt->SetMinimum(0); gt->SetMaximum(chi2Max); gt->Draw("AL"); gs->Draw("L");
  TLine *lb = new TLine(R.baked, 0, R.baked, chi2Max); lb->SetLineColor(kBlue); lb->SetLineStyle(2); lb->Draw();
  TLegend *l2 = new TLegend(0.50, 0.72, 0.89, 0.89);
  l2->AddEntry(gt, Form("plain template: %+.4f", R.tpl.delta), "l");
  l2->AddEntry(gs, Form("stretched (k=%.2f): %+.4f", R.str.stretch, R.str.delta), "l");
  l2->AddEntry(lb, Form("baked %+.4f", R.baked), "l");
  l2->Draw();
  c2.Print(pdf.c_str());

  // Page 3: offset for each cut variation (core method), with the best value and its band
  TCanvas c3("c3", "", kPageWidth, kPageHeight);
  const int nv = R.varName.size();
  TH1D *frame = new TH1D(Form("%s_frame", R.name.c_str()), Form("%s: core-Gaussian offset per variation;;offset (GeV^{2})", R.name.c_str()), nv, 0, nv);
  frame->SetDirectory(0);
  for (int i = 0; i < nv; i++) frame->GetXaxis()->SetBinLabel(i + 1, R.varName[i].c_str());
  double ymin = std::min(R.best, R.baked) - 0.03, ymax = std::max(R.best, R.baked) + 0.03;
  frame->SetMinimum(ymin); frame->SetMaximum(ymax); frame->Draw("AXIS");
  TGraphErrors *gv = new TGraphErrors();
  for (int i = 0; i < nv; i++) { gv->SetPoint(i, i + 0.5, R.varOffset[i]); gv->SetPointError(i, 0, R.varError[i]); }
  gv->SetMarkerStyle(20); gv->Draw("P");
  TLine *lbest = new TLine(0, R.best, nv, R.best); lbest->SetLineColor(kRed); lbest->Draw();
  TLine *lbak  = new TLine(0, R.baked, nv, R.baked); lbak->SetLineColor(kBlue); lbak->SetLineStyle(2); lbak->Draw();
  TLegend *l3 = new TLegend(0.60, 0.78, 0.89, 0.89);
  l3->AddEntry(lbest, Form("best %+.4f", R.best), "l"); l3->AddEntry(lbak, Form("baked %+.4f", R.baked), "l"); l3->Draw();
  c3.Print(pdf.c_str());
}

// Summary page: best - baked per setting, with total error.
void DrawSummary(const std::vector<SettingResult> &all, const std::string &pdf) {
  TCanvas c("csum", "", kPageWidth, kPageHeight);
  const int n = all.size();
  TH1D *frame = new TH1D("summary_frame", "Needed minus baked W2_offset (LH2);;#delta_{best} - #delta_{baked} (GeV^{2})", n, 0, n);
  frame->SetDirectory(0);
  for (int i = 0; i < n; i++) frame->GetXaxis()->SetBinLabel(i + 1, all[i].name.c_str());
  frame->SetMinimum(-0.03); frame->SetMaximum(0.03); frame->Draw("AXIS");
  TGraphErrors *g = new TGraphErrors();
  for (int i = 0; i < n; i++) { g->SetPoint(i, i + 0.5, all[i].best - all[i].baked); g->SetPointError(i, 0, std::hypot(all[i].stat, all[i].sys)); }
  g->SetMarkerStyle(20); g->Draw("P");
  TLine *l0 = new TLine(0, 0, n, 0); l0->Draw();
  TLine *lp = new TLine(0,  kFlagThreshold, n,  kFlagThreshold); lp->SetLineStyle(2); lp->Draw();
  TLine *lm = new TLine(0, -kFlagThreshold, n, -kFlagThreshold); lm->SetLineStyle(2); lm->Draw();
  c.Print(pdf.c_str());
}

// ---------------------------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------------------------
int w2_offset_check(const char *configfilename = "configs/conf_w2_offset_check.json", std::string only = "") {
  gErrorIgnoreLevel = kError;
  gStyle->SetOptStat(0);

  // ---- config
  JSONManager *jmgr = new JSONManager(configfilename);
  const std::string dataDir = jmgr->GetValueFromSubKey_str("common", "data_dir");
  const std::string simuDir = jmgr->GetValueFromSubKey_str("common", "simu_dir");
  const std::string outDir  = jmgr->GetValueFromSubKey_str("common", "out_dir");
  std::string outBase       = jmgr->GetValueFromSubKey_str("common", "out_file_base");
  const double thpqCut      = jmgr->GetValueFromSubKey<double>("common", "thpq_p_cut");
  const double vzCut        = jmgr->GetValueFromSubKey<double>("common", "vz_cut");
  const double coinTWindow  = jmgr->GetValueFromSubKey<double>("common", "coinT_window");
  std::vector<std::string> settings; jmgr->GetVectorFromKey_str("settings", settings);
  if (!only.empty()) { settings = {only}; outBase += "_" + only; }

  mkdir(outDir.c_str(), 0755);
  const std::string pdf = outDir + "/" + outBase + ".pdf";
  TFile fout((outDir + "/" + outBase + ".root").c_str(), "RECREATE");
  FILE *csv = fopen((outDir + "/" + outBase + ".csv").c_str(), "w");
  fprintf(csv, "setting,baked,core,core_err,template,template_err,stretched,stretched_err,stretch_k,"
               "best,stat,sys,sys_cut,sys_fit,sys_method,worst_cut,best_minus_baked,sigma_data,sigma_mc,verdict\n");

  // ---- read back the baked offsets first: TTree::Draw with an entry limit, before multithreading
  std::vector<std::vector<double>> baked(settings.size());
  for (size_t s = 0; s < settings.size(); s++) {
    std::vector<std::string> files; jmgr->GetVectorFromSubKey<std::string>(settings[s], "simu_files", files);
    for (auto &f : files) baked[s].push_back(ReadBakedOffset(simuDir + "/" + f));
  }
  ROOT::EnableImplicitMT();

  TCanvas cOpen("cOpen", "", kPageWidth, kPageHeight); cOpen.Print((pdf + "[").c_str());
  std::vector<SettingResult> all;
  for (size_t s = 0; s < settings.size(); s++) {
    SettingResult R; R.name = settings[s];
    const std::string &key = settings[s];
    const std::string dataFile = dataDir + "/" + jmgr->GetValueFromSubKey_str(key, "data_file");
    const std::string runCut   = jmgr->GetValueFromSubKey_str(key, "data_run_cut");
    const double eHCALmin      = jmgr->GetValueFromSubKey<double>(key, "eHCAL_min");
    std::vector<std::string> simuFiles; jmgr->GetVectorFromSubKey<std::string>(key, "simu_files", simuFiles);
    std::vector<double> simuWeights;    jmgr->GetVectorFromSubKey<double>(key, "simu_weights", simuWeights);
    if (simuFiles.size() != simuWeights.size()) { printf("ERROR %s: simu_files and simu_weights differ in length\n", key.c_str()); continue; }
    // All MC files of a mix must carry the same offset, otherwise "baked" is ambiguous.
    R.baked = baked[s][0];
    for (double b : baked[s]) if (std::abs(b - R.baked) > 1e-4) printf("WARNING %s: MC files have different baked offsets\n", key.c_str());
    printf("==== %s: baked offset %+.4f\n", key.c_str(), R.baked);

    // ---- fill histograms for all cut variations (one pass per file)
    ROOT::RDataFrame dfData("Tout", dataFile);
    auto dData = dfData.Filter(runCut);
    std::vector<ROOT::RDF::RResultPtr<TH2D>> hData;
    for (auto &v : kCutVariations)
      hData.push_back(dData.Filter(ElasticCut(v, thpqCut, vzCut, eHCALmin))
                      .Histo2D({Form("%s_data_%s", key.c_str(), v.name.c_str()), "", kNbinsW2, kW2Min, kW2Max, kNbinsCoinT, kCoinTMin, kCoinTMax},
                               "W2", "coinT_ADC_c"));
    std::vector<TH1D*> hMC(kCutVariations.size(), nullptr);
    for (size_t f = 0; f < simuFiles.size(); f++) {
      ROOT::RDataFrame dfMC("Tout", simuDir + "/" + simuFiles[f]);
      // zero-offset W2, weighted with the per-event MC weight times the mix weight
      auto dMC = dfMC.Define("W2_0", Form("W2-(%.10f)", baked[s][f])).Define("w_mix", Form("weight*%.10f", simuWeights[f]));
      std::vector<ROOT::RDF::RResultPtr<TH1D>> h;
      for (auto &v : kCutVariations)
        h.push_back(dMC.Filter(ElasticCut(v, thpqCut, vzCut, eHCALmin))
                    .Histo1D({Form("%s_mc%zu_%s", key.c_str(), f, v.name.c_str()), "", kNbinsW2, kW2Min, kW2Max}, "W2_0", "w_mix"));
      for (size_t i = 0; i < h.size(); i++) {
        if (!hMC[i]) { hMC[i] = (TH1D*)h[i]->Clone(Form("%s_mc_%s", key.c_str(), kCutVariations[i].name.c_str())); hMC[i]->SetDirectory(0); }
        else hMC[i]->Add(h[i].GetPtr());
      }
    }

    // ---- coinT peak of elastic events (nominal cuts)
    TH2D *d2 = hData[0].GetPtr();
    TH1D *hct = d2->ProjectionY(Form("%s_coinT", key.c_str()), d2->GetXaxis()->FindBin(kElasW2Lo), d2->GetXaxis()->FindBin(kElasW2Hi));
    double cmax = hct->GetBinCenter(hct->GetMaximumBin());
    hct->GetXaxis()->SetRangeUser(cmax - kCoinTPeakSearch, cmax + kCoinTPeakSearch);
    const double coinTPeak = hct->GetMean();
    hct->GetXaxis()->SetRange(0, 0);

    // ---- (a) core Gaussian, nominal cuts
    TH1D *hd = ProjectCoinTWindow(d2, coinTPeak, coinTWindow, key + "_data_nominal_W2");
    TH1D *hm = hMC[0];
    R.data = CorePeak(hd); R.mc = CorePeak(hm);
    R.offCore = R.data.mu - R.mc.mu; R.errCore = std::hypot(R.data.emu, R.mc.emu);

    // ---- (b) template scans on coarser data bins, over all fit windows
    TH1D *hdR = (TH1D*)hd->Clone((key + "_data_nominal_W2_rebinned").c_str()); hdR->SetDirectory(0); hdR->Rebin(kDataRebin);
    double strMin = 1e9, strMax = -1e9;
    for (size_t iw = 0; iw < kFitWindows.size(); iw++) {
      auto [lo, hi] = FitWindow(R.data, iw);
      ScanResult plain = ScanDelta(hdR, hm, R.mc.mu, 1, kNominalPolOrder, lo, hi, R.offCore);
      ScanResult str   = ScanDeltaAndStretch(hdR, hm, R.mc.mu, kNominalPolOrder, lo, hi, R.offCore);
      printf("  window [%.3f,%.3f]: template %+.4f +- %.4f (chi2/ndf %.1f/%d) | stretched %+.4f +- %.4f k %.2f\n",
             lo, hi, plain.delta, plain.err, plain.chi2, plain.ndf, str.delta, str.err, str.stretch);
      if (iw == 0) { R.tpl = plain; R.str = str; }
      if (str.err > 0) { strMin = std::min(strMin, str.delta); strMax = std::max(strMax, str.delta); }
    }

    // ---- cut variations (core method) and coinT windows
    for (size_t i = 0; i < kCutVariations.size(); i++) {
      TH1D *a = ProjectCoinTWindow(hData[i].GetPtr(), coinTPeak, coinTWindow, key + "_data_" + kCutVariations[i].name + "_W2");
      Peak pa = CorePeak(a), pb = CorePeak(hMC[i]);
      R.varName.push_back(kCutVariations[i].name); R.varOffset.push_back(pa.mu - pb.mu); R.varError.push_back(std::hypot(pa.emu, pb.emu));
      fout.cd(); a->Write(); hMC[i]->Write();
    }
    for (double w : kCoinTWindowVariations) {
      TH1D *a = ProjectCoinTWindow(d2, coinTPeak, w, key + Form("_data_coinT%.1f_W2", w));
      Peak pa = CorePeak(a);
      R.varName.push_back(Form("coinT #pm%.1f ns", w)); R.varOffset.push_back(pa.mu - R.mc.mu); R.varError.push_back(std::hypot(pa.emu, R.mc.emu));
    }
    for (size_t i = 1; i < R.varName.size(); i++)
      if (std::abs(R.varOffset[i] - R.offCore) > R.sysCut) { R.sysCut = std::abs(R.varOffset[i] - R.offCore); R.worstCut = R.varName[i]; }

    // ---- combine
    R.best = 0.5 * (R.offCore + R.str.delta);
    R.stat = std::max(R.errCore, R.str.err);
    R.sysFit = (strMax > strMin) ? 0.5 * (strMax - strMin) : 0;
    R.sysMethod = 0.5 * (std::max({R.offCore, R.tpl.delta, R.str.delta}) - std::min({R.offCore, R.tpl.delta, R.str.delta}));
    R.sys = std::sqrt(R.sysCut * R.sysCut + R.sysFit * R.sysFit + R.sysMethod * R.sysMethod);
    const double diff = R.best - R.baked;
    const double total = std::hypot(R.stat, R.sys);
    std::string verdict = "ok";
    if (std::abs(diff) > kFlagThreshold) verdict = (std::abs(diff) > kFlagNsigma * total) ? "FLAG" : "check";
    printf("  RESULT %s: core %+.4f | template %+.4f | stretched %+.4f (k %.2f) | best %+.4f +- %.4f (stat) +- %.4f (sys) | baked %+.4f | diff %+.4f -> %s\n",
           key.c_str(), R.offCore, R.tpl.delta, R.str.delta, R.str.stretch, R.best, R.stat, R.sys, R.baked, diff, verdict.c_str());
    fprintf(csv, "%s,%.4f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%.2f,%.5f,%.5f,%.5f,%.5f,%.5f,%.5f,%s,%.5f,%.4f,%.4f,%s\n",
            key.c_str(), R.baked, R.offCore, R.errCore, R.tpl.delta, R.tpl.err, R.str.delta, R.str.err, R.str.stretch,
            R.best, R.stat, R.sys, R.sysCut, R.sysFit, R.sysMethod, R.worstCut.c_str(), diff, R.data.sigma, R.mc.sigma, verdict.c_str());

    DrawSetting(R, hdR, hm, pdf);
    all.push_back(R);
  }
  DrawSummary(all, pdf);
  cOpen.Print((pdf + "]").c_str());
  fclose(csv);
  fout.Close();
  printf("Output: %s/%s.{pdf,root,csv}\n", outDir.c_str(), outBase.c_str());
  return 0;
}
