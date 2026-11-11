/*
   This macro checks how the BigBite (BB) momentum reconstruction of the MC treats the energy
   loss of the scattered e- on its way out of the target (MC consistency audit M2, ../MC).
   E.g. run from scripts/simulation:
     analyzer -b -q $ROOTLOGON 'eloss_check_simu.cpp("configs/sbs14-sbs70p-simu/conf_qelas_ana_simu.json",2,"simc","deeN")'
   It reads the same config. blocks as qelas_ana_simu.cpp / elas_ana_simu.cpp (only the keys
   listed at the bottom of this script are used).
   -----
   Physics:
   - The MC BB optics/momentum fit (scripts/pcalib/Optics_GEN.C) was done on g4sbs e-gun events
     with the target and scattering chamber in the geometry but no BB-side shield, and was fit
     to the TRUE momentum at the vertex (ev.ep). So it should already undo, on average, the
     outgoing loss in the target + chamber (dT), but not the loss in the shield (dS).
   - The e- momentum g4sbs launches at the vertex is MC.mc_ep (for SIMC events = simc.p_e, i.e.
     after internal radiation; SIMC applies no outgoing energy loss for targ%can = 4).
   - So dp = p_recon - mc_ep should peak at:
       Scenario A (optics absorbs dT):    -dS  (shield samples)   or  0     (no-shield samples)
       Scenario B (optics absorbs nothing): -(dT + dS)            or  -dT
     dT and dS are estimated per event with the data energy-loss model (util_pd::GetElossInTgt)
     at the reconstructed vz and e- angle.
   - Histograms are NOT cross-section weighted: we study the detector response, not the yield.
   -----
   P. Datta  Created  11-11-2026
*/

#include <vector>
#include <iostream>

#include "TCut.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../include/gmn_ana.h"
#include "../../dflay/src/JSONManager.cxx"

// Histogram settings (dp in MeV; BB momentum resolution is ~1%, i.e. a few tens of MeV)
static const double kGeVtoMeV = 1E3;
static const int    kNbins_dp = 500;
static const double kdp_min = -250.;          // MeV
static const double kdp_max = 250.;           // MeV
static const int    kNbins_vz = 15;
static const double kvz_min = -0.15;          // m
static const double kvz_max = 0.15;           // m
static const int    kNbins_th = 20;
static const double kth_halfwidth = 8.;       // deg, theta range = BB angle +/- this
// Core Gaussian fit: the low-momentum side has a radiative (bremsstrahlung) tail, so we fit
// only +/- kCoreNsigma around the peak, iterating kCoreNiter times.
static const double kCoreNsigma = 1.5;
static const int    kCoreNiter = 3;
static const int    kMinEntriesSlice = 100;   // min. entries for a vz/theta slice fit

// Fits the core of a peak with a Gaussian; returns {mean, mean error, sigma}
std::vector<double> FitCore(TH1 *h) {
  double mean = h->GetBinCenter(h->GetMaximumBin());
  double sigma = h->GetRMS();
  TF1 *fg = new TF1(Form("fg_%s",h->GetName()),"gaus",kdp_min,kdp_max);
  for (int i=0; i<kCoreNiter; i++) {
    fg->SetRange(mean-kCoreNsigma*sigma, mean+kCoreNsigma*sigma);
    h->Fit(fg,"QR0");
    mean = fg->GetParameter(1);
    sigma = std::abs(fg->GetParameter(2));
  }
  return {mean, fg->GetParError(1), sigma};
}

// Median of a histogram (robust against the radiative tail)
double Median(TH1 *h) {
  double const half = 0.5;
  double median = 0.;
  h->GetQuantiles(1, &median, &half);
  return median;
}

// Fits the core of dp in every x-slice of a 2D histogram; returns a histogram of the peaks vs x
TH1D* FitCoreSlices(TH2F *h2, std::string const &name) {
  TH1D *hpeak = new TH1D(name.c_str(), "", h2->GetNbinsX(), h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax());
  for (int ix=1; ix<=h2->GetNbinsX(); ix++) {
    TH1D *hs = h2->ProjectionY(Form("%s_slice%d",name.c_str(),ix), ix, ix);
    if (hs->GetEntries() < kMinEntriesSlice) continue;
    std::vector<double> core = FitCore(hs);
    hpeak->SetBinContent(ix, core[0]);
    hpeak->SetBinError(ix, core[1]);
  }
  return hpeak;
}

int eloss_check_simu (const char *configfilename, int model=2, std::string generator="simc", std::string process="deeN")
{
  gErrorIgnoreLevel = kError; // Ignore all ROOT warnings

  // Define a clock to check macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  std::string key = generator + "_" + process + "_model" + std::to_string(model);

  // LH2 for the elastic (heep) samples, LD2 otherwise (QE deeN / deep / deen)
  std::string const target = process.compare("heep")==0 ? "LH2" : "LD2";

  int verbosefn = jmgr->GetValueFromSubKey<int>(key,"verbose_function");
  int conf = jmgr->GetValueFromSubKey<int>(key,"SBS_config");
  int sbsmag = jmgr->GetValueFromSubKey<int>(key,"SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  std::cout << sbsconf;

  // building the MC chain, exactly as in qelas_ana_simu.cpp
  std::string rfd = jmgr->GetValueFromSubKey_str(key,"rootfile_dir");
  std::vector<std::string> prefix; jmgr->GetVectorFromSubKey<std::string>(key,"prefix_to_filebase",prefix);
  int njobs = jmgr->GetValueFromSubKey<int>(key,"Njobs_to_ana");
  std::vector<SimuJob> sjobs;
  util_pd::ReadSimuJobSummary(rfd,prefix,conf,sbsmag,generator,process,njobs,verbosefn,sjobs);
  TChain *C = new TChain("T");
  util_pd::LoadSimuROOTTree(sjobs,verbosefn,C);

  // global cut (track quality, vertex, PS energy) from the config
  std::string gcut = jmgr->GetValueFromSubKey_str(key,"global_cut");
  std::cout << "\nGlobal cut: " << gcut << "\n\n";
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",(TCut)gcut.c_str(),C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr = jmgr->GetValueFromSubKey<int>(key,"max_N_tracks");
  C->SetBranchStatus("*",0);
  // reconstructed BB track
  double ntrack, p[maxNtr], px[maxNtr], py[maxNtr], pz[maxNtr], vz[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","vz"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&vz};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);
  // MC truth: e- momentum launched by g4sbs at the vertex (GeV). For SIMC events = simc.p_e.
  double mc_ep;
  std::vector<std::string> mcvar = {"mc_ep"};
  std::vector<void*> mcvar_mem = {&mc_ep};
  setrootvar::setbranch(C,"MC",mcvar,mcvar_mem);
  // turning on the remaining branches used by the global cut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);
  C->SetBranchStatus("bb.ps.e", 1);

  // defining histograms ---------------------------------------
  double const thmin = sbsconf.GetBBtheta() - kth_halfwidth;
  double const thmax = sbsconf.GetBBtheta() + kth_halfwidth;
  TH1F *h_dp = new TH1F("h_dp", "p_{recon} - p_{true} (MC);p_{recon} - p_{true} (MeV)", kNbins_dp, kdp_min, kdp_max);
  TH2F *h2_dp_vz = new TH2F("h2_dp_vz", "p_{recon} - p_{true} vs vz;vz_{recon} (m);p_{recon} - p_{true} (MeV)",
			    kNbins_vz, kvz_min, kvz_max, kNbins_dp, kdp_min, kdp_max);
  TH2F *h2_dp_th = new TH2F("h2_dp_th", "p_{recon} - p_{true} vs #theta_{e};#theta_{e,recon} (deg);p_{recon} - p_{true} (MeV)",
			    kNbins_th, thmin, thmax, kNbins_dp, kdp_min, kdp_max);
  // energy-loss model (data procedure), averaged over the same events
  TProfile *pr_dT_vz = new TProfile("pr_dT_vz", "model: -dT (target+chamber);vz_{recon} (m);MeV", kNbins_vz, kvz_min, kvz_max);
  TProfile *pr_dT_th = new TProfile("pr_dT_th", "model: -dT (target+chamber);#theta_{e,recon} (deg);MeV", kNbins_th, thmin, thmax);
  TH1F *h_dT = new TH1F("h_dT", "model: dT (target+chamber);dT (MeV)", 200, 0., 10.);
  TH1F *h_dS = new TH1F("h_dS", "model: dS (Al shield);dS (MeV)", 200, 0., 10.);

  // looping through the tree ---------------------------------------
  std::vector<double> Eloss_none, Eloss_Al; // [0]=>before, [1]=>after scattering (GeV)
  long nevent = 0, nevents = C->GetEntries(), ngoodevs = 0;
  int treenum = 0, currenttreenum = 0;
  while (C->GetEntry(nevent++)) {
    if (nevent%1000 == 0) std::cout << nevent << "/" << nevents << "\r"; std::cout.flush();

    // apply global cut (update the formula leaves whenever a new file starts)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) { treenum = currenttreenum; GlobalCut->UpdateFormulaLeaves(); }
    if (ntrack < 1 || GlobalCut->EvalInstance(0) == 0) continue;
    ngoodevs++;

    // reconstructed e- angle (as in the data scripts)
    TLorentzVector Peprime(px[0],py[0],pz[0],p[0]);
    double const etheta = kine::etheta(Peprime);

    // dp = reconstructed - true momentum at the vertex
    double const dp = (p[0] - mc_ep)*kGeVtoMeV;

    // outgoing energy loss from the data model, without and with the Al shield
    util_pd::GetElossInTgt(target, conf, vz[0], etheta, expconst::ShieldType::kNone, verbosefn, Eloss_none);
    util_pd::GetElossInTgt(target, conf, vz[0], etheta, expconst::ShieldType::kAl, verbosefn, Eloss_Al);
    double const dT = Eloss_none[1]*kGeVtoMeV;             // target + cell wall
    double const dS = (Eloss_Al[1] - Eloss_none[1])*kGeVtoMeV; // Al shield only

    double const etheta_deg = etheta*TMath::RadToDeg();
    h_dp->Fill(dp);
    h2_dp_vz->Fill(vz[0], dp);
    h2_dp_th->Fill(etheta_deg, dp);
    pr_dT_vz->Fill(vz[0], -dT);
    pr_dT_th->Fill(etheta_deg, -dT);
    h_dT->Fill(dT);
    h_dS->Fill(dS);
  }
  std::cout << "\n";

  // results ---------------------------------------
  std::vector<double> core = FitCore(h_dp);
  double const dT_mean = h_dT->GetMean(), dS_mean = h_dS->GetMean();
  std::cout << "\n------------------------------------------------------------\n"
	    << Form(" SBS%d, %s, %s: %ld events passed the cuts\n", conf, target.c_str(), process.c_str(), ngoodevs)
	    << Form(" Measured p_recon - p_true: core peak = %.2f +/- %.2f MeV (core sigma %.1f MeV)\n", core[0], core[1], core[2])
	    << Form("                            mean = %.2f MeV, median = %.2f MeV\n", h_dp->GetMean(), Median(h_dp))
	    << Form(" Model (data procedure):    dT = %.2f MeV, dS(Al) = %.2f MeV\n", dT_mean, dS_mean)
	    << " Expected peak:\n"
	    << Form("   A (optics absorbs dT):   Al-shield sample %.2f MeV | no-shield sample 0 MeV\n", -dS_mean)
	    << Form("   B (optics absorbs none): Al-shield sample %.2f MeV | no-shield sample %.2f MeV\n", -(dT_mean+dS_mean), -dT_mean)
	    << "------------------------------------------------------------\n";

  // peak vs vz and theta, compared with the model's -dT dependence
  TH1D *h_peak_vz = FitCoreSlices(h2_dp_vz, "h_peak_vz");
  TH1D *h_peak_th = FitCoreSlices(h2_dp_th, "h_peak_th");

  // output ---------------------------------------
  std::string filebase = jmgr->GetValueFromSubKey_str(key,"outfile_prefix");
  filebase = filebase.empty() ? "" : filebase + "_";
  TString outFile = Form("siout/%seloss_check_simu_%s_%s_sbs%d_sbs%dp_model%d.root",filebase.c_str(),generator.c_str(),
			 process.c_str(),conf,sbsmag,model);
  TFile *fout = new TFile(outFile.Data(),"RECREATE");
  TString outPlot = outFile; outPlot.ReplaceAll(".root",".pdf");

  TCanvas *c1 = new TCanvas("c1","dp",1200,500);
  c1->Divide(2,1);
  c1->cd(1); h_dp->Draw();
  c1->cd(2); h_dT->SetLineColor(kBlue); h_dT->Draw(); h_dS->SetLineColor(kRed); h_dS->Draw("same");
  TLegend *leg = new TLegend(0.55,0.7,0.88,0.88);
  leg->AddEntry(h_dT,"dT (target+chamber)","l"); leg->AddEntry(h_dS,"dS (Al shield)","l"); leg->Draw();
  c1->SaveAs(outPlot+"[");
  c1->SaveAs(outPlot);
  c1->Write();

  // peak vs vz / theta (black) with the model -dT shape (blue) for comparison
  TCanvas *c2 = new TCanvas("c2","dp vs vz and theta",1200,500);
  c2->Divide(2,1);
  c2->cd(1); h_peak_vz->SetTitle("core peak of p_{recon} - p_{true};vz_{recon} (m);MeV"); h_peak_vz->SetMarkerStyle(20); h_peak_vz->Draw("E");
  pr_dT_vz->SetLineColor(kBlue); pr_dT_vz->Draw("hist same");
  c2->cd(2); h_peak_th->SetTitle("core peak of p_{recon} - p_{true};#theta_{e,recon} (deg);MeV"); h_peak_th->SetMarkerStyle(20); h_peak_th->Draw("E");
  pr_dT_th->SetLineColor(kBlue); pr_dT_th->Draw("hist same");
  c2->SaveAs(outPlot);
  c2->SaveAs(outPlot+"]");
  c2->Write();

  h_dp->Write(); h2_dp_vz->Write(); h2_dp_th->Write(); h_dT->Write(); h_dS->Write();
  pr_dT_vz->Write(); pr_dT_th->Write(); h_peak_vz->Write(); h_peak_th->Write();
  fout->Close();

  std::cout << " Summary plots  : " << outPlot << "\n"
	    << " Output ROOT file  : " << outFile << "\n\n";
  sw->Stop();
  std::cout << "CPU time = " << sw->CpuTime() << "s. Real time = " << sw->RealTime() << "s.\n\n";

  delete jmgr;
  return 0;
}

/*
  Config. file parameters used (same blocks as qelas_ana_simu.cpp / elas_ana_simu.cpp):
  ** SBS_config         : SBS configuration (4, 7, 8, 9, 11, 14)
  ** SBS_magnet_percent : SBS field setting (%)
  ** rootfile_dir       : directory of the replayed MC files (with simcout/ job summaries)
  ** prefix_to_filebase : file prefixes (deeN: [deep, deen])
  ** Njobs_to_ana       : # jobs to analyze (-1 => all)
  ** global_cut         : event selection
  ** max_N_tracks       : max. # tracks (array size)
  ** outfile_prefix     : prefix for the output file names
  ** verbose_function   : verbosity of the library functions
  NOTE: the shield state of the sample is not read; the printout gives the expected peak for both.
*/
