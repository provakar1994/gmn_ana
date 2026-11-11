/*
  Unit test: util_pd::LoadSimuROOTTree (and its helper util_pd::GetTreeEntries).
  ----
  LoadSimuROOTTree must drop every MC job whose replayed ROOT file cannot be used,
  because the MC normalization (total # tries, charge) is summed over the jobs it
  keeps. A dropped job must not count; a kept job must be in the TChain.
  The test writes small ROOT files to a temporary directory:
    good      : tree T with entries             -> kept
    low       : tree T with few entries         -> kept (only a warning)
    empty     : tree T with 0 entries           -> dropped
    notree    : valid ROOT file, no tree T      -> dropped
    garbage   : not a ROOT file (corrupted)     -> dropped
    missing   : file does not exist             -> dropped
  Checks: the kept job list, the TChain entry count, GetTreeEntries for each case,
  and that the total # tries summed after loading only counts the kept jobs.
*/
#include "Utilities.h"
#include "TestUtils.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TError.h"
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

namespace {

  long const kNGenReq     = 100;  // events requested per job
  long const kNGood       = 90;   // entries in a healthy file (some events have no hits)
  long const kNLow        = 10;   // entries in a suspiciously short file (< 50% of kNGenReq)
  double const kNTriedJob = 1000; // # tries per job, distinct per job below to identify it

  std::string const kTreeName = "T";

  // Writes a ROOT file with tree T holding <nentries> entries
  void WriteTree(std::string const& fname, long nentries) {
    TFile f(fname.c_str(), "RECREATE");
    TTree T(kTreeName.c_str(), "");
    double x = 0;
    T.Branch("x", &x, "x/D");
    for (long i = 0; i < nentries; i++) { x = i; T.Fill(); }
    T.Write();
  }

  // Makes a job pointing to <fname>; ntried tags the job so the sum identifies kept jobs
  SimuJob MakeJob(std::string const& fname, double ntried) {
    SimuJob job;
    job.rfname = fname;
    job.ngenreq = kNGenReq;
    job.ntried = ntried;
    job.is_rfile = true;
    return job;
  }

} //::

int main() {
  gErrorIgnoreLevel = kFatal;  // the bad files are expected to make ROOT complain

  // Write the test files to a fresh temporary directory
  std::filesystem::path const dir = std::filesystem::temp_directory_path() / "gmnana_test_load_simu_tree";
  std::filesystem::remove_all(dir);
  std::filesystem::create_directories(dir);
  std::string const fgood    = (dir / "good.root").string();
  std::string const flow     = (dir / "low.root").string();
  std::string const fempty   = (dir / "empty.root").string();
  std::string const fnotree  = (dir / "notree.root").string();
  std::string const fgarbage = (dir / "garbage.root").string();
  std::string const fmissing = (dir / "missing.root").string();

  WriteTree(fgood, kNGood);
  WriteTree(flow, kNLow);
  WriteTree(fempty, 0);
  { TFile f(fnotree.c_str(), "RECREATE"); TH1F h("h", "", 10, 0, 1); h.Write(); }
  { std::ofstream out(fgarbage); out << "this is not a ROOT file\n"; }

  // GetTreeEntries for each case
  CHECK_TRUE(util_pd::GetTreeEntries(fgood, kTreeName) == kNGood, "GetTreeEntries: good file");
  CHECK_TRUE(util_pd::GetTreeEntries(fempty, kTreeName) == 0, "GetTreeEntries: empty tree -> 0");
  CHECK_TRUE(util_pd::GetTreeEntries(fnotree, kTreeName) == -1, "GetTreeEntries: no tree -> -1");
  CHECK_TRUE(util_pd::GetTreeEntries(fgarbage, kTreeName) == -1, "GetTreeEntries: corrupted -> -1");
  CHECK_TRUE(util_pd::GetTreeEntries(fmissing, kTreeName) == -1, "GetTreeEntries: missing -> -1");

  // Jobs in summary order; each has a different ntried (1x, 2x, ... kNTriedJob)
  std::vector<SimuJob> sjobs = {
    MakeJob(fgood,    1*kNTriedJob),
    MakeJob(fempty,   2*kNTriedJob),
    MakeJob(fnotree,  3*kNTriedJob),
    MakeJob(flow,     4*kNTriedJob),
    MakeJob(fgarbage, 5*kNTriedJob),
    MakeJob(fmissing, 6*kNTriedJob),
  };

  TChain* C = new TChain(kTreeName.c_str());
  util_pd::LoadSimuROOTTree(sjobs, 0, C);

  // Only the good and the low jobs are kept, in their original order
  CHECK_TRUE(sjobs.size() == 2, "two jobs kept");
  if (sjobs.size() == 2) {
    CHECK_TRUE(sjobs[0].rfname == fgood, "first kept job is the good file");
    CHECK_TRUE(sjobs[1].rfname == flow, "second kept job is the low file (warned, not dropped)");
  }
  CHECK_TRUE(C->GetNtrees() == 2, "TChain holds two files");
  CHECK_TRUE(C->GetEntries() == kNGood + kNLow, "TChain entries = good + low");
  // Entry counts are passed to TChain::Add, so check the last entry is really readable
  CHECK_TRUE(C->GetEntry(kNGood + kNLow - 1) > 0, "last chain entry can be read");

  // Normalization summed after loading counts only the kept jobs (1x + 4x)
  std::vector<long double> totNtriesnCh;
  util_pd::GetTotNtriesnCh(sjobs, totNtriesnCh);
  CHECK_REL((double)totNtriesnCh[0], 5*kNTriedJob, 1e-12, "total # tries = kept jobs only");

  // All files bad -> must throw rather than return an empty job list
  std::vector<SimuJob> badjobs = { MakeJob(fempty, kNTriedJob), MakeJob(fmissing, kNTriedJob) };
  TChain* Cbad = new TChain(kTreeName.c_str());
  bool threw = false;
  try { util_pd::LoadSimuROOTTree(badjobs, 0, Cbad); } catch (std::runtime_error const&) { threw = true; }
  CHECK_TRUE(threw, "throws when no file is usable");

  delete C; delete Cbad;
  std::filesystem::remove_all(dir);
  return test::Summary("test_load_simu_tree");
}
