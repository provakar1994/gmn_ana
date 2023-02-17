/* 
   This macro will loop through the EPICS tree of every run in a
   given SBS configuration and then generate a CSV file with the
   following three values: 1. runnum, 2. evnum, 3. HALLA_p
   -----
   P. Datta  Created  11-02-2022 
*/
#include <vector>
#include <iostream>

#include "../../include/gmn-ana.h"
#include "../../dflay/src/JSONManager.cxx"

int get_avgHALLA_p_prun (const char *configfilename) 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  
  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  std::string target = jmgr->GetValueFromKey_str("target");
  int pass = jmgr->GetValueFromKey<int>("replay_pass"); 

  // reading relevant run info from relevant good runlist spreadsheet
  std::string runsheet_dir = jmgr->GetValueFromKey_str("runsheet_dir");
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  vector<CodaRun> crun; util_pd::ReadRunList(runsheet_dir,nruns,conf,target,pass,crun);
 
  // parsing trees depending on target type
  TChain *C = nullptr;
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");

  TString outFile; outFile = Form("epout/get_HALLA_p_prun_SBS%d_%s.csv",conf,target.c_str());
  ofstream outFile_data; outFile_data.open(outFile);
  outFile_data << "runnum," << "evnum," << "HALLA_p" << std::endl;
  
  for (int irun=0; irun<nruns; irun++) {
    std::cout << "Analyzing run " << crun[irun].runnum << std::endl;

    std::string rfname = rootfile_dir + Form("/*%d*",crun[irun].runnum);
    C = new TChain("E");
    C->Add(rfname.c_str());

    // setting up ROOT tree branch addresses ---------------------------------------
    C->SetBranchStatus("*",0);
    // beam energy
    double HALLA_p; setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);
    // global enent number
    long evnum;   setrootvar::setbranch(C, "evnum", "", &evnum);
    
    // looping through EPICS events
    long nevent = 0, nevents = C->GetEntries();  
    while (C->GetEntry(nevent++)) {
      outFile_data << crun[irun].runnum << "," << evnum << "," << HALLA_p << std::endl;
    }

    // getting ready for the next run
    C->Reset();
  }

  return 0;
}
