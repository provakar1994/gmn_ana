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

int get_HALLA_p_prun (const char *configfilename) 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  
  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  std::vector<int> runnums; jmgr->GetVectorFromKey<int>("runnums",runnums);
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  //TChain *C = new TChain("E");
  TChain *C = nullptr;
  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();


  TString outFile; outFile = "test.csv";
  ofstream outFile_data; outFile_data.open(outFile);
  outFile_data << "runnum," << "evnum," << "HALLA_p" << std::endl;
  
  for (int irun=0; irun<nruns; irun++) {
    std::cout << "Analyzing run " << runnums[irun] << std::endl;

    //std::string rfname = rootfile_dir + Form("/*%d_1000k*",runnums[i]);
    std::string rfname = rootfile_dir + Form("/*%d*",runnums[irun]);
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
      outFile_data << runnums[irun] << "," << evnum << "," << HALLA_p << std::endl;
    }

    // getting ready for the next run
    C->Reset();
  }

  return 0;
}
