/* 
   This macro will loop through the scaler tree (TSsbs) of every 
   run in a given SBS configuration and then generate 2 CSV files. 
   One file will contain:
   1)runnum, 2)segnum, 3)sevnum, 4)gevnum 5)dnew.cnt, 6)dnew.current 7)cum. charge
   ** charge calculated with gain factor: 3317.99 +/- 31.69 Hz/uA
   Another will contain:
   1)runnum, 2)tot. charge
   -----
   P. Datta  Created  03-19-2023 
*/
#include <vector>
#include <iostream>

#include "../../include/gmn-ana.h"
#include "../../dflay/src/JSONManager.cxx"

// bcm gain factors from calibration
static const double dnewgain = 3317.99; // +/- 31.69 [Hz/uA]

int get_scalerdata_prun (const char *configfilename) 
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  
  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int pass = jmgr->GetValueFromKey<int>("replay_pass");
  std::string target = jmgr->GetValueFromKey_str("target");

  // reading run info and parsing ROOT trees
  std::string runsheet_dir = jmgr->GetValueFromKey_str("runsheet_dir");
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # of runs to analyze
  vector<CodaRun> crun; util_pd::ReadRunList(runsheet_dir,nruns,conf,target,pass,0,crun);
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  TChain *C = nullptr;

  TString outFile1; outFile1 = Form("epout/get_scalerdata_prun_SBS%d_%s.csv",conf,target.c_str());
  TString outFile2; outFile2 = Form("epout/get_beamcharge_prun_SBS%d_%s.csv",conf,target.c_str());
  ofstream outFile_data1; outFile_data1.open(outFile1);
  ofstream outFile_data2; outFile_data2.open(outFile2);
  outFile_data1 << "runnum," << "segnum," << "sevnum," << "gevnum," << "dnewcnt," << "dnewcurr," << "cumcharge(C)" << std::endl;
  outFile_data2 << "runnum," << "totcharge(C)" << std::endl;

  // looping through runs
  for (int irun=0; irun<nruns; irun++) {
    std::cout << "Analyzing run " << crun[irun].runnum << std::endl;
    C = new TChain("TSsbs"); util_pd::LoadROOTTree(rootfile_dir,crun[irun],1,0,C);

    // setting up ROOT tree branch addresses ----------------------------
    C->SetBranchStatus("*",0);
    // enevt variables
    double sevnum; C->SetBranchAddress("evcount",&sevnum);
    long long gevnum;   C->SetBranchAddress("evnum",&gevnum);

    // dnew variables
    double dnewcnt, dnewcurr;
    std::vector<std::string> dnewvar = {"cnt","current"};
    std::vector<void*> dnewvar_mem = {&dnewcnt,&dnewcurr};
    setrootvar::setbranch(C, "sbs.bcm.dnew", dnewvar, dnewvar_mem);

    // looping through the events ---------------------------------------
    std::cout << std::endl;
    double dnewcharge_cum = 0.;
    long nevent = 0, nevents = C->GetEntries(); 
    int treenum = 0, currenttreenum = 0; int segnum = 0;
    while (C->GetEntry(nevent++)) {
   
      // print progress 
      if( nevent % 1000 == 0 ) std::cout << nevent << "/" << nevents << "\r";
      std::cout.flush();

      // keep track of tree number
      currenttreenum = C->GetTreeNumber();
      if (nevent == 1 || currenttreenum != treenum) {
      	treenum = currenttreenum;
	segnum++;  //keep track of segment number
      } 

      dnewcharge_cum = (dnewcnt / dnewgain) * 1e-6; //C 

      outFile_data1 << crun[irun].runnum << "," << segnum << "," << sevnum << "," << gevnum << "," 
		    << dnewcnt << "," << dnewcurr << "," << dnewcharge_cum << "," << std::endl;
    } //event loop
    
    outFile_data2 << crun[irun].runnum << "," << dnewcharge_cum << std::endl;
    
    // getting ready for next run
    C->Reset();
  } //run loop
  std::cout << std::endl;


  sw->Stop();
  cout << "CPU time elapsed = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s. " << endl << endl;

  sw->Delete();
  delete jmgr;
  return 0;
}
