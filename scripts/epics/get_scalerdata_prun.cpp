/* 
   This macro will loop through the scaler tree (TSsbs) of every 
   run in a given SBS configuration and then generate 2 CSV files. 
   One file will contain:
   1)runnum, 2)segnum, 3)index, 4)sevnum, 5)gevnum 6)dnew.cnt, 7)dnew.current, 8)cum. charge, 9) daq livetime from both bbhi and ts1 scalers
   ** charge calculated with gain factor: 3317.99 +/- 31.69 Hz/uA (https://sbs.jlab.org/DocDB/0001/000164/002/dflay_bcm-ana-update_02-21-22.pdf)
   Another will contain:
   1)runnum, 2)tot. charge, 3) daq livetime from BBHi scaler, 4) daq livetime from TS1 scaler
   -----
   P. Datta  Created  03-19-2023 
*/
#include <vector>
#include <iostream>

#include "../../include/gmn_ana.h"
//#include "../../dflay/src/JSONManager.cxx"

static const int pass = 2; // replay pass
static const double dnewgain = 3317.99; // +/- 31.69 [Hz/uA], bcm gain for dnew source
static const std::string runsheet_dir = "../../DB";

int get_scalerdata_prun (const char *target,        // LH2/LD2
			 int conf,                  // SBS config
			 int sbsmagfield,           // -1=>No restriction
			 int nruns,                 // # runs to analyze
			 const char *data_dir,      // path to data directory containing rootfiles and logs
			 int verbose)               // >0=>Debug
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings
  
  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading runs list ---------------------------------------
  vector<CodaRun> crun; 
  if (sbsmagfield<0) util_pd::ReadRunList(runsheet_dir,nruns,conf,target,pass,verbose,crun);
  else util_pd::ReadRunList(runsheet_dir,nruns,conf,target,pass,sbsmagfield,verbose,crun);
  TChain *C = nullptr;

  TString outFile1, outFile2;
  if (verbose==0) {
    if (sbsmagfield<0) {
      outFile1 = Form("epout/scalerdata_prun_SBS%d_%s.root",conf,target);
      outFile2 = Form("epout/beamCh_n_daqLvTm_prun_SBS%d_%s.csv",conf,target);
    } else {
      outFile1 = Form("epout/scalerdata_prun_SBS%d_%dp_%s.root",conf,sbsmagfield,target);
      outFile2 = Form("epout/beamCh_n_daqLvTm_prun_SBS%d_%dp_%s.csv",conf,sbsmagfield,target);
    }  
  } else {
    if (sbsmagfield<0) {
      outFile1 = Form("epout/test_scalerdata_prun_SBS%d_%s.root",conf,target);
      outFile2 = Form("epout/test_beamCh_n_daqLvTm_prun_SBS%d_%s.csv",conf,target);
    } else {
      outFile1 = Form("epout/test_scalerdata_prun_SBS%d_%dp_%s.root",conf,sbsmagfield,target);
      outFile2 = Form("epout/test_beamCh_n_daqLvTm_prun_SBS%d_%dp_%s.csv",conf,sbsmagfield,target);
    }  
  }
  TFile *fout = new TFile(outFile1.Data(),"RECREATE");
  ofstream outFile_data2; outFile_data2.open(outFile2);
  outFile_data2 << "runnum," << "totcharge_pd(C)," << "daqlvtm_bbhi," << "daqlvtm_ts1" << std::endl;

  // defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  // ev info
  UInt_t T_rnum;        Tout->Branch("rnum", &T_rnum, "rnum/i");
  UInt_t T_segnum;      Tout->Branch("segnum", &T_segnum, "segnum/i");
  ULong64_t T_evindex;  Tout->Branch("evindex", &T_evindex, "evindex/l");
  ULong64_t T_sevnum;   Tout->Branch("sevnum", &T_sevnum, "sevnum/l");
  ULong64_t T_gevnum;   Tout->Branch("gevnum", &T_gevnum, "gevnum/l");
  // bcm info
  double T_dnewcnt;     Tout->Branch("dnewcnt", &T_dnewcnt, "dnewcnt/D");
  double T_dnewcurr;    Tout->Branch("dnewcurr", &T_dnewcurr, "dnewcurr/D");
  double T_dnewcharge;  Tout->Branch("dnewcharge", &T_dnewcharge, "dnewcharge/D");
  // trigger info
  double T_l1ascaler;   Tout->Branch("l1ascaler", &T_l1ascaler, "l1ascaler/D");
  double T_bbhiscaler;   Tout->Branch("bbhiscaler", &T_bbhiscaler, "bbhiscaler/D");
  double T_ts1scaler;   Tout->Branch("ts1scaler", &T_ts1scaler, "ts1scaler/D");

  // looping through runs
  for (int irun=0; irun<nruns; irun++) {
    int runnum = crun[irun].runnum;
    std::cout << "Analyzing run " << crun[irun].runnum << std::endl;
    C = new TChain("TSsbs"); int lr = util_pd::LoadROOTTree(data_dir,crun[irun],1,verbose,1,C);
    if (lr!=0) continue;

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

    // trigger variables
    double l1ascaler; C->SetBranchAddress("sbs.L1A.scaler",&l1ascaler);
    double bbhiscaler; C->SetBranchAddress("sbs.BBCalHi.BBCALTRG.scaler",&bbhiscaler);
    double ts1scaler; C->SetBranchAddress("sbs.TS1_BB.scaler",&ts1scaler);

    // looping through the events ---------------------------------------
    std::cout << std::endl;
    double dnewcharge_cum=0., l1ascaler_cum=0., bbhiscaler_cum=0., ts1scaler_cum=0.;
    long nevent=0, nevents=C->GetEntries(); 
    int treenum=0, currenttreenum=0, segnum=0, index=0;
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

      if (dnewcnt>0) //avoid last segment with zero count 
	dnewcharge_cum = (dnewcnt / dnewgain) * 1e-6; //C 

      if (l1ascaler>0) //avoid last segment with zero count 
	l1ascaler_cum = l1ascaler;

      if (bbhiscaler>0) //avoid last segment with zero count 
	bbhiscaler_cum = bbhiscaler;
      
      if (ts1scaler>0) //avoid last segment with zero count 
	ts1scaler_cum = ts1scaler;

      T_rnum = runnum;
      T_segnum = segnum;
      T_evindex = index;
      T_sevnum = sevnum;
      T_gevnum = gevnum;

      T_dnewcnt = dnewcnt;
      T_dnewcurr = dnewcurr;
      T_dnewcharge = dnewcharge_cum;

      T_l1ascaler = l1ascaler_cum;
      T_bbhiscaler = bbhiscaler_cum;
      T_ts1scaler = ts1scaler_cum;

      index++;  // gets reset at the beginning of every run
      if (verbose>0 && nevent<5)
	std::cout << runnum << "," << segnum << "," << index << "," << sevnum << "," << gevnum << "," 
		  << dnewcnt << "," << dnewcurr << "," << dnewcharge_cum << "," << std::endl;

      Tout->Fill();
    } //event loop

    // DAQ livetime: Accepted BBCal singles trigger
    double daqlvtm_bbhi = (double)crun[irun].BBCalSinglesPassed / (double)bbhiscaler_cum;
    double daqlvtm_ts1 = (double)crun[irun].BBCalSinglesPassed / (double)ts1scaler_cum;
    
    // NOTE: Total charge is nothing but the last cumulative charge entry
    outFile_data2 << runnum << "," << dnewcharge_cum << "," << daqlvtm_bbhi << "," << daqlvtm_ts1 << std::endl;
    if (verbose>0) std::cout << runnum << "," << dnewcharge_cum << "," << daqlvtm_bbhi << "," << daqlvtm_ts1 << std::endl;
    
    // getting ready for next run
    C->Reset();
  } //run loop
  std::cout << std::endl << std::endl;

  cout << "------" << endl;
  cout << " Output CSV file  : " << outFile2 << endl;
  cout << " Output ROOT file : " << outFile1 << endl;
  cout << "------" << endl << endl;

  sw->Stop();
  cout << "CPU time elapsed = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s. " << endl << endl;

  fout->Write();
  sw->Delete();
  // delete jmgr;
  return 0;
}
