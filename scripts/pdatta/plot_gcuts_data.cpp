/*
  This macro plots distributions which are used as global cuts.
  -----
  P. Datta Created 01-18-2023
*/
#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TLatex.h"
#include "TChain.h"
#include "TVector3.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../../include/gmn-ana.h"
#include "../../dflay/src/JSONManager.cxx"

int plot_gcuts_data (const char *configfilename, std::string filebase="pdout/test_plot_gcuts_data")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string target = jmgr->GetValueFromKey_str("target");
  std::string rootfile_dir; //= jmgr->GetValueFromKey_str("rootfile_dir");
  std::vector<int> runnums; //jmgr->GetVectorFromKey<int>("runnums",runnums);
  if (target.compare("lh2") == 0) {
    rootfile_dir =  jmgr->GetValueFromKey_str("rootfile_dir_lh2");
    jmgr->GetVectorFromKey<int>("runnums_lh2",runnums);
  } else if (target.compare("ld2") == 0) {
    rootfile_dir =  jmgr->GetValueFromKey_str("rootfile_dir_ld2");
    jmgr->GetVectorFromKey<int>("runnums_ld2",runnums);
  } else {
    std::cerr << "[Parsing error] Enter valid target type!" << std::endl; throw;
  }
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # runs to analyze
  TChain *C = new TChain("T");
  if (nruns < 1 || nruns > runnums.size()) nruns = runnums.size();
  for (int i=0; i<nruns; i++) {
    //std::string rfname = rootfile_dir + Form("/*%d_1000k*",runnums[i]);
    std::string rfname = rootfile_dir + Form("/*%d*",runnums[i]);
    C->Add(rfname.c_str());
  }
  if (C->GetEntries()==0) {std::cerr << "*!* No ROOT file!" << std::endl; throw;}

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  double sbsfieldfrac = sbsmag / 100.;
  SBSconfig sbsconf(conf, sbsmag);
  sbsconf.Print();

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy - Probably we should take an average over 100 events
  // double HALLA_p;
  // setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  // bbcal clus var
  double eSH;       setrootvar::setbranch(C,"bb.sh","e",&eSH);
  double ePS;       setrootvar::setbranch(C,"bb.ps","e",&ePS);
  
  // hcal clus var
  double eHCAL;     setrootvar::setbranch(C,"sbs.hcal","e",&eHCAL);

  // track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&xTr,&yTr,&thTr,&phTr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  // tdctrig variable (N/A for simulation)
  int tdcElemN;
  double tdcTrig[maxNtr], tdcElem[maxNtr];
  std::vector<std::string> tdcvar = {"tdcelemID","tdcelemID","tdc"};
  std::vector<void*> tdcvar_mem = {&tdcElem,&tdcElemN,&tdcTrig};
  setrootvar::setbranch(C,"bb.tdctrig",tdcvar,tdcvar_mem,1);

  // Other branches
  double e_ov_p;    setrootvar::setbranch(C,"bb","etot_over_p",&e_ov_p);
  double W2;        setrootvar::setbranch(C,"e.kine","W2",&W2);

  // defining the outputfile
  int pass = jmgr->GetValueFromKey<int>("replay_pass"); 
  TString outFile = Form("%s_sbs%d_sbs%dp_pass%d_%s.root", 
			 filebase.c_str(), sbsconf.GetSBSconf(), sbsconf.GetSBSmag(), pass, target.c_str());
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining tree branches
  TTree *Tout = new TTree("Tout", "");
  //kine
  double T_W;           Tout->Branch("W", &T_W, "W/D");
  double T_W2;          Tout->Branch("W2", &T_W2, "W2/D");
  //track
  double T_vz;          Tout->Branch("vz", &T_vz, "vz/D");
  double T_trP;         Tout->Branch("trP", &T_trP, "trP/D");
  double T_trX;         Tout->Branch("trX", &T_trX, "trX/D");
  double T_trY;         Tout->Branch("trY", &T_trY, "trY/D");
  double T_trTh;        Tout->Branch("trTh", &T_trTh, "trTh/D");
  double T_trPh;        Tout->Branch("trPh", &T_trPh, "trPh/D");
  //BBCAL
  double T_ePS;         Tout->Branch("ePS", &T_ePS, "ePS/D"); 
  double T_eSH;         Tout->Branch("eSH", &T_eSH, "eSH/D"); 
  double T_EovP;        Tout->Branch("EovP", &T_EovP, "EovP/D"); 
  //HCAL
  double T_eHCAL;       Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  //coin time trigger
  double T_coinT_trig;  Tout->Branch("coinT_trig", &T_coinT_trig, "coinT_trig/D");

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  long nevent = 0, nevents = C->GetEntries(); 
  while (C->GetEntry(nevent++)) {
    
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    // coin time cut (N/A for simulation) !! Not a reliable cut - loosing a lot of elastics
    double bbcal_time=0., hcal_time=0.;
    for(int ihit=0; ihit<tdcElemN; ihit++){
      if(tdcElem[ihit]==5) bbcal_time=tdcTrig[ihit];
      if(tdcElem[ihit]==0) hcal_time=tdcTrig[ihit];
    }
    double coin_time = hcal_time - bbcal_time;  
    T_coinT_trig = coin_time;

    T_W = max(0., sqrt(W2));
    T_W2 = W2;
    
    T_vz = vz[0];
    T_trP = p[0];
    T_trX = xTr[0];
    T_trY = yTr[0];
    T_trTh = thTr[0];
    T_trPh = phTr[0];

    T_ePS = ePS;
    T_eSH = eSH;
    T_EovP = e_ov_p; 

    T_eHCAL = eHCAL;
    
    Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;

  cout << "------" << endl;
  cout << " Output file : " << outFile << endl;
  cout << "------" << endl << endl;

  sw->Stop();
  cout << "CPU time elapsed = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s. " << endl << endl;

  fout->Write();
  sw->Delete();
  delete jmgr;
  return 0;
}
