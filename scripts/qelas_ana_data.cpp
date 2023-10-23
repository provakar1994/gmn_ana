/* 
   This macro will perform QE analysis for GMn using LD2 data.
   E.g. config. file: sbs14-sbs70p/conf_qelas_ana_data.json
   * A brief description of all the config. file parameters can be
   found at the bottom of this script.
   -----
   P. Datta  Created  11-02-2022 
*/

// TO-DO
// 1. Energy loss calculations - Done (Cell diameter and thickness are guesses)

#include <vector>
#include <iostream>

#include "TCut.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLatex.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TVector3.h"
#include "TEllipse.h"
#include "TStopwatch.h"
#include "TTreeFormula.h"
#include "TLorentzVector.h"

#include "../include/gmn_ana.h"
#include "../dflay/src/JSONManager.cxx"

/* this script will only analyze LD2 data */
static const std::string target = "LD2";

int qelas_ana_data (const char *configfilename,
		    int pass, //replay pass
                    int verbose=-1,  //<-1=>Debug, =-1=>Test
                    int verbosefn=0, //>0=>Debug
		    /* verbose==verbosefn==0 => Production */
                    std::string filebase="pdout/test_qelas_ana")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  std::string key = "pass" + std::to_string(pass);

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromSubKey<int>(key,"SBS_config");
  int sbsmag = jmgr->GetValueFromSubKey<int>(key,"SBS_magnet_percent"); 
  SBSconfig sbsconf(conf, sbsmag);
  std::cout << sbsconf;

  // reading run info and parsing ROOT trees
  std::string runsheet_dir = jmgr->GetValueFromSubKey_str(key,"runsheet_dir");
  int nruns = jmgr->GetValueFromSubKey<int>(key,"Nruns_to_ana"); // # of runs to analyze
  vector<CodaRun> crun; util_pd::ReadRunList(runsheet_dir,nruns,conf,target,pass,sbsmag,verbosefn,crun);
  std::string rootfile_dir = jmgr->GetValueFromSubKey_str(key,"rootfile_dir");
  TChain *C = new TChain("T"); util_pd::LoadROOTTree(rootfile_dir,crun,1,verbosefn,C); 
 
  // reading scaler tree
  TChain *S = new TChain("Tout"); 
  int get_scaler_info = jmgr->GetValueFromSubKey<int>(key,"get_scaler_info");
  if (get_scaler_info) {
    S->Add(Form("epics/epout/scalerdata_prun_SBS%d_%s.root",conf,target.c_str()));
    if (S->GetEntries()==0) throw std::runtime_error("No scaler event found!");
  }

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = jmgr->GetValueFromSubKey<int>(key,"model");
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromSubKey_str(key,"Ntype");

  // setting up global cuts
  std::string gcut = jmgr->GetValueFromSubKey_str(key,"global_cut");
  std::vector<std::string> gCutList; util_pd::SplitString('&',gcut,gCutList);
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",(TCut)gcut.c_str(),C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy 
  // double HALLA_p; setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  // bbsh clus var
  double eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH;
  std::vector<std::string> bbshclvar = {"e","x","y","rowblk","colblk","idblk","atimeblk"};
  std::vector<void*> bbshclvar_mem = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH};
  setrootvar::setbranch(C, "bb.sh", bbshclvar, bbshclvar_mem);

  // bbps clus var
  double ePS, rblkPS, cblkPS, idblkPS, atimePS;
  std::vector<std::string> bbpsclvar = {"e","rowblk","colblk","idblk","atimeblk"};
  std::vector<void*> bbpsclvar_mem = {&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
  setrootvar::setbranch(C, "bb.ps", bbpsclvar, bbpsclvar_mem);
 
  // hcal clus var
  double eHCAL, xHCAL, yHCAL, rblkHCAL, cblkHCAL, idblkHCAL, atimeHCAL, tdcHCAL[maxNtr];
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","atimeblk","clus_blk.tdctime"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&atimeHCAL,&tdcHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // // hcal clus var [2]
  // int maxHCALcl = 50;
  // int idHCAL_clN; double idHCAL_cl[maxHCALcl];
  // double nblkHCAL_cl[maxHCALcl], eblkHCAL_cl[maxHCALcl], atimeblkHCAL_cl[maxHCALcl], tdcblkHCAL_cl[maxHCALcl]; 
  // double eHCAL_cl[maxHCALcl], HCAL_cl[maxHCALcl], yHCAL_cl[maxHCALcl], rblkHCAL_cl[maxHCALcl], cblkHCAL_cl[maxHCALcl];
  // std::vector<std::string> hcalclvar_cl = {"id","id","nblk","eblk","e","x","y","rowblk","colblk","atimeblk","tdctime"};
  // std::vector<void*> hcalclvar_cl_mem = {&idHCAL_cl,&idHCAL_clN,&nblkHCAL_cl,&eblkHCAL_cl,&eHCAL_cl,&xHCAL_cl,&yHCAL_cl,&rblkHCAL_cl,&cblkHCAL_cl,&atimeblkHCAL_cl,&tdcblkHCAL_cl};
  // setrootvar::setbranch(C, "sbs.hcal.clus", hcalclvar_cl, hcalclvar_cl_mem, 1);

  // bbhodo clus var
  int ncltmeanHODO; 
  double cltmeanHODO[maxNtr];
  std::vector<std::string> hodoclvar = {"clus.tmean","clus.tmean"};
  std::vector<void*> hodoclvar_mem = {&ncltmeanHODO,&cltmeanHODO};
  setrootvar::setbranch(C, "bb.hodotdc", hodoclvar, hodoclvar_mem, 0);  

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

  // fEvtHdr variables (N/A for simulation) 
  UInt_t rnum, gevnum, trigbits;
  std::vector<std::string> evhdrvar = {"fRun","fEvtNum","fTrigBits"};
  std::vector<void*> evhdrvar_mem = {&rnum,&gevnum,&trigbits};
  setrootvar::setbranch(C,"fEvtHdr",evhdrvar,evhdrvar_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.etot_over_p", 1);
  C->SetBranchStatus("bb.gem.track.nhits", 1);

  // scaler tree variables
  UInt_t rnumS=0, segnumS;
  double dnewcnt, dnewcurr;
  ULong64_t evindex, gevnumS;
  std::vector<std::string> streevar = {"rnum","segnum","gevnum","evindex","dnewcnt","dnewcurr"};
  std::vector<void*> streevar_mem = {&rnumS,&segnumS,&gevnumS,&evindex,&dnewcnt,&dnewcurr};
  if (get_scaler_info) setrootvar::setbranch(S,"",streevar,streevar_mem);

  // defining the outputfile
  if (verbose==0 && verbosefn==0) filebase = "pdout/qelas_ana";
  TString outFile = Form("%s_data_sbs%d_sbs%dp_model%d_pass%d.root",filebase.c_str(),conf,sbsmag,model,pass);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms
  TH1F *h_W = util_pd::TH1FhW("h_W");
  TH1F *h_W_cut = util_pd::TH1FhW("h_W_cut");
  TH1F *h_W_acut = util_pd::TH1FhW("h_W_acut");
  TH1F *h_dpel = new TH1F("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = util_pd::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_lims",hdx_lim);
  vector<double> hdy_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dyHCAL_lims",hdy_lim);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","W & fiducial cuts;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dxHCAL_nfc = new TH1F("h_dxHCAL_nfc","W cut;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","W & fiducial cuts;y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_dyHCAL_nfc = new TH1F("h_dyHCAL_nfc","W cut;y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_coin_time = new TH1F("h_coin_time","BBCAL-HCAL trigger coincidence time (ns)",200,380,660);

  TH2F *h2_rcHCAL = util_pd::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = util_pd::TH2FdxdyHCAL("h2_dxdyHCAL");

  // defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  Tout->SetMaxTreeSize(4000000000LL);
  //cuts
  bool WCut;            Tout->Branch("WCut", &WCut, "WCut/O");
  bool pCut;            Tout->Branch("pCut", &pCut, "pCut/O");
  bool nCut;            Tout->Branch("nCut", &nCut, "nCut/O");
  // -- a few variations
  bool pCut_1p5sig;     Tout->Branch("pCut_1p5sig", &pCut_1p5sig, "pCut_1p5sig/O");
  bool nCut_1p5sig;     Tout->Branch("nCut_1p5sig", &nCut_1p5sig, "nCut_1p5sig/O");
  bool pCut_2sig;       Tout->Branch("pCut_2sig", &pCut_2sig, "pCut_2sig/O");
  bool nCut_2sig;       Tout->Branch("nCut_2sig", &nCut_2sig, "nCut_2sig/O");
  bool pCut_2p5sig;     Tout->Branch("pCut_2p5sig", &pCut_2p5sig, "pCut_2p5sig/O");
  bool nCut_2p5sig;     Tout->Branch("nCut_2p5sig", &nCut_2p5sig, "nCut_2p5sig/O");
  bool pCut_3sig;       Tout->Branch("pCut_3sig", &pCut_3sig, "pCut_3sig/O");
  bool nCut_3sig;       Tout->Branch("nCut_3sig", &nCut_3sig, "nCut_3sig/O");
  // --
  bool SMCut;           Tout->Branch("SMCut", &SMCut, "SMCut/B");
  bool ARCut;           Tout->Branch("ARCut", &ARCut, "ARCut/B");
  bool fiduCut;         Tout->Branch("fiduCut", &fiduCut, "fiduCut/O");
  //run info
  UInt_t T_rnum;        Tout->Branch("rnum", &T_rnum, "rnum/i");
  UInt_t T_segnum;      Tout->Branch("segnum", &T_segnum, "segnum/i");
  ULong64_t T_gevnum;   Tout->Branch("gevnum", &T_gevnum, "gevnum/l");
  double T_ebeam;       Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  double T_ebeam_corr;  Tout->Branch("ebeam_corr", &T_ebeam_corr, "ebeam_corr/D");
  double T_ebeam_std;   Tout->Branch("ebeam_std", &T_ebeam_std, "ebeam_std/D");
  //bcm/scaler
  //UInt_t T_segnumS;   if (get_scaler_info) Tout->Branch("segnumS", &T_segnumS, "segnumS/i");
  double T_dnewcnt;     if (get_scaler_info) Tout->Branch("dnewcnt", &T_dnewcnt, "dnewcnt/D"); 
  double T_dnewcurr;    if (get_scaler_info) Tout->Branch("dnewcurr", &T_dnewcurr, "dnewcurr/D"); 
  //kine
  double T_nu;          Tout->Branch("nu", &T_nu, "nu/D");
  double T_Q2;          Tout->Branch("Q2", &T_Q2, "Q2/D");
  double T_W2;          Tout->Branch("W2", &T_W2, "W2/D");
  double T_W;           Tout->Branch("W", &T_W, "W/D");
  double T_dpel;        Tout->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi;        Tout->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta;      Tout->Branch("etheta", &T_etheta, "etheta/D");
  double T_pelas;       Tout->Branch("pelas", &T_pelas, "pelas/D");
  double T_epsilon;     Tout->Branch("epsilon", &T_epsilon, "epsilon/D"); // calculated using general eqn.
  double T_epsilon_p;   Tout->Branch("epsilon_p", &T_epsilon_p, "epsilon_p/D");
  double T_epsilon_n;   Tout->Branch("epsilon_n", &T_epsilon_n, "epsilon_n/D");
  double T_thpq_p;      Tout->Branch("thpq_p", &T_thpq_p, "thpq_p/D");
  double T_thpq_n;      Tout->Branch("thpq_n", &T_thpq_n, "thpq_n/D");
  //track
  double T_vz;          Tout->Branch("vz", &T_vz, "vz/D");
  double T_trP;         Tout->Branch("trP", &T_trP, "trP/D");
  double T_trP_corr;    Tout->Branch("trP_corr", &T_trP_corr, "trP_corr/D");
  double T_trX;         Tout->Branch("trX", &T_trX, "trX/D");
  double T_trY;         Tout->Branch("trY", &T_trY, "trY/D");
  double T_trTh;        Tout->Branch("trTh", &T_trTh, "trTh/D");
  double T_trPh;        Tout->Branch("trPh", &T_trPh, "trPh/D");
  double T_tgX;         Tout->Branch("tgX", &T_tgX, "tgX/D");
  double T_tgY;         Tout->Branch("tgY", &T_tgY, "tgY/D");
  double T_tgTh;        Tout->Branch("tgTh", &T_tgTh, "tgTh/D");
  double T_tgPh;        Tout->Branch("tgPh", &T_tgPh, "tgPh/D");
  //cross-section predicted from MC
  double T_sigMott;     Tout->Branch("sigMott", &T_sigMott, "sigMott/D");
  double T_sigRed_p;    Tout->Branch("sigRed_p", &T_sigRed_p, "sigRed_p/D");
  double T_sigRed_n;    Tout->Branch("sigRed_n", &T_sigRed_n, "sigRed_n/D");
  double T_sigBorn_ratio; Tout->Branch("sigBorn_ratio", &T_sigBorn_ratio, "sigBorn_ratio/D");
  //BBCAL
  double T_ePS;         Tout->Branch("ePS", &T_ePS, "ePS/D"); 
  double T_rblkPS;      Tout->Branch("rblkPS", &T_rblkPS, "rblkPS/D"); 
  double T_cblkPS;      Tout->Branch("cblkPS", &T_cblkPS, "cblkPS/D"); 
  double T_idblkPS;     Tout->Branch("idblkPS", &T_idblkPS, "idblkPS/D"); 
  double T_atimePS;     Tout->Branch("atimePS", &T_atimePS, "atimePS/D"); 
  double T_eSH;         Tout->Branch("eSH", &T_eSH, "eSH/D"); 
  double T_xSH;         Tout->Branch("xSH", &T_xSH, "xSH/D"); 
  double T_ySH;         Tout->Branch("ySH", &T_ySH, "ySH/D"); 
  double T_rblkSH;      Tout->Branch("rblkSH", &T_rblkSH, "rblkSH/D"); 
  double T_cblkSH;      Tout->Branch("cblkSH", &T_cblkSH, "cblkSH/D"); 
  double T_idblkSH;     Tout->Branch("idblkSH", &T_idblkSH, "idblkSH/D"); 
  double T_atimeSH;     Tout->Branch("atimeSH", &T_atimeSH, "atimeSH/D"); 
  //HCAL
  double T_eHCAL;       Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  double T_xHCAL;       Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  double T_yHCAL;       Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
  double T_idblkHCAL;   Tout->Branch("idblkHCAL", &T_idblkHCAL, "idblkHCAL/D"); 
  double T_rblkHCAL;    Tout->Branch("rblkHCAL", &T_rblkHCAL, "rblkHCAL/D"); 
  double T_cblkHCAL ;   Tout->Branch("cblkHCAL", &T_cblkHCAL, "cblkHCAL/D"); 
  double T_atimeHCAL;   Tout->Branch("atimeHCAL", &T_atimeHCAL, "atimeHCAL/D"); 
  double T_tdcHCAL;     Tout->Branch("tdcHCAL", &T_tdcHCAL, "tdcHCAL/D"); 
  double T_xHCAL_exp;   Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); 
  double T_yHCAL_exp;   Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D"); 
  double T_dx;          Tout->Branch("dx", &T_dx, "dx/D"); 
  double T_dy;          Tout->Branch("dy", &T_dy, "dy/D");
  double T_ToF_n;       Tout->Branch("ToF_n", &T_ToF_n, "ToF_n/D");
  //HODO
  int T_ncltmeanHODO;   Tout->Branch("ncltmeanHODO", &T_ncltmeanHODO, "ncltmeanHODO/I"); 
  double T_cltmeanHODO; Tout->Branch("cltmeanHODO", &T_cltmeanHODO, "cltmeanHODO/D"); 
  //coin time trigger
  double T_bbT_trig;    Tout->Branch("bbT_trig", &T_bbT_trig, "bbT_trig/D");
  double T_coinT_trig;  Tout->Branch("coinT_trig", &T_coinT_trig, "coinT_trig/D");

  // reading HCAL cut definitions
  vector<double> dx_p_cut; jmgr->GetVectorFromSubKey<double>(key,"dx_p_cut",dx_p_cut);
  double sbs_kick = abs(dx_p_cut[0]);
  vector<double> dy_p_cut; jmgr->GetVectorFromSubKey<double>(key,"dy_p_cut",dy_p_cut);
  vector<double> dx_n_cut; jmgr->GetVectorFromSubKey<double>(key,"dx_n_cut",dx_n_cut);
  vector<double> dy_n_cut; jmgr->GetVectorFromSubKey<double>(key,"dy_n_cut",dy_n_cut);
  vector<double> hcal_active_area = cut::hcal_active_area_data(1,1); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p_cut[1],dx_n_cut[1],dy_p_cut[1],hcal_active_area);
  TH2F *h2_xyHCAL_p = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_p",sbs_kick);
  TH2F *h2_xyHCAL_n = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_n",0);

  // reading W cut limits
  std::vector<double> W_cutR; jmgr->GetVectorFromSubKey<double>(key,"W_cutR",W_cutR);

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = jmgr->GetValueFromSubKey<double>(key,"hcal_voffset");
  double hcal_hoffset = jmgr->GetValueFromSubKey<double>(key,"hcal_hoffset");
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(),HCAL_axes);
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // EMFF fits
  Ye2017 yefit;
  Kelly2004 kellyfit;
  Seamus20XX seamusfit;
  Christy2022 christyfit;

  // looping through the events ---------------------------------------
  std::cout << std::endl;
  std::vector<double> ElossInTgt; // array to hold energy loss correction values per event
  long nevent=0, nevents=C->GetEntries(), neventsS=S->GetEntries(), index=0, tgevnumS, ngoodevs = 0; 
  int treenum=0, currenttreenum=0; UInt_t runnum=0, nseg, tsegnumS;
  double ebeam=sbsconf.GetEbeam(), ebeam_std=0.; 
  double tdnewcurr=0., tdnewcnt=0; 
  while (C->GetEntry(nevent++)) {

    // progress indicator 
    if (nevent % 1000 == 0) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    // reading matching scaler info per event
    if (get_scaler_info) {
      /* finding 1st scaler event for the current run */
      while (rnumS!=rnum) {
	if (index==neventsS) {
	  std::cout << Form("Run %u | GevNum %u | GevNumS %llu",rnum,gevnum,gevnumS) << std::endl;
	  throw std::runtime_error("S tree index out of bounds! *INVESTIGATE!");
	}
	S->GetEntry(index); index++;
	tsegnumS = segnumS; tgevnumS = gevnumS;
	tdnewcnt = dnewcnt; tdnewcurr = dnewcurr;
      }
      /* finding nearest scaler event for the current T event */
      while (gevnum>gevnumS && rnumS==rnum) {
	if (index==neventsS) {
	  std::cout << Form("Run %u | GevNum %u | GevNumS %llu",rnum,gevnum,gevnumS) << std::endl;
	  throw std::runtime_error("S tree index out of bounds! **INVESTIGATE!");
	}
	tsegnumS = segnumS; tgevnumS = gevnumS;
	tdnewcnt = dnewcnt; tdnewcurr = dnewcurr;
	S->GetEntry(index); index++;
	if (verbose==-2) std::cout << rnum << " " << tgevnumS << " " << gevnum << " " << segnumS << "\n";
      }
    }

    // keep track of run number & tree number
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum; nseg++;
      // apply global cuts efficiently (AJRP method)
      GlobalCut->UpdateFormulaLeaves();
      
      // read ebeam once per run
      if (nevent == 1 || rnum != runnum) {
	runnum = rnum; nseg=1;
	/* In search of a faster algorithm */
	auto it = std::find_if(crun.begin(), crun.end(), [&](CodaRun const& cr) {return cr.runnum == runnum;});
	if (it != crun.end()) {
	  ebeam = it->ebeam; 
	  ebeam_std = it->ebeam_std;
	}else 
	  std::cerr << "**!** Run " << runnum << " is not in spreadsheet!" << std::endl;
     }
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (!passedgCut) continue;
    ngoodevs++;
      
    // coin time cut (N/A for simulation)  !! Not a reliable cut - loosing a lot of elastics
    double bbcal_time=0., hcal_time=0.;
    for(int ihit=0; ihit<tdcElemN; ihit++){
      if(tdcElem[ihit]==5) bbcal_time=tdcTrig[ihit];
      if(tdcElem[ihit]==0) hcal_time=tdcTrig[ihit];
    }
    double coin_time = hcal_time - bbcal_time; 
    T_bbT_trig = bbcal_time;
    T_coinT_trig = coin_time; h_coin_time->Fill(coin_time);

    // constructing the 4 vectors
    /* Reaction    : e + e' -> N + N'
       Conservation: Pe + Peprime = PN + PNprime */
    TLorentzVector Peprime_course{px[0],py[0],pz[0],p[0]};
    util_pd::GetElossInTgt(target,rnum,conf,vz[0],kine::etheta(Peprime_course),verbosefn,ElossInTgt);
    double ebeam_corr = ebeam - ElossInTgt[0];
    double trP_corr = p[0] + ElossInTgt[1];
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e- 4-vector
    TLorentzVector Peprime(px[0] * (trP_corr/p[0]), // scattered e- 4-vector
			   py[0] * (trP_corr/p[0]),
			   pz[0] * (trP_corr/p[0]),
			   trP_corr);                 
    TLorentzVector PN;                              // target nucleon 4-vector
    kine::SetPN(Ntype, PN);
    TLorentzVector PNprime;                         // Recoil nucleon 4-vector
    TLorentzVector q = Pe - Peprime;                // 4-momentum of virtual photon

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pelas = kine::pelas(ebeam_corr, etheta, Ntype);
    double thelas = kine::thelas(ebeam_corr, trP_corr, Ntype);

    double nu = 0.;                   // energy of the virtual photon
    double pN_expect = 0.;            // expected recoil nucleon momentum
    double thetaN_expect = 0.;        // expected recoil nucleon theta
    double phiN_expect = ephi + constant::pi; 
    /* Different modes of calculation. Goal is to achieve the best resolution
       model 0 = uses reconstructed p as independent variable
       model 1 = uses reconstructed angles as independent variable 
       model 2 = uses 4-vector calculation */
    TVector3 pNhat;                   // 3-momentum of the recoil nucleon (Unit)
    double Q2recon{0},W2recon{0},epsilon{0},epsilon_p{0},epsilon_n{0};
    if (model == 0) { // p as independent variable
      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((Pe.E() - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      PNprime.SetPxPyPzE(pN_expect*pNhat.X(), pN_expect*pNhat.Y(), pN_expect*pNhat.Z(), nu+PN.E());
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), thelas);
      epsilon = kine::epsilon_general(thelas,Q2recon,nu);
      epsilon_p = kine::epsilon(thelas,Q2recon,"p");
      epsilon_n = kine::epsilon(thelas,Q2recon,"n");
      //W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
      W2recon = kine::W2_general(Pe.E(), Peprime.E(), etheta, Ntype);
    } else if (model == 1) { // angle as independent variable
      nu = Pe.E() - pelas;
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((Pe.E() - pelas*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      PNprime.SetPxPyPzE(pN_expect*pNhat.X(), pN_expect*pNhat.Y(), pN_expect*pNhat.Z(), nu+PN.E());
      Q2recon = kine::Q2(Pe.E(), pelas, etheta);
      epsilon = kine::epsilon_general(etheta,Q2recon,nu);
      epsilon_p = kine::epsilon(etheta,Q2recon,"p");
      epsilon_n = kine::epsilon(etheta,Q2recon,"n");
      //W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
      W2recon = kine::W2_general(Pe.E(), Peprime.E(), etheta, Ntype);
    } else if (model == 2) { // 4-vector calculation
      nu = q.E();
      PNprime = q + PN;
      pNhat = PNprime.Vect().Unit();
      Q2recon = -q.M2();
      epsilon = kine::epsilon_general(etheta,Q2recon,nu);
      epsilon_p = kine::epsilon(etheta,Q2recon,"p");
      epsilon_n = kine::epsilon(etheta,Q2recon,"n");
      W2recon = PNprime.M2();
    }
    h_Q2->Fill(Q2recon); 
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pelas - 1.0; h_dpel->Fill(dpel);

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_W = Wrecon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pelas = pelas;
    T_epsilon = epsilon;
    T_epsilon_p = epsilon_p;
    T_epsilon_n = epsilon_n;

    T_rnum = rnum;
    T_segnum = nseg;
    T_gevnum = gevnum;
    T_ebeam = ebeam;
    T_ebeam_corr = Pe.E();
    T_ebeam_std = ebeam_std;
    if (get_scaler_info) {
      //T_segnumS = tsegnumS;
      T_dnewcnt = tdnewcnt;
      T_dnewcurr = tdnewcurr;
    }

    T_vz = vz[0];
    T_trP = p[0];
    T_trP_corr = trP_corr;
    T_trX = xTr[0];
    T_trY = yTr[0];
    T_trTh = thTr[0];
    T_trPh = phTr[0];

    T_tgX = xtgt[0];
    T_tgY = ytgt[0];
    T_tgTh = thtgt[0];
    T_tgPh = phtgt[0];

    // EMFF extraction using same parametrization used in MC generators
    double GEp_kelly = kellyfit.GetFF(G_t::kGEp,Q2recon);
    double GMp_kelly = kellyfit.GetFF(G_t::kGMp,Q2recon);
    double GEn_seamus = seamusfit.GetFF(G_t::kGEn,Q2recon);
    double GMn_kelly = kellyfit.GetFF(G_t::kGMn,Q2recon);

    T_sigMott = kine::sigmaMott(Pe.E(),Peprime.E(),etheta);
    T_sigRed_p = kine::sigmaReduced(kine::tau(Q2recon,"p"),epsilon_p,GEp_kelly,GMp_kelly);
    T_sigRed_n = kine::sigmaReduced(kine::tau(Q2recon,"n"),epsilon_n,GEn_seamus,GMn_kelly);
    T_sigBorn_ratio = kine::sigmaBorn_ratio(etheta,Q2recon,GEp_kelly,GMp_kelly,GEn_seamus,GMn_kelly);

    T_ePS = ePS;
    T_rblkPS = rblkPS;
    T_cblkPS = cblkPS;
    T_idblkPS = idblkPS;
    T_atimePS = atimePS;

    T_eSH = eSH;
    T_xSH = xSH;
    T_ySH = ySH;
    T_rblkSH = rblkSH;
    T_cblkSH = cblkSH;
    T_idblkSH = idblkSH;
    T_atimeSH = atimeSH;

    T_eHCAL = eHCAL;
    T_xHCAL = xHCAL;
    T_yHCAL = yHCAL;
    T_rblkHCAL = rblkHCAL;
    T_cblkHCAL = cblkHCAL;
    T_idblkHCAL = idblkHCAL;
    T_atimeHCAL = atimeHCAL;
    T_tdcHCAL = tdcHCAL[0];

    T_ncltmeanHODO = ncltmeanHODO;
    T_cltmeanHODO = cltmeanHODO[0];

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL - xyHCAL_exp[0];  
    double dy = yHCAL - xyHCAL_exp[1]; 

    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_dx = dx;
    T_dy = dy;

    /* Calculating thpq (both p & n hypothesis) */
    // n (no deflection)
    TVector3 HCAL_pos = HCAL_origin + xHCAL*HCAL_axes[0] + yHCAL*HCAL_axes[1];
    TVector3 n_dir = (HCAL_pos - vertex);
    T_thpq_n = acos(n_dir.Unit().Dot(pNhat));
    // p 
    double BdL = (sbsmag / 100.) * 2.; //expconst::sbsmaxfield * expconst::sbsdipolegap;
    double proton_thetabend = 0.3 * BdL / PNprime.Vect().Mag();  // p*theta = 0.3*BdL
    double proton_deflection = tan(proton_thetabend)*(sbsconf.GetHCALdist()-(sbsconf.GetSBSdist()+expconst::sbsdipolegap/2.0));
    TVector3 p_dir = (HCAL_pos + proton_deflection*HCAL_axes[0] - vertex);
    T_thpq_p = acos(p_dir.Unit().Dot(pNhat));

    // calculate ToF for neutrons
    double ToF_n = (n_dir.Mag() / constant::c) * sqrt(1. + pow((constant::Mn/PNprime.Vect().Mag()), 2));
    T_ToF_n = ToF_n*1e9; //ns

    // HCAL active area and safety margin cuts [Fiducial region]
    ARCut = cut::inHCAL_activeA(xHCAL,yHCAL,hcal_active_area);
    SMCut = cut::inHCAL_safety_margin(target,xyHCAL_exp[0],xyHCAL_exp[1],sbs_kick,hcal_safety_margin);
    fiduCut = ARCut && SMCut; 
    // defining HCAL cuts
    pCut = pow((dx-dx_p_cut[0])/(dx_p_cut[1]*dx_p_cut[2]),2) + pow((dy-dy_p_cut[0])/(dy_p_cut[1]*dy_p_cut[2]),2) <= 1.;
    nCut = pow((dx-dx_n_cut[0])/(dx_n_cut[1]*dx_n_cut[2]),2) + pow((dy-dy_n_cut[0])/(dy_n_cut[1]*dy_n_cut[2]),2) <= 1.;
    // a few variations
    pCut_1p5sig = pow((dx-dx_p_cut[0])/(dx_p_cut[1]*dx_p_cut[2]*1.5),2) + pow((dy-dy_p_cut[0])/(dy_p_cut[1]*dy_p_cut[2]*1.5),2) <= 1.;
    nCut_1p5sig = pow((dx-dx_n_cut[0])/(dx_n_cut[1]*dx_n_cut[2]*1.5),2) + pow((dy-dy_n_cut[0])/(dy_n_cut[1]*dy_n_cut[2]*1.5),2) <= 1.;
    pCut_2sig = pow((dx-dx_p_cut[0])/(dx_p_cut[1]*dx_p_cut[2]*2.),2) + pow((dy-dy_p_cut[0])/(dy_p_cut[1]*dy_p_cut[2]*2.),2) <= 1.;
    nCut_2sig = pow((dx-dx_n_cut[0])/(dx_n_cut[1]*dx_n_cut[2]*2.),2) + pow((dy-dy_n_cut[0])/(dy_n_cut[1]*dy_n_cut[2]*2.),2) <= 1.;
    pCut_2p5sig = pow((dx-dx_p_cut[0])/(dx_p_cut[1]*dx_p_cut[2]*2.5),2) + pow((dy-dy_p_cut[0])/(dy_p_cut[1]*dy_p_cut[2]*2.5),2) <= 1.;
    nCut_2p5sig = pow((dx-dx_n_cut[0])/(dx_n_cut[1]*dx_n_cut[2]*2.5),2) + pow((dy-dy_n_cut[0])/(dy_n_cut[1]*dy_n_cut[2]*2.5),2) <= 1.;
    pCut_3sig = pow((dx-dx_p_cut[0])/(dx_p_cut[1]*dx_p_cut[2]*3.),2) + pow((dy-dy_p_cut[0])/(dy_p_cut[1]*dy_p_cut[2]*3.),2) <= 1.;
    nCut_3sig = pow((dx-dx_n_cut[0])/(dx_n_cut[1]*dx_n_cut[2]*3.),2) + pow((dy-dy_n_cut[0])/(dy_n_cut[1]*dy_n_cut[2]*3.),2) <= 1.;

    // defining W cut
    WCut = Wrecon >= W_cutR[0] && Wrecon <= W_cutR[1];

    // W cut
    if (WCut) {
      h_dxHCAL_nfc->Fill(dx);
      h_dyHCAL_nfc->Fill(dy);
      // fiducial cut
      if (fiduCut) {
        h_dxHCAL->Fill(dx);
        h_dyHCAL->Fill(dy);
        h2_rcHCAL->Fill(cblkHCAL, rblkHCAL);
        h2_dxdyHCAL->Fill(dy, dx);

        if (pCut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
        if (nCut) h2_xyHCAL_n->Fill(xyHCAL_exp[1], xyHCAL_exp[0]);
      }
    }

    // fiducial cut but no W cut
    if (Wrecon>0) { 
      if (fiduCut) {
	h_W->Fill(Wrecon);
	if (pCut || nCut) { 
	  h_W_cut->Fill(Wrecon);
	} else {
	  h_W_acut->Fill(Wrecon);
	}
      }  
    }

    Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;

  // calculating total charge analyzed
  double totcharge = util_pd::GetTotCharge(crun);

  /////////////////////////////////
  // Generating diagnostic plots //
  /////////////////////////////////
  TString outPlot = outFile; outPlot.ReplaceAll(".root",".pdf");
  /**** Canvas 1 (p&n spots) ****/
  TCanvas *c1 = util_pd::TC("c1",2,2);
  c1->cd(1); //
  h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p_cut[0], dx_p_cut[0], dy_p_cut[2]*dy_p_cut[1], dx_p_cut[2]*dx_p_cut[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n_cut[0], dx_n_cut[0], dy_n_cut[2]*dy_n_cut[1], dx_n_cut[2]*dx_n_cut[1], 0,360,0);
  c1->cd(2); //
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");
  c1->cd(3); //
  h2_xyHCAL_p->Draw("colz");
  util_pd::DrawArea(hcal_active_area,2,4,9);
  util_pd::DrawArea(hcal_safety_margin,4,4,9);
  c1->cd(4); //
  h2_xyHCAL_n->Draw("colz");
  util_pd::DrawArea(hcal_active_area,2,4,9);
  util_pd::DrawArea(hcal_safety_margin,4,4,9);
  c1->SaveAs(Form("%s[",outPlot.Data())); c1->SaveAs(Form("%s",outPlot.Data())); c1->Write();
  //**** -- ***//

  /**** Canvas 2 (dx & dy) ****/
  TCanvas *c2 = util_pd::TC("c2",2,2);
  std::vector<double> hdxp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_p_fitR", hdxp_fitR);
  std::vector<double> hdxn_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_n_fitR", hdxn_fitR);
  std::vector<double> hdy_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dyHCAL_fitR", hdy_fitR);
  gStyle->SetOptFit(1);
  c2->cd(1); //
  TF1 *fdxp = fit::fit_1gs_nbg(hdxp_fitR,h_dxHCAL);
  TF1 *fdxn = fit::fit_1gs_nbg(hdxn_fitR,h_dxHCAL);
  double dxpM = fdxp->GetParameter(1); double dxpS = fdxp->GetParameter(2);
  double dxnM = fdxn->GetParameter(1); double dxnS = fdxn->GetParameter(2);
  fdxp->Draw("same");
  c2->cd(2); //
  TF1 *fdy = fit::fit_1gs_nbg(hdy_fitR,h_dyHCAL);
  double dyM = fdy->GetParameter(1); double dyS = fdy->GetParameter(2);
  c2->cd(3); //
  TF1 *fdxp_nfc = fit::fit_1gs_nbg(hdxp_fitR,h_dxHCAL_nfc);
  TF1 *fdxn_nfc = fit::fit_1gs_nbg(hdxn_fitR,h_dxHCAL_nfc);
  double dxpM_nfc = fdxp_nfc->GetParameter(1); double dxpS_nfc = fdxp_nfc->GetParameter(2);
  double dxnM_nfc = fdxn_nfc->GetParameter(1); double dxnS_nfc = fdxn_nfc->GetParameter(2);
  fdxp_nfc->Draw("same");
  c2->cd(4); //
  TF1 *fdy_nfc = fit::fit_1gs_nbg(hdy_fitR,h_dyHCAL_nfc);
  double dyM_nfc = fdy->GetParameter(1); double dyS_nfc = fdy->GetParameter(2);
  c2->SaveAs(Form("%s",outPlot.Data())); c2->Write();
  //**** -- ***//

  /**** Summary Canvas ****/
  TCanvas *cSummary = new TCanvas("cSummary","Summary");
  cSummary->cd();
  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText(Form(" Date of creation: %s", util_pd::getDate().c_str()));
  pt->AddText(Form("Configfile: %s",configfilename));
  pt->AddText(Form(" Analyzing QE events for SBS%d-SBS%dp settings",conf,sbsmag));
  pt->AddText(Form(" Analysis model: %d",model));
  pt->AddText(Form(" Total # events analyzed: %ld, Total # runs: %d",nevents,nruns));
  pt->AddText(Form(" Total charge: %.7fC",totcharge));
  pt->AddText(Form(" HCAL offsets: v = %.4f, h = %.4f",hcal_voffset,hcal_hoffset));
  pt->AddText(Form(" Global cuts: "));
  std::string tmpstr = "";
  for (std::size_t i=0; i<gCutList.size(); i++) {
    if (i>0 && i%3==0) {pt->AddText(Form(" %s",tmpstr.c_str())); tmpstr="";}
    tmpstr += gCutList[i] + ", "; 
  }
  if (!tmpstr.empty()) pt->AddText(Form(" %s",tmpstr.c_str()));
  pt->AddText(Form(" # events passed global cuts: %ld",ngoodevs));
  pt->AddText(" Elastic cuts: ");
  pt->AddText(Form(" Inbuilt W cut: %.2f #leq W #leq %.2f GeV/c",W_cutR[0],W_cutR[1]));
  pt->AddText(Form(" Inbuilt p cut (#Deltax): Mean = %.4f, %.1f#sigma = %.4f",dx_p_cut[0],dx_p_cut[2],dx_p_cut[1]));
  pt->AddText(Form(" Inbuilt p cut (#Deltay): Mean = %.4f, %.1f#sigma = %.4f",dy_p_cut[0],dy_p_cut[2],dy_p_cut[1]));
  pt->AddText(Form(" Inbuilt n cut (#Deltax): Mean = %.4f, %.1f#sigma = %.4f",dx_n_cut[0],dx_n_cut[2],dx_n_cut[1]));
  pt->AddText(Form(" Inbuilt n cut (#Deltay): Mean = %.4f, %.1f#sigma = %.4f",dy_n_cut[0],dy_n_cut[2],dy_n_cut[1]));
  pt->AddText(" Fit info: ");
  pt->AddText(" p & n peaks, w/ fiducial cut: dxpM,dxpS,dxnM,dxnS,dyM,dyS ");
  pt->AddText(Form(" %.5f,%.5f,%.5f,%.5f,%.5f,%.5f",dxpM,dxpS,dxnM,dxnS,dyM,dyS));
  pt->AddText(" p & n peaks, w/o fiducial cut: dxpM_nfc,dxpS_nfc,dxnM_nfc,dxnS_nfc,dyM_nfc,dyS_nfc ");
  pt->AddText(Form(" %.5f,%.5f,%.5f,%.5f,%.5f,%.5f",dxpM_nfc,dxpS_nfc,dxnM_nfc,dxnS_nfc,dyM_nfc,dyS_nfc));
  sw->Stop();
  pt->AddText(Form("Macro processing time: CPU %.1fs | Real %.1fs",sw->CpuTime(),sw->RealTime()));
  TText *t1 = pt->GetLineWith("Configfile"); t1->SetTextColor(kRed);
  TText *t2 = pt->GetLineWith(" Global"); t2->SetTextColor(kBlue);
  TText *t3 = pt->GetLineWith(" Elastic"); t3->SetTextColor(kBlue);
  TText *t4 = pt->GetLineWith(" Fit info"); t4->SetTextColor(kBlue);
  TText *t5 = pt->GetLineWith("Macro"); t5->SetTextColor(kGreen+3);
  pt->Draw(); 
  cSummary->SaveAs(Form("%s",outPlot.Data())); cSummary->SaveAs(Form("%s]",outPlot.Data())); cSummary->Write();  
  //**** -- ***//

  std::cout << "\n----Fit info----" << "\n";
  std::cout << "dxpM,dxpS,dxnM,dxnS,dyM,dyS,dxpM_nfc,dxpS_nfc,dxnM_nfc,dxnS_nfc,dyM_nfc,dyS_nfc" << "\n";
  std::cout << dxpM<<","<<dxpS<<","<<dxnM<<","<<dxnS<<","<<dyM<<","<<dyS<<","<<dxpM_nfc<<","<<dxpS_nfc<<","<<dxnM_nfc<<","<<dxnS_nfc<<","<<dyM_nfc<<","<<dyS_nfc<<"\n";
  std::cout << "----------------" << "\n\n";

  std::cout << "------" << std::endl;
  std::cout << " Total charge : " << totcharge << " C" << std::endl;
  std::cout << " Summary plots  : " << outPlot << std::endl;
  std::cout << " Output ROOT file  : " << outFile << std::endl;
  std::cout << "------" << std::endl << std::endl;

  std::cout << "CPU time = " << sw->CpuTime() << "s. Real time = " << sw->RealTime() << "s.\n\n";

  Tout->Write("",TObject::kOverwrite);
  h_Q2->Write();
  h_dpel->Write(); h_W->Write();
  h_W_cut->Write(); h_W_acut->Write();
  h_dxHCAL->Write(); h_dyHCAL->Write();
  h2_rcHCAL->Write(); h2_dxdyHCAL->Write();
  h2_xyHCAL_p->Write(); h2_xyHCAL_n->Write();
  h_coin_time->Write();
  sw->Delete();
  delete jmgr;
  return 0;
}

/*
  ////////////////////////////////////////////////////////
  // Brief description of configuration file parameters //
  ////////////////////////////////////////////////////////
  ** Ntype : Struck nucleon type. (Valid options: n, p, np)
     - n(p) => neutrons(protons) 
     - np => Avg. of n and p masses. Use for LD2 data since the struck nucleon is not known apriori.
  ** SBS_config : SBS configuration. (Valid options: 4,7,11,14,8,9)
  ** SBS_magnet_percent : SBS magnet field current as a percentage of 2100A
  ** model : Model of analysis. (Valid options: 0, 1, 2)
     - 0 => uses reconstructed p as independent variable
     - 1 => uses reconstructed angles as independent variable
     - 2 => uses 4-vector calculation
  ** Nruns_to_ana : # CODA runs to analyze
  ** get_scaler_info : If true, writes out matching scaler tree variables in the output ROOT tree
  ** rootfile_dir : Directory name w/ path containing the MC ROOT files to analyze
  ** global_cut : set of global cuts to apply at the start of event processing
  ** SBS_field : 
  ** hcal_v(h)offset : 
  ** dx_p(n)_cut : deltax p(n) peak cut definitions. Needed to constitute p(n) spot cuts. 
     - dx_p(n)_cut[0] => p(n) peak position, 
     - dx_p(n)_cut[1] => p(n) peak RMS, 
     - dx_p(n)_cut[2] => # sigma to include in the cut
  ** dy_p(n)_cut : deltay p(n) peak cut definitions. Needed to constitute p(n) spot cuts.
     - dy_p(n)_cut[0] => p(n) peak position, 
     - dy_p(n)_cut[1] => p(n) peak RMS, 
     - dy_p(n)_cut[2] => # sigma to include in the cut
  ** h_dx(dy)HCAL_lims : h_dx(dy)HCAL histogram limits (Can be found in the output ROOT file)
     - h_dx(dy)HCAL_lims[0] : No. of bins of h_dx(dy)HCAL histograms
     - h_dx(dy)HCAL_lims[1] : xmin
     - h_dx(dy)HCAL_lims[2] : xmax
  ** h_dxHCAL_p(n)_fitR : Fit ranges for the proton(neutron) signal peak in deltax dist. (h_dxHCAL histogram)
                          Algorithm fits the distribution twice for optimization.
     - h_dxHCAL_p(n)_fitR[0] : xmin for 1st fit (crude). Try to avoid any secondary peak.
     - h_dxHCAL_p(n)_fitR[1] : xmax for 1st fit (crude). Try to avoid any secondary peak.
     - h_dxHCAL_p(n)_fitR[2] : # sigma below the peak for 2nd fit (fine)
     - h_dxHCAL_p(n)_fitR[3] : # sigma above the peak for 2nd fit (fine)
  ** h_dyHCAL_fitR : Fit ranges for the signal peak in deltay dist. (h_dxHCAL histogram)
                     Algorithm fits the distribution twice for optimization.
     - h_dyHCAL_fitR[0] : xmin for 1st fit (crude). Try to avoid any secondary peak.
     - h_dyHCAL_fitR[1] : xmax for 1st fit (crude). Try to avoid any secondary peak.
     - h_dyHCAL_fitR[2] : # sigma below the peak for 2nd fit (fine)
     - h_dyHCAL_fitR[3] : # sigma above the peak for 2nd fit (fine)
  ** W_cutR : W cut range.
     - W_cutR[0] : lower limit
     - W_cutR[1] : upper limit
*/
