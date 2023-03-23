/* 
   This macro will perform QE analysis for GMn using LD2 data.
   E.g. config. file: sbs14-sbs70p/conf_qelas_ana_data.json
   -----
   P. Datta  Created  11-02-2022 
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

#include "../include/gmn-ana.h"
#include "../dflay/src/JSONManager.cxx"

/* this script will only analyze LD2 data */
static const std::string target = "LD2";

int qelas_ana_data (const char *configfilename,
                    int verbose=-1,  //<0=>Debug
                    int verbosefn=0, //>0=>Debug
                    std::string filebase="pdout/test_qelas_ana_data")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent"); 
  int pass = jmgr->GetValueFromKey<int>("replay_pass"); 
  SBSconfig sbsconf(conf, sbsmag);
  cout << sbsconf;

  // reading run info and parsing ROOT trees
  std::string runsheet_dir = jmgr->GetValueFromKey_str("runsheet_dir");
  int nruns = jmgr->GetValueFromKey<int>("Nruns_to_ana"); // # of runs to analyze
  vector<CodaRun> crun; util_pd::ReadRunList(runsheet_dir,nruns,conf,target,pass,sbsmag,verbosefn,crun);
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  TChain *C = new TChain("T"); util_pd::LoadROOTTree(rootfile_dir,crun,1,verbosefn,C); 
 
  // reading scaler tree
  TChain *S = new TChain("Tout"); 
  int get_scaler_info = jmgr->GetValueFromKey<int>("get_scaler_info");
  if (get_scaler_info) {
    S->Add(Form("epics/epout/scalerdata_prun_SBS%d_%s.root",conf,target.c_str()));
    if (S->GetEntries()==0) throw std::runtime_error("No scaler event found!");
  }

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = jmgr->GetValueFromKey<int>("model");
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");

  // setting up global cuts
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = gcut.c_str();
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy 
  // double HALLA_p; setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  // bbcal clus var
  double eSH, xSH, ySH, rblkSH, cblkSH, idblkSH, atimeSH, ePS, rblkPS, cblkPS, idblkPS, atimePS;
  std::vector<std::string> bbcalclvar = {"sh.e","sh.x","sh.y","sh.rowblk","sh.colblk","sh.idblk","sh.atimeblk",
					 "ps.e","ps.rowblk","ps.colblk","ps.idblk","ps.atimeblk"};
  std::vector<void*> bbcalclvar_mem = {&eSH,&xSH,&ySH,&rblkSH,&cblkSH,&idblkSH,&atimeSH,&ePS,&rblkPS,&cblkPS,&idblkPS,&atimePS};
  setrootvar::setbranch(C, "bb", bbcalclvar, bbcalclvar_mem);
 
  // hcal clus var
  double eHCAL, xHCAL, yHCAL, rblkHCAL, cblkHCAL, idblkHCAL, atimeHCAL, tdcHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","atimeblk","tdctimeblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&atimeHCAL,&tdcHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

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
  if (verbose==0 && verbosefn==0) filebase = "pdout/qelas_ana_data";
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d_pass%d.root", 
			 filebase.c_str(), sbsconf.GetSBSconf(), sbsconf.GetSBSmag(), model, pass);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms
  TH1F *h_W = util_pd::TH1FhW("h_W");
  TH1F *h_W_cut = util_pd::TH1FhW("h_W_cut");
  TH1F *h_W_acut = util_pd::TH1FhW("h_W_acut");
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = util_pd::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim; jmgr->GetVectorFromKey<double>("h_dxHCAL_lims", hdx_lim);
  vector<double> hdy_lim; jmgr->GetVectorFromKey<double>("h_dyHCAL_lims", hdy_lim);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","W & fiducial cuts;x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","W & fiducial cuts;y_{HCAL} - y_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_coin_time = new TH1F("h_coin_time", "Coincidence time (ns)", 200, 380, 660);

  TH2F *h2_rcHCAL = util_pd::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = util_pd::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_xyHCAL_p = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_p");
  TH2F *h2_xyHCAL_n = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_n");

  // defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  //cuts
  bool WCut;            Tout->Branch("WCut", &WCut, "WCut/O");
  bool pCut;            Tout->Branch("pCut", &pCut, "pCut/O");
  bool nCut;            Tout->Branch("nCut", &nCut, "nCut/O");
  bool fiduCut;         Tout->Branch("fiduCut", &fiduCut, "fiduCut/O");
  //run info
  UInt_t T_rnum;        Tout->Branch("rnum", &T_rnum, "rnum/i");
  UInt_t T_segnum;      Tout->Branch("segnum", &T_segnum, "segnum/i");
  ULong64_t T_gevnum;   Tout->Branch("gevnum", &T_gevnum, "gevnum/l");
  double T_ebeam;       Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
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
  double T_pcentral;    Tout->Branch("pcentral", &T_pcentral, "pcentral/D");
  double T_thetapq_p;   Tout->Branch("thetapq_p", &T_thetapq_p, "thetapq_p/D");
  double T_thetapq_n;   Tout->Branch("thetapq_n", &T_thetapq_n, "thetapq_n/D");
  //track
  double T_vz;          Tout->Branch("vz", &T_vz, "vz/D");
  double T_trP;         Tout->Branch("trP", &T_trP, "trP/D");
  double T_trX;         Tout->Branch("trX", &T_trX, "trX/D");
  double T_trY;         Tout->Branch("trY", &T_trY, "trY/D");
  double T_trTh;        Tout->Branch("trTh", &T_trTh, "trTh/D");
  double T_trPh;        Tout->Branch("trPh", &T_trPh, "trPh/D");
  double T_tgX;         Tout->Branch("tgX", &T_tgX, "tgX/D");
  double T_tgY;         Tout->Branch("tgY", &T_tgY, "tgY/D");
  double T_tgTh;        Tout->Branch("tgTh", &T_tgTh, "tgTh/D");
  double T_tgPh;        Tout->Branch("tgPh", &T_tgPh, "tgPh/D");
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

  // Do the energy loss calculation here ...........

  // reading HCAL cut definitions
  vector<double> dx_p; jmgr->GetVectorFromKey<double>("dx_p", dx_p);
  double sbs_kick = abs(dx_p[0]);
  vector<double> dy_p; jmgr->GetVectorFromKey<double>("dy_p", dy_p);
  double Nsigma_cut_dx_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_p");
  double Nsigma_cut_dy_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_p");
  vector<double> dx_n; jmgr->GetVectorFromKey<double>("dx_n", dx_n);
  vector<double> dy_n; jmgr->GetVectorFromKey<double>("dy_n", dy_n);
  double Nsigma_cut_dx_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_n");
  double Nsigma_cut_dy_n = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_n");
  vector<double> hcal_active_area = cut::hcal_active_area_data(); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p[1], dx_n[1], dy_p[1], hcal_active_area);

  // reading W cut limits
  double Wmin = jmgr->GetValueFromKey<double>("Wmin");
  double Wmax = jmgr->GetValueFromKey<double>("Wmax");

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = jmgr->GetValueFromKey<double>("hcal_voffset");
  double hcal_hoffset = jmgr->GetValueFromKey<double>("hcal_hoffset");
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // looping through the events ---------------------------------------
  std::cout << std::endl;
  long nevent=0, nevents=C->GetEntries(), neventsS=S->GetEntries(), index=0, tgevnumS; 
  int treenum=0, currenttreenum=0; UInt_t runnum=0, nseg, tsegnumS;
  double ebeam=sbsconf.GetEbeam(), ebeam_std=0.; 
  double tdnewcurr=0., tdnewcnt=0; 
  while (C->GetEntry(nevent++)) {

    // progress indicator 
    if (nevent % 1000 == 0) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    // reading matching scaler info per event
    if (get_scaler_info) {
        // finding 1st scaler event for the current run
        if (nevent==1 || rnumS!=rnum) {
          while (rnumS!=rnum && index<neventsS) {
            S->GetEntry(index); index++;
	    tsegnumS = segnumS; tgevnumS = gevnumS;
            tdnewcnt = dnewcnt; tdnewcurr = dnewcurr;
          }
        }   
        // finding nearest scaler event for the current T event
        while (gevnum>gevnumS && rnumS==rnum && index<neventsS) {
	  tsegnumS = segnumS; tgevnumS = gevnumS;
          tdnewcnt = dnewcnt; tdnewcurr = dnewcurr;
          S->GetEntry(index); index++;
          if (verbose==-2) std::cout << tgevnumS << " " << gevnum << " " << segnumS << std::endl;
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
	auto it = std::find_if(crun.begin(), crun.end(), [=](CodaRun const& cr) {return cr.runnum == runnum;});
	if (it != crun.end()) {
	  ebeam = it->ebeam; 
	  ebeam_std = it->ebeam_std;
	}else 
	  std::cerr << "**!** Run " << runnum << " is not in spreadsheet!" << std::endl;
     }
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (!passedgCut) continue;
      
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
    double ebeam_corr = ebeam; //- MeanEloss;
    double precon = p[0]; //+ MeanEloss_outgoing;
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e- 4-vector
    TLorentzVector Peprime(px[0] * (precon/p[0]),   // scattered e- 4-vector
			   py[0] * (precon/p[0]),
			   pz[0] * (precon/p[0]),
			   precon);                 
    TLorentzVector PN;                              // target nucleon 4-vector
    kine::SetPN(Ntype, PN);
    TLorentzVector PNprime;                         // Recoil nucleon 4-vector
    TLorentzVector q = Pe - Peprime;                // 4-momentum of virtual photon

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pcentral = kine::pcentral(ebeam_corr, etheta, Ntype);

    double nu = 0.;                   // energy of the virtual photon
    double pN_expect = 0.;            // expected recoil nucleon momentum
    double thetaN_expect = 0.;        // expected recoil nucleon theta
    double phiN_expect = ephi + constant::pi; 
    /* Different modes of calculation. Goal is to achieve the best resolution
       model 0 = uses reconstructed p as independent variable
       model 1 = uses reconstructed angles as independent variable 
       model 2 = uses 4-vector calculation */
    TVector3 pNhat;                   // 3-momentum of the recoil nucleon (Unit)
    double Q2recon = 0., W2recon = 0.;
    if (model == 0) {
      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((Pe.E() - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      PNprime.SetPxPyPzE(pN_expect*pNhat.X(), pN_expect*pNhat.Y(), pN_expect*pNhat.Z(), nu+PN.E());
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
    } else if (model == 1) {
      nu = Pe.E() - pcentral;
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((Pe.E() - pcentral*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      PNprime.SetPxPyPzE(pN_expect*pNhat.X(), pN_expect*pNhat.Y(), pN_expect*pNhat.Z(), nu+PN.E());
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
    } else if (model == 2) {
      nu = q.E();
      PNprime = q + PN;
      pNhat = PNprime.Vect().Unit();
      Q2recon = -q.M2();
      W2recon = PNprime.M2();
    }
    h_Q2->Fill(Q2recon); 
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pcentral - 1.0; h_dpel->Fill(dpel);

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_W = Wrecon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pcentral = pcentral;

    T_rnum = rnum;
    T_segnum = nseg;
    T_gevnum = gevnum;
    T_ebeam = Pe.E();
    T_ebeam_std = ebeam_std;
    if (get_scaler_info) {
      //T_segnumS = tsegnumS;
      T_dnewcnt = tdnewcnt;
      T_dnewcurr = tdnewcurr;
    }

    T_vz = vz[0];
    T_trP = p[0];
    T_trX = xTr[0];
    T_trY = yTr[0];
    T_trTh = thTr[0];
    T_trPh = phTr[0];

    T_tgX = xtgt[0];
    T_tgY = ytgt[0];
    T_tgTh = thtgt[0];
    T_tgPh = phtgt[0];

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
    T_tdcHCAL = tdcHCAL;

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

    /* Calculating thetapq (both p & n hypothesis) */
    // n (no deflection)
    TVector3 HCAL_pos = HCAL_origin + xHCAL*HCAL_axes[0] + yHCAL*HCAL_axes[1];
    TVector3 n_dir = (HCAL_pos - vertex);
    T_thetapq_n = acos(n_dir.Unit().Dot(pNhat));
    // p 
    double BdL = (sbsmag / 100.) * expconst::sbsmaxfield * expconst::sbsdipolegap;
    double proton_thetabend = 0.3 * BdL / PNprime.Vect().Mag();  // p*theta = 0.3*BdL
    double proton_deflection = tan(proton_thetabend)*(sbsconf.GetHCALdist() - (sbsconf.GetSBSdist() + expconst::sbsdipolegap/2.0));
    TVector3 p_dir = (HCAL_pos + proton_deflection*HCAL_axes[0] - vertex);
    T_thetapq_p = acos(p_dir.Unit().Dot(pNhat));

    // calculate ToF for neutrons
    double ToF_n = (n_dir.Mag() / constant::c) * sqrt(1. + pow((constant::Mn/PNprime.Vect().Mag()), 2));
    T_ToF_n = ToF_n*1e9; //ns

    // HCAL active area and safety margin cuts [Fiducial region]
    bool AR_cut = cut::inHCAL_activeA(xHCAL, yHCAL, hcal_active_area);
    bool FR_cut = cut::inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], sbs_kick, hcal_safety_margin);
    fiduCut = AR_cut && FR_cut;
    // defining HCAL cuts
    pCut = pow((dx-dx_p[0]) / (dx_p[1]*Nsigma_cut_dx_p), 2) + pow((dy-dy_p[0]) / (dy_p[1]*Nsigma_cut_dy_p), 2) <= 1.;
    nCut = pow((dx-dx_n[0]) / (dx_n[1]*Nsigma_cut_dx_n), 2) + pow((dy-dy_n[0]) / (dy_n[1]*Nsigma_cut_dy_n), 2) <= 1.;
    // defining W cut
    WCut = Wrecon >= Wmin && Wrecon <= Wmax;

    // W cut
    if (WCut) {
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
    if (fiduCut) {
      h_W->Fill(Wrecon);
      if (pCut || nCut) { 
        h_W_cut->Fill(Wrecon);
      } else {
        h_W_acut->Fill(Wrecon);
      }
    }  

    Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;

  // calculating total charge analyzed
  double totcharge = util_pd::GetTotCharge(crun);
    
  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
  c1->Divide(2,2);

  c1->cd(1); h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p[0], dx_p[0], Nsigma_cut_dy_p*dy_p[1], Nsigma_cut_dx_p*dx_p[1], 0,360,0);
  TEllipse Ep_n;
  Ep_n.SetFillStyle(0); Ep_n.SetLineColor(3); Ep_n.SetLineWidth(2);
  Ep_n.DrawEllipse(dy_n[0], dx_n[0], Nsigma_cut_dy_n*dy_n[1], Nsigma_cut_dx_n*dx_n[1], 0,360,0);
 
  c1->cd(2);
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");

  c1->cd(3);
  h2_xyHCAL_p->Draw("colz");
  util_pd::DrawArea(hcal_active_area);
  util_pd::DrawArea(hcal_safety_margin,4);

  c1->cd(4); 
  h2_xyHCAL_n->Draw("colz");
  util_pd::DrawArea(hcal_active_area);
  util_pd::DrawArea(hcal_safety_margin,4);

  // let's record the summary
  TCanvas *c2 = new TCanvas("c2","Summary");
  c2->cd();

  TPaveText *pt = new TPaveText(.05,.1,.95,.8);
  pt->AddText(Form("Configfile: %s",configfilename));
  pt->AddText(Form(" Analysis model: %d",model));
  pt->AddText(Form(" Total charge : %f C",totcharge));
  pt->AddText(Form(" Total # runs analyzed: %d",nruns));
  pt->AddText(Form(" Total # events analyzed: %ld",nevents));
  pt->AddText(Form(" First run no.: %d | Last run no.: %d",crun[0].runnum,crun[nruns-1].runnum));
  pt->AddText(Form(" HCAL offsets: v = %.4f, h = %.4f",hcal_voffset,hcal_hoffset));
  pt->AddText(Form(" Global cuts: %s",gcut.c_str()));
  pt->AddText(Form(" Inbuilt W cut: %.2f <= W <= %.2f GeV/c",Wmin,Wmax));
  pt->AddText(Form(" Inbuilt p cut (dx): mean = %.4f, sigma = %.4f, Nsigma = %.1f",dx_p[0],dx_p[1],Nsigma_cut_dx_p));
  pt->AddText(Form(" Inbuilt p cut (dy): mean = %.4f, sigma = %.4f, Nsigma = %.1f",dy_p[0],dy_p[1],Nsigma_cut_dy_p));
  pt->AddText(Form(" Inbuilt n cut (dx): mean = %.4f, sigma = %.4f, Nsigma = %.1f",dx_n[0],dx_n[1],Nsigma_cut_dx_n));
  pt->AddText(Form(" Inbuilt n cut (dy): mean = %.4f, sigma = %.4f, Nsigma = %.1f",dy_n[0],dy_n[1],Nsigma_cut_dy_n));
  TText *t1 = pt->GetLineWith("Configfile");
  t1->SetTextColor(kBlue);
  pt->Draw();

  // outFile.ReplaceAll(".root",".png");
  // c1->Print(outFile.Data(),"png");

  std::cout << "------" << std::endl;
  std::cout << " Total charge : " << totcharge << " C" << std::endl;
  std::cout << " Output file  : " << outFile << std::endl;
  std::cout << "------" << std::endl << std::endl;

  sw->Stop();
  std::cout << "CPU time elapsed = " << sw->CpuTime() 
	    << " s. Real time = " << sw->RealTime() << " s. " << std::endl << std::endl;

  c1->Write();
  c2->Write();
  h_W->Write();
  h_W_cut->Write();
  h_W_acut->Write();
  h_dpel->Write();
  h_Q2->Write();
  h_dxHCAL->Write();
  h_dyHCAL->Write();
  h2_rcHCAL->Write();
  h2_dxdyHCAL->Write();
  h2_xyHCAL_p->Write();
  h2_xyHCAL_n->Write();
  h_coin_time->Write();
  Tout->Write();
  sw->Delete();
  delete jmgr;
  return 0;
}
