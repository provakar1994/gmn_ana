/* 
   This macro will perform elastic analysis for GMn using LH2 data.
   E.g. Config. File: sbs14-sbs70p/conf_elas_ana_data.json
   * A brief description of all the config. file parameters can be
   found at the bottom of this script.
   -----
   P. Datta  Created  02-01-2023 
*/

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

/* this script will only analyze LH2 data */
static const std::string target = "LH2";

int elas_ana_data (const char *configfilename,
		   int pass, //replay pass
		   int model=2) //Analysis model
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  std::string key = "pass" + std::to_string(pass) + "_model" + std::to_string(model);

  // setting verbosity
  int verbose = jmgr->GetValueFromSubKey<int>(key,"verbose");
  int verbosefn = jmgr->GetValueFromSubKey<int>(key,"verbose_function");

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromSubKey<int>(key,"SBS_config");
  int sbsmag = jmgr->GetValueFromSubKey<int>(key,"SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  cout << sbsconf;

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
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

  // Choosing the best HCAL cluster
  /*
    Algorithms:
    1. Default: Cluster w/ highest energy.
    2. In-time: HE cluster passing HCAL/SH ADC coincidence time. 
    3. Smallest thpq: Clusters with smallest thpq_p that passed ADC coin time and sampling fraction cuts
    -----
    How-to access in Tout:
    1. Use HCAL_acl variables with index = 0.
    2. Use any default HCAL variable.
    3. Use HCAL_acl variables with index = idclHCAL_sthpq_p.
    -----
    Flag(s):
    > hcal_acl_ON: False=>Don't loop through all the HCAL cls. Keep using the default HE clusters.
    -----
    NOTE:
    > HE clusters ("Default algo") always have index=0.
    > While using "Atime algo", for events that have no cluster in time we fall back to index=0.
    > While using the 3rd algorithm, for events that have no cls. passing coin. time cut and sampling
    fraction cut we update the index according to "Atime_algo". 
    > For the 3rd algo, one can choose smallest min(thpq_p,thpq_n) instead of smallest thpq_p. But the 
    later gave better results. 
  */
  bool hcal_acl_ON = jmgr->GetValueFromSubKey<int>(key,"hcal_acl_ON");
  double hcal_sF_cutR = jmgr->GetValueFromSubKey<double>(key,"hcal_samp_frac_cutR");

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromSubKey_str(key,"Ntype");

  // setting up global cuts
  std::string gcut = jmgr->GetValueFromSubKey_str(key,"global_cut");
  std::vector<std::string> gCutList; util_pd::SplitString('&',gcut,gCutList);
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut",(TCut)gcut.c_str(),C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr = jmgr->GetValueFromSubKey<int>(key,"max_N_tracks");
  C->SetBranchStatus("*",0);
  // beam energy - Probably we should take an average over 100 events
  // double HALLA_p;
  // setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

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
  double eHCAL, xHCAL, yHCAL, indexHCAL, rblkHCAL, cblkHCAL, idblkHCAL, atimeHCAL, tdcHCAL[maxNtr];
  std::vector<std::string> hcalclvar = {"e","x","y","index","rowblk","colblk","idblk","atimeblk","clus_blk.tdctime"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&indexHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&atimeHCAL,&tdcHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // hcal all clus vars
  int maxNHCALcl = jmgr->GetValueFromSubKey<int>(key,"max_N_HCAL_clusters");
  int idblkHCAL_aclN; double idblkHCAL_acl[maxNHCALcl], rblkHCAL_acl[maxNHCALcl], cblkHCAL_acl[maxNHCALcl];
  double nblkHCAL_acl[maxNHCALcl], eblkHCAL_acl[maxNHCALcl], atimeblkHCAL_acl[maxNHCALcl], tdcblkHCAL_acl[maxNHCALcl]; 
  double eHCAL_acl[maxNHCALcl], xHCAL_acl[maxNHCALcl], yHCAL_acl[maxNHCALcl];
  std::vector<std::string> hcalclvar_acl = {"id","id","nblk","eblk","e","x","y","row","col","atime","tdctime"};
  std::vector<void*> hcalclvar_acl_mem = {&idblkHCAL_acl,&idblkHCAL_aclN,&nblkHCAL_acl,&eblkHCAL_acl,&eHCAL_acl,
					  &xHCAL_acl,&yHCAL_acl,&rblkHCAL_acl,&cblkHCAL_acl,&atimeblkHCAL_acl,&tdcblkHCAL_acl};
  setrootvar::setbranch(C, "sbs.hcal.clus", hcalclvar_acl, hcalclvar_acl_mem, 1);

  // bbhodo clus var
  int ncltmeanHODO; 
  double cltmeanHODO[maxNtr];
  std::vector<std::string> hodoclvar = {"clus.tmean","clus.tmean"};
  std::vector<void*> hodoclvar_mem = {&ncltmeanHODO,&cltmeanHODO};
  setrootvar::setbranch(C, "bb.hodotdc", hodoclvar, hodoclvar_mem, 0);  

  // track var
  double ntrack,p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr],xfp[maxNtr],yfp[maxNtr],thfp[maxNtr],phfp[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph",
				    "r_x","r_y","r_th","r_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&xTr,&yTr,&thTr,&phTr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt,
				  &xfp,&yfp,&thfp,&phfp};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  // GEM variables
  double nhitsGEM[maxNtr], ngoodhitsGEM[maxNtr], trchi2ndf[maxNtr];
  std::vector<std::string> gemvar = {"nhits","ngoodhits","chi2ndf"};
  std::vector<void*> gemvar_mem = {&nhitsGEM,&ngoodhitsGEM,&trchi2ndf};
  setrootvar::setbranch(C,"bb.gem.track",gemvar,gemvar_mem);

  // GRINCH cluster variables
  double clsizeGRINCH, cltmeanGRINCH, cltotmeanGRINCH, cltrindexGRINCH;
  std::vector<std::string> grinchvar = {"size","t_mean","tot_mean","trackindex"};
  std::vector<void*> grinchvar_mem = {&clsizeGRINCH,&cltmeanGRINCH,&cltotmeanGRINCH,&cltrindexGRINCH};
  setrootvar::setbranch(C,"bb.grinch_tdc.clus",grinchvar,grinchvar_mem);

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

  // scaler tree variables
  UInt_t rnumS=0, segnumS;
  double dnewcnt, dnewcurr;
  ULong64_t evindex, gevnumS;
  std::vector<std::string> streevar = {"rnum","segnum","gevnum","evindex","dnewcnt","dnewcurr"};
  std::vector<void*> streevar_mem = {&rnumS,&segnumS,&gevnumS,&evindex,&dnewcnt,&dnewcurr};
  if (get_scaler_info) setrootvar::setbranch(S,"",streevar,streevar_mem);

  // defining the outputfile
  std::string filebase = jmgr->GetValueFromSubKey_str(key,"outfile_prefix");
  filebase = filebase.empty() ? "" : filebase + "_";
  TString outFile = Form("pdout/%selas_ana_data_sbs%d_sbs%dp_model%d_pass%d.root",filebase.c_str(),conf,sbsmag,model,pass);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms
  TH1F *h_W = util_pd::TH1FhW("h_W");
  TH1F *h_W_cut = util_pd::TH1FhW("h_W_cut");
  TH1F *h_W_acut = util_pd::TH1FhW("h_W_acut");
  TH1F *h_W2_cut = new TH1F("h_W2_cut","",200,-1,4);
  TH1F *h_dpel = new TH1F("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = util_pd::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_lims", hdx_lim);
  vector<double> hdy_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dyHCAL_lims", hdy_lim);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","W & fiducial cuts;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dxHCAL_nfc = new TH1F("h_dxHCAL_nfc","W cut;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","W & fiducial cuts;y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_dyHCAL_nfc = new TH1F("h_dyHCAL_nfc","W cut;y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_coinT_trig = new TH1F("h_coinT_trig","BBCAL-HCAL trigger coincidence time (ns)",200,380,660);
  TH1F *h_coinT_ADC = new TH1F("h_coinT_ADC","HCAL-SH ADC time",200,-20,20);
  TH1F *h_coinT_ADC_c = new TH1F("h_coinT_ADC_c","HCAL-SH ADC time (Corrected)",200,-20,20);

  TH1F *h_dx_w_p_def = new TH1F("h_dx_w_p_def",";dx+p_def (m);",200,-1,1);
  
  TH2F *h2_rcHCAL = util_pd::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = util_pd::TH2FdxdyHCAL("h2_dxdyHCAL");

  TH1F *h_hcl_inTime_idcl = new TH1F("h_hcl_inTime_idcl","Index of HCAL cl. in time (HCAL/SH ADC coin)",maxNHCALcl,0,maxNHCALcl);
  TH1F *h_hcl_inTime_idcl_WCut = new TH1F("h_hcl_inTime_idcl_WCut","Index of HCAL cl. in time (HCAL/SH ADC coin) | WCut",maxNHCALcl,0,maxNHCALcl);

  TH2F *h2_hclHE_eng_vs_idcl = new TH2F("h2_hclHE_eng_vs_idcl","Eng. of HE block vs cl. index | HCAL",maxNHCALcl,0,maxNHCALcl,200,0,1);
  TH2F *h2_hclHE_atime_vs_idcl = new TH2F("h2_hclHE_atime_vs_idcl","ADC time of HE block vs cl. index | HCAL",maxNHCALcl,0,maxNHCALcl,250,25,75);
  TH2F *h2_hclHE_tdc_vs_idcl = new TH2F("h2_hclHE_tdc_vs_idcl","TDC of HE block vs cl. index | HCAL",maxNHCALcl,0,maxNHCALcl,350,-100,-30);

  TH2F *h2_hcl_eng_vs_idcl = new TH2F("h2_hcl_eng_vs_idcl","Cl. energy vs cl. index | HCAL",maxNHCALcl,0,maxNHCALcl,200,0,1);
  TH2F *h2_hcl_nblk_vs_idcl = new TH2F("h2_hcl_nblk_vs_idcl","Cl. multiplicity vs cl. index | HCAL",maxNHCALcl,0,maxNHCALcl,10,0,10);
 
  // Defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  Tout->SetMaxTreeSize(4000000000LL);
  //cuts
  bool WCut;              Tout->Branch("WCut", &WCut, "WCut/O");
  bool bbfiduCut;         Tout->Branch("bbfiduCut", &bbfiduCut, "bbfiduCut/O");
  bool pCut;              Tout->Branch("pCut", &pCut, "pCut/O");
  double pdx_nS;          Tout->Branch("pdx_nS", &pdx_nS, "pdx_nS/D"); //# sigma away from p dx peak
  double dy_nS;           Tout->Branch("dy_nS", &dy_nS, "dy_nS/D"); //# sigma away from dy peak
  // --
  bool SMCut;             Tout->Branch("SMCut", &SMCut, "SMCut/O");
  double SMy_nS;          Tout->Branch("SMy_nS", &SMy_nS, "SMy_nS/D"); //min. # sigma away from left & right margins
  double SMx_nS_p;        Tout->Branch("SMx_nS_p", &SMx_nS_p, "SMx_nS_p/D"); //min. # sigma away from top & bot. margins 
  bool ARCut;             Tout->Branch("ARCut", &ARCut, "ARCut/O");
  bool fiduCut;           Tout->Branch("fiduCut", &fiduCut, "fiduCut/O");
  bool coinTADCCut;       Tout->Branch("coinTADCCut", &coinTADCCut, "coinTADCCut/O"); //HCAL/SH ADC coin time cut
  double coinTADC_nS;     Tout->Branch("coinTADC_nS", &coinTADC_nS, "coinTADC_nS/D"); //# sigma away from HCAL-SH good coinT zone
  //run info
  UInt_t T_rnum;          Tout->Branch("rnum", &T_rnum, "rnum/i");
  UInt_t T_segnum;        Tout->Branch("segnum", &T_segnum, "segnum/i");
  ULong64_t T_gevnum;     Tout->Branch("gevnum", &T_gevnum, "gevnum/l");
  double T_ebeam;         Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  double T_ebeam_corr;    Tout->Branch("ebeam_corr", &T_ebeam_corr, "ebeam_corr/D");
  double T_ebeam_std;     Tout->Branch("ebeam_std", &T_ebeam_std, "ebeam_std/D");
  //bcm/scaler
  //UInt_t T_segnumS;     if (get_scaler_info) Tout->Branch("segnumS", &T_segnumS, "segnumS/i");
  double T_dnewcnt;       if (get_scaler_info) Tout->Branch("dnewcnt", &T_dnewcnt, "dnewcnt/D"); 
  double T_dnewcurr;      if (get_scaler_info) Tout->Branch("dnewcurr", &T_dnewcurr, "dnewcurr/D");
  //per run info
  double T_dnewcharge;    Tout->Branch("dnewcharge", &T_dnewcharge, "dnewcharge/D");
  //double T_daqlvtm;       Tout->Branch("daqlvtm", &T_daqlvtm, "daqlvtm/D");
  double T_weight_norm;   Tout->Branch("weight_norm", &T_weight_norm, "weight_norm/D"); // = 1/charge/daqlvtm
  //kine
  double T_nu;            Tout->Branch("nu", &T_nu, "nu/D");
  double T_Q2;            Tout->Branch("Q2", &T_Q2, "Q2/D");
  double T_W2;            Tout->Branch("W2", &T_W2, "W2/D");
  double T_W;             Tout->Branch("W", &T_W, "W/D");
  double T_dpel;          Tout->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi;          Tout->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta;        Tout->Branch("etheta", &T_etheta, "etheta/D");
  double T_pelas;         Tout->Branch("pelas", &T_pelas, "pelas/D");
  double T_ethbend;       Tout->Branch("ethbend", &T_ethbend, "ethbend/D");
  double T_pN_exp;        Tout->Branch("pN_exp", &T_pN_exp, "pN_exp/D"); //exp. nucleon momentum
  double T_thN_exp;       Tout->Branch("thN_exp", &T_thN_exp, "thN_exp/D"); //exp. nucelon theta
  double T_epsilon;       Tout->Branch("epsilon", &T_epsilon, "epsilon/D"); // calculated using general eqn.
  double T_epsilon_p;     Tout->Branch("epsilon_p", &T_epsilon_p, "epsilon_p/D");
  double T_thpq_p;        Tout->Branch("thpq_p", &T_thpq_p, "thpq_p/D");
  double T_thpq_n;        Tout->Branch("thpq_n", &T_thpq_n, "thpq_n/D"); //n=>no deflection
  //track
  double T_vz;            Tout->Branch("vz", &T_vz, "vz/D");
  double T_trP;           Tout->Branch("trP", &T_trP, "trP/D");
  double T_trP_corr;      Tout->Branch("trP_corr", &T_trP_corr, "trP_corr/D");
  double T_trX;           Tout->Branch("trX", &T_trX, "trX/D");
  double T_trY;           Tout->Branch("trY", &T_trY, "trY/D");
  double T_trTh;          Tout->Branch("trTh", &T_trTh, "trTh/D");
  double T_trPh;          Tout->Branch("trPh", &T_trPh, "trPh/D");
  double T_tgX;           Tout->Branch("tgX", &T_tgX, "tgX/D");
  double T_tgY;           Tout->Branch("tgY", &T_tgY, "tgY/D");
  double T_tgTh;          Tout->Branch("tgTh", &T_tgTh, "tgTh/D");
  double T_tgPh;          Tout->Branch("tgPh", &T_tgPh, "tgPh/D");
  double T_fpX;           Tout->Branch("fpX", &T_fpX, "fpX/D");
  double T_fpY;           Tout->Branch("fpY", &T_fpY, "fpY/D");
  double T_fpTh;          Tout->Branch("fpTh", &T_fpTh, "fpTh/D");
  double T_fpPh;          Tout->Branch("fpPh", &T_fpPh, "fpPh/D");
  //BBCAL
  double T_ePS;           Tout->Branch("ePS", &T_ePS, "ePS/D"); 
  double T_rblkPS;        Tout->Branch("rblkPS", &T_rblkPS, "rblkPS/D"); 
  double T_cblkPS;        Tout->Branch("cblkPS", &T_cblkPS, "cblkPS/D"); 
  double T_idblkPS;       Tout->Branch("idblkPS", &T_idblkPS, "idblkPS/D"); 
  double T_atimePS;       Tout->Branch("atimePS", &T_atimePS, "atimePS/D"); 
  double T_eSH;           Tout->Branch("eSH", &T_eSH, "eSH/D"); 
  double T_xSH;           Tout->Branch("xSH", &T_xSH, "xSH/D"); 
  double T_ySH;           Tout->Branch("ySH", &T_ySH, "ySH/D"); 
  double T_rblkSH;        Tout->Branch("rblkSH", &T_rblkSH, "rblkSH/D"); 
  double T_cblkSH;        Tout->Branch("cblkSH", &T_cblkSH, "cblkSH/D"); 
  double T_idblkSH;       Tout->Branch("idblkSH", &T_idblkSH, "idblkSH/D"); 
  double T_atimeSH;       Tout->Branch("atimeSH", &T_atimeSH, "atimeSH/D"); 
  double T_EovP;          Tout->Branch("EovP", &T_EovP, "EovP/D"); 
  //HCAL (HE HCAL cluster in time)
  double T_eHCAL;         Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  double T_xHCAL;         Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  double T_yHCAL;         Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
  double T_indexHCAL;     Tout->Branch("indexHCAL", &T_indexHCAL, "indexHCAL/D"); 
  double T_idblkHCAL;     Tout->Branch("idblkHCAL", &T_idblkHCAL, "idblkHCAL/D"); 
  double T_rblkHCAL;      Tout->Branch("rblkHCAL", &T_rblkHCAL, "rblkHCAL/D"); 
  double T_cblkHCAL ;     Tout->Branch("cblkHCAL", &T_cblkHCAL, "cblkHCAL/D"); 
  double T_atimeHCAL;     Tout->Branch("atimeHCAL", &T_atimeHCAL, "atimeHCAL/D"); 
  double T_tdcHCAL;       Tout->Branch("tdcHCAL", &T_tdcHCAL, "tdcHCAL/D"); 
  double T_xHCAL_exp;     Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); 
  double T_yHCAL_exp;     Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D");
  double T_xHCAL_exp_p;   Tout->Branch("xHCAL_exp_p", &T_xHCAL_exp_p, "xHCAL_exp_p/D"); // average sbs_kick included 
  double T_dx;            Tout->Branch("dx", &T_dx, "dx/D"); 
  double T_dy;            Tout->Branch("dy", &T_dy, "dy/D");
  double T_p_def;         Tout->Branch("p_def", &T_p_def, "p_def/D"); // expected proton deflection
  double T_ToF;           Tout->Branch("ToF", &T_ToF, "ToF/D");
  //HCAL (All clusters)
  int T_idblkHCAL_aclN;   if (hcal_acl_ON) Tout->Branch("idblkHCAL_aclN", &T_idblkHCAL_aclN, "idblkHCAL_aclN/I"); 
  int T_idclHCAL_htote;   if (hcal_acl_ON) Tout->Branch("idclHCAL_htote", &T_idclHCAL_htote, "idclHCAL_htote/I"); //stores HCAL cl. index with highest total energy
  int T_idclHCAL_intime;  if (hcal_acl_ON) Tout->Branch("idclHCAL_intime", &T_idclHCAL_intime, "idclHCAL_intime/I"); //stores HCAL cl. that are in BBCAL/HCAL ADC coin time
  int T_idclHCAL_sthpq_p; if (hcal_acl_ON) Tout->Branch("idclHCAL_sthpq_p", &T_idclHCAL_sthpq_p, "idclHCAL_sthpq_p/I"); //stores HCAL cl. index with smallest thpq value
  double T_idblkHCAL_acl[maxNHCALcl];    if (hcal_acl_ON) Tout->Branch("idblkHCAL_acl", &T_idblkHCAL_acl, "idblkHCAL_acl[idblkHCAL_aclN]/D"); 
  double T_nblkHCAL_acl[maxNHCALcl];     if (hcal_acl_ON) Tout->Branch("nblkHCAL_acl", &T_nblkHCAL_acl, "nblkHCAL_acl[idblkHCAL_aclN]/D"); 
  double T_eblkHCAL_acl[maxNHCALcl];     if (hcal_acl_ON) Tout->Branch("eblkHCAL_acl", &T_eblkHCAL_acl, "eblkHCAL_acl[idblkHCAL_aclN]/D"); 
  double T_atimeblkHCAL_acl[maxNHCALcl]; if (hcal_acl_ON) Tout->Branch("atimeblkHCAL_acl", &T_atimeblkHCAL_acl, "atimeblkHCAL_acl[idblkHCAL_aclN]/D"); 
  double T_tdcblkHCAL_acl[maxNHCALcl];   if (hcal_acl_ON) Tout->Branch("tdcblkHCAL_acl", &T_tdcblkHCAL_acl, "tdcblkHCAL_acl[idblkHCAL_aclN]/D"); 
  double T_eHCAL_acl[maxNHCALcl];        if (hcal_acl_ON) Tout->Branch("eHCAL_acl", &T_eHCAL_acl, "eHCAL_acl[idblkHCAL_aclN]/D"); 
  double T_xHCAL_acl[maxNHCALcl];        if (hcal_acl_ON) Tout->Branch("xHCAL_acl", &T_xHCAL_acl, "xHCAL_acl[idblkHCAL_aclN]/D"); 
  double T_yHCAL_acl[maxNHCALcl];        if (hcal_acl_ON) Tout->Branch("yHCAL_acl", &T_yHCAL_acl, "yHCAL_acl[idblkHCAL_aclN]/D");
  //HODO
  int T_ncltmeanHODO;     Tout->Branch("ncltmeanHODO", &T_ncltmeanHODO, "ncltmeanHODO/I"); 
  double T_cltmeanHODO;   Tout->Branch("cltmeanHODO", &T_cltmeanHODO, "cltmeanHODO/D"); 
  //GEM
  double T_nhitsGEM;      Tout->Branch("nhitsGEM", &T_nhitsGEM, "nhitsGEM/D");
  double T_ngoodhitsGEM;  Tout->Branch("ngoodhitsGEM", &T_ngoodhitsGEM, "ngoodhitsGEM/D");
  double T_trchi2ndf;     Tout->Branch("trchi2ndf", &T_trchi2ndf, "trchi2ndf/D");
  //various coin time (trigger & ADC)
  double T_bbT_trig;      Tout->Branch("bbT_trig", &T_bbT_trig, "bbT_trig/D");
  double T_coinT_trig;    Tout->Branch("coinT_trig", &T_coinT_trig, "coinT_trig/D");
  double T_coinT_ADC_c;   Tout->Branch("coinT_ADC_c", &T_coinT_ADC_c, "coinT_ADC_c/D"); // centered w.r.t the mean
  //GRINCH
  double T_clsizeGRINCH;    if (conf>7) Tout->Branch("clsizeGRINCH", &T_clsizeGRINCH, "clsizeGRINCH/D");
  double T_cltmeanGRINCH;   if (conf>7) Tout->Branch("cltmeanGRINCH", &T_cltmeanGRINCH, "cltmeanGRINCH/D");
  double T_cltotmeanGRINCH; if (conf>7) Tout->Branch("cltotmeanGRINCH", &T_cltotmeanGRINCH, "cltotmeanGRINCH/D");
  double T_cltrindexGRINCH; if (conf>7) Tout->Branch("cltrindexGRINCH", &T_cltrindexGRINCH, "cltrindexGRINCH/D");

  // reading W cut limits
  std::vector<double> W_cutR; jmgr->GetVectorFromSubKey<double>(key,"W_cutR",W_cutR);
  // reading BB fiducial cut limits
  std::vector<double> bbfidu_cutR; jmgr->GetVectorFromSubKey<double>(key,"bbfidu_cutR",bbfidu_cutR);
  // reading HCAL/SH ADC coincidence time cut limits
  std::vector<double> coinTADC_cutR; jmgr->GetVectorFromSubKey<double>(key,"coinT_ADC_cutR",coinTADC_cutR);

  // HCAL cut definitions
  vector<double> dx_p_cut; jmgr->GetVectorFromSubKey<double>(key,"dx_p_cut", dx_p_cut);
  double sbs_kick = abs(dx_p_cut[0]);
  vector<double> dy_p_cut; jmgr->GetVectorFromSubKey<double>(key,"dy_p_cut", dy_p_cut);
  vector<double> hcal_active_area = cut::hcal_active_area_data(1,1,pass); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p_cut[1],dx_p_cut[1],dy_p_cut[1],hcal_active_area);
  // varying safety margin width by +/- 10% in vertical direction
  vector<double> hcal_safety_margin_p10p = cut::hcal_safety_margin(dx_p_cut[1]*1.1,dx_p_cut[1]*1.1,dy_p_cut[1],hcal_active_area);
  vector<double> hcal_safety_margin_m10p = cut::hcal_safety_margin(dx_p_cut[1]*.9,dx_p_cut[1]*.9,dy_p_cut[1],hcal_active_area);
  TH2F *h2_xyHCAL_p = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_p",sbs_kick,pass);
  TH2F *h2_xyHCAL_p_nfc = util_pd::TH2FHCALface_xy_data("h2_xyHCAL_p_nfc",sbs_kick,pass);

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = jmgr->GetValueFromSubKey<double>(key,"hcal_voffset");
  double hcal_hoffset = jmgr->GetValueFromSubKey<double>(key,"hcal_hoffset");
  double hcal_zoffset = jmgr->GetValueFromSubKey<double>(key,"hcal_zoffset");
  //vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(),HCAL_axes);
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetHCALtheta_rad(),HCAL_axes);
  TVector3 HCAL_origin = (sbsconf.GetHCALdist()+hcal_zoffset)*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  std::vector<double> ElossInTgt; // array to hold energy loss correction values per event
  long nevent=0, nevents=C->GetEntries(), neventsS=S->GetEntries(), index=0, tgevnumS, ngoodevs = 0; 
  int treenum=0, currenttreenum=0; UInt_t runnum=0, nseg, tsegnumS;
  double ebeam=sbsconf.GetEbeam(), ebeam_std=0., dnewcharge, daqlvtm; 
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
	  throw std::runtime_error("S tree index out of bounds! INVESTIGATE!");
	}
	S->GetEntry(index); index++;
	tsegnumS = segnumS; tgevnumS = gevnumS;
	tdnewcnt = dnewcnt; tdnewcurr = dnewcurr;
      }
      /* finding nearest scaler event for the current T event */
      while (gevnum>gevnumS && rnumS==rnum) {
	if (index==neventsS) {
	  std::cout << Form("Run %u | GevNum %u | GevNumS %llu",rnum,gevnum,gevnumS) << std::endl;
	  throw std::runtime_error("S tree index out of bounds! INVESTIGATE!");
	}
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

      // code block to read run specific info
      // ebeam, charge, daqlvtm are read
      if (nevent == 1 || rnum != runnum) {
	runnum = rnum; nseg=1;
	auto it = std::find_if(crun.begin(), crun.end(), [=](CodaRun const& cr) {return cr.runnum == runnum;});
	if (it != crun.end()) {
	  ebeam = it->ebeam; ebeam_std = it->ebeam_std;
	  dnewcharge = it->charge; daqlvtm = it->daqlvtm;
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
    T_coinT_trig = coin_time; h_coinT_trig->Fill(coin_time);

    // constructing the 4 vectors
    /* Reaction    : e + e' -> p + p'
       Conservation: Pe + Peprime = Pp + Ppprime */
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
    double Q2recon{0},W2recon{0},epsilon{0},epsilon_p{0};
    if (model == 0) { // p as independent variable
      nu = Pe.E() - Peprime.E();
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((Pe.E() - Peprime.Pz()) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      PNprime.SetPxPyPzE(pN_expect*pNhat.X(), pN_expect*pNhat.Y(), pN_expect*pNhat.Z(), nu+PN.E());
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), thelas);
      epsilon = kine::epsilon_general(thelas,Q2recon,nu);
      epsilon_p = kine::epsilon(thelas,Q2recon,"p");
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
      //W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
      W2recon = kine::W2_general(Pe.E(), Peprime.E(), etheta, Ntype);
    } else if (model == 2) { // 4-vector calculation
      nu = q.E();
      PNprime = q + PN; 
      pNhat = PNprime.Vect().Unit();
      Q2recon = -q.M2();
      epsilon = kine::epsilon_general(etheta,Q2recon,nu);
      epsilon_p = kine::epsilon(etheta,Q2recon,"p");
      W2recon = PNprime.M2();
      //--
      pN_expect = kine::pN_expect(nu, Ntype);
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
    T_pN_exp = pN_expect;
    T_thN_exp = thetaN_expect;
    T_epsilon = epsilon;
    T_epsilon_p = epsilon_p;

    // defining W cut
    WCut = T_W >= W_cutR[0] && T_W <= W_cutR[1];

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
    // per run info
    T_dnewcharge = dnewcharge;
    //T_daqlvtm = daqlvtm;
    T_weight_norm = 1./dnewcharge/daqlvtm;

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

    T_fpX = xfp[0];
    T_fpY = yfp[0];
    T_fpTh = thfp[0];
    T_fpPh = phfp[0];

    // calculating bend angle
    T_ethbend = kine::ethbend(T_tgTh,T_tgPh,T_fpTh,T_fpPh);
    
    // defining BB fiducial cut
    bbfiduCut = abs(T_fpX - 0.9*T_fpTh - bbfidu_cutR[0]) <= bbfidu_cutR[1];

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

    T_EovP = (eSH+ePS)/p[0];

    T_ncltmeanHODO = ncltmeanHODO;
    T_cltmeanHODO = cltmeanHODO[0];

    T_nhitsGEM = nhitsGEM[0];
    T_ngoodhitsGEM = ngoodhitsGEM[0];
    T_trchi2ndf = trchi2ndf[0];

    T_clsizeGRINCH = clsizeGRINCH;
    T_cltmeanGRINCH = cltmeanGRINCH;
    T_cltotmeanGRINCH = cltotmeanGRINCH;
    T_cltrindexGRINCH = cltrindexGRINCH;

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_xHCAL_exp_p = xyHCAL_exp[0]-sbs_kick;

    // Calculating proton deflection angle
    double BdL = (sbsmag / 100.) * expconst::sbsmaxfield_data(conf) * expconst::sbsdipolegap;
    double proton_thetabend = 0.3 * BdL / PNprime.Vect().Mag();  // p*theta = 0.3*BdL
    double proton_deflection = tan(proton_thetabend)*(sbsconf.GetHCALdist()+hcal_zoffset-(sbsconf.GetSBSdist()+expconst::sbsdipolegap/2.0));
    T_p_def = proton_deflection;

    /*
      Implementing HCAL offline clustering algorithm
    */
    int HtotE_idcl = -1;    //HCAL cl. index that highest total energy
    int inTime_idcl = -1;   //HCAL cl. index that is within HCAL/SH ADC coin. time
    int sthpq_p_idcl = -1;  //HCAL cl. index with smallest thpq_p value 
    int highest_idcl = 0;   //Highest HCAL cl. index among all 4 cl. algos
    double temp_thpq = 1e6;
    if (idblkHCAL_aclN>0) { //if HCAL has a cluster
      bool index_mismatch = false; //bool to check HtotE cluster index mismatch
      if (hcal_acl_ON) {
	// Clusters in the *.clus.* variables aren't sorted in terms of their total energy.
	// Instead they are sorted by the energy of the HE blocks in the corresponding cluster.
	// But for offline clustering we want them to be sorted by total energy instead. So, 
	// let's sort the array indices in terms of total cluster energy.
	std::vector<size_t> sortedIndices = util_pd::SortIndices(eHCAL_acl, idblkHCAL_aclN);
	HtotE_idcl = sortedIndices[0]; if (HtotE_idcl>highest_idcl) highest_idcl = HtotE_idcl;
	// Sanity check: check if the sorted HE (tot) matches with indexHCAL
	if (HtotE_idcl!=(int)indexHCAL) {
	  index_mismatch = true;
	  //std::cout << Form("[*INDEX mismatch*] HCAL HtotE cl: rnum %u, segnum %u, gev %u\n",rnum,nseg-1,gevnum);
	  //std::cout << Form(" exp index %d, obs index %d, exp eng %f, obs eng %f\n------\n",(int)indexHCAL,HtotE_idcl,eHCAL,eHCAL_acl[HtotE_idcl]);
	}
	// if index_mismatch and the default cl is in time then accept it
	if (index_mismatch && abs(atimeHCAL-atimeSH-coinTADC_cutR[0])<=coinTADC_cutR[2]*coinTADC_cutR[1]) {
	  inTime_idcl = (int)indexHCAL; sthpq_p_idcl = (int)indexHCAL;	    
	} else {
	  // Looping through ALL HCAL clusters
	  for (int ihcl=0; ihcl<idblkHCAL_aclN; ihcl++) {
	    // getting cl index sorted by the total energy
	    int ihcl_sorted = sortedIndices[ihcl];
      
	    // picking the HE cluster that is in time (i.e. HCAL/SH ADC coin time)
	    bool coinT_cut = abs(atimeblkHCAL_acl[ihcl_sorted]-atimeSH-coinTADC_cutR[0])<=coinTADC_cutR[2]*coinTADC_cutR[1];
	    if (inTime_idcl==-1 && coinT_cut) {
	      inTime_idcl = ihcl_sorted;
	      if (inTime_idcl>highest_idcl) highest_idcl = inTime_idcl;
	    }

	    // // filling HCAL cl vectors
	    // // ** info on the HE block on each cl
	    // T_idblkHCAL_acl[ihcl_sorted] = idblkHCAL_acl[ihcl_sorted];
	    // T_nblkHCAL_acl[ihcl_sorted] = nblkHCAL_acl[ihcl_sorted];
	    // T_eblkHCAL_acl[ihcl_sorted] = eblkHCAL_acl[ihcl_sorted];
	    // T_atimeblkHCAL_acl[ihcl_sorted] = atimeblkHCAL_acl[ihcl_sorted];
	    // T_tdcblkHCAL_acl[ihcl_sorted] = tdcblkHCAL_acl[ihcl_sorted];
	    // // ** info on the params of the cl itself
	    // T_eHCAL_acl[ihcl_sorted] = eHCAL_acl[ihcl_sorted];
	    // T_xHCAL_acl[ihcl_sorted] = xHCAL_acl[ihcl_sorted];
	    // T_yHCAL_acl[ihcl_sorted] = yHCAL_acl[ihcl_sorted];

	    if (inTime_idcl!=-1) h_hcl_inTime_idcl->Fill(inTime_idcl);
	    if (WCut) {
	      if (inTime_idcl!=-1) h_hcl_inTime_idcl_WCut->Fill(inTime_idcl);
      
	      // HE block related variables
	      h2_hclHE_eng_vs_idcl->Fill(ihcl_sorted,eblkHCAL_acl[ihcl_sorted]);
	      h2_hclHE_atime_vs_idcl->Fill(ihcl_sorted,atimeblkHCAL_acl[ihcl_sorted]);
	      h2_hclHE_tdc_vs_idcl->Fill(ihcl_sorted,tdcblkHCAL_acl[ihcl_sorted]);
      
	      // Cluster variables
	      h2_hcl_eng_vs_idcl->Fill(ihcl_sorted,eHCAL_acl[ihcl_sorted]);
	      h2_hcl_nblk_vs_idcl->Fill(ihcl_sorted,nblkHCAL_acl[ihcl_sorted]);
	    }

	    // picking the cluster that has smallest thpq and passes coinT_cut and sFrac_cut
	    bool sFrac_cut = eHCAL_acl[ihcl_sorted]/(ebeam_corr-trP_corr)>hcal_sF_cutR;
	    if (coinT_cut && sFrac_cut) {
	      // Calculating thpq (both w & w/o deflection due to SBS dipole)
	      // assuming no deflection (using "n" for no deflection)
	      TVector3 HCAL_pos = HCAL_origin + xHCAL_acl[ihcl_sorted]*HCAL_axes[0] + yHCAL_acl[ihcl_sorted]*HCAL_axes[1];
	      // TVector3 n_dir = (HCAL_pos - vertex).Unit();
	      // double thpq_n = acos(n_dir.Dot(pNhat));
	      // p
	      TVector3 p_dir = (HCAL_pos + proton_deflection*HCAL_axes[0] - vertex);
	      double thpq_p = acos(p_dir.Unit().Dot(pNhat));
	      // finding the cl. id. with smallest thpq_p value
	      if (thpq_p < temp_thpq) sthpq_p_idcl = ihcl_sorted;
	      temp_thpq = thpq_p;
	      // if (min(thpq_p,thpq_n) < temp_thpq) sthpq_p_idcl = ihcl_sorted;
	      // temp_thpq = min(thpq_p,thpq_n);
	    }
	  }
	}
	if (inTime_idcl==-1) inTime_idcl = HtotE_idcl; // couldn't find any cl. in time, switching to HtotE cl
	if (sthpq_p_idcl==-1) sthpq_p_idcl = inTime_idcl; // couldn't find any cl. w/ sthqp, switching to cl. in time

	// -- filling HCAL cl vectors
	if (!index_mismatch && sthpq_p_idcl>highest_idcl) highest_idcl = sthpq_p_idcl; // last step of determining the highest idcl
	for (int ihcl=0; ihcl<=highest_idcl; ihcl++) {
	  // ** info on the HE block on each cl
	  T_idblkHCAL_acl[ihcl] = idblkHCAL_acl[ihcl];
	  T_nblkHCAL_acl[ihcl] = nblkHCAL_acl[ihcl];
	  T_eblkHCAL_acl[ihcl] = eblkHCAL_acl[ihcl];
	  T_atimeblkHCAL_acl[ihcl] = atimeblkHCAL_acl[ihcl];
	  T_tdcblkHCAL_acl[ihcl] = tdcblkHCAL_acl[ihcl];
	  // ** info on the params of the cl itself
	  T_eHCAL_acl[ihcl] = eHCAL_acl[ihcl];
	  T_xHCAL_acl[ihcl] = xHCAL_acl[ihcl];
	  T_yHCAL_acl[ihcl] = yHCAL_acl[ihcl];
	}
	if (index_mismatch && inTime_idcl==(int)indexHCAL) {
	  // if there is index mismatch and the default cl is in time
	  T_dx = xHCAL - xyHCAL_exp[0];
	  T_dy = yHCAL - xyHCAL_exp[1];
	  T_eHCAL = eHCAL;
	  T_xHCAL = xHCAL;
	  T_yHCAL = yHCAL;
	  T_indexHCAL = indexHCAL;
	  T_rblkHCAL = rblkHCAL;
	  T_cblkHCAL = cblkHCAL;
	  T_idblkHCAL = idblkHCAL;
	  T_atimeHCAL = atimeHCAL;
	  T_tdcHCAL = tdcHCAL[0];
	  T_idblkHCAL_aclN = idblkHCAL_aclN;
	  T_idclHCAL_htote = HtotE_idcl;
	  T_idclHCAL_intime = inTime_idcl;
	  T_idclHCAL_sthpq_p = sthpq_p_idcl;
	} else {
	  // fill basic HCAL variables using inTime clusters
	  T_dx = xHCAL_acl[inTime_idcl] - xyHCAL_exp[0];
	  T_dy = yHCAL_acl[inTime_idcl] - xyHCAL_exp[1];
	  T_eHCAL = eHCAL_acl[inTime_idcl];
	  T_xHCAL = xHCAL_acl[inTime_idcl];
	  T_yHCAL = yHCAL_acl[inTime_idcl];
	  T_indexHCAL = indexHCAL;
	  T_rblkHCAL = rblkHCAL_acl[inTime_idcl];
	  T_cblkHCAL = cblkHCAL_acl[inTime_idcl];
	  T_idblkHCAL = idblkHCAL_acl[inTime_idcl];
	  T_atimeHCAL = atimeblkHCAL_acl[inTime_idcl];
	  T_tdcHCAL = tdcblkHCAL_acl[inTime_idcl];
	  T_idblkHCAL_aclN = idblkHCAL_aclN;
	  T_idclHCAL_htote = HtotE_idcl;
	  T_idclHCAL_intime = inTime_idcl;
	  T_idclHCAL_sthpq_p = sthpq_p_idcl;
	}
      } else { //if HCAL has a cl but secondary clustering is off!
	T_dx = xHCAL - xyHCAL_exp[0];
	T_dy = yHCAL - xyHCAL_exp[1];
	T_eHCAL = eHCAL;
	T_xHCAL = xHCAL;
	T_yHCAL = yHCAL;
	T_indexHCAL = indexHCAL;
	T_rblkHCAL = rblkHCAL;
	T_cblkHCAL = cblkHCAL;
	T_idblkHCAL = idblkHCAL;
	T_atimeHCAL = atimeHCAL;
	T_tdcHCAL = tdcHCAL[0];
      }
    } else { //if HCAL doesn't have a cl
      T_dx = -99;
      T_dy = -99;
      T_eHCAL = -99;
      T_xHCAL = -99;
      T_yHCAL = -99;
      T_indexHCAL = -99;
      T_rblkHCAL = -99;
      T_cblkHCAL = -99;
      T_idblkHCAL = -99;
      T_atimeHCAL = -99;
      T_tdcHCAL = -99;
      T_idblkHCAL_aclN = -99;
      T_idclHCAL_htote = -99;
      T_idclHCAL_intime = -99;
      T_idclHCAL_sthpq_p = -99;
    }

    double dx = T_dx;
    double dy = T_dy;

    // HCAL/SH SDC coincidence time
    double coinTADC_corr_pos = T_atimeHCAL - atimeSH - coinTADC_cutR[0];
    coinTADCCut = abs(coinTADC_corr_pos) <= coinTADC_cutR[2]*coinTADC_cutR[1];
    coinTADC_nS = abs(coinTADC_corr_pos)/coinTADC_cutR[1];
    T_coinT_ADC_c = coinTADC_corr_pos;

    // Calculating thpq (both w & w/o deflection due to SBS dipole)
    // assuming no deflection (using "n" for no deflection)
    TVector3 HCAL_pos = HCAL_origin + T_xHCAL*HCAL_axes[0] + T_yHCAL*HCAL_axes[1];
    TVector3 n_dir = (HCAL_pos - vertex).Unit();
    T_thpq_n = acos(n_dir.Dot(pNhat));
    // p 
    TVector3 p_dir = (HCAL_pos + proton_deflection*HCAL_axes[0] - vertex);
    T_thpq_p = acos(p_dir.Unit().Dot(pNhat));
    // Calculate ToF
    double ToF = (p_dir.Mag() / constant::c) * sqrt(1. + pow((constant::Mp/PNprime.Vect().Mag()), 2));
    T_ToF = ToF*1e9; //ns 

    // HCAL active area and safety margin cuts [Fiducial region]
    ARCut = cut::inHCAL_activeA(T_xHCAL,T_yHCAL,hcal_active_area);
    SMCut = cut::inHCAL_safety_margin(target,xyHCAL_exp[0],xyHCAL_exp[1],sbs_kick,hcal_safety_margin);
    fiduCut = ARCut && SMCut; 
    // *** Implementing Sebastian Seed's idea for flexible SMCut
    // First step is to evalute current xyHCAL_exp values are how many sigmas away from the active area boundaries
    double SMy_nS_l = fabs(T_yHCAL_exp-hcal_active_area[3])/dy_p_cut[1];
    double SMy_nS_r = fabs(T_yHCAL_exp-hcal_active_area[2])/dy_p_cut[1]; 
    double SMx_nS_t_p = fabs(T_xHCAL_exp-sbs_kick-hcal_active_area[0])/dx_p_cut[1];
    double SMx_nS_b_p = fabs(T_xHCAL_exp-sbs_kick-hcal_active_area[1])/dx_p_cut[1];
    // The above calculation is symmetric for the events landing on either sides of AR boundaries. Hence, it is
    // necessary to identify the events landing outside of AR to be able to exclude them from Cut. The easiest I can
    // think of right now is to define a bool that turns 1 if xyHCAL_exp values are within AR
    bool p_xyexp_in_AR = cut::inHCAL_activeA(T_xHCAL_exp-sbs_kick,T_yHCAL_exp,hcal_active_area);
    // Now, its time to define just one variable to be able to cut on both sides of boundaries for a given direction.
    // The way I can think of is to cut on the minimum # sigmas for a given direction to accomplish this. At the same
    // time be sure to assign garbage value for events with xyHCAL_exp values landing outside of AR
    SMy_nS = p_xyexp_in_AR ? min(SMy_nS_l,SMy_nS_r) : -99;
    SMx_nS_p = p_xyexp_in_AR ? min(SMx_nS_t_p,SMx_nS_b_p) : -99;
    // defining HCAL cuts
    pCut = cut::SpotCut(dx,dx_p_cut[0],dx_p_cut[1],dx_p_cut[2],dy,dy_p_cut[0],dy_p_cut[1],dy_p_cut[2]);
    // ***
    pdx_nS = fabs(dx-dx_p_cut[0])/dx_p_cut[1];
    dy_nS = fabs(dy-dy_p_cut[0])/dy_p_cut[1];

    // W cut
    if (WCut&&idblkHCAL_aclN!=0&&bbfiduCut) {
      if (T_eHCAL>0) {
	if (abs(dy)<0.3) h_dxHCAL_nfc->Fill(dx);
	h_dyHCAL_nfc->Fill(dy);
      }
      // fiducial cut (Should we use it for LH2 data?)
      if (fiduCut) {
	if (T_eHCAL>0) {
	  if (abs(dy)<0.3) {
	    h_dxHCAL->Fill(dx);
	    h_dx_w_p_def->Fill(dx+T_p_def);
	  }
	  if (pCut) {
	    h_coinT_ADC->Fill(T_atimeHCAL - atimeSH);
	    h_coinT_ADC_c->Fill(T_coinT_ADC_c);
	  }
	  h_dyHCAL->Fill(dy);
	  h2_rcHCAL->Fill(T_cblkHCAL, T_rblkHCAL);
	  h2_dxdyHCAL->Fill(dy, dx);
	}

	if (pCut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
      }
      if (pCut) h2_xyHCAL_p_nfc->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
    }

    // fiducial cut but no W cut
    if (T_W>0) {      
      if (fiduCut) {
	h_W->Fill(T_W);
	if (pCut) { 
	  h_W_cut->Fill(T_W);
	  h_W2_cut->Fill(T_W2);
	} else {
	  h_W_acut->Fill(T_W);
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
  c1->cd(2); //
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");
  c1->cd(3); //
  h2_xyHCAL_p_nfc->Draw("colz");
  util_pd::PlotFiduCut(pass,hcal_active_area,hcal_safety_margin);
  c1->cd(4); //
  h2_xyHCAL_p->Draw("colz");
  util_pd::PlotFiduCut(pass,hcal_active_area,hcal_safety_margin);
  c1->SaveAs(Form("%s[",outPlot.Data())); c1->SaveAs(Form("%s",outPlot.Data())); c1->Write();
  //**** -- ***//

  /**** Canvas 2 (dx & dy) ****/
  TCanvas *c2 = util_pd::TC("c2",2,2);
  c2->SetGridx();
  std::vector<double> hdxp_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_p_fitR", hdxp_fitR);
  std::vector<double> hdy_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dyHCAL_fitR", hdy_fitR);
  gStyle->SetOptFit(1);
  c2->cd(1); //
  gPad->SetGridx();
  TF1 *fdxp = fit::fit_1gs_nbg(hdxp_fitR,h_dxHCAL);
  double dxpM = fdxp->GetParameter(1); double dxpS = fdxp->GetParameter(2);
  c2->cd(2); //
  gPad->SetGridx();
  TF1 *fdy = fit::fit_1gs_nbg(hdy_fitR,h_dyHCAL);
  double dypM = fdy->GetParameter(1); double dypS = fdy->GetParameter(2);
  c2->cd(3); //
  gPad->SetGridx();
  TF1 *fdxp_nfc = fit::fit_1gs_nbg(hdxp_fitR,h_dxHCAL_nfc);
  double dxpM_nfc = fdxp_nfc->GetParameter(1); double dxpS_nfc = fdxp_nfc->GetParameter(2);
  c2->cd(4); //
  gPad->SetGridx();
  TF1 *fdy_nfc = fit::fit_1gs_nbg(hdy_fitR,h_dyHCAL_nfc);
  double dypM_nfc = fdy->GetParameter(1); double dypS_nfc = fdy->GetParameter(2);
  c2->SaveAs(Form("%s",outPlot.Data())); c2->Write();
  //**** -- ***//

  /**** Canvas 3 (ADC coin time & p Deflection) ****/
  TCanvas *c3 = util_pd::TC("c3",2,2);
  c3->SetGridx();
  gStyle->SetOptFit(1);
  //
  c3->cd(1);
  gPad->SetGridx();
  std::vector<double> hct_fitR{-5,5,1.6,1.6};
  TF1 *fct = fit::fit_1gs_nbg(hct_fitR,h_coinT_ADC);
  double ctM = fct->GetParameter(1); double ctS = fct->GetParameter(2);
  h_coinT_ADC_c->Draw("same");
  h_coinT_ADC_c->SetLineColor(kBlack);
  //
  c3->cd(2);
  gPad->SetGridx();
  std::vector<double> hpdef_fitR{-0.5,0.5,1.2,1.2};
  TF1 *fpdef = fit::fit_1gs_nbg(hpdef_fitR,h_dx_w_p_def);
  double pdefM = fpdef->GetParameter(1); double pdefS = fpdef->GetParameter(2);
  h_dx_w_p_def->Draw("same");
  h_dx_w_p_def->SetLineColor(kBlack);
  //
  c3->cd(3);
  gPad->SetGridx();
  std::vector<double> hW2_fitR{0.6,1.2,1.2,1.2};
  TF1 *fW2 = fit::fit_1gs_nbg(hW2_fitR,h_W2_cut);
  double W2M = fW2->GetParameter(1); double W2S = fW2->GetParameter(2);
  h_W2_cut->Draw("same");
  h_W2_cut->SetLineColor(kBlack); 
  c3->SaveAs(Form("%s",outPlot.Data())); c3->Write();
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
  pt->AddText(Form(" HCAL offsets: v = %.4f, h = %.4f, z = %.4f",hcal_voffset,hcal_hoffset,hcal_zoffset));
  char const * hcl_algo_flag = hcal_acl_ON ? "In-time" : "Default";
  pt->AddText(Form(" Best HCAL cluster choice algorithm: %s",hcl_algo_flag));
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
  pt->AddText(Form(" Inbuilt BB fiducial cut: |fpX-0.9*fpTh-%.2f| #leq %.2f",bbfidu_cutR[0],bbfidu_cutR[1]));
  pt->AddText(Form(" Inbuilt coin. time cut: |atimeHCAL-atimeSH-%.1f| #leq %.1f*%.3f",coinTADC_cutR[0],coinTADC_cutR[2],coinTADC_cutR[1]));
  pt->AddText(" Fit info: ");
  pt->AddText(" p peak, w/ fiducial cut: dxpM,dxpS,dypM,dypS ");
  pt->AddText(Form(" %.5f,%.5f,%.5f,%.5f",dxpM,dxpS,dypM,dypS));
  pt->AddText(" p peak, w/o fiducial cut: dxpM_nfc,dxpS_nfc,dypM_nfc,dypS_nfc ");
  pt->AddText(Form(" %.5f,%.5f,%.5f,%.5f",dxpM_nfc,dxpS_nfc,dypM_nfc,dypS_nfc));
  pt->AddText(" coinT ADC cut: ctM,ctS ");
  pt->AddText(Form(" %.5f,%.5f",ctM,ctS));
  pt->AddText(" W2 distribution: W2M,W2S ");
  pt->AddText(Form(" %.5f,%.5f",W2M,W2S));    
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
  std::cout << "dxpM,dxpS,dypM,dypS,dxpM_nfc,dxpS_nfc,dypM_nfc,dypS_nfc" << "\n";
  std::cout << dxpM<<","<<dxpS<<","<<dypM<<","<<dypS<<","<<dxpM_nfc<<","<<dxpS_nfc<<","<<dypM_nfc<<","<<dypS_nfc<<"\n";
  std::cout << "------------------" << "\n";
  std::cout << "ctM,ctS\n";
  std::cout << ctM<<","<<ctS<<"\n";
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
  h2_xyHCAL_p->Write();
  h2_xyHCAL_p_nfc->Write();
  h_coinT_trig->Write();
  h_coinT_ADC->Write();
  h_dx_w_p_def->Write();
  if (hcal_acl_ON) {
    h_hcl_inTime_idcl->Write(); h_hcl_inTime_idcl_WCut->Write();
    h2_hclHE_eng_vs_idcl->Write(); h2_hclHE_atime_vs_idcl->Write(); h2_hclHE_tdc_vs_idcl->Write();
    h2_hcl_eng_vs_idcl->Write(); h2_hcl_nblk_vs_idcl->Write();
  }
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
** W_cutR : W cut range.
- W_cutR[0] : lower limit
- W_cutR[1] : upper limit
** coinT_ADC_cutR : HCAL/SH ADC coin time cut range.
- coinT_ADC_cutR[0] : Mean
- coinT_ADC_cutR[1] : Sigma
- coinT_ADC_cutR[2] : # sigmas to cut on
** SBS_field : yet to be added!!
** max_N_tracks : # tracks in an event
** max_N_HCAL_clusters : # HCAL clusters being stored per event
** hcal_acl_ON : flag to turn on HCAL best cl. choice algorithms
** hcal_samp_frac_cutR : lower limit of HCAL sampling fraction cut (only used if hcal_acl_ON is true)
** hcal_v(h)offset : vertical(horizontal) offsets of HCAL
** dx_p_cut : deltax p peak cut definitions. Needed to constitute p spot cuts. 
- dx_p_cut[0] => p peak position, 
- dx_p_cut[1] => p peak RMS, 
- dx_p_cut[2] => # sigma to include in the cut
** dy_p_cut : deltay p peak cut definitions. Needed to constitute p spot cuts.
- dy_p_cut[0] => p peak position, 
- dy_p_cut[1] => p peak RMS, 
- dy_p_cut[2] => # sigma to include in the cut
** h_dx(dy)HCAL_lims : h_dx(dy)HCAL histogram limits (Can be found in the output ROOT file)
- h_dx(dy)HCAL_lims[0] : No. of bins of h_dx(dy)HCAL histograms
- h_dx(dy)HCAL_lims[1] : xmin
- h_dx(dy)HCAL_lims[2] : xmax
** h_dxHCAL_p_fitR : Fit ranges for the proton(neutron) signal peak in deltax dist. (h_dxHCAL histogram)
Algorithm fits the distribution twice for optimization.
- h_dxHCAL_p_fitR[0] : xmin for 1st fit (crude). Try to avoid any secondary peak.
- h_dxHCAL_p_fitR[1] : xmax for 1st fit (crude). Try to avoid any secondary peak.
- h_dxHCAL_p_fitR[2] : # sigma below the peak for 2nd fit (fine)
- h_dxHCAL_p_fitR[3] : # sigma above the peak for 2nd fit (fine)
** h_dyHCAL_fitR : Fit ranges for the signal peak in deltay dist. (h_dxHCAL histogram)
Algorithm fits the distribution twice for optimization.
- h_dyHCAL_fitR[0] : xmin for 1st fit (crude). Try to avoid any secondary peak.
- h_dyHCAL_fitR[1] : xmax for 1st fit (crude). Try to avoid any secondary peak.
- h_dyHCAL_fitR[2] : # sigma below the peak for 2nd fit (fine)
- h_dyHCAL_fitR[3] : # sigma above the peak for 2nd fit (fine)
*/
