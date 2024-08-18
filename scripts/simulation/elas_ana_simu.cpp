/* 
   This macro will perform elastic (LH2) analysis for GMn using MC data.
   E.g. Config. File: sbs14-sbs70p-simu/conf_elas_ana_simu.json
   * A brief description of all the config. file parameters can be
   found at the bottom of this script.
   -----
   P. Datta  Created  11-05-2022 
*/

// TO-DO
// 1. Energy loss calculations

#include <vector>
#include <iostream>
#include <unordered_map>

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

#include "../../include/gmn_ana.h"
#include "../../dflay/src/JSONManager.cxx"

/* this script will only analyze LH2 data */
static const std::string target = "LH2";

int elas_ana_simu (const char *configfilename, 
		    int model=1, //Analysis model
		    std::string generator="simc",
		    std::string process="heep")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // clock to keep macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);
  std::string key = generator + "_" + process + "_model" + std::to_string(model);

  // setting verbosity
  int verbose = jmgr->GetValueFromSubKey<int>(key,"verbose");
  int verbosefn = jmgr->GetValueFromSubKey<int>(key,"verbose_function");

  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromSubKey<int>(key,"SBS_config");
  int sbsmag = jmgr->GetValueFromSubKey<int>(key,"SBS_magnet_percent");
  double sbsscalefield = jmgr->GetValueFromSubKey<double>(key,"SBS_scale_field");
  SBSconfig sbsconf(conf, sbsmag);
  std::cout << sbsconf;

  // reading job summary and parsing ROOT trees
  std::string rfd = jmgr->GetValueFromSubKey_str(key,"rootfile_dir");
  std::string prefix = jmgr->GetValueFromSubKey_str(key,"prefix_to_filebase");
  int njobs = jmgr->GetValueFromSubKey<int>(key,"Njobs_to_ana"); // # MC jobs to analyze
  //std::vector<SimuJob> sjobs; util_pd::ReadSimuJobSummary(rfd,prefix,conf,sbsmag,gen,target,njobs,verbosefn,sjobs);
  std::vector<SimuJob> sjobs; util_pd::ReadSimuJobSummary(rfd,prefix,conf,sbsmag,generator,process,njobs,verbosefn,sjobs);
  TChain *C = new TChain("T"); util_pd::LoadSimuROOTTree(sjobs,verbosefn,C);

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
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
  int maxNtr = jmgr->GetValueFromSubKey<int>(key,"max_N_tracks");
  C->SetBranchStatus("*",0);
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
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk","atimeblk","tdctimeblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL,&atimeHCAL,&tdcHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

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

  //MC variables
  double mc_sigma, mc_fnucl, mc_ebeam, mc_np, mc_veE, mc_vetheta;  // mc_sigma => Cross-section weight  std::vector<std::string> mc = {"mc_sigma","mc_fnucl"};  // Default: g4sbs gen.
  std::vector<std::string> mc = {"mc_sigma","mc_fnucl","mc_np"};   // Default: g4sbs gen.
  std::vector<void*> mc_mem = {&mc_sigma,&mc_fnucl}; 
  if (generator.compare("simc")==0) {
    mc = {"simc_Weight","simc_fnucl","simc_Ebeam","simc_veE","simc_vetheta"};  
    mc_mem = {&mc_sigma,&mc_fnucl,&mc_ebeam,&mc_veE,&mc_vetheta}; 
  }
  setrootvar::setbranch(C,"MC",mc,mc_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);

  // defining the outputfile
  std::string filebase = jmgr->GetValueFromSubKey_str(key,"outfile_prefix");
  filebase = filebase.empty() ? "" : filebase + "_";
  if (process.compare("inel")==0) filebase = filebase + "inel_";
  TString outFile = Form("siout/%selas_ana_%s_sbs%d_sbs%dp_model%d.root",filebase.c_str(),generator.c_str(),conf,sbsmag,model);
  TFile *fout = new TFile(outFile.Data(),"RECREATE");

  // defining histograms
  TH1F *h_W = util_pd::TH1FhW("h_W");
  TH1F *h_W_cut = util_pd::TH1FhW("h_W_cut");
  TH1F *h_W_acut = util_pd::TH1FhW("h_W_acut");
  TH1F *h_W2_cut = new TH1F("h_W2_cut","",200,-1,4);
  TH1F *h_dpel = new TH1F("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = util_pd::TH1FhQ2("h_Q2", conf);
  std::vector<double> hdx_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_lims", hdx_lim);
  std::vector<double> hdy_lim; jmgr->GetVectorFromSubKey<double>(key,"h_dyHCAL_lims", hdy_lim);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","W & fiducial cuts;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dxHCAL_nfc = new TH1F("h_dxHCAL_nfc","W cut;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","W & fiducial cuts;y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_dyHCAL_nfc = new TH1F("h_dyHCAL_nfc","W cut;y_{HCAL}^{obs} - y_{HCAL}^{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_dxHCAL_p = new TH1F("h_dxHCAL_p","mc_fnucl = p;x_{HCAL}^{obs} - x_{HCAL}^{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);

  TH1F *h_dx_w_p_def = new TH1F("h_dx_w_p_def",";dx+p_def (m);",200,-1,1);

  TH2F *h2_rcHCAL = util_pd::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = util_pd::TH2FdxdyHCAL("h2_dxdyHCAL");

  // Defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  Tout->SetMaxTreeSize(4000000000LL);
  //cuts
  bool WCut;              Tout->Branch("WCut", &WCut, "WCut/B");
  bool bbfiduCut;         Tout->Branch("bbfiduCut", &bbfiduCut, "bbfiduCut/O");
  bool pCut;              Tout->Branch("pCut", &pCut, "pCut/B");
  double pdx_nS;          Tout->Branch("pdx_nS", &pdx_nS, "pdx_nS/D"); //# sigma away from p dx peak
  double dy_nS;           Tout->Branch("dy_nS", &dy_nS, "dy_nS/D"); //# sigma away from dy peak 
  // --
  bool SMCut;             Tout->Branch("SMCut", &SMCut, "SMCut/B");
  double SMy_nS;          Tout->Branch("SMy_nS", &SMy_nS, "SMy_nS/D"); //min. # sigma away from left & right margins
  double SMx_nS_p;        Tout->Branch("SMx_nS_p", &SMx_nS_p, "SMx_nS_p/D"); //min. # sigma away from top & bot. margins   
  bool ARCut;             Tout->Branch("ARCut", &ARCut, "ARCut/B");
  bool fiduCut;           Tout->Branch("fiduCut", &fiduCut, "fiduCut/B");
  //MC related
  double weight;          Tout->Branch("weight", &weight, "weight/D");
  double weight_norm;     Tout->Branch("weight_norm", &weight_norm, "weight_norm/D");  
  int T_mc_fnucl;         Tout->Branch("mc_fnucl", &T_mc_fnucl, "mc_fnucl/I");
  //
  double T_ebeam;         Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  //kine
  // -- vertex (only for SIMC) ---
  double T_veE;           if (generator.compare("simc")==0) Tout->Branch("veE", &T_veE, "veE/D");
  double T_vetheta;       if (generator.compare("simc")==0) Tout->Branch("vetheta", &T_vetheta, "vetheta/D");
  double T_vQ2;           if (generator.compare("simc")==0) Tout->Branch("vQ2", &T_vQ2, "vQ2/D");
  double T_vepsilon;      if (generator.compare("simc")==0) Tout->Branch("vepsilon", &T_vepsilon, "vepsilon/D"); // calculated using general eqn.
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
  //HCAL
  double T_eHCAL;         Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  double T_xHCAL;         Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  double T_yHCAL;         Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
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
  //GEM
  double T_nhitsGEM;      Tout->Branch("nhitsGEM", &T_nhitsGEM, "nhitsGEM/D");
  double T_ngoodhitsGEM;  Tout->Branch("ngoodhitsGEM", &T_ngoodhitsGEM, "ngoodhitsGEM/D");
  double T_trchi2ndf;     Tout->Branch("trchi2ndf", &T_trchi2ndf, "trchi2ndf/D");

  // Do the energy loss calculation here (only for g4sbs generator)
  double ebeam = sbsconf.GetEbeam(); // gets overwritten in the event loop

  // HCAL cut definitions
  vector<double> dx_p_cut; jmgr->GetVectorFromSubKey<double>(key,"dx_p_cut", dx_p_cut);
  double sbs_kick = abs(dx_p_cut[0]);
  vector<double> dy_p_cut; jmgr->GetVectorFromSubKey<double>(key,"dy_p_cut", dy_p_cut);
  vector<double> hcal_active_area = cut::hcal_active_area_simu(1,1); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p_cut[1], dx_p_cut[1], dy_p_cut[1], hcal_active_area);
  // varying safety margin width by +/- 10% in vertical direction
  vector<double> hcal_safety_margin_p10p = cut::hcal_safety_margin(dx_p_cut[1]*1.1,dx_p_cut[1]*1.1,dy_p_cut[1],hcal_active_area);
  vector<double> hcal_safety_margin_m10p = cut::hcal_safety_margin(dx_p_cut[1]*.9,dx_p_cut[1]*.9,dy_p_cut[1],hcal_active_area);
  TH2F *h2_xyHCAL_p = util_pd::TH2FHCALface_xy_simu("h2_xyHCAL_p",sbs_kick);
  TH2F *h2_xyHCAL_p_nfc = util_pd::TH2FHCALface_xy_simu("h2_xyHCAL_p_nfc",sbs_kick);

  // reading W2 offset and W cut limits
  double dy_offset = jmgr->GetValueFromSubKey<double>(key,"dy_offset");
  double W2_offset = jmgr->GetValueFromSubKey<double>(key,"W2_offset");
  std::vector<double> W_cutR; jmgr->GetVectorFromSubKey<double>(key,"W_cutR",W_cutR);
  // reading BB fiducial cut limits
  std::vector<double> bbfidu_cutR; jmgr->GetVectorFromSubKey<double>(key,"bbfidu_cutR",bbfidu_cutR);

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = jmgr->GetValueFromSubKey<double>(key,"hcal_voffset");
  double hcal_hoffset = jmgr->GetValueFromSubKey<double>(key,"hcal_hoffset");
  double hcal_zoffset = jmgr->GetValueFromSubKey<double>(key,"hcal_zoffset");
  //vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetHCALtheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = (sbsconf.GetHCALdist()+hcal_zoffset)*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // calculating MC normalizetion factors 
  double mc_omega, lumi;
  vector<long double> totNtriesnCh; util_pd::GetTotNtriesnCh(sjobs,totNtriesnCh);
  // -- Rejection sampling (RS) flags and weight
  bool usingRS = sjobs[0].usingRS;
  double maxwtRS = sjobs[0].maxwtRS;
  // report
  if (usingRS) {
    std::cout << "----\n";
    std::cout << "Using RS: " << usingRS << std::endl;
    std::cout << "Max Weight RS: " << maxwtRS << std::endl;
    std::cout << "Total ntries: " << totNtriesnCh[0] << std::endl;
    std::cout << "----\n";
  }

  // implementing hashtable with norm info for efficiency
  std::unordered_map<std::string,SimuJob> mnorm;
  for (auto & sjob: sjobs) mnorm[sjob.rfname] = sjob;

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  double charge;
  long double totntries = totNtriesnCh[0];
  long nevent = 0, nevents = C->GetEntries(), ngoodevs = 0; 
  int treenum = 0, currenttreenum = 0, treeitr = 0;
  while (C->GetEntry(nevent++)) {
   
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      // apply global cuts efficiently (AJRP method)
      GlobalCut->UpdateFormulaLeaves();
      
      // getting normalization factors per run
      const char* rftemp = C->GetFile()->GetName();
      SimuJob sjtemp = mnorm[rftemp];
      ebeam = sjtemp.ebeam; charge = sjtemp.charge;
      lumi = sjtemp.lumi; mc_omega = sjtemp.genvol; 
      if (sjtemp.usingRS) maxwtRS = sjtemp.maxwtRS;
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (!passedgCut) continue;
    ngoodevs++;

    // cross section weighted normalization factor
    weight = usingRS ? maxwtRS*mc_omega*lumi/totntries : mc_sigma*mc_omega*lumi/totntries;
    weight_norm = usingRS ? maxwtRS*mc_omega*lumi/totntries/charge : mc_sigma*mc_omega*lumi/totntries/charge;

    // kinematic parameters
    double ebeam_corr = ebeam; //- MeanEloss;
    if (generator.compare("simc")==0) ebeam_corr = mc_ebeam;
    double precon = p[0]; //+ MeanEloss_outgoing

    // constructing the 4 vectors
    /* Reaction    : e + e' -> N + N'
       Conservation: Pe + Peprime = PN + PNprime */
    TVector3 vertex(0, 0, vz[0]);
    TLorentzVector Pe(0,0,ebeam_corr,ebeam_corr);   // incoming e-
    TLorentzVector Peprime(px[0] * (precon/p[0]),   // scattered e-
			   py[0] * (precon/p[0]),
			   pz[0] * (precon/p[0]),
			   precon);                 
    TLorentzVector PN;                              // target nucleon
    kine::SetPN(Ntype, PN);
    TLorentzVector PNprime;                         // Recoil nucleon 4-vector
    TLorentzVector q = Pe - Peprime;                // 4-momentum of virtual photon

    double etheta = kine::etheta(Peprime);
    double ephi = kine::ephi(Peprime);
    double pelas = kine::pelas(ebeam_corr, etheta, Ntype);
    double thelas = kine::thelas(ebeam_corr, precon, Ntype);

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
    if (model == 0) {
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
    } else if (model == 1) {
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
    } else if (model == 2) {
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

    T_ebeam = Pe.E();

    // at vertex ("true" values)
    T_veE = mc_veE;
    T_vetheta = mc_vetheta;
    T_vQ2 = kine::Q2(Pe.E(),T_veE,T_vetheta);
    T_vepsilon = kine::epsilon_general(T_vetheta,T_vQ2,T_ebeam-T_veE);

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon+W2_offset;
    T_W = Wrecon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pelas = pelas;
    T_pN_exp = pN_expect;
    T_thN_exp = thetaN_expect;
    T_epsilon = epsilon;
    T_epsilon_p = epsilon_p;

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

    T_eHCAL = eHCAL;
    T_xHCAL = xHCAL;
    T_yHCAL = yHCAL;
    T_rblkHCAL = rblkHCAL;
    T_cblkHCAL = cblkHCAL;
    T_idblkHCAL = idblkHCAL;
    T_atimeHCAL = atimeHCAL;
    T_tdcHCAL = tdcHCAL[0];

    T_nhitsGEM = nhitsGEM[0];
    T_ngoodhitsGEM = ngoodhitsGEM[0];
    T_trchi2ndf = trchi2ndf[0];

    T_mc_fnucl = int(mc_fnucl);

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL - xyHCAL_exp[0];  
    double dy = yHCAL - xyHCAL_exp[1]; 

    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_xHCAL_exp_p = xyHCAL_exp[0]-sbs_kick;
    T_dx = dx;
    T_dy = dy+dy_offset;

    /* Calculating thpq (p hypothesis) */
    TVector3 HCAL_pos = HCAL_origin + xHCAL*HCAL_axes[0] + yHCAL*HCAL_axes[1];
    TVector3 n_dir = (HCAL_pos - vertex).Unit();
    T_thpq_n = acos(n_dir.Dot(pNhat));
    // p 
    double BdL = sbsscalefield * expconst::sbsmaxfield_simu * expconst::sbsdipolegap;
    double proton_thetabend = 0.3 * BdL / PNprime.Vect().Mag();  // p*theta = 0.3*BdL
    double proton_deflection = tan(proton_thetabend)*(sbsconf.GetHCALdist()+hcal_zoffset-(sbsconf.GetSBSdist()+expconst::sbsdipolegap/2.0));
    T_p_def = proton_deflection;
    TVector3 p_dir = (HCAL_pos + proton_deflection*HCAL_axes[0] - vertex);
    T_thpq_p = acos(p_dir.Unit().Dot(pNhat));

    // HCAL active area and safety margin cuts [Fiducial region]
    ARCut = cut::inHCAL_activeA(xHCAL,yHCAL,hcal_active_area);
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

    // defining W cut
    WCut = Wrecon >= W_cutR[0] && Wrecon <= W_cutR[1];

    // W cut
    if (WCut&&bbfiduCut&&T_eHCAL>0) {
      if (abs(dy)<0.4) h_dxHCAL_nfc->Fill(dx, weight);
      h_dyHCAL_nfc->Fill(dy, weight);
      if (fiduCut) {
	if (abs(dy)<0.4) {
	  h_dxHCAL->Fill(dx, weight);
	  h_dx_w_p_def->Fill(dx+T_p_def, weight);
	}
	h_dyHCAL->Fill(dy, weight);
	// dx dist. for p & n separately using MC info
	if (int(mc_fnucl)==1) {
	  h_dxHCAL_p->Fill(dx, weight);
	  h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick, weight);
	}else if (process.compare("heep")==0) {
	  std::cerr << "*!* Invalid final state nuclei!" << std::endl; 
	  std::exit(1);
	}

	h2_rcHCAL->Fill(cblkHCAL, rblkHCAL, weight);
	// p & n spots
	h2_dxdyHCAL->Fill(dy, dx, weight);

	// hit map to show p & n in fiducial region
	if (pCut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick, weight);
      }
      if (pCut) h2_xyHCAL_p_nfc->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick, weight);
    }

    // fiducial cut but no W cut
    if (Wrecon>0) { 
      if (fiduCut) { 
	h_W->Fill(Wrecon,weight);
	if (pCut) { 
	  h_W_cut->Fill(Wrecon,weight);
	  h_W2_cut->Fill(W2recon,weight);
	} else {
	  h_W_acut->Fill(Wrecon,weight);
	}
      }
    }
      
    Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;

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
  util_pd::PlotFiduCut(2,hcal_active_area,hcal_safety_margin);
  c1->cd(4); //
  h2_xyHCAL_p->Draw("colz");
  util_pd::PlotFiduCut(2,hcal_active_area,hcal_safety_margin);
  c1->SaveAs(Form("%s[",outPlot.Data())); c1->SaveAs(Form("%s",outPlot.Data())); c1->Write();
  //**** -- ***//

  /**** Canvas 2 (dx & dy) ****/
  TCanvas *c2 = util_pd::TC("c2",2,2);
  std::vector<double> hdx_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dxHCAL_fitR", hdx_fitR);
  std::vector<double> hdy_fitR; jmgr->GetVectorFromSubKey<double>(key,"h_dyHCAL_fitR", hdy_fitR);
  gStyle->SetOptFit(1);
  c2->cd(1); //
  gPad->SetGridx();
  TF1 *fdxp = fit::fit_1gs_nbg(hdx_fitR,h_dxHCAL);
  double dxpM = fdxp->GetParameter(1); double dxpS = fdxp->GetParameter(2);
  c2->cd(2); //
  gPad->SetGridx();
  TF1 *fdyp = fit::fit_1gs_nbg(hdy_fitR,h_dyHCAL);
  double dypM = fdyp->GetParameter(1); double dypS = fdyp->GetParameter(2);
  c2->cd(3); //
  gPad->SetGridx();
  TF1 *fdxp_nfc = fit::fit_1gs_nbg(hdx_fitR,h_dxHCAL_nfc);
  double dxpM_nfc = fdxp_nfc->GetParameter(1); double dxpS_nfc = fdxp_nfc->GetParameter(2);
  c2->cd(4); //
  gPad->SetGridx();
  TF1 *fdyp_nfc = fit::fit_1gs_nbg(hdy_fitR,h_dyHCAL_nfc);
  double dypM_nfc = fdyp_nfc->GetParameter(1); double dypS_nfc = fdyp_nfc->GetParameter(2);
  c2->SaveAs(Form("%s",outPlot.Data())); c2->Write();
  //**** -- ***//

  /**** Canvas 3 (Deflection and W2) ****/
  TCanvas *c3 = util_pd::TC("c3",1,2);
  c3->SetGridx();
  gStyle->SetOptFit(1);
  //
  c3->cd(1);
  gPad->SetGridx();
  std::vector<double> hpdef_fitR{-0.5,0.5,1.2,1.2};
  TF1 *fpdef = fit::fit_1gs_nbg(hpdef_fitR,h_dx_w_p_def);
  double pdefM = fpdef->GetParameter(1); double pdefS = fpdef->GetParameter(2);
  h_dx_w_p_def->Draw("same");
  h_dx_w_p_def->SetLineColor(kBlack);
  //
  c3->cd(2);
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
  pt->AddText(Form(" Analyzing %s generated QE events for SBS%d-SBS%dp settings",generator.c_str(),conf,sbsmag));
  pt->AddText(Form(" Analysis model: %d",model));
  pt->AddText(Form(" HCAL offsets: v = %.4f, h = %.4f, z = %.4f",hcal_voffset,hcal_hoffset,hcal_zoffset));
  pt->AddText(Form(" W2 offset = %.4f, dy offset = %.4f",W2_offset,dy_offset));
  pt->AddText(Form(" Total # events analyzed: %ld",nevents));
  pt->AddText(Form(" Global cuts: "));
    std::string tmpstr = "";
  for (std::size_t i=0; i<gCutList.size(); i++) {
    if (i>0 && i%3==0) {pt->AddText(Form(" %s",tmpstr.c_str())); tmpstr="";}
    tmpstr += gCutList[i] + ", "; 
  }
  if (!tmpstr.empty()) pt->AddText(Form(" %s",tmpstr.c_str()));
  pt->AddText(Form(" # events passed global cuts: %ld", ngoodevs));
  pt->AddText(" Elastic cuts: ");
  pt->AddText(Form(" Inbuilt W cut: %.2f #leq W #leq %.2f GeV/c",W_cutR[0],W_cutR[1]));
  pt->AddText(Form(" Inbuilt p cut (#Deltax): Mean = %.4f, %.1f#sigma = %.4f",dx_p_cut[0],dx_p_cut[2],dx_p_cut[1]));
  pt->AddText(Form(" Inbuilt p cut (#Deltay): Mean = %.4f, %.1f#sigma = %.4f",dy_p_cut[0],dy_p_cut[2],dy_p_cut[1]));
  pt->AddText(Form(" Inbuilt BB fiducial cut: |fpX-0.9*fpTh-%.2f| #leq %.2f",bbfidu_cutR[0],bbfidu_cutR[1]));
  pt->AddText(" Fit info: ");
  pt->AddText(" p peak, w/ fiducial cut: dxpM,dxpS,dypM,dypS ");
  pt->AddText(Form(" %.5f,%.5f,%.5f,%.5f",dxpM,dxpS,dypM,dypS));
  pt->AddText(" p peak, w/o fiducial cut: dxpM_nfc,dxpS_nfc,dypM_nfc,dypS_nfc ");
  pt->AddText(Form(" %.5f,%.5f,%.5f,%.5f",dxpM_nfc,dxpS_nfc,dypM_nfc,dypS_nfc));
  if (usingRS) {
    pt->AddText(Form(" Rejection Sampling (RS) Summary:")); 
    pt->AddText(Form(" Chosen maximum weight: %f",maxwtRS)); 
    //pt->AddText(Form(" Total # tries: %.0Lf",totNtriesNCh[0])); 
    TText *t5 = pt->GetLineWith(" Rejection"); t5->SetTextColor(kMagenta+2);
  }
  sw->Stop();
  pt->AddText(Form("Macro processing time: CPU %.1fs | Real %.1fs",sw->CpuTime(),sw->RealTime()));
  TText *t1 = pt->GetLineWith("Configfile"); t1->SetTextColor(kRed);
  TText *t2 = pt->GetLineWith(" Global"); t2->SetTextColor(kBlue);
  TText *t3 = pt->GetLineWith(" Elastic"); t3->SetTextColor(kBlue);
  TText *t4 = pt->GetLineWith(" Fit info"); t4->SetTextColor(kBlue);
  TText *t6 = pt->GetLineWith("Macro"); t6->SetTextColor(kGreen+3);
  pt->Draw();  
  cSummary->SaveAs(Form("%s",outPlot.Data())); cSummary->SaveAs(Form("%s]",outPlot.Data())); cSummary->Write();  
  //**** -- ***//

  std::cout << "\n----Fit info----" << "\n";
  std::cout << "dxpM,dxpS,dypM,dypS,dxpM_nfc,dxpS_nfc,dypM_nfc,dypS_nfc" << "\n";
  std::cout << dxpM<<","<<dxpS<<","<<dypM<<","<<dypS<<","<<dxpM_nfc<<","<<dxpS_nfc<<","<<dypM_nfc<<","<<dypS_nfc<<"\n";
  std::cout << "----------------" << "\n\n";

  std::cout << "------" << std::endl;
  std::cout << " Summary plots  : " << outPlot << std::endl;
  std::cout << " Output ROOT file  : " << outFile << std::endl;
  std::cout << "------" << std::endl << std::endl;

  std::cout << "CPU time = " << sw->CpuTime() << "s. Real time = " << sw->RealTime() << "s.\n\n";

  Tout->Write("",TObject::kOverwrite);
  h_Q2->Write();
  h_dpel->Write(); h_W->Write();
  h_W_cut->Write(); h_W_acut->Write();
  h_W2_cut->Write();
  h_dxHCAL->Write(); h_dxHCAL_nfc->Write();
  h_dyHCAL->Write(); h_dyHCAL_nfc->Write();
  h_dxHCAL_p->Write(); h2_xyHCAL_p->Write();
  h2_xyHCAL_p_nfc->Write();
  h2_rcHCAL->Write(); h2_dxdyHCAL->Write();
  h_dx_w_p_def->Write();
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
  ** generator : MC enevt generator (Valid options: g4sbs, simc)
  ** prefix_to_filebase : <prefix>_sbs<sbsconfig>_sbs<sbsmagfield>p_<generator>_<process>.root
  ** Njobs_to_ana : # MC jobs to analyze
  ** rootfile_dir : Directory name w/ path containing the MC ROOT files to analyze
  ** global_cut : set of global cuts to apply at the start of event processing
  ** SBS_field : 
  ** hcal_v(h)offset :
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
  ** W_cutR : W cut range.
     - W_cutR[0] : lower limit
     - W_cutR[1] : upper limit
*/
