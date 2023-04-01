/* 
   This macro will perform elastic (LH2) analysis for GMn using MC data.
   E.g. Config. File: sbs14-sbs70p-simu/conf_elas_ana_simu.json
   -----
   P. Datta  Created  11-05-2022 
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

#include "../../include/gmn_ana.h"
#include "../../dflay/src/JSONManager.cxx"

int elas_ana_simu (const char *configfilename, std::string filebase="siout/test_elas_ana_simu")
{
  gErrorIgnoreLevel = kError; // Ignores all ROOT warnings

  // Define a clock to get macro processing time
  TStopwatch *sw = new TStopwatch(); sw->Start();

  // reading input config file ---------------------------------------
  JSONManager *jmgr = new JSONManager(configfilename);

  // parsing trees
  std::string rootfile_dir = jmgr->GetValueFromKey_str("rootfile_dir");
  TChain *C = new TChain("T");
  C->Add(rootfile_dir.c_str());
  if (C->GetEntries()==0) {std::cerr << "*!* No ROOT file!" << std::endl; throw;}
 
  // seting up the desired SBS configuration
  int conf = jmgr->GetValueFromKey<int>("SBS_config");
  int sbsmag = jmgr->GetValueFromKey<int>("SBS_magnet_percent");
  SBSconfig sbsconf(conf, sbsmag);
  cout << sbsconf;

  // Choosing the model of calculation
  // model 0 => uses reconstructed p as independent variable
  // model 1 => uses reconstructed angles as independent variable
  // model 2 => uses 4-vector calculation
  int model = jmgr->GetValueFromKey<int>("model");
  if (model == 0) std::cout << "Using model 0 [recon. p as indep. var.] for analysis.." << std::endl;
  else if (model == 1) std::cout << "Using model 1 [recon. angle as indep. var.] for analysis.." << std::endl;
  else if (model == 2) std::cout << "Using model 2 [4-vector calculation] for analysis.." << std::endl;
  else { std::cerr << "Enter a valid model number! **!**" << std::endl; throw; }

  // parameters for normalization
  std::string targetType = jmgr->GetValueFromKey_str("target_type");
  double I_beam = jmgr->GetValueFromKey<double>("beam_current");
  // total generated simulation events per job
  double ngen_total = jmgr->GetValueFromKey<double>("ngen_total"); 
  double lumi = kine::Luminosity(I_beam, targetType);

  // choosing nucleon type 
  std::string Ntype = jmgr->GetValueFromKey_str("Ntype");

  // setting up global cuts
  std::string gcut = jmgr->GetValueFromKey_str("global_cut");
  TCut globalcut = gcut.c_str();
  TTreeFormula *GlobalCut = new TTreeFormula("GlobalCut", globalcut, C);

  // setting up ROOT tree branch addresses ---------------------------------------
  int maxNtr=1000;
  C->SetBranchStatus("*",0);
  // beam energy - Probably we should take an average over 100 events
  // double HALLA_p;
  // setrootvar::setbranch(C, "HALLA_p", "", &HALLA_p);

  // bbcal clus var
  double eSH; setrootvar::setbranch(C,"bb.sh","e",&eSH);
  double ePS; setrootvar::setbranch(C,"bb.ps","e",&ePS);
  
  // hcal clus var
  double eHCAL, xHCAL, yHCAL, rblkHCAL, cblkHCAL, idblkHCAL;
  std::vector<std::string> hcalclvar = {"e","x","y","rowblk","colblk","idblk"};
  std::vector<void*> hcalclvar_mem = {&eHCAL,&xHCAL,&yHCAL,&rblkHCAL,&cblkHCAL,&idblkHCAL};
  setrootvar::setbranch(C, "sbs.hcal", hcalclvar, hcalclvar_mem);

  // track var
  double ntrack, p[maxNtr],px[maxNtr],py[maxNtr],pz[maxNtr],xTr[maxNtr],yTr[maxNtr],thTr[maxNtr],phTr[maxNtr];
  double vx[maxNtr],vy[maxNtr],vz[maxNtr];
  double xtgt[maxNtr],ytgt[maxNtr],thtgt[maxNtr],phtgt[maxNtr];
  std::vector<std::string> trvar = {"n","p","px","py","pz","x","y","th","ph","vx","vy","vz","tg_x","tg_y","tg_th","tg_ph"};
  std::vector<void*> trvar_mem = {&ntrack,&p,&px,&py,&pz,&xTr,&yTr,&thTr,&phTr,&vx,&vy,&vz,&xtgt,&ytgt,&thtgt,&phtgt};
  setrootvar::setbranch(C,"bb.tr",trvar,trvar_mem);

  //MC variables
  double mc_sigma, mc_omega, mc_fnucl;
  std::vector<std::string> mc = {"mc_sigma","mc_omega","mc_fnucl"};
  std::vector<void*> mc_mem = {&mc_sigma,&mc_omega,&mc_fnucl};
  setrootvar::setbranch(C,"MC",mc,mc_mem);

  //MC variables (SIMC)
  double simc_sigma, simc_Weight;
  std::vector<std::string> simc = {"simc_sigma","simc_Weight"};
  std::vector<void*> simc_mem = {&simc_sigma,&simc_Weight};
  setrootvar::setbranch(C,"MC",simc,simc_mem);

  // turning on the remaining branches we use for the globalcut
  C->SetBranchStatus("bb.gem.track.nhits", 1);
  C->SetBranchStatus("bb.etot_over_p", 1);

  // defining the outputfile
  TString outFile = Form("%s_sbs%d_sbs%dp_model%d.root", 
			 filebase.c_str(), sbsconf.GetSBSconf(), sbsconf.GetSBSmag(), model);
  TFile *fout = new TFile(outFile.Data(), "RECREATE");

  // defining histograms
  TH1F *h_W2 = new TH1F("h_W2","h_W2",300,-3.,3.);
  TH1F *h_W = util_pd::TH1FhW("h_W");
  TH1F *h_W_cut = util_pd::TH1FhW("h_W_cut");
  TH1F *h_W_acut = util_pd::TH1FhW("h_W_acut");
  TH1D *h_dpel = new TH1D("h_dpel",";p/p_{elastic}(#theta)-1;",100,-0.3,0.3);
  
  TH1F *h_Q2 = util_pd::TH1FhQ2("h_Q2", conf);
  vector<double> hdx_lim; jmgr->GetVectorFromKey<double>("h_dxHCAL_lims", hdx_lim);
  vector<double> hdy_lim; jmgr->GetVectorFromKey<double>("h_dyHCAL_lims", hdy_lim);
  TH1F *h_dxHCAL = new TH1F("h_dxHCAL","; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);
  TH1F *h_dyHCAL = new TH1F("h_dyHCAL","; y_{HCAL} - y_{exp} (m);",int(hdy_lim[0]),hdy_lim[1],hdy_lim[2]);
  TH1F *h_dxHCAL_p = new TH1F("h_dxHCAL_p","mc_fnucl = p; x_{HCAL} - x_{exp} (m);",int(hdx_lim[0]),hdx_lim[1],hdx_lim[2]);

  TH2F *h2_rcHCAL = util_pd::TH2FHCALface_rc("h2_rcHCAL");
  TH2F *h2_dxdyHCAL = util_pd::TH2FdxdyHCAL("h2_dxdyHCAL");
  TH2F *h2_xyHCAL_p = util_pd::TH2FHCALface_xy_simu("h2_xyHCAL_p");

  // Defining interesting ROOT tree branches 
  TTree *Tout = new TTree("Tout", "");
  //cuts
  bool WCut;            Tout->Branch("WCut", &WCut, "WCut/B");
  bool pCut;            Tout->Branch("pCut", &pCut, "pCut/B");
  bool fiduCut;         Tout->Branch("fiduCut", &fiduCut, "fiduCut/B");
  //MC related
  double weight;        Tout->Branch("weight", &weight, "weight/D");  
  int T_mc_fnucl;       Tout->Branch("mc_fnucl", &T_mc_fnucl, "mc_fnucl/I");
  //
  double T_ebeam;       Tout->Branch("ebeam", &T_ebeam, "ebeam/D");
  //kine
  double T_nu;          Tout->Branch("nu", &T_nu, "nu/D");
  double T_Q2;          Tout->Branch("Q2", &T_Q2, "Q2/D");
  double T_W2;          Tout->Branch("W2", &T_W2, "W2/D");
  double T_W;           Tout->Branch("W", &T_W, "W/D");
  double T_dpel;        Tout->Branch("dpel", &T_dpel, "dpel/D");
  double T_ephi;        Tout->Branch("ephi", &T_ephi, "ephi/D");
  double T_etheta;      Tout->Branch("etheta", &T_etheta, "etheta/D");
  double T_pcentral;    Tout->Branch("pcentral", &T_pcentral, "pcentral/D");
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
  //HCAL
  double T_eHCAL;       Tout->Branch("eHCAL", &T_eHCAL, "eHCAL/D"); 
  double T_xHCAL;       Tout->Branch("xHCAL", &T_xHCAL, "xHCAL/D"); 
  double T_yHCAL;       Tout->Branch("yHCAL", &T_yHCAL, "yHCAL/D"); 
  double T_xHCAL_exp;   Tout->Branch("xHCAL_exp", &T_xHCAL_exp, "xHCAL_exp/D"); 
  double T_yHCAL_exp;   Tout->Branch("yHCAL_exp", &T_yHCAL_exp, "yHCAL_exp/D"); 
  double T_dx;          Tout->Branch("dx", &T_dx, "dx/D"); 
  double T_dy;          Tout->Branch("dy", &T_dy, "dy/D");  

  // Do the energy loss calculation here ...........

  // HCAL cut definitions
  double sbs_kick = jmgr->GetValueFromKey<double>("sbs_kick");
  vector<double> dx_p; jmgr->GetVectorFromKey<double>("dx_p", dx_p);
  vector<double> dy_p; jmgr->GetVectorFromKey<double>("dy_p", dy_p);
  double Nsigma_cut_dx_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dx_p");
  double Nsigma_cut_dy_p = jmgr->GetValueFromKey<double>("Nsigma_cut_dy_p");
  vector<double> hcal_active_area = cut::hcal_active_area_simu(); // Exc. 1 blk from all 4 sides
  vector<double> hcal_safety_margin = cut::hcal_safety_margin(dx_p[1], dx_p[1], dy_p[1], hcal_active_area);

  // elastic cut limits
  double Wmin = jmgr->GetValueFromKey<double>("Wmin");
  double Wmax = jmgr->GetValueFromKey<double>("Wmax");

  // costruct axes of HCAL CoS in Hall CoS
  double hcal_voffset = jmgr->GetValueFromKey<double>("hcal_voffset");
  double hcal_hoffset = jmgr->GetValueFromKey<double>("hcal_hoffset");
  vector<TVector3> HCAL_axes; kine::SetHCALaxes(sbsconf.GetSBStheta_rad(), HCAL_axes);
  TVector3 HCAL_origin = sbsconf.GetHCALdist()*HCAL_axes[2] + hcal_voffset*HCAL_axes[0] + hcal_hoffset*HCAL_axes[1];

  // looping through the tree ---------------------------------------
  std::cout << std::endl;
  long nevent = 0, nevents = C->GetEntries(); 
  int treenum = 0, currenttreenum = 0;
  while (C->GetEntry(nevent++)) {
   
    // print progress 
    if( nevent % 1000 == 0 ) std::cout << nevent << "/" << nevents << "\r";
    std::cout.flush();

    // apply global cuts efficiently (AJRP method)
    currenttreenum = C->GetTreeNumber();
    if (nevent == 1 || currenttreenum != treenum) {
      treenum = currenttreenum;
      GlobalCut->UpdateFormulaLeaves();
    } 
    bool passedgCut = GlobalCut->EvalInstance(0) != 0;   
    if (!passedgCut) continue;

    // cross section weighted normalization factor
    //weight = mc_sigma*mc_omega*lumi / ngen_total;
    weight = simc_Weight;

    // kinematic parameters
    double ebeam = sbsconf.GetEbeam();       // Expected beam energy (GeV) [Get it from EPICS, eventually]
    double ebeam_corr = ebeam; //- MeanEloss;
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
    TLorentzVector PN;                              // target nucleon [Ntype ??]
    kine::SetPN(Ntype, PN);

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
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
    } else if (model == 1) {
      nu = Pe.E() - pcentral;
      pN_expect = kine::pN_expect(nu, Ntype);
      thetaN_expect = acos((Pe.E() - pcentral*cos(etheta)) / pN_expect);
      pNhat = kine::qVect_unit(thetaN_expect, phiN_expect);
      Q2recon = kine::Q2(Pe.E(), Peprime.E(), etheta);
      W2recon = kine::W2(Pe.E(), Peprime.E(), Q2recon, Ntype);
    } else if (model == 2) {
      TLorentzVector q = Pe - Peprime; // 4-momentum of virtual photon
      nu = q.E();
      pNhat = q.Vect().Unit();
      Q2recon = -q.M2();
      W2recon = (PN + q).M2();
    }
    h_Q2->Fill(Q2recon); 
    double Wrecon = sqrt(max(0., W2recon));
    double dpel = Peprime.E()/pcentral - 1.0; h_dpel->Fill(dpel);

    T_ebeam = Pe.E();

    T_nu = nu;
    T_Q2 = Q2recon;
    T_W2 = W2recon;
    T_W = Wrecon;
    T_dpel = dpel;
    T_ephi = ephi;
    T_etheta = etheta;
    T_pcentral = pcentral;

    T_vz = vz[0];
    T_trP = p[0];
    T_trX = xTr[0];
    T_trY = yTr[0];
    T_trTh = thTr[0];
    T_trPh = phTr[0];

    T_ePS = ePS;
    T_eSH = eSH;

    T_eHCAL = eHCAL;
    T_xHCAL = xHCAL;
    T_yHCAL = yHCAL;

    T_mc_fnucl = int(mc_fnucl);

    // Expected position of the q vector at HCAL
    vector<double> xyHCAL_exp; // xyHCAL_exp[0] = xHCAL_exp & xyHCAL_exp[1] = yHCAL_exp
    kine::GetxyHCALexpect(vertex, pNhat, HCAL_origin, HCAL_axes, xyHCAL_exp);
    double dx = xHCAL - xyHCAL_exp[0];  
    double dy = yHCAL - xyHCAL_exp[1]; 

    T_xHCAL_exp = xyHCAL_exp[0];
    T_yHCAL_exp = xyHCAL_exp[1];
    T_dx = dx;
    T_dy = dy;

    // HCAL active area and safety margin cuts [Fiducial region]
    /* Using fiducut for LH2 data doesn't make sense. Still keeping the branch just 
       in case we need to check something. Won't be using this on any histograms. */
    bool AR_cut = cut::inHCAL_activeA(xHCAL, yHCAL, hcal_active_area);
    bool FR_cut = cut::inHCAL_fiducial(xyHCAL_exp[0], xyHCAL_exp[1], sbs_kick, hcal_safety_margin);
    fiduCut = AR_cut && FR_cut;
    // HCAL cuts
    pCut = pow((dx-dx_p[0]) / (dx_p[1]*Nsigma_cut_dx_p), 2) + pow((dy-dy_p[0]) / (dy_p[1]*Nsigma_cut_dy_p), 2) <= 1.;

    WCut = Wrecon >= Wmin && Wrecon <= Wmax;

    // W cut
    if (WCut) {
      // fiducial cut  (not using it for LH2 data)
      // if (fiduCut) {
      h_dxHCAL->Fill(dx, weight);
      h_dyHCAL->Fill(dy, weight);
      // dx dist. for p & n separately using MC info
      if (int(mc_fnucl)==1) {
	h_dxHCAL_p->Fill(dx, weight);
      } else {
	std::cerr << "*!* Invalid final state nuclei!" << std::endl; 
	throw;
      }

      h2_rcHCAL->Fill(cblkHCAL, rblkHCAL);
      // p & n spots
      h2_dxdyHCAL->Fill(dy, dx);

      // hit map to show p & n in fiducial region
      if (pCut) h2_xyHCAL_p->Fill(xyHCAL_exp[1], xyHCAL_exp[0] - sbs_kick);
      // }
    }

    // fiducial cut but no W cut
    // if (fiduCut) {  (not using it for LH2 data)
    h_W2->Fill(W2recon, weight);
    h_W->Fill(Wrecon, weight);
    if (pCut) { 
      h_W_cut->Fill(Wrecon);
    } else {
      h_W_acut->Fill(Wrecon);
    }
    // }
      
    Tout->Fill();
  } // event loop
  std::cout << std::endl << std::endl;

  TCanvas *c1 = new TCanvas("c1", "c1", 1200, 1000);
  c1->Divide(2,2);

  c1->cd(1); h2_dxdyHCAL->Draw("colz");
  TEllipse Ep_p;
  Ep_p.SetFillStyle(0); Ep_p.SetLineColor(2); Ep_p.SetLineWidth(2);
  Ep_p.DrawEllipse(dy_p[0], dx_p[0], Nsigma_cut_dy_p*dy_p[1], Nsigma_cut_dx_p*dx_p[1], 0,360,0);
 
  c1->cd(2);
  h_W->Draw(); h_W->SetLineColor(1);
  h_W_cut->Draw("same"); h_W_cut->SetLineColor(2);
  h_W_acut->Draw("same");

  c1->cd(3);
  h_dxHCAL->Draw();
  // h2_xyHCAL_p->Draw("colz");
  // util_pd::DrawArea(hcal_active_area);
  // util_pd::DrawArea(hcal_safety_margin,4);

  c1->cd(4);
  h_W2->Draw();
  // h_dyHCAL->Draw();
  // h2_xyHCAL_n->Draw("colz");
  // util_pd::DrawArea(hcal_active_area);
  // util_pd::DrawArea(hcal_safety_margin,4);

  // outFile.ReplaceAll(".root",".png");
  // c1->Print(outFile.Data(),"png");

  cout << "------" << endl;
  cout << " Output file : " << outFile << endl;
  cout << "------" << endl << endl;

  sw->Stop();
  cout << "CPU time elapsed = " << sw->CpuTime() << " s. Real time = " << sw->RealTime() << " s. " << endl << endl;

  c1->Write();
  fout->Write();
  sw->Delete();
  delete jmgr;
  return 0;
}
