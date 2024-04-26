#include "../include/EMFFFits.h"

void EMFFFits::StripGDandMu(G_t const kG, double const Q2, double & value) {
  /* All the fits implemented here returns a ratio of Sachs FF (G) to 
     GDip. Furthermore, the returned magnetic FF values are divided by
     appropriate nucleon anomaleous magnetic moment.
     This fn takes such FF values and returns pure Sachs FF by cancelling
     out GDip and mu contributions appropriately.*/
  value *= EMFFFits::GetGDip(Q2);
  switch (kG) {
  case G_t::kGMp:
    value *= constant::mup;
    break;
  case G_t::kGMn:
    value *= constant::mun;
  }
}

G_t EMFFFits::StrToG_t(const std::string& str) {
  /* Returns G_t object according to input string */
  if (str == "GEp") {
    return G_t::kGEp;
  } else if (str == "GMp") {
    return G_t::kGMp;
  } else if (str == "GEn") {
    return G_t::kGEn;
  } else if (str == "GMn") {
    return G_t::kGMn;
  } else {
    throw std::invalid_argument("Invalid string for conversion to enum");
  }
}

// ########################
// ## Galster Fit (1971) ##
// ########################
int Galster1971::GalsterFit(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err) {
  // This fit has been mentioned in the GMn proposal. 
  // It returns GEn/GEp, but it seems to me that GEp for their case is essentially GDip
  
  // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  if (kID!=3) {
    std::cerr<<"*** ERROR***, Galster fit only supports kID=3 ie GEn"<<std::endl;
    GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
    return -1;
  }

  //// Applying parametrization formula
  double tau = kID<3 ? kine::tau(kQ2,"p") : kine::tau(kQ2,"n");
  double numerator = -constant::mun;
  double denominator = 1. + 5.6*tau;

  GNGD_Fit[0] = numerator / denominator;
  GNGD_Err[0] = -1000.;

  return 0;
  
}

// ######################
// ## Kelly Fit (2004) ##
// ######################
int Kelly2004::KellyFit(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err) {
  // ** Link to the original paper:
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.70.068202

  // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  if (kID<1 || kID>4){
    std::cerr<<"*** ERROR***, kID is not any of [1->GEp, 2->GMp, 3->GEn, 4->GMn]"<<std::endl;
    GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
    return -1;
  }

  ////////////////////////////////////////////////
  //// a_i, b_i, A, & B Parameters for Form Factor Values
  /////////////////////////////////////////////////*{{{*/
  const double GN_Coef_Fit[4][4] ={
    {-0.24, 10.98, 12.82, 21.97}, /*GEp*/
    {0.12, 10.97, 18.86, 6.55}, /*GMp/mu_p*/
    {1.70, 3.30, -1000, -1000}, /*GEn/GD*/ /*(A & B params)*/
    {2.33, 14.72, 24.20, 84.1} /*GMn/mu_n*/
  };/*}}}*/

  //// Applying parametrization formula
  double tau = kID<3 ? kine::tau(kQ2,"p") : kine::tau(kQ2,"n");
  double numerator = 0.;
  double denominator = 1.;
  if (kID-1!=2) { // for all but GEn
    numerator = 1. + GN_Coef_Fit[kID-1][0]*tau;
    for (int i=1; i<4; i++) denominator += GN_Coef_Fit[kID-1][i]*pow(tau,i);
  } else {
    numerator = GN_Coef_Fit[kID-1][0]*tau;
    denominator += GN_Coef_Fit[kID-1][1]*tau;
  }

  if (kID-1!=2) GNGD_Fit[0] = numerator / (denominator*EMFFFits::GetGDip(kQ2));
  else GNGD_Fit[0] = numerator / denominator;
  GNGD_Err[0] = -1000.;

  return 0;
}

// #######################
// ## Seamus Fit (20??) ##
// #######################
int Seamus20XX::SeamusFit(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err) {
  // This is the fit that is implemented in g4sbs, presumably by Seamus Riordan,
  // for GEn parametrization. In the code it has been indicated as "Our fit"! I
  // couldn't find any corresponding paper related to this parametrization. Asking
  // Andrew and Googling didn't help as well. For the sake of consistency, I 
  // implemented the same parametrization in SIMC generator as well.
  
  // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  if (kID!=3) {
    std::cerr<<"*** ERROR***, Seamus fit only supports kID=3 ie GEn"<<std::endl;
    GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
    return -1;
  }
  
  ////////////////////////////////////////////////
  //// a_i & b_i Parameters for Form Factor Values
  /////////////////////////////////////////////////*{{{*/
  const double GN_Coef_Fit[4][6] ={
    {-1000, -1000, -1000, -1000, -1000, -1000}, /*Not available!*/
    {-1000, -1000, -1000, -1000, -1000, -1000}, /*Not available!*/
    {1.52, 2.629, 3.055, 5.222, 0.04, 11.438}, /*GEn/GD*/
    {-1000, -1000, -1000, -1000, -1000, -1000} /*Not available!*/
  };/*}}}*/

  //// Applying parametrization formula
  double tau = kID<3 ? kine::tau(kQ2,"p") : kine::tau(kQ2,"n");
  double numerator = 0.;
  double denominator = 1.;
  for (int i=0; i<3; i++) numerator += GN_Coef_Fit[kID-1][i]*pow(tau,i+1);
  for (int i=3; i<6; i++) denominator += GN_Coef_Fit[kID-1][i]*pow(tau,i-2);

  GNGD_Fit[0] = numerator / denominator;
  GNGD_Err[0] = -1000.;

  return 0;
  
}

// ###################
// ## Ye Fit (2017) ##
// ###################
int Ye2017::YeFit(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err) {
  // ** Link to the original paper:
  // https://www.sciencedirect.com/science/article/pii/S0370269317309152?via%3Dihub

  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  //// Parameterized Form Factor Central Value and Error
  /////////////////////////////////////////////////////////
  //// ID = 1 for GEp, 2 for GMp, 3 for GEn, 4 for GMn, 
  //// Q2 in GeV^2
  ////
  // The parameterization formula returns the uncertainty devided by G(0)*GD, where 
  //  GD(Q2) = 1./(1+Q2/0.71)^2
  // and GEp(0) = 1, GMp(0) = 2.79284356, GEn(0) = 1, GMn(0) = -1.91304272,
  //
  // The parameterization formula for the Form Factor value is:
  //  $$ GN(z) = sum_{i=0}^{N=12}(a_i * z^i) 
  // Note that the return value has been divided by (G(Q2=0)*G_Dip)
  //
  // The parameterization formula for the Form Factor error is:
  // $$ log_{10}\frac{\delta G}{G_D} = (L+c_0)\Theta_a(L_1-L) 
  //                                 +\sum_{i=1}^{N}(c_i+d_i L)[\Theta_a(L_i-L)-\Theta_a(L_{i+1}-L)]
  //                                 +log_{10}(E_{\inf})\Theta_a(L-L_{N+1})$$
  // where $L=log_{10}(Q^2)$, $\Theta_{a}(x)=[1+10^{-ax}]^{-1}$. $a=1$.

  // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  if (kID<1 || kID>4){
    std::cerr<<"*** ERROR***, kID is not any of [1->GEp, 2->GMp, 3->GEn, 4->GMn]"<<std::endl;
    GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
    return -1;
  }
  ////////////////////////////////////////////////
  //// z-Expansion Parameters for Form Factor Values
  /////////////////////////////////////////////////*{{{*/
  const double GN_Coef_Fit[4][13] ={
    {0.239163298067, -1.10985857441, 1.44438081306, 0.479569465603, -2.28689474187,  1.12663298498, 1.25061984354,-3.63102047159, 4.08221702379,  0.504097346499,  -5.08512046051,  3.96774254395,-0.981529071103}, /*GEp*/
    {0.264142994136, -1.09530612212, 1.21855378178, 0.661136493537, -1.40567892503, -1.35641843888, 1.44702915534, 4.2356697359, -5.33404565341, -2.91630052096,    8.70740306757, -5.70699994375, 1.28081437589}, /*GMp/mu_p*/
    {0.048919981379,-0.064525053912,-0.240825897382,0.392108744873, 0.300445258602,-0.661888687179,-0.175639769687, 0.624691724461,-0.077684299367,-0.236003975259, 0.090401973470, 0.0, 0.0}, /*GEn*/
    {0.257758326959,-1.079540642058, 1.182183812195,0.711015085833,-1.348080936796,-1.662444025208, 2.624354426029, 1.751234494568,-4.922300878888, 3.197892727312,-0.712072389946, 0.0, 0.0} /*GMn/mu_n*/
  };/*}}}*/

    ////////////////////////////////////////////////
    //// Parameters for Form Factor Errors
    ////////////////////////////////////////////////   /*{{{*/
  const double parL[4][2] ={
    {-0.97775297,  0.99685273}, //GEp
    {-0.68452707,  0.99709151}, //GMp
    {-2.02311829, 1.00066282}, //GEn
    {-0.20765505, 0.99767103}, //GMn
  };
  const double parM[4][15] = {
    {  -1.97750308e+00,  -4.46566998e-01,   2.94508717e-01,   1.54467525e+00,
       9.05268347e-01,  -6.00008111e-01,  -1.10732394e+00,  -9.85982716e-02,
       4.63035988e-01,   1.37729116e-01,  -7.82991627e-02,  -3.63056932e-02,
       2.64219326e-03,   3.13261383e-03,   3.89593858e-04}, //GEp

    {  -1.76549673e+00,   1.67218457e-01,  -1.20542733e+00,  -4.72244127e-01,
       1.41548871e+00,   6.61320779e-01,  -8.16422909e-01,  -3.73804477e-01,
       2.62223992e-01,   1.28886639e-01,  -3.90901510e-02,  -2.44995181e-02,
       8.34270064e-04,   1.88226433e-03,   2.43073327e-04}, //GMp

    {  -2.07343771e+00,   1.13218347e+00,   1.03946682e+00,  -2.79708561e-01,
       -3.39166129e-01,   1.98498974e-01,  -1.45403679e-01,  -1.21705930e-01,
       1.14234312e-01,   5.69989513e-02,  -2.33664051e-02,  -1.35740738e-02,
       7.84044667e-04,   1.19890550e-03,   1.55012141e-04}, //GEn

    {  -2.07087611e+00,   4.32385770e-02,  -3.28705077e-01,   5.08142662e-01,
       1.89103676e+00,   1.36784324e-01,  -1.47078994e+00,  -3.54336795e-01,
       4.98368396e-01,   1.77178596e-01,  -7.34859451e-02,  -3.72184066e-02,
       1.97024963e-03,   2.88676628e-03,   3.57964735e-04} //GMn
  };
  const double parH[4][3] = {
    {0.78584754,  1.89052183, -0.4104746}, //GEp
    {0.80374002,  1.98005828, -0.69700928}, //GMp
    {0.4553596,  1.95063341,  0.32421279}, //GEn
    {0.50859057, 1.96863291,  0.2321395} //GMn 
  };
  /*}}}*/

  //// Apply the z-expansion formula for form factor/*{{{*/
  const double tcut = 0.0779191396 ;
  const double t0 = -0.7 ;
  double z = (sqrt(tcut+kQ2)-sqrt(tcut-t0))/(sqrt(tcut+kQ2)+sqrt(tcut-t0)) ;
  double GNQ2 = 0.0;
  for (int i=0;i<13;i++) GNQ2 += GN_Coef_Fit[kID-1][i] * pow(z, i);
  double GDip= pow(1./(1. + kQ2/0.71), 2);
  GNGD_Fit[0] = GNQ2 / GDip; //Note that Coef_Fit have been divided by mu_p or mu_n
  /*}}}*/

  //// Apply the parameterization formula for error/*{{{*/
  double lnQ2 = log10(kQ2);
  double lnGNGD_Err=0.0;
  if (kQ2<1e-3)
    lnGNGD_Err = parL[kID-1][0] + parL[kID-1][1]*lnQ2;
  else if (kQ2>1e2)
    lnGNGD_Err = parH[kID-1][0]*sqrt(lnQ2-parH[kID-1][1]) + parH[kID-1][2];
  else{
    for (int i=0; i<15;i++) lnGNGD_Err += parM[kID-1][i] * pow(lnQ2,i);
  }
  GNGD_Err[0] = pow(10.,(lnGNGD_Err));    //LOG10(dG/G(0)/GD);
  /*}}}*/

  return 0;
}

int Ye2017::YeFitNoTPE(const int kID, const double kQ2, double &FF, double &FF_err) {
  // Only returs GEp, GMp, and the associated errors

  /////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////
  // The proton FFs w/o TPE correction are not included
  // in the official supplymentary materials. The associated 
  // lookup table have been collected from the authors via
  // private communications. 

  // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  if (kID>2) {
    std::cerr<<"*** ERROR***, Ye fit w/o TPE only supports kID=1 & 2 ie GEp & GMp"<<std::endl;
    FF = -1000;  FF_err= -1000;
    return -1;
  }

  // Reading the lookup table
  LookUpTableReader reader;
  std::string filename = "data.csv";
  reader.readCSV(filename);

  // Reading the FFs and the errors
  double GEp = reader.GetClosestValueByKey(kQ2,1);
  double GEpErr = reader.GetClosestValueByKey(kQ2,2);
  double GMp = reader.GetClosestValueByKey(kQ2,3); 
  double GMpErr = reader.GetClosestValueByKey(kQ2,4);;

  // Returning EMFF values
  if (kID==1) {
    FF = GEp; FF_err = GEpErr;
  }
  else {
    FF = GMp; FF_err = GMpErr;
  }

  return 0;
}

// ########################
// ## Christy Fit (2022) ##
// ########################
int Christy2022::ChristyFit(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err) {
  // ** Link to the original paper:
  // https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.102002

  // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  if (kID>2) {
    std::cerr<<"*** ERROR***, Christy fit only supports kID=1 & 2 ie GEp & GMp"<<std::endl;
    GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
    return -1;
  }

  ////////////////////////////////////////////////
  //// a_i, b_i, & c_i Parameters for Form Factor Values
  /////////////////////////////////////////////////*{{{*/
  const double GN_Coef_Fit[2][4] ={
    {0.072, 10.73, 19.81, 4.75}, /*GMp/mu_p*/
    {-0.46, 0.12, -1000, -1000}  /*RS ie (mu_p*GEp/GMp)^2*/ /*c1 and c2*/
  };/*}}}*/

  //// Applying parametrization formula
  double tau = kine::tau(kQ2,"p");
  double numerator = 1.;
  double denominator = 1.;
  // Calculating GMp first
  numerator += GN_Coef_Fit[0][0]*tau;
  for (int i=1; i<4; i++) denominator += GN_Coef_Fit[0][i]*pow(tau,i);
  double GMp_ov_mun = numerator / denominator; 
  // Now calculating RS and then GEp
  double RS = 1.;
  for (int i=0; i<2; i++) RS += GN_Coef_Fit[1][i]*pow(tau,i+1);
  double GEp = sqrt(RS)*GMp_ov_mun;

  // Returning EMFF values
  if (kID-1==1) GNGD_Fit[0] = GMp_ov_mun / EMFFFits::GetGDip(kQ2);
  else GNGD_Fit[0] = GEp / EMFFFits::GetGDip(kQ2);
  GNGD_Err[0] = -1000.;

  return 0;
}

// ##########################
// ## Arrington Fit (2007) ##
// ##########################
int Arrington2007::ArringtonFit(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err) {
  // ** Link to the original paper:
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.76.035205 (TABLE I)

  // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  if (kID>2) {
    std::cerr<<"*** ERROR***, Arrington fit only supports kID=1 & 2 ie GEp & GMp"<<std::endl;
    GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
    return -1;
  }

  ////////////////////////////////////////////////
  //// a_i, b_i, & c_i Parameters for True Form Factor Values (with TPE corrections)
  /////////////////////////////////////////////////*{{{*/
  const double GN_Coef_Fit[2][8] ={
    {-1.465, 1.260, 0.262, 9.627, 0.000, 0.000, 11.179, 13.245}, /*GMp/mu_p*/
    {3.439, -1.602, 0.068, 15.055, 48.061, 99.304, 0.012, 8.650}  /*GEp*/
  };/*}}}*/

  //// Applying parametrization formula
  double tau = kine::tau(kQ2,"p");
  double numerator = 1.;
  double denominator = 1.;

  //// Implementing the parametrization
  if (kID-1==1) { // GMp/mu_p
    for (int i=0; i<3; i++) numerator += GN_Coef_Fit[0][i]*pow(tau,i+1);
    for (int i=3; i<8; i++) denominator += GN_Coef_Fit[0][i]*pow(tau,i-2);
    double GMp_ov_mup = numerator/denominator;
    GNGD_Fit[0] = GMp_ov_mup / EMFFFits::GetGDip(kQ2);
  } else { // GEp
    for (int i=0; i<3; i++) numerator += GN_Coef_Fit[1][i]*pow(tau,i+1);
    for (int i=3; i<8; i++) denominator += GN_Coef_Fit[1][i]*pow(tau,i-2);
    double GEp = numerator/denominator;
    GNGD_Fit[0] = GEp / EMFFFits::GetGDip(kQ2);    
  }
  GNGD_Err[0] = -1000.;

  return 0;
}

int Arrington2007::ArringtonFitNoTPE(const int kID, const double kQ2, double *GNGD_Fit, double* GNGD_Err) {
  // ** Link to the original paper:
  // https://journals.aps.org/prc/pdf/10.1103/PhysRevC.76.035205 (TABLE IV)

  // GEp->kID=1, GMp->kID=2, GEn->kID=3, GMn->kID=4
  if (kID>2) {
    std::cerr<<"*** ERROR***, Arrington fit only supports kID=1 & 2 ie GEp & GMp"<<std::endl;
    GNGD_Fit[0] = -1000;  GNGD_Err[0] = -1000;
    return -1;
  }

  ////////////////////////////////////////////////
  //// a_i, b_i, & c_i Parameters for Effective Form Factor Values (No TPE Correction)
  /////////////////////////////////////////////////*{{{*/
  const double GN_Coef_Fit[2][8] ={
    {-2.151, 4.261, 0.159, 8.647, 0.001, 5.245, 82.817, 14.191}, /*Fm/mu_p -> GMp/mu_p in OPE*/
    {-1.651, 1.287, -0.185, 9.531, 0.591, 0.000, 0.000, 4.994}  /*Fe -> GEp in OPE*/
  };/*}}}*/

  //// Applying parametrization formula
  double tau = kine::tau(kQ2,"p");
  double numerator = 1.;
  double denominator = 1.;

  //// Implementing the parametrization
  if (kID-1==1) { // GMp/mu_p
    for (int i=0; i<3; i++) numerator += GN_Coef_Fit[0][i]*pow(tau,i+1);
    for (int i=3; i<8; i++) denominator += GN_Coef_Fit[0][i]*pow(tau,i-2);
    double GMp_ov_mup = numerator/denominator;
    GNGD_Fit[0] = GMp_ov_mup / EMFFFits::GetGDip(kQ2);
  } else { // GEp
    for (int i=0; i<3; i++) numerator += GN_Coef_Fit[1][i]*pow(tau,i+1);
    for (int i=3; i<8; i++) denominator += GN_Coef_Fit[1][i]*pow(tau,i-2);
    double GEp = numerator/denominator;
    GNGD_Fit[0] = GEp / EMFFFits::GetGDip(kQ2);    
  }
  GNGD_Err[0] = -1000.;

  return 0;
}
