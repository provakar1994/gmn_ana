#include "../include/ExpConstants.h"

/* 436 => GEN 8-) */

namespace expconst {

  double ebeam(int config) {
    /* These GMn energy values are for LD2 target from good runlist. For LH2 runs
       these numbers are ever so slightly different! 
    */
    if(config==1)
      return 1.916;
    else if(config==4)
      return 3.7393;
    else if(config==7)
      return 7.9308;
    else if(config==11)
      return 9.889;
    else if(config==14)
      return 5.9828;
    else if(config==8)
      return 5.9826;
    else if(config==9)
      return 4.0268;
    else if(config==4363)
      return 6.373;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double pcentral(int config){
    /* Returns central scattered electron energy per config */
    if(config==1)
      return 1.09;
    else if(config==4)
      return 2.12;
    else if(config==7)
      return 2.66;
    else if(config==11)
      return 2.67;
    else if(config==14)
      return 2.0;
    else if(config==8)
      return 3.58;
    else if(config==9)
      return 1.63;
    else if(config==4363)
      return 2.73;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double bbtheta(int config){
    if(config==1)
      return 51.0;
    else if(config==4)
      return 36.0;
    else if(config==7)
      return 40.0;
    else if(config==11)
      return 42.0;
    else if(config==14)
      return 46.5;
    else if(config==8)
      return 26.5;
    else if(config==9)
      return 49.0;
    else if(config==4363)
      return 36.5;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double bbdist(int config){
    if(config==1)
      return 1.8518;
    else if(config==4)
      return 1.7988;
    else if(config==7)
      return 1.84896;
    else if(config==11)
      return 1.55146;
    else if(config==14)
      return 1.84787;
    else if(config==8)
      return 1.97473;
    else if(config==9)
      return 1.550;
    else if(config==4363)
      return 1.63;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double sbstheta(int config){
    if(config==1)
      return 33.5;
    else if(config==4)
      return 31.9;
    else if(config==7)
      return 16.1;
    else if(config==11)
      return 13.3;
    else if(config==14)
      return 17.3;
    else if(config==8)
      return 29.9;
    else if(config==9)
      return 22.5;
    else if(config==4363)
      return 22.1;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double sbsdist(int config){
    if(config==1||config==4||config==7||config==11
       ||config==14||config==8||config==9)
      return 2.25;
    else if(config==4363)
      return 2.8;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double hcaldist(int config){
    if(config==1)
      return 13.5;
    else if(config==4||config==8||config==9)
      return 11.0;
    else if(config==7||config==14)
      return 14.0;
    else if(config==11)
      return 14.5;
    else if(config==4363)
      return 17.0;
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }
  //--------------------------------------------
  double hcaltheta(int config){
    if(config==1)
      return 33.5;
    else if(config==4)
      return 31.9;
    else if(config==7)
      return 16.1;
    else if(config==11)
      return 13.3;
    else if(config==14)
      return 17.3;
    else if(config==8)
      return 29.4;
    else if(config==9)
      return 22.0;
    else if(config==4363)
      return 21.6; 
    else{
      std::cerr << "Enter a valid SBS configuration!" << std::endl;
      return -1;
    }
  }

  /* ###########################################################
     ##   Functions to get collision stopping powers (dEdx)   ##  
     ########################################################### */
  /* Returns collision stopping power before and after scattering based on SBS configuration. Before scattering
     dE/dx is calculated based on beam energy whereas after scattering dE/dx is calculated based on the energy of
     central elastically scattered e-. See Utilities::GetElossInTgt for more info.
    *. Stopping power: 
       - Useful website by NIST ESTAR: https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
       - Above website don't have data on deuterium but I found a master's thesis from 1952 that shows that the 
         stopping power of deuterium is pretty close to hydrogen. https://open.library.ubc.ca/media/stream/pdf/831/1.0085416/1
  */

  //--------------------------------------------
  double GetdEdxCollH(int const config,                // SBS config
		      bool const is_before_scattering) // 1=>YES, 0=>NO (i.e. after scattering)
  /* Collisional stopping power of hydrogen/deuterium. Read notes above for more info. */
  {
    if (is_before_scattering) {
      std::unordered_map<int,double> m = {{1,  5.713E-3},
					  {4,  5.833E-3},
					  {7,  5.947E-3},
					  {11, 5.981E-3},
					  {14, 5.904E-3},
					  {8,  5.904E-3},
					  {9,  5.844E-3}};
      if (m.find(config)==m.end()) 
	throw std::invalid_argument("[ExpConstants::GetdEdxCollH] ERROR! Invalid SBS config!!");
      return m[config];
      
    }else {
      std::unordered_map<int,double> m = {{1,  5.645E-3},
					  {4,  5.746E-3},
					  {7,  5.781E-3},
					  {11, 5.781E-3},
					  {14, 5.737E-3},
					  {8,  5.826E-3},
					  {9,  5.706E-3}};
      if (m.find(config)==m.end()) 
	throw std::invalid_argument("[ExpConstants::GetdEdxCollH] ERROR! Invalid SBS config!!");
      return m[config];
    }
  }
  //--------------------------------------------
  double GetdEdxCollAl(int const config,                // SBS config
		       bool const is_before_scattering) // 1=>YES, 0=>NO (i.e. after scattering)
  /* Collisional stopping power of Al */
  {
    if (is_before_scattering) {
      std::unordered_map<int,double> m = {{1,  2.069E-3},
					  {4,  2.118E-3},
					  {7,  2.174E-3},
					  {11, 2.190E-3},
					  {14, 2.153E-3},
					  {8,  2.153E-3},
					  {9,  2.124E-3}};
      if (m.find(config)==m.end()) 
	throw std::invalid_argument("[ExpConstants::GetdEdxCollAl] ERROR! Invalid SBS config!!");
      return m[config];
      
    }else {
      std::unordered_map<int,double> m = {{1,  2.027E-3},
					  {4,  2.076E-3},
					  {7,  2.093E-3},
					  {11, 2.093E-3},
					  {14, 2.072E-3},
					  {8,  2.115E-3},
					  {9,  2.057E-3}};
      if (m.find(config)==m.end()) 
	throw std::invalid_argument("[ExpConstants::GetdEdxCollAl] ERROR! Invalid SBS config!!");
      return m[config];
    }
  }
  //--------------------------------------------
  double GetdEdxCollPE(int const config,                // SBS config
		       bool const is_before_scattering) // 1=>YES, 0=>NO (i.e. after scattering)
  /* Collisional stopping power of PE */
  {
    if (is_before_scattering) {
      std::unordered_map<int,double> m = {{1,  2.526E-3},
					  {4,  2.585E-3},
					  {7,  2.651E-3},
					  {11, 2.670E-3},
					  {14, 2.626E-3},
					  {8,  2.626E-3},
					  {9,  2.591E-3}};
      if (m.find(config)==m.end()) 
	throw std::invalid_argument("[ExpConstants::GetdEdxCollPE] ERROR! Invalid SBS config!!");
      return m[config];
      
    }else {
      std::unordered_map<int,double> m = {{1,  2.477E-3},
					  {4,  2.535E-3},
					  {7,  2.555E-3},
					  {11, 2.555E-3},
					  {14, 2.530E-3},
					  {8,  2.581E-3},
					  {9,  2.512E-3}};
      if (m.find(config)==m.end()) 
	throw std::invalid_argument("[ExpConstants::GetdEdxCollPE] ERROR! Invalid SBS config!!");
      return m[config];
    }
  }


  // //--------------------------------------------
  // double dEdx_coll_bs(int config, std::string material)
  // /* Returns collision stopping power before scattering based on SBS configuration. Calculating
  //    dEdx based on beam energy (mimics before scattering situation).
  //   *. Stopping power: 
  //      - Useful website by NIST ESTAR: https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
  //      - Above website don't have data on deuterium but I found a master's thesis from 1952 that shows that the 
  //        stopping power of deuterium is pretty close to hydrogen. https://open.library.ubc.ca/media/stream/pdf/831/1.0085416/1
  // */
  // {
  //   if (material.compare("Al")==0) {
  //     std::unordered_map<int,double> m = {{1,  2.069E-3},
  // 					  {4,  2.118E-3},
  // 					  {7,  2.174E-3},
  // 					  {11, 2.190E-3},
  // 					  {14, 2.153E-3},
  // 					  {8,  2.153E-3},
  // 					  {9,  2.124E-3}};
  //     if (m.find(config)==m.end()) 
  // 	throw std::invalid_argument("[ExpConstants::dEdx_coll_bs] ERROR Invalid SBS config!! Valid options are: 4,7,11,14,8,9");
  //     return m[config];
      
  //   }else if (material.compare("H")==0) {
  //     std::unordered_map<int,double> m = {{1,  5.713E-3},
  // 					  {4,  5.833E-3},
  // 					  {7,  5.947E-3},
  // 					  {11, 5.981E-3},
  // 					  {14, 5.904E-3},
  // 					  {8,  5.904E-3},
  // 					  {9,  5.844E-3}};
  //     if (m.find(config)==m.end()) 
  // 	throw std::invalid_argument("[ExpConstants::dEdx_coll_bs] ERROR Invalid SBS config!! Valid options are: 4,7,11,14,8,9");
  //     return m[config];

  //   }else if (material.compare("PE")==0) { // Polyethylene
  //     std::unordered_map<int,double> m = {{1,  2.526E-3},
  // 					  {4,  2.585E-3},
  // 					  {7,  2.651E-3},
  // 					  {11, 2.670E-3},
  // 					  {14, 2.626E-3},
  // 					  {8,  2.626E-3},
  // 					  {9,  2.591E-3}};
  //     if (m.find(config)==m.end()) 
  // 	throw std::invalid_argument("[ExpConstants::dEdx_coll_bs] ERROR Invalid SBS config!! Valid options are: 4,7,11,14,8,9");
  //     return m[config];

  //   }else {
  //     throw std::invalid_argument("[ExpConstants::dEdx_coll_bs] Enter a valid material type!");
  //   } 
  // }
  // //--------------------------------------------
  // double dEdx_coll_as(int config, std::string material)
  // /* Returns collision stopping power after scattering based on SBS configuration. Calculating
  //    dEdx based on the central scattered elastic e- energy (mimics before scattering situation).
  //   *. Stopping power: 
  //      - Useful website by NIST ESTAR: https://physics.nist.gov/PhysRefData/Star/Text/ESTAR.html
  //      - Above website don't have data on deuterium but I found a master's thesis from 1952 that shows that the 
  //        stopping power of deuterium is pretty close to hydrogen. https://open.library.ubc.ca/media/stream/pdf/831/1.0085416/1
  // */
  // {
  //   if (material.compare("Al")==0) {
  //     std::unordered_map<int,double> m = {{1,  2.027E-3},
  // 					  {4,  2.076E-3},
  // 					  {7,  2.093E-3},
  // 					  {11, 2.093E-3},
  // 					  {14, 2.072E-3},
  // 					  {8,  2.115E-3},
  // 					  {9,  2.057E-3}};
  //     if (m.find(config)==m.end()) 
  // 	throw std::invalid_argument("[ExpConstants::dEdx_coll_as] ERROR Invalid SBS config!! Valid options are: 4,7,11,14,8,9");
  //     return m[config];

  //   }else if (material.compare("H")==0) { // Hydrogen/Deuterium
  //     std::unordered_map<int,double> m = {{1,  5.645E-3},
  // 					  {4,  5.746E-3},
  // 					  {7,  5.781E-3},
  // 					  {11, 5.781E-3},
  // 					  {14, 5.737E-3},
  // 					  {8,  5.826E-3},
  // 					  {9,  5.706E-3}};
  //     if (m.find(config)==m.end()) 
  // 	throw std::invalid_argument("[ExpConstants::dEdx_coll_as] ERROR Invalid SBS config!! Valid options are: 4,7,11,14,8,9");
  //     return m[config];

  //   }else if (material.compare("PE")==0) { // Polyethylene
  //     std::unordered_map<int,double> m = {{1,  2.477E-3},
  // 					  {4,  2.535E-3},
  // 					  {7,  2.555E-3},
  // 					  {11, 2.555E-3},
  // 					  {14, 2.530E-3},
  // 					  {8,  2.581E-3},
  // 					  {9,  2.512E-3}};
  //     if (m.find(config)==m.end()) 
  // 	throw std::invalid_argument("[ExpConstants::dEdx_coll_as] ERROR Invalid SBS config!! Valid options are: 4,7,11,14,8,9");
  //     return m[config];

  //   }else
  //     throw std::invalid_argument("[ExpConstants::dEdx_coll_as] Enter a valid material type!");
  // }

} //::expconst


