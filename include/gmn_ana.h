/* 
   gmn_ana.h will include all the gmn-ana libraries in the most
   efficient way. Just including this file in an analysis script
   will be enough to get access to all the gmn_ana libraries.
   -----
   P. Datta <pdbforce@jlab.org> Created 09-17-2022
*/

#include "Cut.h"
#include "Fit.h"
#include "FitFns.h"
#include "CodaRun.h"
#include "SimuJob.h"
#include "EMFFFits.h"
#include "Constants.h"
#include "Utilities.h"
#include "SetROOTVar.h"
#include "KinematicVar.h"
#include "ExpConstants.h"

/*
#include "Constants.h"                // namespace constant (General constants) 
#include "CodaRun.h"                  // struct CodaRun      
#include "SimuJob.h"                  // struct SimuJob  
#include "FitFns.h"                   // class FitFn        (Various fit fns.) 
#include "../src/Fit.cpp"             // namespace fit      (Various fit methods.) 
#include "../src/Cut.cpp"             // namespace cut      (Various cut defn.) 
#include "../src/Utilities.cpp"       // namespace util_pd  (Various utility fns.) 
#include "../src/ExpConstants.cpp"    // namespace expconst & class SBSconfig (Experimental constants) 
#include "../src/SetROOTVar.cpp"      // namespace setrootvar  
#include "../src/KinematicVar.cpp"    // namespace kine     (fns. to calculate physics variables) 
*/

/* --- List of gmn_ana libraries --- */
// Fit.h          : namespace fit
// Cut.h          : namespace cut
// Constants.h    : namespace constant
// KinematicVar.h : namespace kine
// ExpConstants.h : namespace expconst & class SBSconfig
// SetROOTVar.h   : namespace setrootvar
// Utilities.h    : namespace util_pd
