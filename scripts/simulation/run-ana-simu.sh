#!/bin/bash

#SBATCH --mem-per-cpu=1000

# Location of this repository and of the installed rootlogon.C (loads the gmnana library).
# Set GMN_ANA / ROOTLOGON_ENV in your environment before submitting (sbatch passes the
# environment to the job); the defaults below are the original author's ifarm setup.
export GMN_ANA=${GMN_ANA:-/u/home/pdbforce/gmn_ana}
ROOTLOGON_ENV=${ROOTLOGON_ENV:-/home/pdbforce/install_el9/runenv/rootlogon.C}
if [[ ! -d $GMN_ANA/scripts || ! -r $ROOTLOGON_ENV ]]; then
    echo "ERROR: GMN_ANA ($GMN_ANA) or ROOTLOGON_ENV ($ROOTLOGON_ENV) not found" >&2
    exit 1
fi

process=$1  # elas/qelas
conf=$2
field=$3
gen=$4      # simc/g4sbs
reaction=$5 # heep/deeN/inel
variant=$6  # optional: config block suffix, e.g. noshield (default block if empty)

MODULES=/etc/profile.d/modules.sh 

if [[ $(type -t module) != function && -r ${MODULES} ]]; then 
source ${MODULES} 
fi 

if [ -d /apps/modulefiles ]; then 
module use /apps/modulefiles 
fi 

module use /group/halla/modulefiles
#module load geant4/11.1.2
module load analyzer/1.7.12
#module load cmake/3.23.2

cd $GMN_ANA/scripts/simulation

if [[ "$process" == "qelas" ]]; then
    analyzer -b -q $ROOTLOGON_ENV ''$process'_ana_simu.cpp("configs/sbs'$conf'-sbs'$field'p-simu/conf_'$process'_ana_simu.json",2,"'$gen'","'$reaction'","'$variant'")'
else
    analyzer -b -q $ROOTLOGON_ENV ''$process'_ana_simu.cpp("configs/sbs'$conf'-sbs'$field'p-simu/conf_'$process'_ana_simu.json",1,"'$gen'","'$reaction'","'$variant'")'
fi

#sbatch --mem-per-cpu=2000 run-ana-simu.sh elas 4 30 simc heep
#sbatch --mem-per-cpu=2000 run-ana-simu.sh qelas 4 30 simc deeN
#sbatch --mem-per-cpu=2000 run-ana-simu.sh qelas 11 100 simc deeN noshield
