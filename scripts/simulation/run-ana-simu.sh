#!/bin/bash

#SBATCH --mem-per-cpu=1000

export GMN_ANA=/u/home/pdbforce/gmn_ana

process=$1  # elas/qelas
conf=$2
field=$3
gen=$4      # simc/g4sbs
reaction=$5 # heep/deeN/inel

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
    analyzer -b -q ''$process'_ana_simu.cpp("sbs'$conf'-sbs'$field'p-simu/conf_'$process'_ana_simu.json",2,"'$gen'","'$reaction'")'
else
    analyzer -b -q ''$process'_ana_simu.cpp("sbs'$conf'-sbs'$field'p-simu/conf_'$process'_ana_simu.json",1,"'$gen'","'$reaction'")'
fi

#sbatch --mem-per-cpu=2000 run-ana-simu.sh elas 4 30 simc heep
#sbatch --mem-per-cpu=2000 run-ana-simu.sh qelas 4 30 simc deeN
