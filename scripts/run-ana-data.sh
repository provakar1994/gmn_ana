#!/bin/bash

# Location of this repository and of the installed rootlogon.C (loads the gmnana library).
# Set GMN_ANA / ROOTLOGON_ENV in your environment before submitting (sbatch passes the
# environment to the job); the defaults below are the original author's ifarm setup.
export GMN_ANA=${GMN_ANA:-/u/home/pdbforce/gmn_ana}
ROOTLOGON_ENV=${ROOTLOGON_ENV:-/home/pdbforce/install_el9/runenv/rootlogon.C}
if [[ ! -d $GMN_ANA/scripts || ! -r $ROOTLOGON_ENV ]]; then
    echo "ERROR: GMN_ANA ($GMN_ANA) or ROOTLOGON_ENV ($ROOTLOGON_ENV) not found" >&2
    exit 1
fi

process=$1
conf=$2
field=$3
pass=$4
model=$5
target=$6 #LD2/Dummy

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

cd $GMN_ANA/scripts

if [[ "$process" == "qelas" ]]; then
    analyzer -b -q $ROOTLOGON_ENV ''$process'_ana_data.cpp("configs/sbs'$conf'-sbs'$field'p/conf_'$process'_ana_data.json","'$target'",'$pass','$model')'
else
    analyzer -b -q $ROOTLOGON_ENV ''$process'_ana_data.cpp("configs/sbs'$conf'-sbs'$field'p/conf_'$process'_ana_data.json",'$pass','$model')'
fi

#sbatch --mem-per-cpu=2000 run-ana-data.sh elas 11 100 pass model
#sbatch --mem-per-cpu=2000 run-ana-data.sh qelas 11 100 pass model LD2
