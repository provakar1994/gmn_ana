#!/bin/bash

export GMN_ANA=/u/home/pdbforce/gmn_ana

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
module load geant4/11.1.2
module load analyzer/1.7.12
#module load cmake/3.23.2

cd $GMN_ANA/scripts

if [[ "$process" == "qelas" ]]; then
    analyzer -b -q ''$process'_ana_data.cpp("sbs'$conf'-sbs'$field'p/conf_'$process'_ana_data.json","'$target'",'$pass','$model')'
else
    analyzer -b -q ''$process'_ana_data.cpp("sbs'$conf'-sbs'$field'p/conf_'$process'_ana_data.json",'$pass','$model')'
fi

#sbatch --mem-per-cpu=2000 run-ana-data.sh elas 11 100 pass model
#sbatch --mem-per-cpu=2000 run-ana-data.sh qelas 11 100 pass model LD2
