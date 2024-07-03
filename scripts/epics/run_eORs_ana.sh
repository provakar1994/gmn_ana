#!/bin/bash

###################################################
# This script will run EPICS/Scaler data analysis 
# ------
# P. Datta CREATED 03/19/2023
###################################################

anamode=$1
target=$2
conf=$3
sbsmagfield=$4
nruns=$5
verbose=$6

# Validating arguments
if [[ "$#" -ne 6 ]]; then
    echo -e "\n--!-- Illegal number of arguments!"
    echo -e "This script expects 6 arguments:" 
    echo -e " ---- "
    echo -e " <anamode>      : e=>EPICS, s=>Scaler"
    echo -e " <target>       : LH2/LD2/Optics/Dummy"
    echo -e " <conf>         : SBS config, -1=>All"
    echo -e " <sbsmagfield>  : -1=>No restriction"
    echo -e " <nruns>        : # runs to analyze"
    echo -e " <verbose>      : >0=>Debug"
    echo -e " ---- "
    echo -e "E.g.: ./run_eORs_ana.sh s LH2 4 30 2 1 \n"
    exit;
fi

# List of ROOT file directories 
## -- Pass 0/1
# sbs4dir='/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS4/'$target''
# sbs7dir='/w/halla-scshelf2102/sbs/sbs-gmn/pass0/SBS7/'$target''
# sbs11dir='/lustre19/expphy/volatile/halla/sbs/sbs-gmn/GMN_REPLAYS/pass1/SBS11/'$target''
# sbs14dir='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS14/'$target''
# sbs8dir='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS8/'$target''
# sbs9dir='/w/halla-scshelf2102/sbs/sbs-gmn/pass1/SBS9/'$target''
## -- Pass 2
sbs4dir='/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS4/'$target''
sbs7dir='/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS7/'$target''
sbs11dir='/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS11/'$target''
sbs14dir='/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS14/'$target''
sbs8dir='/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS8/'$target''
sbs9dir='/w/halla-scshelf2102/sbs/pdbforce/gmn-data/pass2/SBS9/'$target''

# Defining arrays for combined analysis
confarr=(4 7 11 14 8 9)
rfdirarr=($sbs4dir $sbs7dir $sbs11dir $sbs14dir $sbs8dir $sbs9dir)

# Validate SBS config value
set +e
invalidConf() {
    local e match="$1"
    shift
    for e; do [[ "$e" == "$match" ]] && return 0; done
    return 1
}

# Executing the analysis
if [[ $anamode = 's' ]]; then
    if [[ $conf -eq -1 ]]; then
	index=0
	sbsmagfield=-1
	for i in ${confarr[@]}; do
	    analyzer -l -q 'get_scalerdata_prun.cpp("'$target'",'$i','$sbsmagfield','$nruns',"'${rfdirarr[$index]}'",'$verbose')'
	    index=$((index+1))
	done
    else
	invalidConf "$conf" "${confarr[@]}"
	if [[ $? -eq 0 ]]; then
	    index=-1
	    if [[ $conf -eq 4 ]]; then
		index=0
	    elif [[ $conf -eq 7 ]]; then
		index=1
	    elif [[ $conf -eq 11 ]]; then
		index=2
	    elif [[ $conf -eq 14 ]]; then
		index=3
	    elif [[ $conf -eq 8 ]]; then
		index=4
	    else
		index=5
	    fi
    	    analyzer -l -q 'get_scalerdata_prun.cpp("'$target'",'$conf','$sbsmagfield','$nruns',"'${rfdirarr[$index]}'",'$verbose')'
	else 
    	    echo ${confarr[2]}
    	    echo -e "\n--!-- [ERROR]\n Invalid SBS configuration! Expected: -1/4/7/11/14/8/9\n"
    	    exit;
	fi
    fi
fi

