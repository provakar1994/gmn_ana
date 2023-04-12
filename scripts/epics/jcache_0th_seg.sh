#!/bin/sh

## Edits
# P. Datta <pdbforce@jlab.org> Created 04-12-2023

# This script jcache the 0th segment of evio files
# of a given GMn run.

# Define arguments
runnum=$1 #run number

DATA_DIR=/mss/halla/sbs/raw
CACHE_DIR=/cache/halla/sbs/raw

# check whether the file is already in CACHE_DIR
if [[ ! -f "$CACHE_DIR/e1209019_"$runnum".evio.0.0" ]]; then
    jcache get "$DATA_DIR/e1209019_"$runnum".evio.0.0"
else
    echo "0th segment of run $runnum is already in CACHE_DIR"
    break;
fi
