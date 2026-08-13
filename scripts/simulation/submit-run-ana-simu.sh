#!/bin/bash

process=$1  # elas/qelas
conf=$2
field=$3
gen=$4      # simc/g4sbs
reaction=$5 # heep/deeN/inel

jobname=simu_${process}_sbs${conf}_f${field}p_${gen}_${reaction}

# optional: log directory
logdir=slurm_out
#mkdir -p "$logdir"

sbatch \
  --job-name="$jobname" \
  --output=${logdir}/%x_%j.out \
  --error=${logdir}/%x_%j.err \
  --mem-per-cpu=2000 \
  run-ana-simu.sh \
  "$process" "$conf" "$field" "$gen" "$reaction"

#Execution: ./submit-run-ana-simu.sh elas conf field pass model
#Example:
#./submit-run-ana-simu.sh elas 4 30 simc heep
#./submit-run-ana-simu.sh qelas 4 30 simc deeN
#./submit-run-ana-simu.sh qelas 4 30 g4sbs inel



