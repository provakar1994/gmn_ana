#!/bin/bash

process=$1  # elas/qelas
conf=$2
field=$3
gen=$4      # simc/g4sbs
reaction=$5 # heep/deeN/inel
variant=$6  # optional: config block suffix, e.g. noshield

jobname=simu_${process}_sbs${conf}_f${field}p_${gen}_${reaction}${variant:+_$variant}

# log directory (sbatch fails if it does not exist)
logdir=slurm_out
mkdir -p "$logdir"

sbatch \
  --job-name="$jobname" \
  --output=${logdir}/%x_%j.out \
  --error=${logdir}/%x_%j.err \
  --mem-per-cpu=2000 \
  run-ana-simu.sh \
  "$process" "$conf" "$field" "$gen" "$reaction" "$variant"

#Execution: ./submit-run-ana-simu.sh elas/qelas conf field gen reaction [variant]
#Example:
#./submit-run-ana-simu.sh elas 4 30 simc heep
#./submit-run-ana-simu.sh qelas 4 30 simc deeN
#./submit-run-ana-simu.sh qelas 4 30 g4sbs inel
#./submit-run-ana-simu.sh qelas 11 100 simc deeN noshield   # non-default block simc_deeN_model2_noshield



