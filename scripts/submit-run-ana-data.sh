#!/bin/bash

process=$1
conf=$2
field=$3
pass=$4
model=$5
target=$6

# build readable job name
if [[ "$process" == "elas" ]]; then
    target=LH2
fi
jobname=data_${process}_sbs${conf}_f${field}p_${target}_pass${pass}_model${model}

# optional: log directory
logdir=slurm_out
#mkdir -p "$logdir"

sbatch \
  --job-name="$jobname" \
  --output=${logdir}/%x_%j.out \
  --error=${logdir}/%x_%j.err \
  --mem-per-cpu=2000 \
  run-ana-data.sh \
  "$process" "$conf" "$field" "$pass" "$model" "$target"

#Execution: ./submit-run-ana-data.sh elas conf field pass model
#Example:
#./submit-run-ana-data.sh elas 4 30 3 1
#./submit-run-ana-data.sh qelas 4 30 3 2 LD2

