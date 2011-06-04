#!/bin/bash

config=hypha.config
ntupledir=/scratch/$USER/htt

for dataset in `cat $config | grep -v ^# | tr -s ' ' | cut -d' ' -f 1`; do
  if [ -s $ntupledir/${dataset}_ntuple.root ]; then
    echo "!!--->merged file ${dataset}_ntuple.root already exists."
  else
    echo "merging..."
    echo $ntupledir/${dataset}_ntuple.root   > merge.txt
    ls $ntupledir/${dataset}_????_ntuple.root >>merge.txt
    eval `scramv1 runtime -sh`
    root -b -q ../macros/MergeNtuples.C+\(\"merge.txt\"\)
    rm merge.txt
  fi
done


