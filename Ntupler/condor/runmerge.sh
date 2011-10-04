#!/bin/bash

config=htt.config
if [ "` hostname | grep '\.mit\.edu'`" ]; then
    ntupledir=/scratch/$USER/htt/nojetid
else
    ntupledir=/data/blue/$USER/htt/nojetid
fi

for dataset in `cat $config | grep -v ^# | tr -s ' ' | cut -d' ' -f 1`; do
  if [ -s $ntupledir/${dataset}_ntuple.root ] || [ -s $ntupledir/${dataset}_ntuple.json ]; then
    echo "!!--->merged files ${dataset}_ntuple.{root,json} already exists."
  else
    echo "merging..."
    echo $ntupledir/${dataset}_ntuple.root    > merge.txt
    ls $ntupledir/${dataset}_????_ntuple.root >>merge.txt
    eval `scramv1 runtime -sh`
    root -b -l -q ../macros/MergeNtuples.C+\(\"merge.txt\"\)
    if [ "`echo $dataset | grep 'r11\|r10'`" ]; then
	./mergeJSON.sh merge.txt
    fi
    rm merge.txt
  fi
done


