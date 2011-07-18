#!/bin/bash

config=hypha.config
if [ "` hostname | grep '\.mit\.edu'`" ]; then
    ntupledir=/scratch/$USER/htt
else
    ntupledir=/data/blue/$USER/htt
fi

for dataset in `cat $config | grep -v ^# | tr -s ' ' | cut -d' ' -f 1`; do

  # remove .root files
  if [ -s $ntupledir/${dataset}_ntuple.root ]; then
    echo "merged .root file exists for $dataset"
    for file in `ls $ntupledir/${dataset}_????_ntuple.root`; do
      echo "          rm $file"
      rm $file
    done
  else
    echo "!!-->no merged root file for $dataset"
  fi

  # remove .json files
  if [ -s $ntupledir/${dataset}_ntuple.json ]; then
   echo "merged .json file exists for $dataset"
   for file in `ls $ntupledir/${dataset}_????_ntuple.json`; do
     echo "          rm $file"
     rm $file
   done
  else
    echo "!!-->no merged .json file for $dataset"
  fi

done
