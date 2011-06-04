#!/bin/bash

config=hypha.config
ntupledir=/scratch/$USER/htt

for dataset in `cat $config | grep -v ^# | tr -s ' ' | cut -d' ' -f 1`; do
  if [ -s $ntupledir/${dataset}_ntuple.root ]; then
    echo "merged file exists for $dataset"
    for file in `ls $ntupledir/${dataset}_????_ntuple.root`; do
      echo "          rm $file"
      rm $file
    done
  else
    echo "!!-->no merged file for $dataset"
  fi
done
