#!/bin/bash

skimname=emu
config=htt.config
if [ "` hostname | grep '\.mit\.edu'`" ]; then
    ntupledir=/scratch/$USER/htt
else
    ntupledir=/data/blue/$USER/htt
fi

for dataset in `cat $config | grep -v ^# | tr -s ' ' | cut -d' ' -f 1`; do
    line=`grep -v ^# $config | grep $dataset`
    usegen=`echo $line 		| tr -s ' ' | cut -d ' ' -f 5`
    if [[ $usegen == 0 ]]; then
	sed -i 's@.*#define _USEGEN_@//#define _USEGEN_@' ../macros/SkimNtuples.C
    else
	sed -i 's@.*#define _USEGEN_@#define _USEGEN_@'   ../macros/SkimNtuples.C
    fi
#  if [ -s $ntupledir/${dataset}_ntuple.root ]; then
#    echo "!!--->merged file ${dataset}_ntuple.root already exists."
#  else
    echo $ntupledir/${dataset}_${skimname}_skim.root    > skim.input
    echo $ntupledir/${dataset}_ntuple.root              >>skim.input
    eval `scramv1 runtime -sh`
    root -b -l -q ../macros/SkimNtuples.C+\(\"skim.input\"\)
    rm skim.input
#  fi
done

exit 0
