#!/bin/bash
   runMacro=$1
  outputDir=$2
  clusterid=$3
    dataset=$4
     infile=$5
    nevents=$6
         id=$7
    workdir=$8
    channel=$9

workDir=$CMSSW_BASE #`pwd`
scramDir=$CMSSW_BASE #`pwd`
echo `hostname`
echo "workDir: $workDir"
echo "workDir: $scramDir"
echo "args:    $*"


cd $scramDir/src/
eval `scramv1 runtime -sh`
cd $workDir

root -l -b -q apply.C+\(${id},${channel},\"${dataset}\",\"${outputDir}/${infile}\",1000\)

status=`echo $?`
echo "Status - $status"

mkdir                                  ${outputDir}/${dataset}
mv ${dataset}_${id}_ntuple.root        ${outputDir}/${dataset}

exit $status
