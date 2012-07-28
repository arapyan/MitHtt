#!/bin/bash

  outputDir=$1
   runMacro=$2
    dataset=$3 
     infile=$4
    nevents=$5
    channel=$6

workDir=$CMSSW_BASE #`pwd`            # should be run from $src/xxx/condor
script=runint.sh

clusterstr=`ls -ltr /tmp/ | grep krb | grep \`id -u\` | sed "s@_@ @g" | tail -1 | awk '{print $11}'`

soFile=`echo $runMacro | sed 's/\./_/'`.so

# check a few things
if [ ! "$CMSSW_BASE" ]; then
  echo "-------> error: define cms environment."
  exit 1
fi
if [ "$runMacro" -nt "$soFile" ]; then
  echo "-------> error: forgot to recompile run macro."
  exit 1
fi

# grab the files needed for the job
cp rootlogon.C             $workDir
cp $runMacro               $workDir
cp $soFile                 $workDir

# set up kerberos credentials
if ! klist -5 &> /dev/null; then echo "-------> warning: check credentials cache."; exit 1; fi
tktcache=`klist -5 | sed -n 's/Ticket cache: FILE:\(.*\)/\1/p'`
cp /tmp/krb5cc_`id -u`*           $workDir

counter="0"
for x in `seq 0 1000`; do
  if [[ $counter -lt $nevents ]] ; then
    counter=$[$counter+1000]
    rFile="0"
    if [ -d "$outputDir/${dataset}" ] ; then
      rFile=`ls $outputDir/${dataset} | grep -c ${dataset}_${x}_ntuple.root  `
    fi
    if [ $rFile -gt 0 ]; then echo "      ---- file exists already ($dataset $x)"; continue; fi
  #echo "$counter"
    echo " $script $runMacro $outputDir $clusterstr $dataset $infile $nevents $x $workDir $channel "
    bsub -q 8nh  $script $runMacro $outputDir $clusterstr $dataset $infile $nevents $x $workDir $channel 
  fi
done

exit 0

