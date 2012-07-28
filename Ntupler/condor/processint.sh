#!/bin/bash

config=/afs/cern.ch/work/a/arapyan/httflat/httint.conf
outputDir=/afs/cern.ch/work/a/arapyan/httflat
runMacro=apply.C
#mkdir $outputDir

for dataset in `cat $config | grep -v ^# | tr -s ' ' | sed "s@\s\s@ @g" | awk '{print $1}' `; do
  line=`grep -v ^# $config | grep ^$dataset `
  infile=`echo $line           | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $1}'`_select.root
  nevents=`echo $line           | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $2}'`
  channel=`echo $line           | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $3}'`

   	
  ./submitint.sh $outputDir $runMacro $dataset $infile $nevents $channel

  sleep 0.5
done

exit 0
