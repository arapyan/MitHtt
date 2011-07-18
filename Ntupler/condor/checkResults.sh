#!/bin/bash
# loops through all the condor output files in the directory arguments,
# and prints all lines in the .err and .out files that contain the
# listed strings.
#
parentdir=$src/MitHtt/Ntupler/condor
if [ ! "$*" ]; then echo "error: provide directories as args."; fi
for dir in $*; do
    outputdir=$CMSSW_BASE/src/MitHtt/Ntupler/condor/$dir
    for dset in `ls $outputdir` #w10-ggww-z2-v8-pu11 
      do
      dsetdir=$outputdir/$dset
      for file in `ls $dsetdir`
	do
	errors=`sed -n '/_condor_stderr/ d
                        /Muon.ptErr/     d
                        /no such\|not be\|[eE][rR][rR]\|[Ff]atal\|Attaching\|[fF]ail\|condor_exec.exe: Status - [^0]/ p' <$dsetdir/$file`
	if [ "$errors" ]
	    then
	    echo $dsetdir/$file
	    echo " "
	    echo "$errors"
	    echo "================="
	fi
      done
    done
done