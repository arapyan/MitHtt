#!/bin/bash
# note: cannot give this two versions of the same dataset
config=hypha.config
if [ "` hostname | grep '\.mit\.edu'`" ]; then
    catalog=/home/cmsprod/catalog
else
    catalog=/home/mitprod/catalog
fi
runMacro=runhypha.C
outputDir=/scratch/$USER/htt
mkdir -p $outputDir

for dataset in `cat $config | grep -v ^# | tr -s ' ' | cut -d' ' -f 1`
do
  line=`grep -v ^# $config | grep $dataset`
  book=`echo 	$line 		| tr -s ' ' | cut -d ' ' -f 2`/filefi
  version=`echo $line 		| tr -s ' ' | cut -d ' ' -f 3`
  book=$book/$version
  isdata=`echo $line  		| tr -s ' ' | cut -d ' ' -f 4`
  usegen=`echo $line 		| tr -s ' ' | cut -d ' ' -f 5`
  nevents=`echo $line 		| tr -s ' ' | cut -d ' ' -f 6`
  skiphltfail=`echo $line 	| tr -s ' ' | cut -d ' ' -f 7`

  ./submitjobs.sh $outputDir $runMacro $dataset $book $catalog $isdata $usegen $nevents $skiphltfail
  sleep 3

done

exit 0
