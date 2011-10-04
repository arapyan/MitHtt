#!/bin/bash
# note: cannot give this two versions of the same dataset
config=htt.config
if [ "` hostname | grep '\.mit\.edu'`" ]; then
    catalog=/home/cmsprod/catalog
    outputDir=/scratch/$USER/htt/nojetid
else
    catalog=/home/mitprod/catalog
    outputDir=/data/blue/$USER/htt/nojetid
fi
runMacro=runHttNtupler.C
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
  json=`echo $line      	| tr -s ' ' | cut -d ' ' -f 8`

  ./submitjobs.sh $outputDir $runMacro $dataset $book $catalog $isdata $usegen $nevents $skiphltfail $json

  sleep 0.5
done

exit 0
