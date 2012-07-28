#!/bin/bash
# note: cannot give this two versions of the same dataset
config=htt.config
if [ "` hostname | grep '\.mit\.edu'`" ]; then
    catalog=/home/cmsprod/catalog
    outputDir=/scratch/$USER/ntuple-prod
else
    catalog=~/mitprod/catalog
    outputDir=/store/user/arapyan/production/2012
fi
runMacro=runHttNtupler.C
cmsMkdir $outputDir

# file to keep track of number of requests to castor stager
echo 0 > /tmp/requests

for dataset in `cat $config | grep -v ^# | tr -s ' ' | sed "s@\s\s@ @g" | awk '{print $1}' `; do
  line=`grep -v ^# $config | grep $dataset `
  book=`echo    $line           | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $2}'`/filefi
  version=`echo $line           | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $3}'`
  book=$book/$version
  isdata=`echo $line            | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $4}'`
  usegen=`echo $line            | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $5}'`
  nevents=`echo $line           | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $6}'`
  skiphltfail=`echo $line       | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $7}'`
  is2012=`echo $line       | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $8}'`
  json=`echo $line              | tr -s ' '  | sed "s@\s\s@ @g" | awk '{print $9}'`

  ./submitjobs.sh $outputDir $runMacro $dataset $book $catalog $isdata $usegen $nevents $skiphltfail $is2012 $json

  sleep 0.5
done

rm /tmp/requests

exit 0
