#!/bin/bash
#
# For each dataset that appears in the histogram directory, loops 
# through the catalog checking that a .root file exists in 
# the histogram directory for each fileset in the catalog
#
#dset=w10-zee-powheg-c10-v8-pu11

config=hypha.config
if [ "` hostname | grep '\.mit\.edu'`" ]; then
    catalog=/home/cmsprod/catalog
    outputdir=/scratch/$USER/htt
else
    catalog=/home/mitprod/catalog
    outputdir=/data/blue/$USER/htt
fi

if [ ! "$CMSSW_BASE" ] || [ ! "$MIT_VERS" ]; then
  echo "-------> error: define cms environment."
  exit 1
fi

for dataset in `cat $config | grep -v ^# | tr -s ' ' | cut -d' ' -f 1`
do
  line=`grep -v ^# $config | grep $dataset`
  book=`echo 	$line 		| tr -s ' ' | cut -d ' ' -f 2`/filefi
  version=`echo $line 		| tr -s ' ' | cut -d ' ' -f 3`
  echo "===================================="
  echo $dataset; echo " "
  filesets=$catalog/$book/$version/$dataset/Filesets
  filenumbers=()
  for file in `ls $outputdir |grep $dataset`; do
#    tmp=`echo $file | cut -d'_' -f4 | cut -d'.' -f1`
    tmp=`echo $file | sed 's/.*\([0-9][0-9][0-9][0-9]\)\_ntuple.root/\1/'`
    filenumbers=("${filenumbers[@]}" "$tmp")
  done
  for fileset in `cat $filesets | cut -d' ' -f1`
  do
    found=0
    for filenum in `echo ${filenumbers[@]}`
    do
#      echo "looking for $filenum and $fileset"; echo " "
      if [ "$filenum" == "$fileset" ]; then
        echo "found $filenum"
        found=1
        break
      fi
    done
    if [ "$found" == "0" ]; then
      echo "---------------> missing $fileset"
    fi
  done
done

exit 0