#!/bin/bash
#===================================================================================================
# Script to process our complete analysis task.
#
# note: cannot give this two versions of the same dataset
#                                                                             Ch.Paus (Aug 15, 2010)
#===================================================================================================

echo " Config: ${MIT_HTT_DIR}/config/${MIT_PROD_CFG}.txt"
for dataset in `cat ${MIT_HTT_DIR}/config/${MIT_PROD_CFG}.txt|grep -v ^#|tr -s ' '|cut -d' ' -f 2`
do
  # find the line to this dataset and do further analysis
  line=`grep $dataset ${MIT_HTT_DIR}/config/${MIT_PROD_CFG}.txt`
  # determine the MIT version of this sample
  export VERSION=`echo $line | tr -s ' ' | cut -d ' ' -f 1`
  # find potential JSON file
  export MIT_PROD_JSON=`echo $line | tr -s ' ' | cut -d ' ' -f 7`
  # find potential cut to remove overlap
  export MIT_PROD_OVERLAP=`echo $line | tr -s ' ' | cut -d ' ' -f 6`
  echo " JSON: $MIT_PROD_JSON  Overlap: $MIT_PROD_OVERLAP"
  # now submit the sucker
  submit.sh $MIT_PROD_MACRO $MIT_CATALOG $MIT_PROD_BOOK/$VERSION $dataset noskim \
            $MIT_PROD_CFG    $MIT_PROD_HIST 1
done

exit 0
