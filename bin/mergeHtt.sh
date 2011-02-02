#!/bin/bash
#===================================================================================================
# Script to merge all historgram files from a complete analysis task such that each sample with
# a well defined cross section is in one 'merged' file.
#
#                                                                             Ch.Paus (Aug 15, 2010)
#===================================================================================================
# FilenameHead has to be MIT_ANA_CFG sometimes, sometimes MIT_PROD_CFG


source envir.sh

mergeHist.py \
     --InputPath=$MIT_PROD_HIST/$MIT_PROD_CFG/$MIT_PROD_BOOK \
    --OutputPath=$MIT_PROD_HIST/$MIT_PROD_CFG/$MIT_PROD_BOOK/merged \
    --FilenameHeader=$MIT_PROD_CFG \
    --DatasetListFile=$MIT_HTT_DIR/config/${MIT_PROD_CFG}.txt
