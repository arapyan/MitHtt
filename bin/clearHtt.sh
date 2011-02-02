#!/bin/bash
#===================================================================================================
# Script to delete all histogram files of the given entire analysis task.
#
#                                                                             Ch.Paus (Aug 15, 2010)
#===================================================================================================

cd $MIT_PROD_HIST/$MIT_PROD_CFG/$MIT_PROD_BOOK/
rm */*.root

exit 0
