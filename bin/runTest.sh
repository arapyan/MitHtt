#!/bin/bash

export MIT_PROD_JSON="~"
#export MIT_PROD_JSON="Cert_136033-149442_7TeV_Dec22ReReco_Collisions10_JSON.txt"
export MIT_PROD_OVERLAP="-1.0"

mkdir -p /home/$USER/cms/hist/test

run.sh $src/MitHtt/macros/runEMU.C $MIT_CATALOG $MIT_PROD_BOOK/$MIT_VERS \
          p11-zttm20-v1g1-pu noskim 0000 $MIT_PROD_CFG /home/$USER/cms/hist/test 1 1335 100

#p11-zttm20-v1g1-pu
#r10a-mu-d22
#w10-ztt-powheg-c10-v8-pu11 
#w10-h160wwlt-gf-z2-v8-pu11
#w10-gghtt120-v8-pu11
#w10-wjetsl-z2-v8-pu11 
