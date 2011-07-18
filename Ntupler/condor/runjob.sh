#!/bin/bash

   runMacro=$1
  outputDir=$2
  clusterid=$3

    fileset=$4
    dataset=$5
       book=$6
    catalog=$7
     isdata=$8
     usegen=$9
    nevents=${10}
skiphltfail=${11}
       json=${12}

scramdir=/home/$USER/cms/cmssw/$MIT_VERS/$CMSSW_VERSION/src
workDir=`pwd`
echo `hostname`
echo "workDir: $workDir"
echo "args:    $*"
ls -l

tar pxf job_files.tar
mkdir -p /tmp/$USER/.krb5
chmod og-rwx /tmp/$USER/.krb5/
mv tmp/krb5cc_`id -u`* /tmp/$USER/.krb5/ticket.$clusterid
export KRB5CCNAME=/tmp/$USER/.krb5/ticket.$clusterid
klist -5

cd $scramdir
SCRAM_ARCH=slc5_amd64_gcc434
eval `scramv1 runtime -sh`
cd $workDir

root -b -l -q ./rootlogon.C \
  $runMacro+\(\"${fileset}\",\"${dataset}\",\"${book}\",\"${catalog}\",${isdata},${usegen},${nevents},${skiphltfail},\"${json}\"\)

status=`echo $?`
echo "Status - $status"

mkdir -p $outputDir
mv ${dataset}_${fileset}_ntuple.root $outputDir/
mv ${dataset}_${fileset}_ntuple.json $outputDir/

rm -vf /tmp/$USER/.krb5/ticket.$clusterid

exit $status
