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
   scramdir=${12}
     is2012=${13}
       json=${14} 

ulimit -v 3000000
workDir=`pwd`
echo `hostname`
echo "workDir: $workDir"
echo "workDir: $scramdir"
echo "args:    $*"

cd ${scramdir}/src/
SCRAM_ARCH=slc5_amd64_gcc462
eval `scramv1 runtime -sh`
cd $workDir

soFile=`echo $runMacro | sed 's/\./_/'`.so
cp ${scramdir}/${soFile}   .
cp ${scramdir}/${runMacro} .
cp ${scramdir}/rootlogon.C .
cp ${scramdir}/$json       .

#mkdir -p /tmp/$USER/.krb5
#chmod og-rwx /tmp/$USER/.krb5/
#cp ${scramdir}/krb5cc_`id -u`_$clusterid /tmp/$USER/.krb5/ticket.$clusterid
#export KRB5CCNAME=/tmp/$USER/.krb5/ticket.$clusterid
klist -5
#cp ${scramdir}/x509up_u`id -u` /tmp/$USER/.krb5/
#export X509_USER_PROXY=/tmp/$USER/.krb5/x509up_u`id -u`

root -b -l -q ./rootlogon.C $runMacro+\(\"${fileset}\",\"${dataset}\",\"${book}\",\"${catalog}\",${isdata},${usegen},${nevents},${skiphltfail},${is2012},\"${json}\"\)

oldsize=`cmsLs ${outputDir}/${dataset}/${dataset}_${fileset}_ntuple.root | awk '{print $2}'`
newsize=`ls -ltr  ${dataset}_${fileset}_ntuple.root | awk '{print $5}'`
echo "old size : "$oldsize
echo "new size : "$newsize


status=`echo $?`
echo "Status - $status"

if [[ $oldsize -lt $newsize ]]; then
    cmsMkdir                                     ${outputDir}/${dataset}
    cmsMkdir                                     ${outputDir}/${dataset}_json
    cmsStage ${dataset}_${fileset}_ntuple.root   ${outputDir}/${dataset}
    cmsStage ${dataset}_${fileset}_ntuple.json   ${outputDir}/${dataset}_json
fi

exit $status
