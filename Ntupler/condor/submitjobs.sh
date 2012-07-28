#!/bin/bash
echo " "
echo "  --> $*"

  outputDir=$1
   runMacro=$2
    dataset=$3
       book=$4
    catalog=$5
     isdata=$6
     usegen=$7
    nevents=$8
skiphltfail=$9
     is2012=${10}
       json=${11}


workDir=$CMSSW_BASE #`pwd`            # should be run from $src/xxx/condor
iscastor=`echo $outputDir | grep -c castor`
script=runjob.sh
if [ $iscastor -gt 0 ]; then
    script=runjobsCastor.sh
else 
    script=runjobsCMST3.sh
fi

clusterstr=`ls -ltr /tmp/ | grep krb | grep \`id -u\` | sed "s@_@ @g" | tail -1 | awk '{print $11}'`
#'$(Cluster)'  # otherwise it doesn't get expanded properly
soFile=`echo $runMacro | sed 's/\./_/'`.so

# check a few things
if [ ! "$CMSSW_BASE" ]; then
  echo "-------> error: define cms environment."
  exit 1
fi
if [ macros/$runMacro -nt macros/$soFile ]; then
  echo "-------> error: forgot to recompile run macro."
  exit 1
fi

# grab the files needed for the job
cp ../macros/rootlogon.C             $workDir
cp ../macros/$runMacro               $workDir
cp ../macros/$soFile                 $workDir

# set up kerberos credentials
if ! klist -5 &> /dev/null; then echo "-------> warning: check credentials cache."; exit 1; fi
tktcache=`klist -5 | sed -n 's/Ticket cache: FILE:\(.*\)/\1/p'`
cp /tmp/krb5cc_`id -u`*           $workDir

#if ! [ -f $X509_USER_PROXY ]
#then 
#  echo " ERROR -- missing grid proxy ($X509_USER_PROXY)."
#fi
#cp /tmp/x509up_u`id -u`          $workDir

# retrieve number of previous castor requests
reqs=`cat /tmp/$LOGNAME/requests`

# Looping through each single fileset and submitting the condor jobs
filesets=$catalog/$book/$dataset/Filesets
for fileset in `cat $filesets | cut -d' ' -f1 `; do
  # check if the output already exists
 if [ $iscastor -gt 0 ]; then
     rFile=`nsls $outputDir/${dataset} | grep -c ${dataset}_${fileset}_ntuple.root  `
 else 
     rFile=`cmsLs $outputDir/${dataset} | grep -c ${dataset}_${fileset}_ntuple.root  `
 fi
 if [ $rFile -gt 0 ]; then echo "      ---- file exists already ($dataset $fileset)"; continue; fi

  # check if job already in condor queue
  if grep $dataset /tmp/$LOGNAME/condor_q-output.txt | grep " $fileset " &> /dev/null; then
      echo "      ---- job already running ($dataset $fileset)"
      continue
  fi

  # check if staged to castor disk cache
  if [ "` hostname | grep '\.cern\.ch'`" ]; then
      allstaged=yes
      dir=`head -n 1 $catalog/$book/$dataset/Filesets | awk '{print $2}'`
      # loop over each file in the fileset
      for file in `grep $fileset $catalog/$book/$dataset/Files | awk '{print $2}'`; do
	  # check current number of castor requests (doesn't count requests from the actual jobs)
	  if [[ $reqs -gt 2500 ]]; then echo -e "\n\n   Reached 2500 requests. Sleeping..."; sleep 501; echo 0 > /tmp/$LOGNAME/requests; reqs=0; fi
	  staged=`stager_qry -M $dir/$file` # staged already?
	  reqs=$[$reqs+1]
	  staged=`echo $staged | awk '{print $3}'`
	  # stage it if necessary
	  if [ "$staged" == "STAGED" ]; then    # already staged
	      echo "$dir/$file staged" > /dev/null
	  elif [ "$staged" == "STAGEIN" ]; then # already requested to stage this file
	      allstaged=no
	  else # request to stage this file
	      echo -n "            -- staging file: ($dir/$file)"
	      get=`stager_get -M $dir/$file`
	      reqs=$[$reqs+1]
	      if echo $get | grep "SUBREQUEST_FAILED" &> /dev/null; then echo " SUBREQUEST_FAILED"; else echo ""; fi
	      allstaged=no
	  fi
      done
      if [ "$allstaged" == "no" ]; then # had to stage some of the files, so we'll have to come back and resubmit this fileset later
	  echo "      ---- staging files for: ($dir $fileset)"
	  continue
      fi
  fi
  echo  "$script  $runMacro $outputDir $clusterstr $fileset $dataset $book $catalog $isdata $usegen $nevents $skiphltfail $workDir $is2012 $json"
  bsub -q 8nh  $script $runMacro $outputDir $clusterstr $fileset $dataset $book $catalog $isdata $usegen $nevents $skiphltfail $workDir $is2012 $json
  #bsub -q 1nh  $script $runMacro $outputDir $clusterstr $fileset $dataset $book $catalog $isdata $usegen $nevents $skiphltfail $workDir $is2012 $json
done
# cache the total number of requests
echo $reqs > /tmp/$LOGNAME/requests
echo $dataset" "$fileset >> /tmp/$LOGNAME/condor_q-output.txt

exit 0

