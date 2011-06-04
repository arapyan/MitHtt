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

workDir=`pwd`            # should be run from $src/xxx/condor
script=runjob.sh
clusterstr='$(Cluster)'  # otherwise it doesn't get expanded properly
soFile=`echo $runMacro | sed 's/\./_/'`.so

if [ ! "$CMSSW_BASE" ] || [ ! "$MIT_VERS" ]; then
  echo "-------> error: define cms environment."
  exit 1
fi
if [ ../macros/$runMacro -nt ../macros/$soFile ]; then
  echo "-------> error: forgot to recompile run macro."
  exit 1
fi

cp ../macros/rootlogon.C             $workDir
cp ../macros/$runMacro               $workDir
cp ../macros/$soFile                 $workDir

klist -5 &> /dev/null
[ "$?" != "0" ] && ( echo "-------> warning: check credentials cache." )
tktcache=`klist -5 | sed -n 's/Ticket cache: FILE:\(.*\)/\1/p'`
tar pcf job_files.tar rootlogon.C $runMacro $soFile $tktcache &> /dev/null
if [ "$?" != "0" ]; then
  echo "-------> error: tar failed."
  exit 1
fi
chmod og-rwx job_files.tar
mkdir -p /tmp/$USER            # just so it's not lying around on /home
mv job_files.tar /tmp/$USER

# Create the directory for the condor results
datestring=`date +-%g%m%d-%H`
condorOutDir=${workDir}/condor${datestring}/$dataset
mkdir -p $condorOutDir

# Looping through each single fileset and submitting the condor jobs
filesets=$catalog/$book/$dataset/Filesets
for fileset in `cat $filesets | cut -d' ' -f1 `; do
#     if [ "$fileset" != "0000" ]; then
#       continue
#     fi
  # check if the output already exists
  rFile="$outputDir/${dataset}_${fileset}_ntuple.root"
  if [ -f "$rFile" ]; then
    echo "      ---- file exists already ($rFile)"
    continue
  fi
  logFile=`echo $book/$dataset/$fileset | tr '/' '+'`
  logFile=/tmp/$USER/$logFile
  mkdir -p /tmp/$USER
  rm    -f $logFile
  echo "      $script $runMacro $outputDir $clusterstr $fileset $dataset $book $catalog $isdata $usegen $nevents $skiphltfail"
  
cat > submit.cmd <<EOF
Universe                = vanilla
Requirements            = (Arch == "X86_64") && (OpSys == "LINUX") && (Disk >= DiskUsage) && ((Memory * 1024) >= ImageSize) && (HasFileTransfer)
Notification            = Error
Executable              = $script
Arguments               = $runMacro $outputDir $clusterstr $fileset $dataset $book $catalog $isdata $usegen $nevents $skiphltfail
Rank                    = Mips
GetEnv                  = True
Initialdir              = $workDir
Input                   = /dev/null
Output                  = $condorOutDir/${fileset}.out
Error                   = $condorOutDir/${fileset}.err
Log                     = $logFile
job_lease_duration      = 36000
should_transfer_files   = YES
transfer_input_files    = /tmp/$USER/job_files.tar
when_to_transfer_output = ON_EXIT
Queue
EOF

  condor_submit submit.cmd >& /dev/null;
  rm submit.cmd
done

exit 0

#+AccountingGroup        = "research.dkralph"
