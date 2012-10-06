#!/bin/bash

# Macro parameters
MACRO="applySVFit.C"

# Default destination directory
if [ -z "$DEST" ]; then
    echo "Destination path must be specified"
    exit 0
fi

# Check CMSSW_BASE
if [ -z ${CMSSW_BASE+x} ]; then
    echo "CMSSW_BASE is not set. Please set the CMSSW environment."
    echo "Aborting..."
    exit 0
fi

# Set CMSSW environment
# (necessary because of a flaw in how lxbatch handles environment variables)
echo "Setting CMSSW environment"
cd $CMSSW_BASE
eval `scram runtime -sh`
echo -n "Working directory: "
cd -

# Echo parameters
echo "MACRO"      $MACRO
echo "INPUT"      $INPUT
echo "OUTPUT"     $OUTPUT
echo "DEST"       $DEST
echo "SKIPEVENTS" $SKIPEVENTS
echo "EVENTS"     $EVENTS
echo "LEPTYPE1"   $LEPTYPE1
echo "LEPTYPE2"   $LEPTYPE2
echo "ZEROONLY"   $ZEROONLY

# Copy macro to current working directory
DIR=`dirname $0`
cp $DIR/rootlogon.C .
cp $DIR/$MACRO .

# Run macro
root -l -b -q ./rootlogon.C ${MACRO}+"(\"${INPUT}\",\"${OUTPUT}\",${SKIPEVENTS},${EVENTS},\"${LEPTYPE1}\",\"${LEPTYPE2}\", $ZEROONLY)"

# Copy output to destination
rsync -avz $OUTPUT $DEST/
