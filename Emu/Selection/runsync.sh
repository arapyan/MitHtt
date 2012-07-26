#!/bin/bash

#enum { kNo, kDown, kUp }; jet and btagging uncertainties

lumi=5000
ntupledir=.
plotdir=.
root -b -q -l selectEmu.C+\(\"sync.conf\",\"$ntupledir\",$lumi,0,0,0,0,0\)  #2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check' # central value
