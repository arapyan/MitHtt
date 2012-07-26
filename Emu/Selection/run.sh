#!/bin/bash

#enum { kNo, kDown, kUp }; jet and btagging uncertainties

lumi=5065
ntupledir=.
plotdir=.
root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,0\) 
root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,0,0\)

#------------------------------------------------------------------------------------------------------------

# apply up/down electron scale/resolution
#ntupledir=e-down
#plotdir=e-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,1\) 
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,1,0\)
#ntupledir=e-up
#plotdir=e-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,2\)
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,2,0\)

