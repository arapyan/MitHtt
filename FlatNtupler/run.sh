#!/bin/bash

#enum { kNo, kDown, kUp }; jet and btagging uncertainties

lumi=1.0
#ntupledir=/data/blue/arapyan/httprod/hpstau/tautau/
#root -b -q -l selectEmu.C+\(\"htt.conf\",\"$ntupledir\",$lumi,1,0,0,0,0\) 

#ntupledir=/data/blue/arapyan/httprod/mutau/
#root -b -q -l selectMuTau.C+\(\"mutau.conf\",\"$ntupledir\",$lumi,1\)  

#ntupledir=/data/blue/arapyan/httprod/etau/
#root -b -q -l selectETau.C+\(\"etau.conf\",\"$ntupledir\",$lumi,1\)  

ntupledir=/data/blue/arapyan/httprod/ttttnew/
root -b -q -l selectTauTau.C+\(\"tautau.conf\",\"$ntupledir\",$lumi,1\) 


