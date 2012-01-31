#enum { kNo, kDown, kUp }; jet and btagging uncertainties

lumi=4587
ntupledir=/build/vdutta/htt/selections/2011
plotdir=/build/vdutta/htt/selections/2011
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,0\)  2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check' # central value
root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,0\) 2>&1 | grep -v 'Info in <TCanvas::'

#------------------------------------------------------------------------------------------------------------

# apply up/down electron scale/resolution
#ntupledir=/build/vdutta/htt/selections/2011/e-down
#plotdir=/build/vdutta/htt/selections/2011
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,1\)  2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check' # electron energy scale down
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,1\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/build/vdutta/htt/selections/2011/e-up
#plotdir=/build/vdutta/htt/selections/2011
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,2\)  2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check' # electron energy scale up
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,2\) 2>&1 | grep -v 'Info in <TCanvas::'

#------------------------------------------------------------------------------------------------------------

# scale btag efficiency up and down
#ntupledir=/build/vdutta/htt/selections/2011/btag-down
#plotdir=/build/vdutta/htt/selections/2011/btag-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,1,0,0,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale btag efficiency down
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/build/vdutta/htt/selections/2011/btag-up
#plotdir=/build/vdutta/htt/selections/2011/btag-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,2,0,0,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale btag efficiency up
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'

#------------------------------------------------------------------------------------------------------------

# scale jet energy up and down
#ntupledir=/build/vdutta/htt/selections/2011/jet-down
#plotdir=/build/vdutta/htt/selections/2011/jet-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,1,0,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale jets down
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/build/vdutta/htt/selections/2011/jet-up
#plotdir=/build/vdutta/htt/selections/2011/jet-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,2,0,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale jets up
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'

#------------------------------------------------------------------------------------------------------------

# scale b mistag rate up and down
#ntupledir=/build/vdutta/htt/selections/2011/mistag-down
#plotdir=/build/vdutta/htt/selections/2011/mistag-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,1,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale mistag rate down
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/build/vdutta/htt/selections/2011/mistag-up
#plotdir=/build/vdutta/htt/selections/2011/mistag-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,2,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale mistage rate up
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
