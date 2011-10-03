#enum { kNo, kDown, kUp }; jet and btagging uncertainties
lumi=1597
ntupledir=/home/vdutta/htt/selections/susy11
plotdir=/home/vdutta/htt/selections/susy11
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0\)  # central value
root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#------------------------------------------------------------------------------------------------------------

# scale btag efficiency up and down
#ntupledir=$HOME/htt/selections/btag-down
#plotdir=$HOME/htt/selections/btag-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,1,0,0\)  # scale jets down
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=$HOME/htt/selections/btag-up
#plotdir=$HOME/htt/selections/btag-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,2,0,0\)  # scale jets up
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'

# scale jet energy up and down
#ntupledir=$HOME/htt/selections/jet-down
#plotdir=$HOME/htt/selections/jet-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,1,0\)  # scale down
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=$HOME/htt/selections/jet-up
#plotdir=$HOME/htt/selections/jet-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,2,0\)  # scale up
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'

# scale b mistag rate up and down
#ntupledir=$HOME/htt/selections/mistag-down
#plotdir=$HOME/htt/selections/mistag-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,1\)  # scale down
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=$HOME/htt/selections/mistag-up
#plotdir=$HOME/htt/selections/mistag-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,2\)  # scale up
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
