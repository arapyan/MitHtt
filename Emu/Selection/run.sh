#enum { kNo, kDown, kUp }; jet and btagging uncertainties

lumi=4587
ntupledir=/home/vdutta/htt/selections/test
plotdir=/home/vdutta/htt/selections/test
root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,0\)  2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check' # central value
root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,0\) 2>&1 | grep -v 'Info in <TCanvas::'

#------------------------------------------------------------------------------------------------------------

# apply up/down electron scale/resolution
#ntupledir=/home/vdutta/htt/selections/test/e-down
#plotdir=/home/vdutta/htt/selections/test
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,1\)  2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check' # electron energy scale down
#root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,1\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/home/vdutta/htt/selections/test/e-up
#plotdir=/home/vdutta/htt/selections/test
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,0,2\)  2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check' # electron energy scale up
#root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi,2\) 2>&1 | grep -v 'Info in <TCanvas::'

#------------------------------------------------------------------------------------------------------------

# scale btag efficiency up and down
#ntupledir=/home/vdutta/htt/selections/test/btag-down
#plotdir=/home/vdutta/htt/selections/test/btag-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,1,0,0,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale btag efficiency down
#root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/home/vdutta/htt/selections/test/btag-up
#plotdir=/home/vdutta/htt/selections/test/btag-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,2,0,0,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale btag efficiency up
#root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'

#------------------------------------------------------------------------------------------------------------

# scale jet energy up and down
#ntupledir=/home/vdutta/htt/selections/test/jet-down
#plotdir=/home/vdutta/htt/selections/test/jet-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,1,0,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale jets down
#root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/home/vdutta/htt/selections/test/jet-up
#plotdir=/home/vdutta/htt/selections/test/jet-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,2,0,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale jets up
#root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'

#------------------------------------------------------------------------------------------------------------

# scale b mistag rate up and down
#ntupledir=/home/vdutta/htt/selections/test/mistag-down
#plotdir=/home/vdutta/htt/selections/test/mistag-down
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,1,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale mistag rate down
#root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/home/vdutta/htt/selections/test/mistag-up
#plotdir=/home/vdutta/htt/selections/test/mistag-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",$lumi,0,0,2,0\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'  # scale mistage rate up
#root -b -q -l rootlogon2.C plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",$lumi\) 2>&1 | grep -v 'Info in <TCanvas::'
