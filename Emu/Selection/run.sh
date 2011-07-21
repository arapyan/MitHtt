#enum { kCenter, kScDown, kScUp, kResDown, kResUp }; escale uncertainties
#enum { kNo, kDown, kUp }; jet uncertainties
ntupledir=/scratch/$USER/htt/selections/test
plotdir=/scratch/$USER/htt/selections/test
root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",1084,0\)  # central value
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",1084\) 2>&1 | grep -v 'Info in <TCanvas::'
#
#------------------------------------------------------------------------------------------------------------
#
#ntupledir=/scratch/dkralph/htt/selections/jet-down
#plotdir=/scratch/dkralph/htt/selections/jet-down
##root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",1149,1\)  # scale jets down
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",1149\) 2>&1 | grep -v 'Info in <TCanvas::'
#ntupledir=/scratch/dkralph/htt/selections/jet-up
#plotdir=/scratch/dkralph/htt/selections/jet-up
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",1149,2\)  # scale jets up
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",1149\) 2>&1 | grep -v 'Info in <TCanvas::'

