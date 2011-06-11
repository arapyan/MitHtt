ntupledir=/scratch/dkralph/htt/selections/test
plotdir=/scratch/dkralph/htt/selection/test
#root -b -q -l selectEmu.C+\(\"emu.conf\",\"$ntupledir\",190\)
root -b -q -l selectMuMu.C+\(\"mumu.conf\",\"$ntupledir\",190\)
#root -l -b -q plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",190\) 2>&1 | grep -v 'Info in <TCanvas::'
