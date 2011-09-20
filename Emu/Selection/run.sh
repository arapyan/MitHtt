#enum { kCenter, kScDown, kScUp, kResDown, kResUp }; escale uncertainties
#enum { kNo, kDown, kUp }; jet uncertainties
ntupledir=/data/blue/dkralph/htt
#plotdir=/scratch/$USER/htt/selections/b-scale-2
root -b -q -l Selector.cc+ selectEmu.C+\(\"emu.conf\",\"$ntupledir\",1597\)  # central value
#root -b -q -l plotEmu.C+\(\"emu.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",1084\) 2>&1 | grep -v 'Info in <TCanvas::'
