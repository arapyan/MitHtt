ntupledir=/scratch/dkralph/htt/selections/ee-test
plotdir=/scratch/dkralph/htt/selections/ee-test
#root -b -q -l selectee.C+\(\"ee.conf\",\"$ntupledir\",865\)
#root -b -q -l plotee.C+\(\"ee.conf\",\"$ntupledir/ntuples\",\"$plotdir\",\"png\",190\) 2>&1 | grep -v 'Info in <TCanvas::'

root -l RooVoigtianShape.cc+ escalefit.C+\(\"fit.conf\"\)
