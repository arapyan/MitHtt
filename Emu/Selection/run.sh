#ntupledir=full
#plotdir=full
#root -b -q selectEmu.C+\(\"emu.conf\",\"/scratch/dkralph/htt/selections/full\",190\)
root -b -q selectEmu.C+\(\"emu.conf\",\"test\",190\)

# to get right binning for limits you have to write to a directory name containing "limit"
#root -b -q plotEmu.C+\(\"emu.conf\",\"/scratch/dkralph/htt/selections/full/ntuples\",\"/scratch/dkralph/htt/selections/full\",\"png\",190\)
#root -b -q plotEmu.C+\(\"emu.conf\",\"test/ntuples\",\"test\",\"png\",190\)
