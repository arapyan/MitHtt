root -b -l -q predictFakesMuEle.C+\(\"full2011\"\) 2>&1 | grep -v 'Info in <Minuit2' | grep -v 'Mass Check'

rm *.so *.d
