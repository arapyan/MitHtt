#root -b -l -q selectEleFO.C+\(\"data_smu.fo\",\"865\",35,1\)
#root -b -l -q computeFakeRate.C+\(\"elefr.input\",\"png\",1\)
root -b -l -q predictFakesMuEle.C+\(\"865/FakeRate/fr.root\",\"07-july/pt10\"\)

rm *.so *.d
