{
  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
  gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/CPlot.cc+");
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMitHttNtupler.so");
}
