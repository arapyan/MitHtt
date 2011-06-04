{
  // keeping from loading this so many times
  TString addedLibs(gSystem->GetLibraries());
  if(!addedLibs.Contains("setRootEnv_C.so")) {
    gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
  }
  
  loadLibraries("libMitHttNtupler.so");
}
