{
  if (gSystem->Getenv("CMSSW_VERSION")) {
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32") == 0 && str.Contains("-m64") == 0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }
  }
  // keeping from loading this so many times
  TString addedLibs(gSystem->GetLibraries());
  if(!addedLibs.Contains("setRootEnv_C.so")) {
    gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
  }

  gROOT->SetMacroPath(TString(gROOT->GetMacroPath()) + TString(":")
		      +TString(gSystem->Getenv("CMSSW_BASE"))+"/src/MitHtt/macros/"+TString(":"));

  loadmylib("MitPhysics","Mods");
  loadmylib("MitPhysics","SelMods");
  loadmylib("MitPlots",  "Style");
  loadmylib("MitPlots",  "Input");
  loadmylib("MitPlots",  "Plot");

  loadmylib("MitHtt",    "Mods");

  loadmylib("MitHtt",    "Mods");
  loadmylib("MitCommon","DataFormats");
//   loadmylib("MitPhysics","Validation");
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
