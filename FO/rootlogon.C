{  
  if (gSystem->Getenv("CMSSW_VERSION")) {
    TString rfitpath("/server/02a/cmsprod/cmssoft/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms3/include/");
    // TString rfitpath("/afs/cern.ch/cms/sw/$SCRAM_ARCH/lcg/roofit/5.26.00-cms5/include");
    TString path = gSystem->GetIncludePath();
    path += "-I. -I$ROOTSYS/src -I";
    path += rfitpath;
    gSystem->SetIncludePath(path.Data());
    
    TString str = gSystem->GetMakeSharedLib();
    if (str.Contains("-m32")==0 && str.Contains("-m64")==0) {
      str.ReplaceAll("g++", "g++ -m32");
      gSystem->SetMakeSharedLib(str);
    }      
    
    gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
    gROOT->Macro("$CMSSW_BASE/src/Common/CPlot.cc+");
    gROOT->Macro("$CMSSW_BASE/src/Common/MitStyleRemix.cc+");
    gROOT->Macro("$CMSSW_BASE/src/Common/CEffUser2D.cc+");
    
    gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMitHttNtupler.so");  
  }

  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
