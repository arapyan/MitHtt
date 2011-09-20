{  
  if(gSystem->Getenv("CMSSW_VERSION")) {
    TString rfitpath("/afs/cern.ch/cms/sw/$SCRAM_ARCH/lcg/roofit/5.26.00-cms5/include");
    // TString rfitpath("/server/02a/cmsprod/cmssoft/slc5_amd64_gcc434/lcg/roofit/5.28.00a-cms3/include/");
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
    gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMitHttNtupleDefs.so");  
    // gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMitPhysicsFakeMods.so");  
  }
  // else {
  //   cout << "initializing with hardcoded CMSSW_BASE: " << endl;
  //   gROOT->Macro("/home/dkralph/cms/cmssw/020/CMSSW_4_1_3/src/Common/CPlot.cc+");
  //   gROOT->Macro("/home/dkralph/cms/cmssw/020/CMSSW_4_1_3/src/Common/MitStyleRemix.cc+");
  //   gROOT->Macro("/home/dkralph/cms/cmssw/020/src/MitAna/macros/setRootEnv.C+");
  //   gSystem->Load("/home/dkralph/cms/cmssw/020/CMSSW_4_1_3/lib/slc5_amd64_gcc434/libMitCommonMathTools.so");
  //   gSystem->Load("/home/dkralph/cms/cmssw/020/CMSSW_4_1_3/lib/slc5_amd64_gcc434/libMitHttNtupler.so");  
  // }

  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
