{    
  if(gSystem->Getenv("CMSSW_VERSION")) {
    TString rfitpath("/afs/cern.ch/cms/sw/slc5_amd64_gcc434/lcg/roofit/5.26.00-cms6/include");
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
    
    gSystem->Load("$CMSSW_BASE/lib/slc5_amd64_gcc434/libMitHttNtupler.so");
    
    gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/CPlot.cc+");
    gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/MitStyleRemix.cc+"); 
    
  } else {

    TString mypath("/home/$USER/Include/");
    TString path = gSystem->GetIncludePath();
    path += " -I";
    path += mypath;
    gSystem->SetIncludePath(path.Data());      
    gROOT->Macro(mypath+TString("/Goodies/RooVoigtianShape.cc+"));
    gROOT->Macro(mypath+TString("/Goodies/RooErf.cc+"));
    gROOT->Macro(mypath+TString("/Goodies/RooCMSShape.cc+"));
    
    gROOT->Macro(mypath+TString("/Common/Efficiency/CPlot.cc+"));
    gROOT->Macro(mypath+TString("/Common/Efficiency/MitStyleRemix.cc+")); 
    gROOT->Macro(mypath+TString("/Common/Efficiency/CEffUser1D.cc+"));
    gROOT->Macro(mypath+TString("/Common/Efficiency/CEffUser2D.cc+")); 
  }
               
  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
