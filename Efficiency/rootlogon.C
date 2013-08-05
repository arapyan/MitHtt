{  
  TString path = gSystem->GetIncludePath();
  path += " -I$ROOFITSYS/include";
  gSystem->SetIncludePath(path.Data());      
  gROOT->Macro("$CMSSW_BASE/src/MitAna/macros/setRootEnv.C+");
  gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/CPlot.cc+");

  gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/RooVoigtianShape.cc+");
  gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/RooErf.cc+");
  gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/RooCMSShape.cc+");
  
  gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/CEffUser1D.cc+");
  gROOT->Macro("$CMSSW_BASE/src/MitHtt/Common/CEffUser2D.cc+");
               
  gSystem->Load("$CMSSW_BASE/lib/$SCRAM_ARCH/libMitHttNtupler.so");

  // Show which process needs debugging
  gInterpreter->ProcessLine(".! ps |grep root.exe");
}
