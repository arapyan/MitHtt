//--------------------------------------------------------------------------------------------------
// 
// ROOT macro to print generator level information from BAMBU
//
//==================================================================================================

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitHtt/Ntupler/interface/BambuGenDumperMod.hh"
#endif

using namespace mithep;

//--------------------------------------------------------------------------------------------------
void runBambuGenDumper(const char *file = "/castor/cern.ch/user/p/paus/filefi/020/p11-wm-v1g1-pu/B4AF460D-FD4F-E011-B619-0015178C4B24.root")
{  
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  BambuGenDumperMod *mod = new BambuGenDumperMod;

  // set up analysis
  Analysis *ana = new Analysis;
  ana->SetSuperModule(mod);
  ana->AddFile(file);
  ana->SetUseHLT(kFALSE);
  ana->SetProcessNEvents(10);
  
  // run the analysis after successful initialisation
  ana->Run(kTRUE);
}
