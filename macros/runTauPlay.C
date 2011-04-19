// $Id:$
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"

#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PFTauIDMod.h"

#include "MitAna/DataTree/interface/CaloJetCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/MetCol.h" 
#include "MitAna/DataTree/interface/CaloMetCol.h"
#include "MitAna/DataTree/interface/PFMetCol.h"

#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/CaloTauCol.h"
#include "MitAna/DataTree/interface/TauCol.h"

#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Mods/interface/EffMod.h"

#include "MitHtt/Mods/interface/TauPlayAnalysis.h"

#endif
//--------------------------------------------------------------------------------------------------
void runTauPlay(const char *fileset    = "0000",
		const char *skim       = "noskim",
		const char *dataset    = "p11-zttm20-v1g1-pu",
		const char *book       = "local/filefi/020",
		const char *catalogDir = "/home/cmsprod/catalog",
		int         nEvents    = 1000)
{
  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace std;
  using namespace mithep;

  //------------------------------------------------------------------------------------------------
  // object id sequence
  //------------------------------------------------------------------------------------------------

  PFTauIDMod          *pfTauId        = new PFTauIDMod;
  pfTauId->SetIsHPSSel(kTRUE);
  pfTauId->SetPtMin(20.0);

  //------------------------------------------------------------------------------------------------
  // analysis modules
  //------------------------------------------------------------------------------------------------
  TauPlayAnalysis *analysisTauPlay = new TauPlayAnalysis;
  analysisTauPlay->SetTausName      (ModNames::gkGoodPFTausName);

  //------------------------------------------------------------------------------------------------
  // making analysis chain -- the list of modules that will run on each event
  //------------------------------------------------------------------------------------------------
  pfTauId           ->Add(analysisTauPlay);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------

  Analysis *ana = new Analysis;
  ana->SetKeepHierarchy(kTRUE);
  ana->SetProcessNEvents(nEvents);
  ana->SetSuperModule(pfTauId);
  ana->SetPrintScale(100);

  //------------------------------------------------------------------------------------------------
  // organize input
  //------------------------------------------------------------------------------------------------
  Catalog *c = new Catalog(catalogDir);
  TString skimdataset = TString(dataset)+TString("/") +TString(skim);
  Dataset *d = NULL;
  if (TString(skim).CompareTo("noskim") == 0)
    d = c->FindDataset(book,dataset,fileset);
  else 
    d = c->FindDataset(book,skimdataset.Data(),fileset);
  ana->AddDataset(d);
  //ana->AddFile("root://castorcms//castor/cern.ch/user/p/paus/filler/011/s09-ttbar-7-mc3/*.root");
  //ana->AddFile("hgg-skim_r10a-eg-pr-v4_noskim_0000_000.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  ana->SetOutputName("tauplay.root");
  ana->SetCacheSize(64*1024*1024);

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}
