// $Id:$
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"

#include "MitPhysics/Init/interface/ModNames.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/PDFProducerMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCorrectionMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/TauCleaningMod.h"
#include "MitPhysics/Mods/interface/PFTauIDMod.h"
#include "MitPhysics/Mods/interface/PFTauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetIDMod.h"
#include "MitPhysics/Mods/interface/ElectronCleaningMod.h"
#include "MitPhysics/Mods/interface/PhotonCleaningMod.h"
#include "MitPhysics/Mods/interface/TauCleaningMod.h"
#include "MitPhysics/Mods/interface/JetCleaningMod.h"
#include "MitPhysics/Mods/interface/MergeLeptonsMod.h"

#include "MitAna/DataTree/interface/JetCol.h"
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

#include "MitHtt/Mods/interface/EMUAnalysis2.h"

#endif
//--------------------------------------------------------------------------------------------------
void runEMU(const char *fileset    = "0000",
	  const char *skim       = "noskim",
	  const char *dataset    = "p10-pj15-v36",
	  const char *book       = "local/filefi/016",
	  const char *catalogDir = "/home/cmsprod/catalog",
	  const char *outputName = "htt",
	  int         nEvents    = 1000)
{
  printf("\n==== Enter macro  ====\n");

  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024], src[1024];
  float overlapCut = -1;


  if (gSystem->Getenv("MIT_PROD_JSON"))
    sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  else {
    printf(" JSON file was not properly defined. EXIT!\n");
    return;
  } 
  TString jsonFile = TString("/home/cmsprod/json/") + TString(json);
  Bool_t  isData   = (jsonFile.CompareTo("/home/cmsprod/json/~") != 0);

  if (gSystem->Getenv("MIT_PROD_OVERLAP")) {
    sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
    if (EOF == sscanf(overlap,"%f",&overlapCut)) {
      printf(" Overlap was not properly defined. EXIT!\n");
      return;
    }
  }
  else {
    printf(" OVERLAP file was not properly defined. EXIT!\n");
    return;
  } 

  if (gSystem->Getenv("src"))
    sprintf(src,   "%s",gSystem->Getenv("src"));
  else {
    printf(" src dir not defined. EXIT!\n");
    return;
  } 

  printf("\n Initialization worked. \n\n");

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace std;
  using namespace mithep;
  gDebugMask  = Debug::kAnalysis;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------ 
  Bool_t applyISRFilter = kFALSE;
  Bool_t applyMllGenCut = kFALSE;

  RunLumiSelectionMod *runLumiSel = new RunLumiSelectionMod;
  runLumiSel->SetAcceptMC(kTRUE);                          // Monte Carlo events are always accepted

  // if([json is not ~] and [json is not -])
  if ((jsonFile.CompareTo("/home/cmsprod/json/~") != 0) &&    //if not MC and 
      (jsonFile.CompareTo("/home/cmsprod/json/-") != 0)   ) { //if json file not absent
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  if ((jsonFile.CompareTo("/home/cmsprod/json/-") == 0)   ) {
    printf("\n WARNING -- Looking at data without JSON file: always accept.\n\n");
    runLumiSel->SetAbortIfNotAccepted(kFALSE);   // accept all events if there is no valid JSON file
  }

  printf("\n Run lumi worked. \n\n");

  //GeneratorMod
  GeneratorMod *generatorMod = new GeneratorMod;
  generatorMod->SetPrintDebug(kFALSE);
  generatorMod->SetPtLeptonMin(0.0);
  generatorMod->SetEtaLeptonMax(2.7);
  generatorMod->SetPtPhotonMin(15.0);
  generatorMod->SetEtaPhotonMax(2.7);
  generatorMod->SetPtRadPhotonMin(10.0);
  generatorMod->SetEtaRadPhotonMax(2.7);
  generatorMod->SetIsData(isData); 
  generatorMod->SetFillHist(!isData);
  if(applyMllGenCut == kTRUE){
    generatorMod->SetPdgIdCut(23);
    generatorMod->SetMassMinCut( 0.);
    generatorMod->SetMassMaxCut(50.);
  }
  generatorMod->SetApplyISRFilter(applyISRFilter);

  //------------------------------------------------------------------------------------------------
  // HLT information
  //------------------------------------------------------------------------------------------------
  HLTMod *hltMod = new HLTMod;
  if (isData) {
//   valentina:
//     hltMod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",150000,161176);
//     hltMod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",150000,161176);
//     hltMod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",161179,999999);
//     hltMod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",161179,999999);
    hltMod->AddTrigger("HLT_Mu9",132440,147119);
    hltMod->AddTrigger("HLT_Mu15",147120,9999999);
  }
  else {
    hltMod->AddTrigger("HLT_Mu9");
    hltMod->AddTrigger("HLT_Mu15");
//     hltMod->AddTrigger("HLT_Ele10_SW_L1R_v2");
  }
  hltMod->SetTrigObjsName("MyHltObjs");
  hltMod->SetAbortIfNotAccepted(kTRUE);
//   hltMod->SetPrintTable(kTRUE);
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(4);  // should be 4
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);

  //------------------------------------------------------------------------------------------------
  // Publisher Modules
  //------------------------------------------------------------------------------------------------
  PublisherMod<PFJet,Jet> *pubPFJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubPFJet->SetInputName("AKt5PFJets");
  pubPFJet->SetOutputName("PubAKt5PFJets");

  PublisherMod<Met,Met> *pubTCMet = new PublisherMod<Met,Met>("MetTCPub");
  pubTCMet->SetInputName("TCMet");
  pubTCMet->SetOutputName("PubTCMet");

  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("MetPFPub");
  pubPFMet->SetInputName("PFMet");
  pubPFMet->SetOutputName("PubPFMet");

  PublisherMod<CaloMet> *pubCaloMet = new PublisherMod<CaloMet>("CaloMetPub");
  //  pubCaloMet->SetName("CaloMetPub");
  pubCaloMet->SetInputName("CorMuonMet");
  pubCaloMet->SetOutputName("pubCaloMet");

  //------------------------------------------------------------------------------------------------
  // Apply Jet Corrections
  //------------------------------------------------------------------------------------------------
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  jetCorr->SetInputName(pubPFJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");//corrected jets are PF jets!
  jetCorr->ApplyL1FastJetCorrection(); // <<==== apply default L1 correction
  string path(string(src)+"/MitPhysics/data/");
  jetCorr->AddCorrectionFromFile(path+string("START38_V13_AK5PF_L2Relative.txt"));
  jetCorr->AddCorrectionFromFile(path+string("START38_V13_AK5PF_L3Absolute.txt"));
  if(isData) {
    jetCorr->AddCorrectionFromFile(path+string("START38_V13_AK5PF_L2L3Residual.txt"));
  }
  //------------------------------------------------------------------------------------------------
  // Apply Met Corrections
  //------------------------------------------------------------------------------------------------
  CaloMetCorrectionMod *metCaloCorr = new CaloMetCorrectionMod;
  metCaloCorr->SetInputName(pubCaloMet->GetOutputName());
  metCaloCorr->SetCorrectedJetsName(jetCorr->GetOutputName());
  metCaloCorr->SetOutputName("pubCaloCorrectedMet");

  //------------------------------------------------------------------------------------------------
  // object id and cleaning sequence
  //------------------------------------------------------------------------------------------------
  MuonIDMod           *muonId        = new MuonIDMod;
  muonId->SetClassType("Global"); // def is "Global"
  muonId->SetIDType("WWMuId");  // def is "WWMuId"
  muonId->SetIsoType("TrackCaloSliding"); // this is def
  muonId->SetApplyD0Cut(kTRUE);  // def is 1
  muonId->SetPtMin(15.0); // def is 10
  muonId->SetEtaCut(2.4); // def is 2.4

  ElectronIDMod       *electronId    = new ElectronIDMod;
  electronId->SetIDType("CustomTight"); // def is "CustomTight"
  electronId->SetIsoType("TrackJuraSliding"); // this is default
  electronId->SetApplyConversionFilterType1(kTRUE); // default is 1
  electronId->SetApplyConversionFilterType2(kFALSE); // default is 0
  electronId->SetChargeFilter(kTRUE); // def is 1
  electronId->SetApplyD0Cut(kTRUE); // def is 1
  electronId->SetNExpectedHitsInnerCut(0); // def is 999
  electronId->SetPtMin(15.0); // def is 10
  electronId->SetEtaMax(2.5); // def is 2.5

  PhotonIDMod         *photonId      = new PhotonIDMod;

  PFTauIDMod          *pfTauId        = new PFTauIDMod;
  pfTauId->SetIsHPSSel(kTRUE);
  pfTauId->SetPtMin(20.0);

  JetIDMod            *pfJetId       = new JetIDMod;
  pfJetId->SetInputName(jetCorr->GetOutputName());
  pfJetId->SetPtCut(20.0); // def is 35
  pfJetId->SetEtaMaxCut(5.0); // def is 5
  pfJetId->SetJetEEMFractionMinCut(0.01); // def is 0.01
  pfJetId->SetOutputName("GoodPFJets");

  ElectronCleaningMod *electronCleaning     = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning       = new PhotonCleaningMod;
  PFTauCleaningMod    *pfTauCleaning        = new PFTauCleaningMod;
  JetCleaningMod      *pfJetCleaning      = new JetCleaningMod;
  pfJetCleaning->SetGoodJetsName(pfJetId->GetOutputName());
  pfJetCleaning->SetCleanJetsName("CleanPFJets");

  //------------------------------------------------------------------------------------------------
  // merge modules
  //------------------------------------------------------------------------------------------------
  MergeLeptonsMod *mergeLeptonsMod = new MergeLeptonsMod;
  mergeLeptonsMod->SetMuonsName    (muonId->GetOutputName());
  mergeLeptonsMod->SetElectronsName(electronCleaning->GetOutputName());

  //------------------------------------------------------------------------------------------------
  // acceptance modules (gen -> reco)
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  // id efficiency modules (reco -> id)
  //------------------------------------------------------------------------------------------------
  //------------------------------------------------------------------------------------------------
  // analysis modules
  //------------------------------------------------------------------------------------------------
  EMUAnalysis2 *analysisModEMU = new EMUAnalysis2;
  analysisModEMU->SetTrigObjsName  (hltMod           -> GetOutputName());
  analysisModEMU->SetMuonName      (muonId           -> GetOutputName());
  analysisModEMU->SetElecName      (electronCleaning -> GetOutputName());
  analysisModEMU->SetJetName       (pfJetCleaning    -> GetOutputName());
//   analysisModEMU->SetHistNamePref  ("hXXXXX");

  //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // just in case we have run lumi sel
  if ((jsonFile.CompareTo("/home/cmsprod/json/~") != 0) &&
      (jsonFile.CompareTo("/home/cmsprod/json/-") != 0)   )
    runLumiSel->Add(generatorMod);

  generatorMod      ->Add(hltMod);
  hltMod            ->Add(goodPVFilterMod);
  goodPVFilterMod   ->Add(muonId);
  muonId            ->Add(electronId);
  electronId        ->Add(photonId);
  photonId          ->Add(pubPFJet);
  pubPFJet          ->Add(pubTCMet); 
  pubTCMet          ->Add(pubPFMet);
  pubPFMet          ->Add(pubCaloMet);
  pubCaloMet        ->Add(pfTauId);
  pfTauId           ->Add(jetCorr);
  jetCorr           ->Add(metCaloCorr);
  metCaloCorr       ->Add(pfJetId);
  pfJetId           ->Add(electronCleaning);
  electronCleaning  ->Add(photonCleaning);
  photonCleaning    ->Add(pfTauCleaning);
  pfTauCleaning     ->Add(pfJetCleaning);
  pfJetCleaning     ->Add(mergeLeptonsMod);

  mergeLeptonsMod   ->Add(analysisModEMU);

  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------

  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);
  if ((jsonFile.CompareTo("/home/cmsprod/json/~") != 0) &&
      (jsonFile.CompareTo("/home/cmsprod/json/-") != 0)   )
    ana->SetSuperModule(runLumiSel);
  else
    ana->SetSuperModule(generatorMod);
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
  //ana->AddFile("hgg-skim_p10-h110gg-gf-v26_noskim_0000_000.root");

  //------------------------------------------------------------------------------------------------
  // organize output
  //------------------------------------------------------------------------------------------------
  TString rootFile = TString(outputName);
  rootFile += TString("_") + TString(dataset) + TString("_") + TString(skim);
  if (TString(fileset) != TString(""))
    rootFile += TString("_") + TString(fileset);
  rootFile += TString(".root");
  ana->SetOutputName(rootFile.Data());
//  use current directory:
//  ana->SetOutputName((TString(outputName)+TString(".root")).Data());
  ana->SetCacheSize(64*1024*1024);

  //------------------------------------------------------------------------------------------------
  // Say what we are doing
  //------------------------------------------------------------------------------------------------
  printf("\n==== PARAMETER SUMMARY FOR THIS JOB ====\n");
  printf("\n JSON file: %s\n  and overlap cut: %f (%s)\n",jsonFile.Data(),overlapCut,overlap);
  printf("\n Rely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Skim: %s  Fileset: %s <-\n",book,dataset,skim,fileset);
  printf("\n Root output: %s\n\n",rootFile.Data());  
  printf("\n========================================\n");

  //------------------------------------------------------------------------------------------------
  // run the analysis after successful initialisation
  //------------------------------------------------------------------------------------------------
  ana->Run(!gROOT->IsBatch());

  return;
}
