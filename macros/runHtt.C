// $Id:$
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TSystem.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitAna/TreeMod/interface/HLTMod.h"
#include "MitAna/PhysicsMod/interface/PublisherMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitPhysics/Mods/interface/GeneratorMod.h"
#include "MitPhysics/Mods/interface/PDFProducerMod.h"
#include "MitPhysics/Mods/interface/HKFactorProducer.h"
#include "MitPhysics/Mods/interface/JetCorrectionMod.h"
#include "MitPhysics/Mods/interface/CaloMetCorrectionMod.h"
#include "MitPhysics/Mods/interface/MuonIDMod.h"
#include "MitPhysics/Mods/interface/ElectronIDMod.h"
#include "MitPhysics/Mods/interface/PhotonIDMod.h"
#include "MitPhysics/Mods/interface/TauIDMod.h"
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
#include "MitPhysics/Mods/interface/GoodPVFilterMod.h"
#include "MitAna/PhysicsMod/interface/RunLumiSelectionMod.h"
#include "MitPhysics/Mods/interface/EffMod.h"

#include "MitHtt/Mods/interface/HttAnalysis.h"
#include "MitHtt/Mods/interface/ZeeAnalysis.h"
#include "MitHtt/Mods/interface/ZmmAnalysis.h"
#include "MitHtt/Mods/interface/ZttAnalysis.h"

#endif
//--------------------------------------------------------------------------------------------------
void runHtt(const char *fileset    = "0000",
	    const char *skim       = "noskim",
	    const char *dataset    = "p10-pj15-v36",
	    const char *book       = "local/filefi/014",
	    const char *catalogDir = "/home/cmsprod/catalog",
	    const char *outputName = "htt",
	    int         nEvents    = 100000)
{
  printf("\n==== Enter macro  ====\n");

  //------------------------------------------------------------------------------------------------
  // some parameters get passed through the environment
  //------------------------------------------------------------------------------------------------
  char json[1024], overlap[1024];
  float overlapCut = -1;
  sprintf(json,   "%s",gSystem->Getenv("MIT_PROD_JSON"));
  TString jsonFile = TString("/home/klute/cms/root/json/") + TString(json);
  sprintf(overlap,"%s",gSystem->Getenv("MIT_PROD_OVERLAP"));
  if (EOF == sscanf(overlap,"%f",&overlapCut)) {
    printf(" Overlap was not properly defined. EXIT!\n");
    return;
  }

  //------------------------------------------------------------------------------------------------
  // some global setups
  //------------------------------------------------------------------------------------------------
  using namespace mithep;
  gDebugMask  = Debug::kGeneral;
  gDebugLevel = 1;

  //------------------------------------------------------------------------------------------------
  // set up information
  //------------------------------------------------------------------------------------------------ 
  Bool_t applyISRFilter = kFALSE;
  Bool_t applyMllGenCut = kFALSE;
  Bool_t isData         = kFALSE;

  RunLumiSelectionMod *runLumiSel = 0;
  // only slect on run-and lumisection numbers whenjson file present
  if ((jsonFile.CompareTo("/home/klute/cms/root/json/~") != 0) &&
      (jsonFile.CompareTo("/home/klute/cms/root/json/-") != 0)   ) {
    runLumiSel = new RunLumiSelectionMod;
    runLumiSel->SetAcceptMC(kTRUE);
    runLumiSel->AddJSONFile(jsonFile.Data());
  }
  
  // determine if we run on data or not # not very elegant
  if (jsonFile.CompareTo("/home/klute/cms/root/json/~") != 0)
    {
      printf("\n==== THIS IS DATA ====\n");
      isData = kTRUE;
    }

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
  // this needs some more work
  HLTMod *hltMod = new HLTMod;
  hltMod->AddTrigger("HLT_Mu9");
  hltMod->AddTrigger("HLT_Mu11");
  hltMod->AddTrigger("HLT_Mu15");
  hltMod->AddTrigger("HLT_L2Mu9");
  hltMod->AddTrigger("HLT_L2Mu11");
  hltMod->AddTrigger("HLT_DoubleMu0");
  hltMod->AddTrigger("HLT_Photon10_L1R",132440,137028);
  hltMod->AddTrigger("HLT_Photon15_Cleaned_L1R",138564,140401);
  hltMod->AddTrigger("HLT_Ele15_SW_CaloEleId_L1R",141956,999999);
  hltMod->AddTrigger("HLT_Ele15_LW_L1R");
  hltMod->AddTrigger("HLT_Photon15_L1R");
  hltMod->AddTrigger("HLT_Photon20_L1R");
  hltMod->AddTrigger("HLT_Photon20_Cleaned_L1R");
  hltMod->AddTrigger("HLT_Ele15_SW_EleId_L1R");
  hltMod->AddTrigger("HLT_Ele17_SW_EleId_L1R");
  hltMod->SetTrigObjsName("MyHltObjs");
  //------------------------------------------------------------------------------------------------
  // select events with a good primary vertex
  //------------------------------------------------------------------------------------------------
  GoodPVFilterMod *goodPVFilterMod = new GoodPVFilterMod;
  goodPVFilterMod->SetMinVertexNTracks(0);
  goodPVFilterMod->SetMinNDof(5);
  goodPVFilterMod->SetMaxAbsZ(24.0);
  goodPVFilterMod->SetMaxRho(2.0);

  //------------------------------------------------------------------------------------------------
  // Publisher Modules
  //------------------------------------------------------------------------------------------------
  PublisherMod<PFJet,Jet> *pubPFJet = new PublisherMod<PFJet,Jet>("JetPub");
  pubPFJet->SetInputName("AKt5PFJets");
  pubPFJet->SetOutputName("PubAKt5PFJets");

  PublisherMod<Met,Met> *pubMet = new PublisherMod<Met,Met>("MetPub");
  pubMet->SetInputName("TCMet");
  pubMet->SetOutputName("PubTCMet");

  PublisherMod<PFMet,Met> *pubPFMet = new PublisherMod<PFMet,Met>("MetPFPub");
  pubPFMet->SetInputName("PFMet");
  pubPFMet->SetOutputName("PubPFMet");

  PublisherMod<CaloMet> *pubCaloMet = new PublisherMod<CaloMet>;
  pubCaloMet->SetName("CaloMetPub");
  pubCaloMet->SetInputName("CorMuonMet");
  pubCaloMet->SetOutputName("pubCaloMet");

  //------------------------------------------------------------------------------------------------
  // Apply Jet Corrections
  //------------------------------------------------------------------------------------------------
  JetCorrectionMod *jetCorr = new JetCorrectionMod;
  //Check to ensure that we are using the right corrections
  jetCorr->AddCorrectionFromRelease("CondFormats/JetMETObjects/data/Spring10_L2Relative_AK5PF.txt"); 
  jetCorr->AddCorrectionFromRelease("CondFormats/JetMETObjects/data/Spring10_L3Absolute_AK5PF.txt");  
  jetCorr->SetInputName(pubPFJet->GetOutputName());
  jetCorr->SetCorrectedName("CorrectedJets");

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
  muonId->SetPtMin(15.0);
  muonId->SetIDType("Minimal"); //Check this Muon criteria
  muonId->SetIsoType("TrackCaloSliding");
  muonId->SetClassType("Global");
  //muonId->SetApplyD0Cut(kTRUE);

  ElectronIDMod       *electronId    = new ElectronIDMod;
  electronId->SetApplyTriggerMatching(kFALSE); //Check for trigger matching in the analysis module
  electronId->SetPtMin(15.0);
  electronId->SetIDType("VBTFWorkingPoint80Id");
  electronId->SetApplyConversionFilterType1(kFALSE);
  electronId->SetApplyConversionFilterType2(kTRUE);
  electronId->SetChargeFilter(kFALSE);
  electronId->SetIsoType("TrackJuraSliding");
  electronId->SetNExpectedHitsInnerCut(0);

  PhotonIDMod         *photonId      = new PhotonIDMod;

  JetIDMod            *caloJetId         = new JetIDMod;
  caloJetId->SetInputName(jetCorr->GetOutputName());
  caloJetId->SetPtCut(20.0);
  caloJetId->SetEtaMaxCut(5.0);
  caloJetId->SetJetEEMFractionMinCut(0.0);
  caloJetId->SetOutputName("GoodCaloJets");

  JetIDMod *pfJetId=new JetIDMod;
  pfJetId->SetInputName(pubPFJet->GetOutputName());
  pfJetId->SetPtCut(20.0);
  pfJetId->SetEtaMaxCut(5.0);
  pfJetId->SetJetEEMFractionMinCut(0.0);
  pfJetId->SetOutputName("GoodPFJets");

  TauIDMod            *tauId          = new TauIDMod;

  PFTauIDMod          *pfTauId        = new PFTauIDMod;
  //Check options on PFTauIDMod; what is the optimal range for signal mass?

  ElectronCleaningMod *electronCleaning     = new ElectronCleaningMod;
  PhotonCleaningMod   *photonCleaning       = new PhotonCleaningMod;
  TauCleaningMod      *tauCleaning          = new TauCleaningMod;
  PFTauCleaningMod    *pfTauCleaning        = new PFTauCleaningMod;
  JetCleaningMod      *caloJetCleaning      = new JetCleaningMod;
  caloJetCleaning->SetGoodJetsName(caloJetId->GetOutputName());
  caloJetCleaning->SetCleanJetsName("CleanCaloJets");

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
  // analyses modules
  //------------------------------------------------------------------------------------------------
  ZeeAnalysis *analysisModZee = new ZeeAnalysis;
  analysisModZee->SetTrigObjsName   (hltMod    ->GetOutputName());
  analysisModZee->SetElecName       (electronId->GetOutputName());
  analysisModZee->SetElecsFromBranch(kFALSE);
  analysisModZee->SetOverlapCut(double(overlapCut));
  analysisModZee->SetIsData(isData);

  ZmmAnalysis *analysisModZmm = new ZmmAnalysis;
  analysisModZmm->SetTrigObjsName  (hltMod  ->GetOutputName());
  analysisModZmm->SetMuonName      (muonId  ->GetOutputName());
  analysisModZmm->SetMuonsFromBranch(kFALSE);
  analysisModZmm->SetOverlapCut(double(overlapCut));
  analysisModZmm->SetIsData(isData);

  ZttAnalysis *analysisModZtt = new ZttAnalysis;
 //------------------------------------------------------------------------------------------------
  // making analysis chain
  //------------------------------------------------------------------------------------------------
  // just in case we have run lumi sel
  if ((jsonFile.CompareTo("/home/klute/cms/root/json/~") != 0) &&
      (jsonFile.CompareTo("/home/klute/cms/root/json/-") != 0)   )
    runLumiSel->Add(generatorMod);

  generatorMod->Add(hltMod);

  // high level trigger is always first
  hltMod          ->Add(goodPVFilterMod);
  goodPVFilterMod ->Add(muonId);
  // simple object id modules
  muonId          ->Add(electronId);
  electronId      ->Add(photonId);
  photonId        ->Add(tauId);
  tauId           ->Add(pfTauId); 
  pfTauId         ->Add(pubPFJet);
  pubPFJet        ->Add(pubMet); 
  pubMet          ->Add(pubPFMet);
  pubPFMet        ->Add(pubCaloMet);
  pubCaloMet      ->Add(jetCorr);
  jetCorr         ->Add(metCaloCorr);
  metCaloCorr     ->Add(caloJetId);
  caloJetId       ->Add(pfJetId);

 // cleaning modules
  pfJetId         ->Add(electronCleaning);
  electronCleaning->Add(photonCleaning);
  photonCleaning  ->Add(tauCleaning);
  tauCleaning     ->Add(pfTauCleaning);
  pfTauCleaning   ->Add(caloJetCleaning);
  caloJetCleaning ->Add(pfJetCleaning);
  pfJetCleaning   ->Add(mergeLeptonsMod);

  mergeLeptonsMod ->Add(analysisModZee);
  analysisModZee  ->Add(analysisModZmm);
  analysisModZmm  ->Add(analysisModZtt);
  //------------------------------------------------------------------------------------------------
  // setup analysis
  //------------------------------------------------------------------------------------------------

  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  ana->SetKeepHierarchy(kTRUE);
  if (nEvents >= 0)
    ana->SetProcessNEvents(nEvents);
  if ((jsonFile.CompareTo("/home/klute/cms/root/json/~") != 0) &&
      (jsonFile.CompareTo("/home/klute/cms/root/json/-") != 0)   )
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
