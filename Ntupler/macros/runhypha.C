#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitHtt/Ntupler/interface/HyphaMod.hh"
#endif

using namespace mithep;

/*
 *  useGen options:
 *  0: ignore GEN info
 *  1: H->WW->2l2nu
 *  2: Z->ll
 *  3: W->lnu
 *  4: WW->2l2nu
 */
   
//==================================================================================================
/*
 * Run on a BAMBU fileset
 *
 * Example usage:
 *   root -b -q runhypha.C+\(\"0000\",\"p11-h160ww2l-gf-v1g1-pu\",\"cern/filefi/020\",\"/home/mitprod/catalog\",0,1,-1,0\)
 *
 * Output file name has standard format: <dataset>_<fileset>_ntuple.root
 *
 */
void runhypha(
    const char   *fileset,      // "4-digit" string that labels a group of files
    const char   *dataset,      // BAMBU dataset name
    const char   *book,         // BAMBU book containing the dataset
    const char   *catalogDir,   // BAMBU catalog directory
    const Bool_t  isData,       // flag to indicate processing of collision data
    const Int_t   useGen,       // which MC process? 
    const Int_t   nevents,      // number of events to process
    const Bool_t  skipHLTFail   // skip events if no HLT accept
  )
{
  gDebugMask  = Debug::kAnalysis;  // debug message category
  gDebugLevel = 1;                 // higher level allows more messages to print
 
  char output[100];
  sprintf(output,"%s_%s_ntuple.root",dataset,fileset); 
  
  // muon kinematics
  const Double_t muPtMin  = 0;
  const Double_t muPtMax  = 7000;
  const Double_t muEtaMin = -3;
  const Double_t muEtaMax =  3;

  // electron kinematics
  const Double_t eleEtMin  = 9;
  const Double_t eleEtMax  = 7000;
  const Double_t eleEtaMin = -3;
  const Double_t eleEtaMax =  3;
  
  // jet requirements
  const Double_t jetPtMin = 10;

  // photon requirements
  const Double_t photonEtMin = 9;
      
  // good PV requirements
  const UInt_t   minNTracksFit = 0;
  const Double_t minNdof       = 4;
  const Double_t maxAbsZ       = 24;
  const Double_t maxRho        = 2;
  
  //
  // setup analysis object
  //
  Analysis *ana = new Analysis;
  ana->SetUseHLT(kTRUE);
  if(nevents>0) 
    ana->SetProcessNEvents(nevents);
  printf("\nRely on Catalog: %s\n",catalogDir);
  printf("  -> Book: %s  Dataset: %s  Fileset: %s <-\n\n",book,dataset,fileset);
  Catalog *c = new Catalog(catalogDir);
  Dataset *d = NULL;
  d = c->FindDataset(book,dataset,fileset);
  ana->AddDataset(d);
    
  //
  // setup ntupler module
  //
  HyphaMod *hymod = new HyphaMod;
  hymod->SetOutputName(output);          // output ntuple file name
  hymod->SetIsData(isData);              // toggle data specific or MC specific procedures
  hymod->SetUseGen(useGen);              // use generator info
  hymod->SetSkipIfHLTFail(skipHLTFail);  // skip to next event if no HLT accept
  hymod->SetMuPtMin(muPtMin);
  hymod->SetMuPtMax(muPtMax);
  hymod->SetMuEtaMin(muEtaMin);
  hymod->SetMuEtaMax(muEtaMax);
  hymod->SetEleEtMin(eleEtMin);
  hymod->SetEleEtMax(eleEtMax);
  hymod->SetEleEtaMin(eleEtaMin);
  hymod->SetEleEtaMax(eleEtaMax);
  hymod->SetJetPtMin(jetPtMin);
  hymod->SetPhotonEtMin(photonEtMin);
  hymod->SetMinNTracksFit(minNTracksFit);
  hymod->SetMinNdof(minNdof);
  hymod->SetMaxAbsZ(maxAbsZ);
  hymod->SetMaxRho(maxRho);

  hymod->AddTrigger("HLT_Mu8_v1",     kHLT_Mu8);
  hymod->AddTrigger("HLT_Mu8_v2",     kHLT_Mu8);
  hymod->AddTrigger("HLT_Mu8_v3",     kHLT_Mu8);
  hymod->AddTrigger("HLT_Mu9",        kHLT_Mu9);
  hymod->AddTrigger("HLT_Mu11",       kHLT_Mu9); // note: same number as mu9
  hymod->AddTrigger("HLT_Mu15_v1",    kHLT_Mu15);
  hymod->AddTrigger("HLT_Mu15_v2",    kHLT_Mu15);
  hymod->AddTrigger("HLT_Mu15_v3",    kHLT_Mu15);
  hymod->AddTrigger("HLT_Mu15_v4",    kHLT_Mu15);
  hymod->AddTrigger("HLT_Mu21_v1",    kHLT_Mu21, 24);
  hymod->AddTrigger("HLT_Mu24_v1",    kHLT_Mu24);
  hymod->AddTrigger("HLT_Mu24_v2",    kHLT_Mu24);
  hymod->AddTrigger("HLT_Mu24_v3",    kHLT_Mu24);
  hymod->AddTrigger("HLT_IsoMu17_v4", kHLT_IsoMu17);
  hymod->AddTrigger("HLT_IsoMu17_v5", kHLT_IsoMu17);
  hymod->AddTrigger("HLT_IsoMu17_v6", kHLT_IsoMu17);
  hymod->AddTrigger("HLT_IsoMu17_v8", kHLT_IsoMu17);

  hymod->AddTrigger("HLT_Ele8_v1",                                               kHLT_Ele8);
  hymod->AddTrigger("HLT_Ele8_v2",                                               kHLT_Ele8);
  hymod->AddTrigger("HLT_Ele8_v3",                                               kHLT_Ele8);
  hymod->AddTrigger("HLT_Ele8_v4",                                               kHLT_Ele8);
  hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v1",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v2",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v3",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v4",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele17_SW_L1R_v2",                                       kHLT_Ele17_SW_L1R);
  hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v1",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v2",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v3",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v4",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele22_SW_L1R_v2",                                       kHLT_Ele22_SW_L1R, 27);
  hymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
  hymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
  hymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);  
  hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
  hymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
  hymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
  hymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v4",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
  hymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL);
  hymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL);
  hymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
  hymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);

  hymod->AddTrigger("HLT_DoubleMu5_v1", kHLT_DoubleMu5, 7);
  hymod->AddTrigger("HLT_DoubleMu7_v1", kHLT_DoubleMu7);
  hymod->AddTrigger("HLT_DoubleMu7_v2", kHLT_DoubleMu7);
  hymod->AddTrigger("HLT_Mu13_Mu8_v2",  kHLT_Mu13_Mu8);
  hymod->AddTrigger("HLT_Mu5_Jet50U_v3",kHLT_Mu5_Jet50U);
  hymod->AddTrigger("HLT_Mu8_Jet40_v2", kHLT_Mu8_Jet40);
  hymod->AddTrigger("HLT_Mu8_Jet40_v3", kHLT_Mu8_Jet40);
  hymod->AddTrigger("HLT_Mu8_Jet40_v4", kHLT_Mu8_Jet40);
  hymod->AddTrigger("HLT_Mu8_Jet40_v5", kHLT_Mu8_Jet40);

  hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
  hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
  hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
  hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
  
  hymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",         kHLT_Mu17_Ele8_CaloIdL, 0,  "hltL1Mu3EG5L3Filtered17"    , kHLT_Mu17_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",   kHLT_Mu17_Ele8_CaloIdL);
  hymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",         kHLT_Mu17_Ele8_CaloIdL, 0,  "hltL1Mu3EG5L3Filtered17"    , kHLT_Mu17_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",   kHLT_Mu17_Ele8_CaloIdL);
  hymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",         kHLT_Mu17_Ele8_CaloIdL, 0,  "hltL1MuOpenEG5L3Filtered17" , kHLT_Mu17_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",   kHLT_Mu17_Ele8_CaloIdL);
  hymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v4",         kHLT_Mu17_Ele8_CaloIdL, 0,  "hltL1MuOpenEG5L3Filtered17" , kHLT_Mu17_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",   kHLT_Mu17_Ele8_CaloIdL);
  hymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",         kHLT_Mu8_Ele17_CaloIdL, 0,  "hltL1Mu3EG5L3Filtered8"     , kHLT_Mu8_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",   kHLT_Mu8_Ele17_CaloIdL);
  hymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",         kHLT_Mu8_Ele17_CaloIdL, 0,  "hltL1Mu3EG5L3Filtered8"     , kHLT_Mu8_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",   kHLT_Mu8_Ele17_CaloIdL);
  hymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",         kHLT_Mu8_Ele17_CaloIdL, 0,  "hltL1MuOpenEG5L3Filtered8"  , kHLT_Mu8_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",   kHLT_Mu8_Ele17_CaloIdL);
  hymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v4",         kHLT_Mu8_Ele17_CaloIdL, 0,  "hltL1MuOpenEG5L3Filtered8"  ,	kHLT_Mu8_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",   kHLT_Mu8_Ele17_CaloIdL);

//   hymod->AddTrigger("HLT_Mu8_Ele8_v1",                  kHLT_Mu8_Ele8);
//   hymod->AddTrigger("HLT_Mu11_Ele8_v1",                 kHLT_Mu11_Ele8);
  hymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v1",kHLT_Mu8_Photon20_CaloIdVT_IsoT);
  hymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2",kHLT_Mu8_Photon20_CaloIdVT_IsoT);
  hymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v3",kHLT_Mu8_Photon20_CaloIdVT_IsoT);
  hymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v4",kHLT_Mu8_Photon20_CaloIdVT_IsoT);
  //----------------------------------------------------------------------------------------
  // si:
  //
  //Main EMu Triggers
//   mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1Mu3EG5L3Filtered17",     kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter"  );
//   mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1Mu3EG5L3Filtered17",     kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter");
//   mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1MuOpenEG5L3Filtered17",  kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter");
//   mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v4",               kHLT_Mu17_Ele8_CaloIdL, kHLTObject_Mu17, "hltL1MuOpenEG5L3Filtered17",  kHLTObject_Ele8_CaloIdL, "hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter");
//   mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1Mu3EG5L3Filtered8",      kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter" );
//   mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1Mu3EG5L3Filtered8",      kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter");
//   mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1MuOpenEG5L3Filtered8",   kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter");
//   mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v4",               kHLT_Mu8_Ele17_CaloIdL,  kHLTObject_Mu8, "hltL1MuOpenEG5L3Filtered8",   kHLTObject_Ele17_CaloIdL, "hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter");

  hymod->AddTrigger("HLT_Jet15U_v3", kHLT_Jet15U);
  hymod->AddTrigger("HLT_Jet30U_v3", kHLT_Jet30U);
  hymod->AddTrigger("HLT_Jet30_v1",  kHLT_Jet30);
  hymod->AddTrigger("HLT_Jet30_v2",  kHLT_Jet30);
  hymod->AddTrigger("HLT_Jet30_v2",  kHLT_Jet30);
  hymod->AddTrigger("HLT_DiJetAve15U_v3",kHLT_DiJetAve15U);
  hymod->AddTrigger("HLT_DiJetAve30U_v3",kHLT_DiJetAve30U);
  hymod->AddTrigger("HLT_DiJetAve30U_v4",kHLT_DiJetAve30U);
  
  hymod->AddTrigger("HLT_Photon30_Cleaned_L1R_v1",  kHLT_Photon30_Cleaned_L1R);
  hymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v1",kHLT_Photon30_CaloIdVL_IsoL);
  hymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v2",kHLT_Photon30_CaloIdVL_IsoL);
  hymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v3",kHLT_Photon30_CaloIdVL_IsoL);
  hymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v4",kHLT_Photon30_CaloIdVL_IsoL);
  hymod->AddTrigger("HLT_Photon50_Cleaned_L1R_v1",  kHLT_Photon50_Cleaned_L1R);
  hymod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v1",kHLT_Photon50_CaloIdVL_IsoL);
  hymod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v2",kHLT_Photon50_CaloIdVL_IsoL);
  hymod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v3",kHLT_Photon50_CaloIdVL_IsoL);
          
  hymod->AddTrigger("HLT_DoubleMu0_Quarkonium_v1",  kHLT_DoubleMu0_Quarkonium);
  hymod->AddTrigger("HLT_DoubleMu3_Quarkonium_v1",  kHLT_DoubleMu3_Quarkonium);
  hymod->AddTrigger("HLT_DoubleMu3_Quarkonium_v2",  kHLT_DoubleMu3_Quarkonium);
  hymod->AddTrigger("HLT_DoubleMu3_Jpsi_v1",        kHLT_DoubleMu3_Jpsi);
  hymod->AddTrigger("HLT_DoubleMu3_Jpsi_v2",        kHLT_DoubleMu3_Jpsi);
  hymod->AddTrigger("HLT_Dimuon0_Barrel_Upsilon_v1",kHLT_Dimuon0_Barrel_Upsilon);
  hymod->AddTrigger("HLT_Dimuon6p5_Barrel_Jpsi_v1", kHLT_Dimuon6p5_Barrel_Jpsi);
  hymod->AddTrigger("HLT_Dimuon6p5_Jpsi_v1",        kHLT_Dimuon6p5_Jpsi);
  
  hymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  
  ana->AddSuperModule(hymod); 
    
  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}

// //==================================================================================================
// /*
//  * Run on a single BAMBU file (mainly for testing purposes)
//  *
//  */
// void runhypha(
//     const char *file   = "/castor/cern.ch/user/p/paus/filefi/020/p11-h160ww2l-gf-v1g1-pu/3EABE6D3-F150-E011-A1E9-00A0D1EE89E0.root",
//     const char *output = "test.root",
//     Bool_t isData      = kFALSE,
//     Int_t  useGen      = 1,
//     Int_t  nevents     = -1,
//     Bool_t skipHLTFail = kFALSE
//   )
// {
//   gDebugMask  = Debug::kAnalysis;  // debug message category
//   gDebugLevel = 1;                 // higher level allows more messages to print
  
//   // muon kinematics
//   const Double_t muPtMin  = 0;
//   const Double_t muPtMax  = 7000;
//   const Double_t muEtaMin = -3;
//   const Double_t muEtaMax =  3;

//   // electron kinematics
//   const Double_t eleEtMin  = 9;
//   const Double_t eleEtMax  = 7000;
//   const Double_t eleEtaMin = -3;
//   const Double_t eleEtaMax =  3;
  
//   // jet requirements
//   const Double_t jetPtMin = 10;

//   // photon requirements
//   const Double_t photonEtMin = 9;
      
//   // good PV requirements
//   const UInt_t   minNTracksFit = 0;
//   const Double_t minNdof       = 4;
//   const Double_t maxAbsZ       = 24;
//   const Double_t maxRho        = 2;
  
//   //
//   // setup analysis object
//   //
//   Analysis *ana = new Analysis;
//   ana->SetUseHLT(kTRUE);
//   if(nevents>0) 
//     ana->SetProcessNEvents(nevents);
//   ana->AddFile(file);
    
//   //
//   // setup ntupler module
//   //
//   HyphaMod *hymod = new HyphaMod;
//   hymod->SetOutputName(output);          // output ntuple file name
//   hymod->SetIsData(isData);              // toggle data specific or MC specific procedures
//   hymod->SetUseGen(useGen);              // use generator info
//   hymod->SetSkipIfHLTFail(skipHLTFail);  // skip to next event if no HLT accept
//   hymod->SetMuPtMin(muPtMin);
//   hymod->SetMuPtMax(muPtMax);
//   hymod->SetMuEtaMin(muEtaMin);
//   hymod->SetMuEtaMax(muEtaMax);
//   hymod->SetEleEtMin(eleEtMin);
//   hymod->SetEleEtMax(eleEtMax);
//   hymod->SetEleEtaMin(eleEtaMin);
//   hymod->SetEleEtaMax(eleEtaMax);
//   hymod->SetJetPtMin(jetPtMin);
//   hymod->SetPhotonEtMin(photonEtMin);
//   hymod->SetMinNTracksFit(minNTracksFit);
//   hymod->SetMinNdof(minNdof);
//   hymod->SetMaxAbsZ(maxAbsZ);
//   hymod->SetMaxRho(maxRho);

//   hymod->AddTrigger("HLT_Mu8_v1",     kHLT_Mu8);
//   hymod->AddTrigger("HLT_Mu8_v2",     kHLT_Mu8);
//   hymod->AddTrigger("HLT_Mu9",        kHLT_Mu9);
//   hymod->AddTrigger("HLT_Mu15_v1",    kHLT_Mu15);
//   hymod->AddTrigger("HLT_Mu15_v2",    kHLT_Mu15);
//   hymod->AddTrigger("HLT_Mu15_v3",    kHLT_Mu15);
//   hymod->AddTrigger("HLT_Mu21_v1",    kHLT_Mu21, 24);
//   hymod->AddTrigger("HLT_Mu24_v1",    kHLT_Mu24);
//   hymod->AddTrigger("HLT_Mu24_v2",    kHLT_Mu24);
//   hymod->AddTrigger("HLT_IsoMu17_v4", kHLT_IsoMu17);
//   hymod->AddTrigger("HLT_IsoMu17_v5", kHLT_IsoMu17);
//   hymod->AddTrigger("HLT_IsoMu17_v6", kHLT_IsoMu17);

//   hymod->AddTrigger("HLT_Ele8_v1",                                               kHLT_Ele8);
//   hymod->AddTrigger("HLT_Ele8_v2",                                               kHLT_Ele8);
//   hymod->AddTrigger("HLT_Ele8_v3",                                               kHLT_Ele8);
//   hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v1",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v2",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v3",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele17_SW_L1R_v2",                                       kHLT_Ele17_SW_L1R);
//   hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v1",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v2",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v3",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele22_SW_L1R_v2",                                       kHLT_Ele22_SW_L1R, 27);
//   hymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
//   hymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
//   hymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
//   hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
//   hymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
//   hymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
//   hymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL);
//   hymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL);
//   hymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
//   hymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);

//   hymod->AddTrigger("HLT_DoubleMu5_v1", kHLT_DoubleMu5, 7);
//   hymod->AddTrigger("HLT_DoubleMu7_v1", kHLT_DoubleMu7);
//   hymod->AddTrigger("HLT_DoubleMu7_v2", kHLT_DoubleMu7);
//   hymod->AddTrigger("HLT_Mu5_Jet50U_v3",kHLT_Mu5_Jet50U);
//   hymod->AddTrigger("HLT_Mu8_Jet40_v2", kHLT_Mu8_Jet40);
//   hymod->AddTrigger("HLT_Mu8_Jet40_v3", kHLT_Mu8_Jet40);
//   hymod->AddTrigger("HLT_Mu8_Jet40_v4", kHLT_Mu8_Jet40);

//   hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
//   hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
//   hymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
  
//   hymod->AddTrigger("HLT_Mu11_Ele8_v1",                 kHLT_Mu11_Ele8);
//   hymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",         kHLT_Mu17_Ele8_CaloIdL);
//   hymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",         kHLT_Mu17_Ele8_CaloIdL);
//   hymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",         kHLT_Mu17_Ele8_CaloIdL);
//   hymod->AddTrigger("HLT_Mu8_Ele8_v1",                  kHLT_Mu8_Ele8);
//   hymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",         kHLT_Mu8_Ele17_CaloIdL);
//   hymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",         kHLT_Mu8_Ele17_CaloIdL);
//   hymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",         kHLT_Mu8_Ele17_CaloIdL);
//   hymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v1",kHLT_Mu8_Photon20_CaloIdVT_IsoT);
//   hymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2",kHLT_Mu8_Photon20_CaloIdVT_IsoT);
//   hymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v3",kHLT_Mu8_Photon20_CaloIdVT_IsoT);

//   hymod->AddTrigger("HLT_Jet15U_v3", kHLT_Jet15U);
//   hymod->AddTrigger("HLT_Jet30U_v3", kHLT_Jet30U);
//   hymod->AddTrigger("HLT_Jet30_v1",  kHLT_Jet30);
//   hymod->AddTrigger("HLT_Jet30_v2",  kHLT_Jet30);
//   hymod->AddTrigger("HLT_DiJetAve15U_v3",kHLT_DiJetAve15U);
//   hymod->AddTrigger("HLT_DiJetAve30U_v3",kHLT_DiJetAve30U);
//   hymod->AddTrigger("HLT_DiJetAve30U_v4",kHLT_DiJetAve30U);
  
//   hymod->AddTrigger("HLT_Photon30_Cleaned_L1R_v1",  kHLT_Photon30_Cleaned_L1R);
//   hymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v1",kHLT_Photon30_CaloIdVL_IsoL);
//   hymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v2",kHLT_Photon30_CaloIdVL_IsoL);
//   hymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v3",kHLT_Photon30_CaloIdVL_IsoL);
//   hymod->AddTrigger("HLT_Photon50_Cleaned_L1R_v1",  kHLT_Photon50_Cleaned_L1R);
//   hymod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v1",kHLT_Photon50_CaloIdVL_IsoL);
//   hymod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v2",kHLT_Photon50_CaloIdVL_IsoL);
          
//   hymod->AddTrigger("HLT_DoubleMu0_Quarkonium_v1",  kHLT_DoubleMu0_Quarkonium);
//   hymod->AddTrigger("HLT_DoubleMu3_Quarkonium_v1",  kHLT_DoubleMu3_Quarkonium);
//   hymod->AddTrigger("HLT_DoubleMu3_Quarkonium_v2",  kHLT_DoubleMu3_Quarkonium);
//   hymod->AddTrigger("HLT_DoubleMu3_Jpsi_v1",        kHLT_DoubleMu3_Jpsi);
//   hymod->AddTrigger("HLT_DoubleMu3_Jpsi_v2",        kHLT_DoubleMu3_Jpsi);
//   hymod->AddTrigger("HLT_Dimuon0_Barrel_Upsilon_v1",kHLT_Dimuon0_Barrel_Upsilon);
//   hymod->AddTrigger("HLT_Dimuon6p5_Barrel_Jpsi_v1", kHLT_Dimuon6p5_Barrel_Jpsi);
//   hymod->AddTrigger("HLT_Dimuon6p5_Jpsi_v1",        kHLT_Dimuon6p5_Jpsi); 

//   hymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  
//   ana->AddSuperModule(hymod); 
  
//   //
//   // run analysis after successful initialisation
//   //
//   ana->Run(!gROOT->IsBatch());
// }
