#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitHtt/Ntupler/interface/HttNtuplerMod.hh"
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
 *   root -b -l -q runHttNtupler.C+\(\"0000\",\"r11a-del-pr-v4\",\"t2mit/filefi/022\",\"/home/cmsprod/catalog\",1,0,100,1\)
 *   root -b -l -q runHttNtupler.C+\(\"0000\",\"s11-h100tt-gf-v1g1-pu\",\"t2mit/filefi/022\",\"/home/cmsprod/catalog\",0,1,10,0\)
 *   root -b -l -q runHttNtupler.C+\(\"0000\",\"r11a-mueg-m10-v1\",\"local/filefi/021\",\"/home/cmsprod/catalog\",1,0,-1,1\)
 *   root -b -l -q runHttNtupler.C+\(\"0000\",\"r11a-mueg-m10-v1\",\"local/filefi/021\",\"/home/cmsprod/catalog\",1,0,100,1,\"foo.json\"\)|grep -v '^\*'         
 * Output file name has standard format: <dataset>_<fileset>_ntuple.root
 *
 */
void runHttNtupler(
    const char   *fileset,      // "4-digit" string that labels a group of files
    const char   *dataset,      // BAMBU dataset name
    const char   *book,         // BAMBU book containing the dataset
    const char   *catalogDir,   // BAMBU catalog directory
    const Bool_t  isData,       // flag to indicate processing of collision data
    const Int_t   useGen,       // which MC process? 
    const Int_t   nevents,      // number of events to process
    const Bool_t  skipHLTFail,  // skip events if no HLT accept
    const char   *json=""       // file with certified runlumis
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
  HttNtuplerMod *mymod = new HttNtuplerMod;
  mymod->SetOutputName(output);          // output ntuple file name
  mymod->SetIsData(isData);              // toggle data specific or MC specific procedures
  mymod->SetUseGen(useGen);              // use generator info
  mymod->SetSkipIfHLTFail(skipHLTFail);  // skip to next event if no HLT accept
  mymod->SetMuPtMin(muPtMin);
  mymod->SetMuPtMax(muPtMax);
  mymod->SetMuEtaMin(muEtaMin);
  mymod->SetMuEtaMax(muEtaMax);
  mymod->SetEleEtMin(eleEtMin);
  mymod->SetEleEtMax(eleEtMax);
  mymod->SetEleEtaMin(eleEtaMin);
  mymod->SetEleEtaMax(eleEtaMax);
  mymod->SetJetPtMin(jetPtMin);
  mymod->SetPhotonEtMin(photonEtMin);
  mymod->SetMinNTracksFit(minNTracksFit);
  mymod->SetMinNdof(minNdof);
  mymod->SetMaxAbsZ(maxAbsZ);
  mymod->SetMaxRho(maxRho);


  // Jet corrections
  TString path("/home/vdutta/cms/cmssw/023/CMSSW_4_2_4_patch1");
  path += "/src/MitPhysics/data/";

  mymod->AddJetCorr(path   + "START42_V12_AK5PF_L1FastJet.txt");
  mymod->AddJetCorr(path   + "START42_V12_AK5PF_L2Relative.txt");
  mymod->AddJetCorr(path   + "START42_V12_AK5PF_L3Absolute.txt");
  if(isData)		      				  
    mymod->AddJetCorr(path + "START42_V12_AK5PF_L2L3Residual.txt");

  if(TString(json).Length() > 0)
    mymod->AddJSON(json);
  
  //
  // MuEG
  // 
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",kHLT_Mu17_Ele8_CaloIdL,"hltL1Mu3EG5L3Filtered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",kHLT_Mu17_Ele8_CaloIdL,"hltL1Mu3EG5L3Filtered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",kHLT_Mu17_Ele8_CaloIdL,"hltL1MuOpenEG5L3Filtered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v4",kHLT_Mu17_Ele8_CaloIdL,"hltL1MuOpenEG5L3Filtered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v5",kHLT_Mu17_Ele8_CaloIdL,"hltL1MuOpenEG5L3Filtered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v6",kHLT_Mu17_Ele8_CaloIdL,"hltL1MuOpenEG5L3Filtered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v8",kHLT_Mu17_Ele8_CaloIdL,"hltL1Mu7EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v1",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"hltL1Mu3EG5L3Filtered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj,0,"hltMu17Ele8CaloIdTPixelMatchFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v3",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"hltL1Mu7EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj,0,"hltMu17Ele8CaloIdTPixelMatchFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",kHLT_Mu8_Ele17_CaloIdL,"hltL1Mu3EG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",kHLT_Mu8_Ele17_CaloIdL,"hltL1Mu3EG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v4",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v5",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v6",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v8",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v1",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2",kHLT_Mu8_Photon20_CaloIdVT_IsoT,"hltSingleMu8EG5L3Filtered8",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v3",kHLT_Mu8_Photon20_CaloIdVT_IsoT,"hltSingleMu8EG5L3Filtered8",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v4",kHLT_Mu8_Photon20_CaloIdVT_IsoT,"hltSingleMu8EG5L3Filtered8",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v5",kHLT_Mu8_Photon20_CaloIdVT_IsoT,"hltSingleMu8EG5L3Filtered8",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v6",kHLT_Mu8_Photon20_CaloIdVT_IsoT,"hltSingleMu8EG5L3Filtered8",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v2",kHLT_Mu15_Photon20_CaloIdL,"hltL1Mu3EG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v3",kHLT_Mu15_Photon20_CaloIdL,"hltL1Mu3EG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v4",kHLT_Mu15_Photon20_CaloIdL,"hltL1MuOpenEG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v5",kHLT_Mu15_Photon20_CaloIdL,"hltL1MuOpenEG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v6",kHLT_Mu15_Photon20_CaloIdL,"hltL1MuOpenEG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v7",kHLT_Mu15_Photon20_CaloIdL,"hltL1MuOpenEG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);

  //
  // DoubleMu
  //
  mymod->AddTrigger("HLT_DoubleMu7_v1",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj); 
  mymod->AddTrigger("HLT_DoubleMu7_v2",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);
  mymod->AddTrigger("HLT_DoubleMu7_v3",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);
  mymod->AddTrigger("HLT_Mu13_Mu8_v1",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu13_Mu8_v2",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu13_Mu8_v3",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu13_Mu8_v4",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);  
  mymod->AddTrigger("HLT_Mu13_Mu8_v6",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v1",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v2",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v3",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v4",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v6",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu8_Jet40_v2",kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8",kHLT_Mu8_Jet40_MuObj,0,"hltJet40",kHLT_Mu8_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Mu8_Jet40_v3",kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8",kHLT_Mu8_Jet40_MuObj,0,"hltJet40",kHLT_Mu8_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Mu8_Jet40_v4",kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8",kHLT_Mu8_Jet40_MuObj,0,"hltJet40",kHLT_Mu8_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Mu8_Jet40_v5",kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8",kHLT_Mu8_Jet40_MuObj,0,"hltJet40",kHLT_Mu8_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Mu8_Jet40_v6",kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8",kHLT_Mu8_Jet40_MuObj,0,"hltJet40",kHLT_Mu8_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Mu8_Jet40_v7",kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8",kHLT_Mu8_Jet40_MuObj,0,"hltJet40",kHLT_Mu8_Jet40_JetObj,0);

  //
  // SingleMu
  //
  mymod->AddTrigger("HLT_Mu8_v1",kHLT_Mu8,"hltSingleMu8L3Filtered8",kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v2",kHLT_Mu8,"hltSingleMu8L3Filtered8",kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v3",kHLT_Mu8,"hltSingleMu8L3Filtered8",kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v4",kHLT_Mu8,"hltSingleMu8L3Filtered8",kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v5",kHLT_Mu8,"hltSingleMu8L3Filtered8",kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu9",kHLT_Mu9,"",kHLT_Mu9_MuObj,15);
  mymod->AddTrigger("HLT_Mu12_v1",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);
  mymod->AddTrigger("HLT_Mu12_v2",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);
  mymod->AddTrigger("HLT_Mu12_v3",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);
  mymod->AddTrigger("HLT_Mu12_v4",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);
  mymod->AddTrigger("HLT_Mu12_v5",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);
  mymod->AddTrigger("HLT_Mu15_v1",kHLT_Mu15,"",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v2",kHLT_Mu15,"hltL3Muon15",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v3",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v4",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v5",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v6",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu24_v1",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  mymod->AddTrigger("HLT_Mu24_v2",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  mymod->AddTrigger("HLT_Mu24_v3",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  mymod->AddTrigger("HLT_Mu24_v4",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  mymod->AddTrigger("HLT_Mu24_v5",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  mymod->AddTrigger("HLT_Mu30_v1",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v2",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v3",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v4",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v5",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_IsoMu17_v5",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);
  mymod->AddTrigger("HLT_IsoMu17_v6",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);
  mymod->AddTrigger("HLT_IsoMu17_v8",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v1",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v2",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v4",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v5",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v6",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v7",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v8",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);


  //
  // DoubleElectron
  //
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v4",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v5",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v6",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v2",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v3",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v4",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v1",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17,"hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj,0,"hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17,"hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj,0,"hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v3",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17,"hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj,0,"hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter",kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v1",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj,0,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v2",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj,0,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v3",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj,0,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele8_v1",kHLT_Ele8,"hltEle8PixelMatchFilter",kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v2",kHLT_Ele8,"hltEle8PixelMatchFilter",kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v3",kHLT_Ele8,"hltEle8PixelMatchFilter",kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v4",kHLT_Ele8,"hltEle8PixelMatchFilter",kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v5",kHLT_Ele8,"hltEle8PixelMatchFilter",kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v6",kHLT_Ele8,"hltEle8PixelMatchFilter",kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v1",kHLT_Ele8_CaloIdL_TrkIdVL,"hltEle8CaloIdLTrkIdVLDphiFilter",kHLT_Ele8_CaloIdL_TrkIdVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v2",kHLT_Ele8_CaloIdL_TrkIdVL,"hltEle8CaloIdLTrkIdVLDphiFilter",kHLT_Ele8_CaloIdL_TrkIdVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v3",kHLT_Ele8_CaloIdL_TrkIdVL,"hltEle8CaloIdLTrkIdVLDphiFilter",kHLT_Ele8_CaloIdL_TrkIdVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v4",kHLT_Ele8_CaloIdL_TrkIdVL,"hltEle8CaloIdLTrkIdVLDphiFilter",kHLT_Ele8_CaloIdL_TrkIdVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v5",kHLT_Ele8_CaloIdL_TrkIdVL,"hltEle8CaloIdLTrkIdVLDphiFilter",kHLT_Ele8_CaloIdL_TrkIdVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_TrkIdVL_v6",kHLT_Ele8_CaloIdL_TrkIdVL,"hltEle8CaloIdLTrkIdVLDphiFilter",kHLT_Ele8_CaloIdL_TrkIdVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v1",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v2",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v3",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v4",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v5",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v5",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3",kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4",kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v1",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v2",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v3",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v4",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v5",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v6",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj,0,"hltJet40Ele8CaloIdLCaloIsoVLRemoved",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj,0,"hltJet40Ele8CaloIdLCaloIsoVLRemoved",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj,0,"hltJet40Ele8CaloIdLCaloIsoVLRemoved",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj,0,"hltJet40Ele8CaloIdLCaloIsoVLRemoved",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v5",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj,0,"hltJet40Ele8CaloIdLCaloIsoVLRemoved",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v6",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj,0,"hltJet40Ele8CaloIdLCaloIsoVLRemoved",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj,0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj,0,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj,0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj,0,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj,0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj,0,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj,0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj,0,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj,0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj,0,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj,0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v6",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj,0,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj,0);

  //
  // SingleElectron
  //
  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);

  
  mymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  
  ana->AddSuperModule(mymod); 
    
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
// void runHttNtupler(
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
//   HyphaMod *mymod = new HyphaMod;
//   mymod->SetOutputName(output);          // output ntuple file name
//   mymod->SetIsData(isData);              // toggle data specific or MC specific procedures
//   mymod->SetUseGen(useGen);              // use generator info
//   mymod->SetSkipIfHLTFail(skipHLTFail);  // skip to next event if no HLT accept
//   mymod->SetMuPtMin(muPtMin);
//   mymod->SetMuPtMax(muPtMax);
//   mymod->SetMuEtaMin(muEtaMin);
//   mymod->SetMuEtaMax(muEtaMax);
//   mymod->SetEleEtMin(eleEtMin);
//   mymod->SetEleEtMax(eleEtMax);
//   mymod->SetEleEtaMin(eleEtaMin);
//   mymod->SetEleEtaMax(eleEtaMax);
//   mymod->SetJetPtMin(jetPtMin);
//   mymod->SetPhotonEtMin(photonEtMin);
//   mymod->SetMinNTracksFit(minNTracksFit);
//   mymod->SetMinNdof(minNdof);
//   mymod->SetMaxAbsZ(maxAbsZ);
//   mymod->SetMaxRho(maxRho);

//   mymod->AddTrigger("HLT_Mu8_v1",     kHLT_Mu8);
//   mymod->AddTrigger("HLT_Mu8_v2",     kHLT_Mu8);
//   mymod->AddTrigger("HLT_Mu9",        kHLT_Mu9);
//   mymod->AddTrigger("HLT_Mu15_v1",    kHLT_Mu15);
//   mymod->AddTrigger("HLT_Mu15_v2",    kHLT_Mu15);
//   mymod->AddTrigger("HLT_Mu15_v3",    kHLT_Mu15);
//   mymod->AddTrigger("HLT_Mu21_v1",    kHLT_Mu21, 24);
//   mymod->AddTrigger("HLT_Mu24_v1",    kHLT_Mu24);
//   mymod->AddTrigger("HLT_Mu24_v2",    kHLT_Mu24);
//   mymod->AddTrigger("HLT_IsoMu17_v4", kHLT_IsoMu17);
//   mymod->AddTrigger("HLT_IsoMu17_v5", kHLT_IsoMu17);
//   mymod->AddTrigger("HLT_IsoMu17_v6", kHLT_IsoMu17);

//   mymod->AddTrigger("HLT_Ele8_v1",                                               kHLT_Ele8);
//   mymod->AddTrigger("HLT_Ele8_v2",                                               kHLT_Ele8);
//   mymod->AddTrigger("HLT_Ele8_v3",                                               kHLT_Ele8);
//   mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v1",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v2",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v3",                             kHLT_Ele8_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele17_SW_L1R_v2",                                       kHLT_Ele17_SW_L1R);
//   mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v1",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v2",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v3",                            kHLT_Ele17_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele22_SW_L1R_v2",                                       kHLT_Ele22_SW_L1R, 27);
//   mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
//   mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
//   mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",             kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT);
//   mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",     kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
//   mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
//   mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30);
//   mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL);
//   mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL);
//   mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);
//   mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3",      kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL);

//   mymod->AddTrigger("HLT_DoubleMu5_v1", kHLT_DoubleMu5, 7);
//   mymod->AddTrigger("HLT_DoubleMu7_v1", kHLT_DoubleMu7);
//   mymod->AddTrigger("HLT_DoubleMu7_v2", kHLT_DoubleMu7);
//   mymod->AddTrigger("HLT_Mu5_Jet50U_v3",kHLT_Mu5_Jet50U);
//   mymod->AddTrigger("HLT_Mu8_Jet40_v2", kHLT_Mu8_Jet40);
//   mymod->AddTrigger("HLT_Mu8_Jet40_v3", kHLT_Mu8_Jet40);
//   mymod->AddTrigger("HLT_Mu8_Jet40_v4", kHLT_Mu8_Jet40);

//   mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
//   mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
//   mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40);
  
//   mymod->AddTrigger("HLT_Mu11_Ele8_v1",                 kHLT_Mu11_Ele8);
//   mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v1",         kHLT_Mu17_Ele8_CaloIdL);
//   mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v2",         kHLT_Mu17_Ele8_CaloIdL);
//   mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v3",         kHLT_Mu17_Ele8_CaloIdL);
//   mymod->AddTrigger("HLT_Mu8_Ele8_v1",                  kHLT_Mu8_Ele8);
//   mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",         kHLT_Mu8_Ele17_CaloIdL);
//   mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",         kHLT_Mu8_Ele17_CaloIdL);
//   mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",         kHLT_Mu8_Ele17_CaloIdL);
//   mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v1",kHLT_Mu8_Photon20_CaloIdVT_IsoT);
//   mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2",kHLT_Mu8_Photon20_CaloIdVT_IsoT);
//   mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v3",kHLT_Mu8_Photon20_CaloIdVT_IsoT);

//   mymod->AddTrigger("HLT_Jet15U_v3", kHLT_Jet15U);
//   mymod->AddTrigger("HLT_Jet30U_v3", kHLT_Jet30U);
//   mymod->AddTrigger("HLT_Jet30_v1",  kHLT_Jet30);
//   mymod->AddTrigger("HLT_Jet30_v2",  kHLT_Jet30);
//   mymod->AddTrigger("HLT_DiJetAve15U_v3",kHLT_DiJetAve15U);
//   mymod->AddTrigger("HLT_DiJetAve30U_v3",kHLT_DiJetAve30U);
//   mymod->AddTrigger("HLT_DiJetAve30U_v4",kHLT_DiJetAve30U);
  
//   mymod->AddTrigger("HLT_Photon30_Cleaned_L1R_v1",  kHLT_Photon30_Cleaned_L1R);
//   mymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v1",kHLT_Photon30_CaloIdVL_IsoL);
//   mymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v2",kHLT_Photon30_CaloIdVL_IsoL);
//   mymod->AddTrigger("HLT_Photon30_CaloIdVL_IsoL_v3",kHLT_Photon30_CaloIdVL_IsoL);
//   mymod->AddTrigger("HLT_Photon50_Cleaned_L1R_v1",  kHLT_Photon50_Cleaned_L1R);
//   mymod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v1",kHLT_Photon50_CaloIdVL_IsoL);
//   mymod->AddTrigger("HLT_Photon50_CaloIdVL_IsoL_v2",kHLT_Photon50_CaloIdVL_IsoL);
          
//   mymod->AddTrigger("HLT_DoubleMu0_Quarkonium_v1",  kHLT_DoubleMu0_Quarkonium);
//   mymod->AddTrigger("HLT_DoubleMu3_Quarkonium_v1",  kHLT_DoubleMu3_Quarkonium);
//   mymod->AddTrigger("HLT_DoubleMu3_Quarkonium_v2",  kHLT_DoubleMu3_Quarkonium);
//   mymod->AddTrigger("HLT_DoubleMu3_Jpsi_v1",        kHLT_DoubleMu3_Jpsi);
//   mymod->AddTrigger("HLT_DoubleMu3_Jpsi_v2",        kHLT_DoubleMu3_Jpsi);
//   mymod->AddTrigger("HLT_Dimuon0_Barrel_Upsilon_v1",kHLT_Dimuon0_Barrel_Upsilon);
//   mymod->AddTrigger("HLT_Dimuon6p5_Barrel_Jpsi_v1", kHLT_Dimuon6p5_Barrel_Jpsi);
//   mymod->AddTrigger("HLT_Dimuon6p5_Jpsi_v1",        kHLT_Dimuon6p5_Jpsi); 

//   mymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  
//   ana->AddSuperModule(mymod); 
  
//   //
//   // run analysis after successful initialisation
//   //
//   ana->Run(!gROOT->IsBatch());
// }
