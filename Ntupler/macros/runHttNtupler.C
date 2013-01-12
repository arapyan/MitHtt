#if !defined(__CINT__) || defined(__MAKECINT__)
#include "TROOT.h"
#include "cstdlib"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/Catalog/interface/Catalog.h"
#include "MitAna/TreeMod/interface/Analysis.h"
#include "MitHtt/Ntupler/interface/HttNtupler.h"
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
 *   root -b -l -q runHttNtupler.C+\(\"0000\",\"s11-h120tt-gf-v11-pu\",\"cern/filefi/025\",\"/home/mitprod/catalog\",0,1,100,0\)
 * Output file name has standard format: <dataset>_<fileset>_ntuple.root
 *
 */
void runHttNtupler(
    const char   *fileset,      // "4-digit" string that labels a group of files
    const char   *dataset,      // BAMBU dataset name
    const char   *book,         // BAMBU book containing the dataset
    const char   *catalogDir,   // BAMBU catalog directory
    const Int_t  isData,       // flag to indicate processing of collision data
    const Int_t   useGen,       // which MC process? 
    const Int_t   nevents,      // number of events to process
    const Bool_t  skipHLTFail,  // skip events if no HLT accept
    const Int_t  is2012,         // 2012 or 2011 samples
    const char   *json=""       // file with certified runlumis
  )
{
  gDebugMask  = Debug::kAnalysis;  // debug message category
  gDebugLevel = 1;                 // higher level allows more messages to print
 
  char output[100];
  sprintf(output,"%s_%s_ntuple.root",dataset,fileset); 
  
  // muon kinematics
  const Double_t muPtMin  = 3.0; //3  
  const Double_t muPtMax  = 7000;
  const Double_t muEtaMin = -2.4;
  const Double_t muEtaMax =  2.4;

  // electron kinematics
  const Double_t eleEtMin  = 7.0; //7
  const Double_t eleEtMax  = 7000;
  const Double_t eleEtaMin = -2.7;
  const Double_t eleEtaMax =  2.7;
  
  //tau kinematics
  const Double_t  tauPtMin = 18; //18
  const Double_t  tauEtaMax = 2.5;//2.5
  
  // jet requirements
  const Double_t jetPtMin = 18; //30

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
  //ana->AddFile("/castor/cern.ch/user/p/paus/filefi/029/r12c-tau-pr-v2/F8ED0695-7BF1-E111-86F5-002481E94C7E.root");
  //
  // setup ntupler module
  //
  HttNtupler *mymod = new HttNtupler;
  mymod->SetOutputName(output);          // output ntuple file name
  mymod->SetIsData(isData);              // toggle data specific or MC specific procedures
  mymod->SetUseGen(useGen);              // use generator info
  mymod->Set2012(is2012);               // 2012 samples
  mymod->SetSkipIfHLTFail(skipHLTFail);  // skip to next event if no HLT accept
  mymod->SetMuPtMin(muPtMin);
  mymod->SetMuPtMax(muPtMax);
  mymod->SetMuEtaMin(muEtaMin);
  mymod->SetMuEtaMax(muEtaMax);
  mymod->SetEleEtMin(eleEtMin);
  mymod->SetEleEtMax(eleEtMax);
  mymod->SetEleEtaMin(eleEtaMin);
  mymod->SetEleEtaMax(eleEtaMax);
  mymod->SetTauPtMin(tauPtMin);
  mymod->SetTauEtaMax(tauEtaMax);
  mymod->SetJetPtMin(jetPtMin);
  mymod->SetPhotonEtMin(photonEtMin);
  mymod->SetMinNTracksFit(minNTracksFit);
  mymod->SetMinNdof(minNdof);
  mymod->SetMaxAbsZ(maxAbsZ);
  mymod->SetMaxRho(maxRho);

  
  // Jet corrections
  char* PATH = getenv("CMSSW_BASE"); assert(PATH);
  TString path(TString::Format("%s/src/MitPhysics/data/", PATH));
  if(is2012)
    {
      if(isData || useGen==ESampleType::kEmbed){
	mymod->AddJetCorrNew(path + "GR_P_V42_AN3_L1FastJet_AK5PF.txt");
	mymod->AddJetCorrNew(path + "GR_P_V42_AN3_L2Relative_AK5PF.txt");
	mymod->AddJetCorrNew(path + "GR_P_V42_AN3_L3Absolute_AK5PF.txt");
	mymod->AddJetCorrNew(path + "GR_P_V42_AN3_L2L3Residual_AK5PF.txt");
	mymod->AddJetCorr(path + "GR_P_V41_AN2_L1FastJet_AK5PF.txt");
	mymod->AddJetCorr(path + "GR_P_V41_AN2_L2Relative_AK5PF.txt");
	mymod->AddJetCorr(path + "GR_P_V41_AN2_L3Absolute_AK5PF.txt");
	mymod->AddJetCorr(path + "GR_P_V41_AN2_L2L3Residual_AK5PF.txt");
      }
      else
	{
	  mymod->AddJetCorrNew(path + "START53_V15_L1FastJet_AK5PF.txt"   );
	  mymod->AddJetCorrNew(path + "START53_V15_L2Relative_AK5PF.txt"  );
	  mymod->AddJetCorrNew(path + "START53_V15_L3Absolute_AK5PF.txt"  );
	  mymod->AddJetCorr(path + "START53_V7F_L1FastJet_AK5PF.txt"   );
	  mymod->AddJetCorr(path + "START53_V7F_L2Relative_AK5PF.txt"  );
	  mymod->AddJetCorr(path + "START53_V7F_L3Absolute_AK5PF.txt"  );
	}
    } 
  else
    {
      mymod->AddJetCorr(path + "START42_V17_AK5PF_L1FastJet.txt"   );
      mymod->AddJetCorr(path + "START42_V17_AK5PF_L2Relative.txt"  );
      mymod->AddJetCorr(path + "START42_V17_AK5PF_L3Absolute.txt"  );
      if(isData || useGen==ESampleType::kEmbed){
	mymod->AddJetCorr(path + "GR_R_42_V23_AK5PF_L2L3Residual.txt");
      }
    }
  

  if(TString(json).Length() > 0){
    mymod->AddJSON(json);
  }
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
  //mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v9",kHLT_Mu17_Ele8_CaloIdL,"hltL1Mu7EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);
  //mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v12",kHLT_Mu17_Ele8_CaloIdL,"hltL1Mu7EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);//auto
  //mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdL_v13",kHLT_Mu17_Ele8_CaloIdL,"hltL1Mu7EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu17Ele8PixelMatchFilter",kHLT_Mu17_Ele8_CaloIdL_EGObj,0);//auto
  //mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v1",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"hltL1Mu3EG5L3Filtered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj,0,"hltMu17Ele8CaloIdTPixelMatchFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj,0);
  //mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v3",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"hltL1Mu7EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj,0,"hltMu17Ele8CaloIdTPixelMatchFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v4",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"hltL1Mu12EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj,0,"hltMu17Ele8CaloIdTPixelMatchFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v7",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"hltL1Mu12EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj,0,"hltMu17Ele8CaloIdTPixelMatchFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_v8",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL,"hltL1Mu12EG5L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj,0,"hltMu17Ele8CaloIdTPixelMatchFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1Mu12EG7L3MuFiltered17",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu17Ele8CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v1",kHLT_Mu8_Ele17_CaloIdL,"hltL1Mu3EG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v2",kHLT_Mu8_Ele17_CaloIdL,"hltL1Mu3EG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v3",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v4",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v5",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v6",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v12",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v13",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v8",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdL_v9",kHLT_Mu8_Ele17_CaloIdL,"hltL1MuOpenEG5L3Filtered8",kHLT_Mu8_Ele17_CaloIdL_MuObj,0,"hltL1NonIsoHLTNonIsoMu8Ele17PixelMatchFilter",kHLT_Mu8_Ele17_CaloIdL_EGObj,0);//auto
  //mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v1",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v3",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v4",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v7",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_v8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLPixelMatchFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltL1sL1Mu3p5EG12ORL1MuOpenEG12L3Filtered8",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj,0,"hltMu8Ele17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj,0);
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
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v10",kHLT_Mu15_Photon20_CaloIdL,"hltL1MuOpenEG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v13",kHLT_Mu15_Photon20_CaloIdL,"hltL1MuOpenEG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v14",kHLT_Mu15_Photon20_CaloIdL,"hltL1MuOpenEG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_Photon20_CaloIdL_v9",kHLT_Mu15_Photon20_CaloIdL,"hltL1MuOpenEG5L3Filtered15",kHLT_Mu15_Photon20_CaloIdL_MuObj,0,"hltMu15Photon20CaloIdLHEFilter",kHLT_Mu15_Photon20_CaloIdL_EGObj,0);//auto

  //
  // DoubleMu
  //
  mymod->AddTrigger("HLT_DoubleMu7_v1",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj); 
  mymod->AddTrigger("HLT_DoubleMu7_v2",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);
  //mymod->AddTrigger("HLT_DoubleMu7_v3",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);
  //mymod->AddTrigger("HLT_DoubleMu7_v11",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);//auto
  //mymod->AddTrigger("HLT_DoubleMu7_v12",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);//auto
  //mymod->AddTrigger("HLT_DoubleMu7_v4",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);//auto
  //mymod->AddTrigger("HLT_DoubleMu7_v5",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);//auto
  //mymod->AddTrigger("HLT_DoubleMu7_v7",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);//auto
  //mymod->AddTrigger("HLT_DoubleMu7_v8",kHLT_DoubleMu7,"hltDiMuonL3PreFiltered7",kHLT_DoubleMu7_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu13_Mu8_v1",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu13_Mu8_v2",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  //mymod->AddTrigger("HLT_Mu13_Mu8_v3",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu13_Mu8_v4",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);  
  mymod->AddTrigger("HLT_Mu13_Mu8_v6",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu13_Mu8_v7",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu13_Mu8_v16",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu13_Mu8_v17",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu13_Mu8_v18",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu13_Mu8_v19",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu13_Mu8_v21",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu13_Mu8_v10",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu13_Mu8_v11",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu13_Mu8_v22",kHLT_Mu13_Mu8,"hltSingleMu13L3Filtered13",kHLT_Mu13_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu13_Mu8_Mu2Obj,0);//auto
  //mymod->AddTrigger("HLT_Mu17_Mu8_v1",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  //mymod->AddTrigger("HLT_Mu17_Mu8_v2",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  //mymod->AddTrigger("HLT_Mu17_Mu8_v3",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  //mymod->AddTrigger("HLT_Mu17_Mu8_v4",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  //mymod->AddTrigger("HLT_Mu17_Mu8_v6",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v10",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v11",kHLT_Mu17_Mu8,"hltSingleMu13L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltDiMuonL3p5PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v13",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v14",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v15",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v16",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);
  mymod->AddTrigger("HLT_Mu17_Mu8_v17",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Mu8_v18",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Mu8_v19",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Mu8_v21",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Mu8_v7",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);//auto
  mymod->AddTrigger("HLT_Mu17_Mu8_v22",kHLT_Mu17_Mu8,"hltL3fL1DoubleMu10MuOpenL1f0L2f10L3Filtered17",kHLT_Mu17_Mu8_Mu1Obj,0,"hltL3pfL1DoubleMu10MuOpenL1f0L2pf0L3PreFiltered8",kHLT_Mu17_Mu8_Mu2Obj,0);//auto
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
  mymod->AddTrigger("HLT_Mu12_v16",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);//auto
  mymod->AddTrigger("HLT_Mu12_v17",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);//auto
  mymod->AddTrigger("HLT_Mu12_v18",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);//auto
  mymod->AddTrigger("HLT_Mu12_v11",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);//auto
  mymod->AddTrigger("HLT_Mu12_v12",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);//auto
  mymod->AddTrigger("HLT_Mu12_v7",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);//auto
  mymod->AddTrigger("HLT_Mu12_v8",kHLT_Mu12,"hltSingleMu12L3Filtered12",kHLT_Mu12_MuObj);//auto
  mymod->AddTrigger("HLT_Mu15_v1",kHLT_Mu15,"",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v2",kHLT_Mu15,"hltL3Muon15",kHLT_Mu15_MuObj);
  //mymod->AddTrigger("HLT_Mu15_v3",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  //mymod->AddTrigger("HLT_Mu15_v4",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  //mymod->AddTrigger("HLT_Mu15_v5",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  //mymod->AddTrigger("HLT_Mu15_v6",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu24_v1",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  mymod->AddTrigger("HLT_Mu24_v2",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  //mymod->AddTrigger("HLT_Mu24_v3",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  //mymod->AddTrigger("HLT_Mu24_v4",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  //mymod->AddTrigger("HLT_Mu24_v5",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);
  //mymod->AddTrigger("HLT_Mu24_v14",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu24_v15",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu24_v16",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu24_v11",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu24_v12",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu24_v7",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu24_v8",kHLT_Mu24,"hltSingleMu24L3Filtered24",kHLT_Mu24_MuObj);//auto
  mymod->AddTrigger("HLT_Mu30_v1",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v2",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v3",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v4",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v5",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);
  mymod->AddTrigger("HLT_Mu30_v14",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);//auto
  mymod->AddTrigger("HLT_Mu30_v15",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);//auto
  mymod->AddTrigger("HLT_Mu30_v16",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);//auto
  mymod->AddTrigger("HLT_Mu30_v11",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);//auto
  mymod->AddTrigger("HLT_Mu30_v12",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);//auto
  mymod->AddTrigger("HLT_Mu30_v7",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);//auto
  mymod->AddTrigger("HLT_Mu30_v8",kHLT_Mu30,"hltSingleMu30L3Filtered30",kHLT_Mu30_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_v5",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);
  //mymod->AddTrigger("HLT_Mu40_v12",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_v13",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_v14",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_v10",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_v2",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_v3",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_v6",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_v9",kHLT_Mu40,"hltSingleMu40L2QualL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_eta2p1_v1",kHLT_Mu40,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40",kHLT_Mu40_MuObj);
  //mymod->AddTrigger("HLT_Mu40_eta2p1_v4",kHLT_Mu40,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40",kHLT_Mu40_MuObj);
  //mymod->AddTrigger("HLT_Mu40_eta2p1_v10",kHLT_Mu40,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_eta2p1_v11",kHLT_Mu40,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_eta2p1_v9",kHLT_Mu40,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40",kHLT_Mu40_MuObj);//auto
  //mymod->AddTrigger("HLT_Mu40_eta2p1_v5",kHLT_Mu40,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered40",kHLT_Mu40_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu17_v5",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);
  mymod->AddTrigger("HLT_IsoMu17_v6",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);
  mymod->AddTrigger("HLT_IsoMu17_v8",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);
  mymod->AddTrigger("HLT_IsoMu17_v9",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);
  mymod->AddTrigger("HLT_IsoMu17_v10",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu17_v11",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu17_v13",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu17_v14",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17",kHLT_IsoMu17_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu17_eta2p1_v1",kHLT_IsoMu17,"hltSingleMuIsoL3IsoFiltered17Eta21",kHLT_IsoMu17_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v1",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v2",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v4",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v5",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v6",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v7",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v8",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_v15",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_v16",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_v17",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_v12",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_v13",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_v9",kHLT_IsoMu24,"hltSingleMuIsoL3IsoFiltered24",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v9",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v10",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v11",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v12",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v13",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v14",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v15",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v3",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v6",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v7",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoFiltered10",kHLT_IsoMu24_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v13",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v14",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15",kHLT_IsoMu24_MuObj);
  mymod->AddTrigger("HLT_IsoMu24_eta2p1_v15",kHLT_IsoMu24,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f24QL3crIsoRhoFiltered0p15",kHLT_IsoMu24_MuObj);
  //mymod->AddTrigger("HLT_IsoMu30_eta2p1_v11",kHLT_IsoMu30,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f30L3IsoFiltered",kHLT_IsoMu30_MuObj);//auto
  //mymod->AddTrigger("HLT_IsoMu30_eta2p1_v12",kHLT_IsoMu30,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f30L3IsoFiltered",kHLT_IsoMu30_MuObj);//auto
  //mymod->AddTrigger("HLT_IsoMu30_eta2p1_v13",kHLT_IsoMu30,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f30L3IsoFiltered",kHLT_IsoMu30_MuObj);//auto
  //mymod->AddTrigger("HLT_IsoMu30_eta2p1_v14",kHLT_IsoMu30,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f30L3IsoFiltered",kHLT_IsoMu30_MuObj);//auto
  //mymod->AddTrigger("HLT_IsoMu30_eta2p1_v15",kHLT_IsoMu30,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f30L3IsoFiltered",kHLT_IsoMu30_MuObj);//auto
  //mymod->AddTrigger("HLT_IsoMu30_eta2p1_v3",kHLT_IsoMu30,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f30L3IsoFiltered",kHLT_IsoMu30_MuObj);
  //mymod->AddTrigger("HLT_IsoMu30_eta2p1_v3",kHLT_IsoMu30,"hltL3IsoL1sMu14Eta2p1L1f0L2f14QL2IsoL3f30L3IsoFiltered",kHLT_IsoMu30_MuObj);

  //
  // DoubleElectron
  //
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v1",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v2",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v3",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v4",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v5",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v6",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v7",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_v8",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj,0,"hltEle17CaloIdIsoEle8CaloIdIsoPixelMatchDoubleFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj,0);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v2",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v5",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v9",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsolDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v16",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v17",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v18",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v19",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj,0,"hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDoubleFilter",kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj,0);//auto

  //
  // SingleElectron
  //
  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle27CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v1",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v2",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v3",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v4",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v5",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v6",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_v7",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT,"hltEle32CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v3",kHLT_Ele52_CaloIdVT_TrkIdT,"hltEle52CaloIdVTTrkIdTDphiFilter",kHLT_Ele52_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v4",kHLT_Ele52_CaloIdVT_TrkIdT,"hltEle52CaloIdVTTrkIdTDphiFilter",kHLT_Ele52_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v1",kHLT_Ele52_CaloIdVT_TrkIdT,"hltEle52CaloIdVTTrkIdTDphiFilter",kHLT_Ele52_CaloIdVT_TrkIdT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele52_CaloIdVT_TrkIdT_v2",kHLT_Ele52_CaloIdVT_TrkIdT,"hltEle52CaloIdVTTrkIdTDphiFilter",kHLT_Ele52_CaloIdVT_TrkIdT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v1",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v2",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v3",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v4",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v5",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v11",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v12",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v13",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele65_CaloIdVT_TrkIdT_v6",kHLT_Ele65_CaloIdVT_TrkIdT,"hltEle65CaloIdVTTrkIdTDphiFilter",kHLT_Ele65_CaloIdVT_TrkIdT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v2",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v5",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v6",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v7",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v10",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v3",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);//auto
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v8",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v9",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele80_CaloIdVT_TrkIdT_v8",kHLT_Ele80_CaloIdVT_TrkIdT,"hltEle80CaloIdVTTrkIdTDphiFilter",kHLT_Ele80_CaloIdVT_TrkIdT_EleObj);
  mymod->AddTrigger("HLT_Ele27_WP80_v5",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);
  mymod->AddTrigger("HLT_Ele27_WP80_v6",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);
  mymod->AddTrigger("HLT_Ele27_WP80_v7",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);
  mymod->AddTrigger("HLT_Ele27_WP80_v8",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);
  mymod->AddTrigger("HLT_Ele27_WP80_v9",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);
  mymod->AddTrigger("HLT_Ele27_WP80_v10",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);//auto
  mymod->AddTrigger("HLT_Ele27_WP80_v11",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);//auto
  mymod->AddTrigger("HLT_Ele27_WP80_v2",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);//auto
  mymod->AddTrigger("HLT_Ele27_WP80_v3",kHLT_Ele27_WP80,"hltEle27WP80TrackIsoFilter",kHLT_Ele27_WP80_EleObj);//auto

  //T&P triggers

  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v1",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, "hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter", kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj, 0, "hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter", kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v2",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, "hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter", kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj, 0, "hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter", kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdL_CaloIsoVL_SC17_v3",                   kHLT_Ele32_CaloIdL_CaloIsoVL_SC17, "hltEle32CaloIdLCaloIsoVLSC17PixelMatchFilter", kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj, 0, "hltEle32CaloIdLCaloIsoVLSC17HEDoubleFilter", kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v1",     kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v2",     kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v3",     kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v4",     kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v5",     kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v6",     kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v7",     kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_v8",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsolFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj,0,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17HEDoubleFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj,0);//auto

  //mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_v1",    kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele1Obj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17PixelMatchDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele2Obj, 0);
  //mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_v2",    kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele1Obj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17PixelMatchDoubleFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele2Obj, 0);
  //mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_v3",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17TrackIsolFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele1Obj,0,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTEle17PixelMatchDoubleFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele2Obj,0);//auto

  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v1", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v2", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v3", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v4", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v5", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v6", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v7", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v8", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v9", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_v10",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v1", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj, 0, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj);
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v2", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj, 0, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj);
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v3", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj, 0, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj);
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v4", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj, 0, "hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter", kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj);
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v5",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50,"hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj,0,"hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v6",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50,"hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj,0,"hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_v7",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50,"hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4TrackIsoFilter",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj,0,"hltEle20CaloIdVTCaloIsoVTTrkIdTTrkIsoVTSC4PMMassFilter",kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj);//auto
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v1", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_SCObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v2", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_SCObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v3", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_SCObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v4", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_EleObj, 0, "hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter", kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_SCObj);
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v5",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_EleObj,0,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_SCObj);//auto
  mymod->AddTrigger("HLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_v6",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_EleObj,0,"hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter",kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_SCObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v2", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v3", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v4", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v5", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v6", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v7", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v8", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_v9",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsolFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj,0);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v1", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v2", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v3", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v4", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj, 0, "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj, 0);
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v5",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj,0);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v6",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj,0,"hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter",kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj,0);//auto

  //Fake Rate triggers

  mymod->AddTrigger("HLT_Mu8_v1",                                kHLT_Mu8, "hltSingleMu8L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v2",                                kHLT_Mu8, "hltSingleMu8L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v3",                                kHLT_Mu8, "hltSingleMu8L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v5",                                kHLT_Mu8, "hltSingleMu8L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v7",                                kHLT_Mu8, "hltSingleMu8L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v8",                                kHLT_Mu8, "hltSingleMu8L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v11",                                kHLT_Mu8, "hltSingleMu8L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v14",                                kHLT_Mu8, "hltL3fL1sMu3L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v15",                                kHLT_Mu8, "hltL3fL1sMu3L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v16",                                kHLT_Mu8, "hltL3fL1sMu3L3Filtered8", kHLT_Mu8_MuObj);
  mymod->AddTrigger("HLT_Mu8_v17",kHLT_Mu8,"hltL3fL1sMu3L3Filtered8",kHLT_Mu8_MuObj);//auto
  mymod->AddTrigger("HLT_Mu8_v18",kHLT_Mu8,"hltL3fL1sMu3L3Filtered8",kHLT_Mu8_MuObj);//auto
  mymod->AddTrigger("HLT_Mu8_v12",kHLT_Mu8,"hltL3fL1sMu3L3Filtered8",kHLT_Mu8_MuObj);//auto
  mymod->AddTrigger("HLT_Mu15_v1",                               kHLT_Mu15, "hltL3Muon15", kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v2",                               kHLT_Mu15, "hltL3Muon15", kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v3",                               kHLT_Mu15, "hltSingleMu15L3Filtered15", kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v4",                               kHLT_Mu15, "hltSingleMu15L3Filtered15", kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v6",                               kHLT_Mu15, "hltSingleMu15L3Filtered15", kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v8",                               kHLT_Mu15, "hltSingleMu15L3Filtered15", kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v9",                               kHLT_Mu15, "hltSingleMu15L3Filtered15", kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v12",                               kHLT_Mu15, "hltSingleMu15L3Filtered15", kHLT_Mu15_MuObj);
  mymod->AddTrigger("HLT_Mu15_v13",kHLT_Mu15,"hltSingleMu15L3Filtered15",kHLT_Mu15_MuObj);//auto
  mymod->AddTrigger("HLT_Mu17_v1",                               kHLT_Mu17, "hltL3fL1sMu12L3Filtered17", kHLT_Mu17_MuObj);
  mymod->AddTrigger("HLT_Mu17_v2",                               kHLT_Mu17, "hltL3fL1sMu12L3Filtered17", kHLT_Mu17_MuObj);
  mymod->AddTrigger("HLT_Mu17_v3",                               kHLT_Mu17, "hltL3fL1sMu12L3Filtered17", kHLT_Mu17_MuObj);
  mymod->AddTrigger("HLT_Mu17_v4",kHLT_Mu17,"hltL3fL1sMu12L3Filtered17",kHLT_Mu17_MuObj);//auto
  mymod->AddTrigger("HLT_Mu17_v5",kHLT_Mu17,"hltL3fL1sMu12L3Filtered17",kHLT_Mu17_MuObj);//auto
  mymod->AddTrigger("HLT_Mu8_Jet40_v1",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8" , kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v2",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v3",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v4",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v5",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v6",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v7",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v9",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v10",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v14",                          kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8", kHLT_Mu8_Jet40_MuObj);
  mymod->AddTrigger("HLT_Mu8_Jet40_v15",kHLT_Mu8_Jet40,"hltL3Mu8Jet20L3Filtered8",kHLT_Mu8_Jet40_MuObj);//auto
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v2",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMu8EG5L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v3",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMu8EG5L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v3",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMu8EG5L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v4",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMu8EG5L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v5",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMu8EG5L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v6",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMu8EG5L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v8",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMuOpenEG12L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v9",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMuOpenEG12L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v12",         kHLT_Mu8_Photon20_CaloIdVT_IsoT ,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMuOpenEG12L3Filtered8", kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);
  mymod->AddTrigger("HLT_Mu8_Photon20_CaloIdVT_IsoT_v13",kHLT_Mu8_Photon20_CaloIdVT_IsoT,"hltPhoton20CaloIdVTIsoTMu8TrackIsoFilter",kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj,0,"hltSingleMuOpenEG12L3Filtered8",kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj,0);//auto

  mymod->AddTrigger("HLT_Ele8_v1",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v2",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v3",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v4",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v5",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v6",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v7",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v8",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v9",                               kHLT_Ele8,"hltEle8PixelMatchFilter", kHLT_Ele8_EleObj);
  mymod->AddTrigger("HLT_Ele8_v10",kHLT_Ele8,"hltEle8PixelMatchFilter",kHLT_Ele8_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v1",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v2",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v3",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v4",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v5",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v6",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v7",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v8",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v9",             kHLT_Ele8_CaloIdL_CaloIsoVL, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v14",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v15",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v16",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v17",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_v10",kHLT_Ele8_CaloIdL_CaloIsoVL,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v3", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsolFilter", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_v4", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsolFilter", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsolFilter", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsolFilter", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsolFilter", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v7", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsolFilter", kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v10",            kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsoFilter", kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v11",            kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsoFilter", kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v12",            kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsoFilter", kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v13",            kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "hltEle8TightIdLooseIsoTrackIsoFilter", kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v14",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v15",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v8",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v1",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v2",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v3",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v4",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v5",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v6",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v7",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v8",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v9",            kHLT_Ele17_CaloIdL_CaloIsoVL, "hltEle17CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v14",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v15",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v16",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v17",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdL_CaloIsoVL_v10",kHLT_Ele17_CaloIdL_CaloIsoVL,"hltEle17CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v1",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v2",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v3",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v4",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v5",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v6",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v7",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v8",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v11",       kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40, "hltEle8CaloIdLCaloIsoVLPixelMatchFilter", kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdL_CaloIsoVL_Jet40_v12",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40,"hltEle8CaloIdLCaloIsoVLPixelMatchFilter",kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v1",       kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, "hltEle8TightIdLooseIsoTrackIsoFilter", kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v2",       kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, "hltEle8TightIdLooseIsoTrackIsoFilter", kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v3",       kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, "hltEle8TightIdLooseIsoTrackIsoFilter", kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v4",       kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, "hltEle8TightIdLooseIsoTrackIsoFilter", kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v5",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30,"hltEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v6",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30,"hltEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);//auto
  mymod->AddTrigger("HLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v7",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30,"hltEle8TightIdLooseIsoTrackIsoFilter",kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v1",       kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter", kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v2",       kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter", kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v3",       kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter", kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v4",       kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter", kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v5",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30,"hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v6",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30,"hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_v7",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30,"hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v1",       kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter", kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v2",       kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter", kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v3",       kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter", kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v4",       kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL, "hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter", kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v5",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_v6",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL,"hltEle17CaloIdTCaloIsoVLTrkIdVLTrkIsoVLTrackIsoFilter",kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj);//auto
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v1", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v2", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v3", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v4", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v5", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v6", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v8", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v9", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v10", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL , "hltPhoton20CaloIdVTIsoTTrackIsoFilter", kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj, 0, "hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter" , kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj, 0);
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v11",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj,0,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj,0);//auto
  mymod->AddTrigger("HLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_v7",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL,"hltPhoton20CaloIdVTIsoTTrackIsoFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj,0,"hltEle8CaloIdLCaloIsoVLNoL1SeedPixelMatchFilter",kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj,0);//auto


  // Double Hadronic 2012
//   mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v1", kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30, "hltDoublePFTau30TrackPt1MediumIsolation",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj);
//   mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v4",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30,"hltDoublePFTau30TrackPt1MediumIsolation",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj);//auto
//   mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v5",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30,"hltDoublePFTau30TrackPt1MediumIsolation",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj);//auto
 
//   mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v2",kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30,"hltDoublePFTau30TrackPt5MediumIsolation",kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj);//auto
//   mymod->AddTrigger("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v2",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30,"hltDoublePFTau25TrackPt5MediumIsolation",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj);//auto
//   mymod->AddTrigger("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v3",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30,"hltDoublePFTau25TrackPt5MediumIsolation",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj);//auto
//   mymod->AddTrigger("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v4",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30,"hltDoublePFTau25TrackPt5MediumIsolation",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj);//auto

mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v1", kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30, "hltDoublePFTau30TrackPt1MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj);
mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v4",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30,"hltDoublePFTau30TrackPt1MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj);//auto
mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30_v5",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30,"hltDoublePFTau30TrackPt1MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj);//auto
  
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v1",kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30,"hltDoublePFTau30TrackPt5MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj);//auto 
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30_v2",kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30,"hltDoublePFTau30TrackPt5MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj);//auto
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v1",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30,"hltDoublePFTau25TrackPt5MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj);
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v2",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30,"hltDoublePFTau25TrackPt5MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj);//auto
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v3",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30,"hltDoublePFTau25TrackPt5MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj);//auto
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30_v4",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30,"hltDoublePFTau25TrackPt5MediumIsolationProng4Dz02",kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj);//auto

  //1prong
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1_v4",kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1,"hltDoublePFTau35TrackPt1MediumIsolationProng2Dz02",kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1Obj);//auto 
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1_v1",kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1,"hltDoublePFTau35TrackPt1MediumIsolationProng2Dz02",kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1Obj);//auto
  mymod->AddTrigger("HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1_v3",kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1,"hltDoublePFTau35TrackPt1MediumIsolationProng2Dz02",kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Prong1Obj);//auto

  //Double Hadronic 2011
 
   mymod->AddTrigger("HLT_DoubleIsoPFTau35_Trk5_eta2p1_v4", kHLT_DoubleIsoPFTau35_Trk5_eta2p1, "hltFilterDoubleIsoPFTau35Trk5LeadTrack5IsolationL1HLTMatched",kHLT_DoubleIsoPFTau35_Trk5_eta2p1Obj);
   mymod->AddTrigger("HLT_DoubleIsoPFTau35_Trk5_eta2p1_v2", kHLT_DoubleIsoPFTau35_Trk5_eta2p1, "hltFilterDoubleIsoPFTau35Trk5LeadTrack5IsolationL1HLTMatched",kHLT_DoubleIsoPFTau35_Trk5_eta2p1Obj);
 mymod->AddTrigger("HLT_DoubleIsoPFTau35_Trk5_eta2p1_v3", kHLT_DoubleIsoPFTau35_Trk5_eta2p1, "hltFilterDoubleIsoPFTau35Trk5LeadTrack5IsolationL1HLTMatched",kHLT_DoubleIsoPFTau35_Trk5_eta2p1Obj);
   mymod->AddTrigger("HLT_DoubleIsoPFTau25_Trk5_eta2p1_v2", kHLT_DoubleIsoPFTau25_Trk5_eta2p1, "hltFilterDoubleIsoPFTau25Trk5LeadTrack5IsolationL1HLTMatched",kHLT_DoubleIsoPFTau25_Trk5_eta2p1Obj);
   mymod->AddTrigger(" HLT_DoubleIsoPFTau20_Trk5_v4", kHLT_DoubleIsoPFTau20_Trk5, "hltFilterDoubleIsoPFTau20Trk5LeadTrack5IsolationL1HLTMatched",kHLT_DoubleIsoPFTau20_Trk5Obj);
   mymod->AddTrigger(" HLT_DoubleIsoPFTau20_Trk5_v2", kHLT_DoubleIsoPFTau20_Trk5, "hltFilterDoubleIsoPFTau20Trk5LeadTrack5IsolationL1HLTMatched",kHLT_DoubleIsoPFTau20_Trk5Obj);
   mymod->AddTrigger(" HLT_DoubleIsoPFTau20_Trk5_v1", kHLT_DoubleIsoPFTau20_Trk5, "hltFilterDoubleIsoPFTau20Trk5LeadTrack5IsolationL1HLTMatched",kHLT_DoubleIsoPFTau20_Trk5Obj);
  
  // Mu + MET
  mymod->AddTrigger("HLT_IsoMu15_L1ETM20_v0", kHLT_IsoMu15_L1ETM20, "", kHLT_IsoMu15_L1ETM20_MuObj);
  mymod->AddTrigger("HLT_IsoMu15_L1ETM20_v3",kHLT_IsoMu15_L1ETM20,"",kHLT_IsoMu15_L1ETM20_MuObj);//auto
  mymod->AddTrigger("HLT_IsoMu15_L1ETM20_v4",kHLT_IsoMu15_L1ETM20,"",kHLT_IsoMu15_L1ETM20_MuObj);//auto

  // Tau + Muon triggers
  mymod->AddTrigger("HLT_IsoMu12_LooseIsoPFTau10_v0",kHLT_IsoMu12_LooseIsoPFTau10,"hltSingleMuIsoL3IsoFiltered12",kHLT_IsoMu12_LooseIsoPFTau10_MuObj, 0,"hltFilterIsoMu12IsoPFTau10LooseIsolation",kHLT_IsoMu12_LooseIsoPFTau10_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu12_LooseIsoPFTau10_v1",kHLT_IsoMu12_LooseIsoPFTau10,"hltSingleMuIsoL3IsoFiltered12",kHLT_IsoMu12_LooseIsoPFTau10_MuObj,0,"hltFilterIsoMu12IsoPFTau10LooseIsolation",kHLT_IsoMu12_LooseIsoPFTau10_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu12_LooseIsoPFTau10_v2",kHLT_IsoMu12_LooseIsoPFTau10,"hltSingleMuIsoL3IsoFiltered12",kHLT_IsoMu12_LooseIsoPFTau10_MuObj,0,"hltFilterIsoMu12IsoPFTau10LooseIsolation",kHLT_IsoMu12_LooseIsoPFTau10_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu12_LooseIsoPFTau10_v4",kHLT_IsoMu12_LooseIsoPFTau10,"hltSingleMuIsoL3IsoFiltered12",kHLT_IsoMu12_LooseIsoPFTau10_MuObj,0,"hltFilterIsoMu12IsoPFTau10LooseIsolation",kHLT_IsoMu12_LooseIsoPFTau10_TauObj,0);//auto

  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v0",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj, 0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj, 0);
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v13",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v14",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v2",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v4",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v5",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v6",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v8",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau15_v9",kHLT_Mu15_LooseIsoPFTau15,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_Mu15_LooseIsoPFTau15_TauObj,0);//auto

  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau20_v0",kHLT_Mu15_LooseIsoPFTau20,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau20_MuObj, 0,"hltPFTau20TrackLooseIso",kHLT_Mu15_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau20_v1",kHLT_Mu15_LooseIsoPFTau20,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_Mu15_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau20_v2",kHLT_Mu15_LooseIsoPFTau20,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_Mu15_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu15_LooseIsoPFTau20_v4",kHLT_Mu15_LooseIsoPFTau20,"hltSingleMu15L3Filtered15",kHLT_Mu15_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_Mu15_LooseIsoPFTau20_TauObj,0);//auto

  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau15_v0",kHLT_IsoMu15_LooseIsoPFTau15,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau15_MuObj, 0,"hltPFTau15TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau15_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau15_v2",kHLT_IsoMu15_LooseIsoPFTau15,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau15_v4",kHLT_IsoMu15_LooseIsoPFTau15,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau15_v5",kHLT_IsoMu15_LooseIsoPFTau15,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau15_v6",kHLT_IsoMu15_LooseIsoPFTau15,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau15_v8",kHLT_IsoMu15_LooseIsoPFTau15,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau15_v9",kHLT_IsoMu15_LooseIsoPFTau15,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau15_MuObj,0,"hltPFTau15TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau15_TauObj,0);//auto

  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau20_v0",kHLT_IsoMu15_LooseIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau20_MuObj, 0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau20_v2",kHLT_IsoMu15_LooseIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau20_v3",kHLT_IsoMu15_LooseIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau20_v4",kHLT_IsoMu15_LooseIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_LooseIsoPFTau20_v6",kHLT_IsoMu15_LooseIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_LooseIsoPFTau20_TauObj,0);//auto

  mymod->AddTrigger("HLT_IsoMu15_TightIsoPFTau20_v0",kHLT_IsoMu15_TightIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_TightIsoPFTau20_MuObj, 0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_TightIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu15_TightIsoPFTau20_v2",kHLT_IsoMu15_TightIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_TightIsoPFTau20_MuObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_TightIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_TightIsoPFTau20_v3",kHLT_IsoMu15_TightIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_TightIsoPFTau20_MuObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_TightIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_TightIsoPFTau20_v4",kHLT_IsoMu15_TightIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_TightIsoPFTau20_MuObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_TightIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_TightIsoPFTau20_v6",kHLT_IsoMu15_TightIsoPFTau20,"hltSingleMuIsoL3IsoFiltered15",kHLT_IsoMu15_TightIsoPFTau20_MuObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_TightIsoPFTau20_TauObj,0);//auto

  mymod->AddTrigger("HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v0",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_MuObj, 0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v1",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v5",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_LooseIsoPFTau20_v6",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20TrackLooseIso",kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_TauObj,0);//auto

  mymod->AddTrigger("HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v0",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_MuObj, 0,"hltPFTauMediumIso20TrackMediumIso",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v1",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_MuObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v5",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_MuObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_MediumIsoPFTau20_v6",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_MuObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj,0);//auto

  mymod->AddTrigger("HLT_IsoMu15_eta2p1_TightIsoPFTau20_v0",kHLT_IsoMu15_eta2p1_TightIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_TightIsoPFTau20_MuObj, 0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_TightIsoPFTau20_v1",kHLT_IsoMu15_eta2p1_TightIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_TightIsoPFTau20_MuObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_TightIsoPFTau20_v5",kHLT_IsoMu15_eta2p1_TightIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_TightIsoPFTau20_MuObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu15_eta2p1_TightIsoPFTau20_v6",kHLT_IsoMu15_eta2p1_TightIsoPFTau20,"hltSingleMuIsoL1s14L3IsoFiltered15eta2p1",kHLT_IsoMu15_eta2p1_TightIsoPFTau20_MuObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj,0);//auto

  mymod->AddTrigger("HLT_IsoMu18_eta2p1_LooseIsoPFTau20",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoFiltered10",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_MuObj, 0,"hltPFTau20IsoMuVertex",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v4",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoFiltered10",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20IsoMuVertex",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v5",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoFiltered10",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20IsoMuVertex",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu18_eta2p1_LooseIsoPFTau20_v6",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu16Eta2p1L1f0L2f16QL3f18QL3crIsoFiltered10",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20IsoMuVertex",kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_TauObj,0);//auto

  mymod->AddTrigger("HLT_Mu18_eta2p1_LooseIsoPFTau20_v0",kHLT_Mu18_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered18Q",kHLT_Mu18_eta2p1_LooseIsoPFTau20_MuObj, 0,"hltPFTau20MuVertex",kHLT_Mu18_eta2p1_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Mu18_eta2p1_LooseIsoPFTau20_v4",kHLT_Mu18_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered18Q",kHLT_Mu18_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20MuVertex",kHLT_Mu18_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu18_eta2p1_LooseIsoPFTau20_v5",kHLT_Mu18_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered18Q",kHLT_Mu18_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20MuVertex",kHLT_Mu18_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu18_eta2p1_LooseIsoPFTau20_v6",kHLT_Mu18_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu16Eta2p1L1f0L2f16QL3Filtered18Q",kHLT_Mu18_eta2p1_LooseIsoPFTau20_MuObj,0,"hltPFTau20MuVertex",kHLT_Mu18_eta2p1_LooseIsoPFTau20_TauObj,0);//auto

  mymod->AddTrigger("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v0",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu14erORMu16erL1f0L2f14QL3f17QL3crIsoRhoFiltered0p15",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_MuObj, 0,"hltIsoMuPFTau20TrackLooseIso",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v2",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu14erORMu16erL1f0L2f14QL3f17QL3crIsoRhoFiltered0p15",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_MuObj,0,"hltIsoMuPFTau20TrackLooseIso",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v3",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu14erORMu16erL1f0L2f14QL3f17QL3crIsoRhoFiltered0p15",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_MuObj,0,"hltIsoMuPFTau20TrackLooseIso",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v6",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu14erORMu16erL1f0L2f14QL3f17QL3crIsoRhoFiltered0p15",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_MuObj,0,"hltIsoMuPFTau20TrackLooseIso",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v7",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20,"hltL3crIsoL1sMu14erORMu16erL1f0L2f14QL3f17QL3crIsoRhoFiltered0p15",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_MuObj,0,"hltIsoMuPFTau20TrackLooseIso",kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  
  mymod->AddTrigger("HLT_Mu17_eta2p1_LooseIsoPFTau20_v0",kHLT_Mu17_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu14erORMu16erL1f0L2f14QL3Filtered17Q",kHLT_Mu17_eta2p1_LooseIsoPFTau20_MuObj, 0,"hltMuPFTau20TrackLooseIso",kHLT_Mu17_eta2p1_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Mu17_eta2p1_LooseIsoPFTau20_v2",kHLT_Mu17_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu14erORMu16erL1f0L2f14QL3Filtered17Q",kHLT_Mu17_eta2p1_LooseIsoPFTau20_MuObj,0,"hltMuPFTau20TrackLooseIso",kHLT_Mu17_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu17_eta2p1_LooseIsoPFTau20_v3",kHLT_Mu17_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu14erORMu16erL1f0L2f14QL3Filtered17Q",kHLT_Mu17_eta2p1_LooseIsoPFTau20_MuObj,0,"hltMuPFTau20TrackLooseIso",kHLT_Mu17_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu17_eta2p1_LooseIsoPFTau20_v6",kHLT_Mu17_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu14erORMu16erL1f0L2f14QL3Filtered17Q",kHLT_Mu17_eta2p1_LooseIsoPFTau20_MuObj,0,"hltMuPFTau20TrackLooseIso",kHLT_Mu17_eta2p1_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Mu17_eta2p1_LooseIsoPFTau20_v7",kHLT_Mu17_eta2p1_LooseIsoPFTau20,"hltL3fL1sMu14erORMu16erL1f0L2f14QL3Filtered17Q",kHLT_Mu17_eta2p1_LooseIsoPFTau20_MuObj,0,"hltMuPFTau20TrackLooseIso",kHLT_Mu17_eta2p1_LooseIsoPFTau20_TauObj,0);//auto


  // Tau + Electron triggers
  // Non-Iso Ele15 + LooseTau15
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v0",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_EleObj, 0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_TauObj, 0);
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v1",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_EleObj,0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v2",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_EleObj,0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v4",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_EleObj,0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_v6",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_EleObj,0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_TauObj,0);//auto

  // Non-Iso Ele15 + LooseTau20
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_v0",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj, 0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_v2",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_v3",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj,0);//auto

  // Non-Iso Ele15 + TightTau20
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_v0",kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_EleObj, 0,"hltPFTauTightIso20TrackTightIso",kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_v2",kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20,"hltEle15CaloIdVTTrkIdTDphiFilter",kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_EleObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_TauObj,0);//auto

  // Non-Iso Ele18 + MediumTau20
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v0",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20,"hltEle18CaloIdVTTrkIdTDphiFilter",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_EleObj, 0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v1",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20,"hltEle18CaloIdVTTrkIdTDphiFilter",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_EleObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v5",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20,"hltEle18CaloIdVTTrkIdTDphiFilter",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_EleObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_v6",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20,"hltEle18CaloIdVTTrkIdTDphiFilter",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_EleObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_TauObj,0);//auto

  // Ele 15 + LoosePFTau15
 
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v0",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15,"hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj, 0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj, 0);

  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v0",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj, 0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj, 0);
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v1",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj,0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v2",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj,0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v4",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj,0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_v6",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj,0,"hltPFTau15TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj,0);//auto

  // Ele 15 + LoosePFTau20
 
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v0",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle15CaloIdVTTrkIdTCaloIsoTTrkIsoTTrackIsolFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj, 0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj, 0);

  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v0",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj, 0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v1",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v4",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v6",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v8",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v9",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto

  // Ele 15 + TightPFTau20
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v0",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_EleObj,  0,"hltPFTauTightIso20TrackTightIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_v2",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20,"hltEle15CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_EleObj,0,"hltPFTauTightIso20TrackTightIso",kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_TauObj,0);//auto

  // Ele 18 + LoosePFTau20
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v0",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,  0,"hltPFTau20TrackLooseIso",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v2",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v3",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTau20TrackLooseIso",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto

  // Ele18_MediumIsoPFTau20
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v0",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20,"hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj,  0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20,"hltEle18CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilter",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj,0);//auto

  
  // Ele20_MediumIsoPFTau20
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v0",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20,"hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj,  0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v1",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20,"hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v5",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20,"hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_v6",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20,"hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1SingleEG18orL1SingleEG20",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj,0,"hltPFTauMediumIso20TrackMediumIso",kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj,0);//auto

  // Ele20 CaloIsoRho + LooseIsoPFTau20
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v0",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1IsoEG18OrEG20",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,  0,"hltPFTauIsoEleVertex20",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v4",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1IsoEG18OrEG20",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTauIsoEleVertex20",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v5",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1IsoEG18OrEG20",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTauIsoEleVertex20",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_v6",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20,"hltEle20CaloIdVTCaloIsoTTrkIdTTrkIsoTTrackIsoFilterL1IsoEG18OrEG20",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj,0,"hltPFTauIsoEleVertex20",kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj,0);//auto


  // Ele22 WP90 + LooseIsoPFTau20
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v0", kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20, "hltEle22WP90RhoTrackIsoFilter", kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_EleObj,  0,"hltIsoElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v2",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20,"hltEle22WP90RhoTrackIsoFilter",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_EleObj,0,"hltIsoElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v3",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20,"hltEle22WP90RhoTrackIsoFilter",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_EleObj,0,"hltIsoElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v6",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20,"hltEle22WP90RhoTrackIsoFilter",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_EleObj,0,"hltIsoElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_v7",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20,"hltEle22WP90RhoTrackIsoFilter",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_EleObj,0,"hltIsoElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_TauObj,0);//auto
                   
  // Non-Iso Ele22 WP90 + LooseIsoPFTau20
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_v0",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20,"hltEle22WP90NoIsoDphiFilter",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_EleObj, 0,"hltElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_v2",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20,"hltEle22WP90NoIsoDphiFilter",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_EleObj,0,"hltElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_v3",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20,"hltEle22WP90NoIsoDphiFilter",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_EleObj,0,"hltElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_v6",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20,"hltEle22WP90NoIsoDphiFilter",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_EleObj,0,"hltElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_v7",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20,"hltEle22WP90NoIsoDphiFilter",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_EleObj,0,"hltElePFTau20TrackLooseIso",kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_TauObj,0);//auto


  // Ele20 + LooseIsoPFTau20
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_v0",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20,"hltEle20CaloIdVTTrkIdTDphiFilter",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj,  0,"hltPFTauEleVertex20",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj, 0);
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_v4",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20,"hltEle20CaloIdVTTrkIdTDphiFilter",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj,0,"hltPFTauEleVertex20",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_v5",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20,"hltEle20CaloIdVTTrkIdTDphiFilter",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj,0,"hltPFTauEleVertex20",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj,0);//auto
  mymod->AddTrigger("HLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_v6",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20,"hltEle20CaloIdVTTrkIdTDphiFilter",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj,0,"hltPFTauEleVertex20",kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj,0);//auto

  // Jet triggers (lazy, didn't lookup L3 module names)
  mymod->AddTrigger("HLT_Jet150_v0",    kHLT_Jet150, "", kHLT_Jet150_JetObj);
  mymod->AddTrigger("HLT_Jet150_v1",kHLT_Jet150,"",kHLT_Jet150_JetObj);//auto
  mymod->AddTrigger("HLT_Jet150_v2",kHLT_Jet150,"",kHLT_Jet150_JetObj);//auto
  mymod->AddTrigger("HLT_Jet150_v3",kHLT_Jet150,"",kHLT_Jet150_JetObj);//auto
  mymod->AddTrigger("HLT_Jet150_v4",kHLT_Jet150,"",kHLT_Jet150_JetObj);//auto
  mymod->AddTrigger("HLT_Jet150_v5",kHLT_Jet150,"",kHLT_Jet150_JetObj);//auto
  mymod->AddTrigger("HLT_Jet150_v6",kHLT_Jet150,"",kHLT_Jet150_JetObj);//auto
  mymod->AddTrigger("HLT_Jet190_v0",    kHLT_Jet190, "", kHLT_Jet190_JetObj);
  mymod->AddTrigger("HLT_Jet190_v1",kHLT_Jet190,"",kHLT_Jet190_JetObj);//auto
  mymod->AddTrigger("HLT_Jet190_v2",kHLT_Jet190,"",kHLT_Jet190_JetObj);//auto
  mymod->AddTrigger("HLT_Jet190_v3",kHLT_Jet190,"",kHLT_Jet190_JetObj);//auto
  mymod->AddTrigger("HLT_Jet190_v4",kHLT_Jet190,"",kHLT_Jet190_JetObj);//auto
  mymod->AddTrigger("HLT_Jet190_v5",kHLT_Jet190,"",kHLT_Jet190_JetObj);//auto
  mymod->AddTrigger("HLT_Jet190_v6",kHLT_Jet190,"",kHLT_Jet190_JetObj);//auto
  mymod->AddTrigger("HLT_Jet190_v9",kHLT_Jet190,"",kHLT_Jet190_JetObj);//auto
  mymod->AddTrigger("HLT_Jet240_v0",    kHLT_Jet240, "", kHLT_Jet240_JetObj);
  mymod->AddTrigger("HLT_Jet240_v1",kHLT_Jet240,"",kHLT_Jet240_JetObj);//auto
  mymod->AddTrigger("HLT_Jet240_v2",kHLT_Jet240,"",kHLT_Jet240_JetObj);//auto
  mymod->AddTrigger("HLT_Jet240_v3",kHLT_Jet240,"",kHLT_Jet240_JetObj);//auto
  mymod->AddTrigger("HLT_Jet240_v4",kHLT_Jet240,"",kHLT_Jet240_JetObj);//auto
  mymod->AddTrigger("HLT_Jet240_v5",kHLT_Jet240,"",kHLT_Jet240_JetObj);//auto
  mymod->AddTrigger("HLT_Jet240_v6",kHLT_Jet240,"",kHLT_Jet240_JetObj);//auto
  mymod->AddTrigger("HLT_Jet240_v9",kHLT_Jet240,"",kHLT_Jet240_JetObj);//auto


  mymod->SetPrintHLT(kFALSE); // print HLT table at start of analysis?
  
  ana->AddSuperModule(mymod); 
    
  //
  // run analysis after successful initialisation
  //
  ana->Run(!gROOT->IsBatch());
}

