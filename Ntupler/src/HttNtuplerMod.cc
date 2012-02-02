#include "cstdlib"

#include "TFile.h"
#include "TTree.h"

#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/GenJetCol.h"
#include "MitAna/DataTree/interface/Track.h"
#include "MitAna/DataTree/interface/Vertex.h"
#include "MitAna/DataTree/interface/BeamSpotCol.h"
#include "MitAna/DataTree/interface/PileupInfoCol.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/Photon.h"
#include "MitAna/DataTree/interface/TriggerTable.h"
#include "MitAna/DataTree/interface/TriggerObjectsTable.h"
#include "MitAna/DataTree/interface/TriggerName.h"
#include "MitAna/DataTree/interface/Met.h"
#include "MitAna/DataTree/interface/PFMet.h"
#include "MitAna/DataTree/interface/DecayParticle.h"
#include "MitAna/DataTree/interface/StableData.h"

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitHtt/Ntupler/interface/HttNtuplerMod.hh"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

using namespace mithep;

ClassImp(mithep::HttNtuplerMod)

HttNtuplerMod::HttNtuplerMod(const char *name, const char *title):
  BaseMod        (name,title),
  fOutputFile    ( 0),
  fOutputName    ("ntuple.root"),
  fPartName      (Names::gkMCPartBrn),
  fMCEvtInfoName (Names::gkMCEvtInfoBrn),
  fMuonName      (Names::gkMuonBrn),
  fElectronName  (Names::gkElectronBrn),
  fPrimVtxName   (Names::gkPVBrn),
  fBeamSpotName  (Names::gkBeamSpotBrn),
  fPFJetName     (Names::gkPFJetBrn),
  fPhotonName    (Names::gkPhotonBrn),
  fTrigMaskName  (Names::gkHltBitBrn),
  fPFMetName     ("PFMet"),
  fConversionName(Names::gkMvfConversionBrn),
  fPileupName    (Names::gkPileupInfoBrn),
  fPUEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPFCandidateName(Names::gkPFCandidatesBrn),
  fEmbedWeightName("EmbedWeight"),
  fParticles      ( 0),
  fMCEvtInfo      ( 0),
  fGenJets        ( 0),
  fMuons          ( 0),
  fElectrons      ( 0),
  fPrimVerts      ( 0),
  fBeamSpot       ( 0),
  fPFJets         ( 0),
  fPhotons        ( 0),  
  fTrigMask       ( 0),
  fPFMet          ( 0),
  fConversions    ( 0),
  fPileup         ( 0),
  fPUEnergyDensity( 0),
  fPFCandidates   ( 0),
  fIsData         ( 0),
  fUseGen         ( 0),
  fPrintTable     (kFALSE),
  fSkipIfHLTFail  (kFALSE),
  fMuPtMin        (15),
  fMuPtMax        (1000),
  fMuEtaMin       (-3),
  fMuEtaMax       ( 3),
  fEleEtMin       (15),
  fEleEtMax       (1000),
  fEleEtaMin      (-3),
  fEleEtaMax      ( 3),
  fJetPtMin       (15),
  fPhotonEtMin    (10),
  fMinNTracksFit  ( 0),
  fMinNdof        ( 4),
  fMaxAbsZ        (24),
  fMaxRho         ( 2),
  fEventTree      ( 0),
  fJetCorrector   ( 0),
  fJetUnc         ( 0)
{
  // don't write TObject part of the objects
  TEventInfo::Class()->IgnoreTObjectStreamer();
  TGenInfo::Class()  ->IgnoreTObjectStreamer();
  TVertex::Class()   ->IgnoreTObjectStreamer();
  TSVfit::Class()    ->IgnoreTObjectStreamer();
  TJet::Class()      ->IgnoreTObjectStreamer();
  TMuon::Class()     ->IgnoreTObjectStreamer();
  TPhoton::Class()   ->IgnoreTObjectStreamer();
  TElectron::Class() ->IgnoreTObjectStreamer();
}

HttNtuplerMod::~HttNtuplerMod()
{
}	

void 
HttNtuplerMod::Begin()
{
}

void 
HttNtuplerMod::SlaveBegin()
{
  // request Bambu branches for input
  ReqBranch( fPartName            , fParticles ); 
  ReqBranch( fMCEvtInfoName       , fMCEvtInfo );
  ReqBranch( Names::gkGenJetBrn   , fGenJets   );
  ReqBranch( fMuonName            , fMuons     );
  ReqBranch( fElectronName        , fElectrons );
  ReqBranch( fPrimVtxName         , fPrimVerts );
  ReqBranch( fBeamSpotName        , fBeamSpot  );
  ReqBranch( fPFJetName           , fPFJets    );
  ReqBranch( fTrigMaskName        , fTrigMask  );
  ReqBranch( fPFMetName           , fPFMet     );
  ReqBranch( fPhotonName          , fPhotons   );
  ReqBranch( fConversionName      , fConversions );
  ReqBranch( fPileupName          , fPileup    );
  ReqBranch( fPUEnergyDensityName , fPUEnergyDensity );
  ReqBranch( fPFCandidateName     , fPFCandidates );
  // add embbeding weights for embedded samples 
  if(fIsData==2){ ReqBranch( fEmbedWeightName , fEmbedWeight ); }

  // pileup and no-pileup collections of PFCandidates for delta beta corrected 
  // isolation including and not including the association of tracks to primary 
  // vertices based on their distance along the z-axis
  fPFPileUp      = new PFCandidateOArr;
  fPFNoPileUp    = new PFCandidateOArr;
  fPFPileUpNoZ   = new PFCandidateOArr;
  fPFNoPileUpNoZ = new PFCandidateOArr;
 
  // set up arrays for output
  fMuonArr     = new TClonesArray( "mithep::TMuon"     ); assert( fMuonArr     );
  fElectronArr = new TClonesArray( "mithep::TElectron" ); assert( fElectronArr );
  fPFJetArr    = new TClonesArray( "mithep::TJet"      ); assert( fPFJetArr    );
  fPhotonArr   = new TClonesArray( "mithep::TPhoton"   ); assert( fPhotonArr   );
  fPVArr       = new TClonesArray( "mithep::TVertex"   ); assert( fPVArr       );
  fSVfitEMuArr = new TClonesArray( "mithep::TSVfit"    ); assert( fSVfitEMuArr );
  
  // open the output file
  fOutputFile = new TFile( fOutputName, "RECREATE" );
  
  // setup output trees and structs 
  fEventTree = new TTree( "Events", "Events" );
  if( (fIsData!=1) && fUseGen ){ fEventTree->Branch( "Gen", &fGenInfo ); }
  fEventTree->Branch( "Info"     , &fEventInfo   );
  fEventTree->Branch( "Muon"     , &fMuonArr     );
  fEventTree->Branch( "Electron" , &fElectronArr );
  fEventTree->Branch( "PFJet"    , &fPFJetArr    );
  fEventTree->Branch( "Photon"   , &fPhotonArr   );
  fEventTree->Branch( "PV"       , &fPVArr       );
  fEventTree->Branch( "SVfitEMu" , &fSVfitEMuArr );

  // setup jet corrections for article flow jets
  std::vector<JetCorrectorParameters> correctionParameters;
  for(unsigned int icorr=0; icorr<fJetCorrParsv.size(); icorr++){
    correctionParameters.push_back(JetCorrectorParameters(fJetCorrParsv[icorr].Data()));
  }
  // initialize jet corrector class
  fJetCorrector = new FactorizedJetCorrector(correctionParameters);
  // initialize jet uncertainties in analogy to the jes corrections 
  // the jec uncertainties are returned as relative UP/DOWN shifts, 
  // so if scaledJet were a reco::Jet the uncertainty shift should 
  // be applied as scaledJet.scaleEnergy(1+/-uncert). For more details 
  // have a look into the class itself in: CondFormats/JetMETObjects
  char* PATH = getenv("CMSSW_BASE"); assert(PATH);
  JetCorrectorParameters param(std::string(TString::Format("%s/src/MitPhysics/data/START42_V12_AK5PF_Uncertainty.txt", PATH)));
  fJetUnc = new JetCorrectionUncertainty(param);
  // initialize MET significance for svfit
  fMetSignificance = new MetSignificance();
  // initialize tools for electron ID
  fEleTools = new ElectronTools();
  // setup selection with JSON file, if necessary
  fhasJSON = fJSONv.size() > 0;
  for(unsigned int idx=0; idx<fJSONv.size(); idx++) {
    frlrm.AddJSONFile(fJSONv[idx].Data());
  }
  fLastRunLumi = RunLumiRangeMap::RunLumiPairType(0,0);
}

void 
HttNtuplerMod::SlaveTerminate()
{
  // save output to file and close
  fEventTree ->Print();
  fOutputFile->Write();
  fOutputFile->Close();

  // cleanup
  delete fPFPileUp;
  delete fPFNoPileUp;
  delete fPFPileUpNoZ;
  delete fPFNoPileUpNoZ;
  delete fMuonArr;
  delete fElectronArr;
  delete fPFJetArr;
  delete fPhotonArr;
  delete fPVArr;
  delete fSVfitEMuArr;
  delete fJetCorrector;
  delete fJetUnc;
  delete fMetSignificance;
  delete fEleTools;

  // dump json file
  TString jsonfname = fOutputName.ReplaceAll("root","json");
  if(fIsData){
    fRunLumiSet.DumpJSONFile(jsonfname.Data());
  }
}

void 
HttNtuplerMod::Terminate()
{
}

void 
HttNtuplerMod::BeginRun()
{
  if(HasHLTInfo() && fPrintTable) { GetHLTTable()->Print(); }
}

void 
HttNtuplerMod::EndRun()
{
}

void 
HttNtuplerMod::Process()
{
  RunLumiRangeMap::RunLumiPairType rl(GetEventHeader()->RunNum(), GetEventHeader()->LumiSec());
  if(fhasJSON && !frlrm.HasRunLumi(rl)) return;  // not certified run? Skip to next event...
  if(rl!=fLastRunLumi) {
    fLastRunLumi = rl;
    fRunLumiSet.Add(rl);
  }
  // load branches from Bambu
  LoadBranch( fMuonName            );
  LoadBranch( fElectronName        );
  LoadBranch( fPrimVtxName         );
  LoadBranch( fBeamSpotName        );
  LoadBranch( fPFJetName           );
  LoadBranch( fTrigMaskName        );
  LoadBranch( fPFMetName           ); 
  LoadBranch( fPhotonName          );
  LoadBranch( fConversionName      );
  LoadBranch( fPUEnergyDensityName );
  LoadBranch( fPFCandidateName     );
  // branches for monte carlo
  if( fIsData==0 ){ 
    LoadBranch( fPartName          );
    LoadBranch( fPileupName        );
    LoadBranch( Names::gkGenJetBrn );
    if(fUseGen){ 
      LoadBranch(fMCEvtInfoName);
      if(fUseGen==1){ FillGenH();  }
      if(fUseGen==2){ FillGenZ();  }
      if(fUseGen==3){ FillGenW();  }
      if(fUseGen==4){ FillGenWW(); }
    }    
  }
  // branches for embedded sample(s)
  if( fIsData==2 ){ 
    LoadBranch(fEmbedWeightName);
    LoadBranch(fPartName);
  }

  // get HLT info: trigger objects can be matched by name to the corresponding trigger 
  // that has been passed. NOTE: TriggerName::Id() is the bambu numbering scheme, while 
  // fTriggerIdsv[itrig] is an internan numbering scheme
  ULong64_t trigbits=0;
  if( HasHLTInfo() ){
    const TriggerTable *hltTable = GetHLTTable(); assert(hltTable);
    for(unsigned int itrig=0; itrig<fTriggerNamesv.size(); ++itrig){
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      // if event passed this trigger, set the trigger's bit in trigbits
      if(fTrigMask->At(trigname->Id())) { trigbits |= fTriggerIdsv[itrig]; }
    }  
  }
  if(fSkipIfHLTFail && (trigbits==0)){ return; }
  
  // increment the number of events that have been processed
  IncNEventsProcessed();

  // clear internal object arrays
  fMuonArr    ->Clear();
  fElectronArr->Clear();
  fPFJetArr   ->Clear();
  fPhotonArr  ->Clear();
  fPVArr      ->Clear();
  fSVfitEMuArr->Clear();
  
  // get beam spot information. If no beam spot information is available, coordinates are set to 99999.
  double bsx=99999., bsy=99999., bsz=99999.;
  if(fBeamSpot){
    if(fBeamSpot->GetEntries()>1){
      std::cout << "********** Found more than 1 beam spot choosing the first one! **********" << std::endl;
    }
    const BeamSpot* bs = fBeamSpot->At(0);
    bsx = bs->X();
    bsy = bs->Y();
    bsz = bs->Z();
  }

  // get primary vertices. Assume that the primary vertices are ordered 
  // by sum-pt**2 (as should be in CMSSW) NOTE: if no PV is found from 
  // fitting tracks, the beamspot is used. bestPV is used for the calcu-
  // lation of d0 wrt to the recnstructed PV. It is taken to be the PV
  // with the highest sum-pt**2 that fullfils the selection requirements
  // if possible and the PV with the highest sum-pt**2 else
  const Vertex* bestPV = 0;
  bool hasGoodPV = false;  
  for(unsigned int i=0; i<fPrimVerts->GetEntries(); ++i){
    const Vertex* pv = fPrimVerts->At(i);
    // select best PV for d0 wrt to reconstructed PV; if no PV is passing cuts, 
    // the first PV in the collection will be used
    //if(!pv->IsValid()) continue;
    if( pv->NTracksFit()     < fMinNTracksFit ){ continue; }
    if( pv->Ndof()	     < fMinNdof       ){ continue; }
    if( fabs(pv->Z())        > fMaxAbsZ       ){ continue; }
    if( pv->Position().Rho() > fMaxRho        ){ continue; }
    hasGoodPV = true;
    FillPV(pv);
    if(!bestPV){ bestPV = pv; }
  }
  if(!bestPV) bestPV = fPrimVerts->At(0);
  fVertex = bestPV;

  // create pileUp and noPileUp collections with association using along z axis 
  separatePileUp(fPFCandidates, fVertex, fPrimVerts, fPFPileUp   , fPFNoPileUp   , true );
  // create pileUp and noPileUp collections w/o association using along z axis 
  separatePileUp(fPFCandidates, fVertex, fPrimVerts, fPFPileUpNoZ, fPFNoPileUpNoZ, false);
  
  // prepare met significance matrix and fill all information needed to do SVfit 
  for(unsigned int i0=0; i0<fMuons->GetEntries(); ++i0){
    const Muon* pMu = fMuons->At(i0);
    if( !looseMuId(pMu) ) continue;
    for(unsigned int i1=0; i1<fElectrons->GetEntries(); ++i1){
      const Electron* pElectron = fElectrons->At(i1);
      if( !looseEleId(pElectron) ) continue;
      if( MathUtils::DeltaR(pElectron->Mom(),pMu->Mom())<0.3 ) continue;
      TMatrixD lMetMatrix = fMetSignificance->getSignificance(fPFJets, fPFCandidates, 0, pMu, pElectron);
      FillSVfit(fSVfitEMuArr, (Particle*)pMu, (Particle*)pElectron, lMetMatrix);
    }
  }

  // loop muons
  assert(fMuons);
  for(unsigned int i=0; i<fMuons->GetEntries(); ++i){
    const Muon* mu=fMuons->At(i); 
    if(!mu->HasTrk()) continue; 
    // use tracker tracks for kinematics when available
    const Track* muTrk=0;
    if(mu->HasTrackerTrk()) { muTrk = mu->TrackerTrk(); }
    else if(mu->HasGlobalTrk()) { muTrk = mu->GlobalTrk(); }
    else { muTrk = mu->StandaloneTrk(); } 
    // kinematic selection      
    if((muTrk->Eta() < fMuEtaMin) || (muTrk->Eta() > fMuEtaMax)) continue;   
    if((muTrk->Pt()  > fMuPtMax)  || (muTrk->Pt()  < fMuPtMin))  continue;
    FillMuon(mu);  
  }
       
  // loop electrons.
  assert(fElectrons);
  for(unsigned int i=0; i<fElectrons->GetEntries(); ++i){
    const Electron* ele=fElectrons->At(i);  
    if((ele->Pt()  < fEleEtMin)  || (ele->Pt()  > fEleEtMax )) continue;
    if((ele->Eta() < fEleEtaMin) || (ele->Eta() > fEleEtaMax)) continue;
    if(!fEleTools->PassSpikeRemovalFilter(ele)) continue;
    if(fIsData==2 && !ele->BestTrk()) continue;
    FillElectron(ele);
  }
  
  // loop jets
  assert(fPFJets);
  for(unsigned int i=0; i<fPFJets->GetEntries(); ++i){
    const PFJet* jet = fPFJets->At(i);
    const FourVectorM rawMom = jet->RawMom();
    // initializa jet corrector
    fJetCorrector->setJetEta( rawMom.Eta() );
    fJetCorrector->setJetPt ( rawMom.Pt()  );
    fJetCorrector->setJetPhi( rawMom.Phi() );
    fJetCorrector->setJetE  ( rawMom.E()   );
    fJetCorrector->setRho   ( fPUEnergyDensity->At(0)->RhoHighEta() );
    fJetCorrector->setJetA  ( jet->JetArea() );
    fJetCorrector->setJetEMF( -99.0 );     
    
    // keep all jets above specified threshold (after energy correction)
    // and all jets with valid b-tag value (Track Counting High Efficiency 
    // method default is -100)
    double correction = fJetCorrector->getCorrection();
    double pt = rawMom.Pt()*correction;
    if(pt > fJetPtMin || (jet->TrackCountingHighEffBJetTagsDisc() != -100)){
      // loose jetId for particle flow jets. NOTE: energy fractions are on the 
      // raw (i.e. uncorrected energy)
      if(jet->E()==0) continue;
      if(jet->NeutralHadronEnergy()/jet->E()   >  0.99) continue;  //   'neutralHadronEnergyFraction < 0.99'
      if(jet->NeutralEmEnergy()/jet->E()       >  0.99) continue;  //   'neutralEmEnergyFraction < 0.99'
      if(jet->NConstituents()                  <  2   ) continue;  //   'nConstituents > 1'
      if(jet->ChargedHadronEnergy()/jet->E() <= 0  && fabs(jet->Eta()) < 2.4 )  continue;
      if(jet->ChargedEmEnergy()/jet->E()  >  0.99  && fabs(jet->Eta()) < 2.4 )  continue;
      if(jet->ChargedMultiplicity()    < 1         && fabs(jet->Eta()) < 2.4 )  continue;
      FillJet(jet);
    }
  }
    
  // loop photons
  assert(fPhotons);
  for(unsigned int i=0; i<fPhotons->GetEntries(); ++i){
    const Photon* pho= fPhotons->At(i);
    if(pho->SCluster()->Et() > fPhotonEtMin) { FillPhoton(pho); }
  } 
  
  // fill pfMET
  TLorentzVector pfmet; pfmet.SetPxPyPzE(fPFMet->At(0)->Mex(),fPFMet->At(0)->Mey(),0,0);

  // fill trkMET from particlew flow candidates
  double trkMetx=0, trkMety=0, trkSumET=0;
  assert(fPFCandidates);
  for(unsigned int i=0; i<fPFCandidates->GetEntries(); ++i){
    const Double_t trkDzCut = 0.1;
    const PFCandidate *pfcand = fPFCandidates->At(i);
    if( (pfcand->HasTrackerTrk() && (fabs(pfcand->TrackerTrk()->DzCorrected(*fVertex))<trkDzCut)) || (pfcand->HasGsfTrk()     && (fabs(pfcand->GsfTrk()->DzCorrected(*fVertex))<trkDzCut)) ) { 
      trkMetx  -= pfcand->Px();
      trkMety  -= pfcand->Py();
      trkSumET += pfcand->Pt();
    }
  }
  TLorentzVector trkmet; trkmet.SetPxPyPzE(trkMetx,trkMety,0,0);

  // get number of PU interactions for MC events (also available 
  // for embedded samples?)
  int npu = -1;
  if( !fIsData ){
    for(unsigned int i=0; i<fPileup->GetEntries(); ++i){
      if(fPileup->At(i)->GetBunchCrossing() == 0)
	npu = fPileup->At(i)->GetPU_NumInteractions();
    }
    assert(npu>=0);
  } 

  // fill event info
  fEventInfo.runNum       = GetEventHeader()->RunNum();
  fEventInfo.evtNum       = GetEventHeader()->EvtNum();
  fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
  fEventInfo.nPU          = fIsData ? 0 : npu;
  fEventInfo.triggerBits  = trigbits;
  fEventInfo.pvx          = fVertex->X();
  fEventInfo.pvy          = fVertex->Y();
  fEventInfo.pvz          = fVertex->Z();
  fEventInfo.bsx          = bsx;
  fEventInfo.bsy          = bsy;
  fEventInfo.bsz          = bsz;
  fEventInfo.pfMET        = pfmet.Pt();
  fEventInfo.pfMETphi     = pfmet.Phi();
  fEventInfo.pfSumET      = fPFMet->At(0)->SumEt();
  fEventInfo.trkMET       = trkmet.Pt();
  fEventInfo.trkMETphi    = trkmet.Phi();
  fEventInfo.trkSumET     = trkSumET;
  fEventInfo.rho          = fPUEnergyDensity->At(0)->RhoHighEta();
  fEventInfo.hasGoodPV    = hasGoodPV;
  fEventInfo.embWeight    = (fIsData==2) ? fEmbedWeight->At(fEmbedWeight->GetEntries()-1)->Weight() : 1;

  // fill the tree
  fEventTree->Fill();
}

void 
HttNtuplerMod::FillMuon(const Muon* mu)
{
  assert(mu);
  TClonesArray& rMuonArr = *fMuonArr; 
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const int index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon* pMuon = (TMuon*)rMuonArr[index];
  
  // use tracker track when available
  const Track* muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk();    }
  else if(mu->HasGlobalTrk())     { muTrk = mu->GlobalTrk();     }
  else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); }
  assert(muTrk);                  
  
  pMuon->pt       = muTrk->Pt();
  pMuon->ptErr    = muTrk->PtErr();
  pMuon->eta      = muTrk->Eta();
  pMuon->phi      = muTrk->Phi();
  pMuon->trkIso03 = mu->IsoR03SumPt();
  pMuon->emIso03  = mu->IsoR03EmEt();
  pMuon->hadIso03 = mu->IsoR03HadEt();
  pMuon->pfIso03  = computePFMuonIso(mu,0.3);
  pMuon->pfIso04  = computePFMuonIso(mu,0.4);  
  pMuon->d0       = muTrk->D0Corrected(*fVertex);
  pMuon->dz       = muTrk->DzCorrected(*fVertex);
  pMuon->tkNchi2  = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->RChi2() : 0;
  if(mu->HasGlobalTrk())          { pMuon->muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { pMuon->muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { pMuon->muNchi2 = mu->TrackerTrk()->RChi2();    }
  pMuon->q          = muTrk->Charge();
  pMuon->nValidHits = mu->NValidHits();
  pMuon->qualityBits = mu->Quality().QualityMask().Mask();

  // NOTE: It is possible for a muon to be TK+SA. The muon reco associates a TK with a SA
  // if chamber matches for the TK and hits for the SA share DetIDs see hypernews thread:
  // https://hypernews.cern.ch/HyperNews/CMS/get/csa08-muons/57/2/1/1/1.html
  //	      
  pMuon->typeBits = 0;
  if(mu->IsGlobalMuon())     { pMuon->typeBits |= kGlobal; }
  if(mu->IsTrackerMuon())    { pMuon->typeBits |= kTracker; }
  if(mu->IsStandaloneMuon()) { pMuon->typeBits |= kStandalone; }
  pMuon->nTkHits      = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->NHits() : 0;
  pMuon->nPixHits     = muTrk->NPixelHits();
  pMuon->nSeg         = mu->NSegments();
  pMuon->nMatch       = mu->NMatches();
  pMuon->hltMatchBits = MatchHLT(muTrk->Pt(),muTrk->Eta(),muTrk->Phi());
  pMuon->staPt  = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Pt()  : -1;  
  pMuon->staEta = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Eta() : -999;
  pMuon->staPhi = (mu->HasStandaloneTrk()) ? mu->StandaloneTrk()->Phi() : -999;
  pMuon->trkID  = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->GetUniqueID() : 0;

  pMuon->pfPx=0; pMuon->pfPy=0;
  assert(fPFCandidates);
  for(unsigned int i=0; i<fPFCandidates->GetEntries(); ++i){    
    const PFCandidate* pfcand = fPFCandidates->At(i);
    if(mu->HasTrackerTrk() && mu->TrackerTrk() == pfcand->TrackerTrk()){
      pMuon->pfPx = pfcand->Px();
      pMuon->pfPy = pfcand->Py();
      break;
    }
  }

  // Mike's PF isolation
  pMuon->pfIsoCharged    = PFIsoNoPileup(mu, fPFNoPileUp    , 0.0, 0.4, 0.0001,  1);
  pMuon->pfIsoChargedNoZ = PFIsoNoPileup(mu, fPFNoPileUpNoZ , 0.0, 0.4, 0.0001,  1);
  pMuon->pfIsoNeutral    = PFIsoNoPileup(mu, fPFNoPileUp    , 0.5, 0.4,   0.01,  2);
  pMuon->pfIsoNeutralNoZ = PFIsoNoPileup(mu, fPFNoPileUpNoZ , 0.5, 0.4,   0.01,  2);
  pMuon->pfIsoGamma      = PFIsoNoPileup(mu, fPFNoPileUp    , 0.5, 0.4,   0.01,  3);
  pMuon->pfIsoGammaNoZ   = PFIsoNoPileup(mu, fPFNoPileUpNoZ , 0.5, 0.4,   0.01,  3);
  pMuon->puIso    = PileupIsolation(mu, fPFPileUp    , 0.5, 0.4, 0.01);
  pMuon->puIsoNoZ = PileupIsolation(mu, fPFPileUpNoZ , 0.5, 0.4, 0.01);
}


void 
HttNtuplerMod::FillElectron(const Electron* ele)
{
  assert(ele);
  TClonesArray& rElectronArr = *fElectronArr;
  assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
  const int index = rElectronArr.GetEntries();  
  new(rElectronArr[index]) TElectron();
  TElectron* pElectron = (TElectron*)rElectronArr[index];
                    
  pElectron->pt              = ele->Pt();
  pElectron->eta             = ele->Eta();
  pElectron->phi             = ele->Phi();
  pElectron->trkIso03        = ele->TrackIsolationDr03();
  pElectron->emIso03         = ele->EcalRecHitIsoDr03();
  pElectron->hadIso03        = ele->HcalTowerSumEtDr03();
  pElectron->pfIso03         = computePFElecIso(ele,0.3); 
  pElectron->pfIso04         = computePFElecIso(ele,0.4);
  pElectron->d0              = ele->BestTrk()->D0Corrected(*fVertex);
  pElectron->dz              = ele->BestTrk()->DzCorrected(*fVertex);  
  pElectron->scEt            = ele->SCluster()->Et();
  pElectron->scEta           = ele->SCluster()->Eta();
  pElectron->scPhi           = ele->SCluster()->Phi();
  pElectron->HoverE          = ele->HadronicOverEm();
  pElectron->EoverP          = ele->ESuperClusterOverP();
  pElectron->fBrem           = ele->FBrem();
  pElectron->deltaEtaIn      = ele->DeltaEtaSuperClusterTrackAtVtx();
  pElectron->deltaPhiIn      = ele->DeltaPhiSuperClusterTrackAtVtx();
  pElectron->sigiEtaiEta     = ele->CoviEtaiEta();
  pElectron->nExpHitsInner   = ele->BestTrk()->NExpectedHitsInner();
  pElectron->partnerDeltaCot = ele->ConvPartnerDCotTheta();
  pElectron->partnerDist     = ele->ConvPartnerDist();
  pElectron->q               = ele->Charge(); 
  pElectron->p               = ele->BestTrk()->P();
  pElectron->ip3d            = ele->Ip3dPV();
  pElectron->ip3dSig         = ele->Ip3dPVSignificance();
  pElectron->ESeedClusterOverPIn   = ele->ESeedClusterOverPIn();
  pElectron->ESeedClusterOverPOut  = ele->ESeedClusterOverPout();
  pElectron->sigiPhiiPhi     = ele->SCluster()->Seed()->CoviPhiiPhi();
  pElectron->nBrem           = ele->NumberOfClusters()-1;  
  pElectron->hltMatchBits    = MatchHLT(ele->SCluster()->Eta(),ele->SCluster()->Phi());
  pElectron->scID            = ele->SCluster()->GetUniqueID();
  pElectron->trkID           = (ele->HasTrackerTrk()) ? ele->TrackerTrk()->GetUniqueID() : 0;
  pElectron->isEcalDriven    = ele->IsEcalDriven();
  pElectron->isConv          = IsConversion(ele);

  pElectron->pfPx=0; pElectron->pfPy=0;
  assert(fPFCandidates);
  for(unsigned int i=0; i<fPFCandidates->GetEntries(); ++i){
    const PFCandidate* pfcand = fPFCandidates->At(i);
    if( (pfcand->HasTrackerTrk() && ele->TrackerTrk() == pfcand->TrackerTrk()) || (pfcand->HasGsfTrk() && ele->GsfTrk() == pfcand->GsfTrk()) ){
      pElectron->pfPx = pfcand->Px();
      pElectron->pfPy = pfcand->Py();
      break;
    }	 
  }
  
  // Mike's PF isolation
  pElectron->pfIsoCharged    = PFIsoNoPileup(ele, fPFNoPileUp    , 0.0, 0.4, 0.015,  1);
  pElectron->pfIsoChargedNoZ = PFIsoNoPileup(ele, fPFNoPileUpNoZ , 0.0, 0.4, 0.015,  1);
  pElectron->pfIsoNeutral    = PFIsoNoPileup(ele, fPFNoPileUp    , 0.5, 0.4,   0.0,  2);
  pElectron->pfIsoNeutralNoZ = PFIsoNoPileup(ele, fPFNoPileUpNoZ , 0.5, 0.4,   0.0,  2);
  pElectron->pfIsoGamma      = PFIsoNoPileup(ele, fPFNoPileUp    , 0.5, 0.4,  0.08,  3);
  pElectron->pfIsoGammaNoZ   = PFIsoNoPileup(ele, fPFNoPileUpNoZ , 0.5, 0.4,  0.08,  3);
  pElectron->puIso    = PileupIsolation(ele, fPFPileUp    , 0.5, 0.4, 0.01);
  pElectron->puIsoNoZ = PileupIsolation(ele, fPFPileUpNoZ , 0.5, 0.4, 0.01);
}

void 
HttNtuplerMod::FillJet(const PFJet* jet)
{
  TClonesArray& rPFJetArr = *fPFJetArr;
  assert(rPFJetArr.GetEntries() < rPFJetArr.GetSize());
  const int index = rPFJetArr.GetEntries();  
  new(rPFJetArr[index]) TJet();
  TJet* pPFJet = (TJet*)rPFJetArr[index]; 
  
  const FourVectorM rawMom = jet->RawMom();
  fJetCorrector->setJetEta( rawMom.Eta() );
  fJetCorrector->setJetPt ( rawMom.Pt()  );
  fJetCorrector->setJetPhi( rawMom.Phi() );
  fJetCorrector->setJetE  ( rawMom.E()   );
  fJetCorrector->setRho   ( fPUEnergyDensity->At(0)->RhoHighEta() );
  fJetCorrector->setJetA  ( jet->JetArea() );
  fJetCorrector->setJetEMF( -99.0 );

  fJetUnc->setJetPt ( rawMom.Pt()  );
  fJetUnc->setJetEta( rawMom.Eta() );
  // this is an up shift, down shift have false as argument
  float hierr  = fJetUnc->getUncertainty(true);

  pPFJet->pt   = (rawMom.Pt())*(fJetCorrector->getCorrection());
  pPFJet->eta  = rawMom.Eta();
  pPFJet->phi  = rawMom.Phi();
  pPFJet->mass = jet->Mass();
  pPFJet->unc  = hierr;
  pPFJet->area = jet->JetArea();
  pPFJet->tche = jet->TrackCountingHighEffBJetTagsDisc();
  pPFJet->tchp = jet->TrackCountingHighPurBJetTagsDisc();
  pPFJet->mcFlavor = jet->MatchedMCFlavor();

  int matchedFlavor = -999;
  if( fGenJets ){
    double dRmin = 0.3;
    for(unsigned int i=0; i<fGenJets->GetEntries(); ++i){
      const GenJet* j = fGenJets->At(i);
      if(MathUtils::DeltaR(*jet,*j) < dRmin) {
        dRmin = MathUtils::DeltaR(*jet,*j);
        matchedFlavor = j->MatchedMCFlavor();
      }
    }
  }
  pPFJet->matchedId = matchedFlavor;
  pPFJet->hltMatchBits = MatchHLT(jet->Eta(),jet->Phi());
}

void 
HttNtuplerMod::FillPhoton(const Photon *pho)
{
  TClonesArray& rPhotonArr = *fPhotonArr;
  assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
  const int index = rPhotonArr.GetEntries();  
  new(rPhotonArr[index]) TPhoton();
  TPhoton* pPhoton = (TPhoton*)rPhotonArr[index];
  
  pPhoton->pt		= pho->Pt(); 
  pPhoton->eta  	= pho->Eta();
  pPhoton->phi  	= pho->Phi();
  pPhoton->scEt		= pho->SCluster()->Et(); 
  pPhoton->scEta  	= pho->SCluster()->Eta();
  pPhoton->scPhi  	= pho->SCluster()->Phi();
  pPhoton->trkIso04     = pho->HollowConeTrkIsoDr04();
  pPhoton->emIso04      = pho->EcalRecHitIsoDr04();
  pPhoton->hadIso04	= pho->HcalTowerSumEtDr04(); 
  pPhoton->HoverE	= pho->HadOverEm();
  pPhoton->R9		= pho->R9();
  pPhoton->sigiEtaiEta  = pho->CoviEtaiEta();
  pPhoton->hltMatchBits = MatchHLT(pho->SCluster()->Eta(),pho->SCluster()->Phi());
  pPhoton->scID         = pho->SCluster()->GetUniqueID();
  pPhoton->hasPixelSeed = pho->HasPixelSeed();
}

void 
HttNtuplerMod::FillPV(const Vertex *pv) 
{
  TClonesArray& rPVArr = *fPVArr;
  assert(rPVArr.GetEntries() < rPVArr.GetSize());
  const int index = rPVArr.GetEntries();  
  new(rPVArr[index]) TVertex();
  TVertex* pVertex = (TVertex*)rPVArr[index];
  
  pVertex->nTracksFit = pv->NTracksFit();
  pVertex->ndof       = pv->Ndof();      
  pVertex->chi2       = pv->Chi2();  
  pVertex->x          = pv->X();
  pVertex->y          = pv->Y();
  pVertex->z          = pv->Z();
 
  pVertex->sumPt=0;
  for(unsigned int itrk=0; itrk<pv->NTracks(); itrk++){
    pVertex->sumPt += pv->Trk(itrk)->Pt();				   
  }
}

ULong64_t HttNtuplerMod::MatchHLT(const double eta, const double phi)
{
  ULong64_t bits = 0;
  const double hltMatchR = 0.2;
  
  if( HasHLTInfo() ){
    const TriggerTable* hltTable = GetHLTTable(); assert(hltTable);
    for(unsigned int itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList* list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue; TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
      
      while( to ){             
        if( to->IsHLT() ){          
	  if(fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0){
	    bool match = true;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=false;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=false;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
	  }
	  if(fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0){
	    bool match = true;
	    if(to->Pt() < fTriggerObjMinPt2v[itrig])                       match=false;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=false;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds2v[itrig];
	  }
	  if(fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0){
	    bool match = true;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=false;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=false;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
          }
        } 
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  return bits;
}

ULong64_t HttNtuplerMod::MatchHLT(const double pt, const double eta, const double phi)
{
  ULong64_t bits = 0;
  const Double_t hltMatchR = 0.2;
  const Double_t hltMatchPtFrac = 1;
  
  if( HasHLTInfo() ){
    const TriggerTable *hltTable = GetHLTTable(); assert(hltTable);
    for(unsigned int itrig=0; itrig<fTriggerNamesv.size(); ++itrig){
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList* list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue; TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while( to ){         
        if( to->IsHLT() ){
	  if(fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0){
	    bool match = true;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=false;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=false;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=false;  // pT matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
	  }
	  if(fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0){
	    bool match = true;
	    if(to->Pt() < fTriggerObjMinPt2v[itrig])                       match=false;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=false;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=false;  // pT matching
	    if(match) bits |= fTriggerObjIds2v[itrig];
	  }
	  if(fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0){
	    bool match = true;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=false;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=false;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=false;  // pT matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
          }
        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  return bits;
}

bool 
HttNtuplerMod::IsConversion(const Electron* ele) 
{
  bool isGoodConversion = false;
  const unsigned int nWrongHitsMax = 0;
  const double probMin = 1e-6;
  const double lxyMin = 2.0;
  const bool matchCkf = true;
  const bool requireArbitratedMerged = false;
  
  for(unsigned int ifc=0; ifc<fConversions->GetEntries(); ++ifc){
    bool ConversionMatchFound = false;
    for(unsigned int d=0; d<fConversions->At(ifc)->NDaughters(); ++d){
      const Track* trk = dynamic_cast<const ChargedParticle*>(fConversions->At(ifc)->Daughter(d))->Trk();
      if( ele->GsfTrk() == trk || (matchCkf && ele->TrackerTrk()==trk) ){
        ConversionMatchFound = true;
        break;
      }
    }
    // if match between the e-track and one of the conversion legs
    if( ConversionMatchFound ){
      isGoodConversion = (fConversions->At(ifc)->Prob() > probMin) && (!requireArbitratedMerged || fConversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) && (fConversions->At(ifc)->LxyCorrected((BaseVertex*)fVertex) > lxyMin);
      if (isGoodConversion == true) {
        for (unsigned int d=0; d<fConversions->At(ifc)->NDaughters(); ++d){
          const Track* trk = dynamic_cast<const ChargedParticle*> (fConversions->At(ifc)->Daughter(d))->Trk();
          if( trk ){
            const StableData* sd = dynamic_cast<const StableData*> (fConversions->At(ifc)->DaughterDat(d));
            if( sd->NWrongHits() > nWrongHitsMax ){
              isGoodConversion = false;
	    }
          } 
	  else {
            isGoodConversion = false;
          }
        }
      }
    }
    if(isGoodConversion == true){
      break;
    }
  }
  return isGoodConversion;
}

float 
HttNtuplerMod::computePFMuonIso(const Muon* muon, const double dRMax)
{
  const double dRMin = 0;
  const double neuPtMin = 1.0;
  const double dzMax = 0.1;
  double zLepton = (muon->BestTrk()) ? muon->BestTrk()->DzCorrected(*fVertex) : 0.0;
  float iso=0;
  for(unsigned int ipf=0; ipf<fPFCandidates->GetEntries(); ++ipf){
    const PFCandidate* pfcand = fPFCandidates->At(ipf);
    if(!pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin)) continue;  // pT cut on neutral particles
    // exclude the muon itself from the list of particle flow candidates
    if(pfcand->TrackerTrk() && muon->TrackerTrk() && (pfcand->TrackerTrk()==muon->TrackerTrk())) continue;
    // dz cut
    double dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(*fVertex) - zLepton) : 0;
    if(dz >= dzMax) continue;
    // consider only those particles within the isolation cone
    double dr = MathUtils::DeltaR(muon->Mom(), pfcand->Mom());
    if(dr<dRMax && dr>=dRMin){
      iso += pfcand->Pt(); 
    }
  }
  return iso;
}

float 
HttNtuplerMod::computePFElecIso(const Electron* electron, const double dRMax)
{
  const double dRMin = 0;
  const double neuPtMin = 1.0;
  const double dzMax = 0.1;
  double zLepton = (electron->BestTrk()) ? electron->BestTrk()->DzCorrected(*fVertex) : 0.0;
  float iso=0;
  for(unsigned int ipf=0; ipf<fPFCandidates->GetEntries(); ++ipf){
    const PFCandidate* pfcand = fPFCandidates->At(ipf);
    if(!pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin)) continue;  // pT cut on neutral particles
    // dz cut
    double dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(*fVertex) - zLepton) : 0;
    if(dz >= dzMax) continue;
    // exclude the electron itself from the list of particle flow candidates
    if(pfcand->TrackerTrk() && electron->TrackerTrk() && (pfcand->TrackerTrk()==electron->TrackerTrk())) continue;
    if(pfcand->GsfTrk()     && electron->GsfTrk()     && (pfcand->GsfTrk()==electron->GsfTrk())) continue;
    // check iso cone
    Double_t dr = MathUtils::DeltaR(electron->Mom(), pfcand->Mom());
    if(dr<dRMax && dr>=dRMin) {
      // eta-strip veto for photons
      if((pfcand->PFType() == PFCandidate::eGamma) && fabs(electron->Eta() - pfcand->Eta()) < 0.025) continue;
      // inner cone (one tower = dR < 0.07) veto for non-photon neutrals
      if(!pfcand->HasTrk() && (pfcand->PFType() == PFCandidate::eNeutralHadron) && (MathUtils::DeltaR(electron->Mom(), pfcand->Mom()) < 0.07)) continue;
      iso += pfcand->Pt();
    }
  }
  return iso;
}

double 
HttNtuplerMod::PFIsoNoPileup(const ChargedParticle* p, PFCandidateCol* pfNoPileUp, double ptMin, double extRadius, double intRadius, int isoType)
{
  assert(p);
  double ptSum = 0.0;
  for(unsigned int i=0; i<pfNoPileUp->GetEntries(); ++i){
    const PFCandidate* pf = pfNoPileUp->At(i); assert(pf);
    PFCandidate::EPFType pfType = pf->PFType();
    if((isoType == 1 && (pfType == PFCandidate::eHadron || pfType == PFCandidate::eElectron || pfType == PFCandidate::eMuon)) || (isoType == 2 && pfType == PFCandidate::eNeutralHadron) || (isoType == 3 && pfType == PFCandidate::eGamma)){
      if(pf->Pt() >= ptMin && !(pf->TrackerTrk() && p->TrackerTrk() && pf->TrackerTrk() == p->TrackerTrk())){
	// add pt to running sum if the particle flow candidate is close enough in deltaR
	double dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
	if(dr < extRadius && dr >= intRadius) ptSum += pf->Pt();
      }
    }
  }
  return ptSum;
}

double 
HttNtuplerMod::PileupIsolation(const Particle* p, PFCandidateCol* pfPileUp, double ptMin, double extRadius, double intRadius)
{   
  assert(pfPileUp);
  double ptSum = 0.0;
  for(unsigned int i=0; i<pfPileUp->GetEntries(); ++i){ 
    const PFCandidate* pf = pfPileUp->At(i); assert(pf);
    if(pf->Pt() >= ptMin){
      // add pt to running sum if particle flow candidate is within the isolation cone
      double dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
      if(dr < extRadius && dr >= intRadius) ptSum += pf->Pt();
    }
  }
  return ptSum;
}

void 
HttNtuplerMod::separatePileUp(const PFCandidateCol* pfCandidates, const Vertex* pv, const VertexCol* primVerts, PFCandidateOArr* pfPileUp, PFCandidateOArr* pfNoPileUp, bool checkClosestZVertex)
{
  assert(pfPileUp); assert(pfNoPileUp);
  pfPileUp->Reset(); pfNoPileUp->Reset();
  for(unsigned int i=0; i<pfCandidates->GetEntries(); ++i){
    const PFCandidate* pf = pfCandidates->At(i); assert(pf);
    if(pf->PFType() == PFCandidate::eHadron){
      if(pf->HasTrackerTrk() && pv->HasTrack(pf->TrackerTrk()) /*&& pv->TrackWeight(pf->TrackerTrk()) > 0*/){
	pfNoPileUp->Add(pf);
      }
      else{
        bool vertexFound = false;
        const Vertex* closestVtx = 0;
        double dzmin = 10000;
        for(unsigned int j=0; j<primVerts->GetEntries(); ++j){
          const Vertex* vtx = primVerts->At(j); assert(vtx);
          if(pf->HasTrackerTrk() && vtx->HasTrack(pf->TrackerTrk()) /*&& vtx->TrackWeight(pf->TrackerTrk()) > 0*/){
            vertexFound = true;
            closestVtx = vtx;
            break;
          }
          double dz = fabs(pf->SourceVertex().Z() - vtx->Z());
          if(dz < dzmin){
            closestVtx = vtx;
            dzmin = dz;
          }
        }
        if(checkClosestZVertex){
          // Fallback: if track is not associated with any vertex
          // associate it with the vertex, which is closest in z
          if(vertexFound || closestVtx != pv){
            pfPileUp->Add(pf);
	  }
          else{
            pfNoPileUp->Add(pf);
	  }
	}
        else{
          if(vertexFound && closestVtx != pv){
            pfPileUp->Add(pf);
	  }
          else{
            pfNoPileUp->Add(pf);
	  }
        }
      }
    }
    else{
      pfNoPileUp->Add(pf);
    }
  }
}

void 
HttNtuplerMod::FillGenH() 
{
  assert(fParticles); assert(fMCEvtInfo);
  const MCParticle* boson=0, *dau1=0, *dau2=0, *tau1=0, *tau2=0;
  int id1=0, id2=0;
  for(unsigned int i=0; i<fParticles->GetEntries(); ++i){
    const MCParticle* p = fParticles->At(i);
    if( (p->Status()==3) && ((p->PdgId()==25) || (p->PdgId()==35) || (p->PdgId()==36)) ) {
      boson = p;
      // loop through daughters and look for taus
      for(unsigned int ii=0; ii<boson->NDaughters(); ii++){
        const MCParticle* tau = boson->Daughter(ii); 
	if(abs(tau->PdgId())==15){  
	  while(tau->HasDaughter(tau->PdgId(), true) && (tau->Status()!=1)){
	    tau = tau->FindDaughter(tau->PdgId(), true);	  
	  }
	  if(tau->PdgId()>0){ tau1 = tau; }
	  if(tau->PdgId()<0){ tau2 = tau; }
          // loop daughters and look for leptons  
          for(unsigned j=0; j<tau->NDaughters(); j++){
            const MCParticle* d = tau->Daughter(j);
	    if(d->PdgId()==22) continue; // skip photons
            if((abs(d->PdgId())==11) || (abs(d->PdgId())==13)){
	      // traverse down daughter muon tree
              while(d->HasDaughter(d->PdgId(), true) && (d->Status()!=1)){
                d = d->FindDaughter(d->PdgId(),kTRUE);	  
	      }
              if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	      if(d->PdgId()==-11) id2 = -EGenType::kTauElectron;
	      if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	      if(d->PdgId()==-13) id2 = -EGenType::kTauMuon;
	      if(d->PdgId()>0) dau1 = d;
              if(d->PdgId()<0) dau2 = d;
            }
	    if(abs(d->PdgId())>16 || abs(d->PdgId())<7) { // hadronic tau
	      if(tau->PdgId()>0) {
		dau1 = d;
		id1 = EGenType::kTauHadr;
	      } 
	      else {
		dau2 = d;
		id2 = EGenType::kTauHadr;
	      }
	    }
          }
	}	
      }      		      
    }                            
  }
  assert(boson);

  if(!dau1 || !dau2) {
    for(unsigned int j=0; j<boson->NDaughters(); ++j){
      const MCParticle* d = boson->Daughter(j);
      if(d->PdgId()==boson->PdgId()) continue;
      if(d->PdgId()>0) dau1 = d;
      if(d->PdgId()<0) dau2 = d;
    }
    cout << "HttNtuplerMod::FillGenH >> Weird error: Boson ID=" 
	 << boson->PdgId() << " --> Dauther1 ID=" << dau1->PdgId() 
	 << " + Daughter2 ID=" << dau2->PdgId() << endl;
  }
  assert(dau1); assert(dau2);
  FourVectorM vDilep = dau1->Mom() + dau2->Mom();
  
  fGenInfo.pid_1  = fMCEvtInfo->Id1();
  fGenInfo.pid_2  = fMCEvtInfo->Id2();
  fGenInfo.x_1    = fMCEvtInfo->X1();
  fGenInfo.x_2    = fMCEvtInfo->X2();
  fGenInfo.weight = fMCEvtInfo->Weight();
  fGenInfo.vmass  = boson->Mass();
  fGenInfo.vpt    = boson->Pt();
  fGenInfo.vy     = boson->Rapidity();
  fGenInfo.vphi   = boson->Phi();
  fGenInfo.mass   = vDilep.M();
  fGenInfo.pt     = vDilep.Pt(); 
  fGenInfo.y      = vDilep.Rapidity(); 
  fGenInfo.phi    = vDilep.Phi(); 
  fGenInfo.id     = EGenType::kHiggs;
  fGenInfo.pt_1   = dau1->Pt(); 
  fGenInfo.eta_1  = dau1->Eta(); 
  fGenInfo.phi_1  = dau1->Phi();
  fGenInfo.id_1   = id1;
  fGenInfo.pt_2   = dau2->Pt();
  fGenInfo.eta_2  = dau2->Eta(); 
  fGenInfo.phi_2  = dau2->Phi(); 
  fGenInfo.id_2   = id2;
  fGenInfo.decx   = boson->DecayVertex().X();
  fGenInfo.decy   = boson->DecayVertex().Y(); 
  fGenInfo.decz   = boson->DecayVertex().Z();
}

void 
HttNtuplerMod::FillGenZ() 
{
  assert(fParticles); assert(fMCEvtInfo);
  const MCParticle* boson=0, *dau1=0, *dau2=0, *tau1=0, *tau2=0;
  int id1=0, id2=0;
  for(unsigned int i=0; i<fParticles->GetEntries(); ++i){
    const MCParticle* p = fParticles->At(i);
    if( (p->PdgId() == 23) && ((p->Status() == 3 && fIsData!=2 ) || (p->Status() == 2 && fIsData==2)) ){
      boson = p;
    }
    // loop daughters and look for leptons  
    for(unsigned int j=0; j<boson->NDaughters(); ++j){
      const MCParticle* tau = boson->Daughter(j);
      if(abs(tau->PdgId())==15){
        while(tau->HasDaughter(tau->PdgId(), true) && (tau->Status()!=1)){
          tau = tau->FindDaughter(tau->PdgId(), true);
	}
	if(tau->PdgId()>0) tau1 = tau;
	if(tau->PdgId()<0) tau2 = tau;
      }
      const MCParticle* d = boson->Daughter(j);
      if((abs(d->PdgId())==11) || (abs(d->PdgId())==13) || (abs(d->PdgId())==15)){
        // traverse down daughter lepton tree
        while(d->HasDaughter(d->PdgId(), true) && (d->Status()!=1)){
          d = d->FindDaughter(d->PdgId(), true);	  
	}
        if(d->PdgId()== 11) id1 =  EGenType::kElectron;
        if(d->PdgId()==-11) id2 = -EGenType::kElectron;
        if(d->PdgId()== 13) id1 =  EGenType::kMuon;
        if(d->PdgId()==-13) id2 = -EGenType::kMuon;
        if(d->PdgId()== 15) id1 =  EGenType::kTau;
        if(d->PdgId()==-15) id2 = -EGenType::kTau;
        if(abs(d->PdgId())==15){
          if(d->HasDaughter(11)){ 
	    d = d->FindDaughter(11); 
	    if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	    if(d->PdgId()==-11) id2 = -EGenType::kTauElectron;
          }
          if(d->HasDaughter(13)) { 
	    d = d->FindDaughter(13); 
	    if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	    if(d->PdgId()==-13) id2 = -EGenType::kTauMuon;
          }
        }
        if(d->PdgId()>0) dau1 = d;
        if(d->PdgId()<0) dau2 = d;
      }
    }
  }                           
  assert(boson); assert(dau1); assert(dau2);
  FourVectorM vDilep = dau1->Mom() + dau2->Mom();
  if(fMCEvtInfo){
    fGenInfo.pid_1  = fMCEvtInfo->Id1();
    fGenInfo.pid_2  = fMCEvtInfo->Id2();
    fGenInfo.x_1    = fMCEvtInfo->X1();
    fGenInfo.x_2    = fMCEvtInfo->X2();
    fGenInfo.weight = fMCEvtInfo->Weight();
  }
  fGenInfo.vmass  = boson->Mass();
  fGenInfo.vpt    = boson->Pt();
  fGenInfo.vy     = boson->Rapidity();
  fGenInfo.vphi   = boson->Phi();
  fGenInfo.mass   = vDilep.M();
  fGenInfo.pt     = vDilep.Pt(); 
  fGenInfo.y      = vDilep.Rapidity(); 
  fGenInfo.phi    = vDilep.Phi(); 
  fGenInfo.id     = EGenType::kZ;
  fGenInfo.pt_1   = dau1->Pt(); 
  fGenInfo.eta_1  = dau1->Eta(); 
  fGenInfo.phi_1  = dau1->Phi();
  fGenInfo.id_1   = id1;
  fGenInfo.pt_2   = dau2->Pt();
  fGenInfo.eta_2  = dau2->Eta(); 
  fGenInfo.phi_2  = dau2->Phi(); 
  fGenInfo.id_2   = id2;
  fGenInfo.decx   = boson->DecayVertex().X();
  fGenInfo.decy   = boson->DecayVertex().Y(); 
  fGenInfo.decz   = boson->DecayVertex().Z();
}

void 
HttNtuplerMod::FillGenW() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  const MCParticle* boson=0, *dau1=0, *dau2=0;
  int id1=0, id2=0;
  for(unsigned int i=0; i<fParticles->GetEntries(); ++i){
    const MCParticle* p = fParticles->At(i);
    if( (abs(p->PdgId()) == 24) && (p->Status() == 3) ){
      boson = p;
      // loop daughters and look for leptons  
      for(unsigned int j=0; j<boson->NDaughters(); j++){
        const MCParticle* d = boson->Daughter(j);
        if(d->PdgId() == boson->PdgId()) continue;  
	// traverse down daughter lepton tree
        while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1)){
          d = d->FindDaughter(d->PdgId(),kTRUE);	  
	}
        if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	if(d->PdgId()== 15) id1 =  EGenType::kTau;
	if(d->PdgId()==-15) id2 = -EGenType::kTau;
	
	if(abs(d->PdgId())==15){
	  if(d->HasDaughter(11)){ 
	    d = d->FindDaughter(11);
	    if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	    if(d->PdgId()==-11) id2 = -EGenType::kTauElectron; 
	  }
	  if(d->HasDaughter(13)){ 
	    d = d->FindDaughter(13);
	    if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	    if(d->PdgId()==-13) id2 = -EGenType::kTauMuon; 
	  }
	}
	if(d->PdgId()>0) dau1 = d;
        if(d->PdgId()<0) dau2 = d; 
      }		      
    }                            
  }
  assert(boson);
  if(!dau1 || !dau2) {
    while(boson->HasDaughter(boson->PdgId(),kTRUE)){
      boson = boson->FindDaughter(boson->PdgId(),kTRUE);
    }
    // loop daughters and look for leptons  
    for(unsigned int j=0; j<boson->NDaughters(); j++) {
      const MCParticle* d = boson->Daughter(j);
      if(d->PdgId() == boson->PdgId()) continue;  
      // traverse down daughter lepton tree
      while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1)){
    	d = d->FindDaughter(d->PdgId(),kTRUE);  	
      }
      if(d->PdgId()== 11) id1 =  EGenType::kElectron;
      if(d->PdgId()==-11) id2 = -EGenType::kElectron;
      if(d->PdgId()== 13) id1 =  EGenType::kMuon;
      if(d->PdgId()==-13) id2 = -EGenType::kMuon;
      if(d->PdgId()== 15) id1 =  EGenType::kTau;
      if(d->PdgId()==-15) id2 = -EGenType::kTau;
      
      if(abs(d->PdgId())==15){
        if(d->HasDaughter(11)){ 
          d = d->FindDaughter(11);
          if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
          if(d->PdgId()==-11) id2 = -EGenType::kTauElectron; 
        }
        if(d->HasDaughter(13)){ 
          d = d->FindDaughter(13);
          if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
          if(d->PdgId()==-13) id2 = -EGenType::kTauMuon; 
        }
      }
      if(d->PdgId()>0) dau1 = d;
      if(d->PdgId()<0) dau2 = d; 
    }
  }
  assert(dau1);
  assert(dau2);
  
  FourVectorM vDilep = dau1->Mom() + dau2->Mom();
  fGenInfo.pid_1  = fMCEvtInfo->Id1();
  fGenInfo.pid_2  = fMCEvtInfo->Id2();
  fGenInfo.x_1    = fMCEvtInfo->X1();
  fGenInfo.x_2    = fMCEvtInfo->X2();
  fGenInfo.weight = fMCEvtInfo->Weight();
  fGenInfo.vmass  = boson->Mass();
  fGenInfo.vpt    = boson->Pt();
  fGenInfo.vy     = boson->Rapidity();
  fGenInfo.vphi   = boson->Phi();
  fGenInfo.mass   = vDilep.M();
  fGenInfo.pt     = vDilep.Pt(); 
  fGenInfo.y      = vDilep.Rapidity(); 
  fGenInfo.phi    = vDilep.Phi(); 
  fGenInfo.id     = EGenType::kW;
  fGenInfo.pt_1   = dau1->Pt(); 
  fGenInfo.eta_1  = dau1->Eta(); 
  fGenInfo.phi_1  = dau1->Phi();
  fGenInfo.id_1   = id1;
  fGenInfo.pt_2   = dau2->Pt();
  fGenInfo.eta_2  = dau2->Eta(); 
  fGenInfo.phi_2  = dau2->Phi(); 
  fGenInfo.id_2   = id2;
  fGenInfo.decx   = boson->DecayVertex().X();
  fGenInfo.decy   = boson->DecayVertex().Y(); 
  fGenInfo.decz   = boson->DecayVertex().Z();
}

void 
HttNtuplerMod::FillGenWW() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  const MCParticle* boson=0, *dau1=0, *dau2=0;
  int id=0, id1=0, id2=0;
  int nW=0;
  for(unsigned int i=0; i<fParticles->GetEntries(); ++i){
    const MCParticle* p = fParticles->At(i);
    //--------------- PYTHIA FSR mode ---------------//
    // a "branching" in the process tree is created
    // for every physical process; need to scan down
    // the lepton branches and search for photons
    //-----------------------------------------------//
    if( (p->PdgId() == 23) && (p->Status() == 3) ) boson = p;
    if( (abs(p->PdgId()) == 24) && (p->Status() == 3) ){
      boson = p;
      ++nW;
      if(nW==2) id = EGenType::kWW;
      // loop daughters and look for leptons  
      for(unsigned int j=0; j<boson->NDaughters(); j++){
        const MCParticle* d = boson->Daughter(j);
        if(d->PdgId() == boson->PdgId()) continue;  
        
	// traverse down daughter lepton tree
        while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1)){
          d = d->FindDaughter(d->PdgId(),kTRUE);	  
	}
        // ignore neutrinos
	if(abs(d->PdgId())==12) continue;
	if(abs(d->PdgId())==14) continue;
	if(abs(d->PdgId())==16) continue;
	
	if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	if(d->PdgId()== 15) id1 =  EGenType::kTau;
	if(d->PdgId()==-15) id2 = -EGenType::kTau;
	
	if(abs(d->PdgId())==15){
	  if(d->HasDaughter(11)){ 
	    d = d->FindDaughter(11); 
	    if(d->PdgId()== 11) id1 =  EGenType::kTauElectron;
	    if(d->PdgId()==-11) id2 = -EGenType::kTauElectron; 
	  }
	  if(d->HasDaughter(13)){ 
	    d = d->FindDaughter(13); 
	    if(d->PdgId()== 13) id1 =  EGenType::kTauMuon;
	    if(d->PdgId()==-13) id2 = -EGenType::kTauMuon;
	  }
	}
	if(d->PdgId()>0) dau1 = d;
        if(d->PdgId()<0) dau2 = d; 
      }		      
    }                            
  }
  
  FourVectorM vDilep;
  if(dau1 && dau2){
    vDilep = dau1->Mom() + dau2->Mom();
  }
  fGenInfo.pid_1  = fMCEvtInfo->Id1();
  fGenInfo.pid_2  = fMCEvtInfo->Id2();
  fGenInfo.x_1    = fMCEvtInfo->X1();
  fGenInfo.x_2    = fMCEvtInfo->X2();
  fGenInfo.weight = fMCEvtInfo->Weight();
  fGenInfo.vmass  = 0;
  fGenInfo.vpt    = 0;
  fGenInfo.vy     = 0;
  fGenInfo.vphi   = 0;
  fGenInfo.mass   = (dau1 && dau2) ? vDilep.M()        : 0;
  fGenInfo.pt     = (dau1 && dau2) ? vDilep.Pt()       : 0; 
  fGenInfo.y      = (dau1 && dau2) ? vDilep.Rapidity() : 0; 
  fGenInfo.phi    = (dau1 && dau2) ? vDilep.Phi()      : 0; 
  fGenInfo.id     = id;
  fGenInfo.pt_1   = dau1 ? dau1->Pt()  : 0; 
  fGenInfo.eta_1  = dau1 ? dau1->Eta() : 0; 
  fGenInfo.phi_1  = dau1 ? dau1->Phi() : 0;
  fGenInfo.id_1   = id1;
  fGenInfo.pt_2   = dau2 ? dau2->Pt()  : 0;
  fGenInfo.eta_2  = dau2 ? dau2->Eta() : 0; 
  fGenInfo.phi_2  = dau2 ? dau2->Phi() : 0; 
  fGenInfo.id_2   = id2;
  fGenInfo.decx   = boson ? boson->DecayVertex().X() : -999;
  fGenInfo.decy   = boson ? boson->DecayVertex().Y() : -999; 
  fGenInfo.decz   = boson ? boson->DecayVertex().Z() : -999;
}

void 
HttNtuplerMod::FillSVfit(TClonesArray* iArr, Particle* iPart0, Particle* iPart1, TMatrixD iMatrix) 
{
  TClonesArray& rSVfitArr = *iArr;
  const int index = rSVfitArr.GetEntries();
  new(rSVfitArr[index]) TSVfit();
  TSVfit* pSVfit = (TSVfit*)rSVfitArr[index];
  pSVfit->daughter1 = iPart0->Mom();
  pSVfit->daughter2 = iPart1->Mom();
  pSVfit->cov_00    = iMatrix(0,0);
  pSVfit->cov_10    = iMatrix(1,0);
  pSVfit->cov_01    = iMatrix(0,1);
  pSVfit->cov_11    = iMatrix(1,1);
}

bool 
HttNtuplerMod::looseEleId(const Electron* iElectron) 
{ 
  if(iElectron->Pt()  < fEleEtMin)  return false;
  if(iElectron->Pt()  > fEleEtMax)  return false;
  if(iElectron->Eta() < fEleEtaMin) return false;
  if(iElectron->Eta() > fEleEtaMax) return false;
  if(!fEleTools->PassSpikeRemovalFilter(iElectron)) return false;
  if(IsConversion(iElectron))       return false;
  if(fIsData==2 && !iElectron->BestTrk())         return false;
  if(iElectron->BestTrk()->NExpectedHitsInner() > 0) return false;
  return true;
}

bool 
HttNtuplerMod::looseMuId(const Muon* iMu) 
{ 
  if(iMu->TrackerTrk() == 0)                  return false;
  if(iMu->BestTrk()->Pt()       < fMuPtMin)   return false;
  if(iMu->BestTrk()->Eta()      < fMuEtaMin)  return false;
  if(iMu->BestTrk()->Eta()      > fMuEtaMax)  return false;
  if(iMu->BestTrk()->Pt()       > fMuPtMax)   return false;
  if(iMu->TrackerTrk()->PtErr()/iMu->Pt()     > 0.1)        return false;
  if(iMu->TrackerTrk()->NPixelHits()          < 1  )        return false;
  if(iMu->TrackerTrk()->NHits() < 11 )        return false;
  if(iMu->NValidHits()          < 1  )        return false;
  if(iMu->NMatches()            < 1  )        return false;
  return true;
}

