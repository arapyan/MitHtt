#include "cstdlib"

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
//#include "MitAna/DataTree/interface/DCASig.h"

#include "MitHtt/Ntupler/interface/HttNtupler.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/IsolationTools.h"

#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"

using namespace mithep;

ClassImp(mithep::BaconNtupler)

HttNtupler::HttNtupler(const char *name, const char *title):
  BaseMod        (name,title),
  fOutputFile    ( 0),
  fEventTree     ( 0),
  fOutputName    ("ntuple.root"),
  fBeamSpotName  (Names::gkBeamSpotBrn),
  fPrimVtxName   (Names::gkPVBrn),
  fTrigMaskName  (Names::gkHltBitBrn),
  fPileupName    (Names::gkPileupInfoBrn),
  fPUEnergyDensityName(Names::gkPileupEnergyDensityBrn),
  fPartName      (Names::gkMCPartBrn),
  fMCEvtInfoName (Names::gkMCEvtInfoBrn),
  fMuonName      (Names::gkMuonBrn),
  fElectronName  (Names::gkElectronBrn),
  fHPSTauName    ("HPSTaus"),
  fPFJetName     (Names::gkPFJetBrn),
  fPFMetName     ("PFMet"),
  fPhotonName    (Names::gkPhotonBrn),
  fConversionName(Names::gkMvfConversionBrn),
  fPFCandidateName(Names::gkPFCandidatesBrn),
  fEmbedWeightName("EmbedWeight"),
  //fDCASigName     ("DCASig"),
  fParticles      ( 0),
  fMCEvtInfo      ( 0),
  fGenJets        ( 0),
  fBeamSpot       ( 0),
  fPrimVerts      ( 0),
  fTrigMask       ( 0),
  fPileup         ( 0),
  fPUEnergyDensity( 0),
  fMuons          ( 0),
  fElectrons      ( 0),
  fPFJets         ( 0),
  fPFMet          ( 0),
  fPhotons        ( 0),  
  fConversions    ( 0),
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
  fPFTauPtMin     (1.), 
  fPFTauEtaMax    (6.),
  fJetPtMin       (15),
  fPhotonEtMin    (10),
  fMinNTracksFit  ( 0),
  fMinNdof        ( 4),
  fMaxAbsZ        (24),
  fMaxRho         ( 2),
  fJetCorrector   ( 0),
  fJetUncertainties( 0)
  //BaseMod(name,title), fOutputName("HttNtuple.root"), fBeamSpotName(Names::gkBeamSpotBrn), fPrimVtxName(Names::gkPVBrn), fTrigMaskName(Names::gkHltBitBrn), 
  //fPileupName(Names::gkPileupInfoBrn), fPUEnergyDensityName(Names::gkPileupEnergyDensityBrn),fMCEvtInfoName(Names::gkMCEvtInfoBrn), fPartName(Names::gkMCPartBrn), 
  //fMuonName(Names::gkMuonBrn), fElectronName(Names::gkElectronBrn), fHPSTauName("HPSTaus"), fPFJetName(Names::gkPFJetBrn), fPFMetName("PFMet"), fPhotonName(Names::gkPhotonBrn), 
  //fConversionName(Names::gkMvfConversionBrn), fPFCandidateName(Names::gkPFCandidatesBrn), fEmbedWeightName("EmbedWeight"), fParticles(0), fMCEvtInfo(0), fGenJets(0), 
  //fBeamSpot(0), fPrimVerts(0), fTrigMask(0), fPileup(0), fPUEnergyDensity(0), fMuons(0), fElectrons(0), fPFTaus(0), fPFJets(0), fPFMet(0), fPhotons(0), fConversions(0), 
  //fPFCandidates(0), fIsData(0), fUseGen(0), fPrintTable(false), fSkipIfHLTFail(false), fMuPtMin(15.), fMuPtMax(9999.), fMuEtaMin(-3.), fMuEtaMax(3.), fEleEtMin(15.),
  //fEleEtMax(9999.), fEleEtaMin(-3.), fEleEtaMax(3.), fPFTauPtMin(1.), fPFTauEtaMax(6.), fJetPtMin(15.), fPhotonEtMin(10.), fMinNTracksFit(0), fMinNdof(4), fMaxAbsZ(24.),
  //fMaxRho(2.), fDzMax(0.2)
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
  //TMet::Class()    ->IgnoreTObjectStreamer();
  TPFTau::Class()    ->IgnoreTObjectStreamer();
}

HttNtupler::~HttNtupler()
{
}	

void 
HttNtupler::Begin()
{
}

void 
HttNtupler::SlaveBegin()
{
  fHLTTool = new HLTTool();
  setupInput(); setupOutput();
  // setup jet corrections for particle flow jets
  std::vector<JetCorrectorParameters> correctionParameters;
  for(unsigned int icorr=0; icorr<fJetCorrParsv.size(); icorr++){ correctionParameters.push_back(JetCorrectorParameters(fJetCorrParsv[icorr].Data())); }
  fJetCorrector = new FactorizedJetCorrector(correctionParameters); 
  // setup jet energy scale uncertainties
  JetCorrectorParameters param(std::string(TString::Format("%s/src/MitPhysics/data/START42_V12_AK5PF_Uncertainty.txt", getenv("CMSSW_BASE"))));
  fJetUncertainties = new JetCorrectionUncertainty(param);
  // initialize tools for electron ID
  fEleTools = new ElectronTools();
  // initialize tools for muon ID
  fMuonTools = new MuonTools();

  metSign = new MetSignificance();
  
  // setup selection with JSON file, if necessary
  for(unsigned int idx=0; idx<fJSONv.size(); ++idx){
    frlrm.AddJSONFile(fJSONv[idx].Data());
  }
  fLastRunLumi = RunLumiRangeMap::RunLumiPairType(0,0);
}

void 
HttNtupler::BeginRun()
{
  if(HasHLTInfo() && fPrintTable) { GetHLTTable()->Print(); }
}

void 
HttNtupler::EndRun()
{
}

void 
HttNtupler::SlaveTerminate()
{
  fEventTree ->Print(); fOutputFile->Write(); fOutputFile->Close(); cleanup();
  delete fJetCorrector; delete fJetUncertainties; delete fEleTools, delete fMuonTools, delete metSign;

  // dump json file
  TString jsonfname = fOutputName.ReplaceAll("root","json");
  if( fIsData ){
    fRunLumiSet.DumpJSONFile(jsonfname.Data());
  }
}

void 
HttNtupler::Terminate()
{
}

void 
HttNtupler::Process()
{
  // check for run and lumi ranges
  RunLumiRangeMap::RunLumiPairType rl(GetEventHeader()->RunNum(), GetEventHeader()->LumiSec());
  if( fJSONv.size()>0 && !frlrm.HasRunLumi(rl) ){ return; } // not certified run? Skip to next event...
  if( rl!=fLastRunLumi ){ fLastRunLumi = rl; fRunLumiSet.Add(rl); }

  // load relevant Bambu branches
  loadBambuBranches();
  // load trigger table
  loadTriggerTable(fTriggerBits);
  // check whether a trigger table is available or not (if required)
  if( fSkipIfHLTFail && fTriggerBits==0 ){ return; }
  // increment the number of events that have been processed
  IncNEventsProcessed();
  // reset object arrays befor filling
  resetOutputArrays();
  // fill generator information
  if( fUseGen>0 ){ fillGenerator(); }
  // fill general event info
  fillCommon(fTriggerBits);
  // loop and fill muons
  fMuonFiller->fillMuons(fMuons,fVertex,fPFNoPU,fPFPU,fPFCandidates);
  //fillMuons();
  // loop and fill electrons
  fillElecs();
  // loop and fill taus
  fillPFTaus();
  // fill information needed for svfit
  fillSVfit();
  // loop and fill jets
  fillJets();
  // loop and fill photons
  //fillPhotons();
  // fill the tree
  fEventTree->Fill();
}

void
HttNtupler::loadTriggerTable(TriggerBits& trigBits) 
{ 
  trigBits=0;
  if( HasHLTInfo() ){
    const TriggerTable* hltTable = GetHLTTable();
    for(unsigned int itrig=0; itrig<fTriggerNamesv.size(); ++itrig){
      const TriggerName* trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if( !trigname ) continue;
      if( fTrigMask->At(trigname->Id()) ) { trigBits[fTriggerIdsv[itrig]] = true; }
    }  
  }
}

void 
HttNtupler::loadBambuBranches()
{
  LoadBranch( fMuonName            );
  LoadBranch( fElectronName        );
  LoadBranch( fHPSTauName          );
  LoadBranch( fPrimVtxName         );
  LoadBranch( fBeamSpotName        );
  LoadBranch( fPFJetName           );
  LoadBranch( fTrigMaskName        );
  LoadBranch( fPFMetName           ); 
  LoadBranch( fPhotonName          );
  LoadBranch( fConversionName      );
  LoadBranch( fPUEnergyDensityName );
  LoadBranch( fPFCandidateName     );
  //LoadBranch( fDCASigName);
  if( !fIsData ){
    LoadBranch( fPartName          );
    if( fUseGen==ESampleType::kEmbed ){ 
      LoadBranch( fEmbedWeightName   );
    } else {
      LoadBranch( fMCEvtInfoName     );
      LoadBranch( fPileupName        );
      LoadBranch( Names::gkGenJetBrn );
    }
  }
}

void 
HttNtupler::resetOutputArrays()
{ 
  fMuonArr       ->Clear();
  fElectronArr   ->Clear();
  fHPSTauArr     ->Clear();
  fPFJetArr      ->Clear();
  //fMetArr      ->Clear();
  fPhotonArr     ->Clear();
  fPVArr         ->Clear();
  fSVfitEMuArr   ->Clear();
  fSVfitETauArr  ->Clear();
  fSVfitMuTauArr ->Clear();
}

void 
HttNtupler::fillGenerator() 
{
  assert(fParticles);
  const MCParticle* boson_a=0, *dau1_a=0, *dau2_a=0; 
  const MCParticle* boson_b=0, *dau1_b=0, *dau2_b=0;
  for(unsigned int i=0; i<fParticles->GetEntries(); ++i){
    const MCParticle* p = fParticles->At(i);
    if((p->Status() != 3 && fUseGen != ESampleType::kEmbed) || (fUseGen == ESampleType::kEmbed && p->Status() != 2)) continue;
    // fill boson 
    int id = p->PdgId();
    if((fUseGen==ESampleType::kH)      && (id==25 || id==35 || id==36)) boson_a = p;
    if((fUseGen==ESampleType::kEmbed)  &&  id== 23)                     boson_a = p;
    if((fUseGen==ESampleType::kZ)      &&  id== 23)                     boson_a = p;
    if((fUseGen==ESampleType::kVV) && (id==23 || abs(id)==24)) { if(!boson_a) { boson_a = p; } else if(!boson_b) boson_b = p; }
  }
  // special treatment for powheg-herwig interface...
  if(!boson_a ) {
    assert(boson_b);
    while(boson_a->NDaughters()==1) boson_a = boson_a->FindDaughter(boson_a->PdgId());
    while(boson_b->NDaughters()==1) boson_b = boson_b->FindDaughter(boson_b->PdgId());
  }
  assert(boson_a);
  // madgraph zz sample has a few one-z events
  if(fUseGen==ESampleType::kVV && !boson_b && boson_a->NDaughters()==5) boson_b = boson_a;
  if(fUseGen==ESampleType::kVV) assert(boson_b);
  // assign daughters for first boson
  if(boson_a) fillMCParticles(boson_a,dau1_a,dau2_a);
  // assign daughters for second boson
  else if(boson_b) fillMCParticles(boson_b,dau1_b,dau2_b);
  // madgraph zz sample...
  if(fUseGen==ESampleType::kVV && boson_a->NDaughters()==5) {
    dau1_a = dau2_a = dau1_b = dau2_b = 0;
    for(unsigned int idau=0; idau<boson_a->NDaughters(); idau++) { 
      // find the two leptons that aren't already assigned to boson_a
      const MCParticle* daughter = boson_a->Daughter(idau);
      if(daughter->PdgId()==boson_a->PdgId())  continue;
      if(daughter->PdgId()==22) continue; // skip photons
      // this assumes that leptons of same flavor are adjacent in the list of daughters...
      int id = daughter->PdgId();
      if(id > 0) { (!dau1_a) ? dau1_a = daughter : dau1_b = daughter;}
      if(id < 0) { (!dau2_a) ? dau2_a = daughter : dau2_b = daughter;}
    }
  }
  // do the filling
  if(fMCEvtInfo){
    // fill general event info
    fGenInfo.pid_1  = fMCEvtInfo->Id1();
    fGenInfo.pid_2  = fMCEvtInfo->Id2();
    fGenInfo.x_1    = fMCEvtInfo->X1();
    fGenInfo.x_2    = fMCEvtInfo->X2();
    fGenInfo.weight = fMCEvtInfo->Weight();
  }
  fGenInfo.id_a     = pdgId(boson_a);
  fGenInfo.vmass_a  = boson_a->Mass();
  fGenInfo.vpt_a    = boson_a->Pt();
  fGenInfo.vy_a     = boson_a->Rapidity();
  fGenInfo.vphi_a   = boson_a->Phi();
  if(dau1_a != 0) { 
    fGenInfo.pt_1_a   = dau1_a->Pt(); 
    fGenInfo.eta_1_a  = dau1_a->Eta(); 
    fGenInfo.phi_1_a  = dau1_a->Phi();
    fGenInfo.id_1_a   = pdgId(dau1_a);
  }
  if(dau2_a != 0) { 
    fGenInfo.pt_2_a   = dau2_a->Pt();
    fGenInfo.eta_2_a  = dau2_a->Eta(); 
    fGenInfo.phi_2_a  = dau2_a->Phi(); 
    fGenInfo.id_2_a   = pdgId(dau2_a);
  }
  fGenInfo.id_b     = boson_b ? pdgId(boson_b)      : 0;
  fGenInfo.vmass_b  = boson_b ? boson_b->Mass()     : 0;
  fGenInfo.vpt_b    = boson_b ? boson_b->Pt()       : 0;
  fGenInfo.vy_b     = boson_b ? boson_b->Rapidity() : 0;
  fGenInfo.vphi_b   = boson_b ? boson_b->Phi()      : 0;
  
  if(dau1_a != 0) { 
    FourVectorM lL1;  int lId1 = 0;
 
    lL1  = visibleMCMomentum(dau1_a); 
    lId1 = pdgId(dau1_a);
    
    fGenInfo.pt_1_b   = lL1.Pt(); 
    fGenInfo.eta_1_b  = lL1.Eta(); 
    fGenInfo.phi_1_b  = lL1.Phi();
    fGenInfo.id_1_b   = lId1;
  } 
  if(dau2_a != 0) { 
    FourVectorM lL2; int lId2 = 0;
    lL2  = visibleMCMomentum(dau2_a);
    lId2 = pdgId(dau2_a);
    
    fGenInfo.pt_2_b   = lL2.Pt();
    fGenInfo.eta_2_b  = lL2.Eta(); 
    fGenInfo.phi_2_b  = lL2.Phi(); 
    fGenInfo.id_2_b   = lId2;
  }
  fGenInfo.decx   = boson_a->DecayVertex().X();
  fGenInfo.decy   = boson_a->DecayVertex().Y(); 
  fGenInfo.decz   = boson_a->DecayVertex().Z();
}

void 
HttNtupler::fillCommon(TriggerBits& trigBits) 
{
  // get pileup information (for MC samples only)
  int npu0 = -1; int npu1 = -1; int npu2 = -1;
  if( !fIsData && fUseGen!=ESampleType::kEmbed ){
    for(unsigned int idx=0; idx<fPileup->GetEntries(); ++idx){
      if( fPileup->At(idx)->GetBunchCrossing() ==  0 ) npu0 = fPileup->At(idx)->GetPU_NumInteractions();
      if( fPileup->At(idx)->GetBunchCrossing() ==  1 ) npu1 = fPileup->At(idx)->GetPU_NumInteractions();
      if( fPileup->At(idx)->GetBunchCrossing() == -1 ) npu2 = fPileup->At(idx)->GetPU_NumInteractions();
    }
  }  
  // get beamspot information
  double bsx=99999, bsy=99999, bsz=99999;
  if( fBeamSpot ){
    const BeamSpot* bs = fBeamSpot->At(0);
    bsx = bs->X(); bsy = bs->Y(); bsz = bs->Z();
  }
  // get primary vertices
  const Vertex* bestPV=0; bool hasGoodPV=false;  
  for(unsigned int idx=0; idx<fPrimVerts->GetEntries(); ++idx){
    const Vertex* pv = fPrimVerts->At(idx);
    if( pv->NTracksFit()     < fMinNTracksFit ){ continue; }
    if( pv->Ndof()	     < fMinNdof       ){ continue; }
    if( fabs(pv->Z())        > fMaxAbsZ       ){ continue; }
    if( pv->Position().Rho() > fMaxRho        ){ continue; }
    TClonesArray& rPVArr = *fPVArr;
    assert(rPVArr.GetEntries() < rPVArr.GetSize());
    const int index = rPVArr.GetEntries();  
    new(rPVArr[index]) TVertex();
    TVertex* pVertex = (TVertex*)rPVArr[index];
    for(unsigned int itrk=0; itrk<pv->NTracks(); ++itrk){
      pVertex->sumPt+=pv->Trk(itrk)->Pt();				   
    }
    pVertex->nTracksFit = pv->NTracksFit();
    pVertex->ndof       = pv->Ndof();      
    pVertex->chi2       = pv->Chi2();  
    pVertex->x          = pv->X();
    pVertex->y          = pv->Y();
    pVertex->z          = pv->Z();
    pVertex->sumPt=0;
    hasGoodPV = true;
    if(!bestPV){ bestPV = pv; }
  }
  if(!bestPV) bestPV = fPrimVerts->At(0); fVertex = bestPV;
  // calculate track MET
  double trkMetx=0, trkMety=0, trkSumET=0;
  for(unsigned int i=0; i<fPFCandidates->GetEntries(); ++i){
    const double trkDzCut  = 0.1;
    const PFCandidate* pfcand = fPFCandidates->At(i);
    if( (pfcand->HasTrackerTrk() && (fabs(pfcand->TrackerTrk()->DzCorrected(*fVertex))<trkDzCut)) || (pfcand->HasGsfTrk() && (fabs(pfcand->GsfTrk()->DzCorrected(*fVertex))<trkDzCut)) ) {
      trkMetx -= pfcand->Px(); trkMety -= pfcand->Py(); trkSumET += pfcand->Pt();
    }
  }
  TLorentzVector trkMet; trkMet.SetPxPyPzE(trkMetx, trkMety, 0, 0);
  // create pileup and noPileup collections with association using cut along z axis 
  separatePileup(fPFCandidates, fVertex, fPrimVerts, fPFPileUp   , fPFNoPileUp   , true );
  // create pileup and noPileup collections w/o association using cut along z axis 
  separatePileup(fPFCandidates, fVertex, fPrimVerts, fPFPileUpNoZ, fPFNoPileUpNoZ, false);
  
  fEventInfo.runNum       = GetEventHeader()->RunNum();
  fEventInfo.evtNum       = GetEventHeader()->EvtNum();
  fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
  fEventInfo.nPU          = fIsData ? 0 : npu0;
  fEventInfo.nPUPlus      = fIsData ? 0 : npu1;
  fEventInfo.nPUMinus     = fIsData ? 0 : npu2;
  fEventInfo.triggerBits  = trigBits;
  fEventInfo.pvx          = fVertex->X();
  fEventInfo.pvy          = fVertex->Y();
  fEventInfo.pvz          = fVertex->Z();
  fEventInfo.bsx          = bsx;
  fEventInfo.bsy          = bsy;
  fEventInfo.bsz          = bsz;
  fEventInfo.rho          = fPUEnergyDensity->At(0)->Rho();
  fEventInfo.rhoHighEta   = fPUEnergyDensity->At(0)->RhoHighEta();
  fEventInfo.hasGoodPV    = hasGoodPV;
  fEventInfo.pfMET        = fPFMet->At(0)->Et();
  fEventInfo.pfMETphi     = fPFMet->At(0)->Phi();
  fEventInfo.pfSumET      = fPFMet->At(0)->SumEt();
  fEventInfo.trkMET       = trkMet.Pt();
  fEventInfo.trkMETphi    = trkMet.Phi();
  fEventInfo.trkSumET     = trkSumET;
  fEventInfo.embWeight    = (fUseGen==ESampleType::kEmbed) ? fEmbedWeight->At(fEmbedWeight->GetEntries()-1)->Weight() : 1;  
}

void 
HttNtupler::fillMuons()
{
  assert(fMuons); assert(fPFCandidates);
  for(unsigned int i=0; i<fMuons->GetEntries(); ++i){
    const Muon* mu=fMuons->At(i); 
    if(!mu->HasTrk()) continue; 

    double Muon_ChargedIso03=0, Muon_ChargedIso03FromOtherVertices=0, Muon_NeutralIso03_01Threshold=0, Muon_NeutralIso03_05Threshold=0, Muon_NeutralIso03_10Threshold=0, Muon_NeutralIso03_15Threshold=0;
    double Muon_ChargedIso04=0, Muon_ChargedIso04FromOtherVertices=0, Muon_NeutralIso04_01Threshold=0, Muon_NeutralIso04_05Threshold=0, Muon_NeutralIso04_10Threshold=0, Muon_NeutralIso04_15Threshold=0;
    double tmpPFIso_NormalizedMeanEta=0, tmpPFIso_NormalizedMeanPhi=0;
    double tmpChargedIso_DR0p0To0p1=0, tmpChargedIso_DR0p1To0p2=0, tmpChargedIso_DR0p2To0p3=0, tmpChargedIso_DR0p3To0p4=0, tmpChargedIso_DR0p4To0p5=0, tmpChargedIso_DR0p5To0p7=0, tmpChargedIso_DR0p7To1p0=0;
    double tmpChargedIso_MeanEta_numerator=0, tmpChargedIso_MeanPhi_numerator=0, tmpChargedIso_SigmaEtaEta_numerator=0, tmpChargedIso_SigmaPhiPhi_numerator=0, tmpChargedIso_SigmaEtaPhi_numerator=0, tmpChargedIso_Covariance_denominator=0;
    double tmpGammaIso_DR0p0To0p1=0, tmpGammaIso_DR0p1To0p2=0, tmpGammaIso_DR0p2To0p3=0, tmpGammaIso_DR0p3To0p4=0, tmpGammaIso_DR0p4To0p5=0, tmpGammaIso_DR0p5To0p7=0, tmpGammaIso_DR0p7To1p0=0;
    double tmpGammaIso_MeanEta_numerator=0, tmpGammaIso_MeanPhi_numerator=0, tmpGammaIso_SigmaEtaEta_numerator=0, tmpGammaIso_SigmaPhiPhi_numerator=0, tmpGammaIso_SigmaEtaPhi_numerator=0, tmpGammaIso_Covariance_denominator=0;
    double tmpNeutralHadronIso_DR0p0To0p1=0, tmpNeutralHadronIso_DR0p1To0p2=0, tmpNeutralHadronIso_DR0p2To0p3=0, tmpNeutralHadronIso_DR0p3To0p4=0, tmpNeutralHadronIso_DR0p4To0p5=0, tmpNeutralHadronIso_DR0p5To0p7=0, tmpNeutralHadronIso_DR0p7To1p0=0;
    double tmpNeutralHadronIso_MeanEta_numerator=0, tmpNeutralHadronIso_MeanPhi_numerator=0, tmpNeutralHadronIso_SigmaEtaEta_numerator=0, tmpNeutralHadronIso_SigmaPhiPhi_numerator=0, tmpNeutralHadronIso_SigmaEtaPhi_numerator=0, tmpNeutralHadronIso_Covariance_denominator=0;


    double zMuon = 0.0;
    if(mu->BestTrk()) zMuon = mu->BestTrk()->DzCorrected(*fVertex);
    const PFCandidate *pfMuonMatched = 0;
    for (unsigned int p=0; p<fPFCandidates->GetEntries();p++) {   
      const PFCandidate *pf = fPFCandidates->At(p);
      if (!pfMuonMatched) {
        if(pf->TrackerTrk() && mu->TrackerTrk() && pf->TrackerTrk() == mu->TrackerTrk()) {
          pfMuonMatched = pf;
        }           
      }

      //exclude the muon itself
      if(pf->TrackerTrk() && mu->TrackerTrk() && pf->TrackerTrk() == mu->TrackerTrk()) continue;      

      double dr = MathUtils::DeltaR(mu->Mom(), pf->Mom());
 
      if(pf->BestTrk()) {
        double deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*fVertex) - zMuon);
        //remove charged particles from other vertices
        if (deltaZ <= 0.1) {
          if (dr < 0.4) Muon_ChargedIso04 += pf->Pt();
          if (dr < 0.3) Muon_ChargedIso03 += pf->Pt();
        } else {
          if (dr < 0.4) Muon_ChargedIso04FromOtherVertices += pf->Pt();          
          if (dr < 0.3) Muon_ChargedIso03FromOtherVertices += pf->Pt();
        }
      } else {
        if (dr < 0.4) {
          if (pf->Pt() > 0.1) Muon_NeutralIso04_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Muon_NeutralIso04_05Threshold += pf->Pt();
          if (pf->Pt() > 1.0) Muon_NeutralIso04_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Muon_NeutralIso04_15Threshold += pf->Pt();
        }
        if (dr < 0.3) {
          if (pf->Pt() > 0.1) Muon_NeutralIso03_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Muon_NeutralIso03_05Threshold += pf->Pt();
          if (pf->Pt() > 1.0) Muon_NeutralIso03_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Muon_NeutralIso03_15Threshold += pf->Pt();
        }
      }

      // New Isolation Calculations
      double deta = (mu->Eta() - pf->Eta());
      double dphi = MathUtils::DeltaPhi(double(mu->Phi()),double(pf->Phi()));
      if (dr < 1.0) {
        bool isLeptonFootprint = findLeptonFootprint(pf);
        if (!isLeptonFootprint) {
	  bool passVeto = true;
	  //Charged
	  if(pf->BestTrk()) {	   
	    if (!(fabs(pf->BestTrk()->DzCorrected(*fVertex) - zMuon) < 0.2)) passVeto = false;
  	    // Veto any PFmuon, or PFEle
  	    if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = false;
  	    // Footprint Veto
  	    if (dr < 0.015) passVeto = false;
  	    if (passVeto) {
  	      if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += pf->Pt();
              if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += pf->Pt();
  	      if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += pf->Pt();
   	      if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->Pt();
   	      if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->Pt();
  	      if (dr >= 0.5 && dr < 0.7) tmpChargedIso_DR0p5To0p7 += pf->Pt();
  	      if (dr >= 0.7 && dr < 1.0) tmpChargedIso_DR0p7To1p0 += pf->Pt();
  	      tmpChargedIso_Covariance_denominator += pf->Pt();
  	      tmpChargedIso_MeanEta_numerator += pf->Pt() * deta;
  	      tmpChargedIso_MeanPhi_numerator += pf->Pt() * dphi;
  	      tmpChargedIso_SigmaEtaEta_numerator += pf->Pt() * deta*deta;
  	      tmpChargedIso_SigmaPhiPhi_numerator += pf->Pt() * dphi*dphi;
  	      tmpChargedIso_SigmaEtaPhi_numerator += pf->Pt() * deta*dphi;
  	    } //pass veto	  
  	  }
  	  //Gamma
  	  else if (pf->PFType() == PFCandidate::eGamma) {	   
  	    if (dr < 0.05) passVeto = false;
  	    if (passVeto) {
  	      if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->Pt();
              if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->Pt();
  	      if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->Pt();
  	      if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->Pt();
  	      if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->Pt();
  	      if (dr >= 0.5 && dr < 0.7) tmpGammaIso_DR0p5To0p7 += pf->Pt();
  	      if (dr >= 0.7 && dr < 1.0) tmpGammaIso_DR0p7To1p0 += pf->Pt();
  	      tmpGammaIso_Covariance_denominator += pf->Pt();
  	      tmpGammaIso_MeanEta_numerator += pf->Pt() * deta;
  	      tmpGammaIso_MeanPhi_numerator += pf->Pt() * dphi;
  	      tmpGammaIso_SigmaEtaEta_numerator += pf->Pt() * deta*deta;
  	      tmpGammaIso_SigmaPhiPhi_numerator += pf->Pt() * dphi*dphi;
  	      tmpGammaIso_SigmaEtaPhi_numerator += pf->Pt() * deta*dphi;	   	   
  	    }
  	  }
  	  //NeutralHadron
  	  else {
  	    if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += pf->Pt();
            if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += pf->Pt();
  	    if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->Pt();
  	    if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->Pt();
  	    if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->Pt();
  	    if (dr >= 0.5 && dr < 0.7) tmpNeutralHadronIso_DR0p5To0p7 += pf->Pt();
  	    if (dr >= 0.7 && dr < 1.0) tmpNeutralHadronIso_DR0p7To1p0 += pf->Pt();
  	    tmpNeutralHadronIso_Covariance_denominator += pf->Pt();
  	    tmpNeutralHadronIso_MeanEta_numerator += pf->Pt() * deta;
  	    tmpNeutralHadronIso_MeanPhi_numerator += pf->Pt() * dphi;
  	    tmpNeutralHadronIso_SigmaEtaEta_numerator += pf->Pt() * deta*deta;
  	    tmpNeutralHadronIso_SigmaPhiPhi_numerator += pf->Pt() * dphi*dphi;
  	    tmpNeutralHadronIso_SigmaEtaPhi_numerator += pf->Pt() * deta*dphi;
  	  }
  
  	  //For directional iso calculation  
  	  if (dr > 0 && passVeto) {
  	    tmpPFIso_NormalizedMeanEta += pf->Pt() * deta / dr;
  	    tmpPFIso_NormalizedMeanPhi += pf->Pt() * dphi / dr;
  	  }
        } //not lepton footprint
      } //in 1.0 dr cone
    }//loop over pfcands

    TVector2 PFIsoCentroid(tmpPFIso_NormalizedMeanEta,tmpPFIso_NormalizedMeanPhi);
    double tmpPFIsoAngleSqrSum=0, tmpChargedIsoAngleSqrSum=0, tmpGammaIsoAngleSqrSum=0, tmpNeutralHadronIsoAngleSqrSum=0;

    for (unsigned int p=0; p<fPFCandidates->GetEntries();p++) {   
      const PFCandidate *pf = fPFCandidates->At(p);
        
      //exclude the muon itself
      if(pf->TrackerTrk() && mu->TrackerTrk() && pf->TrackerTrk() == mu->TrackerTrk()) continue;      
  
      // New Isolation Calculations
      double deta = (mu->Eta() - pf->Eta());
      double dphi = mithep::MathUtils::DeltaPhi(Double_t(mu->Phi()),Double_t(pf->Phi()));
      double dr = mithep::MathUtils::DeltaR(mu->Phi(),mu->Eta(), pf->Phi(), pf->Eta());
      
      if (dr < 1.0) {
        bool isLeptonFootprint = findLeptonFootprint(pf);
        if (!isLeptonFootprint) {
  	TVector2 *tmpPFAngle;
  	if (dr > 0) {
  	  tmpPFAngle = new TVector2(pf->Pt()*deta/dr, pf->Pt()*dphi/dr);
  	} else {
  	  continue;
  	}
  	assert(tmpPFAngle);
  	if (tmpPFAngle->Mod()*PFIsoCentroid.Mod() == 0) continue;
  	if (fabs((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod()) - 1) < 0.001) continue;
  	
  	bool passVeto = true;
  	
  	//Charged
  	if(pf->BestTrk()) {	  
  	  if (!(fabs(pf->BestTrk()->DzCorrected(*fVertex) - zMuon) < 0.2)) passVeto = false;
  	  // Veto any PFmuon, or PFEle
  	  if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = false;
  	  // Footprint Veto
  	  if (dr < 0.015) passVeto = false;
  	  if (passVeto) tmpChargedIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);
  	}
  	//Gamma
  	else if (pf->PFType() == PFCandidate::eGamma) {	  
  	  if (dr < 0.05) passVeto = false;
  	  if (passVeto) tmpGammaIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);
  	}
  	//NeutralHadron
  	else {
  	  tmpNeutralHadronIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);
  	}
  	if (passVeto) tmpPFIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);
  	if (tmpPFAngle) delete tmpPFAngle;
        } //not lepton footprint  
      } //in 1.0 dr cone
    } //loop over PF candidates

    double tmpChargedIso_MeanEta  = (tmpChargedIso_Covariance_denominator>0) ? tmpChargedIso_MeanEta_numerator/tmpChargedIso_Covariance_denominator : 0;
    double tmpChargedIso_MeanPhi  = (tmpChargedIso_Covariance_denominator>0) ? tmpChargedIso_MeanPhi_numerator/tmpChargedIso_Covariance_denominator : 0;
    double tmpChargedIso_SigmaEtaEta  = (tmpChargedIso_Covariance_denominator>0) ? sqrt(tmpChargedIso_SigmaEtaEta_numerator/tmpChargedIso_Covariance_denominator) : 0;
    double tmpChargedIso_SigmaPhiPhi  = (tmpChargedIso_Covariance_denominator>0) ? sqrt(tmpChargedIso_SigmaPhiPhi_numerator/tmpChargedIso_Covariance_denominator) : 0;
    double tmpChargedIso_CovEtaPhi  = (tmpChargedIso_Covariance_denominator>0) ? (tmpChargedIso_SigmaEtaPhi_numerator/tmpChargedIso_Covariance_denominator) : 0;
    double tmpChargedIso_SigmaEtaPhi = 0;
    if (tmpChargedIso_SigmaEtaEta*tmpChargedIso_SigmaPhiPhi != 0) tmpChargedIso_SigmaEtaPhi = tmpChargedIso_CovEtaPhi / (tmpChargedIso_SigmaEtaEta*tmpChargedIso_SigmaPhiPhi);
    else tmpChargedIso_SigmaEtaPhi = (tmpChargedIso_CovEtaPhi < 0) ? -1 : (tmpChargedIso_CovEtaPhi > 0);

    double tmpGammaIso_MeanEta  = (tmpGammaIso_Covariance_denominator>0) ? tmpGammaIso_MeanEta_numerator/tmpGammaIso_Covariance_denominator : 0;
    double tmpGammaIso_MeanPhi  = (tmpGammaIso_Covariance_denominator>0) ? tmpGammaIso_MeanPhi_numerator/tmpGammaIso_Covariance_denominator : 0;
    double tmpGammaIso_SigmaEtaEta  = (tmpGammaIso_Covariance_denominator>0) ? sqrt(tmpGammaIso_SigmaEtaEta_numerator/tmpGammaIso_Covariance_denominator) : 0;
    double tmpGammaIso_SigmaPhiPhi  = (tmpGammaIso_Covariance_denominator>0) ? sqrt(tmpGammaIso_SigmaPhiPhi_numerator/tmpGammaIso_Covariance_denominator) : 0;
    double tmpGammaIso_CovEtaPhi  = (tmpGammaIso_Covariance_denominator>0) ? (tmpGammaIso_SigmaEtaPhi_numerator/tmpGammaIso_Covariance_denominator) : 0;
    double tmpGammaIso_SigmaEtaPhi = 0;
    if (tmpGammaIso_SigmaEtaEta*tmpGammaIso_SigmaPhiPhi != 0) tmpGammaIso_SigmaEtaPhi = tmpGammaIso_CovEtaPhi / (tmpGammaIso_SigmaEtaEta*tmpGammaIso_SigmaPhiPhi);
    else tmpGammaIso_SigmaEtaPhi = (tmpGammaIso_CovEtaPhi < 0) ? -1 : (tmpGammaIso_CovEtaPhi > 0);

    double tmpNeutralHadronIso_MeanEta  = (tmpNeutralHadronIso_Covariance_denominator>0) ? tmpNeutralHadronIso_MeanEta_numerator/tmpNeutralHadronIso_Covariance_denominator : 0;
    double tmpNeutralHadronIso_MeanPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? tmpNeutralHadronIso_MeanPhi_numerator/tmpNeutralHadronIso_Covariance_denominator : 0;
    double tmpNeutralHadronIso_SigmaEtaEta  = (tmpNeutralHadronIso_Covariance_denominator>0) ? sqrt(tmpNeutralHadronIso_SigmaEtaEta_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
    double tmpNeutralHadronIso_SigmaPhiPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? sqrt(tmpNeutralHadronIso_SigmaPhiPhi_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
    double tmpNeutralHadronIso_CovEtaPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? (tmpNeutralHadronIso_SigmaEtaPhi_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
    double tmpNeutralHadronIso_SigmaEtaPhi = 0;
    if (tmpNeutralHadronIso_SigmaEtaEta*tmpNeutralHadronIso_SigmaPhiPhi != 0) tmpNeutralHadronIso_SigmaEtaPhi = tmpNeutralHadronIso_CovEtaPhi / (tmpNeutralHadronIso_SigmaEtaEta*tmpNeutralHadronIso_SigmaPhiPhi);
    else tmpNeutralHadronIso_SigmaEtaPhi = (tmpNeutralHadronIso_CovEtaPhi < 0) ? -1 : (tmpNeutralHadronIso_CovEtaPhi > 0);

    TClonesArray& rMuonArr = *fMuonArr; 
    assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
    const int index = rMuonArr.GetEntries();  
    new(rMuonArr[index]) TMuon();
    TMuon* pMuon = (TMuon*)rMuonArr[index];

    // use tracker tracks for kinematics when available
    const Track* muTrk=0;
    if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk();    }
    else if(mu->HasGlobalTrk())     { muTrk = mu->GlobalTrk();     }
    else if(mu->HasStandaloneTrk()) { muTrk = mu->StandaloneTrk(); }
    if((muTrk->Eta() < fMuEtaMin) || (muTrk->Eta() > fMuEtaMax)) continue;
    if((muTrk->Pt () > fMuPtMax ) || (muTrk->Pt () < fMuPtMin) ) continue;

    pMuon->pt       = muTrk->Pt();
    pMuon->eta      = muTrk->Eta();
    pMuon->phi      = muTrk->Phi();
    pMuon->ptErr    = muTrk->PtErr();
    pMuon->staPt    = mu->HasStandaloneTrk() ? mu->StandaloneTrk()->Pt()  : -9999.;  
    pMuon->staEta   = mu->HasStandaloneTrk() ? mu->StandaloneTrk()->Eta() : -9999.;
    pMuon->staPhi   = mu->HasStandaloneTrk() ? mu->StandaloneTrk()->Phi() : -9999.;
    pMuon->trkIso03 = mu->IsoR03SumPt();
    pMuon->emIso03  = mu->IsoR03EmEt();
    pMuon->hadIso03 = mu->IsoR03HadEt();
    pMuon->hoIso03  = mu->IsoR03HoEt();
    pMuon->trkIso05 = mu->IsoR05SumPt();
    pMuon->emIso05  = mu->IsoR05EmEt();
    pMuon->hadIso05 = mu->IsoR05HadEt();
    pMuon->hoIso05  = mu->IsoR05HoEt();

    pMuon->chargedIso03 = Muon_ChargedIso03;
    pMuon->chargedIso03FromOtherVertices = Muon_ChargedIso03FromOtherVertices;
    pMuon->neutralIso03_05Threshold = Muon_NeutralIso03_05Threshold ;
    pMuon->neutralIso03_10Threshold = Muon_NeutralIso03_10Threshold ;
    pMuon->chargedIso04 = Muon_ChargedIso04;
    pMuon->chargedIso04FromOtherVertices = Muon_ChargedIso04FromOtherVertices;
    pMuon->neutralIso04_05Threshold = Muon_NeutralIso04_05Threshold ;
    pMuon->neutralIso04_10Threshold = Muon_NeutralIso04_10Threshold ;
    pMuon->chargedIso_DR0p0To0p1 = tmpChargedIso_DR0p0To0p1;
    pMuon->chargedIso_DR0p1To0p2 = tmpChargedIso_DR0p1To0p2;
    pMuon->chargedIso_DR0p2To0p3 = tmpChargedIso_DR0p2To0p3;
    pMuon->chargedIso_DR0p3To0p4 = tmpChargedIso_DR0p3To0p4;
    pMuon->chargedIso_DR0p4To0p5 = tmpChargedIso_DR0p4To0p5;
    pMuon->chargedIso_DR0p5To0p7 = tmpChargedIso_DR0p5To0p7;
    pMuon->chargedIso_DR0p7To1p0 = tmpChargedIso_DR0p7To1p0;
    pMuon->gammaIso_DR0p0To0p1 = tmpGammaIso_DR0p0To0p1;
    pMuon->gammaIso_DR0p1To0p2 = tmpGammaIso_DR0p1To0p2;
    pMuon->gammaIso_DR0p2To0p3 = tmpGammaIso_DR0p2To0p3;
    pMuon->gammaIso_DR0p3To0p4 = tmpGammaIso_DR0p3To0p4;
    pMuon->gammaIso_DR0p4To0p5 = tmpGammaIso_DR0p4To0p5;
    pMuon->gammaIso_DR0p5To0p7 = tmpGammaIso_DR0p5To0p7;
    pMuon->gammaIso_DR0p7To1p0 = tmpGammaIso_DR0p7To1p0;
    pMuon->neutralIso_DR0p0To0p1 = tmpNeutralHadronIso_DR0p0To0p1;
    pMuon->neutralIso_DR0p1To0p2 = tmpNeutralHadronIso_DR0p1To0p2;
    pMuon->neutralIso_DR0p2To0p3 = tmpNeutralHadronIso_DR0p2To0p3;
    pMuon->neutralIso_DR0p3To0p4 = tmpNeutralHadronIso_DR0p3To0p4;
    pMuon->neutralIso_DR0p4To0p5 = tmpNeutralHadronIso_DR0p4To0p5;
    pMuon->neutralIso_DR0p5To0p7 = tmpNeutralHadronIso_DR0p5To0p7;
    pMuon->neutralIso_DR0p7To1p0 = tmpNeutralHadronIso_DR0p7To1p0;
    pMuon->chargedIso_MeanEta = tmpChargedIso_MeanEta;
    pMuon->chargedIso_MeanPhi = tmpChargedIso_MeanPhi;
    pMuon->chargedIso_SigEtaEta = tmpChargedIso_SigmaEtaEta;
    pMuon->chargedIso_SigPhiPhi = tmpChargedIso_SigmaPhiPhi;
    pMuon->chargedIso_SigEtaPhi = tmpChargedIso_SigmaEtaPhi;
    pMuon->gammaIso_MeanEta = tmpGammaIso_MeanEta;
    pMuon->gammaIso_MeanPhi = tmpGammaIso_MeanPhi;
    pMuon->gammaIso_SigEtaEta = tmpGammaIso_SigmaEtaEta;
    pMuon->gammaIso_SigPhiPhi = tmpGammaIso_SigmaPhiPhi;
    pMuon->gammaIso_SigEtaPhi = tmpGammaIso_SigmaEtaPhi;
    pMuon->neutralIso_MeanEta = tmpNeutralHadronIso_MeanEta;
    pMuon->neutralIso_MeanPhi = tmpNeutralHadronIso_MeanPhi;
    pMuon->neutralIso_SigEtaEta = tmpNeutralHadronIso_SigmaEtaEta;
    pMuon->neutralIso_SigPhiPhi = tmpNeutralHadronIso_SigmaPhiPhi;
    pMuon->neutralIso_SigEtaPhi = tmpNeutralHadronIso_SigmaEtaPhi;
    pMuon->directionalPFIso = tmpPFIsoAngleSqrSum;
    pMuon->directionalChargedIso = tmpChargedIsoAngleSqrSum;
    pMuon->directionalGammaIso = tmpGammaIsoAngleSqrSum;
    pMuon->directionalNeutralIso = tmpNeutralHadronIsoAngleSqrSum;

    pMuon->pfIso03  = computePFMuonIso(mu,0.3);
    pMuon->pfIso04  = computePFMuonIso(mu,0.4);  
    pMuon->d0       = muTrk->D0Corrected(*fVertex);
    pMuon->dz       = muTrk->DzCorrected(*fVertex);
    pMuon->d0Sig    = mu->D0PVSignificance();
    pMuon->ip3d     = mu->Ip3dPV();
    pMuon->ip3dSig  = mu->Ip3dPVSignificance();
    pMuon->d0Ub     = mu->D0PVUB();
    pMuon->d0UbSig  = mu->D0PVUBSignificance();
    pMuon->ip3dUb   = mu->Ip3dPVUB();
    pMuon->ip3dUbSig= mu->Ip3dPVUBSignificance();
    pMuon->tkNchi2  = muTrk->RChi2();
    if     (mu->HasGlobalTrk()    ) { pMuon->muNchi2 = mu->GlobalTrk()->RChi2();     }
    else if(mu->HasStandaloneTrk()) { pMuon->muNchi2 = mu->StandaloneTrk()->RChi2(); }
    else if(mu->HasTrackerTrk()   ) { pMuon->muNchi2 = mu->TrackerTrk()->RChi2();    }
    pMuon->q           = muTrk->Charge();
    pMuon->nValidHits  = mu->NValidHits();
    pMuon->qualityBits = mu->Quality().QualityMask().Mask();
    // NOTE: It is possible for a muon to be TK+SA. The muon reco associates a TK with a SA
    // if chamber matches for the TK and hits for the SA share DetIDs see hypernews thread:
    // https://hypernews.cern.ch/HyperNews/CMS/get/csa08-muons/57/2/1/1/1.html
    pMuon->typeBits = 0;
    if(mu->IsGlobalMuon())     { pMuon->typeBits |= kGlobal;     }
    if(mu->IsTrackerMuon())    { pMuon->typeBits |= kTracker;    }
    if(mu->IsStandaloneMuon()) { pMuon->typeBits |= kStandalone; }
    pMuon->nTkHits         = muTrk->NHits();
    pMuon->nPixHits        = muTrk->NPixelHits();
    pMuon->nSeg            = mu->NSegments();
    pMuon->nMatch          = mu->NMatches();
    pMuon->hltMatchBits    = matchHLT(muTrk->Eta(),muTrk->Phi(),muTrk->Pt());
    pMuon->trkID           = mu->HasTrackerTrk() ? mu->TrackerTrk()->GetUniqueID() :  0;

    // additional variables used for MVA id
    pMuon->trkKink         = mu->TrkKink();
    pMuon->globalKink      = mu->GlbKink();
    pMuon->segCompatibility  = fMuonTools->GetSegmentCompatability(mu);
    pMuon->caloCompatibility = fMuonTools->GetCaloCompatability(mu, kTRUE, kTRUE);
    pMuon->hadEnergy         = mu->HadEnergy();
    pMuon->hadS9Energy       = mu->HadS9Energy();
    pMuon->hoEnergy          = mu->HoEnergy();
    pMuon->hoS9Energy        = mu->HoS9Energy();
    pMuon->emEnergy          = mu->EmEnergy();
    pMuon->emS9Energy        = mu->EmS9Energy();

    pMuon->pfIsoCharged    = computeCommonIso(mu, fPFNoPileUp    , 0.0, 0.4, 0.0001,  1);
    pMuon->pfIsoChargedNoZ = computeCommonIso(mu, fPFNoPileUpNoZ , 0.0, 0.4, 0.0001,  1);
    pMuon->pfIsoNeutral    = computeCommonIso(mu, fPFNoPileUp    , 0.5, 0.4,   0.01,  2);
    pMuon->pfIsoNeutralNoZ = computeCommonIso(mu, fPFNoPileUpNoZ , 0.5, 0.4,   0.01,  2);
    pMuon->pfIsoGamma      = computeCommonIso(mu, fPFNoPileUp    , 0.5, 0.4,   0.01,  3);
    pMuon->pfIsoGammaNoZ   = computeCommonIso(mu, fPFNoPileUpNoZ , 0.5, 0.4,   0.01,  3);
    pMuon->puIso           = computeNaiveIso(mu, fPFPileUp    , 0.5, 0.4, 0.01);
    pMuon->puIsoNoZ        = computeNaiveIso(mu, fPFPileUpNoZ , 0.5, 0.4, 0.01);
    pMuon->pfPx = 0; pMuon->pfPy = 0;
    for(unsigned int i=0; i<fPFCandidates->GetEntries(); ++i){    
      const PFCandidate* pfcand = fPFCandidates->At(i);
      if(mu->HasTrackerTrk() && mu->TrackerTrk() == pfcand->TrackerTrk()){
	pMuon->pfPx = pfcand->Px(); pMuon->pfPy = pfcand->Py();
	break;
      }
    }
  }
}

void 
HttNtupler::fillElecs()
{
  assert(fElectrons); assert(fPFCandidates);

  for(unsigned int i=0; i<fElectrons->GetEntries(); ++i){
    const Electron* ele=fElectrons->At(i);  
    if((ele->Pt()  < fEleEtMin)  || (ele->Pt()  > fEleEtMax )) continue;
    if((ele->Eta() < fEleEtaMin) || (ele->Eta() > fEleEtaMax)) continue;
    if(!fEleTools->PassSpikeRemovalFilter(ele)) continue;
    if(fUseGen==ESampleType::kEmbed && !ele->BestTrk()) continue;

    double Electron_ChargedIso03=0, Electron_ChargedIso03FromOtherVertices=0, Electron_NeutralHadronIso03_01Threshold=0, Electron_GammaIso03_01Threshold=0, Electron_NeutralHadronIso03_05Threshold=0, Electron_GammaIso03_05Threshold=0, Electron_NeutralHadronIso03_10Threshold=0, Electron_GammaIso03_10Threshold=0, Electron_NeutralHadronIso03_15Threshold=0, Electron_GammaIso03_15Threshold=0;
    double Electron_ChargedIso04=0, Electron_ChargedIso04FromOtherVertices=0, Electron_NeutralHadronIso04_01Threshold=0, Electron_GammaIso04_01Threshold=0, Electron_NeutralHadronIso04_05Threshold=0, Electron_GammaIso04_05Threshold=0, Electron_NeutralHadronIso04_10Threshold=0, Electron_GammaIso04_10Threshold=0, Electron_NeutralHadronIso04_15Threshold=0, Electron_GammaIso04_15Threshold=0;
    double Electron_ChargedEMIsoVetoEtaStrip03=0, Electron_ChargedEMIsoVetoEtaStrip04=0, Electron_NeutralHadronIso007_01Threshold=0, Electron_GammaIsoVetoEtaStrip03_01Threshold=0, Electron_GammaIsoVetoEtaStrip04_01Threshold=0, Electron_NeutralHadronIso007_05Threshold=0, Electron_GammaIsoVetoEtaStrip03_05Threshold=0, Electron_GammaIsoVetoEtaStrip04_05Threshold=0, Electron_NeutralHadronIso007_10Threshold=0, Electron_GammaIsoVetoEtaStrip03_10Threshold=0, Electron_GammaIsoVetoEtaStrip04_10Threshold=0, Electron_NeutralHadronIso007_15Threshold=0, Electron_GammaIsoVetoEtaStrip03_15Threshold=0, Electron_GammaIsoVetoEtaStrip04_15Threshold=0;
    double tmpPFIso_NormalizedMeanEta=0, tmpPFIso_NormalizedMeanPhi=0;
    double tmpChargedIso_DR0p0To0p1=0, tmpChargedIso_DR0p1To0p2=0, tmpChargedIso_DR0p2To0p3=0, tmpChargedIso_DR0p3To0p4=0, tmpChargedIso_DR0p4To0p5=0, tmpChargedIso_DR0p5To0p7=0, tmpChargedIso_DR0p7To1p0=0;
    double tmpChargedIso_MeanEta_numerator=0, tmpChargedIso_MeanPhi_numerator=0, tmpChargedIso_SigmaEtaEta_numerator=0, tmpChargedIso_SigmaPhiPhi_numerator=0, tmpChargedIso_SigmaEtaPhi_numerator=0, tmpChargedIso_Covariance_denominator=0;
    double tmpGammaIso_DR0p0To0p1=0, tmpGammaIso_DR0p1To0p2=0, tmpGammaIso_DR0p2To0p3=0, tmpGammaIso_DR0p3To0p4=0, tmpGammaIso_DR0p4To0p5=0, tmpGammaIso_DR0p5To0p7=0, tmpGammaIso_DR0p7To1p0=0;
    double tmpGammaIso_MeanEta_numerator=0, tmpGammaIso_MeanPhi_numerator=0, tmpGammaIso_SigmaEtaEta_numerator=0, tmpGammaIso_SigmaPhiPhi_numerator=0, tmpGammaIso_SigmaEtaPhi_numerator=0, tmpGammaIso_Covariance_denominator=0;
    double tmpNeutralHadronIso_DR0p0To0p1=0, tmpNeutralHadronIso_DR0p1To0p2=0, tmpNeutralHadronIso_DR0p2To0p3=0, tmpNeutralHadronIso_DR0p3To0p4=0, tmpNeutralHadronIso_DR0p4To0p5=0, tmpNeutralHadronIso_DR0p5To0p7=0, tmpNeutralHadronIso_DR0p7To1p0=0;
    double tmpNeutralHadronIso_MeanEta_numerator=0, tmpNeutralHadronIso_MeanPhi_numerator=0, tmpNeutralHadronIso_SigmaEtaEta_numerator=0, tmpNeutralHadronIso_SigmaPhiPhi_numerator=0, tmpNeutralHadronIso_SigmaEtaPhi_numerator=0, tmpNeutralHadronIso_Covariance_denominator=0;
  
    const PFCandidate *pfEleMatch = 0;
    double zElectron = 0.0;
    if(ele->BestTrk()) zElectron = ele->BestTrk()->DzCorrected(*fVertex);
    for (unsigned int p=0; p<fPFCandidates->GetEntries();p++) {   
      const PFCandidate *pf = fPFCandidates->At(p);
      //find matching PFElectron
      if (!pfEleMatch) {
        if((pf->GsfTrk() && ele->GsfTrk() && pf->GsfTrk() == ele->GsfTrk()) || (pf->TrackerTrk() && ele->TrackerTrk() && pf->TrackerTrk() == ele->TrackerTrk())) {
          pfEleMatch = pf;
        }
      }
      //exclude the electron itself
      if(pf->GsfTrk() && ele->GsfTrk() && pf->GsfTrk() == ele->GsfTrk()) continue;
      if(pf->TrackerTrk() && ele->TrackerTrk() && pf->TrackerTrk() == ele->TrackerTrk()) continue;      
  
      double dr = MathUtils::DeltaR(ele->Mom(), pf->Mom());
      
      if(pf->BestTrk()) {
        //remove charged particles from other vertices
        double deltaZ = TMath::Abs(pf->BestTrk()->DzCorrected(*fVertex) - zElectron);
  
        if (deltaZ <= 0.1) {
          if (dr < 0.4) Electron_ChargedIso04 += pf->Pt();          
          if (dr < 0.3) Electron_ChargedIso03 += pf->Pt();
          if (pf->PFType() == PFCandidate::eElectron) {
            if (fabs(ele->Eta() - pf->Eta()) < 0.025 && dr < 0.3) Electron_ChargedEMIsoVetoEtaStrip03 += pf->Pt();
            if (fabs(ele->Eta() - pf->Eta()) < 0.025 && dr < 0.4) Electron_ChargedEMIsoVetoEtaStrip04 += pf->Pt();
          } 
        } else {
          if (dr < 0.4) Electron_ChargedIso04FromOtherVertices += pf->Pt();          
          if (dr < 0.3) Electron_ChargedIso03FromOtherVertices += pf->Pt();
        }
      } else if (pf->PFType() == PFCandidate::eGamma) {
        if (dr < 0.4) {
          if (pf->Pt() > 0.1) Electron_GammaIso04_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Electron_GammaIso04_05Threshold += pf->Pt();
          if (pf->Pt() > 1.0) Electron_GammaIso04_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Electron_GammaIso04_15Threshold += pf->Pt();
        }
        if (dr < 0.3) {
          if (pf->Pt() > 0.1) Electron_GammaIso03_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Electron_GammaIso03_05Threshold += pf->Pt();
          if (pf->Pt() > 1.0) Electron_GammaIso03_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Electron_GammaIso03_15Threshold += pf->Pt();
        }
        if (fabs(ele->Eta() - pf->Eta()) < 0.025 && dr < 0.3) {
          if (pf->Pt() > 0.1) Electron_GammaIsoVetoEtaStrip03_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Electron_GammaIsoVetoEtaStrip03_05Threshold += pf->Pt();
          if (pf->Pt() > 1.0) Electron_GammaIsoVetoEtaStrip03_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Electron_GammaIsoVetoEtaStrip03_15Threshold += pf->Pt();
        }
        if (fabs(ele->Eta() - pf->Eta()) < 0.025 && dr < 0.4) {
          if (pf->Pt() > 0.1) Electron_GammaIsoVetoEtaStrip04_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Electron_GammaIsoVetoEtaStrip04_05Threshold += pf->Pt();
          if (pf->Pt() > 1.0) Electron_GammaIsoVetoEtaStrip04_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Electron_GammaIsoVetoEtaStrip04_15Threshold += pf->Pt();
        }
      } else {
        if (dr < 0.4) {
          if (pf->Pt() > 0.1) Electron_NeutralHadronIso04_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Electron_NeutralHadronIso04_05Threshold += pf->Pt();
          if (pf->Pt() > 1.0) Electron_NeutralHadronIso04_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Electron_NeutralHadronIso04_15Threshold += pf->Pt();
        }
        if (dr < 0.3) {
          if (pf->Pt() > 0.1) Electron_NeutralHadronIso03_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Electron_NeutralHadronIso03_05Threshold += pf->Pt();        
          if (pf->Pt() > 1.0) Electron_NeutralHadronIso03_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Electron_NeutralHadronIso03_15Threshold += pf->Pt();
        }
        if (dr < 0.07) {
          if (pf->Pt() > 0.1) Electron_NeutralHadronIso007_01Threshold += pf->Pt();
          if (pf->Pt() > 0.5) Electron_NeutralHadronIso007_05Threshold += pf->Pt();
          if (pf->Pt() > 1.0) Electron_NeutralHadronIso007_10Threshold += pf->Pt();
          if (pf->Pt() > 1.5) Electron_NeutralHadronIso007_15Threshold += pf->Pt();
        }
      }
      
      // New Isolation Calculations
      double deta = (ele->Eta() - pf->Eta());
      double dphi = MathUtils::DeltaPhi(double(ele->Phi()),double(pf->Phi()));
  
      if (dr < 1.0) {
        bool isLeptonFootprint = findLeptonFootprint(pf);
        if (!isLeptonFootprint) {
  	bool passVeto = true;
  	//Charged
  	 if(pf->BestTrk()) {	  	   
  	   if (!(fabs(pf->BestTrk()->DzCorrected(*fVertex) - zElectron) < 0.2)) passVeto = false;
  	   // Veto any PFmuon, or PFEle
  	   if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = false;
  	   // Footprint Veto
  	   if (dr < 0.015) passVeto = false;
  	   if (passVeto) {
  	     if (dr < 0.1) tmpChargedIso_DR0p0To0p1 += pf->Pt();
             if (dr >= 0.1 && dr < 0.2) tmpChargedIso_DR0p1To0p2 += pf->Pt();
  	     if (dr >= 0.2 && dr < 0.3) tmpChargedIso_DR0p2To0p3 += pf->Pt();
  	     if (dr >= 0.3 && dr < 0.4) tmpChargedIso_DR0p3To0p4 += pf->Pt();
  	     if (dr >= 0.4 && dr < 0.5) tmpChargedIso_DR0p4To0p5 += pf->Pt();
  	     if (dr >= 0.5 && dr < 0.7) tmpChargedIso_DR0p5To0p7 += pf->Pt();
  	     if (dr >= 0.7 && dr < 1.0) tmpChargedIso_DR0p7To1p0 += pf->Pt();
  	     tmpChargedIso_Covariance_denominator += pf->Pt();
  	     tmpChargedIso_MeanEta_numerator += pf->Pt() * deta;
  	     tmpChargedIso_MeanPhi_numerator += pf->Pt() * dphi;
  	     tmpChargedIso_SigmaEtaEta_numerator += pf->Pt() * deta*deta;
  	     tmpChargedIso_SigmaPhiPhi_numerator += pf->Pt() * dphi*dphi;
  	     tmpChargedIso_SigmaEtaPhi_numerator += pf->Pt() * deta*dphi;
  	   } //pass veto
  	 }
  	 //Gamma
  	 else if (pf->PFType() == PFCandidate::eGamma) {
  	   // Footprint Veto
  	   if (fabs(ele->Eta()) < 1.479) {
  	     if (fabs(deta) < 0.025) passVeto = false;
  	   } else {
  	     if (dr < 0.08) passVeto = false;
  	   }
  	   if (passVeto) {
  	     if (dr < 0.1) tmpGammaIso_DR0p0To0p1 += pf->Pt();
             if (dr >= 0.1 && dr < 0.2) tmpGammaIso_DR0p1To0p2 += pf->Pt();
  	     if (dr >= 0.2 && dr < 0.3) tmpGammaIso_DR0p2To0p3 += pf->Pt();
  	     if (dr >= 0.3 && dr < 0.4) tmpGammaIso_DR0p3To0p4 += pf->Pt();
  	     if (dr >= 0.4 && dr < 0.5) tmpGammaIso_DR0p4To0p5 += pf->Pt();
  	     if (dr >= 0.5 && dr < 0.7) tmpGammaIso_DR0p5To0p7 += pf->Pt();
  	     if (dr >= 0.7 && dr < 1.0) tmpGammaIso_DR0p7To1p0 += pf->Pt();
  	     tmpGammaIso_Covariance_denominator += pf->Pt();
  	     tmpGammaIso_MeanEta_numerator += pf->Pt() * deta;
  	     tmpGammaIso_MeanPhi_numerator += pf->Pt() * dphi;
  	     tmpGammaIso_SigmaEtaEta_numerator += pf->Pt() * deta*deta;
  	     tmpGammaIso_SigmaPhiPhi_numerator += pf->Pt() * dphi*dphi;
  	     tmpGammaIso_SigmaEtaPhi_numerator += pf->Pt() * deta*dphi;	   	   
  	   }
  	 }
  	 //NeutralHadron
  	 else {
  	   if (dr < 0.1) tmpNeutralHadronIso_DR0p0To0p1 += pf->Pt();
           if (dr >= 0.1 && dr < 0.2) tmpNeutralHadronIso_DR0p1To0p2 += pf->Pt();
  	   if (dr >= 0.2 && dr < 0.3) tmpNeutralHadronIso_DR0p2To0p3 += pf->Pt();
  	   if (dr >= 0.3 && dr < 0.4) tmpNeutralHadronIso_DR0p3To0p4 += pf->Pt();
  	   if (dr >= 0.4 && dr < 0.5) tmpNeutralHadronIso_DR0p4To0p5 += pf->Pt();
  	   if (dr >= 0.5 && dr < 0.7) tmpNeutralHadronIso_DR0p5To0p7 += pf->Pt();
  	   if (dr >= 0.7 && dr < 1.0) tmpNeutralHadronIso_DR0p7To1p0 += pf->Pt();
  	   tmpNeutralHadronIso_Covariance_denominator += pf->Pt();
  	   tmpNeutralHadronIso_MeanEta_numerator += pf->Pt() * deta;
  	   tmpNeutralHadronIso_MeanPhi_numerator += pf->Pt() * dphi;
  	   tmpNeutralHadronIso_SigmaEtaEta_numerator += pf->Pt() * deta*deta;
  	   tmpNeutralHadronIso_SigmaPhiPhi_numerator += pf->Pt() * dphi*dphi;
  	   tmpNeutralHadronIso_SigmaEtaPhi_numerator += pf->Pt() * deta*dphi;
  	 }
  
  	 //For directional iso calculation  
  	 if (dr > 0 && passVeto) {
  	   tmpPFIso_NormalizedMeanEta += pf->Pt() * deta / dr;
  	   tmpPFIso_NormalizedMeanPhi += pf->Pt() * dphi / dr;
  	 }
        } //not lepton footprint
      } //in 1.0 dr cone
    } //loop over PF candidates
  
    TVector2 PFIsoCentroid(tmpPFIso_NormalizedMeanEta,tmpPFIso_NormalizedMeanPhi);
    Double_t tmpPFIsoAngleSqrSum = 0;
    Double_t tmpChargedIsoAngleSqrSum = 0;
    Double_t tmpGammaIsoAngleSqrSum = 0;
    Double_t tmpNeutralHadronIsoAngleSqrSum = 0;
    
  
    for (unsigned int p=0; p<fPFCandidates->GetEntries();p++) {   
      const PFCandidate *pf = fPFCandidates->At(p);
        
      //exclude the electron itself
      if(pf->GsfTrk() && ele->GsfTrk() && pf->GsfTrk() == ele->GsfTrk()) continue;
      if(pf->TrackerTrk() && ele->TrackerTrk() && pf->TrackerTrk() == ele->TrackerTrk()) continue;      
      
      // New Isolation Calculations
      double deta = (ele->Eta() - pf->Eta());
      double dphi = mithep::MathUtils::DeltaPhi(double(ele->Phi()),double(pf->Phi()));
      double dr = mithep::MathUtils::DeltaR(ele->Phi(),ele->Eta(), pf->Phi(), pf->Eta());
      
      if (dr < 1.0) {
        bool isLeptonFootprint = findLeptonFootprint(pf);
        if (!isLeptonFootprint) {
  	TVector2 *tmpPFAngle;
  	if (dr > 0) {
  	  tmpPFAngle = new TVector2(pf->Pt()*deta/dr, pf->Pt()*dphi/dr);
  	} else {
  	  continue;
  	}
  	assert(tmpPFAngle);
  	if (tmpPFAngle->Mod()*PFIsoCentroid.Mod() == 0) continue;
  	if (fabs((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod()) - 1) < 0.001) continue;
  	
  	bool passVeto = true;	
  	//Charged
  	if(pf->BestTrk()) {	  
  	  if (!(fabs(pf->BestTrk()->DzCorrected(*fVertex) - zElectron) < 0.2)) passVeto = false;
  	  // Veto any PFmuon, or PFEle
  	  if (pf->PFType() == PFCandidate::eElectron || pf->PFType() == PFCandidate::eMuon) passVeto = false;
  	  // Footprint Veto
  	  if (dr < 0.015) passVeto = false;
  	  if (passVeto) tmpChargedIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);
  	}
  	//Gamma
  	else if (pf->PFType() == PFCandidate::eGamma) {
  	  // Footprint Veto
 	  if (fabs(ele->Eta()) < 1.479) {
  	    if (fabs(deta) < 0.025) passVeto = false;
  	  } else {
  	    if (dr < 0.08) passVeto = false;
  	  }
  	  if (passVeto) tmpGammaIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);
  	}
  	//NeutralHadron
  	else {
  	  tmpNeutralHadronIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);
  	}
  	if (passVeto) tmpPFIsoAngleSqrSum += pow(acos((tmpPFAngle->Px()*PFIsoCentroid.Px() + tmpPFAngle->Py()*PFIsoCentroid.Py())/(tmpPFAngle->Mod()*PFIsoCentroid.Mod())),2);   
  	if(tmpPFAngle) delete tmpPFAngle;
        } //not lepton footprint  
      } //in 1.0 dr cone
    } //loop over PF candidates
  
    double tmpChargedIso_MeanEta  = (tmpChargedIso_Covariance_denominator>0) ? tmpChargedIso_MeanEta_numerator/tmpChargedIso_Covariance_denominator : 0;
    double tmpChargedIso_MeanPhi  = (tmpChargedIso_Covariance_denominator>0) ? tmpChargedIso_MeanPhi_numerator/tmpChargedIso_Covariance_denominator : 0;
    double tmpChargedIso_SigmaEtaEta  = (tmpChargedIso_Covariance_denominator>0) ? sqrt(tmpChargedIso_SigmaEtaEta_numerator/tmpChargedIso_Covariance_denominator) : 0;
    double tmpChargedIso_SigmaPhiPhi  = (tmpChargedIso_Covariance_denominator>0) ? sqrt(tmpChargedIso_SigmaPhiPhi_numerator/tmpChargedIso_Covariance_denominator) : 0;
    double tmpChargedIso_CovEtaPhi  = (tmpChargedIso_Covariance_denominator>0) ? (tmpChargedIso_SigmaEtaPhi_numerator/tmpChargedIso_Covariance_denominator) : 0;
    double tmpChargedIso_SigmaEtaPhi = 0;
    if (tmpChargedIso_SigmaEtaEta*tmpChargedIso_SigmaPhiPhi != 0) tmpChargedIso_SigmaEtaPhi = tmpChargedIso_CovEtaPhi / (tmpChargedIso_SigmaEtaEta*tmpChargedIso_SigmaPhiPhi);
    else tmpChargedIso_SigmaEtaPhi = (tmpChargedIso_CovEtaPhi < 0) ? -1 : (tmpChargedIso_CovEtaPhi > 0);
    double tmpGammaIso_MeanEta  = (tmpGammaIso_Covariance_denominator>0) ? tmpGammaIso_MeanEta_numerator/tmpGammaIso_Covariance_denominator : 0;
    double tmpGammaIso_MeanPhi  = (tmpGammaIso_Covariance_denominator>0) ? tmpGammaIso_MeanPhi_numerator/tmpGammaIso_Covariance_denominator : 0;
    double tmpGammaIso_SigmaEtaEta  = (tmpGammaIso_Covariance_denominator>0) ? sqrt(tmpGammaIso_SigmaEtaEta_numerator/tmpGammaIso_Covariance_denominator) : 0;
    double tmpGammaIso_SigmaPhiPhi  = (tmpGammaIso_Covariance_denominator>0) ? sqrt(tmpGammaIso_SigmaPhiPhi_numerator/tmpGammaIso_Covariance_denominator) : 0;
    double tmpGammaIso_CovEtaPhi  = (tmpGammaIso_Covariance_denominator>0) ? (tmpGammaIso_SigmaEtaPhi_numerator/tmpGammaIso_Covariance_denominator) : 0;
    double tmpGammaIso_SigmaEtaPhi = 0;
    if (tmpGammaIso_SigmaEtaEta*tmpGammaIso_SigmaPhiPhi != 0) tmpGammaIso_SigmaEtaPhi = tmpGammaIso_CovEtaPhi / (tmpGammaIso_SigmaEtaEta*tmpGammaIso_SigmaPhiPhi);
    else tmpGammaIso_SigmaEtaPhi = (tmpGammaIso_CovEtaPhi < 0) ? -1 : (tmpGammaIso_CovEtaPhi > 0);
    double tmpNeutralHadronIso_MeanEta  = (tmpNeutralHadronIso_Covariance_denominator>0) ? tmpNeutralHadronIso_MeanEta_numerator/tmpNeutralHadronIso_Covariance_denominator : 0;
    double tmpNeutralHadronIso_MeanPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? tmpNeutralHadronIso_MeanPhi_numerator/tmpNeutralHadronIso_Covariance_denominator : 0;
    double tmpNeutralHadronIso_SigmaEtaEta  = (tmpNeutralHadronIso_Covariance_denominator>0) ? sqrt(tmpNeutralHadronIso_SigmaEtaEta_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
    double tmpNeutralHadronIso_SigmaPhiPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? sqrt(tmpNeutralHadronIso_SigmaPhiPhi_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
    double tmpNeutralHadronIso_CovEtaPhi  = (tmpNeutralHadronIso_Covariance_denominator>0) ? (tmpNeutralHadronIso_SigmaEtaPhi_numerator/tmpNeutralHadronIso_Covariance_denominator) : 0;
    double tmpNeutralHadronIso_SigmaEtaPhi = 0;
    if (tmpNeutralHadronIso_SigmaEtaEta*tmpNeutralHadronIso_SigmaPhiPhi != 0) tmpNeutralHadronIso_SigmaEtaPhi = tmpNeutralHadronIso_CovEtaPhi / (tmpNeutralHadronIso_SigmaEtaEta*tmpNeutralHadronIso_SigmaPhiPhi);
    else tmpNeutralHadronIso_SigmaEtaPhi = (tmpNeutralHadronIso_CovEtaPhi < 0) ? -1 : (tmpNeutralHadronIso_CovEtaPhi > 0);
  
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
    pElectron->trkIso04        = ele->TrackIsolationDr04();
    pElectron->emIso04         = ele->EcalRecHitIsoDr04();
    pElectron->hadIso04        = ele->HcalTowerSumEtDr04();
    pElectron->pfIso03         = computePFElecIso(ele,0.3); 
    pElectron->pfIso04         = computePFElecIso(ele,0.4);

    pElectron->chargedIso03 = Electron_ChargedIso03;
    pElectron->chargedIso03FromOtherVertices = Electron_ChargedIso03FromOtherVertices;
    pElectron->neutralHadronIso03_01Threshold = Electron_NeutralHadronIso03_01Threshold;
    pElectron->gammaIso03_01Threshold = Electron_GammaIso03_01Threshold;
    pElectron->neutralHadronIso03_05Threshold = Electron_NeutralHadronIso03_05Threshold;
    pElectron->gammaIso03_05Threshold = Electron_GammaIso03_05Threshold;
    pElectron->neutralHadronIso03_10Threshold = Electron_NeutralHadronIso03_10Threshold;
    pElectron->gammaIso03_10Threshold = Electron_GammaIso03_10Threshold;
    pElectron->neutralHadronIso03_15Threshold = Electron_NeutralHadronIso03_15Threshold;
    pElectron->gammaIso03_15Threshold = Electron_GammaIso03_15Threshold;
    pElectron->chargedIso04 = Electron_ChargedIso04;
    pElectron->chargedIso04FromOtherVertices = Electron_ChargedIso04FromOtherVertices;
    pElectron->neutralHadronIso04_01Threshold = Electron_NeutralHadronIso04_01Threshold;
    pElectron->gammaIso04_01Threshold = Electron_GammaIso04_01Threshold;
    pElectron->neutralHadronIso04_05Threshold = Electron_NeutralHadronIso04_05Threshold;
    pElectron->gammaIso04_05Threshold = Electron_GammaIso04_05Threshold;
    pElectron->neutralHadronIso04_10Threshold = Electron_NeutralHadronIso04_10Threshold;
    pElectron->gammaIso04_10Threshold = Electron_GammaIso04_10Threshold;
    pElectron->neutralHadronIso04_15Threshold = Electron_NeutralHadronIso04_15Threshold;
    pElectron->gammaIso04_15Threshold = Electron_GammaIso04_15Threshold;
    pElectron->chargedEMIsoVetoEtaStrip03 = Electron_ChargedEMIsoVetoEtaStrip03;
    pElectron->chargedEMIsoVetoEtaStrip04 = Electron_ChargedEMIsoVetoEtaStrip04;
    pElectron->neutralHadronIso007_01Threshold = Electron_NeutralHadronIso007_01Threshold ;
    pElectron->gammaIsoVetoEtaStrip03_01Threshold = Electron_GammaIsoVetoEtaStrip03_01Threshold;
    pElectron->gammaIsoVetoEtaStrip04_01Threshold = Electron_GammaIsoVetoEtaStrip04_01Threshold ;
    pElectron->neutralHadronIso007_05Threshold = Electron_NeutralHadronIso007_05Threshold ;
    pElectron->gammaIsoVetoEtaStrip03_05Threshold = Electron_GammaIsoVetoEtaStrip03_05Threshold;
    pElectron->gammaIsoVetoEtaStrip04_05Threshold = Electron_GammaIsoVetoEtaStrip04_05Threshold ;
    pElectron->neutralHadronIso007_10Threshold = Electron_NeutralHadronIso007_10Threshold ;
    pElectron->gammaIsoVetoEtaStrip03_10Threshold = Electron_GammaIsoVetoEtaStrip03_10Threshold;
    pElectron->gammaIsoVetoEtaStrip04_10Threshold = Electron_GammaIsoVetoEtaStrip04_10Threshold ;
    pElectron->neutralHadronIso007_15Threshold = Electron_NeutralHadronIso007_15Threshold ;
    pElectron->gammaIsoVetoEtaStrip03_15Threshold = Electron_GammaIsoVetoEtaStrip03_15Threshold;
    pElectron->gammaIsoVetoEtaStrip04_15Threshold = Electron_GammaIsoVetoEtaStrip04_15Threshold;

    pElectron->chargedIso_DR0p0To0p1 = tmpChargedIso_DR0p0To0p1;
    pElectron->chargedIso_DR0p1To0p2 = tmpChargedIso_DR0p1To0p2;
    pElectron->chargedIso_DR0p2To0p3 = tmpChargedIso_DR0p2To0p3;
    pElectron->chargedIso_DR0p3To0p4 = tmpChargedIso_DR0p3To0p4;
    pElectron->chargedIso_DR0p4To0p5 = tmpChargedIso_DR0p4To0p5;
    pElectron->chargedIso_DR0p5To0p7 = tmpChargedIso_DR0p5To0p7;
    pElectron->chargedIso_DR0p7To1p0 = tmpChargedIso_DR0p7To1p0;
    pElectron->gammaIso_DR0p0To0p1 = tmpGammaIso_DR0p0To0p1;
    pElectron->gammaIso_DR0p1To0p2 = tmpGammaIso_DR0p1To0p2;
    pElectron->gammaIso_DR0p2To0p3 = tmpGammaIso_DR0p2To0p3;
    pElectron->gammaIso_DR0p3To0p4 = tmpGammaIso_DR0p3To0p4;
    pElectron->gammaIso_DR0p4To0p5 = tmpGammaIso_DR0p4To0p5;
    pElectron->gammaIso_DR0p5To0p7 = tmpGammaIso_DR0p5To0p7;
    pElectron->gammaIso_DR0p7To1p0 = tmpGammaIso_DR0p7To1p0;
    pElectron->neutralIso_DR0p0To0p1 = tmpNeutralHadronIso_DR0p0To0p1;
    pElectron->neutralIso_DR0p1To0p2 = tmpNeutralHadronIso_DR0p1To0p2;
    pElectron->neutralIso_DR0p2To0p3 = tmpNeutralHadronIso_DR0p2To0p3;
    pElectron->neutralIso_DR0p3To0p4 = tmpNeutralHadronIso_DR0p3To0p4;
    pElectron->neutralIso_DR0p4To0p5 = tmpNeutralHadronIso_DR0p4To0p5;
    pElectron->neutralIso_DR0p5To0p7 = tmpNeutralHadronIso_DR0p5To0p7;
    pElectron->neutralIso_DR0p7To1p0 = tmpNeutralHadronIso_DR0p7To1p0;
    pElectron->chargedIso_MeanEta = tmpChargedIso_MeanEta;
    pElectron->chargedIso_MeanPhi = tmpChargedIso_MeanPhi;
    pElectron->chargedIso_SigEtaEta = tmpChargedIso_SigmaEtaEta;
    pElectron->chargedIso_SigPhiPhi = tmpChargedIso_SigmaPhiPhi;
    pElectron->chargedIso_SigEtaPhi = tmpChargedIso_SigmaEtaPhi;
    pElectron->gammaIso_MeanEta = tmpGammaIso_MeanEta;
    pElectron->gammaIso_MeanPhi = tmpGammaIso_MeanPhi;
    pElectron->gammaIso_SigEtaEta = tmpGammaIso_SigmaEtaEta;
    pElectron->gammaIso_SigPhiPhi = tmpGammaIso_SigmaPhiPhi;
    pElectron->gammaIso_SigEtaPhi = tmpGammaIso_SigmaEtaPhi;
    pElectron->neutralIso_MeanEta = tmpNeutralHadronIso_MeanEta;
    pElectron->neutralIso_MeanPhi = tmpNeutralHadronIso_MeanPhi;
    pElectron->neutralIso_SigEtaEta = tmpNeutralHadronIso_SigmaEtaEta;
    pElectron->neutralIso_SigPhiPhi = tmpNeutralHadronIso_SigmaPhiPhi;
    pElectron->neutralIso_SigEtaPhi = tmpNeutralHadronIso_SigmaEtaPhi;
    pElectron->directionalPFIso = tmpPFIsoAngleSqrSum;
    pElectron->directionalChargedIso = tmpChargedIsoAngleSqrSum;
    pElectron->directionalGammaIso = tmpGammaIsoAngleSqrSum;
    pElectron->directionalNeutralIso = tmpNeutralHadronIsoAngleSqrSum;

    pElectron->d0              = ele->BestTrk()->D0Corrected(*fVertex);
    pElectron->dz              = ele->BestTrk()->DzCorrected(*fVertex);  
    pElectron->d0Sig           = ele->D0PVSignificance();
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
    pElectron->E               = ele->SCluster()->Energy();
    pElectron->p               = ele->BestTrk()->P();
    pElectron->ip3d            = ele->Ip3dPV();
    pElectron->ip3dSig         = ele->Ip3dPVSignificance();
    pElectron->ESeedClusterOverPIn   = ele->ESeedClusterOverPIn();
    pElectron->ESeedClusterOverPOut  = ele->ESeedClusterOverPout();
    pElectron->sigiPhiiPhi     = ele->SCluster()->Seed()->CoviPhiiPhi();
    pElectron->nBrem           = ele->NumberOfClusters()-1;  
    pElectron->hltMatchBits    = matchHLT(ele->SCluster()->Eta(), ele->SCluster()->Phi(), ele->SCluster()->Et());
    pElectron->scID            = ele->SCluster()->GetUniqueID();
    pElectron->trkID           = (ele->HasTrackerTrk()) ? ele->TrackerTrk()->GetUniqueID() : 0;
    pElectron->isEcalDriven    = ele->IsEcalDriven();
    pElectron->isConv          = isConversion(ele);
    pElectron->isEB            = ele->IsEB();
    pElectron->ip3d            = ele->Ip3dPV();
    pElectron->ip3dSig         = ele->Ip3dPVSignificance();
    pElectron->d0Ub            = ele->D0PVUB();
    pElectron->d0UbSig         = ele->D0PVUBSignificance();
    pElectron->ip3dUb          = ele->Ip3dPVUB();
    pElectron->ip3dUbSig       = ele->Ip3dPVUBSignificance();

    // additional variables used for mva id
    pElectron->gsfTrackChi2OverNdof = ele->BestTrk()->Chi2() / ele->BestTrk()->Ndof();
    pElectron->hcalDepth1OverEcal = ele->HcalDepth1OverEcal();
    pElectron->hcalDepth2OverEcal = ele->HcalDepth2OverEcal();
    pElectron->deltaEtaCalo    = ele->DeltaEtaSeedClusterTrackAtCalo();
    pElectron->deltaPhiCalo    = ele->DeltaPhiSeedClusterTrackAtCalo();
    pElectron->R9              = ele->SCluster()->R9();
    pElectron->scEtaWidth      = ele->SCluster()->EtaWidth();
    pElectron->scPhiWidth      = ele->SCluster()->PhiWidth();
    pElectron->coviEtaiPhi     = ele->SCluster()->Seed()->CoviEtaiPhi();
    pElectron->psOverRaw       = ele->SCluster()->PreshowerEnergy() / ele->SCluster()->RawEnergy();
    pElectron->seedEMaxOverE = ele->SCluster()->Seed()->EMax() / ele->SCluster()->Seed()->Energy();
    pElectron->seedETopOverE = ele->SCluster()->Seed()->ETop() / ele->SCluster()->Seed()->Energy();
    pElectron->seedEBottomOverE = ele->SCluster()->Seed()->EBottom() / ele->SCluster()->Seed()->Energy();
    pElectron->seedELeftOverE = ele->SCluster()->Seed()->ELeft() / ele->SCluster()->Seed()->Energy();
    pElectron->seedERightOverE = ele->SCluster()->Seed()->ERight() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE2ndOverE = ele->SCluster()->Seed()->E2nd() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE2x5RightOverE = ele->SCluster()->Seed()->E2x5Right() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE2x5LeftOverE = ele->SCluster()->Seed()->E2x5Left() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE2x5TopOverE = ele->SCluster()->Seed()->E2x5Top() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE2x5BottomOverE = ele->SCluster()->Seed()->E2x5Bottom() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE2x5MaxOverE = ele->SCluster()->Seed()->E2x5Max() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE1x3OverE = ele->SCluster()->Seed()->E1x3() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE3x1OverE = ele->SCluster()->Seed()->E3x1() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE1x5OverE = ele->SCluster()->Seed()->E1x5() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE2x2OverE = ele->SCluster()->Seed()->E2x2() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE3x2OverE = ele->SCluster()->Seed()->E3x2() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE3x3OverE = ele->SCluster()->Seed()->E3x3() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE4x4OverE = ele->SCluster()->Seed()->E4x4() / ele->SCluster()->Seed()->Energy();
    pElectron->seedE5x5OverE = ele->SCluster()->Seed()->E5x5() / ele->SCluster()->Seed()->Energy();

    pElectron->pfIsoCharged    = computeCommonIso(ele, fPFNoPileUp    , 0.0, 0.4, 0.015,  1);
    pElectron->pfIsoChargedNoZ = computeCommonIso(ele, fPFNoPileUpNoZ , 0.0, 0.4, 0.015,  1);
    pElectron->pfIsoNeutral    = computeCommonIso(ele, fPFNoPileUp    , 0.5, 0.4,   0.0,  2);
    pElectron->pfIsoNeutralNoZ = computeCommonIso(ele, fPFNoPileUpNoZ , 0.5, 0.4,   0.0,  2);
    pElectron->pfIsoGamma      = computeCommonIso(ele, fPFNoPileUp    , 0.5, 0.4,  0.08,  3);
    pElectron->pfIsoGammaNoZ   = computeCommonIso(ele, fPFNoPileUpNoZ , 0.5, 0.4,  0.08,  3);
    pElectron->puIso           = computeNaiveIso(ele, fPFPileUp    , 0.5, 0.4,  0.01);
    pElectron->puIsoNoZ        = computeNaiveIso(ele, fPFPileUpNoZ , 0.5, 0.4,  0.01);
    pElectron->pfPx=0; pElectron->pfPy=0;
    for(unsigned int i=0; i<fPFCandidates->GetEntries(); ++i){
      const PFCandidate* pfcand = fPFCandidates->At(i);
      if( (pfcand->HasTrackerTrk() && ele->TrackerTrk() == pfcand->TrackerTrk()) || (pfcand->HasGsfTrk() && ele->GsfTrk() == pfcand->GsfTrk()) ){
	pElectron->pfPx = pfcand->Px(); pElectron->pfPy = pfcand->Py();
	break;
      }	 
    }
  }
}

void 
HttNtupler::fillPFTaus() 
{
  assert(fPFTaus);
  for(unsigned idx=0; idx<fPFTaus->GetEntries(); ++idx){
    const PFTau* pftau = fPFTaus->At(idx);
    if(pftau->Pt() > fPFTauPtMin && fabs(pftau->Eta()) < fPFTauEtaMax ){
      TClonesArray &rPFTauArr = *fHPSTauArr;
      assert(rPFTauArr.GetEntries() < rPFTauArr.GetSize());
      const Int_t index = rPFTauArr.GetEntries();
      new(rPFTauArr[index]) TPFTau();
      TPFTau* pPFTau = (TPFTau*)rPFTauArr[index]; 
      pPFTau->pt                    = pftau->Pt();
      pPFTau->eta                   = pftau->Eta();
      pPFTau->phi                   = pftau->Phi();
      pPFTau->m                     = pftau->Mass();
      pPFTau->e                     = pftau->E();
      pPFTau->q                     = pftau->Charge();//charge(pftau);//->Charge();
      pPFTau->leadPFCandSignD0Sig   = pftau->LeadPFCandSignD0Sig();
      pPFTau->isoChargedHadronPtSum = pftau->IsoChargedHadronPtSum();
      pPFTau->isoGammaEtSum         = pftau->IsoGammaEtSum();
      pPFTau->hcal3x3EOverP         = pftau->HCal3x3EOverP();
      pPFTau->elePreIDOutput        = pftau->ElectronPreIDOutput();
      pPFTau->emFraction            = pftau->EMFraction();
      const PFCandidate* leadPF     = pftau->LeadChargedHadronPFCand();
      if( leadPF != 0 ){
	pPFTau->hasLeadChargedHadronPFCand       = true;
	pPFTau->hcalOverP                        = leadPF->EHCal()/leadPF->P();
	pPFTau->ecalOverP                        = leadPF->EECal()/leadPF->P();
	pPFTau->hasGsf                           = leadPF->HasGsfTrk();
	pPFTau->leadChargedHadronPFCand.pt       = leadPF->Pt();
	pPFTau->leadChargedHadronPFCand.eta      = leadPF->Eta();
	pPFTau->leadChargedHadronPFCand.phi      = leadPF->Phi();
	pPFTau->leadChargedHadronPFCand.m        = leadPF->Mass();
	pPFTau->leadChargedHadronPFCand.e        = leadPF->E();
	pPFTau->leadChargedHadronPFCand.q        = leadPF->Charge();
	pPFTau->leadChargedHadronPFCand.mvaEPi   = leadPF->MvaEPi();
	pPFTau->leadChargedHadronPFCand.pfType   = leadPF->PFType();
	pPFTau->leadChargedHadronPFCand.hasMuon  = leadPF->Mu() ? kTRUE : kFALSE;
	pPFTau->leadChargedHadronPFCand.nSeg     = pPFTau->leadChargedHadronPFCand.hasMuon ? leadPF->Mu()->NSegments() : 0;
	pPFTau->leadChargedHadronPFCand.nMatches = pPFTau->leadChargedHadronPFCand.hasMuon ? leadPF->Mu()->NMatches()  : 0;
	if( leadPF->BestTrk() ){
	  pPFTau->leadChargedHadronPFCand.hasTrack = kTRUE;
	  pPFTau->leadChargedHadronPFCand.d0       = leadPF->BestTrk()->D0Corrected(*fVertex);
	  pPFTau->leadChargedHadronPFCand.dz       = leadPF->BestTrk()->DzCorrected(*fVertex);
	}
      }
      pPFTau->nSignalPFChargedHadrCands = pftau->NSignalPFChargedHadrCands(); 
      pPFTau->nSignalPFGammaCands       = pftau->NSignalPFGammaCands(); 
      pPFTau->gammaPtR                  = 0;
      pPFTau->gammaDEta2                = 0;
      pPFTau->gammaDPhi2                = 0;
      for(unsigned int icand=0; icand<pftau->NSignalPFGammaCands(); ++icand){
	const PFCandidate* pf = pftau->SignalPFGammaCand(icand);
	double lDEta = fabs(pf->Eta() - pftau->Eta());
	double lDPhi = fabs(pf->Phi() - pftau->Phi());
	if(leadPF) {
	  lDEta =  fabs(pf->Eta() - leadPF->Eta());
	  lDPhi =  fabs(pf->Phi() - leadPF->Phi());
	}
	if(lDPhi > 2.*TMath::Pi()-lDPhi) lDPhi = 2.*TMath::Pi()-lDPhi;
	if(icand == 0){
	  pPFTau->gammaPt    = pf->Pt();
	  pPFTau->gammaDEta  = lDEta;
	  pPFTau->gammaDPhi  = lDPhi;
	}
	if(icand == 1){
	  pPFTau->gamma2Pt   = pf->Pt();
	  pPFTau->gamma2DEta = lDEta;
	  pPFTau->gamma2DPhi = lDPhi;
	}
	pPFTau->gammaPtR    += pf->Pt();
	pPFTau->gammaDEta2  += pf->Pt()*lDEta*lDEta;
	pPFTau->gammaDPhi2  += pf->Pt()*lDPhi*lDPhi;
      }
      pPFTau->gammaPtR   /= pftau->Pt();
      pPFTau->gammaDEta2  = TMath::Sqrt(pPFTau->gammaDEta2)*TMath::Sqrt(pPFTau->gammaPtR)*pftau->Pt();
      pPFTau->gammaDPhi2  = TMath::Sqrt(pPFTau->gammaDPhi2)*TMath::Sqrt(pPFTau->gammaPtR)*pftau->Pt();
      pPFTau->Iso         = computeOfficialPFTauIso(pftau);
      pPFTau->puIso       = computeCommonIso((const ChargedParticle*) pftau,                   fPFPileUp    , 0.5, 0.3, 0.0001, 0);
      pPFTau->puIsoNoZ    = computeCommonIso((const ChargedParticle*) pftau, (PFCandidateCol*)fPFCandidates, 0.5, 0.5, 0.0001, 0);
      pPFTau->puIsoNoPt   = computeCommonIso((const ChargedParticle*) pftau,                   fPFPileUp    ,  0., 0.5,     0., 0);
      if(pftau->LeadChargedHadronPFCand() && pftau->LeadChargedHadronPFCand()->Trk()) {
	pPFTau->isoEtPU = computePFTauIso(pftau, pftau->LeadChargedHadronPFCand()->Trk());
      }
      if( pftau->DiscriminationByLooseElectronRejection()             ) pPFTau->hpsDiscriminators |= TPFTau::kLooseEle;
      if( pftau->DiscriminationByMediumElectronRejection()            ) pPFTau->hpsDiscriminators |= TPFTau::kMediumEle;
      if( pftau->DiscriminationByTightElectronRejection()             ) pPFTau->hpsDiscriminators |= TPFTau::kTightEle;
      if( pftau->DiscriminationByLooseMuonRejection()                 ) pPFTau->hpsDiscriminators |= TPFTau::kLooseMu;
      if( pftau->DiscriminationByTightMuonRejection()                 ) pPFTau->hpsDiscriminators |= TPFTau::kTightMu;
      if( pftau->DiscriminationByDecayModeFinding()                   ) pPFTau->hpsDiscriminators |= TPFTau::kDecayMode;
      if( pftau->DiscriminationByVLooseIsolation()                    ) pPFTau->hpsDiscriminators |= TPFTau::kVLooseIso;
      if( pftau->DiscriminationByLooseIsolation()                     ) pPFTau->hpsDiscriminators |= TPFTau::kLooseIso;
      if( pftau->DiscriminationByMediumIsolation()                    ) pPFTau->hpsDiscriminators |= TPFTau::kMediumIso;
      if( pftau->DiscriminationByTightIsolation()                     ) pPFTau->hpsDiscriminators |= TPFTau::kTightIso;
      if( pftau->DiscriminationByVLooseCombinedIsolationDBSumPtCorr() ) pPFTau->hpsDiscriminators |= TPFTau::kVLooseCombIso;
      if( pftau->DiscriminationByLooseCombinedIsolationDBSumPtCorr()  ) pPFTau->hpsDiscriminators |= TPFTau::kLooseCombIso;
      if( pftau->DiscriminationByMediumCombinedIsolationDBSumPtCorr() ) pPFTau->hpsDiscriminators |= TPFTau::kMediumCombIso;
      if( pftau->DiscriminationByTightCombinedIsolationDBSumPtCorr()  ) pPFTau->hpsDiscriminators |= TPFTau::kTightCombIso;
      // HLT matching
      pPFTau->hltMatchBits = matchHLT(pftau->Eta(), pftau->Phi(), pftau->Pt());
      for(unsigned int itrig=0; itrig<fTriggerNamesv.size(); ++itrig){
	if( pPFTau->hltMatchBits[fTriggerObjIds1v[itrig]] && !pPFTau->hltMatchBits[fTriggerObjIds2v[itrig]]) pPFTau->hltMatchBits[fTriggerIdsv[itrig]] = false;
      }
    }
  }
}

void 
HttNtupler::fillJets()
{
  assert(fPFJets);
  for(unsigned int i=0; i<fPFJets->GetEntries(); ++i){
    const PFJet* jet = fPFJets->At(i);
    const FourVectorM rawMom = jet->RawMom();
    fJetCorrector->setJetEta( rawMom.Eta() );
    fJetCorrector->setJetPt ( rawMom.Pt()  );
    fJetCorrector->setJetPhi( rawMom.Phi() );
    fJetCorrector->setJetE  ( rawMom.E()   );
    fJetCorrector->setRho   ( fPUEnergyDensity->At(0)->RhoHighEta() );
    fJetCorrector->setJetA  ( jet->JetArea() );
    fJetCorrector->setJetEMF( -99.0 );        
    double correction = fJetCorrector->getCorrection();
    double pt = rawMom.Pt()*correction;
    if(pt > fJetPtMin || (jet->TrackCountingHighEffBJetTagsDisc() != -100)){
      // loose jetId for particle flow jets. NOTE: energy fractions are on the raw (i.e. uncorrected energy)
      if(jet->E()==0) continue;
      if(jet->NeutralHadronEnergy()/jet->E()   >  0.99                       ) continue;
      if(jet->NeutralEmEnergy()/jet->E()       >  0.99                       ) continue;
      if(jet->NConstituents()                  <  2                          ) continue;
      if(jet->ChargedHadronEnergy()/jet->E() <= 0  && fabs(jet->Eta()) < 2.4 ) continue;
      if(jet->ChargedEmEnergy()/jet->E()  >  0.99  && fabs(jet->Eta()) < 2.4 ) continue;
      if(jet->ChargedMultiplicity()    < 1         && fabs(jet->Eta()) < 2.4 ) continue;
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
      fJetUncertainties->setJetPt ( rawMom.Pt()  );
      fJetUncertainties->setJetEta( rawMom.Eta() );
      pPFJet->pt          = rawMom.Pt()*fJetCorrector->getCorrection();
      pPFJet->eta         = rawMom.Eta();
      pPFJet->phi         = rawMom.Phi();
      pPFJet->mass        = jet->Mass();
      pPFJet->unc         = fJetUncertainties->getUncertainty(true);
      pPFJet->ptraw       = rawMom.Pt();
      pPFJet->beta        = jet->Beta();
      pPFJet->area        = jet->JetArea();
      pPFJet->nCharged    = jet->ChargedMultiplicity();
      pPFJet->chgEMfrac   = jet->ChargedEmEnergy()/jet->E();
      pPFJet->neuEMfrac   = jet->NeutralEmEnergy()/jet->E();
      pPFJet->chgHadrfrac = jet->ChargedHadronEnergy()/jet->E();
      pPFJet->neuHadrfrac = jet->NeutralHadronEnergy()/jet->E();
      pPFJet->tche        = jet->TrackCountingHighEffBJetTagsDisc();
      pPFJet->tchp        = jet->TrackCountingHighPurBJetTagsDisc();
      pPFJet->csv         = jet->CombinedSecondaryVertexBJetTagsDisc();
      pPFJet->csvMva      = jet->CombinedSecondaryVertexMVABJetTagsDisc();
      pPFJet->mcFlavor    = jet->MatchedMCFlavor();
      int matchedFlavor   = -999;
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
      pPFJet->matchedFlavor = matchedFlavor;

      int matchedId       = -999;
      if (fParticles) {
        Double_t dRmin = 0.3;
        for (UInt_t i=0; i<fParticles->GetEntries(); ++i) {
          const MCParticle *p = fParticles->At(i);
          if((p->Status()==3 || (p->Status()==2 && !p->HasDaughter(1) && !p->HasDaughter(2) && !p->HasDaughter(3) && !p->HasDaughter(4) && !p->HasDaughter(5) && !p->HasDaughter(6) && !p->HasDaughter(21))) && p->IsParton()) {
            if(MathUtils::DeltaR(*jet,*p) < dRmin) {
              dRmin = MathUtils::DeltaR(*jet,*p);
              matchedId = p->PdgId();
            }
          }
        }
      }

      pPFJet->matchedId    = matchedId;
      pPFJet->hltMatchBits = matchHLT(jet->Eta(), jet->Phi(), jet->Pt());
    }
  }
}

void 
HttNtupler::fillPhotons()
{
  assert(fPhotons);
  for(unsigned int i=0; i<fPhotons->GetEntries(); ++i){
    const Photon* pho= fPhotons->At(i);
    if(pho->SCluster()->Et() > fPhotonEtMin) {
      TClonesArray& rPhotonArr = *fPhotonArr;
      assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
      const int index = rPhotonArr.GetEntries();  
      new(rPhotonArr[index]) TPhoton();
      TPhoton* pPhoton = (TPhoton*)rPhotonArr[index];
      pPhoton->pt	    = pho->Pt(); 
      pPhoton->eta  	    = pho->Eta();
      pPhoton->phi  	    = pho->Phi();
      pPhoton->scEt	    = pho->SCluster()->Et(); 
      pPhoton->scEta  	    = pho->SCluster()->Eta();
      pPhoton->scPhi  	    = pho->SCluster()->Phi();
      pPhoton->trkIso04     = pho->HollowConeTrkIsoDr04();
      pPhoton->emIso04      = pho->EcalRecHitIsoDr04();
      pPhoton->hadIso04	    = pho->HcalTowerSumEtDr04(); 
      pPhoton->HoverE	    = pho->HadOverEm();
      pPhoton->R9	    = pho->R9();
      pPhoton->sigiEtaiEta  = pho->CoviEtaiEta();
      pPhoton->hltMatchBits = matchHLT(pho->SCluster()->Eta(), pho->SCluster()->Phi(), pho->SCluster()->Et());
      pPhoton->scID         = pho->SCluster()->GetUniqueID();
      pPhoton->hasPixelSeed = pho->HasPixelSeed();
    }
  }
}

void 
HttNtupler::fillSVfit() 
{ 
  for(unsigned int idx=0; idx<fMuons->GetEntries(); ++idx){ 
    const Muon* pMu = fMuons->At(idx); if( !looseMuId(pMu) ){ continue; }
    for(unsigned int jdx=0; jdx<fElectrons->GetEntries(); ++jdx){ 
      const Electron* pElectron = fElectrons->At(jdx); if( !looseEleId(pElectron) ){ continue; }
      if( MathUtils::DeltaR(pElectron->Mom(), pMu->Mom()) < 0.3 ){ continue; }
      TMatrixD lMetMatrix = metSign->getSignificance(fPFJets, fPFCandidates, 0, pMu, pElectron);	
      double dcaSig3D = -999, dcaSig2D = -999, dca3DErr=-999, dca2DErr=-999;
      /*for(unsigned int  i=0; i<fDCASigs->GetEntries(); ++i){
	const DCASig* dca=fDCASigs->At(i);
	if(dca->Type()!=DCASig::eEMu) continue;
	if(MathUtils::DeltaR(pElectron->Mom(), dca->GetElectron()->Mom()) > 0.1) continue;
        if(MathUtils::DeltaR(pMu->Mom(), dca->GetMuon()->Mom()) > 0.1) continue;
	dcaSig3D = dca->DCASig3D();
        dcaSig2D = dca->DCASig2D();
	dca3DErr = dca->DCA3DErr();
        dca2DErr = dca->DCA2DErr();
      }*/
      fillSVfit(fSVfitEMuArr, (Particle*)pMu, EGenType::kMuon, (Particle*)pElectron, EGenType::kElectron, lMetMatrix, dcaSig3D, dcaSig2D, dca3DErr, dca2DErr);
    }
  }
  for(unsigned int idx=0; idx<fPFTaus->GetEntries(); ++idx){ 
    const PFTau* pPFTau = fPFTaus->At(idx); if( !looseTauId(pPFTau) ){ continue; }
    for(unsigned int jdx=0; jdx<fMuons->GetEntries(); ++jdx){
      const Muon* pMu = fMuons->At(jdx); if( !looseMuId(pMu) ){ continue; }
      if( MathUtils::DeltaR(pPFTau->Mom(), pMu->Mom())<0.3 ){ continue; }
      TMatrixD lMetMatrix = metSign->getSignificance(fPFJets, fPFCandidates, pPFTau, pMu, 0);
      //std::cout << "fill" << std::endl;	
      fillSVfit(fSVfitMuTauArr, (Particle*)pMu, EGenType::kMuon, (Particle*)pPFTau, EGenType::kTau, lMetMatrix);
    }
    for(unsigned int jdx=0; jdx<fElectrons->GetEntries(); ++jdx){ 
      const Electron* pElectron = fElectrons->At(jdx); if( !looseEleId(pElectron) ){ continue; }
      if( MathUtils::DeltaR(pPFTau->Mom(), pElectron->Mom())<0.3 ){ continue; }
      TMatrixD lMetMatrix = metSign->getSignificance(fPFJets, fPFCandidates, pPFTau, 0, pElectron);	
      fillSVfit(fSVfitETauArr, (Particle*)pElectron, EGenType::kElectron, (Particle*)pPFTau, EGenType::kTau, lMetMatrix);
    }
  }
}

void 
HttNtupler::fillSVfit(TClonesArray*& iArr, Particle* lep1, unsigned int lepId1, Particle* lep2, unsigned int lepId2, TMatrixD iMatrix, double dcaSig3D, double dcaSig2D, double dca3DErr, double dca2DErr) 
{
  TClonesArray& rSVfitArr = *iArr;
  const int index = rSVfitArr.GetEntries();
  new(rSVfitArr[index]) TSVfit();
  TSVfit* pSVfit = (TSVfit*)rSVfitArr[index];
  pSVfit->daughter1 = lep1->Mom(); pSVfit->daughterId1 = lepId1;
  pSVfit->daughter2 = lep2->Mom(); pSVfit->daughterId2 = lepId2;
  pSVfit->cov_00 = iMatrix(0,0)  ; pSVfit->cov_10 = iMatrix(1,0);
  pSVfit->cov_01 = iMatrix(0,1)  ; pSVfit->cov_11 = iMatrix(1,1);
  pSVfit->dcaSig3D = dcaSig3D;
  pSVfit->dcaSig2D = dcaSig2D;
  pSVfit->dcaSig3DErr = dca3DErr;
  pSVfit->dcaSig2DErr = dca2DErr;
}

void
HttNtupler::fillSVfit(TClonesArray*& iArr, Particle* lep1, unsigned int lepId1, Particle* lep2, unsigned int lepId2, TMatrixD iMatrix)
{
  TClonesArray& rSVfitArr = *iArr;
  const int index = rSVfitArr.GetEntries();
  new(rSVfitArr[index]) TSVfit();
  TSVfit* pSVfit = (TSVfit*)rSVfitArr[index];
  pSVfit->daughter1 = lep1->Mom(); pSVfit->daughterId1 = lepId1;
  pSVfit->daughter2 = lep2->Mom(); pSVfit->daughterId2 = lepId2;
  pSVfit->cov_00 = iMatrix(0,0)  ; pSVfit->cov_10 = iMatrix(1,0);
  pSVfit->cov_01 = iMatrix(0,1)  ; pSVfit->cov_11 = iMatrix(1,1);
}

void 
HttNtupler::fillMCParticles(const MCParticle*& boson, const MCParticle*& daughterA, const MCParticle*& daughterB) 
{
  // descend as long as boson is its own daughter in the listing
  while(boson->NDaughters()==1){ boson = boson->FindDaughter(boson->PdgId()); }
  for(unsigned int idau=0; idau<boson->NDaughters(); ++idau){
    const MCParticle* d = boson->Daughter(idau);
    if(d->PdgId()==boson->PdgId()){ continue; }
    if(d->PdgId()==22){ continue; }
    if(d->PdgId() > 0) daughterA = d;
    if(d->PdgId() < 0) daughterB = d;
  }
}

FourVectorM
HttNtupler::visibleMCMomentum(const MCParticle* tauLepton) 
{ 
  FourVectorM visMomentum;
  // descend as long as particle is its own daughter in the listing
  while(tauLepton->NDaughters()==1){ tauLepton = tauLepton->FindDaughter(tauLepton->PdgId()); }
  visMomentum+=tauLepton->Mom();
  // loop tau daughters and subract neutrino momenta from the original tau momentum
  for(unsigned int idx=0; idx<tauLepton->NDaughters(); ++idx){ 
    const MCParticle* daughter = tauLepton->Daughter(idx);
    if( !(daughter->AbsPdgId() == EGenType::kTNeutrino || daughter->AbsPdgId() == EGenType::kMNeutrino || daughter->AbsPdgId() == EGenType::kENeutrino) ){ continue; }
    visMomentum-=daughter->Mom();
  }
  return visMomentum;
}

void 
HttNtupler::stableMCParticle(const MCParticle*& part, int status) 
{
  if(status==0){ while(part->HasDaughter(part->PdgId(), true)){ part = part->FindDaughter(part->PdgId(), true); } } 
  else{ while(part->HasDaughter(part->PdgId(), true) && (part->Status()!=status)){ part = part->FindDaughter(part->PdgId(), true); } }
}

int 
HttNtupler::pdgId(const MCParticle* part) 
{ 
  int id = part->PdgId();
  // for everything but taus the Id is the actual pdgId
  if(fabs(id) !=  EGenType::kTau){ return id; }
  const MCParticle* tau = part; stableMCParticle(tau, 1);
  for(unsigned int idx=0; idx<tau->NDaughters(); ++idx){ 
    const MCParticle* daughter = tau->Daughter(idx);
    if( daughter->AbsPdgId() == EGenType::kTNeutrino || daughter->AbsPdgId() == EGenType::kMNeutrino || daughter->AbsPdgId() == EGenType::kENeutrino ){ continue; }
    if( daughter->PdgId() == EGenType::kPhoton ){ continue; }
    if( fabs(daughter->PdgId()) == EGenType::kElectron ){ return daughter->Charge() * EGenType::kTauElectron; }
    if( fabs(daughter->PdgId()) == EGenType::kMuon ){ return daughter->Charge() * EGenType::kTauMuon; }
    return daughter->Charge() * EGenType::kTauHadr;
  }
  return 0;
}

TriggerObjects 
HttNtupler::matchHLT(const double eta, const double phi, const double pt) 
{
  TriggerObjects toBits=0;
  if( !HasHLTInfo() ) return toBits;  
  const double hltMatchR      = 0.5;
  const double hltMatchPtFrac = 0.5;
  const TriggerTable* hltTable = GetHLTTable(); assert(hltTable);
  for(unsigned int itrig=0; itrig<fTriggerNamesv.size(); ++itrig){
    const TriggerName* trigname = hltTable->Get(fTriggerNamesv[itrig].Data()); if(!trigname) continue;
    const TList* list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data()); if(!list) continue;
    TIter iter(list->MakeIterator());
    const TriggerObject* to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    while( to ){         
      if( to->IsHLT() ){
	if( fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0 ){
	  bool match = true;
	  if( to->Pt() < fTriggerObjMinPt1v[itrig] ) match=false;
	  if( MathUtils::DeltaR(phi, eta, to->Phi(), to->Eta()) > hltMatchR ) match=false;
	  if( hltMatchPtFrac>0 && (fabs(pt-to->Pt() )>hltMatchPtFrac*(to->Pt())) ) match=false;
	  if( match ) {toBits[fTriggerObjIds1v[itrig]] = true; }
	}
	if( fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0 ){
	  bool match = true;
	  if( to->Pt() < fTriggerObjMinPt2v[itrig] ) match=false;
	  if( MathUtils::DeltaR(phi, eta, to->Phi(), to->Eta()) > hltMatchR ) match=false;
	  if( hltMatchPtFrac>0 && (fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt())) ) match=false;
	  if( match ) {toBits[fTriggerObjIds2v[itrig]] = true; }
	}
	if( fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0 ){
	  bool match = true;
	  if( to->Pt() < fTriggerObjMinPt1v[itrig] ) match=false;
	  if( MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR ) match=false;
	  if( hltMatchPtFrac>0 && (fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt())) ) match=false;
	  if( match ) {toBits[fTriggerObjIds1v[itrig]] = true; }
	}
      }
      to = dynamic_cast<const TriggerObject*>(iter.Next());
    }    
  }
  return toBits;
}

bool 
HttNtupler::isConversion(const Electron* ele) 
{
  bool isGoodConversion              = false;
  const unsigned int nWrongHitsMax   = 0;
  const double probMin               = 1e-6;
  const double lxyMin                =  2.0;
  const bool matchCkf                = true;
  const bool requireArbitratedMerged = false;
  for(unsigned int ifc=0; ifc<fConversions->GetEntries(); ++ifc){
    bool conversionMatchFound = false;
    for(unsigned int d=0; d<fConversions->At(ifc)->NDaughters(); ++d){
      const Track* trk = dynamic_cast<const ChargedParticle*>(fConversions->At(ifc)->Daughter(d))->Trk();
      if( ele->GsfTrk() == trk || (matchCkf && ele->TrackerTrk()==trk) ){ conversionMatchFound = true; break; }
    }
    if( conversionMatchFound ){
      isGoodConversion = (fConversions->At(ifc)->Prob() > probMin) && (!requireArbitratedMerged || fConversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) && (fConversions->At(ifc)->LxyCorrected((BaseVertex*)fVertex) > lxyMin);
      if( isGoodConversion == true ){
        for(unsigned int d=0; d<fConversions->At(ifc)->NDaughters(); ++d){
          const Track* trk = dynamic_cast<const ChargedParticle*> (fConversions->At(ifc)->Daughter(d))->Trk();
          if( trk ){
            const StableData* sd = dynamic_cast<const StableData*> (fConversions->At(ifc)->DaughterDat(d));
            if( sd->NWrongHits() > nWrongHitsMax ){ isGoodConversion = false; }
          } 
	  else { isGoodConversion = false; }
        }
      }
    }
    if(isGoodConversion == true){ break; }
  }
  return isGoodConversion;
}

float 
HttNtupler::computePFMuonIso(const Muon* muon, const double dRMax)
{
  const double dRMin    = 0.01;
  const double neuPtMin =  1.0;
  const double dzMax    =  0.1;
  double zLepton = (muon->BestTrk()) ? muon->BestTrk()->DzCorrected(*fVertex) : 0.0;
  float iso=0;
  for(unsigned int ipf=0; ipf<fPFCandidates->GetEntries(); ++ipf){
    const PFCandidate* pfcand = fPFCandidates->At(ipf);
    if( !pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin) ){ continue; }
    if( pfcand->TrackerTrk() && muon->TrackerTrk() && (pfcand->TrackerTrk()==muon->TrackerTrk()) ){ continue; }
    double dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(*fVertex) - zLepton) : 0;
    if( dz >= dzMax ){ continue; }
    double dr = MathUtils::DeltaR(muon->Mom(), pfcand->Mom());
    if(dr<dRMax && dr>=dRMin){ iso += pfcand->Pt(); }
  }
  return iso;
}

float 
HttNtupler::computePFElecIso(const Electron* electron, const double dRMax)
{
  const double dRMin    = 0.01;
  const double neuPtMin =  1.0;
  const double dzMax    =  0.1;
  double zLepton = (electron->BestTrk()) ? electron->BestTrk()->DzCorrected(*fVertex) : 0.0;
  float iso=0;
  for(unsigned int ipf=0; ipf<fPFCandidates->GetEntries(); ++ipf){
    const PFCandidate* pfcand = fPFCandidates->At(ipf);
    if( !pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin) ){ continue; }
    double dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(*fVertex) - zLepton) : 0;
    if( dz >= dzMax ){ continue; }
    if( pfcand->TrackerTrk() && electron->TrackerTrk() && (pfcand->TrackerTrk()==electron->TrackerTrk()) ){ continue; }
    if( pfcand->GsfTrk() && electron->GsfTrk() && (pfcand->GsfTrk()==electron->GsfTrk()) ){ continue; }
    double dr = MathUtils::DeltaR(electron->Mom(), pfcand->Mom());
    if( dr<dRMax && dr>=dRMin ){
      // eta-strip veto for photons
      if( pfcand->PFType()==PFCandidate::eGamma && fabs(electron->Eta()-pfcand->Eta()) < 0.025 ){ continue; }
      // inner cone (one tower = dR < 0.07) veto for non-photon neutrals
      if( !pfcand->HasTrk() && (pfcand->PFType()==PFCandidate::eNeutralHadron) && (MathUtils::DeltaR(electron->Mom(), pfcand->Mom()) < 0.07) ){ continue; }
      iso += pfcand->Pt();
    }
  }
  return iso;
}

bool
HttNtupler::findLeptonFootprint(const PFCandidate* pfcand)
{
  bool isLeptonFootprint = false;
  for (UInt_t q=0; q < fElectrons->GetEntries() ; ++q) {
    //consider only electrons that pass ID cuts
    if(fabs(fElectrons->At(q)->BestTrk()->D0Corrected(*fVertex)) > 0.04) continue;
    if(fabs(fElectrons->At(q)->BestTrk()->DzCorrected(*fVertex)) > 0.2)  continue;
    if(fabs(fElectrons->At(q)->SCluster()->Eta())<1.479) {
      if(fElectrons->At(q)->Pt() > 20) {
        if(fElectrons->At(q)->CoviEtaiEta()       > 0.01)  continue;
        if(fabs(fElectrons->At(q)->DeltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
        if(fabs(fElectrons->At(q)->DeltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
        if(fElectrons->At(q)->HadronicOverEm()            > 0.15)  continue;
      } else {
        if(fElectrons->At(q)->CoviEtaiEta()       > 0.012)  continue;
        if(fabs(fElectrons->At(q)->DeltaEtaSuperClusterTrackAtVtx()) > 0.007) continue;
        if(fabs(fElectrons->At(q)->DeltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
        if(fElectrons->At(q)->HadronicOverEm()            > 0.15) continue;
      }
    } else {
      if(fElectrons->At(q)->Pt() > 20) {
        if(fElectrons->At(q)->CoviEtaiEta()       > 0.03)  continue;
        if(fabs(fElectrons->At(q)->DeltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
        if(fabs(fElectrons->At(q)->DeltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
      } else {
        if(fElectrons->At(q)->CoviEtaiEta()       > 0.032)  continue;
        if(fabs(fElectrons->At(q)->DeltaEtaSuperClusterTrackAtVtx()) > 0.010) continue;
        if(fabs(fElectrons->At(q)->DeltaPhiSuperClusterTrackAtVtx()) > 0.8)  continue;
      }
    }
    //if pf candidate matches an electron passing ID cuts, then veto it
    if(pfcand->GsfTrk() && fElectrons->At(q)->GsfTrk() && pfcand->GsfTrk() == fElectrons->At(q)->GsfTrk()) isLeptonFootprint = true;
    if(pfcand->TrackerTrk() && fElectrons->At(q)->TrackerTrk() && pfcand->TrackerTrk() == fElectrons->At(q)->TrackerTrk()) isLeptonFootprint = true;
    //if pf candidate lies in veto regions of electron passing ID cuts, then veto it
    if(pfcand->BestTrk() && MathUtils::DeltaR(fElectrons->At(q)->Mom(), pfcand->Mom()) < 0.015) isLeptonFootprint = true;
    if(pfcand->PFType() == PFCandidate::eGamma && fabs(fElectrons->At(q)->Eta()) < 1.479 && fabs(fElectrons->At(q)->Eta() - pfcand->Eta()) < 0.025) isLeptonFootprint = true;
    if(pfcand->PFType() == PFCandidate::eGamma && fabs(fElectrons->At(q)->Eta()) >= 1.479 && MathUtils::DeltaR(fElectrons->At(q)->Mom(), pfcand->Mom()) < 0.08) isLeptonFootprint = true;
  }
  for (UInt_t q=0; q < fMuons->GetEntries() ; ++q) {
    // Use tracker track when available
    const Track *tmpMuTrk=0;
    if(fMuons->At(q)->HasTrackerTrk())         { tmpMuTrk = fMuons->At(q)->TrackerTrk();    }
    else if(fMuons->At(q)->HasStandaloneTrk()) { tmpMuTrk = fMuons->At(q)->StandaloneTrk(); }
    //consider only muons that pass ID cuts
    if(!(fMuons->At(q)->IsGlobalMuon() || fMuons->At(q)->IsTrackerMuon())) continue;
    if(!tmpMuTrk || tmpMuTrk->NHits() < 11 ) continue;
    if( fMuons->At(q)->Ip3dPVErr() <= 0 || fMuons->At(q)->Ip3dPVSignificance() > 4 ) continue;
    //if pf candidate matches an muon passing ID cuts, then veto it
    if(pfcand->TrackerTrk() && fMuons->At(q)->TrackerTrk() && pfcand->TrackerTrk() == fMuons->At(q)->TrackerTrk()) isLeptonFootprint = true;
    //if pf candidate lies in veto regions of muon passing ID cuts, then veto it
    if(pfcand->BestTrk() && MathUtils::DeltaR(fMuons->At(q)->Mom(), pfcand->Mom()) < 0.015) isLeptonFootprint = true;
    if(pfcand->PFType() == PFCandidate::eGamma && MathUtils::DeltaR(fMuons->At(q)->Mom(), pfcand->Mom()) < 0.05) isLeptonFootprint = true;
  }
  return isLeptonFootprint;
}

float 
HttNtupler::computePFTauIso(const PFTau* tau, const Track* track) 
{ 
  double lPtMin =  0.;
  double lDRMin =  0.;
  double lDRMax = 0.5;
  double lDZMin = 0.2;
  float  iso    = 0.0;
  for(unsigned int idx=0; idx<fPFCandidates->GetEntries(); ++idx){
    const PFCandidate* pf = fPFCandidates->At(idx);
    if( !pf->Trk() ){ continue; }
    if( pf->PFType() != PFCandidate::eHadron ){ continue; }
    if( pf->Pt() < lPtMin ){ continue; }
    if( fabs(pf->Trk()->DzCorrected(*fVertex)-track->DzCorrected(*fVertex)) < lDZMin ){ continue; }
    if( MathUtils::DeltaR(tau->Mom(), pf->Mom())<lDRMax && MathUtils::DeltaR(tau->Mom(), pf->Mom())>=lDRMin ){ iso += pf->Pt(); }
  }
  return iso;
}

float 
HttNtupler::computeOfficialPFTauIso(const PFTau* tau) 
{ 
  double lChargePt=0;  double lPuPt=0;  double lNeutPt=0;
  const Vertex* lVtx = officialTauIsoAssociatedVertex(tau);
  for(unsigned int idx=0; idx<tau->NIsoPFCandS(); ++idx){
    const PFCandidate* pCand = tau->IsoPFCand(idx);
    if( !officialTauIsoPFCandidateQuality(pCand, lVtx) ){ continue; }
    if( pCand->PFType() != PFCandidate::eHadron && pCand->PFType() != PFCandidate::eElectron && pCand->PFType() != PFCandidate::eMuon   && pCand->PFType() != PFCandidate::eGamma ){ continue; }
    pCand->Charge()==0 ? lNeutPt+=pCand->Pt() : lChargePt+=pCand->Pt();
  }
  for(unsigned int idx=0; idx<fPFCandidates->GetEntries(); ++idx){ 
    const PFCandidate* pCand = fPFCandidates->At(idx);
    if( pCand->Charge() == 0 ){ continue; }
    double pDR =  MathUtils::DeltaR(tau->Mom(), pCand->Mom());
    if( pDR > 0.8 ){ continue; }
    if( !officialTauIsoPFCandidateQuality(pCand, lVtx, true) ){ continue; }
    lPuPt += pCand->Pt();
  }
  return lChargePt + max(lNeutPt-0.4576*lPuPt,0.);
}

const Vertex*
HttNtupler::officialTauIsoAssociatedVertex(const PFTau* iTau) 
{
  const PFJet* lJet=0;
  double minDR = 100.0;
  for(unsigned int ijet=0; ijet<fPFJets->GetEntries(); ++ijet){
    const PFJet* pJet = fPFJets->At(ijet);
    double pDR = MathUtils::DeltaR(iTau->Phi(), iTau->Eta(), pJet->Phi(), pJet->Eta());
    if( pDR<minDR ){ 
      minDR = pDR;
      lJet = pJet;
    }
  }
  const Track* leadTrack=0;
  const PFCandidate* leadPF=0;
  for(unsigned int icand=0; icand<lJet->NPFCands(); ++icand){
    const PFCandidate* pCand = lJet->PFCand(icand);
    if( pCand->PFType() != PFCandidate::eHadron && pCand->PFType() != PFCandidate::eElectron && pCand->PFType() != PFCandidate::eMuon ){ continue; }
    if( !leadPF || pCand->Pt()>leadPF->Pt() ){ leadPF=pCand; leadPF->HasTrackerTrk() ? leadTrack=leadPF->TrackerTrk() : leadTrack=leadPF->GsfTrk(); }
  }
  int iId = 0;
  //if(!leadTrack) return fPrimVerts->At(iId);
  //  for(unsigned int ivtx = 0; ivtx < fPrimVerts->GetEntries(); ++ivtx) {
  //    const Vertex *pVtx = fPrimVerts->At(ivtx);
  //    iId = ivtx;
  //    break;
  //  }
  return fPrimVerts->At(iId);
}

bool 
HttNtupler::officialTauIsoPFCandidateQuality(const PFCandidate* cand, const Vertex *vtx, Bool_t invertDz, Bool_t invertGeneral) 
{
  if( cand->PFType() == PFCandidate::eHadron || cand->PFType() == PFCandidate::eElectron || cand->PFType() == PFCandidate::eMuon ){
    bool passDzCut       = false;
    bool passGeneralCuts = false;
    const Track* lTrk=0;
    if( cand->HasTrackerTrk() ) lTrk = cand->TrackerTrk();
    else if( cand->HasGsfTrk() ) lTrk = cand->GsfTrk();
    if( lTrk && fabs(lTrk->DzCorrected(*vtx))<=0.2 ){ passDzCut = true; }
    if( lTrk && lTrk->Pt()>0.5 && lTrk->RChi2()<=100. && lTrk->NPixelHits()>=0. && lTrk->NHits()>=8 && fabs(lTrk->D0Corrected(*vtx))<= 0.03 ){ passGeneralCuts = true; }
    return (passDzCut^invertDz && passGeneralCuts^invertGeneral);
  } 
  else if( cand->PFType() == PFCandidate::eNeutralHadron ){
    bool passNeutralHadronCuts = false;
    if( cand->Et()>0. ) passNeutralHadronCuts = true;
    return passNeutralHadronCuts^invertGeneral;
  }  
  else if( cand->PFType() == PFCandidate::eGamma ){
    bool passGammaCuts = false;
    if( cand->Et()>0.5 ) passGammaCuts = true;
    return passGammaCuts^invertGeneral;
  }
  return false;
}

double 
HttNtupler::computeCommonIso(const ChargedParticle* p, PFCandidateCol* pfNoPileUp, double ptMin, double extRadius, double intRadius, int isoType, double dzLep, double dzMax)
{
  double ptSum = 0.0;
  for(unsigned int i=0; i<pfNoPileUp->GetEntries(); ++i){
    const PFCandidate* pf = pfNoPileUp->At(i); assert(pf);
    PFCandidate::EPFType pfType = pf->PFType();
    if( isoType == 1 && !(pfType == PFCandidate::eHadron || pfType == PFCandidate::eElectron || pfType == PFCandidate::eMuon) ){ continue; }
    if( isoType == 2 && !(pfType == PFCandidate::eNeutralHadron) ){ continue; }
    if( isoType == 3 && !(pfType == PFCandidate::eGamma) ){ continue; }
    if(pf->Pt() >= ptMin && !(pf->TrackerTrk() && p->TrackerTrk() && pf->TrackerTrk() == p->TrackerTrk())){
      if( dzMax>0 ){
	double dz = pf->BestTrk()==0 ? 0 : fabs(pf->BestTrk()->DzCorrected(*fVertex)-dzLep);
	if( dz>dzMax && dz != 0 ){ continue; }
      }
      // add pt to running sum if the particle flow candidate is close enough in deltaR
      double dr = MathUtils::DeltaR(p->Mom(), pf->Mom());
      if(dr < extRadius && dr >= intRadius) ptSum += pf->Pt();
    }
  }
  return ptSum;
}

double 
HttNtupler::computeNaiveIso(const Particle* p, PFCandidateCol* pfPileUp, double ptMin, double extRadius, double intRadius)
{   
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
HttNtupler::separatePileup(const PFCandidateCol* pfCandidates, const Vertex* pv, const VertexCol* primVerts, PFCandidateOArr* pfPileUp, PFCandidateOArr* pfNoPileUp, bool checkClosestZVertex)
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
            vertexFound = true; closestVtx = vtx; break;
          }
          double dz = fabs(pf->SourceVertex().Z() - vtx->Z());
          if(dz < dzmin){
            closestVtx = vtx; dzmin = dz;
          }
        }
        if(checkClosestZVertex){
          // Fallback: if track is not associated with any vertex
          // associate it with the vertex, which is closest in z
          if(vertexFound || closestVtx != pv){ pfPileUp->Add(pf); }
          else{ pfNoPileUp->Add(pf); }
	}
        else{
          if(vertexFound && closestVtx != pv){ pfPileUp->Add(pf); }
          else{ pfNoPileUp->Add(pf); }
        }
      }
    }
    else{ pfNoPileUp->Add(pf); }
  }
}

bool 
HttNtupler::looseTauId(const PFTau* tau) 
{ 
  if( tau->Pt() < fPFTauPtMin                        ) return false;
  if( fabs(tau->Eta()) > fPFTauEtaMax                ) return false;
  if( !tau->DiscriminationByLooseElectronRejection() ) return false;
  if( !tau->DiscriminationByLooseMuonRejection()     ) return false;
  if( !tau->DiscriminationByDecayModeFinding()       ) return false;
  return true;
} 
bool 
HttNtupler::looseEleId(const Electron* elec) 
{ 
  if( elec->Pt () < fEleEtMin                        ) return false;
  if( elec->Pt () > fEleEtMax                        ) return false;
  if( elec->Eta() < fEleEtaMin                       ) return false;
  if( elec->Eta() > fEleEtaMax                       ) return false;
  if( !fEleTools->PassSpikeRemovalFilter(elec)       ) return false;
  if( isConversion(elec)                             ) return false;
  if( fUseGen==ESampleType::kEmbed && !elec->BestTrk()                 ) return false;
  if( elec->BestTrk()->NExpectedHitsInner() > 0      ) return false;
  return true;
}

bool 
HttNtupler::looseMuId(const Muon* muon) 
{ 
  if( muon->TrackerTrk() == 0                        ) return false;
  if( muon->BestTrk()->Pt () < fMuPtMin              ) return false;
  if( muon->BestTrk()->Eta() < fMuEtaMin             ) return false;
  if( muon->BestTrk()->Eta() > fMuEtaMax             ) return false;
  if( muon->BestTrk()->Pt () > fMuPtMax              ) return false;
  if( muon->TrackerTrk()->PtErr()/muon->Pt() > 0.1   ) return false;
  if( muon->TrackerTrk()->NPixelHits() < 1           ) return false;
  if( muon->TrackerTrk()->NHits() < 11               ) return false;
  if( muon->NValidHits()          < 1                ) return false;
  if( muon->NMatches()            < 1                ) return false;
  return true;
}

