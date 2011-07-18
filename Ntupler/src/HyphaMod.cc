#include "MitHtt/Ntupler/interface/HyphaMod.hh"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
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
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TLorentzVector.h>
#include <vector>

using namespace mithep;

ClassImp(mithep::HyphaMod)

HyphaMod::HyphaMod(const char *name, const char *title):
  BaseMod        (name,title),
  fOutputFile    (0),
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
  fParticles     (0),
  fMCEvtInfo     (0),
  fMuons         (0),
  fElectrons     (0),
  fPrimVerts     (0),
  fBeamSpot      (0),
  fPFJets        (0),
  fPhotons       (0),  
  fTrigMask      (0),
  fPFMet         (0),
  fConversions   (0),
  fPileup        (0),
  fPUEnergyDensity(0),
  fPFCandidates  (0),
  fIsData        (kFALSE),
  fUseGen        (0),
  fPrintTable    (kFALSE),
  fSkipIfHLTFail (kFALSE),
  fMuPtMin       (15),
  fMuPtMax       (1000),
  fMuEtaMin      (-3),
  fMuEtaMax      (3),
  fEleEtMin      (15),
  fEleEtMax      (1000),
  fEleEtaMin     (-3),
  fEleEtaMax     (3),
  fJetPtMin      (15),
  fPhotonEtMin   (10),
  fMinNTracksFit (0),
  fMinNdof       (4),
  fMaxAbsZ       (24),
  fMaxRho        (2),
  fEventTree     (0),
  fJetCorrector  (0),
  fJetUnc       (0)
{
  // Constructor
  
  // Don't write TObject part of the objects
  TEventInfo::Class()->IgnoreTObjectStreamer();
  TGenInfo::Class()->IgnoreTObjectStreamer();
  TMuon::Class()->IgnoreTObjectStreamer();
  TElectron::Class()->IgnoreTObjectStreamer();
  TJet::Class()->IgnoreTObjectStreamer();
  TPhoton::Class()->IgnoreTObjectStreamer();
  TVertex::Class()->IgnoreTObjectStreamer();
}

//--------------------------------------------------------------------------------------------------
HyphaMod::~HyphaMod()
{
  // Destructor
}	

//--------------------------------------------------------------------------------------------------      
void HyphaMod::Begin()
{
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::SlaveBegin()
{
  //
  // Request BAMBU branches
  //
  ReqBranch(fPartName,            fParticles); 
  ReqBranch(fMCEvtInfoName,       fMCEvtInfo);
  ReqBranch(fMuonName,            fMuons);
  ReqBranch(fElectronName,        fElectrons);
  ReqBranch(fPrimVtxName,         fPrimVerts);
  ReqBranch(fBeamSpotName,        fBeamSpot);
  ReqBranch(fPFJetName,           fPFJets);
  ReqBranch(fTrigMaskName,        fTrigMask);
  ReqBranch(fPFMetName,           fPFMet);
  ReqBranch(fPhotonName,          fPhotons);
  ReqBranch(fConversionName,      fConversions);
  ReqBranch(fPileupName,          fPileup);
  ReqBranch(fPUEnergyDensityName, fPUEnergyDensity);
  ReqBranch(fPFCandidateName,     fPFCandidates);  
  
  //
  // Set up arrays
  //
  fMuonArr     = new TClonesArray("mithep::TMuon");	assert(fMuonArr);
  fElectronArr = new TClonesArray("mithep::TElectron"); assert(fElectronArr);
  fPFJetArr    = new TClonesArray("mithep::TJet");	assert(fPFJetArr);
  fPhotonArr   = new TClonesArray("mithep::TPhoton");	assert(fPhotonArr);
  fPVArr       = new TClonesArray("mithep::TVertex");	assert(fPVArr);
  
  //
  // Create output file
  //
  fOutputFile = new TFile(fOutputName, "RECREATE");
  
  //
  // Initialize data trees and structs 
  // 
  fEventTree = new TTree("Events","Events");

  fEventTree->Branch("Info",&fEventInfo);
  if(!fIsData && fUseGen)
    fEventTree->Branch("Gen",&fGenInfo);
  
  fEventTree->Branch("Muon",    &fMuonArr);
  fEventTree->Branch("Electron",&fElectronArr);
  fEventTree->Branch("PFJet",   &fPFJetArr);
  fEventTree->Branch("Photon",  &fPhotonArr);
  fEventTree->Branch("PV",      &fPVArr);

  // fTauTree = new TTree("Taus","Taus");
  // fTauTree->Branch("tau1",&fTau1);
  // fTauTree->Branch("tau2",&fTau2);
  
  //
  // Set up jet corrections for PF jets
  //
  std::vector<JetCorrectorParameters> correctionParameters;
  for(UInt_t icorr=0; icorr<fJetCorrParsv.size(); icorr++)
    correctionParameters.push_back(JetCorrectorParameters(fJetCorrParsv[icorr].Data()));
    
  // initialize jet corrector class
  fJetCorrector = new FactorizedJetCorrector(correctionParameters);

  // initialize jet uncertainties
  JetCorrectorParameters param(string("/home/dkralph/cms/cmssw/022/CMSSW_4_2_4_patch1/src/MitPhysics/data/START42_V12_AK5PF_Uncertainty.txt"));
  fJetUnc = new JetCorrectionUncertainty(param);

  // setup selecting with JSON file, if necessary
  fhasJSON = fJSONv.size() > 0;
  for(UInt_t i=0; i<fJSONv.size(); i++) {
    frlrm.AddJSONFile(fJSONv[i].Data());
  }

  fLastRunLumi = RunLumiRangeMap::RunLumiPairType(0,0);
}

// instantiate the jec uncertainty object param are the jet corrector parameters
// that you usually get form the jes correction file (in analogy to the factorized jet corrector)
// scaledJet in this case was a reco::Jet, you just multiply the p4 I guess
// scaledJet.scaleEnergy( 1+jetMet);


// The JetUncertainties class itself:
// http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h?revision=1.5&view=markup




//--------------------------------------------------------------------------------------------------
void HyphaMod::SlaveTerminate()
{
  //
  // Save to ROOT file
  //
  fEventTree->Print();

  fOutputFile->Write();
  fOutputFile->Close();
  
  delete fMuonArr;
  delete fElectronArr;
  delete fPFJetArr;
  delete fPhotonArr;
  delete fPVArr;
  
  delete fJetCorrector;
  delete fJetUnc;

  //
  // dump json file
  //
  TString jsonfname = fOutputName.ReplaceAll("root","json");
  if(fIsData)
    fRunLumiSet.DumpJSONFile(jsonfname.Data());
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::Terminate()
{
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::BeginRun()
{
  if(HasHLTInfo() && fPrintTable) { GetHLTTable()->Print(); }
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::EndRun()
{
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::Process()
{
  RunLumiRangeMap::RunLumiPairType rl(GetEventHeader()->RunNum(), GetEventHeader()->LumiSec());
  if(fhasJSON && !frlrm.HasRunLumi(rl)) return;  // not certified run? Skip to next event...
  if(rl!=fLastRunLumi) {
    fLastRunLumi = rl;
    fRunLumiSet.Add(rl);
  }

  //
  // Load branches
  //
  if(!fIsData && fUseGen) LoadBranch(fPartName);
  if(!fIsData && fUseGen) LoadBranch(fMCEvtInfoName);
  LoadBranch(fMuonName);
  LoadBranch(fElectronName);
  LoadBranch(fPrimVtxName);
  LoadBranch(fBeamSpotName);
  LoadBranch(fPFJetName);
  LoadBranch(fTrigMaskName);
  LoadBranch(fPFMetName); 
  LoadBranch(fPhotonName);
  LoadBranch(fConversionName);
  LoadBranch(fPUEnergyDensityName);
  LoadBranch(fPFCandidateName);
  if(!fIsData)
    LoadBranch(fPileupName);

  //
  // Scan generator info
  //
  if(fUseGen) {
    if(fUseGen==1) FillGenH();
    if(fUseGen==2) FillGenZ();
    if(fUseGen==3) FillGenW();
    if(fUseGen==4) FillGenWW();
  }
  
  //
  // Get HLT info. Trigger objects can be matched by name to the corresponding trigger that passed.
  // note: TriggerName::Id() is bambu numbering scheme, fTriggerIdsv[itrig] is kevin's
  //
  ULong_t trigbits=0;
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      // if event passed this trigger, set the trigger's bit in trigbits
      if(fTrigMask->At(trigname->Id())) { trigbits |= fTriggerIdsv[itrig]; }
    }  
  }
  if(fSkipIfHLTFail && (trigbits==0))
    return;
  
  IncNEventsProcessed();
  
  fMuonArr->Clear();
  fElectronArr->Clear();
  fPFJetArr->Clear();
  fPhotonArr->Clear();
  fPVArr->Clear();


  //
  // Get beam spot. If no beam spot information is available, default the coordinates to 99999
  //
  Double_t bsx=99999, bsy=99999, bsz=99999;
  if(fBeamSpot) {
    if(fBeamSpot->GetEntries() > 1) 
      std::cout << "********** More than 1 beam spot! **********" << std::endl;
    const BeamSpot *bs = fBeamSpot->At(0);
    bsx = bs->X();
    bsy = bs->Y();
    bsz = bs->Z();
  }

  //
  // Get primary vertices
  // Assumes primary vertices are ordered by sum-pT^2 (as should be in CMSSW)
  // NOTE: if no PV is found from fitting tracks, the beamspot is used
  //
  const Vertex *bestPV = 0;
  Bool_t hasGoodPV = kFALSE;  
  for(UInt_t i=0; i<fPrimVerts->GetEntries(); ++i) {
    const Vertex *pv = fPrimVerts->At(i);
    
    // Select best PV for corrected d0; if no PV passing cuts, the first PV in the collection will be used
    //if(!pv->IsValid()) continue;
    if(pv->NTracksFit()     < fMinNTracksFit) continue;
    if(pv->Ndof()	    < fMinNdof)	      continue;
    if(fabs(pv->Z())	    > fMaxAbsZ)	      continue;
    if(pv->Position().Rho() > fMaxRho)	      continue;    
    hasGoodPV = kTRUE;
    
    FillPV(pv);
    
    if(!bestPV) bestPV = pv;
  }
  if(!bestPV) bestPV = fPrimVerts->At(0);
  fVertex.SetPosition(bestPV->X(),bestPV->Y(),bestPV->Z());
  fVertex.SetErrors(bestPV->XErr(),bestPV->YErr(),bestPV->ZErr());
  

  //
  // Loop through muons (and general tracks if desired).
  //
  vector<const Muon*> muonv;    // array of pointers to preselected muons ... 
  
  assert(fMuons);
  for(UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *mu = fMuons->At(i); 
    if(!mu->HasTrk()) continue; 
    
    // Use tracker tracks for kinematics when available
    const Track *muTrk=0;
    if(mu->HasTrackerTrk())     { muTrk = mu->TrackerTrk(); }
    else if(mu->HasGlobalTrk()) { muTrk = mu->GlobalTrk(); }
    else                        { muTrk = mu->StandaloneTrk(); } 
              
    if((muTrk->Eta() < fMuEtaMin) || (muTrk->Eta() > fMuEtaMax)) continue;   
    if((muTrk->Pt()  > fMuPtMax)  || (muTrk->Pt()  < fMuPtMin))  continue;
    
    FillMuon(mu);  
  }
       

  //
  // Loop through electrons.
  //
  vector<const Electron*> elev;  // array of pointers to preselected electrons ... 
  ElectronTools eleTools;        // helper class for electron ID decisions
  
  assert(fElectrons);                
  for(UInt_t i=0; i<fElectrons->GetEntries(); ++i) {
    const Electron *ele = fElectrons->At(i);  
   
    if((ele->Pt()  < fEleEtMin)  || (ele->Pt()  > fEleEtMax))  continue;  // electron pT cut
    if((ele->Eta() < fEleEtaMin) || (ele->Eta() > fEleEtaMax)) continue;  // electron eta cut
    if(!eleTools.PassSpikeRemovalFilter(ele))                  continue;  // spike cleaning
        
    FillElectron(ele);  // fill electron data object
  }
  
  //
  // Loop through jets
  //
  assert(fPFJets);
  for(UInt_t i=0; i<fPFJets->GetEntries(); ++i) {
    const PFJet *jet = fPFJets->At(i);

    const FourVectorM rawMom = jet->RawMom();

    fJetCorrector->setJetEta(rawMom.Eta());
    fJetCorrector->setJetPt(rawMom.Pt());
    fJetCorrector->setJetPhi(rawMom.Phi());
    fJetCorrector->setJetE(rawMom.E());
    fJetCorrector->setRho(fPUEnergyDensity->At(0)->RhoHighEta());
    fJetCorrector->setJetA(jet->JetArea());
    fJetCorrector->setJetEMF(-99.0);     
    
    // keep all jets above specified threshold (after energy correction)
    // and all jets with valid b-tag value (Track Counting High Efficiency method default is -100)
    Double_t correction = fJetCorrector->getCorrection();
    Double_t pt = rawMom.Pt()*correction;
    if(pt > fJetPtMin || (jet->TrackCountingHighEffBJetTagsDisc() != -100)) {

      if(jet->E()==0) continue;
      if(jet->ChargedHadronEnergy()/jet->E()  <=  0)	continue;  //   'chargedHadronEnergyFraction > 0.0 &'	  
      if(jet->NeutralHadronEnergy()/jet->E()   >  0.99)	continue;  //	'neutralHadronEnergyFraction < 0.99 &'	  
      if(jet->ChargedEmEnergy()/jet->E()       >  0.99)	continue;  //	'chargedEmEnergyFraction < 0.99 &'	  
      if(jet->NeutralEmEnergy()/jet->E()       >  0.99)	continue;  //	'neutralEmEnergyFraction < 0.99 &'	  
      if(jet->ChargedMultiplicity()           ==  0)	continue;  //	'chargedMultiplicity > 0 &'		  
      if(jet->NConstituents()                  <  2)	continue;  //	'nConstituents > 1'

      FillJet(jet);
    }
  }
    
  //
  // Loop through photons
  //
  assert(fPhotons);
  for(UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    const Photon *pho= fPhotons->At(i);
    if(pho->SCluster()->Et() > fPhotonEtMin) { FillPhoton(pho); }
  } 
  
  //
  // Compute MET
  //
  TLorentzVector pfmet; pfmet.SetPxPyPzE(fPFMet->At(0)->Mex(),fPFMet->At(0)->Mey(),0,0);
  Double_t trkMetx=0, trkMety=0, trkSumET=0;
   
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const Double_t trkDzCut  = 0.1;
    
    const PFCandidate *pfcand = fPFCandidates->At(i);
    if( (pfcand->HasTrackerTrk() && (fabs(pfcand->TrackerTrk()->DzCorrected(fVertex))<trkDzCut)) ||
        (pfcand->HasGsfTrk()     && (fabs(pfcand->GsfTrk()->DzCorrected(fVertex))<trkDzCut)) ) {
      
      trkMetx  -= pfcand->Px();
      trkMety  -= pfcand->Py();
      trkSumET += pfcand->Pt();
    }
  }
  TLorentzVector trkmet; trkmet.SetPxPyPzE(trkMetx,trkMety,0,0);

  Int_t npu = -1;
  if(!fIsData) {
    for(UInt_t i=0;i<fPileup->GetEntries();i++) {
      if(fPileup->At(i)->GetBunchCrossing() == 0)
	npu = fPileup->At(i)->GetPU_NumInteractions();
    }
    assert(npu>=0);
  }  
  //
  // Fill event info tree
  //    
  fEventInfo.runNum       = GetEventHeader()->RunNum();
  fEventInfo.evtNum       = GetEventHeader()->EvtNum();
  fEventInfo.lumiSec      = GetEventHeader()->LumiSec();
  fEventInfo.nPU          = fIsData ? 0 : npu;
  fEventInfo.triggerBits  = trigbits;
  fEventInfo.pvx          = fVertex.X();
  fEventInfo.pvy          = fVertex.Y();
  fEventInfo.pvz          = fVertex.Z();
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
  
  // Fill the tree
  fEventTree->Fill();
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::FillMuon(const Muon *mu)
{
  assert(mu);

  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];
  
  // Use tracker track when available
  const Track *muTrk=0;
  if(mu->HasTrackerTrk())         { muTrk = mu->TrackerTrk(); }
  else if(mu->HasGlobalTrk())     { muTrk = mu->GlobalTrk(); }
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
  pMuon->d0       = muTrk->D0Corrected(fVertex);
  pMuon->dz       = muTrk->DzCorrected(fVertex);
  pMuon->tkNchi2  = (mu->HasTrackerTrk()) ? mu->TrackerTrk()->RChi2() : 0;
  
  if(mu->HasGlobalTrk())          { pMuon->muNchi2 = mu->GlobalTrk()->RChi2();     }
  else if(mu->HasStandaloneTrk()) { pMuon->muNchi2 = mu->StandaloneTrk()->RChi2(); }
  else if(mu->HasTrackerTrk())    { pMuon->muNchi2 = mu->TrackerTrk()->RChi2();    }
  
  pMuon->q          = muTrk->Charge();
  pMuon->nValidHits = mu->NValidHits();
  
  pMuon->qualityBits = mu->Quality().QualityMask().Mask();

  //
  // NOTE:
  // It is possible for a muon to be TK+SA. The muon reco associates a TK with a SA if
  // chamber matches for the TK and hits for the SA share DetIDs
  //   (see hypernews thread: https://hypernews.cern.ch/HyperNews/CMS/get/csa08-muons/57/2/1/1/1.html)
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

  pMuon->pfPx=0;
  pMuon->pfPy=0;
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {    
    const PFCandidate *pfcand = fPFCandidates->At(i);
    if(mu->HasTrackerTrk() && mu->TrackerTrk() == pfcand->TrackerTrk()) {
      pMuon->pfPx = pfcand->Px();
      pMuon->pfPy = pfcand->Py();
      break;
    }
  }
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::FillElectron(const Electron *ele)
{
  assert(ele);
 
  ElectronTools eleTools;
  
  TClonesArray &rElectronArr = *fElectronArr;
  assert(rElectronArr.GetEntries() < rElectronArr.GetSize());
  const Int_t index = rElectronArr.GetEntries();  
  new(rElectronArr[index]) TElectron();
  TElectron *pElectron = (TElectron*)rElectronArr[index];
                    
  pElectron->pt              = ele->Pt();
  pElectron->eta             = ele->Eta();
  pElectron->phi             = ele->Phi();
  pElectron->trkIso03        = ele->TrackIsolationDr03();
  pElectron->emIso03         = ele->EcalRecHitIsoDr03();
  pElectron->hadIso03        = ele->HcalTowerSumEtDr03();
  pElectron->pfIso03         = computePFEleIso(ele,0.3); 
  pElectron->pfIso04         = computePFEleIso(ele,0.4);
  pElectron->d0              = ele->BestTrk()->D0Corrected(fVertex);
  pElectron->dz              = ele->BestTrk()->DzCorrected(fVertex);  
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
  
  pElectron->hltMatchBits    = MatchHLT(ele->SCluster()->Eta(),ele->SCluster()->Phi());
  pElectron->scID            = ele->SCluster()->GetUniqueID();
  pElectron->trkID           = (ele->HasTrackerTrk()) ? ele->TrackerTrk()->GetUniqueID() : 0;
  pElectron->isEcalDriven    = ele->IsEcalDriven();
  pElectron->isConv          = IsConversion(ele);
  
  pElectron->pfPx=0;
  pElectron->pfPy=0;
  assert(fPFCandidates);
  for(UInt_t i=0; i<fPFCandidates->GetEntries(); ++i) {
    const PFCandidate *pfcand = fPFCandidates->At(i);
 
    if( (pfcand->HasTrackerTrk() && ele->TrackerTrk() == pfcand->TrackerTrk()) ||
        (pfcand->HasGsfTrk() && ele->GsfTrk() == pfcand->GsfTrk()) ) {
          pElectron->pfPx = pfcand->Px();
          pElectron->pfPy = pfcand->Py();
          break;
    }	 
  }
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::FillJet(const PFJet *jet)
{
  TClonesArray &rPFJetArr = *fPFJetArr;
  assert(rPFJetArr.GetEntries() < rPFJetArr.GetSize());
  const Int_t index = rPFJetArr.GetEntries();  
  new(rPFJetArr[index]) TJet();
  TJet *pPFJet = (TJet*)rPFJetArr[index]; 
  
  const FourVectorM rawMom = jet->RawMom();
  fJetCorrector->setJetEta(rawMom.Eta());
  fJetCorrector->setJetPt(rawMom.Pt());
  fJetCorrector->setJetPhi(rawMom.Phi());
  fJetCorrector->setJetE(rawMom.E());
  fJetCorrector->setRho(fPUEnergyDensity->At(0)->RhoHighEta());
  fJetCorrector->setJetA(jet->JetArea());
  fJetCorrector->setJetEMF(-99.0);

  fJetUnc->setJetPt(rawMom.Pt());
  fJetUnc->setJetEta(rawMom.Eta());
  // this is an up shift, down shift have false as argument
  Float_t hierr  = fJetUnc->getUncertainty(true);
  // fJetUnc->setJetPt(rawMom.Pt());
  // fJetUnc->setJetEta(rawMom.Eta());
  // Float_t lowerr = fJetUnc->getUncertainty(false);

  pPFJet->pt   = (rawMom.Pt())*(fJetCorrector->getCorrection());
  pPFJet->eta  = rawMom.Eta();
  pPFJet->phi  = rawMom.Phi();
  pPFJet->mass = jet->Mass();
  pPFJet->unc  = hierr;
  pPFJet->area = jet->JetArea();
  
  pPFJet->tche = jet->TrackCountingHighEffBJetTagsDisc();
  pPFJet->tchp = jet->TrackCountingHighPurBJetTagsDisc();
  
  pPFJet->mcFlavor = jet->MatchedMCFlavor();

  pPFJet->hltMatchBits = MatchHLT(jet->Eta(),jet->Phi());
}

      
//--------------------------------------------------------------------------------------------------
void HyphaMod::FillPhoton(const Photon *pho)
{
  TClonesArray &rPhotonArr = *fPhotonArr;
  assert(rPhotonArr.GetEntries() < rPhotonArr.GetSize());
  const Int_t index = rPhotonArr.GetEntries();  
  new(rPhotonArr[index]) TPhoton();
  TPhoton *pPhoton = (TPhoton*)rPhotonArr[index];
  
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

//--------------------------------------------------------------------------------------------------
void HyphaMod::FillPV(const Vertex *pv) 
{
  TClonesArray &rPVArr = *fPVArr;
  assert(rPVArr.GetEntries() < rPVArr.GetSize());
  const Int_t index = rPVArr.GetEntries();  
  new(rPVArr[index]) TVertex();
  TVertex *pVertex = (TVertex*)rPVArr[index];
  
  pVertex->nTracksFit = pv->NTracksFit();
  pVertex->ndof       = pv->Ndof();      
  pVertex->chi2       = pv->Chi2();  
  pVertex->x          = pv->X();
  pVertex->y          = pv->Y();
  pVertex->z          = pv->Z();
  
  pVertex->sumPt=0;
  for(UInt_t itrk=0; itrk<pv->NTracks(); itrk++)
    pVertex->sumPt += pv->Trk(itrk)->Pt();				   
}
      
      
//--------------------------------------------------------------------------------------------------
ULong_t HyphaMod::MatchHLT(const Double_t eta, const Double_t phi)
{
  ULong_t bits = 0;
  
  const Double_t hltMatchR = 0.2;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
   
      while(to) {             
        if(to->IsHLT()) {          
	  
	  if(fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
	  }
	  
	  if(fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt2v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds2v[itrig];
	  }
	  
	  if(fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
          }
        } 
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}

ULong_t HyphaMod::MatchHLT(const Double_t pt, const Double_t eta, const Double_t phi)
{
  ULong_t bits = 0;
  
  const Double_t hltMatchR = 0.2;
  const Double_t hltMatchPtFrac = 1;
  
  if(HasHLTInfo()) {
    const TriggerTable *hltTable = GetHLTTable();
    assert(hltTable);
    for(UInt_t itrig=0; itrig<fTriggerNamesv.size(); itrig++) {
      const TriggerName *trigname = hltTable->Get(fTriggerNamesv[itrig].Data());
      if(!trigname) continue;
      
      const TList *list = GetHLTObjectsTable()->GetList(fTriggerNamesv[itrig].Data());
      if(!list) continue;
      TIter iter(list->MakeIterator());
      const TriggerObject *to = dynamic_cast<const TriggerObject*>(iter.Next()); 
    
      while(to) {         
        if(to->IsHLT()) {
	  
	  if(fTriggerObjNames1v[itrig].Length()>0 && fTriggerObjNames1v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
	  }
	  
	  if(fTriggerObjNames2v[itrig].Length()>0 && fTriggerObjNames2v[itrig].CompareTo(to->ModuleName())==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt2v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds2v[itrig];
	  }
	  
	  if(fTriggerObjNames1v[itrig].Length()==0 && fTriggerObjNames2v[itrig].Length()==0) {
	    Bool_t match = kTRUE;
	    if(to->Pt() < fTriggerObjMinPt1v[itrig])                       match=kFALSE;  // minimum pT threshold on trigger object
	    if(MathUtils::DeltaR(phi,eta,to->Phi(),to->Eta()) > hltMatchR) match=kFALSE;  // eta-phi matching
	    if(fabs(pt - to->Pt())>hltMatchPtFrac*(to->Pt()))              match=kFALSE;  // pT matching
	    if(match) bits |= fTriggerObjIds1v[itrig];
          }
        }
        to = dynamic_cast<const TriggerObject*>(iter.Next());
      }    
    }
  }
  
  return bits;
}
      
//--------------------------------------------------------------------------------------------------
Bool_t HyphaMod::IsConversion(const Electron *ele) 
{
  Bool_t isGoodConversion = kFALSE;
  
  const UInt_t   nWrongHitsMax = 0;
  const Double_t probMin       = 1e-6;
  const Double_t lxyMin        = 2.0;
  const Bool_t   matchCkf      = kTRUE;
  const Bool_t   requireArbitratedMerged = kFALSE;
  
  for (UInt_t ifc=0; ifc<fConversions->GetEntries(); ifc++) {
    Bool_t ConversionMatchFound = kFALSE;
    for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
      const Track *trk = dynamic_cast<const ChargedParticle*>
        (fConversions->At(ifc)->Daughter(d))->Trk();
      if(ele->GsfTrk() == trk || (matchCkf && ele->TrackerTrk()==trk)) {
        ConversionMatchFound = kTRUE;
        break;
      }
    }

    // if match between the e-track and one of the conversion legs
    if (ConversionMatchFound == kTRUE){
      isGoodConversion =  (fConversions->At(ifc)->Prob() > probMin) &&
        (!requireArbitratedMerged || fConversions->At(ifc)->Quality().Quality(ConversionQuality::arbitratedMerged)) &&
        (fConversions->At(ifc)->LxyCorrected(&fVertex) > lxyMin);

      if (isGoodConversion == kTRUE) {
        for (UInt_t d=0; d<fConversions->At(ifc)->NDaughters(); d++) {
          const Track *trk = dynamic_cast<const ChargedParticle*>
            (fConversions->At(ifc)->Daughter(d))->Trk();
          if (trk) {
            const StableData *sd = dynamic_cast<const StableData*>
              (fConversions->At(ifc)->DaughterDat(d));
            if (sd->NWrongHits() > nWrongHitsMax)
              isGoodConversion = kFALSE;
          } else {
            isGoodConversion = kFALSE;
          }
        }
      }
    }

    if(isGoodConversion == kTRUE) break;
    
  } // loop over all fConversions 

  return isGoodConversion;
}

//--------------------------------------------------------------------------------------------------
Float_t HyphaMod::computePFMuonIso(const Muon *muon, const Double_t dRMax)
{
  const Double_t dRMin    = 0;
  const Double_t neuPtMin = 1.0;
  const Double_t dzMax    = 0.1;
    
  Double_t zLepton = (muon->BestTrk()) ? muon->BestTrk()->DzCorrected(fVertex) : 0.0;
  
  Float_t iso=0;
  for(UInt_t ipf=0; ipf<fPFCandidates->GetEntries(); ipf++) {
    const PFCandidate *pfcand = fPFCandidates->At(ipf);
    
    if(!pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin)) continue;  // pT cut on neutral particles
    
    // exclude THE muon
    if(pfcand->TrackerTrk() && muon->TrackerTrk() && (pfcand->TrackerTrk()==muon->TrackerTrk())) continue;
    
    // dz cut
    Double_t dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(fVertex) - zLepton) : 0;
    if(dz >= dzMax) continue;
    
    // check iso cone
    Double_t dr = MathUtils::DeltaR(muon->Mom(), pfcand->Mom());
    if(dr<dRMax && dr>=dRMin)
      iso += pfcand->Pt(); 
  }
  
  return iso;
}

Float_t HyphaMod::computePFEleIso(const Electron *electron, const Double_t dRMax)
{
  const Double_t dRMin    = 0;
  const Double_t neuPtMin = 1.0;
  const Double_t dzMax    = 0.1;
    
  Double_t zLepton = (electron->BestTrk()) ? electron->BestTrk()->DzCorrected(fVertex) : 0.0;
  
  Float_t iso=0;
  for(UInt_t ipf=0; ipf<fPFCandidates->GetEntries(); ipf++) {
    const PFCandidate *pfcand = fPFCandidates->At(ipf);
    
    if(!pfcand->HasTrk() && (pfcand->Pt()<=neuPtMin)) continue;  // pT cut on neutral particles
    
    // dz cut
    Double_t dz = (pfcand->BestTrk()) ? fabs(pfcand->BestTrk()->DzCorrected(fVertex) - zLepton) : 0;
    if(dz >= dzMax) continue;
    
    // remove THE electron
    if(pfcand->TrackerTrk() && electron->TrackerTrk() && (pfcand->TrackerTrk()==electron->TrackerTrk())) continue;
    if(pfcand->GsfTrk()     && electron->GsfTrk()     && (pfcand->GsfTrk()==electron->GsfTrk()))         continue;
    
    // check iso cone
    Double_t dr = MathUtils::DeltaR(electron->Mom(), pfcand->Mom());
    if(dr<dRMax && dr>=dRMin) {
      // eta-strip veto for photons
      if((pfcand->PFType() == PFCandidate::eGamma) && fabs(electron->Eta() - pfcand->Eta()) < 0.025) continue;
      
      // Inner cone (one tower = dR < 0.07) veto for non-photon neutrals
      if(!pfcand->HasTrk() && (pfcand->PFType() == PFCandidate::eNeutralHadron) && 
         (MathUtils::DeltaR(electron->Mom(), pfcand->Mom()) < 0.07)) continue;
      
      iso += pfcand->Pt();
    }
  }
  
  return iso;
}
      
//--------------------------------------------------------------------------------------------------
void HyphaMod::FillGenH() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
 
  const MCParticle *boson=0, *dau1=0, *dau2=0, *tau1=0, *tau2=0;
  
  Int_t id1=0, id2=0;

  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
    if( (p->Status()==3) && ((p->PdgId()==25) || (p->PdgId()==35) || (p->PdgId()==36)) ) {
      boson = p;
      
      // loop through daughters and look for taus
      for(UInt_t ii=0; ii<boson->NDaughters(); ii++) {
        const MCParticle *tau = boson->Daughter(ii); 
	if(abs(tau->PdgId())==15) {  
	
	  while(tau->HasDaughter(tau->PdgId(),kTRUE) && (tau->Status()!=1))
	    tau = tau->FindDaughter(tau->PdgId(),kTRUE);	  

	  if(tau->PdgId()>0) tau1 = tau;
	  if(tau->PdgId()<0) tau2 = tau;
	  
          // Loop through daughters and look for leptons  
          for(UInt_t j=0; j<tau->NDaughters(); j++) {
            const MCParticle *d = tau->Daughter(j);
	    if(d->PdgId()==22) continue; // skip photons
            if((abs(d->PdgId())==11) || (abs(d->PdgId())==13)) {
          
	      // traverse down daughter muon tree
              while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
                d = d->FindDaughter(d->PdgId(),kTRUE);	  
	      
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
	      } else {
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
    for(UInt_t j=0; j<boson->NDaughters(); j++) {
      const MCParticle *d = boson->Daughter(j);
      if(d->PdgId()==boson->PdgId()) continue;
      if(d->PdgId()>0) dau1 = d;
      if(d->PdgId()<0) dau2 = d;
    }
    cout << "Weird error: " << boson->PdgId() << " --> " << dau1->PdgId() << " + " << dau2->PdgId() << endl;
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

  // Bool_t filltau = (tau1!=0) && (tau2!=0);
  // if(filltau) {
  //   fTau1.SetPtEtaPhiM(tau1->Pt(),tau1->Eta(),tau1->Phi(),tau1->Mass());
  //   fTau2.SetPtEtaPhiM(tau2->Pt(),tau2->Eta(),tau2->Phi(),tau2->Mass());
  //   fTauTree->Fill();
  // }
    
}

//--------------------------------------------------------------------------------------------------
void HyphaMod::FillGenZ() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0, *tau1=0, *tau2=0;
  
  Int_t id1=0, id2=0;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
      
    if( (p->PdgId() == 23) && (p->Status() == 3) ) {
      boson = p;

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<boson->NDaughters(); j++) {
        const MCParticle *tau = boson->Daughter(j);
	if(abs(tau->PdgId())==15) {
          while(tau->HasDaughter(tau->PdgId(),kTRUE) && (tau->Status()!=1))
            tau = tau->FindDaughter(tau->PdgId(),kTRUE);
	  if(tau->PdgId()>0) tau1 = tau;
	  if(tau->PdgId()<0) tau2 = tau;
	}
	
        const MCParticle *d = boson->Daughter(j);
        if((abs(d->PdgId())==11) || (abs(d->PdgId())==13) || (abs(d->PdgId())==15)) {
          
          // traverse down daughter lepton tree
          while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
            d = d->FindDaughter(d->PdgId(),kTRUE);	  

          if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	  if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	  if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	  if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	  if(d->PdgId()== 15) id1 =  EGenType::kTau;
	  if(d->PdgId()==-15) id2 = -EGenType::kTau;
	  
	  if(abs(d->PdgId())==15) {
	    if(d->HasDaughter(11)) { 
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
  }
  
  assert(boson);
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

  // Bool_t filltau= (tau1!=0) && (tau2!=0);
  // if(filltau) {
  //   fTau1.SetPtEtaPhiM(tau1->Pt(),tau1->Eta(),tau1->Phi(),tau1->Mass());
  //   fTau2.SetPtEtaPhiM(tau2->Pt(),tau2->Eta(),tau2->Phi(),tau2->Mass());
  //   fTauTree->Fill();
  // }

}

//--------------------------------------------------------------------------------------------------
void HyphaMod::FillGenW() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id1=0, id2=0;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
      
    if( (abs(p->PdgId()) == 24) && (p->Status() == 3) ) {
      boson = p;

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<boson->NDaughters(); j++) {
        const MCParticle *d = boson->Daughter(j);
        if(d->PdgId() == boson->PdgId()) continue;  
        
	// traverse down daughter lepton tree
        while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
          d = d->FindDaughter(d->PdgId(),kTRUE);	  

        if(d->PdgId()== 11) id1 =  EGenType::kElectron;
	if(d->PdgId()==-11) id2 = -EGenType::kElectron;
	if(d->PdgId()== 13) id1 =  EGenType::kMuon;
	if(d->PdgId()==-13) id2 = -EGenType::kMuon;
	if(d->PdgId()== 15) id1 =  EGenType::kTau;
	if(d->PdgId()==-15) id2 = -EGenType::kTau;
	
	if(abs(d->PdgId())==15) {
	  if(d->HasDaughter(11)) { 
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
  
  assert(boson);
  if(!dau1 || !dau2) {
    while(boson->HasDaughter(boson->PdgId(),kTRUE))
      boson = boson->FindDaughter(boson->PdgId(),kTRUE);
    
    // Loop through daughters and look for leptons  
    for(UInt_t j=0; j<boson->NDaughters(); j++) {
      const MCParticle *d = boson->Daughter(j);
      if(d->PdgId() == boson->PdgId()) continue;  
      
      // traverse down daughter lepton tree
      while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
    	d = d->FindDaughter(d->PdgId(),kTRUE);  	

      if(d->PdgId()== 11) id1 =  EGenType::kElectron;
      if(d->PdgId()==-11) id2 = -EGenType::kElectron;
      if(d->PdgId()== 13) id1 =  EGenType::kMuon;
      if(d->PdgId()==-13) id2 = -EGenType::kMuon;
      if(d->PdgId()== 15) id1 =  EGenType::kTau;
      if(d->PdgId()==-15) id2 = -EGenType::kTau;

      if(abs(d->PdgId())==15) {
        if(d->HasDaughter(11)) { 
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

//--------------------------------------------------------------------------------------------------
void HyphaMod::FillGenWW() 
{
  assert(fParticles);
  assert(fMCEvtInfo);
  
  const MCParticle *boson=0, *dau1=0, *dau2=0;
  
  Int_t id=0, id1=0, id2=0;
  Int_t nW=0;
  
  for(UInt_t i=0; i<fParticles->GetEntries(); ++i) {
    const MCParticle *p = fParticles->At(i);
      
    //--------------- PYTHIA FSR mode ---------------//
    // a "branching" in the process tree is created
    // for every physical process; need to scan down
    // the lepton branches and search for photons
    //
    if( (p->PdgId() == 23) && (p->Status() == 3) ) boson = p;
    
    if( (abs(p->PdgId()) == 24) && (p->Status() == 3) ) {
      boson = p;
      nW++;
      if(nW==2) id = EGenType::kWW;

      // Loop through daughters and look for leptons  
      for(UInt_t j=0; j<boson->NDaughters(); j++) {
        const MCParticle *d = boson->Daughter(j);
        if(d->PdgId() == boson->PdgId()) continue;  
        
	// traverse down daughter lepton tree
        while(d->HasDaughter(d->PdgId(),kTRUE) && (d->Status()!=1))
          d = d->FindDaughter(d->PdgId(),kTRUE);	  

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
	
	if(abs(d->PdgId())==15) {
	  if(d->HasDaughter(11)) { 
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
  
  FourVectorM vDilep;
  if(dau1 && dau2)
    vDilep = dau1->Mom() + dau2->Mom();
  
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
