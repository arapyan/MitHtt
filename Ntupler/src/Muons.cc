void Muons::setupTreeBranch(TTree *iTree,HLTTool *iHLT,double iPtMin,double iEtaMin) { 
  TMuon::Class()->IngnoreTObjectStreamer();
  fMuonArr = new TClonesArray("mithep::TMuon"); assert( fMuonArr );
  iTree->Branch("Muon", &fMuonArr,);
  fEtaMax = iEtaMax;
  fPtMin  = iPtMin;
  fHLT    = iHLT;
}
void Muons::fillMuons(const muonCol *iMuons, const Vertex *iPV,const PFCandidateCol* iPFNoPU,const PFCandidateCol* iPFPU,const PFCandidates *iCands) { 
  fMuonArr->Clear();
  for(UInt_t i0=0; i0 < iMuons->GetEntries(); i0++) {
    const Muon *pMu = fMuons->At(i0); 
    if(!pMu->HasTrk()) continue; 
    if((fabs(pMu->BestTrk()->Eta())                                  > fMuEtaMax)) continue;   
    if((pMu->BestTrk()->Pt()  > fMuPtMax)  || (pMu->BestTrk()->Pt()  < fMuPtMin))  continue;
    fillMuon(pMu,iPV,iPFNoPU,iPFPU);  
  }
}
bool Muons::isPFMuon(const Muon *iMu,const PFCandidates *iCands) {
  bool lPFMuon = false;
  for(int i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const PFCandidate *pCand = iCands->At(i0); 
    if(pCand-PFType() != PFCandidate::eMuon) continue;
    if(MathUtils::DeltaR(iMu->Mom(),pCand->Mom()) > 0.01) continue;
    lPFMuon = true;
  }
  return lPFMuon;
}
void Muons::fillMuon(const Muon *iMu,const Vertex *iPV,const PFCandidateCol *iPFNoPU,const PFCandidateCol *iPFPU,const PFCandidates *iCands) {
  //Make  a New Muon
  TClonesArray &rMuonArr = *fMuonArr;
  assert(rMuonArr.GetEntries() < rMuonArr.GetSize());
  const Int_t index = rMuonArr.GetEntries();  
  new(rMuonArr[index]) TMuon();
  TMuon *pMuon = (TMuon*)rMuonArr[index];
  
  //Fill
  const Track *muTrk=iMu->BestTrk();
  if(muTrk == 0) return;
  
  pMuon->pt              = muTrk->Pt();
  pMuon->ptErr           = muTrk->PtErr();
  pMuon->eta             = muTrk->Eta();
  pMuon->phi             = muTrk->Phi();
  pMuon->d0	         = muTrk->D0Corrected(*fVertex);
  pMuon->ip3d            = iMu->Ip3dPV();
  pMuon->ip3dSig         = iMu->Ip3dPVSignificance();
  pMuon->dz              = muTrk->DzCorrected(*fVertex);
  pMuon->tkNchi2         = (iMu->HasTrackerTrk()) ? iMu->TrackerTrk()->RChi2() : 0;
  pMuon->trkIso03        = iMu->IsoR03SumPt();
  pMuon->emIso03         = iMu->IsoR03EmEt();
  pMuon->hadIso03        = iMu->IsoR03HadEt();
  pMuon->pfIso03         = IsoTools::computePFMuonIso(iMu,0.3);
  pMuon->pfIso04         = IsoTools::computePFMuonIso(iMu,0.4);  
  pMuon->pfIsoCharged    = IsoTools::PFIsoNoPileUp((const ChargedParticle*)iMu, (const PFCandidateCol*) iPFNoPU,    0.0, 0.4, 0.0001, 1);
  pMuon->pfIsoNeutral    = IsoTools::PFIsoNoPileUp((const ChargedParticle*)iMu, (const PFCandidateCol*) iPFNoPU,    0.5, 0.4, 0.01  , 2);
  pMuon->pfIsoGamma      = IsoTools::PFIsoNoPileUp((const ChargedParticle*)iMu, (const PFCandidateCol*) iPFNoPU,    0.5, 0.4, 0.01  , 3);
  pMuon->puIso           = IsoTools::PFIsoNoPileUp((const ChargedParticle*)iMu, (const PFCandidateCol*) iPFPU,      0.5, 0.4, 0.01  , -1);

  if(iMu->HasGlobalTrk())          { pMuon->muNchi2 = iMu->GlobalTrk()->RChi2();     }
  else if(iMu->HasStandaloneTrk()) { pMuon->muNchi2 = iMu->StandaloneTrk()->RChi2(); }
  else if(iMu->HasTrackerTrk())    { pMuon->muNchi2 = iMu->TrackerTrk()->RChi2();    }
  
  pMuon->q               = muTrk->Charge();
  pMuon->nValidHits      = iMu->NValidHits();
  pMuon->qualityBits     = iMu->Quality().QualityMask().Mask();

  pMuon->typeBits        = 0;
  if(iMu->IsGlobalMuon())     { pMuon->typeBits |= kGlobal; }
  if(iMu->IsTrackerMuon())    { pMuon->typeBits |= kTracker; }
  if(iMu->IsStandaloneMuon()) { pMuon->typeBits |= kStandalone; }
  if(isPFMuon(iMu))           { pMuon->typeBits |= kPFMuon; }
  
  pMuon->nTkHits        = (iMu->HasTrackerTrk()) ? iMu->TrackerTrk()->NHits() : 0;
  pMuon->nPixHits       = muTrk->NPixelHits();
  pMuon->nSeg           = iMu->NSegments();
  pMuon->nMatch         = iMu->NMatches();
  pMuon->hltMatchBits   = 0;
  pMuon->hltTOMatchBits = fHLT->matchHLT(muTrk->Pt(),muTrk->Eta(),muTrk->Phi(),pMuon->hltMatchBits);
  
  pMuon->staPt  = (iMu->HasStandaloneTrk()) ? iMu->StandaloneTrk()->Pt()  : -1;  
  pMuon->staEta = (iMu->HasStandaloneTrk()) ? iMu->StandaloneTrk()->Eta() : -999;
  pMuon->staPhi = (iMu->HasStandaloneTrk()) ? iMu->StandaloneTrk()->Phi() : -999;
  pMuon->trkID  = (iMu->HasTrackerTrk()) ? iMu->TrackerTrk()->GetUniqueID() : 0;

}
