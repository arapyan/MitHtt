#include "iostream"

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "MitHtt/Ntupler/interface/MetSignificance.h"
#include "RecoMET/METAlgorithms/interface/significanceAlgo.h"

using namespace mithep;

MetSignificance::MetSignificance()  
{
  loadResolutions();
  fDRLepJetMin = 0.3;
  fDRCandMin   = 0.1;
}

bool 
MetSignificance::filter(const mithep::PFJet* iJet,const mithep::Particle* iCandidate) 
{ 
  if(iCandidate == 0) return true;
  for(unsigned int i1 = 0; i1 < iJet->NConstituents(); i1++) {
    const mithep::PFCandidate* pJet = iJet->PFCand(i1);
    if(MathUtils::DeltaR(pJet->Mom(),iCandidate->Mom()) < fDRCandMin) return true;
  }
  return false;
}

bool 
MetSignificance::filter(const mithep::PFCandidate* iPart,const mithep::PFJetCol* iPFJetCol) 
{ 
  for(unsigned int i0 = 0; i0 < iPFJetCol->GetEntries(); i0++) { 
    for(unsigned int i1 = 0; i1 < iPFJetCol->At(i0)->NConstituents(); i1++) {
      const mithep::PFCandidate* pPart = iPFJetCol->At(i0)->PFCand(i1);
      if(MathUtils::DeltaR(iPart->Mom(),pPart->Mom()) < fDRCandMin) return true;
    }
  }
   return false;
}

void 
MetSignificance::addJets(std::vector<metsig::SigInputObj>& fSig, const mithep::PFJetCol* iJets, const mithep::Particle* iCan1, const mithep::Particle* iCan2) 
{ 
  for(unsigned int i0 = 0; i0 < iJets->GetEntries(); i0++) { 
    const mithep::PFJet *pJet = iJets->At(i0);
    if(iCan1 != 0) {if(MathUtils::DeltaR(pJet->Mom(),iCan1->Mom()) < fDRLepJetMin) continue;}
    if(iCan2 != 0) {if(MathUtils::DeltaR(pJet->Mom(),iCan2->Mom()) < fDRLepJetMin) continue;}
    
    if(filter(pJet,iCan1)) continue;
    if(filter(pJet,iCan2)) continue;
    add(pJet,fSig);
  }
}

void 
MetSignificance::addCandidates(const mithep::PFCandidateCol* iCands, const mithep::PFJetCol* iJets) 
{ 
  fCandSig.clear();
  for(unsigned int i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const mithep::PFCandidate* pCand = iCands->At(i0);
    //if(filter(pCand,iMuons    ,fDRCandMin)) continue;
    //if(filter(pCand,iElectrons,fDRCandMin)) continue;
    if(filter(pCand           ,iJets))      continue;
    std::pair<metsig::SigInputObj,mithep::FourVector> pCandSig;
    pCandSig.first  = fMetRes->evalPF(pCand);
    pCandSig.second = pCand->Mom();
    fCandSig.push_back(pCandSig);
  }
}

void 
MetSignificance::subtractCandidates(std::vector<metsig::SigInputObj>& fSig,const mithep::PFCandidateCol* iCands, const mithep::PFJetCol* iJets, const mithep::Particle* iCan1, const mithep::Particle* iCan2) 
{ 
  for(unsigned int i0 = 0; i0 < fCandSig.size(); i0++) { 
    //const mithep::PFCandidate* pCand = iCands->At(i0);
    if(iCan1 != 0) {if(MathUtils::DeltaR(fCandSig[i0].second,iCan1->Mom()) < fDRLepJetMin) continue;}
    if(iCan2 != 0) {if(MathUtils::DeltaR(fCandSig[i0].second,iCan2->Mom()) < fDRLepJetMin) continue;}
    //if(filter(pCand,iMuons    ,fDRCandMin)) continue;
    //if(filter(pCand,iElectrons,fDRCandMin)) continue;
    //if(filter(pCand           ,iJets))      continue;
    //add(pCand,fSig);
    fSig.push_back(fCandSig[i0].first);
  }
}

void 
MetSignificance::addTau(std::vector<metsig::SigInputObj> &fSig,const mithep::PFTau* iTau) { //const PFTauOArr *iTaus) { 
  if(iTau->SourcePFJet() == 0) return;
  add(iTau->SourcePFJet(),fSig);
}

void 
MetSignificance::add(const mithep::PFJet *iJet,std::vector<metsig::SigInputObj> &fSig) 
{ 
  fSig.push_back(fMetRes->evalPFJet(iJet));
}

void 
MetSignificance::add(const mithep::PFCandidate *iCandidate,std::vector<metsig::SigInputObj> &fSig) 
{ 
  fSig.push_back(fMetRes->evalPF(iCandidate));
}

void 
MetSignificance::add(const mithep::Muon *iMuon,std::vector<metsig::SigInputObj> &fSig) 
{ 
  double lPtE  = 0; double lPhiE = 0;       std::string lParticleType = "muon";
  double lPt   = iMuon->Pt(); double lPhi = iMuon->Phi(); double lEta = iMuon->Eta();
  //if(iMuon->track().isNonnull()) {
  if (iMuon->TrackerTrk() != 0) {
    lPtE  =              iMuon->TrackerTrk()->PtErr();
    lPhiE = iMuon->Pt()*(iMuon->TrackerTrk()->Phi0Err()); // XXXXX check CV: pt*dphi is indeed correct
  } else {
    lPtE  = fMetRes->eval(mithep::PFtype3, mithep::ET,  lPt, lPhi, lEta);
    lPhiE = fMetRes->eval(mithep::PFtype3, mithep::PHI, lPt, lPhi, lEta);
  }
  fSig.push_back(metsig::SigInputObj(lParticleType,lPt,lPhi,lPtE,lPhiE));
}

void 
MetSignificance::add(const mithep::Electron *iElectron,std::vector<metsig::SigInputObj> &fSig) 
{ 
  std::string lParticleType = "electron";
  double lPt   = iElectron->Pt(); double lPhi = iElectron->Phi(); double lEta = iElectron->Eta();
  double lPtE  = fMetRes->eval(mithep::PFtype2, mithep::ET,  lPt,lPhi,lEta);
  double lPhiE = fMetRes->eval(mithep::PFtype2, mithep::PHI, lPt,lPhi,lEta);
  fSig.push_back(metsig::SigInputObj(lParticleType,lPt,lPhi,lPtE,lPhiE));
}

TMatrixD MetSignificance::getSignificance(const PFJetCol *iJets, const PFCandidateCol *iCands, const mithep::PFTau *iTau, const mithep::PFTau *jTau, const mithep::Muon *iMuon, const mithep::Electron *iElectron)
{
  if(fComputeBase) { 
    addCandidates(iCands,iJets);
    fComputeBase = false;
  }
  std::vector<metsig::SigInputObj> fSig;
  //addJets      (fSig,iJets       ,0,0);//(mithep::Particle*) iMuon,(mithep::Particle*) iElectron);
  addJets           (fSig,iJets       ,(mithep::Particle*) iMuon,(mithep::Particle*) iElectron);
  subtractCandidates(fSig,iCands,iJets,(mithep::Particle*) iMuon,(mithep::Particle*) iElectron);
  if(iTau      != 0) addTau        (fSig,iTau);
  if(jTau      != 0) addTau        (fSig,iTau);
  if(iMuon     != 0) add (iMuon,    fSig);
  if(iElectron != 0) add (iElectron,fSig);
  
  metsig::significanceAlgo lMEtSignAlgorithm;
  lMEtSignAlgorithm.addObjects(fSig);
  TMatrixD lPFMETCov = lMEtSignAlgorithm.getSignifMatrix();

  const double lEpsilon = 1.e-4;
  if ( TMath::Abs(lPFMETCov.Determinant()) < lEpsilon ) {
    lPFMETCov(0,0) = 100;
    lPFMETCov(0,1) = 0.;
    lPFMETCov(1,0) = 0.;
    lPFMETCov(1,1) = 100;
  }
  return lPFMETCov;
}

void 
MetSignificance::loadResolutions() 
{ 
  fMetRes = new mithep::SignAlgoResolutions();
}
