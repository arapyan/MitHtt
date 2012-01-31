#include "iostream"

#include "MitCommon/MathTools/interface/MathUtils.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "MitHtt/Ntupler/interface/MetSignificance.h"
#include "RecoMET/METAlgorithms/interface/significanceAlgo.h"

using namespace mithep;


MetSignificance::MetSignificance()  {
  loadResolutions();
  fDRLepJetMin = 0.3;
  fDRCandMin   = 0.1;
}
/*
bool MetSignificance::filter(const PFJet *iPart,MuonOArr *iParticleCol,double iDRMin) { 
  for(unsigned int i0 = 0; i0 < iParticleCol->GetEntries(); i0++) { 
    const mithep::Muon *pPart = iParticleCol->At(i0);
    if(MathUtils::DeltaR(iPart->Mom(),pPart->Mom()) < iDRMin) return true;
  }
  return false;
}
bool MetSignificance::filter(const PFCandidate *iPart,MuonOArr *iParticleCol,double iDRMin) { 
  for(unsigned int i0 = 0; i0 < iParticleCol->GetEntries(); i0++) { 
    const mithep::Muon *pPart = iParticleCol->At(i0);
    if(MathUtils::DeltaR(iPart->Mom(),pPart->Mom()) < iDRMin) return true;
  }
  return false;
}
bool MetSignificance::filter(const PFJet *iPart,ElectronOArr *iParticleCol,double iDRMin) { 
  for(unsigned int i0 = 0; i0 < iParticleCol->GetEntries(); i0++) { 
    const mithep::Electron *pPart = iParticleCol->At(i0);
    if(MathUtils::DeltaR(iPart->Mom(),pPart->Mom()) < iDRMin) return true;
  }
  return false;
}
bool MetSignificance::filter(const PFCandidate *iPart,ElectronOArr *iParticleCol,double iDRMin) { 
  for(unsigned int i0 = 0; i0 < iParticleCol->GetEntries(); i0++) { 
    const mithep::Electron *pPart = iParticleCol->At(i0);
    if(MathUtils::DeltaR(iPart->Mom(),pPart->Mom()) < iDRMin) return true;
  }
  return false;
}
*/
bool MetSignificance::filter(const mithep::PFJet *iJet,const mithep::Particle *iCandidate) { 
  if(iCandidate == 0) return true;
  //for(unsigned int i0 = 0; i0 < iParticleCol->GetEntries(); i0++) { 
  //const mithep::Electron *pPart = iParticleCol->At(i0);
  for(unsigned int i1 = 0; i1 < iJet->NConstituents(); i1++) {
    const mithep::PFCandidate *pJet = iJet->PFCand(i1);
    if(MathUtils::DeltaR(pJet->Mom(),iCandidate->Mom()) < fDRCandMin) return true;
    //if(deltaR(pJet.p4(),iCandidate->p4()) < fDRCandMin) return true;
  }
  return false;
}
/*
bool MetSignificance::filter(const PFJet *iJet,ElectronOArr *iParticleCol) { 
  for(unsigned int i0 = 0; i0 < iParticleCol->GetEntries(); i0++) { 
    const mithep::Electron *pPart = iParticleCol->At(i0);
    for(unsigned int i1 = 0; i1 < iJet->NConstituents(); i1++) {
      const mithep::PFCandidate *pJet = iJet->PFCand(i1);
      if(MathUtils::DeltaR(pJet->Mom(),pPart->Mom()) < fDRCandMin) return true;
    }
  }
   return false;
}
bool MetSignificance::filter(const PFJet *iJet,MuonOArr *iParticleCol) { 
  for(unsigned int i0 = 0; i0 < iParticleCol->GetEntries(); i0++) { 
    const mithep::Muon *pPart = iParticleCol->At(i0);
    for(unsigned int i1 = 0; i1 < iJet->NConstituents(); i1++) {
      const mithep::PFCandidate *pJet = iJet->PFCand(i1);
      if(MathUtils::DeltaR(pJet->Mom(),pPart->Mom()) < fDRCandMin) return true;
    }
  }
   return false;
}
bool MetSignificance::filter(const PFCandidate *iPart,const PFJetCol *iPFJetCol) { 
  for(unsigned int i0 = 0; i0 < iPFJetCol->GetEntries(); i0++) { 
    for(unsigned int i1 = 0; i1 < iPFJetCol->At(i0)->NConstituents(); i1++) {
      const mithep::PFCandidate *pPart = iPFJetCol->At(i0)->PFCand(i1);
      if(MathUtils::DeltaR(iPart->Mom(),pPart->Mom()) < fDRCandMin) return true;
    }
  }
   return false;
}

*/
bool MetSignificance::filter(const mithep::PFCandidate *iPart,const mithep::PFJetCol *iPFJetCol) { 
  for(unsigned int i0 = 0; i0 < iPFJetCol->GetEntries(); i0++) { 
    for(unsigned int i1 = 0; i1 < iPFJetCol->At(i0)->NConstituents(); i1++) {
      const mithep::PFCandidate  *pPart = iPFJetCol->At(i0)->PFCand(i1);
      if(MathUtils::DeltaR(iPart->Mom(),pPart->Mom()) < fDRCandMin) return true;
      //if(deltaR(iPart->p4(),pPart.p4()) < fDRCandMin) return true;
    }
  }
   return false;
}

void MetSignificance::addJets(std::vector<metsig::SigInputObj> &fSig,const mithep::PFJetCol *iJets,
			      const mithep::Particle *iCan1,const mithep::Particle *iCan2) { 
  for(unsigned int i0 = 0; i0 < iJets->GetEntries(); i0++) { 
    const mithep::PFJet *pJet = iJets->At(i0);
    if(iCan1 != 0) {if(MathUtils::DeltaR(pJet->Mom(),iCan1->Mom()) < fDRLepJetMin) continue;}
    if(iCan2 != 0) {if(MathUtils::DeltaR(pJet->Mom(),iCan2->Mom()) < fDRLepJetMin) continue;}
    
    if(filter(pJet,iCan1)) continue;
    if(filter(pJet,iCan2)) continue;
    add(pJet,fSig);
  }
}
void MetSignificance::addCandidates(std::vector<metsig::SigInputObj> &fSig,const mithep::PFCandidateCol *iCands,
				    const mithep::PFJetCol *iJets,const mithep::Particle *iCan1,const mithep::Particle *iCan2) { 
  for(unsigned int i0 = 0; i0 < iCands->GetEntries(); i0++) { 
    const mithep::PFCandidate* pCand = iCands->At(i0);
    if(iCan1 != 0) {if(MathUtils::DeltaR(pCand->Mom(),iCan1->Mom()) < fDRLepJetMin) continue;}
    if(iCan2 != 0) {if(MathUtils::DeltaR(pCand->Mom(),iCan2->Mom()) < fDRLepJetMin) continue;}
    //if(filter(pCand,iMuons    ,fDRCandMin)) continue;
    //if(filter(pCand,iElectrons,fDRCandMin)) continue;
    if(filter(pCand           ,iJets))      continue;
    add(pCand,fSig);
  }
}
void MetSignificance::addTau(std::vector<metsig::SigInputObj> &fSig,const mithep::PFTau* iTau) { //const PFTauOArr *iTaus) { 
  //for(unsigned int i0 = 0; i0 < iTaus->GetEntries(); i0++) { 
  //const PFTau *pTau = iTaus->At(i0);
  if(iTau->SourcePFJet() == 0) return;
  add(iTau->SourcePFJet(),fSig);
  //if(iTau->pfJetRef().get() == 0) return;
  //add(iTau->pfJetRef().get(),fSig);
}
/*
void MetSignificance::addMuons(std::vector<metsig::SigInputObj> &fSig,const MuonOArr *iMuons) {
  for(unsigned int i0 = 0; i0 < iMuons->GetEntries(); i0++) { 
    const Muon *pMuon = iMuons->At(i0);
    add(pMuon,fSig);
  }
}
void MetSignificance::addElectrons(std::vector<metsig::SigInputObj> &fSig,const ElectronOArr *iElectrons) {
  for(unsigned int i0 = 0; i0 < iElectrons->GetEntries(); i0++) { 
    const Electron *pElectron = iElectrons->At(i0);
    add(pElectron,fSig);
  }
}
*/
void MetSignificance::add(const mithep::PFJet *iJet,std::vector<metsig::SigInputObj> &fSig) { 
  fSig.push_back(fMetRes->evalPFJet(iJet));
}
void MetSignificance::add(const mithep::PFCandidate *iCandidate,std::vector<metsig::SigInputObj> &fSig) { 
  fSig.push_back(fMetRes->evalPF(iCandidate));
}
void MetSignificance::add(const mithep::Muon *iMuon,std::vector<metsig::SigInputObj> &fSig) { 
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

void MetSignificance::add(const mithep::Electron *iElectron,std::vector<metsig::SigInputObj> &fSig) { 
  std::string lParticleType = "electron";
  double lPt   = iElectron->Pt(); double lPhi = iElectron->Phi(); double lEta = iElectron->Eta();
  double lPtE  = fMetRes->eval(mithep::PFtype2, mithep::ET,  lPt,lPhi,lEta);
  double lPhiE = fMetRes->eval(mithep::PFtype2, mithep::PHI, lPt,lPhi,lEta);
  fSig.push_back(metsig::SigInputObj(lParticleType,lPt,lPhi,lPtE,lPhiE));
}
TMatrixD MetSignificance::getSignificance(const PFJetCol *iJets,const PFCandidateCol *iCands,
					  const mithep::PFTau *iTau,const mithep::Muon *iMuon,const mithep::Electron *iElectron) {
					  //PFTauOArr *iTaus,MuonOArr *iMuons,ElectronOArr *iElectrons) {
  std::vector<metsig::SigInputObj> fSig;
  addJets      (fSig,iJets       ,(mithep::Particle*) iMuon,(mithep::Particle*) iElectron);
  addCandidates(fSig,iCands,iJets,(mithep::Particle*) iMuon,(mithep::Particle*) iElectron);//iMuons,iElectrons);
  if(iTau      != 0) addTau        (fSig,iTau);
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

void MetSignificance::loadResolutions() { 
  fMetRes = new mithep::SignAlgoResolutions();
}
//void MetSignificance::setup(const PFJetCol *iJets,const PFCandidateCol *iCands,const PFTauCol *iTaus,
//			    const MuonCol *iMuons,const ElectronCol *iElectrons) {
//  std::cout << "--> Check" << std::endl;
//  std::cout << "===> " << fJets << std::endl;
//  fJets      = iJets;
//  fPFCands   = iCands;
//  std::cout << "===> Again" << std::endl;
//  fTaus      = iTaus;
//  fMuons     = iMuons;
//  std::cout << "===> WTF" << std::endl;
//  fElectrons = iElectrons;
//}
