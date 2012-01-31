#ifndef MITHTT_NTUPLER_METSIGNIFICANCE_HH
#define MITHTT_NTUPLER_METSIGNIFICANCE_HH
#include <TVectorD.h>

//#include "DataFormats/Candidate/interface/Candidate.h"
//#include "DataFormats/JetReco/interface/PFJet.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"

//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/Muon.h"
//#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "DataFormats/PatCandidates/interface/Tau.h"

#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"

#include "RecoMET/METAlgorithms/interface/SigInputObj.h"
#include "MitHtt/Ntupler/interface/SignAlgoResolutions.h"


namespace mithep { 
  class MetSignificance { 
  public:
    MetSignificance();
    /*
    bool filter(const PFJet *iPart,   MuonOArr *iParticleCol,double iDRMin);
    bool filter(const PFJet *iPart,   ElectronOArr *iParticleCol,double iDRMin);
    bool filter(const PFCandidate *iPart,   MuonOArr *iParticleCol,double iDRMin);
    bool filter(const PFCandidate *iPart,   ElectronOArr *iParticleCol,double iDRMin);
    bool filter(const PFJet *iJet,       MuonOArr     *iParticleol);
    bool filter(const PFJet *iJet,       ElectronOArr *iParticleol);
    bool filter(const PFCandidate *iPart,const PFJetCol    *iPFJetCol);
    void addJets      (std::vector<metsig::SigInputObj> &fSig,const PFJetCol *iJets,MuonOArr *iMuons,ElectronOArr *iElectrons);
    void addCandidates(std::vector<metsig::SigInputObj> &fSig,const PFCandidateCol *iPFCands,const PFJetCol *iJets,
		       MuonOArr *iMuons,ElectronOArr *iElectrons);
    void addMuons     (std::vector<metsig::SigInputObj> &fSig,const MuonOArr  *iMuons);
    void addElectrons (std::vector<metsig::SigInputObj> &fSig,const ElectronOArr *iElectrons);
    */
    bool filter(const mithep::PFJet *iJet,       const mithep::Particle *iCandidate);
    bool filter(const mithep::PFCandidate *iPart,const mithep::PFJetCol  *iPFJetCol);
    void addJets(std::vector<metsig::SigInputObj> &fSig,const mithep::PFJetCol *iJets,
		 const mithep::Particle *iCan1,const mithep::Particle *iCan2);
    void addCandidates(std::vector<metsig::SigInputObj> &fSig,const mithep::PFCandidateCol *iCands,
		       const mithep::PFJetCol *iJets,const mithep::Particle *iCan1,const mithep::Particle *iCan2);
    void addTau       (std::vector<metsig::SigInputObj> &fSig,const mithep::PFTau *iTau);
    void add(const mithep::PFJet       *iJet,      std::vector<metsig::SigInputObj> &fSig);
    void add(const mithep::PFCandidate *iCandidate,std::vector<metsig::SigInputObj> &fSig);
    void add(const mithep::Muon        *iMuon,     std::vector<metsig::SigInputObj> &fSig);
    void add(const mithep::Electron    *iElectron, std::vector<metsig::SigInputObj> &fSig);
    TMatrixD getSignificance(const mithep::PFJetCol *iJets,const mithep::PFCandidateCol *iCands,
			     const mithep::PFTau    *iTau ,const mithep::Muon *iMuon,const mithep::Electron *iElectron);
    void loadResolutions();
    //void setup(const PFJetCol *iJets,const PFCandidateCol *iCands,const PFTauCol *iTaus,
    //const MuonCol *iMuons,const ElectronCol *iElectrons);

  private:
    double                           fDRLepJetMin;
    double                           fDRCandMin;
    mithep::SignAlgoResolutions   *fMetRes;
    //std::vector<metsig::SigInputObj> *fSig;
    //const MuonCol                   *fMuons;
    //const ElectronCol               *fElectrons;
    //const PFTauCol                  *fTaus;
    //const PFJetCol                  *fJets;
    //const PFCandidateCol            *fPFCands;
  };
}
#endif
