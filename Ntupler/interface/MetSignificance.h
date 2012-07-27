#ifndef MITHTT_NTUPLER_METSIGNIFICANCE_HH
#define MITHTT_NTUPLER_METSIGNIFICANCE_HH

#include <vector>
#include <utility>

#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFJetCol.h"
#include "MitAna/DataTree/interface/PFCandidateCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"

#include "TVectorD.h"
#include "RecoMET/METAlgorithms/interface/SigInputObj.h"
#include "MitHtt/Ntupler/interface/SignAlgoResolutions.h"

/**
   \class MetSignificance MetSignificance.h MitHtt/Ntupler/include/MetSignificance.h

   \brief Description: Class to calculate the MET significance as input to the svfit

   This is a class to calculate the MET significance from particle flow candidates as
   input to the svfit. Input objects are all selected particle flow jets, all selected 
   particle flow candidates outside jets (excluding the selected electron, muon and/or 
   tau) and explicitely the selected electron, muon and/or tau. The actual significance 
   matrix (2d) is calculatec from metsig::significanceAlgo in RecoMET/METAlgorithms. 

   This class is a mirror of the same class in RecoMET/METAlgorithms but replacing 
   official dataformats by Bambu dataformats and specializing it to the needs for Higgs 
   to Tau Tau analyses.
*/

namespace mithep { 
  class MetSignificance { 
  public:

    /// default constructor
    MetSignificance();
    /// returns true if the particle flow candidate is located within a radius of fDRCandMin in the vicinity of the particle flow jet
    bool filter(const mithep::PFJet* iJet, const mithep::Particle* iCandidate);
    /// returns true if the particle flow candidate is located within a radius of fDRCandMin in the vicinity of any particle flow jet in the collection
    bool filter(const mithep::PFCandidate* iPart, const mithep::PFJetCol* iPFJetCol);
    /// add jets to the input for the significance calculation if they do not contain iCan1 or iCan2
    void addJets(std::vector<metsig::SigInputObj>& fSig, const mithep::PFJetCol* iJets, const mithep::Particle* iCan1, const mithep::Particle* iCan2);
    /// add particle flow candidates to the input for the significance calculation if they do not belong to the jets in iJets and if they are not iCan1 or iCan2
    void addCandidates(const mithep::PFCandidateCol *iCands, const mithep::PFJetCol *iJets);
    /// remove pf candidates associated to lepton and compuate pf cand part of significance
    void subtractCandidates(std::vector<metsig::SigInputObj>& fSig,const mithep::PFCandidateCol* iCands, const mithep::PFJetCol* iJets, const mithep::Particle* iCan1, const mithep::Particle* iCan2); 
    /// add particle flow jet corresponding to particle flow tau to the input for the significance calculation
    void addTau(std::vector<metsig::SigInputObj> &fSig,const mithep::PFTau *iTau);
    /// add particle flow jet to input for the significance calculation
    void add(const mithep::PFJet       *iJet,      std::vector<metsig::SigInputObj> &fSig);
    /// add particle flow candidate to input for the significance calculation
    void add(const mithep::PFCandidate *iCandidate,std::vector<metsig::SigInputObj> &fSig);
    /// add muon to input for the significance calculation
    void add(const mithep::Muon        *iMuon,     std::vector<metsig::SigInputObj> &fSig);
    /// add electron to input for the significance calculation
    void add(const mithep::Electron    *iElectron, std::vector<metsig::SigInputObj> &fSig);
    /// get the MET significance
    TMatrixD getSignificance(const mithep::PFJetCol *iJets,const mithep::PFCandidateCol *iCands, const mithep::PFTau    *iTau , const mithep::PFTau    *jTau, const mithep::Muon *iMuon,const mithep::Electron *iElectron);
    /// load resolutions and instantiate SgnAlgoResolutions object
    void loadResolutions();
    /// reset compuate base
    void reset() { fComputeBase = true;}
  private:
    /// minimal distance between jet and leptons (for electorn and muon)
    double fDRLepJetMin;
    /// minimal distance between particle flow candidate and particle flow jet
    double fDRCandMin;
    /// contains object resolutions for MET significance calculation
    mithep::SignAlgoResolutions* fMetRes;
    /// container of cleaned pf candidates (speed up processing)
    std::vector<std::pair<metsig::SigInputObj,mithep::FourVector> > fCandSig;
    /// force compuatation of pfcandidate calculation 
    bool fComputeBase;
  };
}
#endif
