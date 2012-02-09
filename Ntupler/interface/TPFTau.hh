#ifndef MITHTT_NTUPLEDEF_TPFTAU_HH
#define MITHTT_NTUPLEDEF_TPFTAU_HH

#include "TObject.h"
#include "TClonesArray.h"

#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TPFCandidate.hh"

/**
   \class TPFTau TPFTau.h MitHtt/Ntupler/include/TPFTau.h

   \brief Description: Bacon tau

   All information that is available on a Bacon tau.
*/

namespace mithep 
{
  class TPFTau : public TObject
  {
  public:
    /// flags for HPS reconstruction
    enum hpsFlags{
      kLooseEle      = 0x00000001,
      kMediumEle     = 0x00000002,
      kTightEle      = 0x00000004,
      kLooseMu	     = 0x00000008,
      kTightMu	     = 0x00000010,
      kDecayMode     = 0x00000020,
      kVLooseIso     = 0x00000040,
      kLooseIso      = 0x00000080,
      kMediumIso     = 0x00000100,
      kTightIso	     = 0x00000200,
      kVLooseCombIso = 0x00000400,
      kLooseCombIso  = 0x00000800,
      kMediumCombIso = 0x00001000,
      kTightCombIso  = 0x00002000
    };

    /// default constructor
    TPFTau(){}
    /// default destructor
    ~TPFTau(){} 
    
    /// tau kinematics
    float pt, eta, phi, m, e;
    /// tau charge
    float q;
    /// signed d0 significance of the leading particle flow candidate
    float leadPFCandSignD0Sig;
    /// sum ot of the selected charged hadron particle flow candidate in the isolation cone
    float isoChargedHadronPtSum;
    /// sum et of sel. photon pfcands in iso cone
    float isoGammaEtSum; 
    /// hadronic energy over measured track momentum in 3x3 cluster of calo towers        
    float hcal3x3EOverP;
    /// hadronic energy over measured momentum for leading particle flow charged track (if it exists)
    float hcalOverP;
    /// electromagnetic energy over measured momentum for leading particle flow charged track (if it exists)
    float ecalOverP;
    /// tau electromagnetic fraction
    float emFraction;
    /// electron pre ID output (from Bambu ElectronPreIDOutput())
    float elePreIDOutput;
    /// absolute difference in eta btw leading photon candidate and leading track in tau
    float gammaDEta;
    /// absolute difference in phi btw leading photon candidate and leading track in tau
    float gammaDPhi;
    /// pt of leading photon candidate in tau
    float gammaPt;
    /// absolute difference in eta btw second leading photon candidate and leading track in tau
    float gamma2DEta;
    /// absolute difference in phi btw second leading photon candidate and leading track in tau
    float gamma2DPhi;
    /// pt of second leading photon candidate in tau
    float gamma2Pt;
    /// squared sum of differences in eta btw photon candidates and leading track in tau divided by linear sum pt of photon candidates  
    float gammaDEta2;
    /// squared sum of differences in phi btw photon candidates and leading track in tau divided by linear sum pt of photon candidates   
    float gammaDPhi2;
    /// linear sum of pt of photon candidates in tau divided by tau pt
    float gammaPtR;
    /// indicates whether the leading track in the tau candidate has a GSFTrack or not
    bool  hasGsf;
    /// HPS discriminators (according to hpsFlags)
    unsigned int  hpsDiscriminators;
    /// number of charges hadrons in tau candidate
    unsigned int nSignalPFChargedHadrCands;
    /// number of photons in tau candidate
    unsigned int nSignalPFGammaCands;
    /// indicated whether the tau candidate has a charged leading track or not
    bool hasLeadChargedHadronPFCand;
    /// leading charged hadron (if it exists, empty else)
    TPFCandidate leadChargedHadronPFCand;
    /// HLT bits for which the offline reconstructed tau could be matched on trigger level
    TriggerBits hltMatchBits;
    /// isolation using computePFTauIso 
    float isoEtPU;
    /// official tau POG isolation  
    float Iso;
    /// common isolation using pfPileup candidates, ptMin=0.5, dRMax=0.3, dRMin=0.0001
    float puIso;
    /// common isolation using all particle flow candidates, ptMin=0.5, dRMax=0.5, dRMin=0.0001
    float puIsoNoZ;
    /// common isolation using pfPileup candidates, ptMin=0., dRMax=0.5, dRMin=0.
    float puIsoNoPt;

    ClassDef(TPFTau, 1)
  };  
}
#endif
