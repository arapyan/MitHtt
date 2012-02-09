#ifndef MITHTT_NTUPLER_TJET_HH
#define MITHTT_NTUPLER_TJET_HH

#include <TObject.h>
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"

/**
   \class TJet TJet.h MitHtt/Ntupler/include/TJet.h

   \brief Description: Bacon jet

   All information that is available on a Bacon jet.
*/

namespace mithep 
{
  class TJet : public TObject
  {
  public:
    /// default constructor
    TJet(){}
    /// default destructor
    ~TJet(){}
    
    /// jet kinematics (pt is corrected up to L3/L2L3, ptraw is uncorrected)
    float pt, ptraw, eta, phi, mass;
    /// beta 
    float beta;
    /// jet uncertainty as a unsigned relative uncertainty on the corrected jet pt
    float unc;
    /// jet area (from Fastjet)
    float area;
    /// track counting high efficiency btag discriminator
    float tche;
    /// track counting high purity btag discriminator
    float tchp;
    /// combined secondary vertex btag discriminator
    float csv;
    /// combined secondary vertex MVA discriminator
    float csvMva;
    /// number of charged hardons in jet
    unsigned int nCharged;
    /// charged electromagnetic energy over uncorrected jet energy
    float chgEMfrac;
    /// neutral electromagnetic energy over uncorrected jet energy
    float neuEMfrac;
    /// charged hadronic energy over uncorrected jet energy
    float chgHadrfrac;
    /// neutral hadronic energy over uncorrected jet energy
    float  neuHadrfrac;
    /// pdgId of matched parton flavour
    int mcFlavor;
    /// pdgId of matched hen jet
    int matchedId;
    /// HLT bits for which the offline reconstructed jet could be matched on trigger level
    TriggerBits  hltMatchBits;
    
    ClassDef(TJet,4)
  };
}
#endif
