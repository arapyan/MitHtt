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
    
    /// jet kinematics (pt is corrected up to L3/L2L3, ptraw is uncBorrected)
    float pt, ptraw, eta, phi, mass;
    /// beta 
    float beta;
    /// jet uncertainty as a unsigned relative uncertainty on the corrected jet pt
    float unc;
    /// jet area (from Fastjet)
    float area;
    /// combined secondary vertex btag discriminator
    float csv;
    /// jet MVA output
    float mva;
    /// moriond jet id
    float mvaold;
    /// pass or fail MVA id
    unsigned int id;
    /// pass or tail moriond mva id 
    unsigned int idold;
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
    /// pdgId of matched MC particle
    int matchedId;
    /// pdgId of matched GenJet parton flavor
    int matchedFlavor;
    /// pdgId of matched GenJet parton flavor
    float genpt,geneta,genphi;
    /// HLT bits for which the offline reconstructed jet could be matched on trigger level
    TriggerObjects  hltMatchBits;
    //l1 match
    bool l1match;
    /// Quark Likelihood
    float  quark;
    /// Gluon Likelhiood
    float  gluon;
    /// PU Likelhiood
    float  pu;
    
    ClassDef(TJet,9)
  };
}
#endif
