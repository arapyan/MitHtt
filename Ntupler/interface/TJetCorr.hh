#ifndef MITHTT_NTUPLER_TJETCORR_HH
#define MITHTT_NTUPLER_TJETCORR_HH

#include <TObject.h>
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"

/**
   \class TJet TJet.h MitHtt/Ntupler/include/TJet.h

   \brief Description: Bacon jet

   All information that is available on a Bacon jet.
*/

namespace mithep 
{
  class TJetCorr : public TObject
  {
  public:
    /// default constructor
    TJetCorr(){}
    /// default destructor
    ~TJetCorr(){}
    
    /// jet kinematics (pt is corrected up to L3/L2L3, ptraw is uncBorrected)
    float pt, ptraw, eta, phi, mass,area;
    /// jet kinematics for pruning type1 : 0.1, 0.5
    float pt_p1, ptraw_p1, eta_p1, phi_p1, mass_p1, area_p1;
    /// jet kinematics for pruning type2 : 0.1, 0.2
    float pt_p2, ptraw_p2, eta_p2, phi_p2, mass_p2, area_p2;
    /// jet kinematics for filter type1 : 0.2, 3
    float pt_f1, ptraw_f1, eta_f1, phi_f1, mass_f1, area_f1;
    /// jet kinematics for filter type2 : 0.3, 3
    float pt_f2, ptraw_f2, eta_f2, phi_f2, mass_f2, area_f2;
    /// jet kinematics for trimmer type1 : 0.2, 0.05
    float pt_t1, ptraw_t1, eta_t1, phi_t1, mass_t1, area_t1;
    /// jet kinematics for trimmer type2 : 0.2, 0.03
    float pt_t2, ptraw_t2, eta_t2, phi_t2, mass_t2, area_t2;
    /// jet kinematics for massdrop type1 : 0.33
    float pt_m1, ptraw_m1, eta_m1, phi_m1, mass_m1, area_m1;
    /// jet kinematics for massdrop type2 : 0.67
    float pt_m2, ptraw_m2, eta_m2, phi_m2, mass_m2, area_m2;
    /// N subjettiness variables
    float tau1,tau2,tau3,tau4;

    /// beta 
    //float beta;
    /// jet uncertainty as a unsigned relative uncertainty on the corrected jet pt
    float unc;
    /// combined secondary vertex btag discriminator
    //float csv;
    /// jet MVA output
    //float mva;
    /// moriond jet id
    //float mvaold;
    /// pass or fail MVA id
    // unsigned int id;
    /// pass or tail moriond mva id 
    //unsigned int idold;
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
    //float  quark;
    /// Gluon Likelhiood
    //float  gluon;
    /// PU Likelhiood
    //float  pu;
    
    ClassDef(TJetCorr,9)
  };
}
#endif
