#ifndef MITHTT_NTUPLER_TJET_HH
#define MITHTT_NTUPLER_TJET_HH

#include <TObject.h>

namespace mithep 
{
  class TJet : public TObject
  {
    public:
      TJet(){}
      ~TJet(){}

      Float_t pt, eta, phi, mass;  // kinematics
      Float_t unc;                 // energy scale uncertainty
      Float_t area;                // jet area
      Float_t tche;                // TrackCountingHighEfficiency b-tag discriminator
      Float_t tchp;                // TrackCountingHighPurity b-tag discriminator
      Int_t   mcFlavor;            // PDG ID of matched parton flavor
      UInt_t  hltMatchBits;        // bits from matching with HLT primitives

    ClassDef(TJet,2)
  };
}
#endif