#ifndef MITHTT_NTUPLER_TMUONLITE_HH
#define MITHTT_NTUPLER_TMUONLITE_HH

#include <TObject.h>

namespace mithep 
{
  class TMuonLite : public TObject
  {
    public:
      TMuonLite(){}
      ~TMuonLite(){} 
  
      Float_t pt, eta, phi;    // kinematics
      Float_t trkIso03;	       // track isolation
      Float_t emIso03;	       // ECAL-based isolation
      Float_t hadIso03;	       // HCAL-based isolation
      Float_t pfIso;           // Particle Flow isolation
      Float_t d0, dz;          // impact parameter
      Float_t ip2d, ip2dSig;
      Float_t ip3d, ip3dSig;      
      Int_t   q;	       // charge
      UInt_t  hltMatchBits;    // bits for matching with HLT primitives
      Bool_t  passId;          // pass muon ID 

    ClassDef(TMuonLite,1)
  };  
}
#endif
