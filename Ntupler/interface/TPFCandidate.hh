#ifndef MITHTT_NTUPLER_TPFCANDIDATE_HH
#define MITHTT_NTUPLER_TPFCANDIDATE_HH

#include "TObject.h"

namespace mithep { 
  class TPFCandidate : public TObject {
  public:
    enum EPFType { 
      eX = 0,          //unidentified
      eHadron,         //charged hadron
      eElectron,       //electron
    eMuon,           //muon
      eGamma,          //photon
      eNeutralHadron,  //neutral hadron
      eHadronHF,       //hadron in HF
      eEGammaHF        //EM object in HF
    };
    
    TPFCandidate() :
      pt(0),
      eta(0),
      phi(0),
      m(0),
      e(0),
      q(0),
      mvaEPi(0),
      pfType(0),
      hasMuon(kFALSE),
      nSeg(0),
      nMatches(0),
      hasTrack(kFALSE),
    d0(0), dz(0) {}
    
    ~TPFCandidate(){} 
    
    Float_t pt, eta, phi, m, e;    // kinematics
    Float_t   q;	           // charge
    
    Float_t mvaEPi;           // electron-pion discriminant
    UInt_t  pfType;           // particle flow type
    Bool_t  hasMuon;          // has a corresponding muon
    UInt_t  nSeg;             // number of muon segments
    UInt_t  nMatches;         // number of muon chambers matched to segments
    Bool_t  hasTrack;
    Float_t d0, dz;
    
    ClassDef(TPFCandidate,2)
  };
}
#endif
