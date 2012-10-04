#ifndef MITHTT_NTUPLER_TMUON_HH
#define MITHTT_NTUPLER_TMUON_HH

#include "TObject.h"
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"

/**
   \class TMuon TMuon.h MitHtt/Ntupler/include/TMuon.h

   \brief Description: Bacon muon

   All information that is available on a Bacon muon.
*/

namespace mithep 
{
  class TMuon : public TObject
  {
  public:
    /// default constructor
    TMuon(){}
    /// default destructor
    ~TMuon(){} 
    
    /// kinematics of the muon
    float pt, eta, phi, ptErr;
    /// classic detector based track isolation with isolation cone of 0.3
    float trkIso03;
    /// classic detector based ECAL isolation with isolation cone of 0.3
    float emIso03;
    /// classic detector based HCAL isolation with isolation cone of 0.3
    float hadIso03;
    /// classic detector based HO isolation with isolation cone of 0.3
    float hoIso03;
    /// impact parameter wrt the selected primary vertex along the z-axis and perpendicular to to the z-axis, impact parameter significance
    float d0, dz, d0Sig;
    /// 3d impact parameter and 3d impact parameter significance
    float ip3d, ip3dSig;
    // impact parameters calculated with respect to unbiased vertex
    float d0Ub, d0UbSig, ip3dUb, ip3dUbSig;
    // impact parameters calculated with beam spot constraint
    float d0Bs, d0BsSig, ip3dBs, ip3dBsSig;
    /// track chi**2/ndof
    float tkNchi2;
    /// muon fit chi**2/ndof (in order global, standalone, tracker); first come first fill
    float muNchi2;
    /// muon charge
    int q;
    /// number of valid hits in muon chambers (filled directly from Bambu NValidHits)
    int nValidHits;
    /// muon quality bit mask (filled from Bambu Quality().QualityMask().Mask())
    unsigned int qualityBits;
    /// muon type bits (kGlobal, kTracker, kStandalone)
    unsigned int  typeBits;	              // global muon, tracker muon, or standalone muon
    /// number of tracker hits (0 if the muon does not have a tracker track)
    unsigned int nTkHits;
    /// number of pixel hits
    unsigned int nPixHits;
    /// number of hit segments in the muon system 
    unsigned int nSeg;
    /// number of muon chambers that match to track segments 
    unsigned int nMatch;
    /// HLT bits for which the offline reconstructed muon could be matched on trigger level
    TriggerObjects hltMatchBits;
    /// unique track ID (filled from Bambu TrackerTrk()->GetUniqueID())    
    unsigned int trkID;
    /// common isolation from pfNoPileup, ptMin=0.0, dRMax=0.4, dRMin=0.0001, type charged hadron  
    float pfIsoCharged;
    /// common isolation from pfNoPileup, ptMin=0.0, dRMax=0.4, dRMin=0.0001, type charged hadron (no Z restriction)
    float pfIsoChargedNoZ;
    /// common isolation from pfNoPileup, ptMin=0.5, dRMax=0.4, dRMin=0.01  , type neutral hadron  
    float pfIsoNeutral;
    /// common isolation from pfNoPileup, ptMin=0.5, dRMax=0.4, dRMin=0.01  , type neutral hadron (no Z restriction)
    float pfIsoNeutralNoZ;
    /// common isolation from pfNoPileup, ptMin=0.5, dRMax=0.4, dRMin=0.01  , type photon
    float pfIsoGamma;
    /// common isolation from pfNoPileup, ptMin=0.5, dRMax=0.4, dRMin=0.01  , type photon (no Z restriction)
    float pfIsoGammaNoZ;
    /// naive isolation from pfPileup ptMin=0.5, dRMax=0.4, dRMin=0.01
    float puIso;
    /// naive isolation from pfPileup ptMin=0.5, dRMax=0.4, dRMin=0.01 (no Z restriction)
    float puIsoNoZ;
    /// px and py of the matching particle flow candidate
    float pfPx, pfPy;
    /// check if the muon is matched to a PF candidate
    bool matchesPFCand;
    // total delta beta corrected isolation
    float pfDeltaBetaIso;
    /// check type of matched PF candidate
    float nTkLayersHit;
    unsigned int matchedPFType;
    ClassDef(TMuon, 9)
  };  
}
#endif
