#ifndef MITHTT_NTUPLER_TPHOTON_HH
#define MITHTT_NTUPLER_TPHOTON_HH

#include "TObject.h"
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"

/**
   \class TPhoton TPhoton.h MitHtt/Ntupler/include/TPhoton.h

   \brief Description: Bacon photon

   All information that is available on a Bacon photon.
*/

namespace mithep 
{
  class TPhoton : public TObject
  {
    public:
    /// default constructor
    TPhoton(){}
    /// default destructor
    ~TPhoton(){}
    
    /// kinematics of photon
    float pt, eta, phi;
    /// kinematics from super cluster
    float scEt, scEta, scPhi;
    /// filled from Bambu HollowConeTrkIsoDr04()
    float trkIso04;
    /// ECAL rec hit based isolation with isolation cone of 0.4 
    float emIso04;
    /// HCAL calo tower isolation with isolation cone of 0.4
    float hadIso04;
    /// hadronic over electromagnetic energy in the calorimeters
    float HoverE;
    /// ratio of energy od super cluster over energy in 3x3 matrix of ECAL crystals
    float R9;
    /// eta-width of shower in number of crystals
    float sigiEtaiEta;
    /// super cluster ID (filled from Bambu SCluster()->GetUniqueID())
    unsigned int  scID;
    /// indicated whether the super cluster has a seed in the pixel or not
    bool hasPixelSeed;
    /// HLT bits for which the offline reconstructed photon could be matched on trigger level
    TriggerBits hltMatchBits;
    
    ClassDef(TPhoton,1)
  };  
}
#endif
