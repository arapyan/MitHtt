#ifndef MITHTT_NTUPLER_TELECTRON_HH
#define MITHTT_NTUPLER_TELECTRON_HH

#include "TObject.h"
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"

/**
   \class TElectron TElectron.h MitHtt/Ntupler/include/TElectron.h

   \brief Description: Bacon electron

   All information that is available on a Bacon electron.
*/

namespace mithep
{
  class TElectron : public TObject
  {
  public:
    /// default constructor
    TElectron(){}
    /// default destructor
    ~TElectron(){}
    
    /// kinematics of the electron
    float pt, eta, phi;
    /// classic detector based track isolation with isolation cone of 0.3
    float trkIso03;
    /// classic detector based ECAL isolation with isolation cone of 0.3
    float emIso03;
    /// classic detector based HCAL isolation with isolation cone of 0.3
    float hadIso03;
    /// impact parameter wrt the selected primary vertex along the z-axis and perpendicular to to the z-axis, impact parameter significance
    float d0, dz, d0Sig;
    /// kinematics from super cluster 
    float scEt, scEta, scPhi;
    /// hadronic energy over electromagnetic energy
    float HoverE;
    /// hadronic energy over momentum
    float EoverP;
    /// fbrem
    float fBrem;
    /// deltaEta btw super position and track at vertex
    float deltaEtaIn;
    /// deltaPhi btw super position and track at vertex
    float deltaPhiIn;
    /// eta width of the cluster in number of crystals
    float sigiEtaiEta;
    /// phi width of the cluster in number of crystals
    float sigiPhiiPhi;
    /// difference in cot(theta) for conversion partner tracks
    float partnerDeltaCot;
    /// distance in x-y plane to nearest conversion partner track
    float partnerDist; 
    /// electron charge  
    int q;
    /// electron supercluster energy
    float E;
    /// electron momentum
    float p, P;
    // 3d impact parameter and 3d impact paramter significance
    float ip3d, ip3dSig;       
    // impact parameters calculated with respect to unbiased vertex
    float d0Ub, d0UbSig, ip3dUb, ip3dUbSig;
    // impact parameters calculated with beam spot constraint
    float d0Bs, d0BsSig, ip3dBs, ip3dBsSig;
    /// number of hits expected before the first hit has been observed
    unsigned int nExpHitsInner;
    /// number of bremsstahlung cluster (determined from NumberOfClusters-1) 
    float nBrem;
    /// HLT bits for which the offline reconstructed electron could be matched on trigger level
    TriggerObjects  hltMatchBits;
    // supercluster ID (for matching to photon superclusters)
    unsigned int scID;
    /// unique track ID (filled from Bambu TrackerTrk()->GetUniqueID()) 
    unsigned int trkID;
    /// indicates whether measurement is ECAL driven
    bool isEcalDriven;
    /// indicates whether this come from a conversion (as a result of HttNtupler::isConversion())
    bool isConv;
    /// indicates whether this is an ECAL barrel electron
    bool isEB;
    /// mva values
    float mvaValID;
    /// common isolation from pfNoPileup, ptMin=0.0, dRMax=0.4, dRMin=0.015, type charged hadron  
    float pfIsoCharged;
    /// common isolation from pfNoPileup, ptMin=0.0, dRMax=0.4, dRMin=0.015, type charged hadron (no Z restriction)
    float pfIsoChargedNoZ;
    /// common isolation from pfNoPileup, ptMin=0.5, dRMax=0.4, dRMin=0.0   , type neutral hadron 
    float pfIsoNeutral;
    /// common isolation from pfNoPileup, ptMin=0.5, dRMax=0.4, dRMin=0.0   , type neutral hadron (no Z restriction) 
    float pfIsoNeutralNoZ;
    /// common isolation from pfNoPileup, ptMin=0.5, dRMax=0.4, dRMin=0.08  , type photon
    float pfIsoGamma;
    /// common isolation from pfNoPileup, ptMin=0.5, dRMax=0.4, dRMin=0.08  , type photon (no Z restriction) 
    float pfIsoGammaNoZ;
    /// naive isolation from pfPileup ptMin=0.5, dRMax=0.4, dRMin=0.01
    float puIso;
    /// naive isolation from pfPileup ptMin=0.5, dRMax=0.4, dRMin=0.01 (no Z restriction)
    float puIsoNoZ;
    // total delta beta corrected isolation
    float pfDeltaBetaIso;
    // number of tracker layers from associated ctf trk
    int kfhits;
    // additional variables for MVA ID
    float kfChi2, gsfChi2, etaWidth, phiWidth, deltaEtaCalo, E1x5E5x5, R9, eleEoPout, psOverRaw;

    ClassDef(TElectron, 7)
  };
}
#endif
