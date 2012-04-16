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
    /// classic detector based track isolation with isolation cone of 0.4
    float trkIso04;
    /// classic detector based ECAL isolation with isolation cone of 0.4
    float emIso04;
    /// classic detector based HCAL isolation with isolation cone of 0.4
    float hadIso04;
    /// particle flow isolation with charged component restricted to the hard interaction vertex with isolation cone 0.3 and 0.4
    float pfIso03, pfIso04;    // Particle Flow isolation
    /// miscellaneous PF isolation variables used for MVA id with isolation cone of 0.3
    float chargedIso03, chargedIso03FromOtherVertices, neutralHadronIso03_01Threshold, gammaIso03_01Threshold, neutralHadronIso03_05Threshold, gammaIso03_05Threshold, neutralHadronIso03_10Threshold, gammaIso03_10Threshold, neutralHadronIso03_15Threshold, gammaIso03_15Threshold;
    /// miscellaneous PF isolation variables used for MVA id with isolation cone of 0.4
    float chargedIso04, chargedIso04FromOtherVertices, neutralHadronIso04_01Threshold, gammaIso04_01Threshold, neutralHadronIso04_05Threshold, gammaIso04_05Threshold, neutralHadronIso04_10Threshold, gammaIso04_10Threshold, neutralHadronIso04_15Threshold, gammaIso04_15Threshold;
    /// miscellaneous PF isolation variables with eta strip vetoes used for MVA id
    float chargedEMIsoVetoEtaStrip03, chargedEMIsoVetoEtaStrip04, neutralHadronIso007_01Threshold, gammaIsoVetoEtaStrip03_01Threshold, gammaIsoVetoEtaStrip04_01Threshold, neutralHadronIso007_05Threshold, gammaIsoVetoEtaStrip03_05Threshold, gammaIsoVetoEtaStrip04_05Threshold, neutralHadronIso007_10Threshold, gammaIsoVetoEtaStrip03_10Threshold, gammaIsoVetoEtaStrip04_10Threshold, neutralHadronIso007_15Threshold, gammaIsoVetoEtaStrip03_15Threshold, gammaIsoVetoEtaStrip04_15Threshold;
    /// charged PF isolation computed in rings of DR
    float chargedIso_DR0p0To0p1, chargedIso_DR0p1To0p2, chargedIso_DR0p2To0p3, chargedIso_DR0p3To0p4, chargedIso_DR0p4To0p5, chargedIso_DR0p5To0p7, chargedIso_DR0p7To1p0;
    /// gamma PF isolation computed in rings of DR
    float gammaIso_DR0p0To0p1, gammaIso_DR0p1To0p2, gammaIso_DR0p2To0p3, gammaIso_DR0p3To0p4, gammaIso_DR0p4To0p5, gammaIso_DR0p5To0p7, gammaIso_DR0p7To1p0;
    /// neutral PF isolation computed in rings of DR
    float neutralIso_DR0p0To0p1, neutralIso_DR0p1To0p2, neutralIso_DR0p2To0p3, neutralIso_DR0p3To0p4, neutralIso_DR0p4To0p5, neutralIso_DR0p5To0p7, neutralIso_DR0p7To1p0;
    /// shower shape variables from charged PF candidates
    float chargedIso_MeanEta, chargedIso_MeanPhi, chargedIso_SigEtaEta, chargedIso_SigEtaPhi, chargedIso_SigPhiPhi;
    /// shower shape variables from neutral PF candidates
    float neutralIso_MeanEta, neutralIso_MeanPhi, neutralIso_SigEtaEta, neutralIso_SigEtaPhi, neutralIso_SigPhiPhi;
    /// shower shape variables from gamma PF candidates
    float gammaIso_MeanEta, gammaIso_MeanPhi, gammaIso_SigEtaEta, gammaIso_SigEtaPhi, gammaIso_SigPhiPhi;
    /// directional isolation variables
    float directionalChargedIso, directionalNeutralIso, directionalGammaIso, directionalPFIso;
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
    float p;
    // 3d impact parameter and 3d impact paramter significance
    float ip3d, ip3dSig;       
    // impact parameters calculated with respect to unbiased vertex
    float d0Ub, d0UbSig, ip3dUb, ip3dUbSig;
    /// number of hits expected before the first hit has been observed
    unsigned int nExpHitsInner;
    /// number of bremsstahlung cluster (determined from NumberOfClusters-1) 
    float nBrem;
    /// energy of the seed cluster over pt at vertex and at calorimeter surface 
    float ESeedClusterOverPIn, ESeedClusterOverPOut;
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
    /// additional variables used for MVA id
    float gsfTrackChi2OverNdof, hcalDepth1OverEcal, hcalDepth2OverEcal, deltaEtaCalo, deltaPhiCalo, R9, scEtaWidth, scPhiWidth, coviEtaiPhi, psOverRaw, seedEMaxOverE, seedETopOverE, seedEBottomOverE, seedELeftOverE, seedERightOverE, seedE2ndOverE, seedE2x5RightOverE, seedE2x5LeftOverE, seedE2x5TopOverE, seedE2x5BottomOverE, seedE2x5MaxOverE, seedE1x3OverE, seedE3x1OverE, seedE1x5OverE, seedE2x2OverE, seedE3x2OverE, seedE3x3OverE, seedE4x4OverE, seedE5x5OverE;
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
    /// px and py of the matching particle flow candidate
    float pfPx, pfPy;

    ClassDef(TElectron, 4)
  };
}
#endif
