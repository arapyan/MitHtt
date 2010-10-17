//--------------------------------------------------------------------------------------------------
// $Id: $
//
//--------------------------------------------------------------------------------------------------

#ifndef MITHTT_MODS_HTTANALYSIS_H
#define MITHTT_MODS_HTTANALYSIS_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class MCEventInfo;

namespace mithep 
{
  class HttAnalysis : public BaseMod
  {
  public:
    HttAnalysis(const char *name  = "HggAnalysis", 
		const char *title = "Higgs To Tau Tau Analysis");
    
  protected:
    void                     Begin();
    void                     Process();
    void                     SlaveBegin();
    void                     SlaveTerminate();
    void                     Terminate();

    void                     SetTrigObjsName(const char *name) { fTrigObjsName = name; }
    void                     SetPhotonName(TString name)       { fPhotonName = name; }
    void                     SetPhotonsFromBranch(Bool_t b)    { fPhotonsFromBranch = b; }
    void                     SetIsData(Bool_t b)               { fIsData = b; }

    void                     ResetSelectedPhotons();    // reset stored photon indices
    void                     MatchPhotonsToTrigger();   // create trig photon indices from sel.
    Bool_t                   PhotonTriggered(UInt_t i); // did this photon trigger
    
    void                     SetOverlapCut(double c)           { fOverlapCut = c; }

    TString                  fTrigObjsName;             // name of trigger objects
    TString                  fPhotonName;               // name of photon collection
    Bool_t                   fPhotonsFromBranch;        // photons are loaded from a branch
    TString                  fMcEventInfoName;          // name of MC Event information
    const PhotonCol         *fPhotons;                  //! photon from data stream
    const MCEventInfo       *fMcEventInfo;              //! MC event information branch
    Bool_t                   fIsData;                   // looking at Data (or MC)

    // general selection cuts
    double                   fOverlapCut;               // cut to allow for rejection of overlap
    // photon selection cuts
    double                   fPhotonEtMin;              // minimal transverse energy for photon
    Bool_t                   fNoPixelSeedCut;           // request that there be no pixel seed?
    double                   fPhotonSigEeBarMax;        // maximal sigma(eta,eta), barrel
    double                   fPhotonSigEeEcpMax;        // maximal sigma(eta,eta), endcap
    double                   fPhotonIsoMax;             // maximal combined isolation
    // high mass di-photon selection cuts
    double                   f2PhotonMassMin;           // minimum di-photon mass, high mass sample
    // R9 cleaned photon selection cuts
    double                   fPhotonR9Min;              // minimum R9 for cleaner photons
 
    // keep track of our photons
    std::vector<UInt_t>      fSelectedPhotons;          // photon indices, selected
    std::vector<UInt_t>      fTriggeredPhotons;         // photon indices, selected and triggered

    //----------------------------------------------------------------------------------------------
    // histograms
    //----------------------------------------------------------------------------------------------
    // all photons
    TH1D                    *hAllNPhotons;         //! all number of photons
    TH1D                    *hAllPhotonPt;         //! all photon pt
    TH1D                    *hAllPhotonEta;        //! all photon eta
    TH1D                    *hAllPhotonPixelSeed;  //! pixel seed or not
    TH1D                    *hAllPhotonSigEeBar;   //! sigma(eta,eta) of photons, barrel
    TH1D                    *hAllPhotonSigEeEcp;   //! sigma(eta,eta) of photons, endcap
    TH1D                    *hAllPhotonEcHcTrIso;  //! combined isolation
    // trigger monitoring
    TH1D                    *hTrigDeltaEta;        //! trigger match (delta Eta)
    TH1D                    *hTrigDeltaPhi;        //! trigger match (delta Phi)
    TH1D                    *hTrigDeltaEt;         //! trigger match (delta Et)
    // selected photons and di-photons
    TH1D                    *hNPhotons;            //! number of photons, selected
    TH1D                    *hPhotonEtMx;          //! max Et photon
    TH1D                    *hPhotonSigEtaEta;     //! sigma_eta,eta for photon
    TH1D                    *hPhotonEcHcTrIso;     //! combined ECAL, HCAL, Tracker isolation
    TH1D                    *hN2Photons;           //! number of di-photons
    TH1D                    *hPhotonEt1;           //! Et 1 photon
    TH1D                    *hPhotonEt2;           //! Et 2 photon
    TH1D                    *hPhotonEta1;          //! Eta 1 photon
    TH1D                    *hPhotonEta2;          //! Eta 2 photon
    TH1D                    *hPhotonPhi1;          //! Phi 1 photon
    TH1D                    *hPhotonPhi2;          //! Phi 2 photon
    TH1D                    *hPhotonR91;           //! R9  1 photon
    TH1D                    *hPhotonR92;           //! R9  2 photon
    TH1D                    *hPhotonDelR;          //! delta R(photon1,photon2)
    TH1D                    *h2PhotonEtSum;        //! di-photon Et sum
    TH1D                    *h2PhotonPt;           //! di-photon Pt
    TH1D                    *h2PhotonMass;         //! di-photon mass (at least one trig)
    // high mass di-photons
    TH1D                    *hHmPhotonR91;         //! R9  1 photon
    TH1D                    *hHmPhotonR92;         //! R9  2 photon
    // R9 cleaned photons
    TH1D                    *h2R9PhotonMass;       //! di-photon mass (after R9 cleanup)
    // selected and triggered
    TH1D                    *hNTrigPhotons;        //! number of photons, selected and triggered
    TH1D                    *h2TrigPhotonMass;     //! di-photon mass (2 triggered photons)
    
    ClassDef(HttAnalysis, 1) // Higgs to tau tau analysis
  };
}
#endif
