//--------------------------------------------------------------------------------------------------
// $Id: $
//
// ZmmAnalysis
//
// This module uses Muons to look for Z->mm events and higher mass resonances.
//
// Authors: C.Paus
//--------------------------------------------------------------------------------------------------

#ifndef MITHGG_MODS_ZMMANALYSIS_H
#define MITHGG_MODS_ZMMANALYSIS_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class TNtuple;
class MCEventInfo;

namespace mithep 
{
  class ZmmAnalysis : public BaseMod
  {
  public:
    ZmmAnalysis(const char *name  = "ZmmAnalysis", 
		const char *title = "Z To Muon Muon Analysis");
    
  protected:
    void                     Begin();
    void                     Process();
    void                     SlaveBegin();
    void                     SlaveTerminate();
    void                     Terminate();

    void                     SetTrigObjsName(const char *name) { fTrigObjsName = name; }
    void                     SetMuonName(TString name)         { fMuonName = name; }
    void                     SetMuonsFromBranch(Bool_t b)      { fMuonsFromBranch = b; }
    void                     SetIsData(Bool_t b)               { fIsData = b; }
    void                     SetSkimming(Bool_t b)             { fSkimming = b; }

    void                     ResetSelectedMuons();      // reset stored muon indices
    void                     MatchMuonsToTrigger();     // create trig muon indices from sel.
    Bool_t                   MuonTriggered(UInt_t i);   // did this muon trigger
    
    void                     SetOverlapCut(double c)           { fOverlapCut = c; }

    TString                  fTrigObjsName;             // name of trigger objects
    TString                  fMuonName;                 // name of muon collection
    Bool_t                   fMuonsFromBranch;          // muons are loaded from a branch
    TString                  fMcEventInfoName;          // name of MC Event information
    const MuonCol           *fMuons;                    //! muon from data stream
    const MCEventInfo       *fMcEventInfo;              //! MC event information branch
    Bool_t                   fIsData;                   // looking at Data (or MC)
    Bool_t                   fSkimming;                 // are we skimming (skipping event)

    // general selection cuts
    double                   fOverlapCut;               // cut to allow for rejection of overlap
    // muon selection cuts
    double                   fMuonEtMin;                // minimal transverse energy for muon
    //Bool_t                   fNoPixelSeedCut;           // request that there be no pixel seed?
    //double                   fMuonSigEeBarMax;          // maximal sigma(eta,eta), barrel
    //double                   fMuonSigEeEcpMax;          // maximal sigma(eta,eta), endcap
    double                   fMuonIsoMax;               // maximal combined isolation
    // high mass di-muon selection cuts
    double                   f2MuonMassMin;             // minimum di-muon mass, high mass sample
    //// R9 cleaned muon selection cuts
    //double                   fMuonR9Min;              // minimum R9 for cleaner muons
 
    // keep track of our muons
    std::vector<UInt_t>      fSelectedMuons;            // muon indices, selected
    std::vector<UInt_t>      fTriggeredMuons;           // muon indices, selected and triggered

    //----------------------------------------------------------------------------------------------
    // histograms
    //----------------------------------------------------------------------------------------------
    // all muons
    TH1D                    *hAllNMuons;         //! all number of muons
    TH1D                    *hAllMuonPt;         //! all muon pt
    TH1D                    *hAllMuonEta;        //! all muon eta
    //TH1D                    *hAllMuonPixelSeed;  //! pixel seed or not
    //TH1D                    *hAllMuonSigEeBar;   //! sigma(eta,eta) of muons, barrel
    //TH1D                    *hAllMuonSigEeEcp;   //! sigma(eta,eta) of muons, endcap
    TH1D                    *hAllMuonEcHcTrIso;  //! combined isolation
    // trigger monitoring
    TH1D                    *hTrigDeltaEta;      //! trigger match (delta Eta)
    TH1D                    *hTrigDeltaPhi;      //! trigger match (delta Phi)
    TH1D                    *hTrigDeltaEt;       //! trigger match (delta Et)
    // selected muons and di-muons
    TH1D                    *h2SeleMuonMass;     //! di-muon mass (2 selected muons)
    TH1D                    *hNMuons;            //! number of muons, selected
    TH1D                    *hMuonEtMx;          //! max Et muon
    TH1D                    *hMuonSigEtaEta;     //! sigma_eta,eta for muon
    TH1D                    *hMuonEcHcTrIso;     //! combined ECAL, HCAL, Tracker isolation
    TH1D                    *hN2Muons;           //! number of di-muons
    TH1D                    *hMuonEt1;           //! Et 1 muon
    TH1D                    *hMuonEt2;           //! Et 2 muon
    TH1D                    *hMuonEta1;          //! Eta 1 muon
    TH1D                    *hMuonEta2;          //! Eta 2 muon
    TH1D                    *hMuonPhi1;          //! Phi 1 muon
    TH1D                    *hMuonPhi2;          //! Phi 2 muon
    //TH1D                    *hMuonR91;           //! R9  1 muon
    //TH1D                    *hMuonR92;           //! R9  2 muon
    TH1D                    *hMuonDelR;          //! delta R(muon1,muon2)
    TH1D                    *h2MuonEtSum;        //! di-muon Et sum
    TH1D                    *h2MuonPt;           //! di-muon Pt
    TH1D                    *h2MuonMass;         //! di-muon mass (at least one trig)
    //// high mass di-muons
    //TH1D                    *hHmMuonR91;         //! R9  1 muon
    //TH1D                    *hHmMuonR92;         //! R9  2 muon
    // R9 cleaned muons
    TH1D                    *h2R9MuonMass;       //! di-muon mass (after R9 cleanup)
    // selected and triggered
    TH1D                    *hNTrigMuons;        //! number of muons, selected and triggered
    TH1D                    *h2TrigMuonMass;     //! di-muon mass (2 triggered muons)
    
    TNtuple     *fNt;
    ClassDef(ZmmAnalysis, 1) // Z to Muon Muon analysis
  };
}
#endif
