//--------------------------------------------------------------------------------------------------
// $Id: $
//
// ZeeAnalysis
//
// This module uses Electrons to look for Z->ee events and higher mass resonances.
//
// Authors: C.Paus
//--------------------------------------------------------------------------------------------------

#ifndef MITHGG_MODS_ZEEANALYSIS_H
#define MITHGG_MODS_ZEEANALYSIS_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class MCEventInfo;

namespace mithep 
{
  class ZeeAnalysis : public BaseMod
  {
  public:
    ZeeAnalysis(const char *name  = "ZeeAnalysis", 
		const char *title = "Z To Electron Electron Analysis");
    
  protected:
    void                     Begin();
    void                     Process();
    void                     SlaveBegin();
    void                     SlaveTerminate();
    void                     Terminate();

    void                     SetTrigObjsName(const char *name) { fTrigObjsName = name; }
    void                     SetElecName(TString name)         { fElecName = name; }
    void                     SetElecsFromBranch(Bool_t b)      { fElecsFromBranch = b; }
    void                     SetIsData(Bool_t b)               { fIsData = b; }
    void                     SetSkimming(Bool_t b)             { fSkimming = b; }

    void                     ResetSelectedElecs();      // reset stored elec indices
    void                     MatchElecsToTrigger();     // create trig elec indices from sel.
    Bool_t                   ElecTriggered(UInt_t i);   // did this elec trigger
    
    void                     SetOverlapCut(double c)           { fOverlapCut = c; }

    TString                  fTrigObjsName;             // name of trigger objects
    TString                  fElecName;                 // name of elec collection
    Bool_t                   fElecsFromBranch;          // elecs are loaded from a branch
    TString                  fMcEventInfoName;          // name of MC Event information
    const ElectronCol       *fElecs;                    //! elec from data stream
    const MCEventInfo       *fMcEventInfo;              //! MC event information branch
    Bool_t                   fIsData;                   // looking at Data (or MC)
    Bool_t                   fSkimming;                 // are we skimming (skipping event)

    // general selection cuts
    double                   fOverlapCut;               // cut to allow for rejection of overlap
    // elec selection cuts
    double                   fElecEtMin;                // minimal transverse energy for elec
    //Bool_t                   fNoPixelSeedCut;           // request that there be no pixel seed?
    double                   fElecSigEeBarMax;          // maximal sigma(eta,eta), barrel
    double                   fElecSigEeEcpMax;          // maximal sigma(eta,eta), endcap
    double                   fElecIsoMax;               // maximal combined isolation
    // high mass di-elec selection cuts
    double                   f2ElecMassMin;             // minimum di-elec mass, high mass sample
    //// R9 cleaned elec selection cuts
    //double                   fElecR9Min;              // minimum R9 for cleaner elecs
 
    // keep track of our elecs
    std::vector<UInt_t>      fSelectedElecs;            // elec indices, selected
    std::vector<UInt_t>      fTriggeredElecs;           // elec indices, selected and triggered

    //----------------------------------------------------------------------------------------------
    // histograms
    //----------------------------------------------------------------------------------------------
    // all elecs
    TH1D                    *hAllNElecs;         //! all number of elecs
    TH1D                    *hAllElecPt;         //! all elec pt
    TH1D                    *hAllElecEta;        //! all elec eta
    //TH1D                    *hAllElecPixelSeed;  //! pixel seed or not
    TH1D                    *hAllElecSigEeBar;   //! sigma(eta,eta) of elecs, barrel
    TH1D                    *hAllElecSigEeEcp;   //! sigma(eta,eta) of elecs, endcap
    TH1D                    *hAllElecEcHcTrIso;  //! combined isolation
    // trigger monitoring
    TH1D                    *hTrigDeltaEta;      //! trigger match (delta Eta)
    TH1D                    *hTrigDeltaPhi;      //! trigger match (delta Phi)
    TH1D                    *hTrigDeltaEt;       //! trigger match (delta Et)
    // selected elecs and di-elecs
    TH1D                    *h2SeleElecMass;     //! di-elec mass (2 selected elecs)
    TH1D                    *hNElecs;            //! number of elecs, selected
    TH1D                    *hElecEtMx;          //! max Et elec
    TH1D                    *hElecSigEtaEta;     //! sigma_eta,eta for elec
    TH1D                    *hElecEcHcTrIso;     //! combined ECAL, HCAL, Tracker isolation
    TH1D                    *hN2Elecs;           //! number of di-elecs
    TH1D                    *hElecEt1;           //! Et 1 elec
    TH1D                    *hElecEt2;           //! Et 2 elec
    TH1D                    *hElecEta1;          //! Eta 1 elec
    TH1D                    *hElecEta2;          //! Eta 2 elec
    TH1D                    *hElecPhi1;          //! Phi 1 elec
    TH1D                    *hElecPhi2;          //! Phi 2 elec
    //TH1D                    *hElecR91;           //! R9  1 elec
    //TH1D                    *hElecR92;           //! R9  2 elec
    TH1D                    *hElecDelR;          //! delta R(elec1,elec2)
    TH1D                    *h2ElecEtSum;        //! di-elec Et sum
    TH1D                    *h2ElecPt;           //! di-elec Pt
    TH1D                    *h2ElecMass;         //! di-elec mass (at least one trig)
    //// high mass di-elecs
    //TH1D                    *hHmElecR91;         //! R9  1 elec
    //TH1D                    *hHmElecR92;         //! R9  2 elec
    // R9 cleaned elecs
    TH1D                    *h2R9ElecMass;       //! di-elec mass (after R9 cleanup)
    // selected and triggered
    TH1D                    *hNTrigElecs;        //! number of elecs, selected and triggered
    TH1D                    *h2TrigElecMass;     //! di-elec mass (2 triggered elecs)
    
    ClassDef(ZeeAnalysis, 1) // Z to Electron Electron analysis
  };
}
#endif
