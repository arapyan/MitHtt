//--------------------------------------------------------------------------------------------------
// $Id: $
//
// ZttAnalysis
//
// This module to analysis Z to tau tau events
//
// Authors: M.Klute
//--------------------------------------------------------------------------------------------------

#ifndef MITHTT_MODS_ZTTANALYSIS_H
#define MITHTT_MODS_ZTTANALYSIS_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"

class TH1D;
class MCEventInfo;

namespace mithep 
{
  class ZttAnalysis : public BaseMod
  {
  public:
    ZttAnalysis(const char *name  = "ZttAnalysis", 
		const char *title = "Z To Tau Tau Analysis");
    
    void         SetCleanLeptonsName(const char *n)  { fCleanLeptonsName = n; }
    void         SetMetName(const char *n)           { fMetName = n;          }
    void         SetIgnoreElCharge(Bool_t b)         { fIgnoreElCharge = b;   }
    void         SetMaxZMass(Double_t m)             { fMaxZMass = m;         }
    void         SetMinDilMass(Double_t m)           { fDilMinMass = m;       }
    void         SetMinPt(Double_t pt)               { fMinPt      = pt;      }
    void         SetMinZMass(Double_t m)             { fMinZMass = m;         }
    void         SetElecName(TString name)           { fElecName = name; }
    void         SetMuonName(TString name)           { fMuonName = name; }
    void         SetTrigObjsName(const char *name)   { fTrigObjsName = name; }
    
  protected:
    void         Begin();
    void         Process();
    void         SlaveBegin();
    void         SlaveTerminate();
    void         Terminate();
    
    void         MatchCollections(Bool_t *, const char *, const char *);
    Bool_t       *MatchTriggerCollection(const char *, const char *);
    //----------------------------------------------------------------------------------------------
    // input collections
    TString      fCleanLeptonsName;     //clean leptons name (input)
    TString      fMetName;              //met name (input)
    TString      fTrigObjsName;         // name of trigger objects
    TString                  fMuonName;                 // name of muon collection
    const MuonCol           *fMuons;                    //! muon from data stream
    TString                  fElecName;                 // name of electron collection
    const ElectronCol       *fElecs;                    //! elec from data stream
    const MCEventInfo       *fMcEventInfo;              //! MC event information branch
 
    //----------------------------------------------------------------------------------------------
    // module setup
    TString      fFinalState;           // used to set two leptons and name of histograms
    UInt_t        fFinalStateId;        //       
    //----------------------------------------------------------------------------------------------
    // selection cuts
    Double_t     fMinPt;                //minimum pt for leptons
    Double_t     fDilMinMass;           //minimum dilepton mass
    Double_t     fMinZMass;             //minimum Z mass
    Double_t     fMaxZMass;             //maximum Z mass
    Bool_t       fIgnoreElCharge;       //=true then ignore electron charge for z mass cut

    //----------------------------------------------------------------------------------------------
    // histograms
    TH1D        *fNAccCounters;         //!history of cuts
    TH1D        *fLLMass;               //!dilepton mass for all dilepton pairs
    TH1D        *fNLeptons;             //!number of leptons
    TH1D        *fNMuons;         
    TH1D        *fNElecs;         
    TH1D        *fNAccMuons;      
    TH1D        *fNAccElecs;      
    TH1D        *fAllMuonPt;
    TH1D        *fAllMuonEta;
    TH1D        *fAllElecPt;
    TH1D        *fAllElecEta;

    TH1D        *fNVertex;              //!number of vertices
    TH1D        *fNGPairs;              //!number of good pairs
    TH1D        *fNZPairs;              //!number of bad (Z) pairs

    TH1D        *fmet;                  //!met for all dilepton pairs
    TH1D        *fpt1;                  //!pt1 for all dilepton pairs
    TH1D        *fpt2;                  //!pt2 for all dilepton pairs
    TH1D        *feta1;                 //!eta1 for all dilepton pairs
    TH1D        *feta2;                 //!eta2 for all dilepton pairs
    TH1D        *fphi1;                 //!eta1 for all dilepton pairs
    TH1D        *fphi2;                 //!eta2 for all dilepton pairs
    TH1D        *fdphi;                 //!dphi for all dilepton pairs
    TH1D        *fdeta;                 //!deta for all dilepton pairs
    TH1D        *fdphiMet1;             //!dphi for all dilepton pairs between first lepton and met
    TH1D        *fdphiMet2;             //!dphi for all dilepton pairs between second lepton and met
    TH1D        *fmt1;                  //!transverse mass for 1st lepton
    TH1D        *fmt2;                  //!transverse mass for 2nd lepton

    TH1D        *frecoMass;
    TH1D        *ftransMass;
    TH1D        *ftransEll;
    TH1D        *ftransEnn;
    TH1D        *fvisMass;
    TH1D        *fxtau1;
    TH1D        *fxtau2;         
    //----------------------------------------------------------------------------------------------
    ClassDef(ZttAnalysis, 1) // Z to Tau Tau analysis
  };
}
#endif
