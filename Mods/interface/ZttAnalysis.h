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
    void         SetMetName(const char *n)           { fMetName = n; }
    void         SetIgnoreElCharge(Bool_t b)         { fIgnoreElCharge = b;   }
    void         SetMaxZMass(Double_t m)             { fMaxZMass = m;         }
    void         SetMinDilMass(Double_t m)           { fDilMinMass = m;       }
    void         SetMinPt(Double_t pt)               { fMinPt      = pt;      }
    void         SetMinZMass(Double_t m)             { fMinZMass = m;         }
    
  protected:
    void                     Begin();
    void                     Process();
    void                     SlaveBegin();
    void                     SlaveTerminate();
    void                     Terminate();

    //----------------------------------------------------------------------------------------------
    // input collections
    TString      fCleanLeptonsName;     //clean leptons name (input)
    TString      fMetName;              //met name (input)
 
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
    TH1D        *fAllDiLepMass;         //!dilepton mass for all dilepton pairs
    TH1D        *fDiElMass;             //!dielectron mass 
    TH1D        *fDiMuMass;             //!dimuon mass
    TH1D        *fElMuMass;             //!electron-muon mass 
    TH1D        *fAllDiLepMassAcc;      //!accepted dilepton mass for all dilepton pairs
    TH1D        *fDiElMassAcc;          //!accepted dielectron mass 
    TH1D        *fDiMuMassAcc;          //!accepted dimuon mass
    TH1D        *fElMuMassAcc;          //!accepted electron-muon mass 
    TH1D        *fNLeptons;             //!number of leptons
    TH1D        *fNGPairs;              //!number of good pairs
    TH1D        *fNZPairs;              //!number of bad (Z) pairs

    TH1D        *fAllDiLepMet;          //!met for all dilepton pairs
    TH1D        *fAllDiLepPt1;          //!pt1 for all dilepton pairs
    TH1D        *fAllDiLepPt2;          //!pt2 for all dilepton pairs
    TH1D        *fAllDiLepDPhi;         //!dphi for all dilepton pairs
    TH1D        *fAllDiLepDEta;         //!deta for all dilepton pairs
    TH1D        *fAllDiLepDPhiMetOne;   //!dphi for all dilepton pairs between first lepton and met
    TH1D        *fAllDiLepDPhiMetTwo;   //!dphi for all dilepton pairs between second lepton and met
    TH1D        *fAllDiLepMtOne;        //!transverse mass for 1st lepton
    TH1D        *fAllDiLepMtTwo;        //!transverse mass for 2nd lepton

    TH1D        *fDiTauRecoMass;
    TH1D        *fDiTauTransverseMass;
    TH1D        *fDiTauTransverseEll;
    TH1D        *fDiTauTransverseEnn;
    TH1D        *fDiTauVisMass;
    TH1D        *fDiTauXTau1;
    TH1D        *fDiTauXTau2;         
    TH1D        *fEmuRecoMass;
    TH1D        *fEmuTransverseMass;
    TH1D        *fEmuTransverseEll;
    TH1D        *fEmuTransverseEnn;
    TH1D        *fEmuVisMass;
    TH1D        *fEmuXTau1;
    TH1D        *fEmuXTau2;         
    TH1D        *fMumuRecoMass;
    TH1D        *fMumuTransverseMass;
    TH1D        *fMumuTransverseEll;
    TH1D        *fMumuTransverseEnn;
    TH1D        *fMumuVisMass;
    TH1D        *fMumuXTau1;
    TH1D        *fMumuXTau2;         
    TH1D        *fEeRecoMass;
    TH1D        *fEeTransverseMass;
    TH1D        *fEeTransverseEll;
    TH1D        *fEeTransverseEnn;
    TH1D        *fEeVisMass;
    TH1D        *fEeXTau1;
    TH1D        *fEeXTau2;         
    
    //----------------------------------------------------------------------------------------------
    ClassDef(ZttAnalysis, 1) // Z to Tau Tau analysis
  };
}
#endif
