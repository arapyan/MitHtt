//--------------------------------------------------------------------------------------------------
// $Id: $
//
// TauPlayAnalysis
//
// This module doesn't do much yet
//
// Authors:
//--------------------------------------------------------------------------------------------------

#ifndef MITHTT_MODS_TAUPLAYANALYSIS_H
#define MITHTT_MODS_TAUPLAYANALYSIS_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitHtt/Mods/interface/PlotInfo.h"

class TH1D;
class MCEventInfo;

namespace mithep 
{
  class TauPlayAnalysis : public BaseMod
  {
  public:
    TauPlayAnalysis(const char *name  = "TauPlayAnalysis", 
		const char *title = "Messing around with taus");
    ~TauPlayAnalysis();
    
    void         SetCleanLeptonsName(const char *n)  { fCleanLeptonsName = n; }
    void         SetTrigObjsName(const char *n)      { fTrigObjsName = n;     }
    void         SetElecName(TString name)           { fElectronsName = name; }
    void         SetMuonName(TString name)           { fMuonsName = name;     }
    void         SetTauName(TString name)            { fTausName = name;      }
    void         SetHPSTauCandidatesName(TString name)         { fHPSTauCandidatesName = name; }
    void         SetMetName(const char *n)           { fMetColName = n;          }
    void         SetJetName(const char *n)           { fJetsName = n;          }
    
  protected:
    void         Begin();
    void         Process();
    void         SlaveBegin();
    void         SlaveTerminate();
    void         Terminate();
    
    //----------------------------------------------------------------------------------------------
    // input collections
    TString      fCleanLeptonsName;     
    TString      fTrigObjsName;         
    TString      fMuonsName;            
    TString      fElectronsName;
    TString      fTausName;
    TString      fHPSTauCandidatesName;
    TString      fMetColName;
    TString      fJetsName;
    TString              fMCPartName;         //name of MCParticle branch                                             
    const MCEventInfo   *fMcEventInfo;
    const ParticleCol   *fLeptons;
    const PFTauCol      *fTaus;
    const PFTauCol      *fHPSTauCandidates;
    const MCParticleCol *fParticles;          //!MCParticle branch                                                    
 
    //----------------------------------------------------------------------------------------------
    // module setup
    vector<PlotInfo*> finfov;
    vector<TH1D*> fHistsv;
    UInt_t fnhists;
    map<string,int> fmap;
    //----------------------------------------------------------------------------------------------
    // selection cuts
    //----------------------------------------------------------------------------------------------
    // histograms
    TH1D        *fNAccCounters;         //!history of cuts

    //----------------------------------------------------------------------------------------------

    //----------------------------------------------------------------------------------------------
    ClassDef(TauPlayAnalysis, 1) // Messing around with Taus
  };
}
#endif
