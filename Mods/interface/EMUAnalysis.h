//--------------------------------------------------------------------------------------------------
// $Id: $
//
// EMUAnalysis
//
// This module to analysis EMU events
//
// Authors: M.Klute
//--------------------------------------------------------------------------------------------------

#ifndef MITHTT_MODS_EMUANALYSIS_H
#define MITHTT_MODS_EMUANALYSIS_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"

class TH1D;
class TNtuple;
class MCEventInfo;

namespace mithep 
{
  class EMUAnalysis : public BaseMod
  {
  public:
    EMUAnalysis(const char *name  = "EMUAnalysis", 
		const char *title = "Z To Tau Tau To E MU Analysis");
    
    void         SetCleanLeptonsName(const char *n)  { fCleanLeptonsName = n; }
    void         SetTrigObjsName(const char *n)      { fTrigObjsName = n; }
    void         SetElecName(TString name)           { fElectronsName = name; }
    void         SetMuonName(TString name)           { fMuonsName = name; }
    void         SetJetName(TString name)            { fJetsName = name; }
    void         SetCaloJetName(TString name)        { fCaloJetsName = name; }
    void         SetMetName(const char *n)           { fMetName = n; }
    
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
    TString      fJetsName;             
    TString      fCaloJetsName;             
    TString      fMetName;              
    TString      fVertexName;     
    const MCEventInfo       *fMcEventInfo;
    const CaloJetCol        *fCaloJet;             
 
    //----------------------------------------------------------------------------------------------
    // module setup

    //----------------------------------------------------------------------------------------------
    // selection cuts
    Double_t     cutPtLeadingMuon;
    Double_t     cutPtLeadingElec;
    Double_t     cutPtSecondMuon;
    Double_t     cutPtSecondElec;
    Double_t     cutPtTriggerMuon;
    Double_t     cutPtTriggerElec;
    UInt_t       cutTriggerMuon;
    UInt_t       cutTriggerElec;
    //----------------------------------------------------------------------------------------------
    // histograms
    TH1D        *fNAccCounters;         //!history of cuts

    // lepton histograms
    TH1D        *fptmtrig;                 
    TH1D        *fptmnotrig;                 
    TH1D        *fptm;                 
    TH1D        *fpte;                 
    TH1D        *fetam;                
    TH1D        *fetae;                
    TH1D        *fphim;                
    TH1D        *fphie;                
    TH1D        *fdcam;                
    TH1D        *fdcae;                
    TH1D        *fchargem;                
    TH1D        *fchargee;                

    // di-tau histograms
    TH1D        *frecoMass;
    TH1D        *ftransMass;
    TH1D        *ftransEll;
    TH1D        *ftransEnn;
    TH1D        *fproj;
    TH1D        *fprojVis;
    TH1D        *fprojMet;
    TH1D        *fprojPhi;
    TH1D        *fhT;
    TH1D        *fvisMass;
    TH1D        *fscaledVisMass;
    TH1D        *fxtaum;
    TH1D        *fxtaue;    
    TH1D        *ftransMnmu;
    TH1D        *ftransMnel;
    TH1D        *fdphinmu;
    TH1D        *fdphinel;

    TH1D        *fdphi;                 
    TH1D        *fdeta;                 
    TH1D        *fcharge;                 
    TH1D        *ftype;                 
        
    // event histograms
    TH1D        *fnleptons;  
    TH1D        *fnmuons;  
    TH1D        *fnelecs;  
    TH1D        *fnjets;  
    TH1D        *fnjetspt15;  
    TH1D        *fnjetspt20;  
    TH1D        *fnbjetspt20sv2;
    TH1D        *fnbjetspt20sv05;
    TH1D        *fnjetspt25;  
    TH1D        *fmet;  
    TH1D        *fbtag;  
    TH1D        *feemass;
    TH1D        *fmmmass;
    TH1D        *fdzlepton;
    TH1D        *fdzelec;
    TH1D        *fdzmuon;

    // higgs selection
    TH1D        *fvisMassCut1;
    TH1D        *fvisMassCut2;
    TH1D        *fvisMassCut3;
    TH1D        *fvisMassCut4;
    TH1D        *fvisMassCut5;
    TH1D        *fvisMassCut6;
    TH1D        *fvisMassCut7;
    TH1D        *fvisMassCut8;
    
    TNtuple     *fNt;
    //----------------------------------------------------------------------------------------------
    ClassDef(EMUAnalysis, 1) // electron muon analysis
  };
}
#endif
