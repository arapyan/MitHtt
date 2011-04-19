//--------------------------------------------------------------------------------------------------
//
//
//
//--------------------------------------------------------------------------------------------------

#ifndef MITHTT_MODS_EMUANALYSIS2_H
#define MITHTT_MODS_EMUANALYSIS2_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitHtt/Mods/interface/PlotInfo.h"

class TH1D;
class MCEventInfo;
class TTree;

namespace mithep 
{
  class EMUAnalysis2 : public BaseMod
  {
  public:
    EMUAnalysis2(const char *name  = "EMUAnalysis2", 
		const char *title = "tau tau to e mu analysis");
    ~EMUAnalysis2();
    
    void         SetMuonName(TString s)         { fMuonsName     	= s; }
    void         SetElecName(TString s)         { fElectronsName 	= s; }
    void         SetTrigObjsName(TString s)     { fTrigObjsName  	= s; }
    void         SetJetName(TString s)          { fJetsName      	= s; }
    void         SetMetName(TString s)          { fMetsName      	= s; }
    void         SetHistNamePref(TString s) 	{ fHistNamePref 	= s; }

    void         Fill(const char *hname, double val);
    void         Fill(const char *hname, UInt_t val)   { Fill(hname,double(val)); }
    void         Fill(const char *hname, int val)      { Fill(hname,double(val)); }
    
  protected:
    void         Process();
    void         SlaveBegin();
    // find out if event passes all but the named cut
    bool         NMinusOnePass(const char *cut);
    bool         NMinusThreePass(const char *cut1,const char *cut2,const char *cut3);
    
    //----------------------------------------------------------------------------------------------
    // input collections
    TString      	 fTrigObjsName;         
    TString     	 fMuonsName;            
    TString  		 fElectronsName;        
    TString     	 fJetsName;             
    TString      	 fCaloJetsName;             
    TString      	 fMetsName;
    const MuonCol     	*fMuons;
    const ElectronCol   *fElectrons;
    const JetCol        *fJets;
    const CaloJetCol    *fCaloJets;
    const MetCol        *fMets;
    const MCEventInfo   *fMcEventInfo;
 
/*     TTree * fTr; */
/*     Double_t fLepPt,fTauPt,fLepEta,fTauEta,fMet,fDPhi,fMll,fMt,fD0,fm,fisoL,fisoM,fisoT; */
/*     UInt_t fNJets,fNVCJets; */

    //----------------------------------------------------------------------------------------------
    // selection cuts
    Double_t     cutPtLeadingMuon;
    Double_t     cutPtLeadingElec;
    Double_t     cutPtSecondMuon;
    Double_t     cutPtSecondElec;
    Double_t     cutPtTriggerMuon;
    //----------------------------------------------------------------------------------------------
    // plot info
    UInt_t fNHists;
    vector<PlotInfo*> fInfoV;
    vector<TH1D*> fHistV;
    map<string,int> fHistMap;
    TString fHistNamePref;      // prefix for histogram names

    // cut info
    UInt_t fNCuts;
    vector<string> fCutNames;
    map<string,UInt_t> fCutMap;
    vector<bool> fPassCuts;

    //----------------------------------------------------------------------------------------------
    ClassDef(EMUAnalysis2, 1) //
  };
}
#endif
