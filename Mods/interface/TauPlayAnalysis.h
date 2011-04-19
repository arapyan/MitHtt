//--------------------------------------------------------------------------------------------------
// $Id: TauPlayAnalysis.h,v 1.1 2011/02/18 14:56:34 dkralph Exp $
//
// TauPlayAnalysis
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
    
    void         SetTausName(TString name)            { fTausName = name;      }
    
  protected:
    void         SlaveBegin();
    void         Process();
    
    //----------------------------------------------------------------------------------------------
    // input collections
    TString      fTausName;
    const PFTauCol      *fTaus;
 
    // histograms
    TH1D *tauPt;
    TH1D *tauEta;

    //----------------------------------------------------------------------------------------------
    ClassDef(TauPlayAnalysis, 1) // Messing around with Taus
  };
}
#endif
