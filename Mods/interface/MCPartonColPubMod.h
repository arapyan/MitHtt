//--------------------------------------------------------------------------------------------------
// $Id:$
//
// MCPartonColPubMod
//
// This module calculates reconstruction efficiency (and fake rate) between reconstructed
// objects and the MC truth by simple matching in DeltaR.
//
// Authors: Z.Hynes
//--------------------------------------------------------------------------------------------------

#ifndef MITMODS_MODS_MCPARTONCOLPUBMOD_H
#define MITMODS_MODS_MCPARTONCOLPUBMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCParticle.h"

//testing...

class TH1I;

namespace mithep 
{
  class MCPartonColPubMod : public BaseMod
  {
    public:
      MCPartonColPubMod(const char *name="MCPartonColPubMod", 
             const char *title="MCPartonColPub module");

      void                     SetMCParticleColName(const char *n)       { fMCParticleColName = n;  }
      void                     SetOutputName(const char *n)       { fOutputColName = n;  }
      void                     SetDaughterPdgId(UInt_t t) { fDaughterPdgId=t; }
      void                     SetMotherPdgId(Int_t t) { fMotherPdgId=t; }
      void                     AddGaugeBosonIntermediatePdgId(Int_t a) { fIntermediates->push_back(a) ;}
      Int_t                    GetMotherPdgId() {return fMotherPdgId;}
      UInt_t                   GetDaughterPdgId() {return fDaughterPdgId;}

    protected:
      void                     Process();
      void                     Begin();
      void                     SlaveBegin();
      void                     SlaveTerminate();
      void                     Terminate();
      const char*              GetMCParticleColName() {return fMCParticleColName.Data();}
      const char*              GetOutputName() {return fOutputColName.Data();}
    
      TString                  fMCParticleColName;           //first  collection name (input)
      TString                  fOutputColName;           //second collection name (input)
      Bool_t                   fPubPerEvent;

      const Collection<MCParticle> *fMCParticleCol;
      ObjArray<MCParticle>*  fOutputCol;
      Int_t   fDaughterPdgId;           //particle type (default = kNone)
      Int_t fMotherPdgId;
      std::vector<Int_t>*      fIntermediates;  
    ClassDef(MCPartonColPubMod, 1) // Collection Publisher Module
  };
}
#endif
