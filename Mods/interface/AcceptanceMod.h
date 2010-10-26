//--------------------------------------------------------------------------------------------------
// $Id: AcceptanceMod.h,v 1.3 2009/06/15 15:00:21 loizides Exp $
//
// AcceptanceMod
//
// This module calculates reconstruction efficiency (and fake rate) between reconstructed
// objects and the MC truth by simple matching in DeltaR.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITMODS_MODS_ACCEPTANCEMOD_H
#define MITMODS_MODS_ACCEPTANCEMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCParticle.h"
#include "MitAna/DataTree/interface/Particle.h"

class TH1D;

namespace mithep 
{
  class AcceptanceMod : public BaseMod
  {
    public:
      AcceptanceMod(const char *name="AcceptanceMod", 
             const char *title="Efficiency analysis module");
      
      void                     SetSuffix(const char* t) {fSuffix=TString(t);}
      const char*              AppendSuffix(const char* input) {return (TString(input)+fSuffix).Data();}
      void                     AddFinalStateInfo(const char* a, const char* b) 
	{string histName=string(a);
	string partType=string(b);
	std::pair<string,string> ab=make_pair(histName,partType);
	fParticlesVec->push_back(ab);
	}

    
      
	
    protected:
      void                     Process();
      void                     SlaveBegin();
      void                     Terminate();
      void                     AddAllowedFinalState(Int_t NElec,Int_t NMuon,Int_t NOther) {
	std::vector<Int_t> fs;
	fs.push_back(NElec);
	fs.push_back(NMuon);
	fs.push_back(NOther);
	fFinalStatesVec->push_back(fs);
	Int_t index=0+fFinalStatesVec->size();
	fFinalStatesNumer[index-1]=0;
	fFinalStatesDenom[index-1]=0;
	++fNFinalStates;
      }
      std::vector< std::pair<string,string> >* fParticlesVec;
      std::vector< std::vector<Int_t> >* fFinalStatesVec;
      Int_t*                   fFinalStatesDenom;
      Int_t*                   fFinalStatesNumer;
      Int_t                    fNFinalStates;
      Int_t                    fNEvents;
      TH1D*                    fAcceptanceHistByFinalState;
      TH1D*                    fNumeratorHistByFinalState;
      TH1D*                    fDenominatorHistByFinalState;
      TH1D*                    fNEventsHist;
      
      TString                  fSuffix;             //Suffix to histogram names
     
       

    ClassDef(AcceptanceMod, 1) // Efficiency analysis module
  };
}
#endif
