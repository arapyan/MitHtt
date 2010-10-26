//--------------------------------------------------------------------------------------------------
// $Id:$
//
// MatchingMod
//
// This module matches objects in Col1 to objects in a SuperCluster collection and outputs a collection of matched objects. Output collection to be used in an event-level acceptance calculation.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITMODS_MODS_MATCHINGMOD_H
#define MITMODS_MODS_MATCHINGMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/MCParticle.h"
#include "MitAna/DataTree/interface/SuperCluster.h"
#include "MitAna/DataTree/interface/ParticleFwd.h"
#include "MitAna/DataTree/interface/SuperClusterFwd.h"
#include "MitAna/DataTree/interface/MCParticleFwd.h"

class TH1D;

namespace mithep 
{ 
  class MatchingMod : public BaseMod
  {
    public:
      MatchingMod(const char *name="MatchingMod", 
             const char *title="Efficiency analysis module");

      void                     SetMCParticleColName(const char *n)       { fMCParticleColName = n;  }
      void                     SetEndcapSCColName(const char *n)       { fEndcapSCColName = n;  }
      void                     SetBarrelSCColName(const char* n) { fBarrelSCColName=n;}
      void                     SetMinPt1(Double_t pt)            { fMinPt1    = pt; }
      // void                     SetMinEt1(Double_t e)            { fMinEt1   = e; }
      
      void                     SetMinEt2(Double_t et) {fMinEt2=et;}

      void                     AddAllowedAbsEta(Double_t a,Double_t b ) { fAllowedAbsEta->push_back(make_pair(a,b));  }
 
      void                     SetRadius(Double_t r)            { fRadius = r; }
      void                     SetSuffix(const char* t) {fSuffix=TString(t);}
      void                     SetMatchToSC(Bool_t a) {kMatchToSC=a;}
      void                     SetMatchToTrack(Bool_t b) {kMatchToTrack=b;}
      const char*              AppendSuffix(const char* input) {return (TString(input)+fSuffix).Data();}
      const char*              GetPublicName() {return fPublicName.Data();}
      void                     SetPublicName() {fPublicName=AppendSuffix("fAcceptHist");}
    protected:
      void                     Process();
      void                     SlaveBegin();
      void                     Terminate();

      const mithep::SuperClusterCol* fBarrelSCCol;
      const mithep::SuperClusterCol* fEndcapSCCol;
      const mithep::MCParticleOArr*  fMCParticleCol;
      mithep::MCParticleOArr*    fOutputCol;

      TString                  fMCParticleColName;           //first  collection name (input)
      TString                  fEndcapSCColName;           //second collection name (input)
      TString                  fBarrelSCColName;
      TString                  fOutputColName;    //output collection name
      TString                  fSuffix;             //Suffix to histogram names
      TString                  fPublicName;
      Double_t                 fMinPt1;              //minimum pt for MCParticleCol
      Double_t                 fMinEt1;             //minimum Et for MCParticleCol
      Double_t                 fMinEt2;             //minimum Et for SCCol
      std::vector<std::pair<Double_t,Double_t> >* fAllowedAbsEta; //Allowed Abs Eta ranges for col1
      Double_t                 fRadius;             //radius used for matching
      Bool_t kMatchToSC;      //flags for matching to supercluster
      Bool_t kMatchToTrack;   //flag for matching to track
      ClassDef(MatchingMod, 1) // Efficiency analysis module
  };
}
#endif
