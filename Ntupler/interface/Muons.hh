#ifndef MITHTT_NTUPLER_MUONS_HH

#include "MitAna/DataTree/CollectionsFwd.h"
#include "MitPhysics/Utils/interface/HLTTool.h"

namespace mithep {
  class Muons 
  {
  public:
    Muons()  {fPtMin = 0; fEtaMax=0; } 
    ~Muons() {delete fMuonArr; }
    //Filler Functions
    void setupTreeBranch(TTree *iTree,double iPtMin=10.,double iEtaMax=2.5);
    void fillMuons(const muonCol *iMuons, const Vertex *iPV,const PFCandidateCol* iPFNoPU,const PFCandidateCol* iPFPU,const PFCandidates *iCands);
    void fillMuon (const Muon    *iMu,    const Vertex *iPV,const PFCandidateCol *iPFNoPU,const PFCandidateCol *iPFPU,const PFCandidates *iCands);
    //Helpers
    bool isPFMuon(const Muon *iMu,const PFCandidates *iCands);
  
  private:
    double fPtMin;
    dobule fEtaMax;
    TClonseArray *fMuonArr;
    HLTTool      *fHLT;
  };
}
