// $Id:$

#include "MitHtt/Mods/interface/MatchingMod.h"
#include "MitAna/DataCont/interface/BaseCollection.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TH1D.h>
#include <TH1I.h>

using namespace mithep;

ClassImp(mithep::MatchingMod)

//--------------------------------------------------------------------------------------------------
MatchingMod::MatchingMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMCParticleColName(TString("")), 
  fBarrelSCColName(TString("")),
  fEndcapSCColName(TString("")),
  fAllowedAbsEta(new std::vector<std::pair<Double_t,Double_t> >),
  fOutputCol(new mithep::ObjArray<MCParticle>),
  fOutputColName(TString("")),
  fEndcapSCCol(0),
  fBarrelSCCol(0),
  fMCParticleCol(0),
  fMinPt1(10),
  fMinEt1(2.4),
  fMinEt2(2.4),		 
  fRadius(0.1),
  fSuffix(TString("")),
  fPublicName(AppendSuffix("fAcceptHist")),
  kMatchToSC(kFALSE),
  kMatchToTrack(kFALSE) 
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
void MatchingMod::Process()
{
  // Process entries of the tree.

  LoadEventObject(fMCParticleColName,fMCParticleCol,kFALSE);
  if (kMatchToSC) {LoadBranch(fBarrelSCColName); LoadBranch(fEndcapSCColName);}

  TH1D* fAcceptHist= new TH1D(AppendSuffix("fAcceptHist"),"hAcceptHist",2,-0.5,1.5);
  
 
  UInt_t nMCParticles = 0;
  if (fMCParticleCol)
    nMCParticles = fMCParticleCol->GetEntries();

  // std::cout<<"Number of entries in input collection for matchingmod: " << nMCParticles << std::endl;

  UInt_t nBarrelSC = 0;
  if (fBarrelSCCol)
    nBarrelSC = fBarrelSCCol->GetEntries();

  UInt_t nEndcapSC=0;
  if (fEndcapSCCol)
    nEndcapSC = fEndcapSCCol->GetEntries();

  UInt_t nSC=nEndcapSC+nBarrelSC;
  Bool_t *found = new Bool_t[nSC];
  for (UInt_t i=0; i<nSC; ++i) 
    found[i] = kFALSE;
  
  // find matches
  for (UInt_t j=0; j<nMCParticles; ++j) {
    Bool_t kAccept=kTRUE;
    const MCParticle *p1 = fMCParticleCol->At(j);
    if (!p1)
      continue;
    Double_t absEta=p1->AbsEta();

    Double_t pt = p1->Pt();
   
    if (pt < fMinPt1) {
      kAccept=kFALSE;
      //std::cout<<"Reject!"<<std::endl;
      //std::cout<<"Pt: "<< pt << "is too low" <<std::endl;
      fAcceptHist->Fill(0);
      continue;
    }
    //Should checks to ensure that this eta is the one that we want for matching particles. There's the implication of a vector based on the momentum of the generated particle, which may or may not be the same as the location; depends on the expediency of the decay
    
    Bool_t kAllowedEta=kFALSE;
    for (UInt_t etaInd=0;etaInd<fAllowedAbsEta->size();++etaInd) {
      std::pair<Double_t,Double_t> currPair=fAllowedAbsEta->at(etaInd);
      Double_t etaMin=currPair.first;
      Double_t etaMax=currPair.second;
      if (absEta>etaMin && absEta<etaMax) kAllowedEta=kTRUE;
    }
    
    if (!kAllowedEta) {kAccept=kFALSE; fAcceptHist->Fill(0); //std::cout<<"Reject!"<<std::endl;
      continue;}

    Double_t phi = p1->Phi();
    Double_t eta = p1->Eta();

    if (kMatchToSC) {
      UInt_t foundInd  = nSC;
      Double_t foundD  = 1e12;
      
      const SuperClusterCol *fSCCol;
      UInt_t scCandidates;
      UInt_t initIndex;UInt_t finalIndex;
      if (absEta<1.5) {fSCCol=fBarrelSCCol; scCandidates=nSC;initIndex=0;finalIndex=nBarrelSC;}
      else {fSCCol=fEndcapSCCol;scCandidates=nSC;initIndex=nBarrelSC;finalIndex=nSC;}
      
      for (UInt_t k=initIndex; k<finalIndex; ++k) {
	if (found[k]) continue;
	const SuperCluster *sc1 = fSCCol->At(k-initIndex);
	if (!sc1) //Ensure sc1 is not null
	  continue;
	if (!sc1->Et()>fMinEt2) //Ensure sc1 meets kinematic requirements
	  continue;

	if(MathUtils::DeltaR(phi, eta, sc1->Phi(), sc1->Eta()) < fRadius) {
	  Double_t newDiff = TMath::Abs(pt - sc1->Et());
	  if (newDiff < foundD) {
	    foundInd = k;
	    foundD = newDiff;
	  } 
	}
      }
      
      if (foundInd < finalIndex) {
	found[foundInd] = 1;
	fOutputCol->Add(fMCParticleCol->At(j));
	kAccept=kTRUE;
      }
      else {
	kAccept=kFALSE;
	//std::cout<<"Reject!"<<std::endl;
	//std::cout<<"No SC Match found"<<std::endl;
	fAcceptHist->Fill(0);
	continue;
      }
    }
    if (kMatchToTrack) {//std::cout<<"matchtotrack"<<std::endl
      ;}
    if (kAccept) {
	fAcceptHist->Fill(1);
	//std::cout<<"Accept!"<<std::endl;
 }
    else {fAcceptHist->Fill(0); //std::cout<<"Reject!"<<std::endl;
    }
  }

  // std::cout<<"AddingObjThisEvt"<<std::endl;
  AddObjThisEvt(fAcceptHist,AppendSuffix("fAcceptHist"));

}

//--------------------------------------------------------------------------------------------------
void MatchingMod::SlaveBegin()
{
  // Create and add histograms to the output list.
  ReqBranch(fBarrelSCColName ,fBarrelSCCol);
  ReqBranch(fEndcapSCColName ,fEndcapSCCol);
    
}

//--------------------------------------------------------------------------------------------------
void MatchingMod::Terminate()
{
  // Create and add ratio histograms to the output list.

  
}
