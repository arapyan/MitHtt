// $Id:$

#include "MitHtt/Mods/interface/AcceptanceMod.h"
#include "MitAna/DataCont/interface/BaseCollection.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/TreeMod/interface/BaseMod.h"
#include <TH1D.h>

using namespace mithep;

ClassImp(mithep::AcceptanceMod)

//--------------------------------------------------------------------------------------------------
AcceptanceMod::AcceptanceMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fParticlesVec(new std::vector< std::pair<string,string> >),
  fFinalStatesVec(new std::vector< std::vector<Int_t> >),
  fFinalStatesNumer(new Int_t[20]),
  fFinalStatesDenom(new Int_t[20]),
  fNFinalStates(0)
{
  // Constructor.
}


//--------------------------------------------------------------------------------------------------
void AcceptanceMod::Process()
{
  // Process entries of the tree.
  fNEventsHist->Fill(0);
  TH1D* elecHist=new TH1D("elecHist","hElecHist",2,-0.5,1.5);
  TH1D* muonHist=new TH1D("muonHist","hMuonHist",2,-0.5,1.5);
  TH1D* otherHist=new TH1D("otherHist","hOtherHist",2,-0.5,1.5);

  for (UInt_t i=0;i<fParticlesVec->size();++i) {
    const TH1D* currHist=GetObjThisEvt<TH1D>(((fParticlesVec->at(i)).first).c_str());
    //    std::cout<<"Entries in current histogram: "<<currHist->GetEntries()<<std::endl;
    if (!currHist) continue;
    if (TString(fParticlesVec->at(i).second)==string("Elec"))
      elecHist->Add(currHist);
    else {
      if (TString(fParticlesVec->at(i).second)==string("Muon"))
      muonHist->Add(currHist);
      
      else otherHist->Add(currHist);
    }
  }
  if (elecHist->GetEntries()>0) {
  std::cout<<"accepted electrons: "<<elecHist->GetBinContent(2)<<std::endl;
  std::cout<<"rejected electrons: "<<elecHist->GetBinContent(1)<<std::endl;
  }

  if (muonHist->GetEntries()>0) {
  std::cout<<"accepted muons: "<<muonHist->GetBinContent(2)<<std::endl;
  std::cout<<"rejected muons: "<<muonHist->GetBinContent(1)<<std::endl;
  }
  Bool_t kReject=kFALSE;
  Int_t NElec=elecHist->GetEntries();
  if (elecHist->GetBinContent(1)>0) kReject=kTRUE;      
  Int_t NMuon=muonHist->GetEntries();
  if (muonHist->GetBinContent(1)>0) kReject=kTRUE;
  Int_t NOther=otherHist->GetEntries();
  //if (otherHist->GetBinContent(1)>(10^(-10))) kReject=kTRUE;

  //std::cout<<"kReject: "<<kReject<<std::endl;

  //Int_t finalStateIndex;
  for (UInt_t fsi=0;fsi<fFinalStatesVec->size();++fsi) {
    if (NElec==(fFinalStatesVec->at(fsi)).at(0) and NMuon==fFinalStatesVec->at(fsi).at(1) and NOther==(fFinalStatesVec->at(fsi)).at(2)) {
      if (!kReject) {
        Int_t currentNumerator=fFinalStatesNumer[fsi]+1;
	fFinalStatesNumer[fsi]=currentNumerator;
	fNumeratorHistByFinalState->Fill(fsi);
      }
      Int_t currentDenominator=fFinalStatesDenom[fsi]+1;
      fFinalStatesDenom[fsi]=currentDenominator;
      fDenominatorHistByFinalState->Fill(fsi);
    }
  }
  delete elecHist;
  delete muonHist;
  delete otherHist;
}


//--------------------------------------------------------------------------------------------------
void AcceptanceMod::SlaveBegin()
{
  // Create and add histograms to the output list.
  
  AddTH1(fAcceptanceHistByFinalState,AppendSuffix("hAcceptanceHist"),";Final State Index;Acceptance",6,-0.5,5.5);
  AddTH1(fDenominatorHistByFinalState,AppendSuffix("hDenominatorHist"),";Final State Index;Events",6,-0.5,5.5);
  AddTH1(fNumeratorHistByFinalState,AppendSuffix("hNumeratorHist"),";Final State Index;Events",6,-0.5,5.5);
  AddTH1(fNEventsHist,"hNEventsProcessedAcceptanceModHist",";;NEventsProcessed",1,-0.5,0.5);
}

//--------------------------------------------------------------------------------------------------
void AcceptanceMod::Terminate()
{
  // Create and add ratio histograms to the output list.

 if (gROOT->IsBatch()) {
   for (Int_t k=0;k<fNFinalStates;++k) {
     Double_t Acceptance;
     if (!(fFinalStatesDenom[k]>0)) Acceptance=0.0; 
     else Acceptance=(fFinalStatesNumer[k]+0.0)/(fFinalStatesDenom[k]+0.0);
     fAcceptanceHistByFinalState->Fill(k,Acceptance);
     std::cout<<"Denominator is: "<<fFinalStatesDenom[k];
     std::cout<<"Numerator is: "<<fFinalStatesNumer[k];
     std::cout<<"*******Acceptance is: "<<Acceptance<<std::endl;
   }
   std::cout<<"NEvents Processed: "<<fNEventsHist->GetEntries()<<std::endl;
 }
}
