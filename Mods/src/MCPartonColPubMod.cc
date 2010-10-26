
// $Id:$

#include "MitHtt/Mods/interface/MCPartonColPubMod.h"
#include <TH1D.h>
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"

using namespace mithep;

ClassImp(mithep::MCPartonColPubMod)

//--------------------------------------------------------------------------------------------------
MCPartonColPubMod::MCPartonColPubMod(const char *name, const char *title) : 
  BaseMod(name,title),
  fMCParticleColName("MCParticles"),
  fOutputColName("OutputColName"),
  fPubPerEvent(kTRUE),
  fDaughterPdgId(0),
  fMotherPdgId(0),
  fOutputCol(0),
  fMCParticleCol(0),
  fIntermediates(new std::vector<Int_t>)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void MCPartonColPubMod::Begin()
{
  // Run startup code on the client machine. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void MCPartonColPubMod::Process()
{
  // Process entries of the tree. For this module, we just load the branches and
  // fill the histograms.
 
  LoadEventObject(GetMCParticleColName(),fMCParticleCol);

  UInt_t entries=fMCParticleCol->GetEntries();
  
  
  //std::cout <<"Through with initial histograms" <<std::endl;
  //std::cout<<"Creating output collection to store selected MCParticles"<<std::endl;
  
  fOutputCol=new mithep::ObjArray<MCParticle>(entries,GetOutputName());
  
  std::cout << "Created the new object collection" <<std::endl;

  for (UInt_t i=0;i<entries;++i) {
    const MCParticle* partA=fMCParticleCol->At(i);
    //fOutputCol->Add(fMCParticleCol->At(i));
    //Particle needs to be simulated and have the appropriate PdgId to merit further consideration
    Bool_t AGen= partA->IsGenerated();
    if (!AGen) continue;
    if (partA->AbsPdgId()!= fMotherPdgId) continue;
    std::cout<<"We have identified a mother particle"<<std::endl;
    const MCParticle *dp;
    
    for (UInt_t daug=0;daug<partA->NDaughters();++daug) {
      
      dp=partA->Daughter(daug);

      if (dp->AbsPdgId()==fMotherPdgId) continue;

      //std::cout<<"DaughterPdgId is: "<<dp->PdgId()<<std::endl;
      for (UInt_t gaugeIn=0;gaugeIn<fIntermediates->size();++gaugeIn) {
	if (dp->AbsPdgId()==fIntermediates->at(gaugeIn)) {
	  const MCParticle* dsec=dp->FindDaughter(fDaughterPdgId);
	  if (!dsec) continue;
	  dp=dsec;
	  break;
	}
      }
      std::cout << "DaughterPdgId is: "<<dp->PdgId()<<std::endl;

      if (dp->AbsPdgId()!=fDaughterPdgId) continue;
      
      while (dp->HasDaughter(fDaughterPdgId)) {
	const MCParticle* dn=dp->FindDaughter(fDaughterPdgId);
	if (!dn) break;
	if (dn->IsGenerated()) dp=dn;
	else break;
      }
      
      fOutputCol->Add(dp);
    }
  }

	/*for (UInt_t DPdgIndex=0;DPdgIndex<fDaughterPdgId->size();++DPdgIndex) {
      std::cout<<"In viable daughter particles loop"<<std::endl;
      const MCParticle *dp=partA->FindDaughter(fDaughterPdgId->at(DPdgIndex));
      const MCParticle* permanentFirstDp=dp;
      if (!dp) continue;
      
      //Search for final generated particle
      while (dp->HasDaughter(fDaughterPdgId->at(DPdgIndex))) {
	std::cout<<"In FSR while loop" <<std::endl;
        const MCParticle* daughterCandidate=dp->FindDaughter(fDaughterPdgId->at(DPdgIndex));
	if (!daughterCandidate) break;
	if (daughterCandidate->IsGenerated())
	  dp=daughterCandidate;
        else break;
      }

      fOutputCol->Add(dp);
      
      const MCParticle *dp2=partA->FindDaughter(fDaughterPdgId->at(DPdgIndex),kFALSE,permanentFirstDp);

      std::cout << "Checking for existence of second particle" <<std::endl;
      if (!dp2) continue;
      std::cout<<"Found a second particle!**************"<<std::endl;
      const MCParticle* permanentSecondDp=dp2;
       //Search for final generated particle
      while (dp2->HasDaughter(fDaughterPdgId->at(DPdgIndex))) {
	std::cout<<"In FSR while loop" <<std::endl;
        const MCParticle* daughterCandidate=dp2->FindDaughter(fDaughterPdgId->at(DPdgIndex)); //We assume that there's only one radiation product of the same kind, probably a reasonable assumption
	if (daughterCandidate->IsGenerated())
	  dp2=daughterCandidate;
        else break;
      }

      fOutputCol->Add(dp2);
	*/

  std::cout<<"Collection for AbsPdgId of : " <<fDaughterPdgId<<" has "<<fOutputCol->Entries()<<" entries"<<std::endl;
  if (fPubPerEvent)
    AddObjThisEvt(fOutputCol);


}

//--------------------------------------------------------------------------------------------------
void MCPartonColPubMod::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches. For this module, we request a branch of the MitTree.

  ReqEventObject(GetMCParticleColName(),fMCParticleCol,kTRUE);
 
}
//--------------------------------------------------------------------------------------------------
void MCPartonColPubMod::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
  std::cout<<"In SlaveTerminate..."<<std::endl;
 
  
}

//--------------------------------------------------------------------------------------------------
void MCPartonColPubMod::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
  std::cout<<"In terminate..."<<std::endl;
}
