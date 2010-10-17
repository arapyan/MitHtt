// $Id: $

#include <TMath.h>
#include <TH1D.h>
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitHtt/Mods/interface/ZttAnalysis.h"

using namespace mithep;

ClassImp(mithep::ZttAnalysis)

//--------------------------------------------------------------------------------------------------
ZttAnalysis::ZttAnalysis(const char *name, const char *title) :
  BaseMod             (name,title),
  // initialize bambu objects
  fCleanLeptonsName(ModNames::gkMergedLeptonsName),
  fMetName("PubPFMet"),
  // cuts for selection
  fMinPt(15),
  fDilMinMass(0),
  fMinZMass(60),
  fMaxZMass(120),
  fIgnoreElCharge(kTRUE),
  // initialize the histograms
  fNAccCounters(0),
  fAllDiLepMass(0),
  fDiElMass(0),
  fDiMuMass(0),
  fElMuMass(0),
  fAllDiLepMassAcc(0),
  fDiElMassAcc(0),
  fDiMuMassAcc(0),
  fElMuMassAcc(0),
  fNLeptons(0),
  fNGPairs(0),
  fNZPairs(0),
  fAllDiLepMet(0),        
  fAllDiLepPt1(0),        
  fAllDiLepPt2(0),        
  fAllDiLepDPhi(0),
  fAllDiLepDEta(0),    
  fAllDiLepDPhiMetOne(0), 
  fAllDiLepDPhiMetTwo(0), 
  fAllDiLepMtOne(0),  
  fAllDiLepMtTwo(0),      		      
  fDiTauRecoMass(0),      
  fDiTauTransverseMass(0),
  fDiTauTransverseEll(0), 
  fDiTauTransverseEnn(0), 
  fDiTauVisMass(0),   
  fDiTauXTau1(0),      
  fDiTauXTau2(0), 
  fEmuRecoMass(0),      
  fEmuTransverseMass(0),
  fEmuTransverseEll(0), 
  fEmuTransverseEnn(0), 
  fEmuVisMass(0),   
  fEmuXTau1(0),      
  fEmuXTau2(0), 
  fMumuRecoMass(0),      
  fMumuTransverseMass(0),
  fMumuTransverseEll(0), 
  fMumuTransverseEnn(0), 
  fMumuVisMass(0),   
  fMumuXTau1(0),      
  fMumuXTau2(0), 
  fEeRecoMass(0),      
  fEeTransverseMass(0),
  fEeTransverseEll(0), 
  fEeTransverseEnn(0), 
  fEeVisMass(0),   
  fEeXTau1(0),      
  fEeXTau2(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ZttAnalysis::Begin()
{
  // Run startup code on the client machine. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void ZttAnalysis::Process()
{
  // count the events we have processed
  IncNEventsProcessed();

  // get good leptons
  const ParticleCol *leptons = GetObjThisEvt<ParticleCol>(fCleanLeptonsName);

  if (!leptons) {
    SkipEvent();
    return;
  }

  // get missing Et
  const MetCol *metCol = dynamic_cast<ObjArray<Met>* >(FindObjThisEvt(fMetName));
 
  if (!metCol)
    {
      SkipEvent();
      return;
    }

  const Met *met = metCol->At(0);

  if (!met)
    {
      SkipEvent();
      return;
    }


  fNAccCounters->Fill(1);

  // make sure have found at least 2 leptons
  if (leptons->GetEntries()<2) {
    SkipEvent();
    return;
  }

  // create di-tau system from 1st and 2nd lepton
  DiTauSystem *diTau = new DiTauSystem(leptons->At(0),leptons->At(1),met);
 
  fDiTauRecoMass       ->Fill(diTau->RecoMass());      
  fDiTauTransverseMass ->Fill(diTau->TransverseMass());    
  fDiTauTransverseEll  ->Fill(diTau->TransverseEll());    
  fDiTauTransverseEnn  ->Fill(diTau->TransverseEnn());    
  fDiTauVisMass        ->Fill(diTau->VisMass());    
  fDiTauXTau1          ->Fill(diTau->XTau1());        
  fDiTauXTau2          ->Fill(diTau->XTau2());    

  delete diTau;

  UInt_t nLeps = leptons->GetEntries();
  UInt_t nZPairs    = 0;
  UInt_t nGoodPairs = 0;

  for (UInt_t i=0; i<nLeps; ++i) {
    const Particle *li = leptons->At(i);

    if (li->Pt()<fMinPt)
      continue;

    for (UInt_t j=0; j<i; ++j) {
      const Particle *lj = leptons->At(j);

      if (lj->Pt()<fMinPt)
        continue;

      // found two leptons with min Pt      
      CompositeParticle dil;
      dil.AddDaughter(li);
      dil.AddDaughter(lj);
      Double_t mass = dil.Mass();
      if (mass<fDilMinMass)
        continue;

      // some kinematics
      double deltaPhiLeptons      = MathUtils::DeltaPhi(li->Phi(),lj->Phi())* 180.0 / TMath::Pi();
      
      double deltaEtaLeptons      = TMath::Abs(li->Eta() - lj->Eta());

      double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(met->Phi(), li->Phi()),
				     MathUtils::DeltaPhi(met->Phi(), lj->Phi())};

      double transMass[2]         = {TMath::Sqrt(2.0*li->Pt()*met->Pt()*
						 (1.0 - cos(deltaPhiMetLepton[0]))),
				     TMath::Sqrt(2.0*lj->Pt()*met->Pt()*
						 (1.0 - cos(deltaPhiMetLepton[1])))};

      fAllDiLepMass->Fill(mass);
      fAllDiLepMet->Fill(met->Pt());          
      fAllDiLepPt1->Fill(li->Pt());          
      fAllDiLepPt2->Fill(lj->Pt());          
      fAllDiLepDPhi->Fill(deltaPhiLeptons);        
      fAllDiLepDEta->Fill(deltaEtaLeptons);         
      fAllDiLepDPhiMetOne->Fill(deltaPhiMetLepton[0]);  
      fAllDiLepDPhiMetTwo->Fill(deltaPhiMetLepton[1]);  
      fAllDiLepMtOne->Fill(transMass[0]);        
      fAllDiLepMtTwo->Fill(transMass[1]);      

      if (li->ObjType()!=lj->ObjType()) {
        fElMuMass->Fill(mass);
	DiTauSystem *diTau = new DiTauSystem(li,lj,met);
	fEmuRecoMass       ->Fill(diTau->RecoMass());      
	fEmuTransverseMass ->Fill(diTau->TransverseMass());    
	fEmuTransverseEll  ->Fill(diTau->TransverseEll());    
	fEmuTransverseEnn  ->Fill(diTau->TransverseEnn());    
	fEmuVisMass        ->Fill(diTau->VisMass());    
	fEmuXTau1          ->Fill(diTau->XTau1());        
	fEmuXTau2          ->Fill(diTau->XTau2());    
	delete diTau;

        ++nGoodPairs;
        continue;
      }

      if (li->Is(kMuon)) {
        if (li->Charge()!=lj->Charge()) {
          fDiMuMass->Fill(mass);
	  DiTauSystem *diTau = new DiTauSystem(li,lj,met);
	  fMumuRecoMass       ->Fill(diTau->RecoMass());      
	  fMumuTransverseMass ->Fill(diTau->TransverseMass());    
	  fMumuTransverseEll  ->Fill(diTau->TransverseEll());    
	  fMumuTransverseEnn  ->Fill(diTau->TransverseEnn());    
	  fMumuVisMass        ->Fill(diTau->VisMass());    
	  fMumuXTau1          ->Fill(diTau->XTau1());        
	  fMumuXTau2          ->Fill(diTau->XTau2());    
	  delete diTau;
          if ((mass>fMinZMass) && (mass<fMaxZMass)) {
            ++nZPairs;
            continue;
          }
        }
        ++nGoodPairs;
        continue;
      }

      if (li->Is(kElectron)) {
        if (fIgnoreElCharge || (li->Charge()!=lj->Charge())) {
          fDiElMass->Fill(mass);
	  DiTauSystem *diTau = new DiTauSystem(li,lj,met);
	  fEeRecoMass       ->Fill(diTau->RecoMass());      
	  fEeTransverseMass ->Fill(diTau->TransverseMass());    
	  fEeTransverseEll  ->Fill(diTau->TransverseEll());    
	  fEeTransverseEnn  ->Fill(diTau->TransverseEnn());    
	  fEeVisMass        ->Fill(diTau->VisMass());    
	  fEeXTau1          ->Fill(diTau->XTau1());        
	  fEeXTau2          ->Fill(diTau->XTau2());    
	  delete diTau;
          if ((mass>fMinZMass) && (mass<fMaxZMass)) {
            ++nZPairs;
            continue;
          }
        }
        ++nGoodPairs;
        continue;
      }
    }
  }

  fNLeptons->Fill(nLeps);
  fNGPairs->Fill(nGoodPairs);
  fNZPairs->Fill(nZPairs);
  fNAccCounters->Fill(2);

  // cut on number of Z pairs
  if (nZPairs>=1) {
    SkipEvent();
    return;
  }

  fNAccCounters->Fill(3);

  // cut on number of good pairs
  if (nGoodPairs<1) {
    SkipEvent();
    return;
  }

  fNAccCounters->Fill(4);
  for (UInt_t i=0; i<nLeps; ++i) {
    const Particle *li = leptons->At(i);

    if (li->Pt()<fMinPt)
      continue;

    for (UInt_t j=0; j<i; ++j) {
      const Particle *lj = leptons->At(j);

      if (lj->Pt()<fMinPt)
        continue;

      CompositeParticle dil;
      dil.AddDaughter(li);
      dil.AddDaughter(lj);
      Double_t mass = dil.Mass();
      if (mass<fDilMinMass)
        continue;

      fAllDiLepMassAcc->Fill(mass);

      if (li->ObjType()!=lj->ObjType()) {
        fElMuMassAcc->Fill(mass);
        continue;
      }

      if (li->Is(kMuon)) {
        if (li->Charge()!=lj->Charge()) {
          fDiMuMassAcc->Fill(mass);
        continue;
        }
      }

      if (li->Is(kElectron)) {
        if (fIgnoreElCharge || (li->Charge()!=lj->Charge())) {
          fDiElMassAcc->Fill(mass);
        }
        continue;
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------------------
void ZttAnalysis::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.

}

//--------------------------------------------------------------------------------------------------
void ZttAnalysis::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches or objects created by earlier modules.

  // initialize some variables


  // request the following Event Objects (often branches, but could be locally created ones)

  // for MC only to adjust potential overlaps from generation

  // book all histograms

  AddTH1(fNAccCounters,"hNAccCounters",";cut;#",6,-0.5,5.5);
  if (1) {
    TAxis *xa = fNAccCounters->GetXaxis();
    for(Int_t i=1;i<=fNAccCounters->GetNbinsX();++i)
      xa->SetBinLabel(i,"unused");
    xa->SetBinLabel(1,"Enter");
    xa->SetBinLabel(2,"Objs");
    xa->SetBinLabel(3,"2Lep");
    xa->SetBinLabel(4,"ZPair");
    xa->SetBinLabel(5,"GPair");
    xa->SetRangeUser(0,4);
  }
  AddTH1(fAllDiLepMass,"hAllDiLepMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fDiElMass,"hDiElMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fDiMuMass,"hDiMuMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fElMuMass,"hElMuMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fAllDiLepMassAcc,"hAllDiLepMassAcc",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fDiElMassAcc,"hDiElMassAcc",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fDiMuMassAcc,"hDiMuMassAcc",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fElMuMassAcc,"hElMuMassAcc",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fNLeptons,"hNLeptons",";leptons;#",10,-0.5,9.5);
  AddTH1(fNGPairs,"hNGoodPairs",";leptons;#",10,-0.5,9.5);
  AddTH1(fNZPairs,"hNZPairs",";leptons;#",10,-0.5,9.5);
  AddTH1(fAllDiLepMet,"hAllDiLepMet","",200,0,200);
  AddTH1(fAllDiLepPt1,"hAllDiLepPt1","",200,0,200);
  AddTH1(fAllDiLepPt2,"hAllDiLepPt2","",200,0,200);
  AddTH1(fAllDiLepDPhi,"hAllDiLepDPhi",";#Delta #phi;#",360,0,180);
  AddTH1(fAllDiLepDEta,"hAllDiLepDEta",";#Delta #eta;#",400,0,5);   
  AddTH1(fAllDiLepDPhiMetOne,"hAllDiLepDPhiMetOne",";#Delta #phi;#",400,0,TMath::Pi());
  AddTH1(fAllDiLepDPhiMetTwo,"hAllDiLepDPhiMetTwo",";#Delta #phi;#",400,0,TMath::Pi());
  AddTH1(fAllDiLepMtOne,"hAllDiLepMtOne"," transverse Mass [GeV];#",400,0,200);
  AddTH1(fAllDiLepMtTwo,"hAllDiLepMtTwo","; transverse Mass [GeV];#",400,0,200);      		      
  AddTH1(fDiTauRecoMass,"hDiTauRecoMass",";Reco Mass [GeV];#",400,0,400);     
  AddTH1(fDiTauTransverseMass,"hDiTauTransverseMass",";transverse Mass [GeV];#",400,0,400);
  AddTH1(fDiTauTransverseEll,"hDiTauTransverseEll",";transverse Ell [GeV];#",400,0,400); 
  AddTH1(fDiTauTransverseEnn,"hDiTauTransverseEnn",";transverse Enn [GeV];#",400,0,400);
  AddTH1(fDiTauVisMass,"hDiTauVisMass",";visible mass [GeV];#",400,0,400); 
  AddTH1(fDiTauXTau1,"hDiTauXTau1","x1",400,-4,4);    
  AddTH1(fDiTauXTau2,"hDiTauXTau2","x2",400,-4,4);
  AddTH1(fEmuRecoMass,"hEmuRecoMass",";Reco Mass [GeV];#",400,0,400);     
  AddTH1(fEmuTransverseMass,"hEmuTransverseMass",";transverse Mass [GeV];#",400,0,400);
  AddTH1(fEmuTransverseEll,"hEmuTransverseEll",";transverse Ell [GeV];#",400,0,400); 
  AddTH1(fEmuTransverseEnn,"hEmuTransverseEnn",";transverse Enn [GeV];#",400,0,400);
  AddTH1(fEmuVisMass,"hEmuVisMass",";visible mass [GeV];#",400,0,400); 
  AddTH1(fEmuXTau1,"hEmuXTau1","x1",400,-4,4);    
  AddTH1(fEmuXTau2,"hEmuXTau2","x2",400,-4,4);
  AddTH1(fMumuRecoMass,"hMumuRecoMass",";Reco Mass [GeV];#",400,0,400);     
  AddTH1(fMumuTransverseMass,"hMumuTransverseMass",";transverse Mass [GeV];#",400,0,400);
  AddTH1(fMumuTransverseEll,"hMumuTransverseEll",";transverse Ell [GeV];#",400,0,400); 
  AddTH1(fMumuTransverseEnn,"hMumuTransverseEnn",";transverse Enn [GeV];#",400,0,400);
  AddTH1(fMumuVisMass,"hMumuVisMass",";visible mass [GeV];#",400,0,400); 
  AddTH1(fMumuXTau1,"hMumuXTau1","x1",400,-4,4);    
  AddTH1(fMumuXTau2,"hMumuXTau2","x2",400,-4,4);
  AddTH1(fEeRecoMass,"hEeRecoMass",";Reco Mass [GeV];#",400,0,400);     
  AddTH1(fEeTransverseMass,"hEeTransverseMass",";transverse Mass [GeV];#",400,0,400);
  AddTH1(fEeTransverseEll,"hEeTransverseEll",";transverse Ell [GeV];#",400,0,400); 
  AddTH1(fEeTransverseEnn,"hEeTransverseEnn",";transverse Enn [GeV];#",400,0,400);
  AddTH1(fEeVisMass,"hEeVisMass",";visible mass [GeV];#",400,0,400); 
  AddTH1(fEeXTau1,"hEeXTau1","x1",400,-4,4);    
  AddTH1(fEeXTau2,"hEeXTau2","x2",400,-4,4);
}

//--------------------------------------------------------------------------------------------------
void ZttAnalysis::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
