// $Id: $

#include <TMath.h>
#include <TH1D.h>
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
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
  fTrigObjsName       (Names::gkHltObjBrn),
  fMuonName           (Names::gkMuonBrn),
  fMuons              (0),
  fElecName           (Names::gkElectronBrn),
  fElecs              (0),
  // cuts for selection
  fMinPt(15),
  fDilMinMass(0),
  fMinZMass(60),
  fMaxZMass(120),
  fIgnoreElCharge(kTRUE),
  // initialize the histograms
  fNAccCounters(0),  
  fLLMass(0),        
  fNLeptons(0),      
  fNMuons(0),         
  fNElecs(0),       
  fNAccMuons(0),      
  fNAccElecs(0),      
  fAllMuonPt(0),
  fAllMuonEta(0),
  fAllElecPt(0),
  fAllElecEta(0),
  fNVertex(0),      
  fNGPairs(0),       
  fNZPairs(0),        
  fmet(0),           
  fpt1(0),           
  fpt2(0),           
  feta1(0),           
  feta2(0),           
  fphi1(0),           
  fphi2(0),           
  fdphi(0),          
  fdeta(0),          
  fdphiMet1(0),      
  fdphiMet2(0),      
  fmt1(0),           
  fmt2(0),      		 
  frecoMass(0),	 
  ftransMass(0),	 
  ftransEll(0),	 
  ftransEnn(0),	 
  fvisMass(0),	 
  fxtau1(0),	 
  fxtau2(0)
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
  UInt_t iAccCounter = 0;
  fNAccCounters->Fill(iAccCounter++);

  // get good leptons
  const ParticleCol *leptons = GetObjThisEvt<ParticleCol>(fCleanLeptonsName);
  if (!leptons) {
    // SkipEvent();
    return;
  }
  UInt_t nLeps = leptons->GetEntries();

  // get missing Et
  const MetCol *metCol = dynamic_cast<ObjArray<Met>* >(FindObjThisEvt(fMetName));
  if (!metCol) {
    // SkipEvent();
    return;
  }
  const Met *met = metCol->At(0);
  if (!met) {
    // SkipEvent();
    return;
  }

  // access the muons
  LoadEventObject(fMuonName,fMuons);

  // access the elecs
  LoadEventObject(fElecName,fElecs);

  fNMuons -> Fill(fMuons->GetEntries());
  UInt_t accMuons = 0;
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *p = fMuons->At(i);
    fAllMuonPt ->Fill(p->Pt());
    fAllMuonEta->Fill(p->Eta());
    
    if (fabs(p->Eta()) < 2.1 && p->Pt() < fMinPt)
      {
	accMuons++;
      }
      
  }

  fNElecs -> Fill(fElecs->GetEntries());
  UInt_t accElecs = 0;
  for (UInt_t i=0; i<fElecs->GetEntries(); ++i) {
    const Electron *p = fElecs->At(i);
    fAllElecPt ->Fill(p->Pt());
    fAllElecEta->Fill(p->Eta());
    
    if (fabs(p->Eta()) < 2.1 && p->Pt() < fMinPt)
      {
	accElecs++;
      }
  }
   
  fNAccMuons -> Fill(accMuons);
  fNAccElecs -> Fill(accElecs);

  // !!! has at least one muon in acceptance
  if (accMuons < 1) {     
    // SkipEvent();
    return;
  }
  fNAccCounters->Fill(iAccCounter++);
  
  // !!! has at least one good muon
  fNAccCounters->Fill(iAccCounter++);
  
  // !!! has at least one electron 
  fNAccCounters->Fill(iAccCounter++);
  if (accElecs < 1) {     
    // SkipEvent();
    return;
  }

  // !!! has at least one good electron 
  fNAccCounters->Fill(iAccCounter++);

  // get trigger objects
  const TriggerObjectCol *trigger = GetHLTObjects(fTrigObjsName);
  if (!trigger) {             // this can only happen if HLTMod::SetAbortIfNotAccepted(kFALSE) was called
    // SkipEvent();
    return;
  }


  // !!! at least one good lepton has a trigger match
  fNAccCounters->Fill(iAccCounter++);
  Bool_t *leptonTriggerMatch = MatchTriggerCollection(fTrigObjsName, fCleanLeptonsName);
  Bool_t hasTriggerMatch = kFALSE;
  for (UInt_t i=0; i<nLeps; ++i) {
    if (leptonTriggerMatch[i]) 
      hasTriggerMatch = kTRUE;
  }
  if (!hasTriggerMatch) {       
    // SkipEvent();
    return;
  }

  // select objects for various final states
  const Particle *li = 0;
  const Particle *lj = 0;

  // make sure have found at least 2 leptons
  if (leptons->GetEntries()<2) {
    // SkipEvent();
    return;
  }

  UInt_t nZPairs    = 0;
  UInt_t nGoodPairs = 0;

  for (UInt_t i=0; i<nLeps; ++i) {
    li = leptons->At(i);

    if (li->Pt()<fMinPt)
      continue;

    for (UInt_t j=0; j<i; ++j) {
      lj = leptons->At(j);

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


      if (li->ObjType()!=lj->ObjType()) {
	DiTauSystem *diTau = new DiTauSystem(li,lj,met);
	frecoMass       ->Fill(diTau->RecoMass());      
	ftransMass      ->Fill(diTau->TransverseMass());    
	ftransEll       ->Fill(diTau->TransverseEll());    
	ftransEnn       ->Fill(diTau->TransverseEnn());    
	fvisMass        ->Fill(diTau->VisMass());    
	fxtau1          ->Fill(diTau->XTau1());        
	fxtau2          ->Fill(diTau->XTau2());    
	delete diTau;
	fLLMass->Fill(mass);
	fmet->Fill(met->Pt());          
	fpt1->Fill(li->Pt());          
	fpt2->Fill(lj->Pt());          
	feta1->Fill(li->Eta());          
	feta2->Fill(lj->Eta());          
	fphi1->Fill(li->Phi());          
	fphi2->Fill(lj->Phi());          
	fdphi->Fill(deltaPhiLeptons);        
	fdeta->Fill(deltaEtaLeptons);         
	fdphiMet1->Fill(deltaPhiMetLepton[0]);  
	fdphiMet2->Fill(deltaPhiMetLepton[1]);  
	fmt1->Fill(transMass[0]);        
	fmt2->Fill(transMass[1]);      

        ++nGoodPairs;
        continue;
      }

      if (li->Is(kMuon)) {
        if (li->Charge()!=lj->Charge()) {
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
  fNAccCounters->Fill(iAccCounter++);

  // cut on number of Z pairs
  if (nZPairs>=1) {
    // SkipEvent();
    return;
  }

  fNAccCounters->Fill(iAccCounter++);

  // cut on number of good pairs
  if (nGoodPairs<1) {
    // SkipEvent();
    return;
  }

  fNAccCounters->Fill(iAccCounter++);
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
  AddTH1(fLLMass,"hLLMass",";m_{ll} [GeV];#",150,0,300);
  AddTH1(fNLeptons,"hNLeptons",";leptons;#",10,-0.5,9.5);
  AddTH1(fNMuons,"hNMuons",";leptons;#",10,-0.5,9.5);
  AddTH1(fNElecs,"hNElecs",";leptons;#",10,-0.5,9.5);
  AddTH1(fNAccMuons,"hNAccMuons",";leptons;#",10,-0.5,9.5);
  AddTH1(fNAccElecs,"hNAccElecs",";leptons;#",10,-0.5,9.5);
  AddTH1(fAllMuonPt,"hAllMuonPt","",200,0,200);
  AddTH1(fAllMuonEta,"hAllMuonEta","",400,-4,4);
  AddTH1(fAllElecPt,"hAllElecPt","",200,0,200);
  AddTH1(fAllElecEta,"hAllElecEta","",400,-4,4);
  AddTH1(fNVertex,"hNVertex",";vertices;#",10,-0.5,9.5);
  AddTH1(fNGPairs,"hNGoodPairs",";leptons;#",10,-0.5,9.5);
  AddTH1(fNZPairs,"hNZPairs",";leptons;#",10,-0.5,9.5);

  AddTH1(fmet,"hmet","",200,0,200);
  AddTH1(fpt1,"hpt1","",200,0,200);
  AddTH1(fpt2,"hpt2","",200,0,200);
  AddTH1(feta1,"heta1","",400,-4,200);
  AddTH1(feta2,"heta2","",400,-4,200);
  AddTH1(fphi1,"hphi1","",400,-4,8);
  AddTH1(fphi2,"hphi2","",400,-4,8);
  AddTH1(fdphi,"hdphi",";#Delta #phi;#",360,0,180);
  AddTH1(fdeta,"hdeta",";#Delta #eta;#",400,0,5);   
  AddTH1(fdphiMet1,"hdphiMet1",";#Delta #phi;#",400,0,TMath::Pi());
  AddTH1(fdphiMet2,"hdphiMet2",";#Delta #phi;#",400,0,TMath::Pi());
  AddTH1(fmt1,"hmt1"," transverse Mass [GeV];#",400,0,200);
  AddTH1(fmt2,"hmt2","; transverse Mass [GeV];#",400,0,200);      		      
  AddTH1(frecoMass,"hrecoMass",";Reco Mass [GeV];#",400,0,400);     
  AddTH1(ftransMass,"htransMass",";transverse Mass [GeV];#",400,0,400);
  AddTH1(ftransEll,"htransEll",";transverse Ell [GeV];#",400,0,400); 
  AddTH1(ftransEnn,"htransEnn",";transverse Enn [GeV];#",400,0,400);
  AddTH1(fvisMass,"hvisMass",";visible mass [GeV];#",400,0,400); 
  AddTH1(fxtau1,"hxtau1","x1",400,-4,4);    
  AddTH1(fxtau2,"hxtau2","x2",400,-4,4);
  
}

//--------------------------------------------------------------------------------------------------
void ZttAnalysis::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
Bool_t *ZttAnalysis::MatchTriggerCollection(const char *fCol1Name, const char *fCol2Name)
{
  // specify some matching criteria
  Double_t fMaxEta = 2.1;  // only match up to here
  Double_t fRadius = 0.3;

  const BaseCollection *col1 = GetObjThisEvt<BaseCollection>(fCol1Name);
  const BaseCollection *col2 = GetObjThisEvt<BaseCollection>(fCol2Name);

  UInt_t ents1 = 0;
  if (col1)
    ents1 = col1->GetEntries();
  
  UInt_t ents2 = 0;
  if (col2)
    ents2 = col2->GetEntries();
  
  // set all to false
  Bool_t *found = new Bool_t[ents2];
  for (UInt_t i=0; i<ents2; ++i) 
    found[i] = kFALSE;
  
  // find matches
  for (UInt_t j=0; j<ents1; ++j) {
    const Particle *p1 = dynamic_cast<const Particle*>(col1->ObjAt(j));
    if (!p1)
      continue;

    // include later to select specific MC types or flavors
    // if (fPartType != MCParticle::kUnknown) {
    // const MCParticle *mc = dynamic_cast<const MCParticle*>(p1);
    // if (!mc || !mc->Is(fPartType))
    // continue;
    //}

    Double_t pt = p1->Pt();
    if (pt < fMinPt) 
      continue;
    if (p1->AbsEta() > fMaxEta) 
      continue;
    
    Double_t phi = p1->Phi();
    Double_t eta = p1->Eta();
    
    UInt_t foundInd  = ents2;
    Double_t foundD  = 1e12;

    for (UInt_t i=0; i<ents2; ++i) {
      if (found[i]) continue;
      const Particle *p2 = dynamic_cast<const Particle*>(col2->ObjAt(i));
      if (!p2)
        continue;
      
      if(MathUtils::DeltaR(phi, eta, p2->Phi(), p2->Eta()) < fRadius) {
        Double_t newDiff = TMath::Abs(pt - p2->Pt());
        if (newDiff < foundD) {
          foundInd = i;
          foundD = newDiff;
        } 
      }
    }
    
    if (foundInd < ents2) {
      found[foundInd] = 1;
    }
  }
  return found;
} 

//--------------------------------------------------------------------------------------------------
void ZttAnalysis::MatchCollections(Bool_t *found, const char *fCol1Name, const char *fCol2Name)
{
  // specify some matching criteria
  Double_t fMaxEta = 2.1;
  Double_t fRadius = 0.3;

  const BaseCollection *col1 = GetObjThisEvt<BaseCollection>(fCol1Name);
  const BaseCollection *col2 = GetObjThisEvt<BaseCollection>(fCol2Name);

  UInt_t ents1 = 0;
  if (col1)
    ents1 = col1->GetEntries();
  
  UInt_t ents2 = 0;
  if (col2)
    ents2 = col2->GetEntries();
  
  // set all to false
  found = new Bool_t[ents2];
  for (UInt_t i=0; i<ents2; ++i) 
    found[i] = kFALSE;
  
  // find matches
  for (UInt_t j=0; j<ents1; ++j) {
    const Particle *p1 = dynamic_cast<const Particle*>(col1->ObjAt(j));
    if (!p1)
      continue;

    // include later to select specific MC types or flavors
    // if (fPartType != MCParticle::kUnknown) {
    // const MCParticle *mc = dynamic_cast<const MCParticle*>(p1);
    // if (!mc || !mc->Is(fPartType))
    // continue;
    //}

    Double_t pt = p1->Pt();
    if (pt < fMinPt) 
      continue;
    if (p1->AbsEta() > fMaxEta) 
      continue;
    
    Double_t phi = p1->Phi();
    Double_t eta = p1->Eta();
    
    UInt_t foundInd  = ents2;
    Double_t foundD  = 1e12;

    for (UInt_t i=0; i<ents2; ++i) {
      if (found[i]) continue;
      const Particle *p2 = dynamic_cast<const Particle*>(col2->ObjAt(i));
      if (!p2)
        continue;
      
      if(MathUtils::DeltaR(phi, eta, p2->Phi(), p2->Eta()) < fRadius) {
        Double_t newDiff = TMath::Abs(pt - p2->Pt());
        if (newDiff < foundD) {
          foundInd = i;
          foundD = newDiff;
        } 
      }
    }
    
    if (foundInd < ents2) {
      found[foundInd] = 1;
    }
  }
  return;
}

