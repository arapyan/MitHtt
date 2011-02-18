// $Id: $

#include <TMath.h>
#include <TH1D.h>
#include <TNtuple.h>
#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"
#include "MitHtt/Mods/interface/EMUAnalysis.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/JetCol.h"

using namespace mithep;

ClassImp(mithep::EMUAnalysis)

//--------------------------------------------------------------------------------------------------
EMUAnalysis::EMUAnalysis(const char *name, const char *title) :
  BaseMod             (name,title),
  // initialize bambu objects
  // fEventHeader(Names::gkEvtTreeName),
  fCleanLeptonsName(ModNames::gkMergedLeptonsName),
  fTrigObjsName(Names::gkHltObjBrn),
  fMuonsName(ModNames::gkCleanMuonsName),        
  fElectronsName(ModNames::gkCleanElectronsName),        
  fJetsName(ModNames::gkPubJetsName),
  fCaloJetsName("AKt5Jets"),
  fMetName("PubPFMet"),
  // initialize the histograms
  cutPtLeadingMuon(15.),
  cutPtLeadingElec(15.),
  cutPtSecondMuon(15.),
  cutPtSecondElec(15.),
  cutPtTriggerMuon(15.),
  cutPtTriggerElec(15.),
  cutTriggerMuon(1),
  cutTriggerElec(0),
  fNAccCounters(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void EMUAnalysis::Begin()
{
  // Run startup code on the client machine. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void EMUAnalysis::Process()
{
  // count the events we have processed

  IncNEventsProcessed();
  UInt_t iAccCounter = 0; 
  fNAccCounters->Fill(iAccCounter++);

  const MuonCol     *muons     = GetObjThisEvt<MuonCol>(fMuonsName);
  const ElectronCol *electrons = GetObjThisEvt<ElectronCol>(fElectronsName);
  const JetCol      *jets      = GetObjThisEvt<JetCol>(fJetsName);//Corrected, ID'd ("good"), Cleaned PF jets 
  const MetCol      *metCol    = dynamic_cast<ObjArray<Met>* >(FindObjThisEvt(fMetName));
  
  LoadBranch(fCaloJetsName);//ungood, unclean calojets (?)

  if (!muons) { SkipEvent(); return; }

  if (!electrons) {
    SkipEvent();
    return;
  }

  if (!jets) {
    SkipEvent();
    return;
  }

  if (!metCol) {
    SkipEvent();
    return;
  }

  const Met *met = metCol->At(0);
  if (!met) {
    SkipEvent();
    return;
  }

  fNAccCounters->Fill(iAccCounter++);
  // get trigger
  const TriggerObjectCol *tos = GetHLTObjects(fTrigObjsName);
  if (! tos)             // this can only happen if HLTMod::SetAbortIfNotAccepted(kFALSE) was called
    return;
  UInt_t nEnts = tos->GetEntries();
  fNAccCounters->Fill(iAccCounter++);


  // get started with the analysis
  UInt_t nMuons   = muons->GetEntries();
  UInt_t nElecs   = electrons->GetEntries();
  UInt_t nJets    = jets->GetEntries();
  UInt_t nJetsPt15= 0;
  UInt_t nJetsPt20= 0;
  UInt_t nJetsPt25= 0;
  UInt_t nBJetsPt20SV2= 0;
  UInt_t nBJetsPt20SV05= 0;
  Double_t btag = 0;//highest b-tag disc of any of the calojets 

  // dimuon and dielectron mass if there is two ee or two mu
  if (nMuons > 1) {
    if (muons->At(0)->Pt() > cutPtLeadingMuon && muons->At(1)->Pt() > cutPtSecondMuon) {
      CompositeParticle dil;
      dil.AddDaughter(muons->At(0));
      dil.AddDaughter(muons->At(1));
      fmmmass -> Fill(dil.Mass());
    }
  }

  if (nElecs > 1) {
    if (electrons->At(0)->Pt() > cutPtLeadingElec && electrons->At(1)->Pt() > cutPtSecondElec) {
      CompositeParticle dil;
      dil.AddDaughter(electrons->At(0));
      dil.AddDaughter(electrons->At(1));
      feemass -> Fill(dil.Mass());
    }
  }

  for (UInt_t i=0; i<nJets; ++i) {//jets: good, clean calojets. fCaloJet: ungood, unclean calojets?
    if (jets->At(i)->Pt() > 15.) nJetsPt15++;
    if (jets->At(i)->Pt() > 20.) 
      {
	nJetsPt20++;
	// match pf and calo jets and extract b-tagging info
	int nCloseStdJet = -1;
	double deltaRMin = 999.;
	for(UInt_t njet=0; njet<fCaloJet->GetEntries(); njet++){//Loop through the calojets to find
	  const CaloJet *caloJet = fCaloJet->At(njet);          //the closest one
	  Double_t deltaR = MathUtils::DeltaR(jets->At(i)->Mom(),caloJet->Mom());
	  if(deltaR < deltaRMin) {
	    nCloseStdJet = njet;
	    deltaRMin = deltaR;
	  }
	}
	
	if(nCloseStdJet >= 0 && deltaRMin < 0.5)//If we found any jets, and if deltaR of closest
	  {                                     //is below 0.5
	    Double_t csvbtag = fCaloJet->At(nCloseStdJet)->TrackCountingHighEffBJetTagsDisc();
	    if (csvbtag > btag ) btag = csvbtag; //Looking for calojet with highest b-tag disc.
	    if (csvbtag > 2.)  nBJetsPt20SV2++;  //just counting the number of b-tagged calojets
	    if (csvbtag > 0.5) nBJetsPt20SV05++;
	  }
      }

    if (jets->At(i)->Pt() > 25.) nJetsPt25++;
  }

  // fill the ntuple for events with at least one electron and muon and cuts as defined in the ID modules.
  if (nMuons > 0 && nElecs > 0)
    {
      const Muon *muon = muons->At(0);
      const Electron *elec = electrons->At(0);
      DiTauSystem *diTau = new DiTauSystem(muon,elec,met);
          
      // !!! leading muon fired trigger 15 GeV 
      UInt_t matchedTrigMuons = 0;
      UInt_t matchedTrigElecs = 0;
      //Just counting how many are matched
      if (1) {
	for (UInt_t j=0; j<nEnts; ++j) {//nEnts: number of trigger objects
	  const TriggerObject *to = tos->At(j);
	  if (to->Pt() > cutPtTriggerMuon) {
	    if (fabs(muon->Eta()-to->Eta()) < 0.20 && fabs(muon->Phi()-to->Phi()) < 0.10)   
	      matchedTrigMuons++;
	    if (to->Pt() > cutPtTriggerElec) {
	      if (fabs(elec->Eta()-to->Eta()) < 0.20 && fabs(elec->Phi()-to->Phi()) < 0.10)   
		matchedTrigElecs++;
	    } 
	  }
	}
      }	

      // write in ntuple
      int lEvt   = GetEventHeader()->EvtNum();
      int lLumi  = GetEventHeader()->LumiSec();
      int lRun   = GetEventHeader()->RunNum();
      if (1)
	{
	  Float_t vals[17] = {lEvt, lLumi, lRun, muon->Pt(), elec->Pt(), 
			      muon->Eta(),elec->Eta(),  muon->Phi(), elec->Phi(),
			      muon->D0PV(), elec->D0PV(), muon->Charge(), elec->Charge(), 
			      diTau->VisMass(), 
			      met->Pt(), matchedTrigMuons, matchedTrigMuons};
	  fNt->Fill(vals);
	}
      delete diTau; 
    }
  // Now start cutting: 

  // !!! has at least one good muon
  if (nMuons < 1) return;
  if (muons->At(0)->Pt() < cutPtLeadingMuon) return;
  fNAccCounters->Fill(iAccCounter++);
  
  // !!! has at least one good electron 
  if (nElecs < 1) return;

  double e_scale = 1.;
  if (electrons->At(0)->IsEE())
    e_scale = 1.;

  if (e_scale*electrons->At(0)->Pt() < cutPtLeadingElec) return;
  fNAccCounters->Fill(iAccCounter++);

  // select the first two leptons
  const Muon *muon = muons->At(0);
  const Electron *elec = electrons->At(0);

  // !!! leading muon fired trigger 15 GeV 
  //Just count the number of matched elecs and muons:
  UInt_t matchedTrigMuons = 0; //(This is the second time we've done this)
  UInt_t matchedTrigElecs = 0;
  if (1) {
    for (UInt_t j=0; j<nEnts; ++j) {
      const TriggerObject *to = tos->At(j);
      if (to->Pt() > cutPtTriggerMuon) {
	if (fabs(muon->Eta()-to->Eta()) < 0.20 && fabs(muon->Phi()-to->Phi()) < 0.10)   
	  matchedTrigMuons++;
	
	if (to->Pt() > cutPtTriggerElec) {
	  if (fabs(elec->Eta()-to->Eta()) < 0.20 && fabs(elec->Phi()-to->Phi()) < 0.10)   
	    matchedTrigElecs++;
	} 
      }
    }
  }

  fNAccCounters->Fill(iAccCounter++);

  // !!! has at least two good leptons of opposite charge
  if (muon->Charge() * elec->Charge() > 0 ) {
    return;
  }
  fNAccCounters->Fill(iAccCounter++);

  fptmnotrig->Fill(muon->Pt());          
  if (matchedTrigMuons < cutTriggerMuon) { //default: 1 for mu
    return;
  }
  fptmtrig->Fill(muon->Pt());          

  if (matchedTrigElecs < cutTriggerElec) { //default: 0 for elec
    return;
  }

  fNAccCounters->Fill(iAccCounter++);

  // muon plots
  fptm->Fill(muon->Pt());          
  fetam->Fill(muon->Eta());          
  fphim->Fill(muon->Phi());          
  fdcam->Fill(muon->D0PV());          
  fchargem->Fill(muon->Charge()); 
  // electron plots
  fpte->Fill(elec->Pt());          
  fetae->Fill(elec->Eta());          
  fphie->Fill(elec->Phi());          
  fdcae->Fill(elec->D0PV());          
  fchargee->Fill(elec->Charge()); 
 
  // !!! di lepton pair has minimal mass
  fNAccCounters->Fill(iAccCounter++);
    
  // di-tau kinematics
  DiTauSystem *diTau = new DiTauSystem(muon,elec,met);
  frecoMass       ->Fill(diTau->RecoMass());      
  ftransMass      ->Fill(diTau->TransverseMass());    
  ftransEll       ->Fill(diTau->TransverseEll());    
  ftransEnn       ->Fill(diTau->TransverseEnn());    
//   fproj           ->Fill(diTau->Projected());  
//   fprojVis        ->Fill(diTau->ProjectedVis()); 
//   fprojMet        ->Fill(diTau->ProjectedMet()); 
//   fprojPhi        ->Fill(diTau->ProjectedPhi()); 
//   fhT             ->Fill(diTau->Ht());    

  fvisMass        ->Fill(e_scale*diTau->VisMass());    
  fscaledVisMass  ->Fill(e_scale*1.8*diTau->VisMass());    
  fxtaum          ->Fill(diTau->XTau1());        
  fxtaue          ->Fill(diTau->XTau2());  

  double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(met->Phi(), muon->Phi()),
				 MathUtils::DeltaPhi(met->Phi(), elec->Phi())};
  
  double mTW[2] = {TMath::Sqrt(2.0*muon->Pt()*met->Pt()*
			       (1.0 - cos(deltaPhiMetLepton[0]))),
		   TMath::Sqrt(2.0*elec->Pt()*met->Pt()*
			       (1.0 - cos(deltaPhiMetLepton[1])))};

  ftransMnmu -> Fill(mTW[0]);
  ftransMnel -> Fill(mTW[1]);
  fdphinmu   -> Fill(deltaPhiMetLepton[0]* 180.0 / TMath::Pi());
  fdphinel   -> Fill(deltaPhiMetLepton[1]* 180.0 / TMath::Pi());
  
  double deltaPhiLeptons      = MathUtils::DeltaPhi(muon->Phi(),elec->Phi())* 180.0 / TMath::Pi();  
  double deltaEtaLeptons      = TMath::Abs(muon->Eta() - elec->Eta());

  fdphi           ->Fill(deltaPhiLeptons);  
  fdeta           ->Fill(deltaEtaLeptons);  

  fcharge         ->Fill(muon->Charge() * elec->Charge());
  ftype           ->Fill(0);

  fnleptons       ->Fill(nMuons+nElecs);
  fnmuons         ->Fill(nMuons);
  fnelecs         ->Fill(nElecs);
  fnjets          ->Fill(nJets);
  fnjetspt15      ->Fill(nJetsPt15);
  fnjetspt20      ->Fill(nJetsPt20);
  fnbjetspt20sv2  ->Fill(nBJetsPt20SV2);
  fnbjetspt20sv05 ->Fill(nBJetsPt20SV05);
  fnjetspt25      ->Fill(nJetsPt25);
  fmet            ->Fill(met->Pt());
  fbtag           ->Fill(btag);
  //cout << "Btag = " << btag << endl;

  // Higgs selection cuts
  if (mTW[0] < 60. && mTW[1] < 60.)
    fvisMassCut1 ->Fill(diTau->VisMass());
  //
  if (mTW[0] < 50. && mTW[1] < 50.)
    fvisMassCut2 ->Fill(diTau->VisMass());
//   if (diTau->ProjectedVis() < 20.)
//     fvisMassCut3 ->Fill(diTau->VisMass());  

  // require just one b-tag
  if (btag > 2.1) //btag is highest discr. value among the jets 
    fvisMassCut4 ->Fill(diTau->VisMass());

  // require transverse mass and btag or w/o btag
  if (mTW[0] < 50. && mTW[1] < 50. && btag > 2.1)
    fvisMassCut5 ->Fill(diTau->VisMass());
  if (mTW[0] < 50. && mTW[1] < 50. && btag <= 2.1)
    fvisMassCut6 ->Fill(diTau->VisMass());

  // ttbar control region
  if (mTW[0] > 50. && mTW[1] > 50.)
    fvisMassCut7 ->Fill(diTau->VisMass());
  if (mTW[0] > 50. && mTW[1] > 50.&& btag > 2.1)
    fvisMassCut8 ->Fill(diTau->VisMass());

  // !!! transverse mass cut
  fNAccCounters->Fill(iAccCounter++);

  // !!! Pzeta - 1.5 Pzeta cut?
  fNAccCounters->Fill(iAccCounter++);

  // !!! Z removal
  fNAccCounters->Fill(iAccCounter++);
  
  delete diTau;

  return;
  }

//--------------------------------------------------------------------------------------------------
void EMUAnalysis::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.

}

//--------------------------------------------------------------------------------------------------
void EMUAnalysis::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches or objects created by earlier modules.
  ReqBranch(fCaloJetsName,    fCaloJet);

  if (1)
    {
      string s;
      s = "run:evt:lumi:ptm:pte:etam:etae:phim:phie:dcam:dcae:chm:che:mvis:met:trigm:trige";
      
      fNt = new TNtuple("nt","EMUAnalysisMod",s.c_str());

      AddOutput(fNt);
    }
  // initialize some variables
  

  // request the following Event Objects (often branches, but could be locally created ones)

  // for MC only to adjust potential overlaps from generation

  // book all histograms

  AddTH1(fNAccCounters,"hEMUNAccCounters",";cut;#",15,-0.5,15.5);
  if (1) {
    TAxis *xa = fNAccCounters->GetXaxis();
    for(Int_t i=1;i<=fNAccCounters->GetNbinsX();++i)
      xa->SetBinLabel(i,"unused");
    xa->SetBinLabel(1,"a");
    xa->SetBinLabel(2,"b");
    xa->SetBinLabel(3,"c");
    xa->SetBinLabel(4,"d");
    xa->SetBinLabel(5,"e");
    xa->SetBinLabel(6,"f");
    xa->SetBinLabel(7,"g");
    xa->SetBinLabel(8,"h");
    xa->SetBinLabel(9,"i");
    xa->SetBinLabel(10,"j");
    xa->SetBinLabel(11,"k");
    xa->SetBinLabel(12,"l");
    xa->SetBinLabel(13,"m");
    xa->SetBinLabel(14,"n");
    xa->SetRangeUser(0,14);
  }

    // lepton histograms
  AddTH1(fptmtrig,"hEMUptmtrig","",400,0,200);
  AddTH1(fptmnotrig,"hEMUptmnotrig","",400,0,200);
  AddTH1(fptm,"hEMUptm","",400,0,200);
  AddTH1(fpte,"hEMUpte","",400,0,200);
  AddTH1(fetam,"hEMUetam","",400,-4,4);
  AddTH1(fetae,"hEMUetae","",400,-4,4);
  AddTH1(fphim,"hEMUphim","",400,(-1)*TMath::Pi(),TMath::Pi());
  AddTH1(fphie,"hEMUphie","",400,(-1)*TMath::Pi(),TMath::Pi());
  AddTH1(fdcam,"hEMUdcam","",600,-0.02,0.1);
  AddTH1(fdcae,"hEMUdcae","",600,-0.02,0.1);
  AddTH1(fchargem,"hEMUchargem","",3,-1.5,1.5);
  AddTH1(fchargee,"hEMUchargee","",3,-1.5,1.5);
  
  // di-tau histograms
  AddTH1(frecoMass,"hEMUrecoMass",";Reco Mass [GeV];#",400,0,200);     
  AddTH1(ftransMass,"hEMUtransMass",";transverse Mass [GeV];#",400,0,200);
  AddTH1(ftransEll,"hEMUtransEll",";transverse Ell [GeV];#",400,0,200); 
  AddTH1(ftransEnn,"hEMUtransEnn",";transverse Enn [GeV];#",400,0,200);
  AddTH1(ftransMnmu,"hEMUtransMnmu",";transverse nmu [GeV];#",400,0,200); 
  AddTH1(ftransMnel,"hEMUtransMnel",";transverse nel [GeV];#",400,0,200);

  AddTH1(fproj,   "hEMUproj",   "projected higgs pT [GeV];#",400,-50,50);
  AddTH1(fprojVis,"hEMUprojVis","projected vis pT [GeV];#",400,-50,50);
  AddTH1(fprojMet,"hEMUprojMet","projected met [GeV];#",400,-50,50);
  AddTH1(fprojPhi,"hEMUprojPhi","projection",400,(-2)*TMath::Pi(),2*TMath::Pi());
  AddTH1(fhT,     "hEMUhT",     "hT [GeV];#",600,0,300);

  AddTH1(fvisMass,"hEMUvisMass",";visible mass [GeV];#",600,0,300); 
  AddTH1(fscaledVisMass,"hEMUscaledVisMass",";scaled visible mass [GeV];#",600,0,300); 
  AddTH1(fxtaue,"hEMUxtaum","xe",400,-4,4);    
  AddTH1(fxtaum,"hEMUxtaue","xm",400,-4,4);
  AddTH1(fdphi,"hEMUdphi",";#Delta #phi;#",360,0,180);
  AddTH1(fdphinmu,"hEMUdphinmu",";#Delta #phi;#",360,0,180);
  AddTH1(fdphinel,"hEMUdphinel",";#Delta #phi;#",360,0,180);
  AddTH1(fdeta,"hEMUdeta",";#Delta #eta;#",400,0,5);   
  AddTH1(fcharge,"hEMUcharge","charge",3,-1.5,1.5);   
  AddTH1(ftype,"hEMUtype","type",10,-0.5,9.5);   

    // event histograms
  AddTH1(fnleptons,"hEMUnleptons",";leptons;#",10,-0.5,9.5);
  AddTH1(fnmuons,"hEMUnmuons",";leptons;#",10,-0.5,9.5);
  AddTH1(fnelecs,"hEMUnelecs",";leptons;#",10,-0.5,9.5);
  AddTH1(fnjets,"hEMUnjets",";jets;#",10,-0.5,9.5);
  AddTH1(fnjetspt15,"hEMUnjetspt15",";jets;#",10,-0.5,9.5);
  AddTH1(fnjetspt20,"hEMUnjetspt20",";jets;#",10,-0.5,9.5);
  AddTH1(fnbjetspt20sv2,"hEMUnbjetspt20sv2",";jets;#",10,-0.5,9.5);
  AddTH1(fnbjetspt20sv05,"hEMUnbjetspt20sv05",";jets;#",10,-0.5,9.5);
  AddTH1(fnjetspt25,"hEMUnjetspt25",";jets;#",10,-0.5,9.5);
  AddTH1(fmet,"hEMUmet","",200,0,200);
  AddTH1(fbtag,"hEMUbtag","",400,0.,20.);
  AddTH1(fmmmass,"hEMUmmmass",";visible mass [GeV];#",600,0,300); 
  AddTH1(feemass,"hEMUeemass",";visible mass [GeV];#",600,0,300); 

  AddTH1(fvisMassCut1,"hHTTvisMassCut1",";visible mass [GeV];#",600,0,300); 
  AddTH1(fvisMassCut2,"hHTTvisMassCut2",";visible mass [GeV];#",600,0,300); 
  AddTH1(fvisMassCut3,"hHTTvisMassCut3",";visible mass [GeV];#",600,0,300); 
  AddTH1(fvisMassCut4,"hHTTvisMassCut4",";visible mass [GeV];#",600,0,300); 
  AddTH1(fvisMassCut5,"hHTTvisMassCut5",";visible mass [GeV];#",600,0,300); 
  AddTH1(fvisMassCut6,"hHTTvisMassCut6",";visible mass [GeV];#",600,0,300); 
  AddTH1(fvisMassCut7,"hHTTvisMassCut7",";visible mass [GeV];#",600,0,300); 
  AddTH1(fvisMassCut8,"hHTTvisMassCut8",";visible mass [GeV];#",600,0,300); 
 
}

//--------------------------------------------------------------------------------------------------
void EMUAnalysis::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}

