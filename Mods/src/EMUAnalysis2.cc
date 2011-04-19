// $Id: $

#include <TMath.h>
#include <TH1D.h>
#include <TString.h>
#include <TTree.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/PFTauFwd.h"
#include "MitAna/DataTree/interface/MetCol.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitPhysics/Utils/interface/DiTauSystem.h"

#include "MitHtt/Mods/interface/EMUAnalysis2.h"

using namespace mithep;

ClassImp(mithep::EMUAnalysis2)

//--------------------------------------------------------------------------------------------------
EMUAnalysis2::EMUAnalysis2(const char *name, const char *title) :
  BaseMod             (name,title),
  fTrigObjsName(0),
  fMuonsName(0),
  fElectronsName(0),
  fJetsName(0),
  fCaloJetsName("AKt5Jets"),
  fMetsName("PubPFMet"),
  fMuons(0),
  fElectrons(0),
  fJets(0),
  fCaloJets(0),
  fMets(0),
  fMcEventInfo(0),
  cutPtLeadingMuon(15.),
  cutPtLeadingElec(15.),
  cutPtSecondMuon(15.),
  cutPtSecondElec(15.),
  cutPtTriggerMuon(15.),
  fHistNamePref("hEMU_")
{
  // Constructor.
  fCutNames.push_back("one mu");
  fCutNames.push_back("one elec");
  fCutNames.push_back("OS");
  fCutNames.push_back("matchTrigMu");

  fNCuts = fCutNames.size();  
  fPassCuts = vector<bool> (fNCuts);
  for(UInt_t i=0;i<fNCuts;i++) fCutMap[fCutNames[i]] = i;

  //----------------------------------------------------------------------------------------
  fInfoV.push_back(new PlotInfo("AccCount","",                         fNCuts+1,-1.5,fNCuts-0.5));
  // lepton histograms
  fInfoV.push_back(new PlotInfo("ptm","",                    400,0,200));
  fInfoV.push_back(new PlotInfo("pte","",                    400,0,200));
  fInfoV.push_back(new PlotInfo("etam","",                   400,-4,4));
  fInfoV.push_back(new PlotInfo("etae","",                   400,-4,4));
  fInfoV.push_back(new PlotInfo("phim","",                   400,(-1)*TMath::Pi(),TMath::Pi()));
  fInfoV.push_back(new PlotInfo("phie","",                   400,(-1)*TMath::Pi(),TMath::Pi()));
  fInfoV.push_back(new PlotInfo("dcam","",                   600,-0.02,0.1));
  fInfoV.push_back(new PlotInfo("dcae","",                   600,-0.02,0.1));
  fInfoV.push_back(new PlotInfo("chargem","",                3,-1.5,1.5));
  fInfoV.push_back(new PlotInfo("chargee","",                3,-1.5,1.5));
  // di-tau histograms
  fInfoV.push_back(new PlotInfo("recoMass","Reco Mass [GeV]",       400,0,300));
  fInfoV.push_back(new PlotInfo("transMass","m_t [GeV]",            400,0,300));
  fInfoV.push_back(new PlotInfo("transEll","transverse Ell [GeV]",  400,0,300)); 
  fInfoV.push_back(new PlotInfo("transEnn","transverse Enn [GeV]",  400,0,300));
  fInfoV.push_back(new PlotInfo("transMnmu","transverse nmu [GeV]", 400,0,300)); 
  fInfoV.push_back(new PlotInfo("transMnel","transverse nel [GeV]", 400,0,300));
//   fInfoV.push_back(new PlotInfo("proj","projected higgs pT [GeV]",  400,-50,50));
//   fInfoV.push_back(new PlotInfo("projVis","projected vis pT [GeV]", 400,-50,50));
//   fInfoV.push_back(new PlotInfo("projMet","projected met [GeV]",    400,-50,50));
//   fInfoV.push_back(new PlotInfo("projPhi","projection",             400,(-2)*TMath::Pi(),2*TMath::Pi()));
//   fInfoV.push_back(new PlotInfo("hT","hT [GeV]",                    600,0,300));
  fInfoV.push_back(new PlotInfo("visMass","visible mass [GeV]",     600,0,300)); 
//   fInfoV.push_back(new PlotInfo("scaledVisMass","scaled visible mass [GeV]",600,0,300)); 
  fInfoV.push_back(new PlotInfo("xtaum","xe",                       400,-4,4));
  fInfoV.push_back(new PlotInfo("xtaue","xm",                       400,-4,4));
  fInfoV.push_back(new PlotInfo("dphi","#Delta #phi",               360,0,180));
  fInfoV.push_back(new PlotInfo("dphinmu","#Delta #phi;#",          360,0,180));
  fInfoV.push_back(new PlotInfo("dphinel","#Delta #phi",            360,0,180));
  fInfoV.push_back(new PlotInfo("deta","#Delta #eta",               400,0,5));
  fInfoV.push_back(new PlotInfo("charge","charge",                  3,-1.5,1.5));
  // event histograms
  fInfoV.push_back(new PlotInfo("nleptons","leptons",               10,-0.5,9.5));
  fInfoV.push_back(new PlotInfo("nmuons","leptons",                 10,-0.5,9.5));
  fInfoV.push_back(new PlotInfo("nelecs","leptons",                 10,-0.5,9.5));
  fInfoV.push_back(new PlotInfo("njets","jets",                     10,-0.5,9.5));
  fInfoV.push_back(new PlotInfo("met","",                           200,0,200));
  fInfoV.push_back(new PlotInfo("mmmass","visible mass [GeV]",      600,0,300));
  fInfoV.push_back(new PlotInfo("eemass","visible mass [GeV]",      600,0,300)); 

  fInfoV.push_back(new PlotInfo("visMassCut1","m_vis [GeV]",    600,0,300)); 
  fInfoV.push_back(new PlotInfo("visMassCut2","m_vis [GeV]",    600,0,300)); 
  fInfoV.push_back(new PlotInfo("visMassCut3","m_vis [GeV]",    600,0,300)); 
  fInfoV.push_back(new PlotInfo("ttbarcontrol","m_vis [GeV]",   600,0,300)); 

  fNHists = fInfoV.size();
  fHistV = vector<TH1D*> (fNHists);
  for(UInt_t i=0;i<fNHists;i++) fHistMap[fInfoV[i]->fname] = i;
}
//--------------------------------------------------------------------------------------------------
EMUAnalysis2::~EMUAnalysis2()
{
  for(UInt_t i=0;i<fNHists;i++) delete fInfoV[i];
}
//--------------------------------------------------------------------------------------------------
void EMUAnalysis2::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
  ReqEventObject(fMuonsName,     fMuons, 	kFALSE);
  ReqEventObject(fElectronsName, fElectrons, 	kFALSE);
  ReqEventObject(fJetsName,      fJets, 	kFALSE);
  ReqEventObject(fCaloJetsName,  fCaloJets, 	kTRUE);
  ReqEventObject(fMetsName,      fMets, 	kFALSE);

//   fTr = new TTree("WWTauTree","WWTauTree");
//   fTr->Branch("lepPt",&fLepPt);

//   AddOutput(fTr);
  //----------------------------------------------------------------------------------------
  for(UInt_t i=0;i<fNHists;i++) {
    AddTH1(fHistV[i],(fHistNamePref + TString(fInfoV[i]->fname)).Data(),
	   fInfoV[i]->ftitle.c_str(),fInfoV[i]->fnbins,fInfoV[i]->fxmin,fInfoV[i]->fxmax);
  }

}

//--------------------------------------------------------------------------------------------------
void EMUAnalysis2::Process()
{
  fHistV[fHistMap["AccCount"]]->GetXaxis()->SetBinLabel(1,"start");
  for(UInt_t i=0;i<fNCuts;i++) {
    fPassCuts[i] = false;
    fHistV[fHistMap["AccCount"]]->GetXaxis()->SetBinLabel(i+2,fCutNames[i].c_str());
  }

  // count the events we have processed
  IncNEventsProcessed();
  Fill("AccCount",-1);

  LoadEventObject( fMuonsName,            fMuons);
  LoadEventObject( fElectronsName,        fElectrons);
  LoadEventObject( fJetsName,             fJets);
  LoadEventObject( fCaloJetsName,         fCaloJets);
  LoadEventObject( fMetsName,             fMets);
  const TriggerObjectCol *tos = GetHLTObjects(fTrigObjsName);

  if(!fMuons || !fElectrons || !fJets || !fCaloJets || !fMets || !tos) {
    SendError(kAbortAnalysis,"Process","null pointer to collection!");
    return;
  }

  UInt_t nMuons   = fMuons->GetEntries();
  UInt_t nElecs   = fElectrons->GetEntries();
  UInt_t nJets    = fJets->GetEntries();
  UInt_t nTos = tos->GetEntries();

  // dimuon and dielectron mass if there is two el or two mu
  if(nMuons>1) {
    if(fMuons->At(0)->Pt()>cutPtLeadingMuon && fMuons->At(1)->Pt()>cutPtSecondMuon) {
      CompositeParticle dil;
      dil.AddDaughter(fMuons->At(0));
      dil.AddDaughter(fMuons->At(1));
      Fill("mmmass",dil.Mass());
    }
  }
  if(nElecs>1) {
    if(fElectrons->At(0)->Pt()>cutPtLeadingElec && fElectrons->At(1)->Pt()>cutPtSecondElec) {
      CompositeParticle dil;
      dil.AddDaughter(fElectrons->At(0));
      dil.AddDaughter(fElectrons->At(1));
      Fill("eemass",dil.Mass());
    }
  }

  // preselection
  UInt_t nPreCuts=2;
  if(nMuons>0)    fPassCuts[fCutMap["one mu"]]   = true;
  if(nElecs>0)    fPassCuts[fCutMap["one elec"]] = true;

  bool passPreCuts = true;
  for(UInt_t i=0;i<nPreCuts;i++) {
    passPreCuts = passPreCuts && fPassCuts[i];
    if(passPreCuts) Fill("AccCount",i);
  }
  if(!passPreCuts)
    return;

  //-------------------------------------------------------------------------------------------
  // note: code from here down will not execute if preseletions fail
  //-------------------------------------------------------------------------------------------
  // from now on, assume we have a mu and an el
  const Muon *mu = fMuons->At(0);
  const Electron *el = fElectrons->At(0);
  const Met *met = fMets->At(0);

  // leading muon fired trigger 15 GeV 
  UInt_t nMatchedTrigMuons = 0; // N tos with pt>pt0 that match the lead mu
  
  vector<string> trigModNames;
  trigModNames.push_back(string("hltSingleMu9L3Filtered9"));
  trigModNames.push_back(string("hltSingleMu15L3Filtered15"));
  for(UInt_t j=0; j<nTos; ++j) {
    const TriggerObject *to = tos->At(j);
    bool match=false;
    // make sure t.o. has right module name
    for(UInt_t in=0;in<trigModNames.size();in++) {
      if(trigModNames[in].compare(to->ModuleName()) == 0) match=true;
    }
    if(!match) continue;
    // match to muons
    if(to->Pt() > cutPtTriggerMuon) {
      double dEta = fabs(mu->Eta()-to->Eta());
      double dPhi = fabs(mu->Phi()-to->Phi());
      if(dEta<0.20 && dPhi<0.10) nMatchedTrigMuons++;
    }
  }
  if(nMatchedTrigMuons<1) {
    mu->Print();
    for(UInt_t i=0;i<nTos;i++) {
      tos->At(i)->Print();
    }
    cout << "=======================" << endl;
  }

  //----------------------------------------------------------------------------------------
  // cuts
  if(mu->Charge()*el->Charge() < 0.)    fPassCuts[fCutMap["OS"]]          = true;
  if(nMatchedTrigMuons > 0)             fPassCuts[fCutMap["matchTrigMu"]] = true;

  //----------------------------------------------------------------------------------------
  // "N-minus-one" plots 

  //----------------------------------------------------------------------------------------
  // "after" plots
  bool passAllCuts = passPreCuts;
  for(UInt_t i=nPreCuts;i<fNCuts;i++) {
    passAllCuts = passAllCuts && fPassCuts[i];
    if(passAllCuts) Fill("AccCount",i);
  }

  if(passAllCuts) {
    Fill("ptm",         mu->Pt());
    Fill("pte",         el->Pt());
    Fill("etam",        mu->Eta());
    Fill("etae",        el->Eta());
    Fill("phim",        mu->Phi());
    Fill("phie",        el->Phi());
    Fill("dcam",        mu->D0PV());
    Fill("dcae",        el->D0PV());
    Fill("chargem",     mu->Charge());
    Fill("chargee",     el->Charge());

    DiTauSystem diTau(mu,el,met);
    Fill("recoMass",    diTau.RecoMass());
    Fill("transMass",  	diTau.TransverseMass());
    Fill("transEll",    diTau.TransverseEll());
    Fill("transEnn",    diTau.TransverseEnn());
//     Fill("proj",        diTau.Projected());
//     Fill("projVis",     diTau.ProjectedVis());
//     Fill("projMet",     diTau.ProjectedMet());
//     Fill("projPhi",     diTau.ProjectedPhi());
//     Fill("hT",          diTau.Ht());
    Fill("visMass",     diTau.VisMass());
    Fill("xtaum",       diTau.XTau1());
    Fill("xtaue", 	diTau.XTau2());

    double deltaPhiMetLepton[2] = {MathUtils::DeltaPhi(met->Phi(), mu->Phi()),
				   MathUtils::DeltaPhi(met->Phi(), el->Phi())};
    double mTW[2] = {sqrt( 2.0*mu->Pt()*met->Pt()*(1.0 - cos(deltaPhiMetLepton[0])) ),
		     sqrt( 2.0*el->Pt()*met->Pt()*(1.0 - cos(deltaPhiMetLepton[1])) )};
    Fill("transMnmu",   mTW[0]);
    Fill("transMnel",   mTW[1]);
    Fill("dphinmu",     deltaPhiMetLepton[0]*180./TMath::Pi());
    Fill("dphinel",     deltaPhiMetLepton[1]*180./TMath::Pi());

    double deltaPhiLeptons  = MathUtils::DeltaPhi(mu->Phi(),el->Phi()) * 180./TMath::Pi();  
    double deltaEtaLeptons  = fabs(mu->Eta() - el->Eta());
    Fill("dphi",           deltaPhiLeptons);
    Fill("deta",           deltaEtaLeptons);
    Fill("charge",         mu->Charge()*el->Charge());
    Fill("nleptons",       nMuons+nElecs);
    Fill("nmuons",         nMuons);
    Fill("nelecs",         nElecs);
    Fill("njets",          nJets);        // total pf jets
    Fill("met",            met->Pt());

    // Higgs selection cuts
    if(mTW[0]<60. && mTW[1]<60.)     		Fill("visMassCut1",  diTau.VisMass());
    if(mTW[0]<50. && mTW[1]<50.)     		Fill("visMassCut2",  diTau.VisMass());
//     if(diTau.ProjectedVis() < 20.)  		Fill("visMassCut3",  diTau.VisMass());
    // ttbar control region
    if(mTW[0]>50. && mTW[1]>50.)                Fill("ttbarcontrol", diTau.VisMass());
    // !!! Pzeta - 1.5 Pzeta cut?
    // !!! Z removal
  }

  return;
}
//--------------------------------------------------------------------------------------------------
bool EMUAnalysis2::NMinusOnePass(const char *cut)
{
  bool pass=true;
  for(UInt_t i=0;i<fCutNames.size();i++) {
    if(TString(cut) == TString(fCutNames[i])) {
      for(UInt_t j=0;j<fNCuts;j++) {
	if(j != fCutMap[cut])
	  pass = pass && fPassCuts[j];
      }
      return pass;
    }
  }

  // shouldn't get here
  SendError(kAbortAnalysis,"Process","Error: bad cut name given.");
  return false;
  
}
//--------------------------------------------------------------------------------------------------
bool EMUAnalysis2::NMinusThreePass(const char *cut1,const char *cut2,const char *cut3)
{
  bool pass=true;
  for(UInt_t j=0;j<fNCuts;j++) {
    if(j!=fCutMap[cut1] && j!=fCutMap[cut2] && j!=fCutMap[cut3])
      pass = pass && fPassCuts[j];
  }

  return pass;
}
//--------------------------------------------------------------------------------------------------
void EMUAnalysis2::Fill(const char *hname, double val)
{
  for(UInt_t i=0;i<fNHists;i++) {
    if(TString(fInfoV[i]->fname) == TString(hname)) {
      fHistV[fHistMap[hname]]->Fill(val);
      return;
    }
  }

  SendError(kAbortAnalysis,"Process",Form("Error: histogram %s not found!",hname));
  return;
}
