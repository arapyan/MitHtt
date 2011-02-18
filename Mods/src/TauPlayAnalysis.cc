// $Id: $

#include <TMath.h>
#include <TH1D.h>
#include <TString.h>
#include "MitPhysics/Init/interface/ModNames.h"
#include "MitCommon/MathTools/interface/MathUtils.h"

#include "MitAna/DataCont/interface/ObjArray.h"
#include "MitAna/DataTree/interface/ParticleCol.h"
#include "MitAna/DataTree/interface/CompositeParticle.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PFTauCol.h"
#include "MitAna/DataTree/interface/MetCol.h"

#include "MitHtt/Mods/interface/TauPlayAnalysis.h"

using namespace mithep;

ClassImp(mithep::TauPlayAnalysis)

//--------------------------------------------------------------------------------------------------
TauPlayAnalysis::TauPlayAnalysis(const char *name, const char *title) :
  BaseMod             (name,title),
  // initialize bambu objects
  // fEventHeader(Names::gkEvtTreeName),
  fCleanLeptonsName(ModNames::gkMergedLeptonsName),
  fTrigObjsName(Names::gkHltObjBrn),
  fMuonsName(ModNames::gkCleanMuonsName),        
  fElectronsName(ModNames::gkCleanElectronsName),
  fTausName(ModNames::gkCleanPFTausName),
  fHPSTauCandidatesName("HPSTaus"),
  fMetColName("PubPFMet"),
  fJetsName("CleanPFJets"),
  fMCPartName(Names::gkMCPartBrn),
  fMcEventInfo(0),
  fLeptons(0),
  fHPSTauCandidates(0),
  fParticles(0),
  fNAccCounters(0)

{
  // Constructor
  finfov.push_back(new PlotInfo("NGenTaus",  "NGenTaus;N Gen Taus;", 50,  0,      4));
  finfov.push_back(new PlotInfo("NPFTaus",   "NPFTaus;N PF Taus;",        50,  0,      4));
  finfov.push_back(new PlotInfo("NHPSTaus",   "NHPSTaus;N HPS Taus;",        50,  0,      4));
  finfov.push_back(new PlotInfo("PtGenTau",  "PtGenTau;Gen Tau Pt [GeV];",                 400, 0,    200));
  finfov.push_back(new PlotInfo("PtPFTau",   "PtPFTau;PF Tau Pt [GeV];",                 400, 0,    200));
  finfov.push_back(new PlotInfo("PtHPSTau",   "PtHPSTau;HPS Tau Pt [GeV];",                 400, 0,    200));
  finfov.push_back(new PlotInfo("PtDiff",    "PtDiff;#Delta Pt HPS-Gen Taus;",                 400, 0,    100));
  finfov.push_back(new PlotInfo("Eff",       "Efficiency;NHPSTaus/NGenTaus;",                 200, -0.5,   2));

  fnhists = finfov.size();
  fHistsv = vector<TH1D*> (fnhists);
  for(UInt_t i=0;i<fnhists;i++) fmap[finfov[i]->fname] = i;
}
//--------------------------------------------------------------------------------------------------
TauPlayAnalysis::~TauPlayAnalysis()
{
  for(UInt_t i=0;i<fnhists;i++) delete finfov[i];
}
//--------------------------------------------------------------------------------------------------
void TauPlayAnalysis::Begin()
{
  // Run startup code on the client machine. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void TauPlayAnalysis::Process()
{
  // count the events we have processed
  IncNEventsProcessed();
  UInt_t iAccCounter = 0; 
  fNAccCounters->Fill(iAccCounter++);
  LoadEventObject(fMCPartName, fParticles);
  LoadEventObject(fTausName, fTaus);
  LoadEventObject(fHPSTauCandidatesName, fHPSTauCandidates);

  if (!fParticles || !fTaus || !fHPSTauCandidates) {
    cout << "null pointer to collection!" << endl;
    SkipEvent();
    return;
  }
  fNAccCounters->Fill(iAccCounter++);

  // find the taus that passed the hps discriminant
  PFTauOArr *HPSTaus = new PFTauOArr;
  for(UInt_t i=0;i<fHPSTauCandidates->GetEntries();i++) {
    double discr = fHPSTauCandidates->At(i)->DiscriminationByDecayModeFinding();
    if(discr>0.5) {
      HPSTaus->Add(fHPSTauCandidates->At(i));
    }
  }

  MCParticleOArr *hadrGenTaus = new MCParticleOArr;
  UInt_t nGenTaus=0;
  for(UInt_t i=0;i<fParticles->GetEntries();i++) {
    const MCParticle *p = fParticles->At(i);
    // hadronic taus:
    if(p->Is(MCParticle::kTau) && p->Status()==2) {
      if(!p->HasDaughter(MCParticle::kEl) && !p->HasDaughter(MCParticle::kMu) &&
	 p->HasDaughter(MCParticle::kTauNu)) {
	if(p->Pt()>15 && p->Eta()<2.5) {
	  nGenTaus++;
	  hadrGenTaus->Add(p);
	  fHistsv[fmap["PtGenTau"]]->Fill(p->Pt());
	}
      }
    }
  }
  UInt_t nPFTaus=fTaus->GetEntries();
  UInt_t nHPSTaus=HPSTaus->GetEntries();

  for(UInt_t i=0;i<nHPSTaus;i++) {
    const PFTau *p = HPSTaus->At(i);
    double minPtDiff=9999.;
    for(UInt_t j=0;j<hadrGenTaus->GetEntries();j++) {
      double difftmp = fabs(hadrGenTaus->At(j)->Pt() - p->Pt());
      if(difftmp<minPtDiff) {
	minPtDiff = difftmp;
      }
    }
    fHistsv[fmap["PtDiff"]]->Fill(minPtDiff);
    fHistsv[fmap["PtHPSTau"]]->Fill(p->Pt());
  }
  if(nGenTaus != 0) fHistsv[fmap["Eff"]]->Fill(double(nHPSTaus)/nGenTaus);

  for(UInt_t i=0;i<nPFTaus;i++) {
    const PFTau *p = fTaus->At(i);
    fHistsv[fmap["PtPFTau"]]->Fill(p->Pt());
  }

  fHistsv[fmap["NGenTaus"]]->Fill(nGenTaus);
  fHistsv[fmap["NPFTaus"]]->Fill(nPFTaus);
  fHistsv[fmap["NHPSTaus"]]->Fill(nHPSTaus);

  delete hadrGenTaus;
  delete HPSTaus;
  return;
}

//--------------------------------------------------------------------------------------------------
void TauPlayAnalysis::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis.
}

//--------------------------------------------------------------------------------------------------
void TauPlayAnalysis::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
  bool fIsData = kFALSE;
  if(fIsData == kFALSE){
    ReqEventObject(fMCPartName, fParticles, kTRUE);
  }
  ReqEventObject(fTausName, fTaus, kFALSE);
  ReqEventObject(fHPSTauCandidatesName, fHPSTauCandidates, kFALSE);

  AddTH1(fNAccCounters,"hTauPlay_AccCounters",";cut;#",15,-0.5,15.5);
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

  //----------------------------------------------------------------------------------------
  for(UInt_t i=0;i<fnhists;i++) {
    AddTH1(fHistsv[i],(TString("hTauPlay_") + TString(finfov[i]->fname)).Data(),finfov[i]->ftitle.c_str(),
	   finfov[i]->fnbins,finfov[i]->fxmin,finfov[i]->fxmax);
  }

}
//--------------------------------------------------------------------------------------------------
void TauPlayAnalysis::Terminate()
{
  // Run finishing code on the client computer.
}

