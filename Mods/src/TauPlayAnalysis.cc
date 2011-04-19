// $Id: TauPlayAnalysis.cc,v 1.1 2011/02/18 14:56:34 dkralph Exp $

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
  fTausName("HPSTaus"),
  tauPt(0),
  tauEta(0)
{
  // Constructor
}
//--------------------------------------------------------------------------------------------------
TauPlayAnalysis::~TauPlayAnalysis()
{
}
//--------------------------------------------------------------------------------------------------
void TauPlayAnalysis::Process()
{
  // This is the function that gets called on every event in the file

  LoadEventObject( fTausName,             fTaus);

  assert(fTaus);

  // make a cut on the number of taus
  if(fTaus->GetEntries() < 1)
    return;

  // fill a histogram with the tau pt and eta
  for(UInt_t i=0;i<fTaus->GetEntries();i++) {
    tauPt->Fill(fTaus->At(i)->Pt());
    tauEta->Fill(fTaus->At(i)->Eta());
  }

  return;
}
//--------------------------------------------------------------------------------------------------
void TauPlayAnalysis::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis.
  AddTH1(tauPt,"hTauPlay_TauPt","",400,0,200);
  AddTH1(tauEta,"hTauPlay_TauEta","",400,-5,5);
}
