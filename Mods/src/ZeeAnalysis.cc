// $Id: $

#include <TMath.h>
#include <TH1D.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"

#include "MitHtt/Mods/interface/ZeeAnalysis.h"

using namespace mithep;

ClassImp(mithep::ZeeAnalysis)

//--------------------------------------------------------------------------------------------------
ZeeAnalysis::ZeeAnalysis(const char *name, const char *title) :
  BaseMod             (name,title),
  // initialize bambu objects
  fTrigObjsName       (Names::gkHltObjBrn),
  fElecName           (Names::gkElectronBrn),
  fElecsFromBranch    (kTRUE),
  fMcEventInfoName    (Names::gkMCEvtInfoBrn),
  fElecs              (0),
  fMcEventInfo        (0),
  fIsData             (kTRUE),
  fSkimming           (kFALSE),
  // scale removal
  fOverlapCut         (-1.0),
  // cuts for selection
  fElecEtMin          (20.0),
  //fNoPixelSeedCut     (kTRUE),
  fElecSigEeBarMax    ( 0.0105),
  fElecSigEeEcpMax    ( 0.0300),
  fElecIsoMax         ( 0.30),
  // cuts for selection
  f2ElecMassMin       (140.0),
  //fElecR9Min          (0.9),
  // initialize the histograms
  hAllNElecs          (0),
  hAllElecPt          (0),
  hAllElecEta         (0),
  //hAllElecPixelSeed   (0),
  hAllElecSigEeBar    (0),
  hAllElecSigEeEcp    (0),
  hAllElecEcHcTrIso   (0),
  hTrigDeltaEta       (0),
  hTrigDeltaPhi       (0),
  hTrigDeltaEt        (0),
  h2SeleElecMass      (0),
  hNElecs             (0),
  hElecEtMx           (0),
  hN2Elecs            (0),
  hElecEt1            (0),
  hElecEt2            (0),
  hElecEta1           (0),
  hElecEta2           (0),
  hElecPhi1           (0),
  hElecPhi2           (0),
  //hElecR91            (0),
  //hElecR92            (0),
  hElecDelR           (0),
  h2ElecEtSum         (0),
  h2ElecPt            (0),
  h2ElecMass          (0),
  //hHmElecR91          (0),
  //hHmElecR92          (0),
  //h2R9ElecMass        (0),
  hNTrigElecs         (0),
  h2TrigElecMass      (0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ZeeAnalysis::Begin()
{
  // Run startup code on the client machine. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void ZeeAnalysis::Process()
{
  // Process entries of the tree. For this module, we just load the branches and fill the
  // histograms.

  // count the events we have processed
  IncNEventsProcessed();

  // access the Mc Information if needed
  if (! fIsData) {
    LoadBranch(fMcEventInfoName);
    if (fOverlapCut > 0. && fMcEventInfo->Scale() > fOverlapCut) {
      MDB(kModules,1)
	printf(" Reject event with scale %f  above cut: %f\n",fMcEventInfo->Scale(),fOverlapCut);
      return;
    }
  }

  // access the elecs
  LoadEventObject(fElecName,fElecs);

  // monitor how many elecs have been found in this event (no additional selection)
  hAllNElecs->Fill(fElecs->GetEntries());

  // loop through them and make some pictures and select the ones relevant to our analysis
  double lastEt = 99999.0;
  ResetSelectedElecs();
  for (UInt_t i=0; i<fElecs->GetEntries(); ++i) {
    const Electron *p = fElecs->At(i);
    hAllElecPt ->Fill(p->Pt());
    hAllElecEta->Fill(p->Eta());

    double sigEtaEta = p->CoviEtaiEta();
    if      (p->IsEB())
      hAllElecSigEeBar->Fill(sigEtaEta);
    else if (p->IsEE())
      hAllElecSigEeEcp->Fill(sigEtaEta);
    double isoTot = (p->TrackIsolationDr03() + p->EcalRecHitIsoDr03() + p->HcalTowerSumEtDr03())
                    / p->Pt();
    hAllElecEcHcTrIso->Fill(isoTot);
    //hAllElecPixelSeed->Fill((int) p->HasPixelSeed());

    // check the elecs are really order according to pt
    if (lastEt <= p->Et()) {
      printf("\n WARNING -- elecs are out of order (%f larger  %f)\n\n",p->Et(),lastEt);
    }

    // apply elec selection (just an Et cut)
    if (p->Et() > fElecEtMin                   &&         // minimum Et
	//	(! p->HasPixelSeed() && fNoPixelSeedCut) &&         // no pixel seed
	isoTot  < fElecIsoMax                  &&         // maximum isolation
	((p->IsEB() && sigEtaEta < fElecSigEeBarMax) ||   // shower shape (barrel/encap separate)
	 (p->IsEE() && sigEtaEta < fElecSigEeEcpMax) 	 ))
      fSelectedElecs.push_back(i);

    // make sure to keep the Et accounting up to speed
    lastEt = p->Et();
  }
  hNElecs->Fill(fSelectedElecs.size());

  MatchElecsToTrigger();
  hNTrigElecs->Fill(fTriggeredElecs.size());

  // loop through our selected elecs and make histograms only if one was matched to trigger
  UInt_t n2Elecs = 0;
  for (UInt_t i=0; i<fSelectedElecs.size(); i++) {
    const Electron *p1 = fElecs->At(fSelectedElecs[i]);

    // histogram maximum elec Et
    if (i==0)
      hElecEtMx->Fill(p1->Et());

    // is this a triggered elec?
    Bool_t isP1Trig = ElecTriggered(fSelectedElecs[i]);

    // loop through all pairs
    const FourVectorM m1 = p1->Mom();
    for (UInt_t j=i+1; j<fSelectedElecs.size(); ++j) {
      const Electron *p2 = fElecs->At(fSelectedElecs[j]);

      // is this a triggered elec?
      Bool_t isP2Trig = ElecTriggered(fSelectedElecs[j]);

	const FourVectorM m2 = p2->Mom();
	double mass = (m1+m2).M();

      h2SeleElecMass ->Fill(mass);

      // only consider pairs, where at least one elec triggered
      if (isP1Trig || isP2Trig) {
        hElecEt1   ->Fill(p1->Et());
        hElecEt2   ->Fill(p2->Et());
        hElecEta1  ->Fill(p1->Eta());
        hElecEta2  ->Fill(p2->Eta());
        hElecPhi1  ->Fill(p1->Phi());
        hElecPhi2  ->Fill(p2->Phi());
        //hElecR91   ->Fill(p1->R9());
        //hElecR92   ->Fill(p2->R9());
        hElecDelR  ->Fill(MathUtils::DeltaR(m1,m2));
        h2ElecEtSum->Fill(p1->Et()+p2->Et());
        h2ElecPt   ->Fill((m1+m2).Pt());
        h2ElecMass ->Fill(mass);
	  //// our high mass selection
	  //if (mass > f2ElecMassMin) {
	  //  hHmElecR91->Fill(p1->R9());
	  //  hHmElecR92->Fill(p2->R9());
	  //}
	  //// our R9 cleaned elec selection
	  //if (p1->R9() > fElecR9Min && p2->R9() > fElecR9Min) {
	  //  h2R9ElecMass->Fill(mass);
	  //}
        // increment the number of good di-elec candidates
        n2Elecs++;
      }
    }
  }
  hN2Elecs->Fill(n2Elecs);

  // work only with triggered elecs
  for (UInt_t i=0; i<fTriggeredElecs.size(); i++) {
    const Electron *p1 = fElecs->At(fTriggeredElecs[i]);

    // loop through all pairs
    const FourVectorM m1 = p1->Mom();
    for (UInt_t j=i+1; j<fTriggeredElecs.size(); ++j) {
      const Electron *p2 = fElecs->At(fTriggeredElecs[j]);
      const FourVectorM m2 = p2->Mom();
      h2TrigElecMass ->Fill((m1+m2).M());
    }
  }

  if (fTriggeredElecs.size()<=0) {
    if (fSkimming)
      SkipEvent();
  }

  return;
}

//--------------------------------------------------------------------------------------------------
void ZeeAnalysis::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches or objects created by earlier modules.

  // initialize some variables
  double pi = TMath::Pi();
  double mPi = -1.0 * TMath::Pi();

  // request the following Event Objects (often branches, but could be locally created ones)
  ReqEventObject(fElecName,fElecs,fElecsFromBranch);

  // for MC only to adjust potential overlaps from generation
  if (! fIsData) {
    ReqBranch(Names::gkMCEvtInfoBrn,fMcEventInfo);
    printf(" Monte Carlo Information block %p\n",(void*) fMcEventInfo);
      printf(" --> this is no data. Access the McEventInfo.\n\n");
  }
  else
    printf(" --> this is data. Drop the McEventInfo.\n\n");

  // book all histograms
  AddTH1(hAllNElecs,       "hAllNElecs",";# of elecs;#",                         21,-0.5, 20.5);
  AddTH1(hAllElecPt,       "hAllElecPt",";p_{T} [GeV];#",                       100, 0.0,250.0);
  AddTH1(hAllElecEta,      "hAllElecEta",";#eta;#",                             160,-8.0,  8.0);
  AddTH1(hAllElecSigEeBar, "hElecSigEeBar", ";#sigma_{#eta#eta}, barrel;#",     200, 0.0,  0.4);
  AddTH1(hAllElecSigEeEcp, "hElecSigEeEcp", ";#sigma_{#eta#eta}, endcap;#",     200, 0.0,  0.4);
  AddTH1(hAllElecEcHcTrIso,"hElecEcHcTrIso", ";#isolation_{ECAL,HCAL,Trk};#",   100, 0.0,  1.0);
  //AddTH1(hAllElecPixelSeed,"hElecPixelSeed", ";#pixel seed;#",                    2,-0.5,  1.5);

  AddTH1(hTrigDeltaEta,    "hTrigDeltaEta",";#Delta #eta;#",                    100,-0.1,  0.1);
  AddTH1(hTrigDeltaPhi,    "hTrigDeltaPhi",";#Delta #phi;#",                    100,-0.05,0.05);
  AddTH1(hTrigDeltaEt,     "hTrigDeltaEt",";#Delta E_{T}/E_{T};#",              100,-0.1,  0.1);

  AddTH1(h2SeleElecMass,   "h2SeleElecMass",";Dielec Mass [GeV/c^{2}];#",       400, 0.0,400.0);

  AddTH1(hNElecs,          "hNElecs",";# of sel elecs;#",                        21,-0.5, 20.5);
  AddTH1(hElecEtMx,        "hElecEtMx", ";E_{T,max};#",                         200, 0.0,400.0);

  AddTH1(hN2Elecs,         "hN2Elecs",";# of di-elecs;#",                        21,-0.5, 20.5);
  AddTH1(hElecEt1,         "hElecEt1", ";E_{T,1};#",                            200, 0.0,400.0);
  AddTH1(hElecEt2,         "hElecEt2", ";E_{T,2};#",                            200, 0.0,400.0);
  AddTH1(hElecEta1,        "hElecEta1", ";#eta_{1};#",                          104,-2.6,  2.6);
  AddTH1(hElecEta2,        "hElecEta2", ";#eta_{2};#",                          104,-2.6,  2.6);
  AddTH1(hElecPhi1,        "hElecPhi1", ";#phi_{1};#",                          100, mPi,   pi);
  AddTH1(hElecPhi2,        "hElecPhi2", ";#phi_{2};#",                          100, mPi,   pi);
  //AddTH1(hElecR91,         "hElecR91", ";#R9_{1};#",                            100, 0.1,  1.1);
  //AddTH1(hElecR92,         "hElecR92", ";#R9_{2};#",                            100, 0.1,  1.1);
  AddTH1(hElecDelR,        "hElecDelR", ";#Delta R_{1,2};#",                    104, 0.0,  5.2);
  AddTH1(h2ElecEtSum,      "h2ElecEtSum", ";E_{T,sum};#",                       400, 0.0,400.0);
  AddTH1(h2ElecPt,         "h2ElecPt", ";pair P_{T};#",                         400, 0.0,400.0);
  AddTH1(h2ElecMass,       "h2ElecMass",";Dielec Mass [GeV/c^{2}];#",           400, 0.0,400.0);

  //AddTH1(hHmElecR91,       "hHmElecR91", ";#R9_{1};#",                          100, 0.1,  1.1);
  //AddTH1(hHmElecR92,       "hHmElecR92", ";#R9_{2};#",                          100, 0.1,  1.1);

  //AddTH1(h2R9ElecMass,     "h2R9ElecMass", ";Dielec Mass [GeV/c^{2}];#",        400, 0.0,400.0);

  AddTH1(hNTrigElecs,      "hNTrigElecs",";# of sel+trig elecs;#",               21,-0.5, 20.5);
  AddTH1(h2TrigElecMass,   "h2TrigElecMass",";Dielec Mass [GeV/c^{2}];#",       400, 0.0,400.0);

}

//--------------------------------------------------------------------------------------------------
void ZeeAnalysis::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void ZeeAnalysis::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ZeeAnalysis::ResetSelectedElecs()
{
  // Reset the selected and trigger matched elecs

  fSelectedElecs .clear();
  fTriggeredElecs.clear();

  return;
}

//--------------------------------------------------------------------------------------------------
void ZeeAnalysis::MatchElecsToTrigger()
{
  // Reset the selected and trigger matched elecs

  // access all trigger objects
  const TriggerObjectCol *tos = GetHLTObjects(fTrigObjsName);
  if (! tos)             // this can only happen if HLTMod::SetAbortIfNotAccepted(kFALSE) was called
    return;
  UInt_t nEnts = tos->GetEntries();

  // loop through our selected elecs
  for (UInt_t i=0; i<fSelectedElecs.size(); i++) {
    const Electron *p = fElecs->At(fSelectedElecs[i]);
    UInt_t matched = 0;
    // loop through all trigger objects and try to find a match
    for (UInt_t j=0; j<nEnts; ++j) {
      const TriggerObject *to = tos->At(j);
      if (fabs(p->Eta()-to->Eta())          < 0.10 &&
          fabs(p->Phi()-to->Phi())          < 0.05 &&
          fabs((p->Pt() -to->Pt())/p->Pt()) < 0.05   ) {
        matched++;
        hTrigDeltaEta->Fill(p->Eta()-to->Eta());
        hTrigDeltaPhi->Fill(p->Phi()-to->Phi());
        hTrigDeltaEt ->Fill((p->Pt()-to->Pt())/p->Pt());
      }
    }
    // add to our trigger collection if it was matched
    if (matched > 0)
      fTriggeredElecs.push_back(fSelectedElecs[i]);
    // print warning if desired
    if (matched > 1) {
      MDB(kGeneral,2) {
        printf(" Ambiguos match found: %d\n",matched);
        printf(" Elec:  %8.3f %8.3f %8.3f\n",p->Eta(), p->Phi(), p->Pt());
        for (UInt_t j=0; j<nEnts; ++j) {
          const TriggerObject *to = tos->At(j);
          printf(" TrigObj: %8.3f %8.3f %8.3f\n",to->Eta(),to->Phi(),to->Pt());
        }
      }
    }
  }

  return;
}

//--------------------------------------------------------------------------------------------------
Bool_t ZeeAnalysis::ElecTriggered(UInt_t iElec)
{
  // Determine whether the requested elec index is in the triggered elec indices

  // per default the elec did not trigger
  Bool_t triggered = kFALSE;

  // loop through triggered elec indices
  for (UInt_t i=0; i<fTriggeredElecs.size(); i++) {
    if (fTriggeredElecs[i] == iElec) {
      // index was found -> elec triggered
      triggered = kTRUE;
      break;
    }
  }

  return triggered;
}
