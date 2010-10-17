// $Id: $

#include <TMath.h>
#include <TH1D.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/ElectronCol.h"
#include "MitAna/DataTree/interface/PhotonCol.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"

#include "MitHtt/Mods/interface/HttAnalysis.h"

using namespace mithep;

ClassImp(mithep::HttAnalysis)

//--------------------------------------------------------------------------------------------------
HttAnalysis::HttAnalysis(const char *name, const char *title) :
  BaseMod             (name,title),
  // initialize bambu objects
  fTrigObjsName       (Names::gkHltObjBrn),
  fPhotonName         (Names::gkPhotonBrn),
  fPhotonsFromBranch  (kTRUE),
  fMcEventInfoName    (Names::gkMCEvtInfoBrn),
  fPhotons            (0),
  fMcEventInfo        (0),
  fIsData             (kTRUE),
  // scale removal
  fOverlapCut         (-1.0),
  // cuts for selection
  fPhotonEtMin        (20.0),
  fNoPixelSeedCut     (kTRUE),
  fPhotonSigEeBarMax  ( 0.0105),
  fPhotonSigEeEcpMax  ( 0.0300),
  fPhotonIsoMax       ( 0.30),
  // cuts for selection
  f2PhotonMassMin     (140.0),
  fPhotonR9Min        (0.9),
  // initialize the histograms
  hAllNPhotons        (0),
  hAllPhotonPt        (0),
  hAllPhotonEta       (0),
  hAllPhotonPixelSeed (0),
  hAllPhotonSigEeBar  (0),
  hAllPhotonSigEeEcp  (0),
  hAllPhotonEcHcTrIso (0),
  hTrigDeltaEta       (0),
  hTrigDeltaPhi       (0),
  hTrigDeltaEt        (0),
  hNPhotons           (0),
  hPhotonEtMx         (0),
  hN2Photons          (0),
  hPhotonEt1          (0),
  hPhotonEt2          (0),
  hPhotonEta1         (0),
  hPhotonEta2         (0),
  hPhotonPhi1         (0),
  hPhotonPhi2         (0),
  hPhotonR91          (0),
  hPhotonR92          (0),
  hPhotonDelR         (0),
  h2PhotonEtSum       (0),
  h2PhotonPt          (0),
  h2PhotonMass        (0),
  hHmPhotonR91        (0),
  hHmPhotonR92        (0),
  h2R9PhotonMass      (0),
  hNTrigPhotons       (0),
  h2TrigPhotonMass    (0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void HttAnalysis::Begin()
{
  // Run startup code on the client machine. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void HttAnalysis::Process()
{
  // Process entries of the tree. For this module, we just load the branches and fill the
  // histograms.

  // count the events we have processed
  IncNEventsProcessed();
  
  return;
   
  // access the Mc Information if needed
  if (! fIsData) {
    LoadBranch(fMcEventInfoName);
    if (fOverlapCut > 0. && fMcEventInfo->Scale() > fOverlapCut) {
      MDB(kModules,1)
	printf(" Reject event with scale %f  above cut: %f\n",fMcEventInfo->Scale(),fOverlapCut);
      return;
    }
  }

  // access the photons
  LoadEventObject(fPhotonName,fPhotons);

  // monitor how many photons have been found in this event (no additional selection)
  hAllNPhotons->Fill(fPhotons->GetEntries());

  // loop through them and make some pictures and select the ones relevant to our analysis
  double lastEt = 99999.0;
  ResetSelectedPhotons();
  for (UInt_t i=0; i<fPhotons->GetEntries(); ++i) {
    const Photon *p = fPhotons->At(i);
    hAllPhotonPt ->Fill(p->Pt());
    hAllPhotonEta->Fill(p->Eta());

    double sigEtaEta = p->SCluster()->Seed()->CoviEtaiEta();
    if      (p->IsEB())
      hAllPhotonSigEeBar->Fill(sigEtaEta);
    else if (p->IsEE())
      hAllPhotonSigEeEcp->Fill(sigEtaEta);
    double isoTot = (p->HollowConeTrkIsoDr04() + p->EcalRecHitIsoDr04() + p->HcalTowerSumEtDr04())
                    / p->Pt();
    hAllPhotonEcHcTrIso->Fill(isoTot);
    hAllPhotonPixelSeed->Fill((int) p->HasPixelSeed());

    // check the photons are really order according to pt
    if (lastEt <= p->Et()) {
      printf("\n WARNING -- photons are out of order (%f larger  %f)\n\n",p->Et(),lastEt);
    }

    // apply photon selection (just an Et cut)
    if (p->Et() > fPhotonEtMin                   &&         // minimum Et
	(! p->HasPixelSeed() && fNoPixelSeedCut) &&         // no pixel seed
	isoTot  < fPhotonIsoMax                  &&         // maximum isolation
	((p->IsEB() && sigEtaEta < fPhotonSigEeBarMax) ||   // shower shape (barrel/encap separate)
	 (p->IsEE() && sigEtaEta < fPhotonSigEeEcpMax) 	 ))
      fSelectedPhotons.push_back(i);

    // make sure to keep the Et accounting up to speed
    lastEt = p->Et();
  }
  hNPhotons->Fill(fSelectedPhotons.size());

  MatchPhotonsToTrigger();
  hNTrigPhotons->Fill(fTriggeredPhotons.size());

  // loop through our selected photons and make histograms only if one was matched to trigger
  if (fTriggeredPhotons.size()>0) {
    UInt_t n2Photons = 0;
    for (UInt_t i=0; i<fSelectedPhotons.size(); i++) {
      const Photon *p1 = fPhotons->At(fSelectedPhotons[i]);

      // histogram maximum photon Et
      if (i==0)
        hPhotonEtMx->Fill(p1->Et());

      // is this a triggered photon?
      Bool_t isP1Trig = PhotonTriggered(fSelectedPhotons[i]);

      // loop through all pairs
      const FourVectorM m1 = p1->Mom();
      for (UInt_t j=i+1; j<fSelectedPhotons.size(); ++j) {
        const Photon *p2 = fPhotons->At(fSelectedPhotons[j]);

        // is this a triggered photon?
        Bool_t isP2Trig = PhotonTriggered(fSelectedPhotons[j]);

        // only consider pairs, where at least one photon triggered
        if (isP1Trig || isP2Trig) {
          const FourVectorM m2 = p2->Mom();
	  double mass = (m1+m2).M();
          hPhotonEt1   ->Fill(p1->Et());
          hPhotonEt2   ->Fill(p2->Et());
          hPhotonEta1  ->Fill(p1->Eta());
          hPhotonEta2  ->Fill(p2->Eta());
          hPhotonPhi1  ->Fill(p1->Phi());
          hPhotonPhi2  ->Fill(p2->Phi());
          hPhotonR91   ->Fill(p1->R9());
          hPhotonR92   ->Fill(p2->R9());
          hPhotonDelR  ->Fill(MathUtils::DeltaR(m1,m2));
          h2PhotonEtSum->Fill(p1->Et()+p2->Et());
          h2PhotonPt   ->Fill((m1+m2).Pt());
          h2PhotonMass ->Fill(mass);
	  // our high mass selection
	  if (mass > f2PhotonMassMin) {
	    hHmPhotonR91->Fill(p1->R9());
	    hHmPhotonR92->Fill(p2->R9());
	  }
	  // our R9 cleaned photon selection
	  if (p1->R9() > fPhotonR9Min && p2->R9() > fPhotonR9Min) {
	    h2R9PhotonMass->Fill(mass);
	  }
          // increment the number of good di-photon candidates
          n2Photons++;
        }
      }
    }
    hN2Photons->Fill(n2Photons);

    // work only with triggered photons
    for (UInt_t i=0; i<fTriggeredPhotons.size(); i++) {
      const Photon *p1 = fPhotons->At(fTriggeredPhotons[i]);

      // loop through all pairs
      const FourVectorM m1 = p1->Mom();
      for (UInt_t j=i+1; j<fTriggeredPhotons.size(); ++j) {
        const Photon *p2 = fPhotons->At(fTriggeredPhotons[j]);
        const FourVectorM m2 = p2->Mom();
        h2TrigPhotonMass ->Fill((m1+m2).M());
      }
    }
  }
  else
    SkipEvent();

  return;
}

//--------------------------------------------------------------------------------------------------
void HttAnalysis::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches or objects created by earlier modules.

  // initialize some variables
  double pi = TMath::Pi();
  double mPi = -1.0 * TMath::Pi();

  // request the following Event Objects (often branches, but could be locally created ones)
  ReqEventObject(fPhotonName,fPhotons,fPhotonsFromBranch);

  // for MC only to adjust potential overlaps from generation
  if (! fIsData) {
    ReqBranch(Names::gkMCEvtInfoBrn,fMcEventInfo);
    printf(" Monte Carlo Information block %p\n",(void*) fMcEventInfo);
      printf(" --> this is no data. Access the McEventInfo.\n\n");
  }
  else
    printf(" --> this is data. Drop the McEventInfo.\n\n");

  // book all histograms
  AddTH1(hAllNPhotons,       "hAllNPhotons",";# of photons;#",                       21,-0.5, 20.5);
  AddTH1(hAllPhotonPt,       "hAllPhotonPt",";p_{T} [GeV];#",                       100, 0.0,250.0);
  AddTH1(hAllPhotonEta,      "hAllPhotonEta",";#eta;#",                             160,-8.0,  8.0);
  AddTH1(hAllPhotonSigEeBar, "hPhotonSigEeBar", ";#sigma_{#eta#eta}, barrel;#",     200, 0.0,  0.4);
  AddTH1(hAllPhotonSigEeEcp, "hPhotonSigEeEcp", ";#sigma_{#eta#eta}, endcap;#",     200, 0.0,  0.4);
  AddTH1(hAllPhotonEcHcTrIso,"hPhotonEcHcTrIso", ";#isolation_{ECAL,HCAL,Trk};#",   100, 0.0,  1.0);
  AddTH1(hAllPhotonPixelSeed,"hPhotonPixelSeed", ";#pixel seed;#",                    2,-0.5,  1.5);

  AddTH1(hTrigDeltaEta,      "hTrigDeltaEta",";#Delta #eta;#",                      100,-0.1,  0.1);
  AddTH1(hTrigDeltaPhi,      "hTrigDeltaPhi",";#Delta #phi;#",                      100,-0.05,0.05);
  AddTH1(hTrigDeltaEt,       "hTrigDeltaEt",";#Delta E_{T}/E_{T};#",                100,-0.1,  0.1);

  AddTH1(hNPhotons,          "hNPhotons",";# of sel photons;#",                      21,-0.5, 20.5);
  AddTH1(hPhotonEtMx,        "hPhotonEtMx", ";E_{T,max};#",                         200, 0.0,400.0);

  AddTH1(hN2Photons,         "hN2Photons",";# of di-photons;#",                      21,-0.5, 20.5);
  AddTH1(hPhotonEt1,         "hPhotonEt1", ";E_{T,1};#",                            200, 0.0,400.0);
  AddTH1(hPhotonEt2,         "hPhotonEt2", ";E_{T,2};#",                            200, 0.0,400.0);
  AddTH1(hPhotonEta1,        "hPhotonEta1", ";#eta_{1};#",                          104,-2.6,  2.6);
  AddTH1(hPhotonEta2,        "hPhotonEta2", ";#eta_{2};#",                          104,-2.6,  2.6);
  AddTH1(hPhotonPhi1,        "hPhotonPhi1", ";#phi_{1};#",                          100, mPi,   pi);
  AddTH1(hPhotonPhi2,        "hPhotonPhi2", ";#phi_{2};#",                          100, mPi,   pi);
  AddTH1(hPhotonR91,         "hPhotonR91", ";#R9_{1};#",                            100, 0.1,  1.1);
  AddTH1(hPhotonR92,         "hPhotonR92", ";#R9_{2};#",                            100, 0.1,  1.1);
  AddTH1(hPhotonDelR,        "hPhotonDelR", ";#Delta R_{1,2};#",                    104, 0.0,  5.2);
  AddTH1(h2PhotonEtSum,      "h2PhotonEtSum", ";E_{T,sum};#",                       400, 0.0,400.0);
  AddTH1(h2PhotonPt,         "h2PhotonPt", ";pair P_{T};#",                         400, 0.0,400.0);
  AddTH1(h2PhotonMass,       "h2PhotonMass",";Diphoton Mass [GeV/c^{2}];#",         400, 0.0,400.0);

  AddTH1(hHmPhotonR91,       "hHmPhotonR91", ";#R9_{1};#",                          100, 0.1,  1.1);
  AddTH1(hHmPhotonR92,       "hHmPhotonR92", ";#R9_{2};#",                          100, 0.1,  1.1);

  AddTH1(h2R9PhotonMass,     "h2R9PhotonMass", ";Diphoton Mass [GeV/c^{2}];#",      400, 0.0,400.0);

  AddTH1(hNTrigPhotons,      "hNTrigPhotons",";# of sel+trig photons;#",             21,-0.5, 20.5);
  AddTH1(h2TrigPhotonMass,   "h2TrigPhotonMass",";Diphoton Mass [GeV/c^{2}];#",     400, 0.0,400.0);

}

//--------------------------------------------------------------------------------------------------
void HttAnalysis::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void HttAnalysis::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void HttAnalysis::ResetSelectedPhotons()
{
  // Reset the selected and trigger matched photons

  fSelectedPhotons .clear();
  fTriggeredPhotons.clear();

  return;
}

//--------------------------------------------------------------------------------------------------
void HttAnalysis::MatchPhotonsToTrigger()
{
  // Reset the selected and trigger matched photons

  // access all trigger objects
  const TriggerObjectCol *tos = GetHLTObjects(fTrigObjsName);
  if (! tos)             // this can only happen if HLTMod::SetAbortIfNotAccepted(kFALSE) was called
    return;
  UInt_t nEnts = tos->GetEntries();

  // loop through our selected photons
  for (UInt_t i=0; i<fSelectedPhotons.size(); i++) {
    const Photon *p = fPhotons->At(fSelectedPhotons[i]);
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
      fTriggeredPhotons.push_back(fSelectedPhotons[i]);
    // print warning if desired
    if (matched > 1) {
      MDB(kGeneral,1) {
        printf(" Ambiguos match found: %d\n",matched);
        printf(" Photon:  %8.3f %8.3f %8.3f\n",p->Eta(), p->Phi(), p->Pt());
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
Bool_t HttAnalysis::PhotonTriggered(UInt_t iPhoton)
{
  // Determine whether the requested photon index is in the triggered photon indices

  // per default the photon did not trigger
  Bool_t triggered = kFALSE;

  // loop through triggered photon indices
  for (UInt_t i=0; i<fTriggeredPhotons.size(); i++) {
    if (fTriggeredPhotons[i] == iPhoton) {
      // index was found -> photon triggered
      triggered = kTRUE;
      break;
    }
  }

  return triggered;
}
