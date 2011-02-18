// $Id: $

#include <TMath.h>
#include <TH1D.h>
#include <TNtuple.h>
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitAna/DataUtil/interface/Debug.h"
#include "MitAna/DataTree/interface/Names.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/MCParticleCol.h"
#include "MitAna/DataTree/interface/MuonCol.h"
#include "MitAna/DataTree/interface/TrackCol.h"

#include "MitHtt/Mods/interface/ZmmAnalysis.h"

using namespace mithep;

ClassImp(mithep::ZmmAnalysis)

//--------------------------------------------------------------------------------------------------
ZmmAnalysis::ZmmAnalysis(const char *name, const char *title) :
  BaseMod             (name,title),
  // initialize bambu objects
  fTrigObjsName       (Names::gkHltObjBrn),
  fMuonName           (Names::gkMuonBrn),
  fMuonsFromBranch    (kTRUE),
  fMcEventInfoName    (Names::gkMCEvtInfoBrn),
  fMuons              (0),
  fMcEventInfo        (0),
  fIsData             (kTRUE),
  fSkimming           (kFALSE),
  // scale removal
  fOverlapCut         (-1.0),
  // cuts for selection
  fMuonEtMin          (20.0),
  //fNoPixelSeedCut     (kTRUE),
  //fMuonSigEeBarMax    ( 0.0105),
  //fMuonSigEeEcpMax    ( 0.0300),
  fMuonIsoMax         ( 0.30),
  // cuts for selection
  f2MuonMassMin       (140.0),
  //fMuonR9Min          (0.9),
  // initialize the histograms
  hAllNMuons          (0),
  hAllMuonPt          (0),
  hAllMuonEta         (0),
  //hAllMuonPixelSeed   (0),
  //hAllMuonSigEeBar    (0),
  //hAllMuonSigEeEcp    (0),
  hAllMuonEcHcTrIso   (0),
  hTrigDeltaEta       (0),
  hTrigDeltaPhi       (0),
  hTrigDeltaEt        (0),
  h2SeleMuonMass      (0),
  hNMuons             (0),
  hMuonEtMx           (0),
  hN2Muons            (0),
  hMuonEt1            (0),
  hMuonEt2            (0),
  hMuonEta1           (0),
  hMuonEta2           (0),
  hMuonPhi1           (0),
  hMuonPhi2           (0),
  //hMuonR91            (0),
  //hMuonR92            (0),
  hMuonDelR           (0),
  h2MuonEtSum         (0),
  h2MuonPt            (0),
  h2MuonMass          (0),
  //hHmMuonR91          (0),
  //hHmMuonR92          (0),
  //h2R9MuonMass        (0),
  hNTrigMuons         (0),
  h2TrigMuonMass      (0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void ZmmAnalysis::Begin()
{
  // Run startup code on the client machine. For this module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void ZmmAnalysis::Process()
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

  // access the muons
  LoadEventObject(fMuonName,fMuons);

  // monitor how many muons have been found in this event (no additional selection)
  hAllNMuons->Fill(fMuons->GetEntries());

  // loop through them and make some pictures and select the ones relevant to our analysis
  double lastEt = 99999.0;
  ResetSelectedMuons();
  for (UInt_t i=0; i<fMuons->GetEntries(); ++i) {
    const Muon *p = fMuons->At(i);
    hAllMuonPt ->Fill(p->Pt());
    hAllMuonEta->Fill(p->Eta());

    //double sigEtaEta = p->CoviEtaiEta();
    //if      (p->IsEB())
    //  hAllMuonSigEeBar->Fill(sigEtaEta);
    //else if (p->IsEE())
    //  hAllMuonSigEeEcp->Fill(sigEtaEta);
    double isoTot = (p->IsoR03SumPt() + p->IsoR03EmEt() + p->IsoR03HadEt())
                    / p->Pt();
    hAllMuonEcHcTrIso->Fill(isoTot);
    //hAllMuonPixelSeed->Fill((int) p->HasPixelSeed());

    // check the muons are really order according to pt
    if (lastEt <= p->Et()) {
      printf("\n WARNING -- muons are out of order (%f larger  %f)\n\n",p->Et(),lastEt);
    }

    // apply muon selection (just an Et cut)
    if (p->Et() > fMuonEtMin                   &&         // minimum Et
	//	(! p->HasPixelSeed() && fNoPixelSeedCut) &&         // no pixel seed
	isoTot  < fMuonIsoMax)                 //&&         // maximum isolation
	//      ((p->IsEB() && sigEtaEta < fMuonSigEeBarMax) ||   // shower shape (bar/ecp separate)
	//      (p->IsEE() && sigEtaEta < fMuonSigEeEcpMax) 	 ))
      fSelectedMuons.push_back(i);

    // make sure to keep the Et accounting up to speed
    lastEt = p->Et();
  }
  hNMuons->Fill(fSelectedMuons.size());

  MatchMuonsToTrigger();
  hNTrigMuons->Fill(fTriggeredMuons.size());

  // loop through our selected muons and make histograms only if one was matched to trigger
  UInt_t n2Muons = 0;
  for (UInt_t i=0; i<fSelectedMuons.size(); i++) {
    const Muon *p1 = fMuons->At(fSelectedMuons[i]);

    // histogram maximum muon Et
    if (i==0)
      hMuonEtMx->Fill(p1->Et());

    // is this a triggered muon?
    Bool_t isP1Trig = MuonTriggered(fSelectedMuons[i]);

    // loop through all pairs
    const FourVectorM m1 = p1->Mom();
    for (UInt_t j=i+1; j<fSelectedMuons.size(); ++j) {
      const Muon *p2 = fMuons->At(fSelectedMuons[j]);

      // is this a triggered muon?
      Bool_t isP2Trig = MuonTriggered(fSelectedMuons[j]);

	const FourVectorM m2 = p2->Mom();
	double mass = (m1+m2).M();

	h2SeleMuonMass->Fill(mass);

      // only consider pairs, where at least one muon triggered
      if (isP1Trig || isP2Trig) {
	

	// write in ntuple
	int lEvt   = GetEventHeader()->EvtNum();
	int lLumi  = GetEventHeader()->LumiSec();
	int lRun   = GetEventHeader()->RunNum();
	if (1)
	  {
	    Float_t vals[17] = {lEvt, lLumi, lRun, p1->Pt(), p2->Pt(), 
				p1->Eta(),p2->Eta(), p1->Phi(), p2->Phi(),
				p1->D0PV(), p2->D0PV(), p1->Charge(), p2->Charge(), 
				mass, 
				0, 0, 0};
	    fNt->Fill(vals);
	  }
	
        hMuonEt1   ->Fill(p1->Et());
        hMuonEt2   ->Fill(p2->Et());
        hMuonEta1  ->Fill(p1->Eta());
        hMuonEta2  ->Fill(p2->Eta());
        hMuonPhi1  ->Fill(p1->Phi());
        hMuonPhi2  ->Fill(p2->Phi());
        //hMuonR91   ->Fill(p1->R9());
        //hMuonR92   ->Fill(p2->R9());
        hMuonDelR  ->Fill(MathUtils::DeltaR(m1,m2));
        h2MuonEtSum->Fill(p1->Et()+p2->Et());
        h2MuonPt   ->Fill((m1+m2).Pt());
        h2MuonMass ->Fill(mass);
	  //// our high mass selection
	  //if (mass > f2MuonMassMin) {
	  //  hHmMuonR91->Fill(p1->R9());
	  //  hHmMuonR92->Fill(p2->R9());
	  //}
	  //// our R9 cleaned muon selection
	  //if (p1->R9() > fMuonR9Min && p2->R9() > fMuonR9Min) {
	  //  h2R9MuonMass->Fill(mass);
	  //}
        // increment the number of good di-muon candidates
        n2Muons++;
      }
    }
  }
  hN2Muons->Fill(n2Muons);

  // work only with triggered muons
  for (UInt_t i=0; i<fTriggeredMuons.size(); i++) {
    const Muon *p1 = fMuons->At(fTriggeredMuons[i]);

    // loop through all pairs
    const FourVectorM m1 = p1->Mom();
    for (UInt_t j=i+1; j<fTriggeredMuons.size(); ++j) {
      const Muon *p2 = fMuons->At(fTriggeredMuons[j]);
      const FourVectorM m2 = p2->Mom();
      h2TrigMuonMass ->Fill((m1+m2).M());
    }
  }

  // if (fTriggeredMuons.size()<=0) {
  // if (fSkimming)
  // SkipEvent();
  // }

  return;
}

//--------------------------------------------------------------------------------------------------
void ZmmAnalysis::SlaveBegin()
{
  // Run startup code on the computer (slave) doing the actual analysis. Here,
  // we typically initialize histograms and other analysis objects and request
  // branches or objects created by earlier modules.
  if (1)
    {
      string s;
      s = "run:evt:lumi:pt1:pt2:eta1:eta2:phi1:phi2:dca1:dca2:ch1:ch2:mass:met:trigm:trigm";
      
      fNt = new TNtuple("ntZmumu","ZmmmAnalysisMod",s.c_str());
      
      AddOutput(fNt);
    }
  
  // initialize some variables
  double pi = TMath::Pi();
  double mPi = -1.0 * TMath::Pi();

  // request the following Event Objects (often branches, but could be locally created ones)
  ReqEventObject(fMuonName,fMuons,fMuonsFromBranch);

  // for MC only to adjust potential overlaps from generation
  if (! fIsData) {
    ReqBranch(Names::gkMCEvtInfoBrn,fMcEventInfo);
    printf(" Monte Carlo Information block %p\n",(void*) fMcEventInfo);
      printf(" --> this is no data. Access the McEventInfo.\n\n");
  }
  else
    printf(" --> this is data. Drop the McEventInfo.\n\n");

  // book all histograms
  AddTH1(hAllNMuons,       "hAllNMuons",";# of muons;#",                         21,-0.5, 20.5);
  AddTH1(hAllMuonPt,       "hAllMuonPt",";p_{T} [GeV];#",                       100, 0.0,250.0);
  AddTH1(hAllMuonEta,      "hAllMuonEta",";#eta;#",                             160,-8.0,  8.0);
  //AddTH1(hAllMuonSigEeBar, "hMuonSigEeBar", ";#sigma_{#eta#eta}, barrel;#",     200, 0.0,  0.4);
  //AddTH1(hAllMuonSigEeEcp, "hMuonSigEeEcp", ";#sigma_{#eta#eta}, endcap;#",     200, 0.0,  0.4);
  AddTH1(hAllMuonEcHcTrIso,"hMuonEcHcTrIso", ";#isolation_{ECAL,HCAL,Trk};#",   100, 0.0,  1.0);
  //AddTH1(hAllMuonPixelSeed,"hMuonPixelSeed", ";#pixel seed;#",                    2,-0.5,  1.5);

  AddTH1(hTrigDeltaEta,    "hTrigDeltaEta",";#Delta #eta;#",                    100,-0.1,  0.1);
  AddTH1(hTrigDeltaPhi,    "hTrigDeltaPhi",";#Delta #phi;#",                    100,-0.05,0.05);
  AddTH1(hTrigDeltaEt,     "hTrigDeltaEt",";#Delta E_{T}/E_{T};#",              100,-0.1,  0.1);

  AddTH1(h2SeleMuonMass,   "h2SeleMuonMass",";Dimuon Mass [GeV/c^{2}];#",       400, 0.0,400.0);

  AddTH1(hNMuons,          "hNMuons",";# of sel muons;#",                        21,-0.5, 20.5);
  AddTH1(hMuonEtMx,        "hMuonEtMx", ";E_{T,max};#",                         200, 0.0,400.0);

  AddTH1(hN2Muons,         "hN2Muons",";# of di-muons;#",                        21,-0.5, 20.5);
  AddTH1(hMuonEt1,         "hMuonEt1", ";E_{T,1};#",                            200, 0.0,400.0);
  AddTH1(hMuonEt2,         "hMuonEt2", ";E_{T,2};#",                            200, 0.0,400.0);
  AddTH1(hMuonEta1,        "hMuonEta1", ";#eta_{1};#",                          104,-2.6,  2.6);
  AddTH1(hMuonEta2,        "hMuonEta2", ";#eta_{2};#",                          104,-2.6,  2.6);
  AddTH1(hMuonPhi1,        "hMuonPhi1", ";#phi_{1};#",                          100, mPi,   pi);
  AddTH1(hMuonPhi2,        "hMuonPhi2", ";#phi_{2};#",                          100, mPi,   pi);
  //AddTH1(hMuonR91,         "hMuonR91", ";#R9_{1};#",                            100, 0.1,  1.1);
  //AddTH1(hMuonR92,         "hMuonR92", ";#R9_{2};#",                            100, 0.1,  1.1);
  AddTH1(hMuonDelR,        "hMuonDelR", ";#Delta R_{1,2};#",                    104, 0.0,  5.2);
  AddTH1(h2MuonEtSum,      "h2MuonEtSum", ";E_{T,sum};#",                       400, 0.0,400.0);
  AddTH1(h2MuonPt,         "h2MuonPt", ";pair P_{T};#",                         400, 0.0,400.0);
  AddTH1(h2MuonMass,       "h2MuonMass",";Dimuon Mass [GeV/c^{2}];#",           400, 0.0,400.0);

  //AddTH1(hHmMuonR91,       "hHmMuonR91", ";#R9_{1};#",                          100, 0.1,  1.1);
  //AddTH1(hHmMuonR92,       "hHmMuonR92", ";#R9_{2};#",                          100, 0.1,  1.1);

  //AddTH1(h2R9MuonMass,     "h2R9MuonMass", ";Dimuon Mass [GeV/c^{2}];#",        400, 0.0,400.0);

  AddTH1(hNTrigMuons,      "hNTrigMuons",";# of sel+trig muons;#",               21,-0.5, 20.5);
  AddTH1(h2TrigMuonMass,   "h2TrigMuonMass",";Dimuon Mass [GeV/c^{2}];#",       400, 0.0,400.0);

}

//--------------------------------------------------------------------------------------------------
void ZmmAnalysis::SlaveTerminate()
{
  // Run finishing code on the computer (slave) that did the analysis. For this
  // module, we dont do anything here.
}

//--------------------------------------------------------------------------------------------------
void ZmmAnalysis::Terminate()
{
  // Run finishing code on the client computer. For this module, we dont do
  // anything here.
}

//--------------------------------------------------------------------------------------------------
void ZmmAnalysis::ResetSelectedMuons()
{
  // Reset the seleted and trigger matched muons

  fSelectedMuons .clear();
  fTriggeredMuons.clear();

  return;
}

//--------------------------------------------------------------------------------------------------
void ZmmAnalysis::MatchMuonsToTrigger()
{
  // Reset the selected and trigger matched muons

  // access all trigger objects
  const TriggerObjectCol *tos = GetHLTObjects(fTrigObjsName);
  if (! tos)             // this can only happen if HLTMod::SetAbortIfNotAccepted(kFALSE) was called
    return;
  UInt_t nEnts = tos->GetEntries();

  // loop through our selected muons
  for (UInt_t i=0; i<fSelectedMuons.size(); i++) {
    const Muon *p = fMuons->At(fSelectedMuons[i]);
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
      fTriggeredMuons.push_back(fSelectedMuons[i]);
    // print warning if desired
    if (matched > 1) {
      MDB(kGeneral,2) {
        printf(" Ambiguos match found: %d\n",matched);
        printf(" Muon:  %8.3f %8.3f %8.3f\n",p->Eta(), p->Phi(), p->Pt());
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
Bool_t ZmmAnalysis::MuonTriggered(UInt_t iMuon)
{
  // Determine whether the requested muon index is in the triggered muon indices

  // per default the muon did not trigger
  Bool_t triggered = kFALSE;

  // loop through triggered muon indices
  for (UInt_t i=0; i<fTriggeredMuons.size(); i++) {
    if (fTriggeredMuons[i] == iMuon) {
      // index was found -> muon triggered
      triggered = kTRUE;
      break;
    }
  }

  return triggered;
}
