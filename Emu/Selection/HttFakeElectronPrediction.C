//root -l MitHtt/Hww/FakeLeptonBkg/HwwFakeElectronPrediction.C+\(\"HwwNtuple_r11a-dmu-pr-v1_noskim_0000.root\",130,\"\"\) |& tee debugLog
//root -l MitHtt/Hww/FakeLeptonBkg/HwwFakeElectronPrediction.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-full-pr_TightPlusRecoNoTriggerSkim.root\",130,\"\"\)


//root -l MitHtt/Hww/FakeLeptonBkg/HwwFakeElectronPrediction.C+\(\"/home/sixie/hist/HwwAnalysis/2011Data/WWAnalysisSkimmed_r11a-full-pr_TightPlusRecoTriggerSkim.root\",130,\"\"\)
//root -l MitHtt/Hww/FakeRate/HwwFakeElectronPrediction.C+\(\"WWAnalysisSkimmed_r10b-mu-d22_TwoRecoLeptonNoTriggerSkim.root\",130,\"\"\)
//root -l MitHtt/Hww/FakeRate/HwwFakeElectronPrediction.C+\(\"WWAnalysisSkimmed_full-d22_TwoRecoLeptonWithTriggerPlusOneTightLeptonSkim.root\",130,\"\"\)
//================================================================================================
//
// Htt selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
// Debug:
// cat diff | grep ">" | awk '{print "|| (info->runNum == " $2 " && info->evtNum == " $4 ")"}'
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TMath.h>                  // 
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "EmuData.hh"
#include "Common/MyTools.hh"        // miscellaneous helper functions                                                                            

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Utils/LeptonIDCuts.hh"
#include "MitAna/DataTree/interface/BaseVertex.h"

// lumi section selection with JSON files
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitPhysics/FakeMods/interface/FakeRate.h"

#endif

using namespace std;
using namespace mithep;


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passElectronCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passOldElectronDenominatorCuts(const mithep::TElectron *ele);
Bool_t passMuonCuts(const mithep::TMuon *mu);
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HttNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}



//=== MAIN MACRO =================================================================================================

void HttFakeElectronPrediction()
{  
  gBenchmark->Start("HttFakeElectronPrediction");


  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 20;

  //--------------------------------------------------------------------------------------------------------------
  // Set up Fake Rate
  //==============================================================================================================
  Bool_t use2DFakeRate = kTRUE;
  Bool_t useFitFunction = kFALSE;
  TString fakeinfile("data/ElectronFakeRate.root");
  ifstream fakechk; fakechk.open(fakeinfile.Data()); assert(fakechk.is_open()); fakechk.close();
  mithep::FakeRate *fFakeRate = new mithep::FakeRate(fakeinfile.Data(),
						     fakeinfile.Data(),
						     "", "",
						     "ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta",
						     "ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta",
						     use2DFakeRate, useFitFunction );

  vector<TString> fnamev;
  fnamev.push_back("/tmp/dkralph/htt/r11a-mueg-pr-v1_ntuple.root");
  fnamev.push_back("/tmp/dkralph/htt/r11a-mueg-pr-v2_ntuple.root");

  // Data structures to store info from TTrees
  mithep::TEventInfo	*info		=	new mithep::TEventInfo();
  TClonesArray		*electronArr	=	new TClonesArray("mithep::TElectron");
  TClonesArray		*muonArr	=	new TClonesArray("mithep::TMuon");
  TClonesArray		*jetArr		=	new TClonesArray("mithep::TJet");
  TClonesArray		*pvArr		=	new TClonesArray("mithep::TVertex");
  

  // make output ntuple file
  TFile fakefile("fakes_select.root","recreate");
  TTree outtree("Events","Events");

  EmuData data;
  Double_t trigeff=1, rawMet, rawprojvar;
  UInt_t npt20jets;
  const UInt_t kMaxPt20Jets=15;
  TArrayF btagArray; btagArray.Set(kMaxPt20Jets); // array to hold b-tag values for pt-20 jets
  outtree.Branch("Events",&data.runNum,
"runNum/i:evtNum:lumiSec:nPV:njets:nbjets:met/F:metphi:mass:dphi:mt:pt:phi:pmet:pvis:lpt1:leta1:lphi1:lpt2:leta2:lphi2:jpt1:jeta1:jphi1:jpt2:jeta2:jphi2:bjpt:bjeta:bjphi:mjj:weight:state/I");
  // extra branches
  outtree.Branch("npt20jets",&npt20jets);
  outtree.Branch("btagArray",&btagArray);
  outtree.Branch("trigeff",&trigeff);
  outtree.Branch("rawMet",&rawMet);
  outtree.Branch("rawprojvar",&rawprojvar);

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    //
    // Access samples and fill histograms
    TFile *inputFile=0;
    TTree *eventTree=0;  
   
    inputFile = new TFile(fnamev[ifile]);
    assert(inputFile);

    //********************************************************
    // Good RunLumi Selection
    //********************************************************
    Bool_t hasJSON = kTRUE;
    mithep::RunLumiRangeMap rlrm;
    TString jsonfname("/home/$USER/json/Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON.txt");
    ifstream jsonchk; jsonchk.open(jsonfname); assert(jsonchk.is_open()); jsonchk.close();
    rlrm.AddJSONFile(jsonfname.Data()); 

    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(fnamev[ifile],"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *jetBr;
    TBranch *pvBr;

    //*****************************************************************************************
    //Loop over muon Data Tree
    //*****************************************************************************************
    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",		&info);		infoBr	    = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron",	&electronArr);	electronBr  = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("Muon",		&muonArr);      muonBr	    = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("PFJet",		&jetArr);       jetBr	    = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("PV",		&pvArr);        pvBr	    = eventTree->GetBranch("PV");
  
    vector<UInt_t> runs;
    vector<UInt_t> events;

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
      UInt_t trigger = kHLT_Mu17_Ele8_CaloIdL | kHLT_Mu8_Ele17_CaloIdL;
      if(!(info->triggerBits & trigger)) continue;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      jetArr->Clear(); 
      pvArr->Clear();
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);
      pvBr->GetEntry(ientry);

      Int_t nGoodPV = 0;
      for(Int_t i=0; i<pvArr->GetEntries(); i++) {
	const mithep::TVertex *pv = (mithep::TVertex*)((*pvArr)[i]);
	if(pv->ndof < 4 && fabs(pv->z) > 24.0) continue;
	mithep::BaseVertex *PV = new mithep::BaseVertex(pv->x,pv->y,pv->z);
	if (PV->Position().Rho() > 2) continue;
	nGoodPV++;
      }

      // add muons and electrons to vectors of good lepton properties
      vector<Int_t> leptonType;   // leptons passing lepton selection
      vector<Int_t> leptonIndex;
      vector<Double_t> leptonPt;
      vector<Double_t> leptonEta;
      vector<Double_t> leptonPhi;
      vector<Int_t> leptonCharge;

      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
	const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);      
	if ( (0==0)
	     &&
	     passMuonCuts(mu)
	     &&
	     fabs(mu->eta) < 2.1
	     && 
	     mu->pt > 15.0
	     ) {
	  if((info->triggerBits & kHLT_Mu17_Ele8_CaloIdL) && mu->pt < 20.0 ) continue;
	  leptonPt.push_back(mu->pt);
	  leptonEta.push_back(mu->eta);
	  leptonPhi.push_back(mu->phi);
	  leptonType.push_back(13);
	  leptonIndex.push_back(i);  
	  leptonCharge.push_back(mu->q);
	}
      }
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
	const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
	Bool_t isMuonOverlap = kFALSE;
	for (UInt_t k=0; k<leptonPt.size(); ++k) {
	  if ( leptonType[k] == 13 
	       && mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[k],leptonEta[k]) < 0.1
	       ) {
	    isMuonOverlap = kTRUE; 
	    break;
	  }        
	}
      
	if ( (0==0)
	     &&
	     passEleID(ele) // this function is coming from LeptonIDCuts.hh
	     // passElectronCuts(ele)
	     &&
	     fabs(ele->eta) < 2.5
	     && 
	     ele->pt > 15.0
	     &&
	     !isMuonOverlap
	     ) {
	  leptonPt.push_back(ele->pt);
	  leptonEta.push_back(ele->eta);
	  leptonPhi.push_back(ele->phi);
	  leptonType.push_back(11);
	  leptonIndex.push_back(i);
	  leptonCharge.push_back(ele->q);
	}
      }

      // sort leptons
      Int_t tempType;
      Int_t tempIndex;
      Double_t tempPt;
      Double_t tempEta;
      Double_t tempPhi;
      Int_t tempCharge;
      for (UInt_t l=0; l<leptonIndex.size(); l++) {
	for (UInt_t k=0; k < leptonIndex.size() - 1; k++) {
	  if (leptonPt[k+1] > leptonPt[k]) {
	    tempType = leptonType[k];
	    tempIndex = leptonIndex[k];
	    tempPt = leptonPt[k];
	    tempEta = leptonEta[k];
	    tempPhi = leptonPhi[k];
	    tempCharge = leptonCharge[k];
          
	    leptonType[k] = leptonType[k+1];
	    leptonIndex[k] = leptonIndex[k+1];
	    leptonPt[k] = leptonPt[k+1];
	    leptonEta[k] = leptonEta[k+1];
	    leptonPhi[k] = leptonPhi[k+1];
	    leptonCharge[k] = leptonCharge[k+1];

	    leptonType[k+1] = tempType;
	    leptonIndex[k+1] = tempIndex;
	    leptonPt[k+1] = tempPt;
	    leptonEta[k+1] = tempEta;
	    leptonPhi[k+1] = tempPhi;
	    leptonCharge[k+1] = tempCharge;
          
	  }
	}
      }
    
      //******************************************************************************
      // select events with 1 lepton only
      //******************************************************************************
      if (leptonPt.size() < 1) continue;
      if (!(leptonPt[0] > 15.0)) continue;
      if (!(leptonPt.size() < 2 || leptonPt[1] < 15.0)) continue;
 
      //******************************************************************************
      // Check duplicate events
      //******************************************************************************
      Bool_t foundRunAndEvent = kFALSE;
      for (UInt_t q=0; q<runs.size(); ++q) {
	if (runs[q] == info->runNum && events[q] == info->evtNum) foundRunAndEvent = kTRUE;
      }
      if (!foundRunAndEvent) {
	runs.push_back(info->runNum);
	events.push_back(info->evtNum);
      } else {
	continue;
      }

      // loop over electron denominators
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
	const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
   
	if(leptonCharge[0] == ele->q)								continue;
	if (!(ele->pt > 15.0 && fabs(ele->eta) < 2.5))						continue;
	if ( mithep::MathUtils::DeltaR(ele->phi, ele->eta, leptonPhi[0],leptonEta[0]) < 0.1 )	continue;
	if((info->triggerBits & kHLT_Mu8_Ele17_CaloIdL) && ele->pt < 20.0 )			continue;

	// select denominators that fail final selection
	if (!passElectronDenominatorCuts(ele))	continue;
	// if (passElectronCuts(ele))		continue;
	if (passEleID(ele))	        	continue;

	//******************************************************************************
	//Calculate Fake Rate
	//******************************************************************************
	double fakeRate          = fFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi)              /     (1 - fFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi));
	double fakeRateErrorLow  = fFakeRate->ElectronFakeRateStatErrorLow(ele->pt, ele->eta, ele->phi)  / pow((1 - fFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi)),2);
	double fakeRateErrorHigh = fFakeRate->ElectronFakeRateStatErrorHigh(ele->pt, ele->eta, ele->phi) / pow((1 - fFakeRate->ElectronFakeRate(ele->pt, ele->eta, ele->phi)),2);

	Double_t eventWeight = fakeRate; 
	Double_t eventWeightError = (  fakeRateErrorLow +   fakeRateErrorHigh) / 2;


	enum { kMuMu, kEleEle, kEleMu, kMuEle };      // final state type
	Int_t finalState = -1;
	TLorentzVector lep1, lep2, dilep;  // lepton 4-vectors
	if(leptonType[0] != 13) continue;
	const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[leptonIndex[0]]);

	if(leptonPt[0] > ele->pt) {
	  finalState = kMuEle;
	  lep1.SetPtEtaPhiM(leptonPt[0], leptonEta[0], leptonPhi[0], 105.658369e-3 );
	  lep2.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
	} else {
	  finalState = kEleMu;
	  lep2.SetPtEtaPhiM(leptonPt[0], leptonEta[0], leptonPhi[0], 105.658369e-3 );
	  lep1.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
	}
	dilep = lep1+lep2;

	if(lep1.Pt() < 20 || lep2.Pt() < 15) continue;

	// loop through jets      
	UInt_t njets = 0, nbjets = 0;
	const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	btagArray.Reset(); npt20jets=0;
	for(Int_t ijet=0; ijet<jetArr->GetEntriesFast(); ijet++) {
	  const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[ijet]);

	  if(toolbox::deltaR(jet->eta,jet->phi,lep1.Eta(),lep1.Phi()) < 0.3) continue;
	  if(toolbox::deltaR(jet->eta,jet->phi,lep2.Eta(),lep2.Phi()) < 0.3) continue;

	  if(fabs(jet->eta) > 5) continue;

	  // look for b-jets
	  if((jet->pt > kBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
	    assert(npt20jets<kMaxPt20Jets);
	    btagArray.AddAt(jet->tche,npt20jets); npt20jets++;
	    if(jet->tche > 3.3) {
	      nbjets++;
	      if(!bjet || jet->pt > bjet->pt)
		bjet = jet; // leading b-jet
	    }
	  }

	  // look for vbf jets
	  if(jet->pt > kJetPtMin) {
	    njets++;
	    if(!jet1 || jet->pt > jet1->pt) { // jet1 is highest pt, jet2 next highest
	      jet2 = jet1;
	      jet1 = jet;
	    } else if(!jet2 || jet->pt > jet2->pt) {
	      jet2 = jet;
	    }
	  }		    
	}

	TLorentzVector jv1, jv2, dijet;
	if(njets>1) {
	  jv1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
	  jv2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
	  dijet = jv1+jv2;
	} 

	// calculate projection variables
	TVector3 m,e,metv;
	m.SetPtEtaPhi(mu->pt,0,mu->phi);
	e.SetPtEtaPhi(ele->pt,0,ele->phi);
	metv.SetPtEtaPhi(info->pfMET,0,info->pfMETphi); // uncorrected met
	TVector3 bisector(m.Unit() + e.Unit());
	bisector = bisector.Unit();
	Double_t projVis  = (m+e).Dot(bisector);
	Double_t projMet  =   metv.Dot(bisector);
	rawMet = info->pfMET;
	rawprojvar  = 0.85*projVis - projMet;

	data.runNum  = info->runNum;
	data.evtNum  = info->evtNum;
	data.lumiSec = info->lumiSec;
	data.nPV     = nGoodPV;
	data.njets   = njets;
	data.nbjets  = nbjets;
	data.met     = info->pfMET;
	data.metphi  = info->pfMETphi;
	data.mass    = dilep.M();
	data.dphi    = toolbox::deltaPhi(ele->phi,mu->phi);
	data.mt      = sqrt( 2.0 * (dilep.Pt()) * info->pfMET * (1.0-cos(toolbox::deltaPhi(dilep.Phi(),info->pfMETphi))) );
	data.pt      = dilep.Pt();
	data.phi     = dilep.Phi();
	data.pmet    = projMet;
	data.pvis    = projVis;
	data.lpt1    = lep1.Pt();
	data.leta1   = lep1.Eta();
	data.lphi1   = lep1.Phi();
	data.lpt2    = lep2.Pt();
	data.leta2   = lep2.Eta();
	data.lphi2   = lep2.Phi();
	data.jpt1    = (jet1) ? jet1->pt  : 0;
	data.jeta1   = (jet1) ? jet1->eta : 0;
	data.jphi1   = (jet1) ? jet1->phi : 0;
	data.jpt2    = (jet2) ? jet2->pt  : 0;
	data.jeta2   = (jet2) ? jet2->eta : 0;
	data.jphi2   = (jet2) ? jet2->phi : 0;
	data.bjpt    = (bjet) ? bjet->pt  : 0;
	data.bjeta   = (bjet) ? bjet->eta : 0;
	data.bjphi   = (bjet) ? bjet->phi : 0;
	data.mjj     = (njets>1) ? dijet.M() : 0;
	data.weight  = eventWeight;//(isam==0) ? 1 : weight*kf*npvWgt*trigeff/lumi;
	data.state   = finalState;

	outtree.Fill();

      } // end loop over electron denominators

    } // end event loop

    delete inputFile;
    inputFile = 0;
  }

  fakefile.Print();
  fakefile.Write();
  fakefile.Close();

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;
  delete pvArr;

  gBenchmark->Show("HttFakeElectronPrediction");       
} 


Bool_t passElectronCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  if (fabs(ele->eta) >= 2.5) pass = kFALSE;

  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            //&& (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
            && ele->pfIso04 / ele->pt < 0.13
            && ele->nExpHitsInner <= 0
            && !ele->isConv
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.2
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
             && ele->HoverE < 0.10
             //&& (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
             && ele->pfIso04 / ele->pt < 0.09
             && ele->nExpHitsInner <= 0
             && !ele->isConv
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.2
          )
      ) {
      pass = kFALSE;
    }
  } 

  if (ele->pt < 20) {
    //Barrel 
    if (fabs(ele->scEta < 1.479)) {
      if (! ( (0==0)
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.03
              && ele->HoverE < 0.025
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else  {
      if (! (  (0==0)
               && fabs(ele->deltaEtaIn) < 0.005
               && fabs(ele->deltaPhiIn) < 0.02
            )
        ) {
        pass = kFALSE;
      }
    } 

    if (ele->fBrem <= 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EoverP > 0.95 )) pass = kFALSE;
      }
    }
  }



  return pass;
}


Bool_t passMuonCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (mu->pt < 10) pass = kFALSE;
  if (fabs(mu->eta) > 2.1) pass = kFALSE;

  Double_t isoCutValue = 0.0;

  if (fabs(mu->eta) < 1.479) {
    if (mu->pt > 20) {
      isoCutValue = 0.13;
    } else {
      isoCutValue = 0.06;
    }
  } else {
    if (mu->pt > 20) {
      isoCutValue = 0.09;
    } else {
      isoCutValue = 0.05;
    }
  }


  if (! 
      ( mu->typeBits & kGlobal
        && mu->typeBits & kTracker
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.02
        && fabs(mu->dz) < 0.2
        //&& (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 0.15
        && mu->pfIso03 / mu->pt < isoCutValue
        && (mu->nMatch > 1 )
        && (mu->nPixHits > 0)
        && (mu->ptErr / mu->pt < 0.1)
        )
    ) pass = kFALSE;

  return pass;
}

Bool_t passElectronDenominatorCuts(const mithep::TElectron *ele) {

  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 15.0 && fabs(ele->eta) < 2.5)) pass = kFALSE;

  //Barrel

  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && ele->nExpHitsInner <= 0
            && !ele->isConv
            && fabs(ele->dz) < 0.1
            && (ele->trkIso03 ) / ele->pt < 0.20
            && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
            && (ele->hadIso03) / ele->pt < 0.20

	    )
	) {
      pass = kFALSE;
    }
  }
  //Endcap    
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.009
             && fabs(ele->deltaPhiIn) < 0.10
             && ele->HoverE < 0.10
             && ele->nExpHitsInner <= 0
             && !ele->isConv
             && fabs(ele->dz) < 0.1
             && (ele->trkIso03) / ele->pt < 0.20
             && (ele->emIso03) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20
	     )
	) {
      pass = kFALSE;
    }
  }

  return pass;
}

Bool_t passOldElectronDenominatorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) < 2.5)) pass = kFALSE;


  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.007
            && fabs(ele->deltaPhiIn) < 0.15
            && ele->HoverE < 0.12
            && ele->nExpHitsInner <= 0
            && !ele->isConv
            && fabs(ele->dz) < 0.1
            && (ele->trkIso03 ) / ele->pt < 0.20
            && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
            && (ele->hadIso03) / ele->pt < 0.20
              
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.009
             && fabs(ele->deltaPhiIn) < 0.10
	     && ele->HoverE < 0.10
             && ele->nExpHitsInner <= 0
             && !ele->isConv
             && fabs(ele->dz) < 0.1
             && (ele->trkIso03) / ele->pt < 0.20
             && (ele->emIso03) / ele->pt < 0.20
             && (ele->hadIso03) / ele->pt < 0.20
          )
      ) {
      pass = kFALSE;
    }
  } 

  return pass;
}



Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu) {
  
  Bool_t pass = kTRUE;

  if (! 
      ( mu->typeBits & kGlobal
        && mu->nTkHits > 10
        && mu->muNchi2 < 10.0
        && (mu->qualityBits & kGlobalMuonPromptTight)
        && fabs(mu->d0) < 0.2
        && fabs(mu->dz) < 0.1
        && (mu->trkIso03 + mu->emIso03 + mu->hadIso03) / mu->pt < 1.0
        && (mu->ptErr / mu->pt < 0.1)
        )
    ) pass = kFALSE;

  return pass;
}

//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2)
{
  ofs << "Run:" << runNum;
  ofs << "  Lumi:" << lumiSec;
  ofs << "  Event:" << evtNum;
  ofs << "  mass: " << mass;
  
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
//   ofs << "    pt    |    eta    |    phi    |   iso    |    d0      | ntk | npx | nseg | nval | chi^2/ndf | TM | HLT" << endl;
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
  ofs << " " ;
  ofs << setw(9) << pt1 << " |";
  ofs << setw(10) << eta1 << " |";
  ofs << setw(10) << phi1 << " |";
  ofs << setw(10) << leptonCharge1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << setw(10) << leptonCharge2 << " |";
  ofs << endl;
  
}
