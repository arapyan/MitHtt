#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                   // access to gROOT, entry point to ROOT system
#include <TSystem.h>                 // interface to OS
#include <TFile.h>                   // file handle class
#include <TTree.h>                   // class to access ntuples
#include <TClonesArray.h>            // ROOT array class
#include <TBenchmark.h>              // class to track macro running statistics
#include <TLorentzVector.h>          // class for 4-vector calculations
#include <TH2D.h>                    // 2D histogram class
#include <vector>                    // vector class
#include <iostream>                  // standard I/O
#include <iomanip>                   // functions to format standard I/O
#include <fstream>                   // functions for file I/O

#include "MitHtt/Common/MyTools.hh"         // custom helper functions
#include "MitHtt/Common/CEffUser2D.hh"

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Ntupler/interface/TSVfit.h"
#include "MitHtt/Ntupler/interface/TSVfitter.h"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// lepton ID helper functions
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitPhysics/FakeMods/interface/FakeRate.h"
#include "MitHtt/Utils/LeptonIDCuts.hh" 

#include "MitHtt/Emu/Selection/EmuData.hh"

#endif


//=== MAIN MACRO =================================================================================================

void predictFakesMuEle(const TString outputDir   // output directory
) {
  gBenchmark->Start("predictFakesMuEle");


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };
  
  vector<TString> fnamev;
  vector<TString> datatypev;
  vector<TString> jsonv;

  fnamev.push_back("/data/blue/vdutta/htt/025/mva/new/r11a-mueg-m10-v1_emu_skim.root"); datatypev.push_back(eMuEl);
  jsonv.push_back("/home/vdutta/cms/root/json/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt");

  fnamev.push_back("/data/blue/vdutta/htt/025/mva/new/r11a-mueg-pr-v4_emu_skim.root");  datatypev.push_back(eMuEl);
  jsonv.push_back("/home/vdutta/cms/root/json/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");

  fnamev.push_back("/data/blue/vdutta/htt/025/mva/new/r11a-mueg-a05-v1_emu_skim.root");  datatypev.push_back(eMuEl);
  jsonv.push_back("/home/vdutta/cms/root/json/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v3.txt");

  fnamev.push_back("/data/blue/vdutta/htt/025/mva/new/r11a-mueg-o03-v1_emu_skim.root");  datatypev.push_back(eMuEl);
  jsonv.push_back("/home/vdutta/cms/root/json/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");

  fnamev.push_back("/data/blue/vdutta/htt/025/mva/new/r11b-mueg-pr-v1_emu_skim.root");  datatypev.push_back(eMuEl);
  jsonv.push_back("/home/vdutta/cms/root/json/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");

  const Double_t kMuonPt1Min = 20;
  const Double_t kMuonPt2Min = 10;
  
  const Double_t kElePt1Min  = 20;
  const Double_t kElePt2Min  = 10;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 15;
  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };      // final state type

  // setup MVA
  mithep::ElectronIDMVA *electronIDMVANoIPInfo = new mithep::ElectronIDMVA();
  electronIDMVANoIPInfo->Initialize("BDTG method",
                              TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")), 
                              TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                              TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                              TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                              TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                              TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                              mithep::ElectronIDMVA::kNoIPInfo );

    
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  //
  // parse fake rate data
  //

  Bool_t use2DFakeRate = kTRUE;
  Bool_t useFitFunction = kFALSE;
  mithep::FakeRate *fFakeRate = new mithep::FakeRate("ElectronFakeRate.root",
			     "MuonFakeRate.root",
			     "", "",
			     "ElectronFakeRateDenominatorV4_Ele8CaloIdLCaloIsoVLCombinedSample_ptThreshold35_PtEta",
			     "MuonFakeRateDenominatorV2_Mu8PtCombinedSample_ptThreshold15_PtEta",
			     use2DFakeRate, useFitFunction );

  
  // cout << fr.getEff(0.3,60) << endl; return;
    
  vector<Double_t> nEventsv;     // number of events per sample
  vector<Double_t> nEventsVarv;  // error^2 on number of events per sample
  Double_t nEleDenomTL  = 0;        // number of tight+loose candidates
  Double_t nMuDenomTL  = 0;        // number of tight+loose candidates
  Double_t nTight    = 0;        // number of tight objects
  
  // electron + fake yield prediction
  Double_t nEleFOTL=0, nEleFOTLVarl=0, nEleFOTLVarh=0, nMuFOTL=0, nMuFOTLVarl=0, nMuFOTLVarh=0;

  // FILE *prs; prs = fopen("fake-print.txt","w");

  //
  // Set up output ntuple file for the sample
  //
  gSystem->mkdir(outputDir,kTRUE);
  TString outName = outputDir + TString("/") + TString("fakes_select.root");
  TFile *outFile = new TFile(outName,"RECREATE");
  TTree outtree("Events","Events");

  EmuData data;
  Double_t rawMet,rawprojvar, trigeff=1;
  UInt_t npt15jets;
  const UInt_t kMaxPt15Jets=50;
  TArrayF btagArray; btagArray.Set(kMaxPt15Jets); // array to hold b-tag values for pt-15 jets
  TArrayF jptArray; jptArray.Set(kMaxPt15Jets); // array to hold b-tag values for pt-15 jets
  TArrayF jetaArray; jetaArray.Set(kMaxPt15Jets); // array to hold b-tag values for pt-15 jets

  Float_t varl,varh;
    outtree.Branch("Events",&data.runNum,
"runNum/i:evtNum:lumiSec:nPV:njets:nbjets:vpt/F:vphi:rawmet:rawmetphi:met:metphi:mass:dphi:mt:pt:phi:pmet:pvis:eleiso:muiso:eled0:eled0sig:eleip3d:eleip3dsig:mud0:mud0sig:muip3d:muip3dsig:lpt1:leta1:lphi1:lpt2:leta2:lphi2:jpt1:jeta1:jphi1:jpt2:jeta2:jphi2:bjpt:bjeta:bjphi:mjj:svfmass:svfmassunc:genlpt1:genleta1:genlphi1:genlpt2:genleta2:genlphi2:weight:state/I");

  // extra branches
  outtree.Branch("npt15jets",&npt15jets);
  outtree.Branch("btagArray",&btagArray);
  outtree.Branch("jptArray",&jptArray);
  outtree.Branch("jetaArray",&jetaArray);
  outtree.Branch("trigeff",&trigeff);
  outtree.Branch("rawMet",&rawMet);
  outtree.Branch("rawprojvar",&rawprojvar);
  outtree.Branch("varl",&varl);
  outtree.Branch("varh",&varh);

  TFile *infile=0;
  TTree *eventTree=0;
    
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVfit");

  mithep::TSVfitter *fitter = new mithep::TSVfitter();

  
  //
  // loop over data samples
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  
  
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);

    const TString json = jsonv[ifile];
    Bool_t hasJSON = kFALSE;
    mithep::RunLumiRangeMap rlrm;
    if(json.CompareTo("NONE")!=0) { 
      hasJSON = kTRUE;
      rlrm.AddJSONFile(json.Data()); 
    }
        
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",    &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Muon",    &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Electron",&electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFJet",   &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("PV",      &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
    eventTree->SetBranchAddress("SVfitEMu",&svfitArr);    TBranch *svfitBr    = eventTree->GetBranch("SVfitEMu");
   
    Double_t weight=1;

    Double_t counter[30];         for(Int_t i=0; i<30; i++) { counter[i] = 0; }
   
    // perform fakeable objects extrapolation
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      
      infoBr->GetEntry(ientry);

      // check for certified runs
      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  
      
      // trigger requirement
      if(info->runNum <= 170053 && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdL])) continue;
      else if(info->runNum >  170053 && info->runNum <= 173199 && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdL])) continue;
      else if(info->runNum >  173199 && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL])) continue;

      // No good primary vertex? Skip to next event...
      if(!(info->hasGoodPV)) continue;
      
      pvArr->Clear();
      pvBr->GetEntry(ientry);
      
      // check for denominator objects and those passing tight selection

      vector<const mithep::TMuon*> goodMuonsv;
      vector<const mithep::TMuon*> looseMuonsv;

      vector<Double_t> muratev;
      vector<Double_t> murateErrlv;
      vector<Double_t> murateErrhv;

      muonArr->Clear();
      muonBr->GetEntry(ientry);

      UInt_t ntightmu=0;

      // loop over muons
      for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
        const mithep::TMuon* muon = (mithep::TMuon*)((*muonArr)[i]);

        Bool_t trigmatch = muon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdL_MuObj] || muon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdL_MuObj] || muon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj] || muon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj];
        //if(!trigmatch)                     continue;

	if(muon->pt < kMuonPt2Min)    continue;
	if(fabs(muon->eta) > 2.1)     continue;
	
	if(passMuonID(muon) && passMuonIsoPU(muon)) {
	  ntightmu++;
	  goodMuonsv.push_back(muon);
	  continue;
	}

        if(isMuonFO(muon,2)) {
          looseMuonsv.push_back(muon);
          Double_t fopt = (muon->pt < 35) ? muon->pt : 34.99;

          muratev.push_back(fFakeRate->MuonFakeRate(fopt, muon->eta, muon->phi) / (1-fFakeRate->MuonFakeRate(fopt, muon->eta, muon->phi)));
          murateErrlv.push_back(fFakeRate->MuonFakeRateStatErrorLow(fopt, muon->eta, muon->phi) / pow((1- fFakeRate->MuonFakeRate(fopt, muon->eta, muon->phi)),2));
          murateErrhv.push_back(fFakeRate->MuonFakeRateStatErrorHigh(fopt, muon->eta, muon->phi) / pow((1- fFakeRate->MuonFakeRate(fopt, muon->eta, muon->phi)),2));

        }

      }

      vector<const mithep::TElectron*> goodElectronsv;
      vector<const mithep::TElectron*> looseElectronsv;

      vector<Double_t> eleratev;
      vector<Double_t> elerateErrlv;
      vector<Double_t> elerateErrhv;
      
      electronArr->Clear();
      electronBr->GetEntry(ientry);

      UInt_t ntightele=0;

      // loop over electrons
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const mithep::TElectron* electron = (mithep::TElectron*)((*electronArr)[i]);
	
	if(electron->pt        < kElePt2Min)  continue;
	if(fabs(electron->eta) > 2.5)         continue;

        Bool_t trigmatch = electron->hltMatchBits[kHLT_Mu17_Ele8_CaloIdL_EGObj] || electron->hltMatchBits[kHLT_Mu8_Ele17_CaloIdL_EGObj] || electron->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj] || electron->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj];
	//if(!trigmatch)                          continue;

        Double_t mvaValue =  electronIDMVANoIPInfo->MVAValue(
                                          electron->pt,electron->scEta,
                                          electron->sigiEtaiEta,
                                          electron->deltaEtaIn,
                                          electron->deltaPhiIn,
                                          electron->HoverE,
                                          electron->d0,
                                          electron->dz,
                                          electron->fBrem,
                                          electron->EoverP,
                                          electron->ESeedClusterOverPOut,
                                          TMath::Sqrt(electron->sigiPhiiPhi),
                                          electron->nBrem,
                                          (1.0/(electron->scEt * TMath::CosH(electron->scEta)) - 1/electron->p),
                                          electron->ESeedClusterOverPIn,
                                          electron->ip3d,
                                          electron->ip3dSig );

	if(passEleMVAID(electron,mvaValue) && passEleIsoPU(electron)) {
	  ntightele++;
	  goodElectronsv.push_back(electron);
	  continue;
	}
	
	if(isEleFO(electron)) {
	  looseElectronsv.push_back(electron);
	  Double_t fopt = (electron->pt < 35) ? electron->pt : 34.99;

	  eleratev.push_back(fFakeRate->ElectronFakeRate(fopt, electron->eta, electron->phi) / (1-fFakeRate->ElectronFakeRate(fopt, electron->eta, electron->phi)));
	  // assert(ratev.back() > 0);
	  elerateErrlv.push_back(fFakeRate->ElectronFakeRateStatErrorLow(fopt, electron->eta, electron->phi) / pow((1- fFakeRate->ElectronFakeRate(fopt, electron->eta, electron->phi)),2));
	  elerateErrhv.push_back(fFakeRate->ElectronFakeRateStatErrorHigh(fopt, electron->eta, electron->phi) / pow((1- fFakeRate->ElectronFakeRate(fopt, electron->eta, electron->phi)),2));
	}
      }
          
      // Tight+Loose
      if(ntightmu>0 && ntightele==0) {

	const mithep::TMuon* mu = goodMuonsv[0];

	for(UInt_t iel=0; iel<looseElectronsv.size(); iel++) {
	  const mithep::TElectron* ele = looseElectronsv[iel];

          if(mu->q == ele->q) continue;                    

          if(mu->pt < kMuonPt2Min  || ele->pt < kElePt2Min) continue;
	  if(mu->pt < kMuonPt1Min  && ele->pt < kElePt1Min) continue;

	  if(mu->pt < kMuonPt1Min) {
            if(!(info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] || info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL])) continue; // if failed trig1
          }
          else if(ele->pt < kElePt1Min) {
            if(!(info->triggerBits[kHLT_Mu17_Ele8_CaloIdL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL])) continue; // if failed trig2
	  }

          Double_t muiso = muonIsoPU(mu);
          Double_t eleiso = eleIsoPU(ele);
	  
	  TLorentzVector lep1, lep2, dilep;  // lepton 4-vectors
	  Int_t finalState=-1;	           // final state type
          Double_t svfmass = -999;
          Double_t svfmassunc = -999;

	  if(mu->pt > ele->pt) {
	    lep1.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
	    lep2.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
	  } else {
	    lep1.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
	    lep2.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
	  }
	  dilep = lep1+lep2;

	  if(ele->pt > mu->pt) finalState=kEleMu; 
	  else                 finalState=kMuEle;

          Double_t met=info->pfMET,metphi=info->pfMETphi;

          // SVFit info
          svfitArr->Clear();
          svfitBr->GetEntry(ientry);

          for(Int_t i = 0; i < svfitArr->GetEntriesFast(); i++) {
            mithep::TSVfit *svfit = (mithep::TSVfit*) svfitArr->At(i);
            Int_t id = 0;
            if(toolbox::deltaR(lep1.Eta(),lep1.Phi(),svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.01            ) id = 1;
            if(toolbox::deltaR(lep2.Eta(),lep2.Phi(),svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.01 && id == 0) id = 2;
            if(id == 0) continue;
            if(toolbox::deltaR(lep1.Eta(),lep1.Phi(),svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.01 && id == 2) id = 3;
            if(toolbox::deltaR(lep2.Eta(),lep2.Phi(),svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.01 && id == 1) id = 4;
            if(id < 3) continue;
            TLorentzVector svf = fitter->fit(svfit,met,metphi);
            svfmass = svf.M();
            svfmassunc = fitter->massUnc();
          }


	  // loop through jets      
	  jetArr->Clear();
	  jetBr->GetEntry(ientry);
	  UInt_t njets = 0, nbjets = 0;
	  const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	  btagArray.Reset();          jptArray.Reset();        jetaArray.Reset();  npt15jets=0;
	  for(Int_t ijet=0; ijet<jetArr->GetEntriesFast(); ijet++) {
	    const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[ijet]);

	    if(toolbox::deltaR(jet->eta,jet->phi,lep1.Eta(),lep1.Phi()) < 0.3) continue;
	    if(toolbox::deltaR(jet->eta,jet->phi,lep2.Eta(),lep2.Phi()) < 0.3) continue;

	    if(fabs(jet->eta) > 5) continue;

	    // look for b-jets
	    if((jet->pt > kBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
	      assert(npt15jets<kMaxPt15Jets);
	      btagArray.AddAt(jet->csv,npt15jets); npt15jets++;
	      if(jet->tche > 3.3) {
		nbjets++;
		if(!bjet || jet->pt > bjet->pt)
		  bjet = jet; // leading b-jet
	      }
	    }

	    // look for vbf jets
	    if(jet->pt > kJetPtMin) {
              assert(njets<kMaxPt15Jets);
              jptArray.AddAt(jet->pt,njets);
              jetaArray.AddAt(jet->eta,njets);
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
	  rawprojvar  = 0.85*projVis - projMet;

	  rawMet = info->pfMET;
//	  Double_t met=info->pfMET,metphi=info->pfMETphi;
	  metv.SetPtEtaPhi(met,0,metphi);
	  projMet  =   metv.Dot(bisector);

	  const Double_t r    = eleratev[iel];
	  if(eleratev[iel]<0) cout << "error: " << ele->pt << " " << eleratev[iel] << endl;
          const Double_t errl = elerateErrlv[iel];
	  const Double_t errh = elerateErrhv[iel];
	  
	  nEleFOTL     += weight*r;
	  nEleFOTLVarl += weight*weight*(r*r+errl*errl);
	  nEleFOTLVarh += weight*weight*(r*r+errh*errh);
	  
	  nEleDenomTL+=weight;
	  
	  varl    = weight*weight*(r*r+errl*errl);
	  varh    = weight*weight*(r*r+errh*errh);

	  // fprintf(prs,"%14d%14d%8.2f\n",info->runNum,info->evtNum,r);

          data.runNum   = info->runNum;
          data.evtNum   = info->evtNum;
          data.lumiSec  = info->lumiSec;
          data.nPV      = pvArr->GetEntriesFast();
          data.njets    = njets;
          data.nbjets   = nbjets;
          data.vpt      = 0;
          data.vphi     = 0;
          data.rawmet   = met;
          data.rawmetphi= metphi;
          data.met      = met;
          data.metphi   = metphi;
          data.mass     = dilep.M();
          data.dphi     = toolbox::deltaPhi(lep1.Phi(),lep2.Phi());
          data.mt       = sqrt( 2.0 * (dilep.Pt()) * met * (1.0-cos(toolbox::deltaPhi(dilep.Phi(),metphi))) );
          data.pt       = dilep.Pt();
          data.phi      = dilep.Phi();
          data.pmet     = projMet;
          data.pvis     = projVis;
          data.eleiso   = eleiso;
          data.muiso    = muiso;
          data.eled0    = ele->d0;
          data.eled0sig = ele->d0Sig;
          data.eleip3d  = ele->ip3d;
          data.eleip3dsig = ele->ip3dSig;
          data.mud0     = mu->d0;
          data.mud0sig  = mu->d0Sig;
          data.muip3d   = mu->ip3d;
          data.muip3dsig  = mu->ip3dSig;
          data.lpt1     = lep1.Pt();
          data.leta1    = lep1.Eta();
          data.lphi1    = lep1.Phi();
          data.lpt2     = lep2.Pt();
          data.leta2    = lep2.Eta();
          data.lphi2    = lep2.Phi();
          data.jpt1     = (jet1) ? jet1->pt  : 0;
          data.jeta1    = (jet1) ? jet1->eta : 0;
          data.jphi1    = (jet1) ? jet1->phi : 0;
          data.jpt2     = (jet2) ? jet2->pt  : 0;
          data.jeta2    = (jet2) ? jet2->eta : 0;
          data.jphi2    = (jet2) ? jet2->phi : 0;
          data.bjpt     = (bjet) ? bjet->pt  : 0;
          data.bjeta    = (bjet) ? bjet->eta : 0;
          data.bjphi    = (bjet) ? bjet->phi : 0;
          data.mjj      = (njets>1) ? dijet.M() : 0;
          data.svfmass  = svfmass;
          data.svfmassunc = svfmassunc;
          data.genlpt1  = 0;
          data.genleta1 = 0;
          data.genlphi1 = 0;
          data.genlpt2  = 0;
          data.genleta2 = 0;
          data.genlphi2 = 0;
          data.weight   = weight*r;
          data.state    = finalState;

          outtree.Fill();
	}
      }

      else if(ntightele>0 && ntightmu==0) {

        const mithep::TElectron* ele = goodElectronsv[0];

        for(UInt_t imu=0; imu<looseMuonsv.size(); imu++) {
          const mithep::TMuon* mu = looseMuonsv[imu];

          if(mu->q == ele->q) continue;

          if(mu->pt < kMuonPt2Min  || ele->pt < kElePt2Min) continue;
          if(mu->pt < kMuonPt1Min  && ele->pt < kElePt1Min) continue;

          if(mu->pt  < kMuonPt1Min) {
            if(!(info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] || info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL])) continue; // if failed trig1
          }
          else if(ele->pt < kElePt1Min) {
            if(!(info->triggerBits[kHLT_Mu17_Ele8_CaloIdL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL])) continue; // if failed trig2
          }

          Double_t muiso = muonIsoPU(mu);
          Double_t eleiso = eleIsoPU(ele);

          TLorentzVector lep1, lep2, dilep;  // lepton 4-vectors
          Int_t finalState=-1;             // final state type
          Double_t svfmass = -999;
          Double_t svfmassunc = -999;

          if(mu->pt > ele->pt) {
            lep1.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
            lep2.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
          } else {
            lep1.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
            lep2.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
          }
          dilep = lep1+lep2;

          if(ele->pt > mu->pt) finalState=kEleMu;
          else                 finalState=kMuEle;

          Double_t met=info->pfMET,metphi=info->pfMETphi;

          // SVFit info
          svfitArr->Clear();
          svfitBr->GetEntry(ientry);

          for(Int_t i = 0; i < svfitArr->GetEntriesFast(); i++) {
            mithep::TSVfit *svfit = (mithep::TSVfit*) svfitArr->At(i);
            Int_t id = 0;
            if(toolbox::deltaR(lep1.Eta(),lep1.Phi(),svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.01            ) id = 1;
            if(toolbox::deltaR(lep2.Eta(),lep2.Phi(),svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.01 && id == 0) id = 2;
            if(id == 0) continue;
            if(toolbox::deltaR(lep1.Eta(),lep1.Phi(),svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.01 && id == 2) id = 3;
            if(toolbox::deltaR(lep2.Eta(),lep2.Phi(),svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.01 && id == 1) id = 4;
            if(id < 3) continue;
            TLorentzVector svf = fitter->fit(svfit,met,metphi);
            svfmass = svf.M();
            svfmassunc = fitter->massUnc();
          }


          // loop through jets      
          jetArr->Clear();
          jetBr->GetEntry(ientry);
          UInt_t njets = 0, nbjets = 0;
          const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
          btagArray.Reset();          jptArray.Reset();        jetaArray.Reset();  npt15jets=0;
          for(Int_t ijet=0; ijet<jetArr->GetEntriesFast(); ijet++) {
            const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[ijet]);

            if(toolbox::deltaR(jet->eta,jet->phi,lep1.Eta(),lep1.Phi()) < 0.3) continue;
            if(toolbox::deltaR(jet->eta,jet->phi,lep2.Eta(),lep2.Phi()) < 0.3) continue;

            if(fabs(jet->eta) > 5) continue;

            // look for b-jets
            if((jet->pt > kBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
              assert(npt15jets<kMaxPt15Jets);
              btagArray.AddAt(jet->csv,npt15jets); npt15jets++;
              if(jet->csv > 0.679) {
                nbjets++;
                if(!bjet || jet->pt > bjet->pt)
                  bjet = jet; // leading b-jet
              }
            }

            // look for vbf jets
            if(jet->pt > kJetPtMin) {
              assert(njets<kMaxPt15Jets);
              jptArray.AddAt(jet->pt,njets);
              jetaArray.AddAt(jet->eta,njets);
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
          rawprojvar  = 0.85*projVis - projMet;

          rawMet = info->pfMET;
//        Double_t met=info->pfMET,metphi=info->pfMETphi;
          metv.SetPtEtaPhi(met,0,metphi);
          projMet  =   metv.Dot(bisector);

          const Double_t r    = muratev[imu];
          if(muratev[imu]<0) cout << "error: " << mu->pt << " " << muratev[imu] << endl;
          const Double_t errl = murateErrlv[imu];
          const Double_t errh = murateErrhv[imu];

          nMuFOTL     += weight*r;
          nMuFOTLVarl += weight*weight*(r*r+errl*errl);
          nMuFOTLVarh += weight*weight*(r*r+errh*errh);

          nMuDenomTL+=weight;

          varl    = weight*weight*(r*r+errl*errl);
          varh    = weight*weight*(r*r+errh*errh);

          // fprintf(prs,"%14d%14d%8.2f\n",info->runNum,info->evtNum,r);

          data.runNum   = info->runNum;
          data.evtNum   = info->evtNum;
          data.lumiSec  = info->lumiSec;
          data.nPV      = pvArr->GetEntriesFast();
          data.njets    = njets;
          data.nbjets   = nbjets;
          data.vpt      = 0;
          data.vphi     = 0;
          data.rawmet   = met;
          data.rawmetphi= metphi;
          data.met      = met;
          data.metphi   = metphi;
          data.mass     = dilep.M();
          data.dphi     = toolbox::deltaPhi(lep1.Phi(),lep2.Phi());
          data.mt       = sqrt( 2.0 * (dilep.Pt()) * met * (1.0-cos(toolbox::deltaPhi(dilep.Phi(),metphi))) );
          data.pt       = dilep.Pt();
          data.phi      = dilep.Phi();
          data.pmet     = projMet;
          data.pvis     = projVis;
          data.eleiso   = eleiso;
          data.muiso    = muiso;
          data.eled0    = ele->d0;
          data.eled0sig = ele->d0Sig;
          data.eleip3d  = ele->ip3d;
          data.eleip3dsig = ele->ip3dSig;
          data.mud0     = mu->d0;
          data.mud0sig  = mu->d0Sig;
          data.muip3d   = mu->ip3d;
          data.muip3dsig  = mu->ip3dSig;
          data.lpt1     = lep1.Pt();
          data.leta1    = lep1.Eta();
          data.lphi1    = lep1.Phi();
          data.lpt2     = lep2.Pt();
          data.leta2    = lep2.Eta();
          data.lphi2    = lep2.Phi();
          data.jpt1     = (jet1) ? jet1->pt  : 0;
          data.jeta1    = (jet1) ? jet1->eta : 0;
          data.jphi1    = (jet1) ? jet1->phi : 0;
          data.jpt2     = (jet2) ? jet2->pt  : 0;
          data.jeta2    = (jet2) ? jet2->eta : 0;
          data.jphi2    = (jet2) ? jet2->phi : 0;
          data.bjpt     = (bjet) ? bjet->pt  : 0;
          data.bjeta    = (bjet) ? bjet->eta : 0;
          data.bjphi    = (bjet) ? bjet->phi : 0;
          data.mjj      = (njets>1) ? dijet.M() : 0;
          data.svfmass  = svfmass;
          data.svfmassunc = svfmassunc;
          data.genlpt1  = 0;
          data.genleta1 = 0;
          data.genlphi1 = 0;
          data.genlpt2  = 0;
          data.genleta2 = 0;
          data.genlphi2 = 0;
          data.weight   = weight*r;
          data.state    = finalState;

	  outtree.Fill();
	}
      }      
    }    
    delete infile;
    infile=0, eventTree=0;
  }
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  delete info;
  delete muonArr;
  delete electronArr;
  delete jetArr;
  delete pvArr;
  delete svfitArr;

  // fclose(prs);
  
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  cout << "   data samples:" << endl;
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "     " << fnamev[ifile] << endl;
  }
  cout << endl;
  cout << "  >>>   Events from TL: " << nEleFOTL << " (+" << sqrt(nEleFOTLVarh) << "/-" << sqrt(nEleFOTLVarl) << ")" << endl;
  cout << "  >>>   Electron FOs from TL: " << nEleDenomTL << endl;  
  cout << endl;
  cout << endl;
  cout << "  >>>   Events from TL: " << nMuFOTL << " (+" << sqrt(nMuFOTLVarh) << "/-" << sqrt(nMuFOTLVarl) << ")" << endl;
  cout << "  >>>   Muon FOs from TL: " << nMuDenomTL << endl;  
  cout << endl;
  
  ofstream txtfile;
  char txtfname[100];
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << endl;

  txtfile << "   data samples:" << endl;
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "     " << fnamev[ifile] << endl;
  }
  txtfile << endl;
  txtfile << "  >>>   Events from TL: " << nEleFOTL << " (+" << sqrt(nEleFOTLVarh) << "/-" << sqrt(nEleFOTLVarl) << ")" << endl;
  txtfile << "  >>>   Electron FOs from TL: " << nEleDenomTL << endl;
  txtfile << "  >>>      Tight muons: " << nTight << endl;
  txtfile << endl;
  txtfile << endl;
  txtfile << "  >>>   Events from TL: " << nMuFOTL << " (+" << sqrt(nMuFOTLVarh) << "/-" << sqrt(nMuFOTLVarl) << ")" << endl;
  txtfile << "  >>>   Muon FOs from TL: " << nMuDenomTL << endl;  
  txtfile << "  >>>      Tight muons: " << nTight << endl;
  txtfile << endl;
  txtfile.close();
  
  cout << " <> Output saved in " << outputDir << "/" << endl;    
  cout << endl; 


  
  gBenchmark->Show("predictFakesMuEle"); 
} 
