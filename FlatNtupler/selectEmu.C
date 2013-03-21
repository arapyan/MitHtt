#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TH1.h>                    // histogram base class
#include <TH2.h>                    // histogram base class
#include <TNtuple.h>                // class to access ntuples
#include <TTree.h>                  // class to access trees
#include <TRegexp.h>                // ROOT regexp class
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TMath.h>                  // ROOT math functions
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "MitHtt/Common/CSample.hh"        // helper class for organizing input ntuple files
#include "MitHtt/Common/MyTools.hh"        // miscellaneous helper functions

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh" 
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"   
#include "MitHtt/Ntupler/interface/TVertex.hh"   
#include "MitHtt/Ntupler/interface/TSVfit.h"
#include "MitHtt/Ntupler/interface/TSVfitter.h"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh"

// scale factros
#include "MitHtt/Utils/DataMC.hh"

// B-tag scale factors
#include "BtagSF.hh"

// format output
#include "Output.hh"

#endif

//=== FUNCTION DECLARATIONS ======================================================================================
const Double_t pi = 3.14159265358979;

//=== MAIN MACRO =================================================================================================

void selectEmu(const TString conf,         // input config file
               const TString outputDir,         // output directory
	       const Double_t lumi,             // luminosity pb^-1
               const UInt_t btageff=0,          // b-tag efficiency scale factor uncertainty
               const UInt_t jetunc=0,           // jet energy uncertainties
               const UInt_t mistag=0,           // b mistag rate scale factor uncertainty
	       const UInt_t elescale=0,         // electron energy scale/resolution uncertainty
	       const Bool_t is2012=kTRUE
) {
  gBenchmark->Start("selectEmu");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  vector<TString>  snamev;      // sample name (for output file)  
  vector<CSample*> samplev;     // data/MC samples
  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    if(line[0]=='$') {
      // fakes come from a separate macro
      if((TString(line).Contains("fake")) && (state>0)) continue;

      samplev.push_back(new CSample());
      stringstream ss(line);
      string chr;
      string sname;
      Int_t color;
      ss >> chr >> sname >> color;
      string label = line.substr(line.find('@')+1);
      snamev.push_back(sname);
      samplev.back()->label = label;
      samplev.back()->color = color;
      continue;
    }
    
    if(state==0) {  // define data sample
      string fname;
      string json;
      Int_t type;
      stringstream ss(line);
      ss >> fname >> type >> json;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(type);
      samplev.back()->xsecv.push_back(0);
      samplev.back()->jsonv.push_back(json);
    
    } else if(state==1) {  // define MC samples
      string fname;
      Double_t xsec;
      stringstream ss(line);
      ss >> fname >> xsec;
      if(TString(fname).Contains("dummy",TString::kIgnoreCase)) continue;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(0);
      samplev.back()->xsecv.push_back(xsec);
    }
  }
  ifs.close();

  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);

  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };  // dataset type  
  
  const Double_t kMuonPt1Min = 20;
  const Double_t kMuonPt2Min = 10;
  
  const Double_t kElePt1Min  = 20;
  const Double_t kElePt2Min  = 10;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 20;
  
  Bool_t doKFactors = !is2012;      // not needed in Summer12

  Bool_t doNpuRwgt = kTRUE;

  // Access samples and fill histograms
  TFile *infile=0;
  TTree *eventTree=0;  
  TTree *nEventTree=0;  
  
  BtagSF* btsf = new BtagSF(12345);

  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVfit");

  Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // loop over samples
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if(isam==0 && !hasData) continue;
  
    CSample* samp = samplev[isam];

    Double_t nSelEvents=0;

    // Set up output ntuple file for the sample
    TString outfname = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    Output* out = new Output(outfname);  

    // loop through files
    cout <<  "processing " << snamev[isam] << ":" << endl;
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      printf("        %-55s",(samp->fnamev[ifile]+"...").Data()); fflush(stdout);
      infile = new TFile(samp->fnamev[ifile]); 
      assert(infile);

      TString sfname    = samp->fnamev[ifile];
      TString basename = sfname(sfname.Last('/')+1,sfname.Last('.') - sfname.Last('/') - 1);

      // which corrections to apply where
      Bool_t isdata     = !(samp->typev[ifile]==eMC);
      Bool_t is52mc     = sfname.Contains("s12-") && sfname.Contains("v9");
      Bool_t is53mc     = sfname.Contains("s12-") && (sfname.Contains("v7a") || sfname.Contains("v7c"));
      Bool_t isemb      = snamev[isam].Contains("emb");
      Bool_t isfall11   = sfname.Contains("f11");
      Bool_t issamesign = snamev[isam].Contains("ss-fakes");
      Bool_t doRecoil   = (sfname.Contains("ztt") || sfname.Contains("-zll") || sfname.Contains("zjets") || snamev[isam].Contains("_sm_") || snamev[isam].Contains("_mssm_")) && !isemb;
      Bool_t reallyDoKf = doKFactors && sfname.Contains("-gf-");
      Bool_t ismadz     = sfname.Contains("-zll") || sfname.Contains("-zjets"); // madgraph z samples
      Bool_t ismadzll   = snamev[isam].Contains("zll") && (sfname.Contains("-zll") || sfname.Contains("-zjets")); // madgraph z samples
      Bool_t ismssm     = sfname.Contains("-ggh-") || sfname.Contains("-bbh-");
      Bool_t issm       = sfname.Contains("-gf-") || sfname.Contains("-vbf-") || sfname.Contains("-vtth-");
      Bool_t doIdScale  = !isdata || isemb;
      Bool_t doTrigScale= !isdata || isemb;
      Bool_t getGen     = doRecoil || reallyDoKf || ismadz ||isemb || ismssm || issm;
      Bool_t doJetUnc   = (jetunc!=kNo);

      int recoiltype = 0;
      if(issm || ismssm) recoiltype = 1;
      else if(doRecoil) recoiltype = 2;
      out->setupRecoil(recoiltype, is2012, 1);

      // PU reweighting
      TString pileupReweightFile;
      if(isfall11) {
	cout << "Fall11 sample!" << endl;
	pileupReweightFile = "/data/blue/vdutta/htt/pileup/PUWeights_F11To2011.root";
      } else if (is52mc) {
	cout << "52X sample!" << endl;
	pileupReweightFile = "/data/blue/vdutta/htt/pileup/PUWeights_S12To2012_5088ipb_true_noLowPU.root";
      } else if (is53mc) {
	cout << "53X sample!" << endl;
	pileupReweightFile = "/data/blue/vdutta/htt/pileup/PUWeights_S1253XTo2012ABCD.root";
      } else {
	cout << "data" << endl;
      }

      TH1D *puWeights = 0;
      TFile *pufile = new TFile(pileupReweightFile.Data());
      puWeights = (TH1D*)pufile->Get("puWeights");

      // setup selecting with JSON file, if necessary
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if(isdata && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	ifstream jsonchk; jsonchk.open(samp->jsonv[ifile].Data()); assert(jsonchk.is_open()); jsonchk.close();
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }

      // k-factors
      TH1D *hKFactors = (reallyDoKf) ? kfFHPInit(higgsmass(basename)) : 0;

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
      nEventTree = (TTree*)infile->Get("hEvents"); assert(nEventTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");      
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
      eventTree->SetBranchAddress("SVfitEMu", &svfitArr);    TBranch *svfitBr    = eventTree->GetBranch("SVfitEMu");
      TBranch *genBr=0;
      if(getGen) {
        eventTree->SetBranchAddress("Gen", &gen);
        genBr = eventTree->GetBranch("Gen");
      }

      // get weights for MC
      Double_t weight=1,treeEntries=-1; // (weight is only initialized for each *file*)
      if(!isdata) {
	if(sfname.Contains("_skim.root")) treeEntries = unskimmedEntries(sfname); // get entries from unskimmed file
	else                              treeEntries = (Double_t)nEventTree->GetEntries();
	assert(treeEntries>0);
        weight = lumi*(samp->xsecv[ifile])/treeEntries;                           // (assumes you've merged filesets)
	if(isemb)  weight=1.0;
      }
      samp->weightv.push_back(weight);

      // counters
      Double_t nsel=0, nselvar=0; 

      cout << eventTree->GetEntries() << " events" << endl;

      // loop over events
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	if(ientry%1000000 == 0) cout << "processing " << ientry << endl;
        infoBr->GetEntry(ientry);

	if(getGen)  genBr->GetEntry(ientry);

	// skip non-tau events in madgraph sample
	if(ismadz && !ismadzll && (fabs(gen->id_1_a)<15 || fabs(gen->id_1_a)>19)) continue;

        // skip non-ee/mumu events in madgraph sample for zll
        if(ismadzll && (fabs(gen->id_1_a)>14 && fabs(gen->id_1_a)<20)) continue;

	// certified run selection
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;

	// trigger
	if(is2012) {
	  if(!isemb && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL])) continue;
	}
	else {
	  if(!isemb && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdL] || info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL])) continue;
	}

        // good primary vertex
        if(!info->hasGoodPV) continue;
	pvArr->Clear();
	pvBr->GetEntry(ientry);	

        // loop through muons
        vector<const mithep::TMuon*> goodMuonsv;
        vector<const mithep::TMuon*> looseMuonsv;
	Int_t nSelMuons = 0;
        muonArr->Clear();
        muonBr->GetEntry(ientry);

        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const mithep::TMuon *muon = (mithep::TMuon*)((*muonArr)[i]);
	  looseMuonsv.push_back(muon);
	  Bool_t passTrigMatch = kTRUE;

	  // trigger matching
	  Bool_t trigmatch = 0;
	  if(is2012) {trigmatch = ((info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && muon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj]) || (info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && muon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj]));}
	  else {trigmatch = ((info->triggerBits[kHLT_Mu17_Ele8_CaloIdL] && muon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdL_MuObj]) || (info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] && muon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdL_MuObj]) || (info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL] && muon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj]) || (info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL] && muon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj]));}
	  if(!isemb && !trigmatch) passTrigMatch = kFALSE;

	  if(!(muon->typeBits & kGlobal))	continue;
          if(muon->pt < kMuonPt2Min)		continue;
	  if(fabs(muon->eta) > 2.4)		continue;
	  if(!issamesign) {
	    if(passTightPFMuonID(muon,1) && passMuonIsoPU(muon,3))  goodMuonsv.push_back(muon);
	  } else {
	    if(passTightPFMuonID(muon,1))  goodMuonsv.push_back(muon);
	  }
          if(passTightPFMuonID(muon,0) && passTrigMatch) {
	    if(!issamesign && !passMuonIsoPU(muon,0)) continue;
	    nSelMuons++;
	  }
        }
	
        // loop through electrons 
        vector<const mithep::TElectron*> goodElectronsv;
        Int_t nSelElectrons = 0;
        electronArr->Clear();
        electronBr->GetEntry(ientry);

        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
	  const mithep::TElectron *electron = (mithep::TElectron*)((*electronArr)[i]);
          Bool_t passTrigMatch = kTRUE;

	  // trigger matching
	  Bool_t trigmatch = 0;
	  if(is2012) {trigmatch = ((info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && electron->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj]) || (info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && electron->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj]));}
	  else {trigmatch = ((info->triggerBits[kHLT_Mu17_Ele8_CaloIdL] && electron->hltMatchBits[kHLT_Mu17_Ele8_CaloIdL_EGObj]) || (info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL] && electron->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj]) || (info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] && electron->hltMatchBits[kHLT_Mu8_Ele17_CaloIdL_EGObj]) || (info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL] && electron->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj]));}
	  if(!isemb && !trigmatch) passTrigMatch = kFALSE;

          if(fabs(electron->eta) > 2.5)		continue;
	  if(elescale==kNo && electron->pt<kElePt2Min) continue;
	  else if(elescale==kDown && electron->pt<0.99*kElePt2Min) continue;
	  else if(elescale==kUp && electron->pt<1.01*kElePt2Min) continue;

	  // track matching against muons
          Bool_t hasMuonTrack=kFALSE;
          for(UInt_t imu=0; imu<goodMuonsv.size(); imu++) {
            if(electron->trkID == goodMuonsv[imu]->trkID) hasMuonTrack=kTRUE;
          }

	  // clean against loose muons
          Bool_t matchLooseMuon=kFALSE;
	  for(UInt_t imu=0;imu<looseMuonsv.size();imu++) {
	    const mithep::TMuon *mu = looseMuonsv[imu];
	    if(toolbox::deltaR(electron->eta,electron->phi,mu->eta,mu->phi) < 0.3) matchLooseMuon=kTRUE;
	  }

	  if(!issamesign) {
	    if(pass2012EleMVAID(electron,kLoose,1) && passEleIsoPU(electron,3))  goodElectronsv.push_back(electron);
	  } else {
	    if(pass2012EleMVAID(electron,kLoose,1))  goodElectronsv.push_back(electron);
	  }
	  if(fabs(electron->eta)<2.3 && pass2012EleMVAID(electron,kLoose,0) && !hasMuonTrack && !matchLooseMuon && passTrigMatch) {
	    if(!issamesign && !passEleIsoPU(electron,0)) continue;
	    nSelElectrons++;
	  }
        }

	//----------------------------------------------------------------------------------------

	if(nSelMuons < 1 || nSelElectrons < 1) continue;

	if(goodMuonsv.size()>1 || goodElectronsv.size()>1) continue;     // veto events with more than 1 electron or muon passing ID

	const mithep::TMuon *mu	   = goodMuonsv[0];
	const mithep::TElectron *ele = goodElectronsv[0];

	Double_t elept = ele->pt;
	if(elescale==kUp) {
	  elept = 1.01*ele->pt;
	}
	if(elescale==kDown) {
          elept = 0.99*ele->pt;
        }

	if(mu->pt < kMuonPt2Min  || elept < kElePt2Min) continue;
	if(mu->pt < kMuonPt1Min  && elept < kElePt1Min) continue;

	// trigger requirements
	if(mu->pt  < kMuonPt1Min) {
	  if(is2012) {
            if(!isemb && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL])) continue; // if failed trig1
	  } else {
	    if(!isemb && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] || info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL])) continue; // if failed trig1
	  }
	}
	else if(elept < kElePt1Min) {
	  if(is2012) {
	    if(!isemb && !(info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL])) continue; // if failed trig2
	  } else {
	    if(!isemb && !(info->triggerBits[kHLT_Mu17_Ele8_CaloIdL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL])) continue; // if failed trig2
	  }
	}

	// same-sign requirements
	if(issamesign) {
	  if(mu->q != ele->q) continue;
	}
	else {
	  if(mu->q == ele->q) continue;
	} 

        // SVFit
        svfitArr->Clear();
        svfitBr->GetEntry(ientry);

        for(Int_t i = 0; i < svfitArr->GetEntriesFast(); i++) {
          mithep::TSVfit *svfit = (mithep::TSVfit*) svfitArr->At(i);
          Int_t id = 0;
          if(toolbox::deltaR(ele->eta,ele->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.01           ) id = 1;
          if(toolbox::deltaR(mu->eta,mu->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi())   < 0.01 && id == 0) id = 2;
          if(id == 0) continue;
          if(toolbox::deltaR(ele->eta,ele->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.01 && id == 2) id = 3;
          if(toolbox::deltaR(mu->eta,mu->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi())   < 0.01 && id == 1) id = 4;
          if(id < 3) continue;
	  out->fillCov(svfit);
        }

	out->fillMuon(mu,0,muonIsoPU(mu),passMuonIsoPU(mu,0));
	out->fillElectron(ele,1,eleIsoPU(ele),passEleIsoPU(ele,0), elescale);

        // loop through jets      
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets = 0, nbjets = 0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet1=0, *bjet2=0;
	out->btagArray.Reset();	out->jptArray.Reset();	out->jetaArray.Reset();	UInt_t npt20jets=0;
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
	  mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	  double jesunc = (is2012 && fabs(jet->eta) > 2.5) ? 2.0*jet->unc : jet->unc;

	  if(doJetUnc) jet->pt *= (jetunc==kDown) ? (1-jesunc) : (1+jesunc);

          if(toolbox::deltaR(jet->eta,jet->phi,ele->eta,ele->phi) < 0.5) continue;
          if(toolbox::deltaR(jet->eta,jet->phi,mu->eta,mu->phi) < 0.5) continue;

          if(fabs(jet->eta) > 4.7) continue;
	  if(!jet->id) continue;

	  // look for b-jets
	  Bool_t btagged = btsf->isbtagged(jet->pt, jet->eta, jet->csv, jet->matchedId, (isdata||isemb) ,btageff, mistag, is2012);
	  if((jet->pt > kBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
	    assert(npt20jets<50);
	    out->btagArray.AddAt(jet->csv,npt20jets);
	    npt20jets++;
	    if(btagged) {
	      nbjets++;
	      if(!bjet1 || jet->pt > bjet1->pt) {
		bjet2 = bjet1; // leading b-jet
		bjet1 = jet;
	      } else  if(!bjet2 || jet->pt > bjet2->pt) {
		bjet2 = jet;
	      }
	    }
	  }

	  // look for jets
          if(jet->pt > kJetPtMin) {
            assert(njets<50);
	    out->jptArray.AddAt(jet->pt,njets);
	    out->jetaArray.AddAt(jet->eta,njets);
  	    njets++;
	    if(!jet1 || jet->pt > jet1->pt) { // jet1 is highest pt, jet2 next highest
	      jet2 = jet1;
	      jet1 = jet;
	    } else if(!jet2 || jet->pt > jet2->pt) {
	      jet2 = jet;
	    }
          }		    
        }

	// dijet system
	Int_t nCentralJets=0;
	if(njets>1) {
          for(Int_t i=2; i<jetArr->GetEntriesFast(); i++) {
            mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	    if(!(jet->pt > kJetPtMin && fabs(jet->eta)<4.7 && jet->id==1)) continue;
            if(toolbox::deltaR(jet->eta,jet->phi,ele->eta,ele->phi) < 0.5) continue;
            if(toolbox::deltaR(jet->eta,jet->phi,mu->eta,mu->phi) < 0.5) continue;
	    if(jet1->eta > jet2->eta && jet->eta > jet2->eta && jet->eta < jet1->eta) nCentralJets++;
	    else if(jet2->eta > jet1->eta && jet->eta > jet1->eta && jet->eta < jet2->eta) nCentralJets++;
	  }
        }

	out->fillJets(jet1,jet2,bjet1,bjet2,njets,nbjets,npt20jets,nCentralJets);

        // get k-factor if necessary
        Double_t kf=1;
        if(reallyDoKf) kf = kfFHPValue(gen->vpt_a, hKFactors);

	// do vertex reweighting
	Double_t npuWgt = 1;
	if(!isdata && !isemb && doNpuRwgt) {
	  assert(puWeights);
	  Int_t npuxbin = puWeights->GetXaxis()->FindFixBin(TMath::Min(double(info->nPUTrue), 59.999));
	  npuWgt = puWeights->GetBinContent(npuxbin);
	}

	// lepton ID corrections
	Double_t idscale = 1;
	if(doIdScale) idscale = muIDscaleEmu(mu->pt,mu->eta,is2012)*eleIDscaleEmu(elept,ele->eta,is2012);

	// trigger scale factor for MC
	Double_t trigscale = 1;
	if(doTrigScale && !isemb) trigscale=muTrigScaleEmu(mu->pt,mu->eta,is2012)*eleTrigScaleEmu(elept,ele->eta,is2012);
	if(doTrigScale && isemb) trigscale=muTrigEffEmu(mu->pt,mu->eta,is2012)*eleTrigEffEmu(elept,ele->eta,is2012);

	// embedding weight for embedded sample
	Double_t embWgt = 1;
	if(!isdata || isemb) out->fillGen(gen);
	if(isemb) embWgt=info->embWeight;

	out->fMCWeight	 = weight*kf*embWgt/lumi;
	out->fPUWeight	 = npuWgt;
	out->fEffWeight	 = trigscale*idscale;
	out->fWeight	 = weight*kf*npuWgt*trigscale*idscale*embWgt/lumi;
	double scalecorr = 0;
	if(elescale!=kNo) scalecorr = ele->pt - elept;
	out->fillEvent(info,0,pvArr->GetEntriesFast(),scalecorr);
	 
	// events passing selection in this file
	nsel    += weight*kf*npuWgt*trigscale*idscale*embWgt;
	nselvar += weight*weight*kf*kf*npuWgt*npuWgt*trigscale*trigscale*idscale*idscale*embWgt*embWgt;

	// passing events in whole sample 
        nSelEvents += weight*kf*npuWgt*trigscale*idscale*embWgt;

      }

      printf("%8.2f +/- %-8.2f\n",nsel,sqrt(nselvar));

      delete infile;
      infile=0, eventTree=0, nEventTree=0;    
    }
    out->save();
    delete out;

    if(samp->typev.size()>0 && samp->typev[0]==eMC)
      printf("    Yields for %1.2f/fb:",lumi/1000.);
    else
      printf("    Yields for data:    ");

    printf("%10.2f\n",nSelEvents);
    cout << endl;
  }

  delete info;
  delete gen;
  delete muonArr;
  delete electronArr;
  delete jetArr;
  delete pvArr;
  delete svfitArr;
  delete btsf;


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  if(!doNpuRwgt) cout << endl << endl << "Not doing npv reweight!" << endl;
  cout << endl; 
  
  gBenchmark->Show("selectEmu");
}
