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

#include "MitHtt/Common/MitStyleRemix.hh"  // style settings for drawing
#include "MitHtt/Common/CSample.hh"        // helper class for organizing input ntuple files
#include "MitHtt/Common/MyTools.hh"        // miscellaneous helper functions

// define structures to read in ntuple-
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TPFTau.hh"
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

// event-based MVA
#include "MitHtt/Utils/HttMVA.hh" 
#include "MitHtt/Utils/HttElectronMVA.hh" 

#include "MitHtt/Utils/TriggerRatio.h"

#include "Output.hh"

// B-tag scale factors
//#include "BtagSF.hh"

#endif

const Double_t pi = 3.14159265358979;

//double prong3(double pt) {
//  return 0.012+0.001*TMath::Min(TMath::Max(pt-32,0.0),18.0);
//}

//double prong1(double pt) {
//  return 0.025+0.001*TMath::Min(TMath::Max(pt-45,0.0),10.0);
//}

double prong3(double pt) {
  return 0.012;
}

double prong1(double pt) {
  return 0.012;
}


//=== MAIN MACRO =================================================================================================

void selectETau(const TString conf="etau.conf",         // input config file
		const TString outputDir="2012/etau",    // output directory
		const Double_t lumi=1.,        // luminosity pb^-1
		const Int_t is2012=true,          //2012 or 2011 data
		const UInt_t btageff=0,     // b-tag efficiency scale factor uncertainty
		const UInt_t jetunc=0,      // jet energy uncertainties
		const UInt_t mistag=0,      // b mistag rate scale factor uncertainty
		const UInt_t elescale=0     // electron energy scale/resolution uncertainty
) {
  gBenchmark->Start("selectETau");
  
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

  const Double_t kTauPtMin = 15; // 15
  Double_t kElePtMin  = 24;
  if(!is2012) kElePtMin = 20;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 20;

  Bool_t doNpuRwgt = kTRUE;

  // Access samples and fill histograms
  TFile *infile=0;
  TTree *eventTree=0;
  TTree *lTree=0;

  //vbf MVA
  HttMVA *vbfMVA = new HttMVA();
  //vbfMVA->Initialize("BDTG method", getenv("CMSSW_BASE")+std::string("/src/MitHtt/data/VBFMVA/MuTau/VBFMVA_BDTG.weights.xml"), HttMVA::kVBF2);   // vbf mva
  if(!is2012) vbfMVA->Initialize("BDTG method", getenv("CMSSW_BASE")+std::string("/src/MitHtt/data/VBFMVA/MuTau/VBFMVA_BDTG_HCP_42X.weights.xml"), HttMVA::kVBF3);   // vbf mva
  if(is2012)  vbfMVA->Initialize("BDTG method", getenv("CMSSW_BASE")+std::string("/src/MitHtt/data/VBFMVA/MuTau/VBFMVA_BDTG_HCP_52X.weights.xml"), HttMVA::kVBF3);   // vbf mva 

  HttElectronMVA *electronIDMVA = new HttElectronMVA();
  vector<TString> weightFilesEleID;
  weightFilesEleID.push_back(TString::Format("%s/src/MitHtt/data/ElectronMVA/Electrons_BDTG_TrigNoIPV0_2012_Cat1.weights.xml", getenv("CMSSW_BASE")));
  weightFilesEleID.push_back(TString::Format("%s/src/MitHtt/data/ElectronMVA/Electrons_BDTG_TrigNoIPV0_2012_Cat2.weights.xml", getenv("CMSSW_BASE")));
  weightFilesEleID.push_back(TString::Format("%s/src/MitHtt/data/ElectronMVA/Electrons_BDTG_TrigNoIPV0_2012_Cat3.weights.xml", getenv("CMSSW_BASE")));
  weightFilesEleID.push_back(TString::Format("%s/src/MitHtt/data/ElectronMVA/Electrons_BDTG_TrigNoIPV0_2012_Cat4.weights.xml", getenv("CMSSW_BASE")));
  weightFilesEleID.push_back(TString::Format("%s/src/MitHtt/data/ElectronMVA/Electrons_BDTG_TrigNoIPV0_2012_Cat5.weights.xml", getenv("CMSSW_BASE")));
  weightFilesEleID.push_back(TString::Format("%s/src/MitHtt/data/ElectronMVA/Electrons_BDTG_TrigNoIPV0_2012_Cat6.weights.xml", getenv("CMSSW_BASE")));
  electronIDMVA->initialize("BDT", weightFilesEleID);


  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *eleArr      = new TClonesArray("mithep::TElectron");
  TClonesArray *muArr       = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVfit");
  TClonesArray *tauArr      = new TClonesArray("mithep::TPFTau");
  
  Bool_t hasData = (samplev[0]->fnamev.size()>0);
  mithep::TrigEffRatio * tautrigscale = mithep::getTauETrigEffR12();
  mithep::TrigEffRatio * eletrigscale  = mithep::getElectronTrigEffR12();
  setupTrigScale(is2012);
  BtagSF* btsf = new BtagSF(12345);
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
   
      // which corrections to apply where
      Bool_t isemb      = snamev[isam].Contains("emb");
      Bool_t isdata     = !(samp->typev[ifile]==eMC || isemb);
      Bool_t ismadz     = snamev[isam].Contains("ztt-mad"); // madgraph z samples
      Bool_t ismadzee   = snamev[isam].Contains("zee"); 
      Bool_t ismssm     = sfname.Contains("-ggh-") || sfname.Contains("-bbh-");
      Bool_t doIdScale  = !isdata;
      Bool_t doTrigScale= !isdata;
      //Bool_t getGen     =  sfname.Contains("wjets") || doRecoil || reallyDoKf || ismadz ||isemb || ismssm;
      Bool_t doJetUnc   = (jetunc!=kNo);
      Int_t  doRecoil   = (snamev[isam].Contains("ztt") || snamev[isam].Contains("zee")) && !isemb;
      if((snamev[isam].Contains("wjets") || snamev[isam].Contains("w1jets") ||  snamev[isam].Contains("w2jets") || snamev[isam].Contains("w3jets") || snamev[isam].Contains("w4jets") ) && !isemb) doRecoil = 2;
      if((snamev[isam].Contains("_sm_") || snamev[isam].Contains("_mssm_")) && !isemb) doRecoil = 3;
      Bool_t getGen     = sfname.Contains("wjets") || (doRecoil > 0)|| ismadz ||isemb || ismssm;
      Bool_t tauescale = (doRecoil==3) || snamev[isam].Contains("ztt-mad")  || isemb;
      if(!is2012)
	tauescale = 0;
      out->setupRecoil(doRecoil,1,0);
      cout << endl << "Doing recoil correction " << doRecoil << endl;
      cout << sfname << endl;
      cout << snamev[isam] << endl;
      cout << ismadzee << endl;
      cout << ismadz << endl;
      if(ismadzee)
	tauescale=0;
      if(tauescale)
	cout << "Applying tau energy scale" << endl;
      // PU reweighting
      TString pileupReweightFile;
      if(!is2012) {
	cout << "Fall11 sample!" << endl;
	pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_Fall11toFull2011_PixelLumi_50bins.root";
      }  else pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_2012.root";
      TH1F *puWeights = 0;
      TFile *pufile = new TFile(pileupReweightFile.Data());
      TString weightname;
      if(is2012)
	weightname="pileup";
      else
	weightname="puWeights";
      puWeights = (TH1F*)pufile->Get(weightname);

      // setup selecting with JSON file, if necessary
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if((isdata || isemb) && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
	cout << "Embedded " << endl;
        hasJSON = kTRUE;
	ifstream jsonchk; jsonchk.open(samp->jsonv[ifile].Data()); assert(jsonchk.is_open()); jsonchk.close();
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
      lTree =  (TTree*)infile->Get("hEvents"); assert(lTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",  &info);           TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("HPSTau", &tauArr);        TBranch *tauBr      = eventTree->GetBranch("HPSTau");
      eventTree->SetBranchAddress("Electron",     &eleArr);  TBranch *eleBr      = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("Muon",     &muArr);       TBranch *muBr       = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");      
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
      eventTree->SetBranchAddress("SVfitETau", &svfitArr);   TBranch *svfitBr    = eventTree->GetBranch("SVfitETau");
      TBranch *genBr=0;
      if(getGen) {
        eventTree->SetBranchAddress("Gen", &gen);
        genBr = eventTree->GetBranch("Gen");
      }

      // get weights for MC
      Double_t weight=1,treeEntries=-1; // (weight is only initialized for each *file*)
      if(!isdata) {
	treeEntries = (Double_t) lTree->GetEntries();
	assert(treeEntries>0);
	weight = lumi*(samp->xsecv[ifile])/treeEntries;                           // (assumes you've merged filesets)
	if(isemb)  weight=1.0;
      }
      samp->weightv.push_back(weight);
      
      // counters
      Double_t nsel=0, nselvar=0; 
      Double_t nlowmass=0; // low mass z events (below 50)

      cout << eventTree->GetEntries() << " events" << endl;
      int lNEvents = 0;
      // loop over events
      // Start loop
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	if(ientry%100000 == 0) cout << "processing " << float(ientry)/float(eventTree->GetEntriesFast()) << endl;
	//cout << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
	if(getGen)  genBr->GetEntry(ientry);	
	if((fabs(gen->id_1_a) + fabs(gen->id_2_a) == 36 && (fabs(gen->id_1_a) == 17 || fabs(gen->id_2_a) == 17)) && gen->pt_1_b > 20 && gen->pt_2_b > 20)  lNEvents++;

	// skip non-tau events in madgraph sample
	if(ismadz && !ismadzee && (fabs(gen->id_1_a)<15 || fabs(gen->id_1_a)>19)) continue;

        // skip non-mumu events in madgraph sample for zmm
	if(ismadzee && (fabs(gen->id_1_a)>14 && fabs(gen->id_1_a)<20)) continue; 
	infoBr->GetEntry(ientry);

	// certified run selection
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
  
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;
	//trigger
	if(is2012)
	  {
	    //if(!isemb && !(info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20] || info->triggerBits[kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20] || info->triggerBits[kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20] || info->triggerBits[kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20])) continue;
	    if(!isemb && !(info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20] || info->triggerBits[kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20] || info->triggerBits[kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20])) continue;
	    if(isemb && !(info->triggerBits[kHLT_Mu17_Mu8])) continue;
	  }
	else
	  //if(!isemb && !(info->triggerBits[kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15] || info->triggerBits[kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20] || info->triggerBits[kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20] || info->triggerBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] || info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20])) continue;
	  if(!isemb && !(info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15] || info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20] || info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20] || info->triggerBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] || info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20])) continue;
	
        // good primary vertex
        if(!info->hasGoodPV) continue;
	pvArr->Clear();
	pvBr->GetEntry(ientry);	

	double nprong, ngamma,lshift=0;
	// loop through electrons
        vector<const mithep::TElectron*> looseElev;
	vector<const mithep::TElectron*> goodElev;
	eleArr->Clear();
        eleBr->GetEntry(ientry);
	double mvaVal = 0;
	const mithep::TElectron *leadEle = NULL;
        for(Int_t i=0; i<eleArr->GetEntriesFast(); i++) {
	  const mithep::TElectron *ele = (mithep::TElectron*)((*eleArr)[i]);

	  double pIso = eleIsoPU(ele);
	  //if(pIso>0.1) continue;
	
	  Bool_t trigmatch = kFALSE;
	  // trigger matching
	  if(is2012)
	    trigmatch = ((info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20] && ele->hltMatchBits[kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj]) || (info->triggerBits[kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20] && ele->hltMatchBits[kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_EleObj]) ||(info->triggerBits[kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20] && ele->hltMatchBits[kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj]) || (info->triggerBits[kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20] && ele->hltMatchBits[kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_EleObj]) );
	  else
	    //trigmatch = ((info->triggerBits[kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15] && ele->hltMatchBits[kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_EleObj]) || (info->triggerBits[kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20] && ele->hltMatchBits[kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_EleObj]) ||(info->triggerBits[kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20] && ele->hltMatchBits[kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_EleObj]) ||  (info->triggerBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] && ele->hltMatchBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj]) ||(info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] && ele->hltMatchBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj]));
	    trigmatch = ((info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15] && ele->hltMatchBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj]) || (info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20] && ele->hltMatchBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_EleObj]) ||(info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20] && ele->hltMatchBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj]) ||  (info->triggerBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] && ele->hltMatchBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj]) ||(info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] && ele->hltMatchBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj]));
	  
	  if(!isemb && !trigmatch)     continue;
	  if(ele->pt < kElePtMin)		continue;
	  if(fabs(ele->eta) > 2.1)		continue;
	  if(!(pass2012EleMVAID(ele,kTight,1))) continue;
	  //cout << "Old id  " << pass2012EleMVAID(ele,kTight,1) << endl;
	  if(!(eleIsoPU(ele)<0.5)) continue;
	  goodElev.push_back(ele);
	  //double pIso = eleIsoPU(ele);
	  mvaVal = electronIDMVA->mvaValue(ele->fBrem,
					ele->kfChi2,
					ele->kfhits,
					ele->gsfChi2,
					ele->deltaEtaIn,
					ele->deltaPhiIn,
					ele->deltaEtaCalo,
					ele->sigiEtaiEta,
					ele->sigiPhiiPhi,
					ele->etaWidth,
					ele->phiWidth,
					ele->E1x5E5x5,
					ele->R9,				    
					ele->HoverE,
					ele->EoverP,
					(1.0/ele->E - 1.0/ele->p),
					ele->eleEoPout,
					info->rho,
					ele->psOverRaw,				    
					ele->scEta,
					ele->pt,
					0);
	  //cout << "newid " << passTrigNoIPEleMVAID(ele,mvaVal,kTight,1) << " " << mvaVal << endl;
	  //if(!(passTrigNoIPEleMVAID(ele,mvaVal,kTight,1))) continue;
	  if(!leadEle || (ele->pt > leadEle->pt && (pIso < eleIsoPU(leadEle) || pIso < 0.1)) || ( pIso  < eleIsoPU(leadEle) && ( eleIsoPU(leadEle)  > 0.1))) 
	    leadEle = ele;
        }
       	if(!leadEle) continue;
	
        // loop through HPSTaus
        
	tauArr->Clear();
	tauBr->GetEntry(ientry);
	const mithep::TPFTau *leadTau = NULL;
	vector<const mithep::TPFTau*> goodHPSTaus;
	for(Int_t i = 0; i < tauArr->GetEntries(); i++)
	  {
	    const mithep::TPFTau *tau = dynamic_cast<mithep::TPFTau *>(tauArr->At(i));
	    assert(tau);
	    if(toolbox::deltaR(tau->eta, tau->phi, leadEle->eta, leadEle->phi) < 0.3) continue;
	    
	    // Tau ID
	    if(!tauIdElectron(tau)) continue;
	    
	    //if(!(tau->passAntiEleMVA3>2.5)) continue;
	    
	    //if(!(fabs(tau->leadChargedHadronPFCand.dz)<0.2)) continue;
	
	    //bool tt = tau->hpsDiscriminators & mithep::TPFTau::kMVAEle;
	   
	    if(tauescale)
	      {
		if(tau->nSignalPFChargedHadrCands==1 && tau->nSignalPFGammaCands==0)
		  lshift = 0.0;
		else if(tau->nSignalPFChargedHadrCands==1 && tau->nSignalPFGammaCands>0)
		  lshift = prong1(tau->pt);
		else
		  lshift = prong3(tau->pt);
	      }
	    
		// Tau Kinematics
	    if(!((1.0+lshift)*tau->pt > kTauPtMin && fabs(tau->eta) < 2.3)) continue;
	    Bool_t trigmatch = kFALSE;
	    // trigger matching
	    if(is2012)
	      trigmatch = ((info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20] && tau->hltMatchBits[kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj]) || (info->triggerBits[kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20] && tau->hltMatchBits[kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_TauObj]) ||(info->triggerBits[kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20] && tau->hltMatchBits[kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj]) || (info->triggerBits[kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20] && tau->hltMatchBits[kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_TauObj]) );
	    else
	      //trigmatch = ((info->triggerBits[kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15] && tau->hltMatchBits[kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_TauObj]) || (info->triggerBits[kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20] && tau->hltMatchBits[kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_TauObj]) ||(info->triggerBits[kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20] && tau->hltMatchBits[kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_TauObj]) ||  (info->triggerBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] && tau->hltMatchBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj]) ||(info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] && tau->hltMatchBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj]));
	      trigmatch = ((info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15] && tau->hltMatchBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj]) || (info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20] && tau->hltMatchBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_TauObj]) ||(info->triggerBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20] && tau->hltMatchBits[kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj]) ||  (info->triggerBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] && tau->hltMatchBits[kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj]) ||(info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20] && tau->hltMatchBits[kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj]));

	    if(!isemb && !trigmatch)     continue;
	   
	    // Tau Isolation
	    //if(!(tau->ringIso > 0.790)) continue;
	    if(!(tau->rawIso3Hits<10)) continue;
	    
	    //if(!(tau->hcalOverP + tau->ecalOverP > 0.2 ||
	    //	 tau->nSignalPFChargedHadrCands > 1 ||
	    //		 tau->nSignalPFGammaCands > 0)) continue;

	    goodHPSTaus.push_back(tau);
	    if(!leadTau || (tau->pt > leadTau->pt && (tau->rawIso3Hits < leadTau->rawIso3Hits || tau->rawIso3Hits < 1.5)) || (tau->rawIso3Hits < leadTau->rawIso3Hits && (leadTau->rawIso3Hits > 1.5)) )
	      leadTau = tau;
          }	
	if(goodHPSTaus.size()<1) continue;

	Bool_t diElectron = kFALSE;
	for(Int_t i = 0; i < eleArr->GetEntries(); i++)
	  {
	    const mithep::TElectron *ele1 = (mithep::TElectron *)(eleArr->At(i));
	    assert(ele1);
	    
	    for(Int_t j = 0; j < eleArr->GetEntries(); j++)
	      {
		const mithep::TElectron *ele2 = (mithep::TElectron *)(eleArr->At(j));
		assert(ele2);
		
		if(ele1 != ele2 &&
		   ele1->q + ele2->q == 0 &&
		   ele1->pt > 15.0 && ele2->pt > 15.0 &&
		   passEleIdVeto(ele1) && passEleIdVeto(ele2) &&
		   (eleIsoPU(ele1)) < 0.3 && 
		   (eleIsoPU(ele2)) < 0.3 &&
		   toolbox::deltaR(ele1->eta, ele1->phi, ele2->eta, ele2->phi) > 0.15)
		  {
		    diElectron = kTRUE;
		    break;
		  }
	      }
	  }
	//cout << " here ----> " << diElectron << " -- "  <<  leadTau  << " -- " << leadEle  << endl;
	if(!leadTau && !leadEle) continue;
	if(toolbox::deltaR(leadTau->eta, leadTau->phi, leadEle->eta, leadEle->phi) < 0.5) continue;

	if(diElectron) continue;
	bool thirdlep=false;                                                                                                                                                                                                              
	for(Int_t i = 0; i < eleArr->GetEntries(); i++)                                                                                                                                                                                   
	  {                                                                                                                                                                                                                    
	    const mithep::TElectron *ele = (mithep::TElectron *)(eleArr->At(i));                                                                                                                                                           
	    if(toolbox::deltaR(leadEle->eta,leadEle->phi,ele->eta,ele->phi) < 0.3)  continue;                                                                             
	    //if(toolbox::deltaR(leadTau->eta,leadTau->phi,ele->eta,ele->phi) < 0.3)  continue;                                                                                                                                
	    if(ele->pt < 10.0) continue;                                                                                                                                                                                                  
	    if(fabs(ele->eta) > 2.5) continue;                                                                                                                                                                                            
	    if(!(pass2012EleMVAID(ele,kLoose,1))) continue;                                                                                                                                                                               
	    if(eleIsoPU(ele) > 0.3) continue;        

	    mvaVal = electronIDMVA->mvaValue(ele->fBrem,
					     ele->kfChi2,
					     ele->kfhits,
					     ele->gsfChi2,
					     ele->deltaEtaIn,
					     ele->deltaPhiIn,
					     ele->deltaEtaCalo,
					     ele->sigiEtaiEta,
					     ele->sigiPhiiPhi,
					     ele->etaWidth,
					     ele->phiWidth,
					     ele->E1x5E5x5,
					     ele->R9,				    
					     ele->HoverE,
					     ele->EoverP,
					     (1.0/ele->E - 1.0/ele->p),
					     ele->eleEoPout,
					     info->rho,
					     ele->psOverRaw,				    
					     ele->scEta,
					     ele->pt,
					     0);
	    //if(!(passTrigNoIPEleMVAID(ele,mvaVal,kLoose,0))) continue;
	    
	    thirdlep=true;                                                                                                                                                                                                 
	  }                
	muArr->Clear();
        muBr->GetEntry(ientry);                                                                                          
	for(Int_t i=0; i<muArr->GetEntriesFast(); i++) {                                                                                                                                              
	  const mithep::TMuon *muon      = (mithep::TMuon*)     (muArr->At(i));                                                                                                                  
	  //if(toolbox::deltaR(leadEle->eta,leadEle->phi,muon->eta,muon->phi) < 0.3)  continue;                                                                                                             
	  //if(toolbox::deltaR(leadTau->eta,leadTau->phi,muon->eta,muon->phi) < 0.3)  continue;                                                                                                        
	  if(muon->pt < 10.0) continue;                                                                                                                                                              
	  if(fabs(muon->eta) > 2.4) continue;                                                                                                                                                        
	  if(!passTightPFMuonID(muon,1)) continue;
	  if(muonIsoPU(muon) > 0.3) continue;                                                                                                                                                       
	  thirdlep=true;                                                                                                                                                                               
	}                                                                                                                                                                                            
	//thrid lepton veto (WlHhadhad)                                                                                                                                                              
	if(thirdlep) continue;   

	out->fillElectron(leadEle,1,eleIsoPU(leadEle),passEleIsoPU(leadEle,1));
	float passele=0;
	if(passAntiEMVA3(leadTau->antiEleMVA3Cat,leadTau->antiEleMVA3,"Loose"))    passele++;
	if(passAntiEMVA3(leadTau->antiEleMVA3Cat,leadTau->antiEleMVA3,"Medium"))    passele++;
	if(passAntiEMVA3(leadTau->antiEleMVA3Cat,leadTau->antiEleMVA3,"Tight"))    passele++;
	if(passAntiEMVA3(leadTau->antiEleMVA3Cat,leadTau->antiEleMVA3,"VeryTight"))    passele++;
	//if(passele<1.5) continue;
	out->fillTau(leadTau,0,leadTau->ringIso > 0.795,tauescale,passele);
	// SVFit
        svfitArr->Clear();
        svfitBr->GetEntry(ientry);
	bool svf=false;

        for(Int_t i = 0; i < svfitArr->GetEntriesFast(); i++) {
          mithep::TSVfit *svfit = (mithep::TSVfit*) svfitArr->At(i);
          Int_t id = 0;
          if(toolbox::deltaR(leadTau->eta,leadTau->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.1           ) id = 1;
          if(toolbox::deltaR(leadEle->eta,leadEle->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.1 && id == 0) id = 2;
          if(id == 0) continue;
	  if(toolbox::deltaR(leadTau->eta,leadTau->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.1 && id == 2) id = 3;
          if(toolbox::deltaR(leadEle->eta,leadEle->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.1 && id == 1) id = 4;
          if(id < 3)  continue;
	  out->fillCov(svfit);
	  svf = true;
	  //cout << "fill sv  " << endl;
        }
	
 	if(!svf)
 	  {
 	    cout << "No svfit my friend" << endl;
 	    //continue; 
	  }
	
        // loop through jets      
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets = 0, njetsclean=0, nbjets = 0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	out->btagArray.Reset();	out->jptArray.Reset();	out->jetaArray.Reset();	UInt_t npt20jets=0;
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
	  mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	  if(doJetUnc) jet->pt *= (jetunc==kDown) ? (1-jet->unc) : (1+jet->unc);
          if(toolbox::deltaR(jet->eta,jet->phi,leadEle->eta,leadEle->phi) < 0.5) continue;

          if(fabs(jet->eta) > 4.7) continue;
	  //if(!jet->id) continue;
	  if(!passJetIDMVA(jet->pt,jet->eta,jet->mva)) continue;

	  if(jet->pt > kJetPtMin) {
	    njets++;
	  }
	  
	  if(toolbox::deltaR(jet->eta,jet->phi,leadTau->eta,leadTau->phi) < 0.5) continue;

	  // look for b-jets
	  Int_t btagopt = 0;
	  if(isdata||isemb) btagopt = 1;
	  else btagopt = 2;
	  Bool_t btagged = btsf->isbtagged(jet->pt,jet->eta,jet->csv,jet->matchedId,(isdata||isemb),0,0,is2012);
	  if((jet->pt > kBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
	    assert(npt20jets<50);
	    out->btagArray.AddAt(jet->csv,npt20jets);
	    npt20jets++;
	    if(btagged) {
	      nbjets++;
	      if(!bjet || jet->pt > bjet->pt)
		bjet = jet; // leading b-jet
	    }
	  }
	  // look for jets
          if(jet->pt > kJetPtMin) {
            assert(njetsclean<50);
            out->jptArray.AddAt(jet->pt,njetsclean);
            out->jetaArray.AddAt(jet->eta,njetsclean);
  	    njetsclean++;
	    if(!jet1 || jet->pt > jet1->pt) { // jet1 is highest pt, jet2 next highest
	      jet2 = jet1;
	      jet1 = jet;
	    } else if(!jet2 || jet->pt > jet2->pt) {
	      jet2 = jet;
	    }
          }		    
        }

	Int_t nCentralJets=0;
	if(njetsclean>1) {
          for(Int_t i=2; i<jetArr->GetEntriesFast(); i++) {
            mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	    if(toolbox::deltaR(jet->eta,jet->phi,leadTau->eta,leadTau->phi) < 0.5) continue;
	    if(toolbox::deltaR(jet->eta,jet->phi,leadEle->eta,leadEle->phi) < 0.5) continue;
	    //if(!(jet->pt > kJetPtMin && fabs(jet->eta)<4.7 && jet->id==1)) continue;
	    if(!(jet->pt > kJetPtMin && fabs(jet->eta)<4.7 && passJetIDMVA(jet->pt,jet->eta,jet->mva))) continue;
	    if(jet1->eta > jet2->eta && jet->eta > jet2->eta && jet->eta < jet1->eta) nCentralJets++;
	    else if(jet2->eta > jet1->eta && jet->eta > jet1->eta && jet->eta < jet2->eta) nCentralJets++;
	  }
        }

	out->fillJets(jet1,jet2,bjet,0,njetsclean,njets,nbjets,npt20jets,nCentralJets);
        // get k-factor if necessary
        Double_t kf=1;
    
	//W+Jets
	if(doRecoil == 2 && is2012 && gen->npartons == 1) kf  *= 0.203*treeEntries/75276487.;
        if(doRecoil == 2 && is2012 && gen->npartons == 2) kf  *= 0.064*treeEntries/75276487.;
        if(doRecoil == 2 && is2012 && gen->npartons == 3) kf  *= 0.041*treeEntries/75276487.;
        if(doRecoil == 2 && is2012 && gen->npartons == 4) kf  *= 0.039*treeEntries/75276487.;
	
	//W+Jets
	if(doRecoil == 2 && !is2012 && gen->npartons == 1) kf  *= 0.1469*treeEntries/81028892.;
        if(doRecoil == 2 && !is2012 && gen->npartons == 2) kf  *= 0.1415*treeEntries/81028892.;
        if(doRecoil == 2 && !is2012 && gen->npartons == 3) kf  *= 0.1056*treeEntries/81028892.;
        if(doRecoil == 2 && !is2012 && gen->npartons == 4) kf  *= 0.037*treeEntries/81028892.;

	//Z+Jets
	if(doRecoil == 1 && is2012 && gen->npartons == 1) kf *= 0.19275*treeEntries/29960737.;
	if(doRecoil == 1 && is2012 && gen->npartons == 2) kf *= 0.0767*treeEntries/29960737.;
	if(doRecoil == 1 && is2012 && gen->npartons == 3) kf *= 0.047146*treeEntries/29960737.;
	if(doRecoil == 1 && is2012 && gen->npartons == 4) kf *= 0.0371707*treeEntries/29960737.;

	// do vertex reweighting
	Double_t npuWgt = 1;
	if(!isdata && !isemb && doNpuRwgt) {
	  assert(puWeights);
	  Int_t npuxbin = puWeights->GetXaxis()->FindBin(TMath::Min(double(info->nPUTrue), 59.499));
	  npuWgt = puWeights->GetBinContent(npuxbin);
	}

	// lepton ID corrections
	Double_t idscale = 1, embidscale=1;
	
	if(doIdScale) 
	  {
	    if(isemb)
	      embidscale = eleIDScaleETau(leadEle->pt,leadEle->eta,is2012);
	    else
	      idscale = eleIDIsoscaleETau(leadEle->pt,leadEle->eta,is2012);
	  }
	
	if(tauescale)
	  {
	    if(leadTau->nSignalPFChargedHadrCands==1 && leadTau->nSignalPFGammaCands==0)
	      lshift = 0.0;
	    else if(leadTau->nSignalPFChargedHadrCands==1 && leadTau->nSignalPFGammaCands>0)
	      lshift = prong1(leadTau->pt);
	    else
	      lshift = prong3(leadTau->pt);
	  }

		  // trigger scale factor for MC
	Double_t trigscale = 1;
	
	if(doTrigScale && !isemb && !is2012) 
	  trigscale= eff2011TrigEle(leadEle->pt,leadEle->eta,0)*eff2011TrigEleTau((1.0+lshift)*leadTau->pt,leadTau->eta,0);
	if(doTrigScale && isemb && !is2012) 
	  trigscale= eff2011TrigEle(leadEle->pt,leadEle->eta,1)*eff2011TrigEleTau((1.0+lshift)*leadTau->pt,leadTau->eta,1);
	if(doTrigScale && !isemb && is2012) 
	  trigscale = tautrigscale->eff((1.0+lshift)*leadTau->pt,leadTau->eta) * eletrigscale->eff(leadEle->pt,leadEle->eta);
	if(doTrigScale && is2012 && isemb)
	  trigscale = tautrigscale->turnOn((1.0+lshift)*leadTau->pt,leadTau->eta) * eletrigscale->turnOn(leadEle->pt,leadEle->eta);


	if(isnan(trigscale) || isinf(trigscale)) trigscale = 1.;
	
	
	// embedding weight for embedded sample
	Double_t embWgt = 1,embgenWgt=1,embspinWgt=1,embmueffWgt=1,embmuradWgt=1,embkinmass=1;
	Double_t embkinWgt = 1;
	if(!isdata) out->fillGen(gen);
	if(isemb) 
	  {
	    embWgt = info->embGenWeight*info->embSpinnerWeight*info->embMuEffWeight*info->embMuRadWeight;
	    embkinWgt = info->embDiTauMassVsGenDiTauPtRec*info->embGenTau2VsGenTau1PtRec*info->embGenTau2VsGenTau1EtaRec;
	    embgenWgt = info->embGenWeight;
	    embspinWgt = info->embSpinnerWeight;
	    embmueffWgt = info->embMuEffWeight;
	    embmuradWgt = info->embMuRadWeight;
	    embkinmass = info->embDiTauMassVsGenDiTauPtRec;
	  }
	out->fEmbWeight = embWgt;
	out->fEmbKinWeight = embkinWgt;
	out->fEmbIdScale = embidscale;
	out->fEmbGenWeight = embgenWgt;
	out->fEmbSpinnerWeight = embspinWgt;
	out->fEmbMuEff = embmueffWgt;
	out->fEmbMuRad =   embmuradWgt;
	out->fEmbKinMassWeight = embkinmass;
	out->fMCWeight	 = weight*kf*embWgt/lumi;
	out->fPUWeight	 = npuWgt;
	out->fEffWeight	 = trigscale*idscale;
	out->fWeight	 = weight*kf*npuWgt*trigscale*idscale*embWgt/lumi;
	out->fillEvent(info,vbfMVA,pvArr->GetEntriesFast());

	// events passing selection in this file
	nsel    += weight*kf*npuWgt*trigscale*idscale*embWgt;
	nselvar += weight*weight*kf*kf*npuWgt*npuWgt*trigscale*trigscale*idscale*idscale*embWgt*embWgt;
	if(doRecoil && (gen->vmass_a < 50)) nlowmass += weight*kf*npuWgt*trigscale*idscale*embWgt;

	// passing events in whole sample 
        nSelEvents += weight*kf*npuWgt*trigscale*idscale*embWgt;
      }
      cout << " total good events " << lNEvents << endl;
      printf("%8.2f +/- %-8.2f\n",nsel,sqrt(nselvar));
      if(nlowmass > 0) printf("           ---> selected events with z mass < 50:  %10.3f (out of %15.3f)\n",nlowmass,nsel);

      delete infile;
      infile=0, eventTree=0, lTree = 0;    
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
  delete eleArr;
  delete tauArr;
  delete jetArr;
  delete pvArr;
  delete svfitArr;
  delete vbfMVA;
  

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  if(!doNpuRwgt) cout << endl << endl << "Not doing npv reweight!" << endl;
  cout << endl; 
  
  gBenchmark->Show("selectTau");
}


