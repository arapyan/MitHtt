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

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TPFTau.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"   
#include "MitHtt/Ntupler/interface/TVertex.hh"   
#include "MitHtt/Ntupler/interface/TSVfit.h"
#include "MitHtt/Ntupler/interface/TSVfitter.h"
#include "MitHtt/Ntupler/interface/TMuon.hh" 
#include "MitHtt/Ntupler/interface/TElectron.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh"

// scale factros
#include "MitHtt/Utils/DataMC.hh"

#include "Output.hh"

// B-tag scale factors
#include "BtagSF.hh"

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

void selectTauTau(const TString conf="tautau.conf",  // input config file
		  const TString outputDir="tmp",    // output directory
		  const Double_t lumi=1.,        // luminosity pb^-1
		  const Int_t is2012=true          //2012 or 2011 data
		  ) {
  gBenchmark->Start("selectTauTau");
  
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

  enum {eMC};  // dataset type  

  const Double_t kTauPtMin = 40; //30;
  const Double_t kJetPtMin   = 30; //30.0;
  const Double_t kBJetPtMin  = 20; //20;
  
  Bool_t doKFactors = kTRUE;
  if(is2012)
    doKFactors = kFALSE;      // not needed in Summer12
  
  Bool_t doNpuRwgt = kTRUE;

  // Access samples and fill histograms
  TFile *infile=0;
  TTree *eventTree=0;  
  TTree *lTree=0;

  BtagSF* btsf = new BtagSF(12345);
 
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVfit");
  TClonesArray *tauArr      = new TClonesArray("mithep::TPFTau");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  
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
  
      // which corrections to apply where
      Bool_t isemb      = snamev[isam].Contains("emb");
      Bool_t isdata     = !(samp->typev[ifile]==eMC || isemb);
      Bool_t reallyDoKf = doKFactors && sfname.Contains("-gf-");
      Bool_t ismadz     = snamev[isam].Contains("ztt-mad"); // madgraph z samples
      Bool_t ismadzmm   =  snamev[isam].Contains("zll"); // madgraph z samples 
      Bool_t ismssm     = sfname.Contains("-ggh-") || sfname.Contains("-bbh-");
      Bool_t doIdScale  = !isdata;
      Bool_t doTrigScale= !isdata;
      Int_t  doRecoil   = (snamev[isam].Contains("ztt") || snamev[isam].Contains("zll") || sfname.Contains("zjets")) && !isemb;
      if(snamev[isam].Contains("wjets") || snamev[isam].Contains("w1jets") ||  snamev[isam].Contains("w2jets") || snamev[isam].Contains("w3jets") || snamev[isam].Contains("w4jets")) doRecoil = 2;
      if(snamev[isam].Contains("_sm_") || snamev[isam].Contains("_mssm_")) doRecoil = 3;
      Bool_t getGen     = doRecoil || ismadz ||isemb || ismssm || ismadzmm;
      Bool_t tauescale = (doRecoil==3) || sfname.Contains("zllm")  ||  sfname.Contains("z1jetsm50") || sfname.Contains("z2jetsm50") || sfname.Contains("z3jetsm50") || sfname.Contains("z4jetsm50") ||  isemb;
      //if(sfname.Contains("vtth")) getGen=0;
   
      //tauescale=0;
      cout << sfname << endl;
      cout << snamev[isam] << endl;
      cout << ismadzmm << endl;
      cout << ismadz << endl;
      cout << isemb << endl;
      out->setupRecoil(doRecoil,1,0);
      if(doRecoil)
	cout << "Doing recoil " << doRecoil << endl;
      if(ismadzmm)
	tauescale=0;
      if(tauescale)
	cout << "Applying tau energy scale" << endl;
      if(tauescale) cout << "Applying tau energy scale corrections" << endl;
     
     
      // PU reweighting
      TString pileupReweightFile;
      if(!is2012) {
	cout << "Fall11 sample!" << endl;
	pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_Fall11toFull2011_PixelLumi_50bins.root";
      } //else pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_S1253XTo2012_12ifb.root";
      else 
	pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_2012.root";
      TH1F *puWeights = 0;
      TFile *pufile = new TFile(pileupReweightFile.Data());
      puWeights = (TH1F*)pufile->Get("pileup");
     
      // setup selecting with JSON file, if necessary
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
     
      if((isdata || isemb) && (samp->jsonv[ifile].CompareTo("NONE")!=0)) {
	cout << "embedded" << endl;
        hasJSON = kTRUE;
	cout << endl;
	cout << samp->jsonv[ifile] << endl;
	ifstream jsonchk; jsonchk.open(samp->jsonv[ifile].Data()); assert(jsonchk.is_open()); jsonchk.close();
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }
     

      // k-factors
      TH1D *hKFactors = (reallyDoKf) ? kfFHPInit(higgsmass(sfname)) : 0;

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
      lTree =  (TTree*)infile->Get("hEvents"); assert(lTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",  &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("HPSTau", &tauArr);   TBranch *tauBr = eventTree->GetBranch("HPSTau");
      eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");      
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
      eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("SVfitTauTau", &svfitArr);    TBranch *svfitBr    = eventTree->GetBranch("SVfitTauTau");
      TBranch *genBr=0;
      if(getGen) {
        eventTree->SetBranchAddress("Gen", &gen);
        genBr = eventTree->GetBranch("Gen");
      }

      // get weights for MC
      Double_t weight=1,treeEntries=-1; // (weight is only initialized for each *file*)
      if(!isdata) {
	treeEntries = (Double_t)lTree->GetEntries();
	assert(treeEntries>0);
        weight = lumi*(samp->xsecv[ifile])/treeEntries;                           // (assumes you've merged filesets)
	if(isemb)  weight=1.0;
      }
      samp->weightv.push_back(weight);

      // counters
      Double_t nsel=0, nselvar=0; 
      Double_t nlowmass=0; // low mass z events (below 50)
      int ev[4]= {3116902,2863619,1430290,3263514};
      cout << eventTree->GetEntries() << " events" << endl;
   
      // loop over events
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	if(ientry%100000 == 0) cout << "processing " << float(ientry)/float(eventTree->GetEntriesFast()) << endl;
	if(getGen)  genBr->GetEntry(ientry);
        infoBr->GetEntry(ientry);
	
	if(ismadz && !ismadzmm && (fabs(gen->id_1_a)<15 || fabs(gen->id_1_a)>19)) continue;
	
        // skip non-mumu events in madgraph sample for zmm
        if(ismadzmm && (fabs(gen->id_1_a)>14 && fabs(gen->id_1_a)<20)) continue;

	//cout << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;

	// certified run selection
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;

	// trigger
 	//if(!isemb  && !(info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1] || info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1])) continue;
	if(!isemb  && !(info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1] || info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1])) continue;
	//if(isemb && !(info->triggerBits[kHLT_Mu17_Mu8])) continue;
	
	//if(!isemb  && !(info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30])) continue;
	
        // good primary vertex
        if(!info->hasGoodPV) continue;
	pvArr->Clear();
	pvBr->GetEntry(ientry);	
	
	double nprong, ngamma,lshift=0;

        // loop through HPSTaus
        
	tauArr->Clear();
	tauBr->GetEntry(ientry);

	const mithep::TPFTau *leadTau = NULL;
	const mithep::TPFTau *subTau  = NULL;
	for(Int_t i = 0; i < tauArr->GetEntries(); i++)
	  {
	    const mithep::TPFTau *tau = dynamic_cast<mithep::TPFTau *>(tauArr->At(i));
	    assert(tau);
        
	    // Tau ID
	    if(!passtauId(tau)) continue;
	    
	    if(tauescale)
	      {
		if(tau->nSignalPFChargedHadrCands==1 && tau->nSignalPFGammaCands==0)
		  lshift = 0.00;
		else if(tau->nSignalPFChargedHadrCands==1 && tau->nSignalPFGammaCands>0)
		  lshift = prong1(tau->pt);
		else
		  lshift = prong3(tau->pt);
	      }
	    
	    // Tau Kinematics
	    if(!((1.0+lshift)*tau->pt > kTauPtMin && fabs(tau->eta) < 2.1)) continue;
	  
	    //Tau HLT
	    //Bool_t trigmatch = ((info->triggerBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1Obj]));
	    //Bool_t trigmatch = ((info->triggerBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj]));
	    Bool_t trigmatch = ((info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1Obj]));
	   
	    if(!isemb && !trigmatch)     continue;
	   
	   
	    // Tau Isolation
	    //if(!(tau->ringIso > 0.5)) continue; //0.921 tight
	    if(!(tau->rawIso3Hits < 1.0)) continue;
	    // if(!(tau->hcalOverP + tau->ecalOverP > 0.2 ||
// 		 tau->nSignalPFChargedHadrCands > 1 ||
// 	    	 tau->nSignalPFGammaCands > 0)) continue;
	 
	    if(!leadTau || tau->rawIso3Hits < leadTau->rawIso3Hits)
	      {
		subTau = leadTau;
		leadTau = tau;
	      } else if(!subTau || tau->rawIso3Hits < subTau->rawIso3Hits) {
	      subTau = tau;
	    }
	  }
	

 	if(!(leadTau && subTau)) continue;

	double lshift1=0,lshift2=0;
 
	if(tauescale)
	  {
	    if(leadTau->nSignalPFChargedHadrCands==1 && leadTau->nSignalPFGammaCands==0)
	      lshift1 = 0.00;
	    else if(leadTau->nSignalPFChargedHadrCands==1 && leadTau->nSignalPFGammaCands>0)
	      lshift1 = prong1(leadTau->pt);
	    else
	      lshift1 = prong3(leadTau->pt);
	    
	    if(subTau->nSignalPFChargedHadrCands==1 && subTau->nSignalPFGammaCands==0)
	      lshift2 = 0.00;
	    else if(subTau->nSignalPFChargedHadrCands==1 && subTau->nSignalPFGammaCands>0)
	      lshift2 = prong1(subTau->pt);
	    else
	      lshift2 = prong3(subTau->pt);
	  }
	
	float passele=0, passele2=0;
	if(passAntiEMVA3(leadTau->antiEleMVA3Cat,leadTau->antiEleMVA3,"Loose"))    passele++;
	if(passAntiEMVA3(leadTau->antiEleMVA3Cat,leadTau->antiEleMVA3,"Medium"))    passele++;
	if(passAntiEMVA3(leadTau->antiEleMVA3Cat,leadTau->antiEleMVA3,"Tight"))    passele++;
	if(passAntiEMVA3(leadTau->antiEleMVA3Cat,leadTau->antiEleMVA3,"VeryTight"))    passele++;

	if(passAntiEMVA3(subTau->antiEleMVA3Cat,subTau->antiEleMVA3,"Loose"))    passele2++;
	if(passAntiEMVA3(subTau->antiEleMVA3Cat,subTau->antiEleMVA3,"Medium"))    passele2++;
	if(passAntiEMVA3(subTau->antiEleMVA3Cat,subTau->antiEleMVA3,"Tight"))    passele2++;
	if(passAntiEMVA3(subTau->antiEleMVA3Cat,subTau->antiEleMVA3,"VeryTight"))    passele2++;

	if((1.0+lshift2)*subTau->pt < (1.0+lshift1)*leadTau->pt)
	  {
	    out->fillTau(leadTau,1,leadTau->ringIso > 0.884,tauescale,passele);
	    out->fillTau(subTau,0,subTau->ringIso > 0.884,tauescale,passele2);
	  }
	else
	  {
	    out->fillTau(leadTau,0,leadTau->ringIso > 0.884,tauescale,passele);
	    out->fillTau(subTau,1,subTau->ringIso > 0.884,tauescale,passele2);
	  }
	

	muonArr->Clear();
	muonBr->GetEntry(ientry);
	electronArr->Clear();
	electronBr->GetEntry(ientry);
	
	bool thirdlep=false;
	for(Int_t i = 0; i < electronArr->GetEntries(); i++)
	  {
	    const mithep::TElectron *ele = (mithep::TElectron *)(electronArr->At(i));
	    assert(ele);
	    if(ele->pt < 1.0) continue;
	    if(fabs(ele->eta) > 2.5) continue;
	    if(!(pass2012EleMVAID(ele,kLoose,1))) continue;
	    if(eleIsoPU(ele) > 0.3) continue;
	    thirdlep=true;
	  }
	for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const mithep::TMuon *muon = (mithep::TMuon*)((*muonArr)[i]);
	  assert(muon);
	  if(muon->pt < 10.0) continue;
	  if(fabs(muon->eta) > 2.4) continue;
	  if(!passTightPFMuonID(muon,1)) continue;
	  if(muonIsoPU(muon) > 0.3) continue;
	  thirdlep=true;
	}
	//thrid lepton veto (WlHhadhad)       
	if(thirdlep) continue;
	
	// SVFit
        svfitArr->Clear();
        svfitBr->GetEntry(ientry);

        for(Int_t i = 0; i < svfitArr->GetEntriesFast(); i++) {
          mithep::TSVfit *svfit = (mithep::TSVfit*) svfitArr->At(i);
          Int_t id = 0;
          if(toolbox::deltaR(leadTau->eta,leadTau->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.01           ) id = 1;
          if(toolbox::deltaR(subTau->eta,subTau->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi())   < 0.01 && id == 0) id = 2;
          if(id == 0) continue;
          if(toolbox::deltaR(leadTau->eta,leadTau->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.01 && id == 2) id = 3;
          if(toolbox::deltaR(subTau->eta,subTau->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi())   < 0.01 && id == 1) id = 4;
          if(id < 3) continue;
          out->fillCov(svfit);
        }
     
        // loop through jets      
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets = 0, njetsclean=0, nbjets = 0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet1=0, *bjet2=0;
	out->btagArray.Reset();	out->jptArray.Reset();	out->jetaArray.Reset();	UInt_t npt20jets=0;
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
	  mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	  
	  if(toolbox::deltaR(jet->eta,jet->phi,leadTau->eta,leadTau->phi) < 0.5) continue;
	  
	  if(fabs(jet->eta) > 4.7) continue;
	  //if(!jet->id) continue;
	  if(!passJetIDMVA(jet->pt,jet->eta,jet->mva)) continue;
	  if(jet->pt > kJetPtMin) {
	    njets++;
	  }
	  	  
 
          if(toolbox::deltaR(jet->eta,jet->phi,subTau->eta,subTau->phi) < 0.5) continue;

         
	  // look for b-jets
	  Int_t btagopt = 0;
	  if(isdata||isemb) btagopt = 1;
	  else btagopt = 2;
	  //Bool_t btagged = isbtagged(jet,btagopt,0,0);
	  Bool_t btagged = btsf->isbtagged(jet->pt,jet->eta,jet->csv,jet->matchedId,(isdata||isemb),0,0,is2012);
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

	// trigger scale factor for MC
	Double_t trigscale = 1;

// 	if(!isemb && (((info->triggerBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30] && leadTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj] && subTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30] && leadTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj] && subTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30] && leadTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj] && subTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj])) && !((info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1] && leadTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1Obj] && subTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau35_Trk1_eta2p1Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1] && leadTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1Obj] && subTau->hltMatchBits[kHLT_DoubleMediumIsoPFTau35_Trk5_eta2p1Obj]))))
// 	{
// 	    if(!jet1) continue;
// 	    if(!(jet1->hltMatchBits[kHLT_hltTripleL2Jets30eta3Obj])) continue;
// 	    if(!(jet1->pt>50.0)) continue;
// 	    if(!(fabs(jet1->eta)<3.0)) continue;
// 	    if(doTrigScale)
// 	      trigscale =  (eff2012IsoTau19fbData((1.0+lshift1)*leadTau->pt,leadTau->eta)/eff2012IsoTau19fbMC((1.0+lshift1)*leadTau->pt,leadTau->eta))*(eff2012IsoTau19fbData((1.0+lshift2)*subTau->pt,subTau->eta)/eff2012IsoTau19fbMC((1.0+lshift2)*subTau->pt,subTau->eta));
// 	}
// 	else
// 	  {
// 	    if(doTrigScale)
// 	      trigscale =  (eff2012IsoParkedTau19fbDataMSSM((1.0+lshift1)*leadTau->pt,leadTau->eta)/eff2012IsoParkedTau19fbMCMSSM((1.0+lshift1)*leadTau->pt,leadTau->eta))*(eff2012IsoParkedTau19fbDataMSSM((1.0+lshift2)*subTau->pt,subTau->eta)/eff2012IsoParkedTau19fbMCMSSM((1.0+lshift2)*subTau->pt,subTau->eta));
// 	    trigscale =  (eff2012IsoParkedTau19fbData((1.0+lshift1)*leadTau->pt,leadTau->eta)/eff2012IsoParkedTau19fbMC((1.0+lshift1)*leadTau->pt,leadTau->eta))*(eff2012IsoParkedTau19fbData((1.0+lshift2)*subTau->pt,subTau->eta)/eff2012IsoParkedTau19fbMC((1.0+lshift2)*subTau->pt,subTau->eta));
// 	  }
	
 // 	if(doTrigScale && !isemb)
// 	{
// 	  trigscale =  (eff2012IsoParkedTau19fbData((1.0+lshift1)*leadTau->pt,leadTau->eta)/eff2012IsoParkedTau19fbMC((1.0+lshift1)*leadTau->pt,leadTau->eta))*(eff2012IsoParkedTau19fbData((1.0+lshift2)*subTau->pt,subTau->eta)/eff2012IsoParkedTau19fbMC((1.0+lshift2)*subTau->pt,subTau->eta));
// 	}

	if(doTrigScale && !isemb)
	{
	  trigscale =  (eff2012IsoParkedTau19fbDataMSSM((1.0+lshift1)*leadTau->pt,leadTau->eta)/eff2012IsoParkedTau19fbMCMSSM((1.0+lshift1)*leadTau->pt,leadTau->eta))*(eff2012IsoParkedTau19fbDataMSSM((1.0+lshift2)*subTau->pt,subTau->eta)/eff2012IsoParkedTau19fbMCMSSM((1.0+lshift2)*subTau->pt,subTau->eta));
	}
	
 	if(doTrigScale && isemb)
	trigscale =  eff2012IsoParkedTau19fbData((1.0+lshift1)*leadTau->pt,leadTau->eta)*eff2012IsoParkedTau19fbData((1.0+lshift2)*subTau->pt,subTau->eta);
					       
	Int_t nCentralJets=0;
	if(njetsclean>1) {
          for(Int_t i=2; i<jetArr->GetEntriesFast(); i++) {
            mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	    if(toolbox::deltaR(jet->eta,jet->phi,leadTau->eta,leadTau->phi) < 0.5) continue;
	    if(toolbox::deltaR(jet->eta,jet->phi,subTau->eta,subTau->phi) < 0.5) continue;
	    //if(!(jet->pt > kJetPtMin && fabs(jet->eta)<4.7 && jet->id==1)) continue;
	    if(!(jet->pt > kJetPtMin && fabs(jet->eta)<4.7 && passJetIDMVA(jet->pt,jet->eta,jet->mva))) continue;
	    if(jet1->eta > jet2->eta && jet->eta > jet2->eta && jet->eta < jet1->eta) nCentralJets++;
	    else if(jet2->eta > jet1->eta && jet->eta > jet1->eta && jet->eta < jet2->eta) nCentralJets++;
	  }
        }
	
	out->fillJets(jet1,jet2,bjet1,bjet2,njetsclean,njets,nbjets,npt20jets,nCentralJets);

        // get k-factor if necessary
        Double_t kf=1;
        if(reallyDoKf) kf = kfFHPValue(gen->vpt_a, hKFactors);

	// do vertex reweighting
	Double_t npuWgt = 1;
	if(!isdata && !isemb && doNpuRwgt) {
	  assert(puWeights);
	  Int_t npuxbin = puWeights->GetXaxis()->FindFixBin(TMath::Min(double(info->nPUTrue), 59.499));
	  npuWgt = puWeights->GetBinContent(npuxbin);
	}
		
// 	//W+Jets
// 	if(doRecoil == 2 && is2012 && gen->npartons == 1) kf  *= 0.308*treeEntries/56883397.;
//         if(doRecoil == 2 && is2012 && gen->npartons == 2) kf  *= 0.090*treeEntries/56883397.;
//         if(doRecoil == 2 && is2012 && gen->npartons == 3) kf  *= 0.061*treeEntries/56883397.;
//         if(doRecoil == 2 && is2012 && gen->npartons == 4) kf  *= 0.030*treeEntries/56883397.;

	
	//W+Jets
	if(doRecoil == 2 && is2012 && gen->npartons == 1) kf  *= 0.203*treeEntries/75276487.;
        if(doRecoil == 2 && is2012 && gen->npartons == 2) kf  *= 0.064*treeEntries/75276487.;
        if(doRecoil == 2 && is2012 && gen->npartons == 3) kf  *= 0.041*treeEntries/75276487.;
        if(doRecoil == 2 && is2012 && gen->npartons == 4) kf  *= 0.039*treeEntries/75276487.;

	//Z+Jets
	if(doRecoil == 1 && is2012 && gen->npartons == 1) kf *= 0.19275*treeEntries/29960737.;
	if(doRecoil == 1 && is2012 && gen->npartons == 2) kf *= 0.0767*treeEntries/29960737.;
	if(doRecoil == 1 && is2012 && gen->npartons == 3) kf *= 0.047146*treeEntries/29960737.;
	if(doRecoil == 1 && is2012 && gen->npartons == 4) kf *= 0.03572*treeEntries/29960737.;


	// lepton ID corrections
	Double_t idscale = 1;
	
	// embedding weight for embedded sample
	//Double_t embWgt = 1;
	Double_t embWgt = 1,embgenWgt=1,embspinWgt=1,embmueffWgt=1,embmuradWgt=1,embkinmass=1;
	Double_t embkinWgt = 1;
      	if(!isdata) out->fillGen(gen);
	if(isemb) 
	  {
	    embWgt = info->embGenWeight*info->embSpinnerWeight*info->embMuEffWeight*info->embMuRadWeight;
	    //embkinWgt = info->embDiTauMassVsGenDiTauPtRec*info->embGenTau2VsGenTau1PtRec*info->embGenTau2VsGenTau1EtaRec;
	    embgenWgt = info->embGenWeight;
	    embspinWgt = info->embSpinnerWeight;
	    embmueffWgt = info->embMuEffWeight;
	    embmuradWgt = info->embMuRadWeight;
	    embkinmass = info->embDiTauMassVsGenDiTauPtRec;
	  }
	out->fEmbWeight = embWgt;
	out->fEmbKinWeight = embkinWgt;
	out->fEmbIdScale = idscale;
	out->fEmbGenWeight = embgenWgt;
	out->fEmbSpinnerWeight = embspinWgt;
	out->fEmbMuEff = embmueffWgt;
	out->fEmbMuRad =   embmuradWgt;
	out->fEmbKinMassWeight = embkinmass;
	out->fMCWeight	 = weight*kf*embWgt/lumi;
	out->fPUWeight	 = npuWgt;
	out->fEffWeight	 = trigscale*idscale;
	out->fWeight	 = weight*kf*npuWgt*trigscale*idscale*embWgt/lumi;
 	out->fillEvent(info,0,pvArr->GetEntriesFast());
        // events passing selection in this file
	nsel    += weight*kf*npuWgt*trigscale*idscale*embWgt;
	nselvar += weight*weight*kf*kf*npuWgt*npuWgt*trigscale*trigscale*idscale*idscale*embWgt*embWgt;
	if(doRecoil && (gen->vmass_a < 50)) nlowmass += weight*kf*npuWgt*trigscale*idscale*embWgt;

	// passing events in whole sample 
        nSelEvents += weight*kf*npuWgt*trigscale*idscale*embWgt;
      }

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
  delete jetArr;
  delete pvArr;
  delete tauArr;
  delete svfitArr;
  delete btsf;
  delete electronArr;
  delete muonArr;

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  if(!doNpuRwgt) cout << endl << endl << "Not doing npv reweight!" << endl;
  cout << endl; 
  gBenchmark->Show("selectTauTau");
}

