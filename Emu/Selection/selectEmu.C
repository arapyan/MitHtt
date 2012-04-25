#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TH1.h>                    // histogram base class
#include <TH2.h>                    // 2d histogram base class
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

// recoil corrections
#include "MitHtt/Utils/RecoilCorrector.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"

// lepton ID helper functions
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitHtt/Utils/LeptonIDCuts.hh"

// define structure for output ntuple
#include "EmuData.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================
const Double_t pi = 3.14159265358979;
TRandom1 randm(0xDEADBEEF);
enum { kNo, kDown, kUp };                     // systematic variations 

// Initialize k-factors
TH1D* kfFHPInit(Int_t mH);

// Get k-factor
Double_t kfFHPValue(Double_t pt, TH1D* hKF);

// Is jet b-tagged?
Bool_t isbtagged(mithep::TJet *jet, Int_t isdata, UInt_t btageff, UInt_t mistag);

// Get higgs mass point from sample name
Int_t higgsmass(TString basename)
{
  stringstream ss(basename(TRegexp("[0-9][0-9]*"),3).Data());
  Int_t mass;
  ss >> mass;
  if(basename.Contains("-gf-")) assert(mass>85 && mass<1200);
  return mass;
}

// Get unfolding weights for embedded
Double_t embUnfoldWgt(Double_t pt1, Double_t eta1, Double_t pt2, Double_t eta2);

// Lepton id scale factors
Double_t eleIDscale(Double_t elept, Double_t eleeta);
Double_t muIDscale(Double_t mupt, Double_t mueta);

// Trigger scale factors/efficiencies
Double_t eleTrigScale(Double_t elept, Double_t eleeta);
Double_t muTrigScale(Double_t mupt, Double_t mueta);
Double_t eleTrigEff(Double_t elept, Double_t eleeta);
Double_t muTrigEff(Double_t mupt, Double_t mueta);

// Get number of entries in unskimmed tree
Double_t unskimmedEntries(TString skimname);

//=== MAIN MACRO =================================================================================================

void selectEmu(const TString conf,         // input config file
               const TString outputDir,    // output directory
	       const Double_t lumi,        // luminosity pb^-1
               const UInt_t btageff=0,     // b-tag efficiency scale factor uncertainty
               const UInt_t jetunc=0,      // jet energy uncertainties
               const UInt_t mistag=0,      // b mistag rate scale factor uncertainty
	       const UInt_t elescale=0     // electron energy scale/resolution uncertainty
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
      // fake-rate fakes come from a separate macro
      if((TString(line).Contains("fake")) && !(TString(line).Contains("ss"))) continue;
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
      if(TString(fname).Contains("dummy",TString::kIgnoreCase)) continue;
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
  enum { kMuMu, kEleEle, kEleMu, kMuEle };      // final state type
  
  const Double_t kMuonPt1Min = 20;
  const Double_t kMuonPt2Min = 10;
  
  const Double_t kElePt1Min  = 20;
  const Double_t kElePt2Min  = 10;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 15;
  
  Bool_t doKFactors = kTRUE;

  Bool_t doNpuRwgt = kTRUE;

  mithep::ElectronIDMVA *electronIDMVANoIPInfo = new mithep::ElectronIDMVA();
  electronIDMVANoIPInfo->Initialize("BDTG method",
                                      TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                                      TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                                      TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                                      TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                                      TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                                      TString::Format("%s/src/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml", getenv("CMSSW_BASE")),
                                      mithep::ElectronIDMVA::kNoIPInfo );

  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0;  
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVfit");

  Bool_t hasData = (samplev[0]->fnamev.size()>0);

  //
  // loop over samples
  //
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if(isam==0 && !hasData) continue;
  
    CSample* samp = samplev[isam];

    Double_t nSelEvents[3]; for(Int_t i=0; i<3; i++) nSelEvents[i]=0; // events in this sample
	  
    //
    // Set up output ntuple file for the sample
    //
    TString outfname = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    TFile outfile(outfname,"RECREATE");
    TTree outtree("Events","Events");

    EmuData data;
    Double_t rawMet,rawprojvar,npuWgt, rawMetphi;
    UInt_t npt15jets;
    const UInt_t kMaxPt15Jets=50;
    TArrayF btagArray; btagArray.Set(kMaxPt15Jets); // array to hold b-tag values for pt-15 jets
    TArrayF jptArray; jptArray.Set(kMaxPt15Jets);   // array to hold jet pt values
    TArrayF jetaArray; jetaArray.Set(kMaxPt15Jets); // array to hold jet eta values
    outtree.Branch("Events",&data.runNum,
"runNum/i:evtNum:lumiSec:nPV:njets:nbjets:vpt/F:vphi:rawmet:rawmetphi:met:metphi:mass:dphi:mt:pt:phi:pmet:pvis:eleiso:muiso:eled0:eled0sig:eleip3d:eleip3dsig:mud0:mud0sig:muip3d:muip3dsig:lpt1:leta1:lphi1:lpt2:leta2:lphi2:jpt1:jeta1:jphi1:jpt2:jeta2:jphi2:bjpt:bjeta:bjphi:mjj:svfmass:svfmassunc:genlpt1:genleta1:genlphi1:genlpt2:genleta2:genlphi2:weight:state/I");

    // extra branches
    outtree.Branch("npt15jets",&npt15jets);
    outtree.Branch("btagArray",&btagArray);
    outtree.Branch("jptArray",&jptArray);
    outtree.Branch("jetaArray",&jetaArray);
    outtree.Branch("rawMet",&rawMet);
    outtree.Branch("rawprojvar",&rawprojvar);
    outtree.Branch("npuWgt",&npuWgt);

    //
    // loop through files
    //
    cout <<  "processing " << snamev[isam] << ":" << endl;
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      printf("        %-55s",(samp->fnamev[ifile]+"...").Data()); fflush(stdout);
      infile = new TFile(samp->fnamev[ifile]); 
      assert(infile);

      TString sfname    = samp->fnamev[ifile];
      TString basename = sfname(sfname.Last('/')+1,sfname.Last('.') - sfname.Last('/') - 1);

      //
      // which corrections to apply where
      //
      Bool_t isdata     = !(samp->typev[ifile]==eMC);
      Bool_t is42mc     = sfname.Contains("s11-") && !sfname.Contains("f11");
      Bool_t isemb      = snamev[isam].Contains("emb");
      Bool_t isfall11   = sfname.Contains("f11");
      Bool_t issamesign = snamev[isam].Contains("ss-fakes");
      Bool_t hasTrigs   = isdata || is42mc;
      Bool_t doRecoil   = (sfname.Contains("ztt") || sfname.Contains("-zll") || sfname.Contains("zjets") || snamev[isam].Contains("_sm_") || snamev[isam].Contains("_mssm_"));
      Bool_t reallyDoKf = doKFactors && sfname.Contains("-gf-");
      Bool_t ismadz     = sfname.Contains("-zll") || sfname.Contains("-zjets"); // madgraph z samples
      Bool_t ismadzmm   = snamev[isam].Contains("zmm") && (sfname.Contains("-zll") || sfname.Contains("-zjets")); // madgraph z samples
      Bool_t doIdScale  = !isdata;
      Bool_t doTrigScale= !isdata;
      Bool_t getGen     = doRecoil || reallyDoKf || ismadz || isemb;
      Bool_t doJetUnc   = (jetunc!=kNo) && (isdata || is42mc);

      // PU reweighting
      TString pileupReweightFile;
      if(sfname.Contains("f11")) {
	pileupReweightFile = "/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root";
      } else pileupReweightFile = "/data/smurf/sixie/Pileup/weights/PileupReweighting.Summer11DYmm_To_Full2011.root";
      TH1F *puWeights = 0;
      TFile *pufile = new TFile(pileupReweightFile.Data());
      puWeights = (TH1F*)pufile->Get("puWeights");

      // setup selecting with JSON file, if necessary
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if(isdata && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	ifstream jsonchk; jsonchk.open(samp->jsonv[ifile].Data()); assert(jsonchk.is_open()); jsonchk.close();
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }

      RecoilCorrector *corrector=0;
      if(doRecoil) {
        corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_datamm_njet.root");
      }

      mithep::TSVfitter *fitter = new mithep::TSVfitter();


      TH1D *hKFactors = (reallyDoKf) ? kfFHPInit(higgsmass(basename)) : 0;

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

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
	else                              treeEntries = (Double_t)eventTree->GetEntries();
	assert(treeEntries>0);
        weight = lumi*(samp->xsecv[ifile])/treeEntries;                           // (assumes you've merged filesets)
	if(isemb)  weight=1.0;
      }
      samp->weightv.push_back(weight);

      // counters
      Double_t nsel[3], nselvar[3]; for(Int_t i=0; i<3; i++)  { nsel[i]    = nselvar[i] = 0; } // events in this file
      Double_t nlowmass=0; // low mass z events (below 50)

      // loop over events
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {

        infoBr->GetEntry(ientry);

	if(getGen)  genBr->GetEntry(ientry);

	// skip non-tau events in madgraph sample
	if(ismadz && !ismadzmm && (fabs(gen->id_1_a)<15 || fabs(gen->id_1_a)>19)) continue;

        // skip non-mumu events in madgraph sample for zmm
        if(ismadzmm && (fabs(gen->id_1_a)>14 && fabs(gen->id_1_a)<20)) continue;

        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

	// require trigger
        if((isemb || isfall11) && !(info->triggerBits[kHLT_Mu17_Ele8_CaloIdL] || info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] || info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL])) continue;
        if(hasTrigs) {
          if(info->runNum <= 170053 && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdL])) continue;
          else if(info->runNum >  170053 && info->runNum <= 173199 && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdL])) continue;
          else if(info->runNum >  173199 && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL])) continue;
        }

        // No good primary vertex? Skip to next event...
        if(!info->hasGoodPV) continue;

	pvArr->Clear();
	pvBr->GetEntry(ientry);	

        // loop through muons
        vector<const mithep::TMuon*> goodMuonsv;
        vector<const mithep::TMuon*> looseMuonsv;
        muonArr->Clear();
        muonBr->GetEntry(ientry);

        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const mithep::TMuon *muon = (mithep::TMuon*)((*muonArr)[i]);

	  looseMuonsv.push_back(muon);

	  Bool_t trigmatch = muon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdL_MuObj] || muon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdL_MuObj] || muon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj] || muon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj];
	  //if(isdata && !trigmatch)                     continue;

          if(muon->pt < kMuonPt2Min)                     continue;
	  if(fabs(muon->eta) > 2.1)                      continue;

	  if(passMuonID(muon) && passMuonIsoPU(muon))  goodMuonsv.push_back(muon);
        }
	
        // loop through electrons 
        vector<mithep::TElectron*> goodElectronsv;   
        electronArr->Clear();
        electronBr->GetEntry(ientry);
        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
	  mithep::TElectron *electron = (mithep::TElectron*)((*electronArr)[i]);

          if(fabs(electron->eta) > 2.5)   	            continue;

	  Bool_t trigmatch = electron->hltMatchBits[kHLT_Mu17_Ele8_CaloIdL_EGObj] || electron->hltMatchBits[kHLT_Mu8_Ele17_CaloIdL_EGObj] || electron->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj] || electron->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj];
	  //if(isdata && !trigmatch)                          continue;

          Bool_t hasMuonTrack=kFALSE;
          for(UInt_t imu=0; imu<goodMuonsv.size(); imu++) {
            if(electron->trkID == goodMuonsv[imu]->trkID) hasMuonTrack=kTRUE;
          }
          if(hasMuonTrack) continue;

	  // clean against loose muons
          Bool_t matchLooseMuon=kFALSE;
	  for(UInt_t imu=0;imu<looseMuonsv.size();imu++) {
	    const mithep::TMuon *mu = looseMuonsv[imu];
	    if(toolbox::deltaR(electron->eta,electron->phi,mu->eta,mu->phi) < 0.3) matchLooseMuon=kTRUE;
	  }
	  if(matchLooseMuon) continue;

	  Double_t mvaValue = electronIDMVANoIPInfo->MVAValue(
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


	  if(passEleMVAID(electron,mvaValue) && passEleIsoPU(electron))    goodElectronsv.push_back(electron);
        }

        TLorentzVector lep1, lep2, dilep;  // lepton 4-vectors
        Int_t finalState=-1;	           // final state type
        Double_t svfmass = -999;
	Double_t svfmassunc = -999;

	//----------------------------------------------------------------------------------------

	if(goodMuonsv.size()<1 || goodElectronsv.size()<1) continue;

	const mithep::TMuon *mu	   = goodMuonsv[0];
	const mithep::TElectron *ele = goodElectronsv[0];

	Double_t mupt, mueta, muphi, muq, elept, eleeta, elephi, eleq;

	mupt = mu->pt; mueta = mu->eta; muphi = mu->phi; muq = mu->q;
	elept = ele->pt; eleeta = ele->eta; elephi = ele->phi; eleq = ele->q;

        rawMet = info->pfMET;
        rawMetphi = info->pfMETphi;
        Double_t met=info->pfMET,metphi=info->pfMETphi;

        if(!isdata) {
          if(elescale==kNo) {
            elept = ele->pt;
          } else if(elescale==kUp) {
            elept = 1.01*ele->pt;
	    rawMet -= 0.01*ele->pt;
	    met -= 0.01*ele->pt;
          } else if(elescale==kDown) {
            elept = 0.99*ele->pt;
            rawMet += 0.01*ele->pt;
	    met += 0.01*ele->pt;
          }
        }

	if(mupt < kMuonPt2Min  || elept < kElePt2Min) continue;
	if(mupt < kMuonPt1Min  && elept < kElePt1Min) continue;

	// trigger requirements
	if(mupt  < kMuonPt1Min) {
          if(!(info->triggerBits[kHLT_Mu8_Ele17_CaloIdL] || info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL])) continue; // if failed trig1
	}
	else if(elept < kElePt1Min) {
          if(!(info->triggerBits[kHLT_Mu17_Ele8_CaloIdL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL])) continue; // if failed trig2
	}

	if(issamesign) {
	  if(muq != eleq) continue;
	}
	else {
	  if(muq == eleq) continue;
	} // skip same-sign events

	Double_t muiso  = muonIsoPU(mu);
	Double_t eleiso = eleIsoPU(ele);

	if(mupt > elept) {
	  lep1.SetPtEtaPhiM(mupt,  mueta,  muphi,  0.105658369);
	  lep2.SetPtEtaPhiM(elept, eleeta, elephi, 0.000511);
	} else {
	  lep1.SetPtEtaPhiM(elept, eleeta, elephi, 0.000511);
	  lep2.SetPtEtaPhiM(mupt,  mueta,  muphi,  0.105658369);
	}
	dilep = lep1+lep2;

	if(elept > mupt) finalState=kEleMu; 
	else             finalState=kMuEle;

        // loop through jets      
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets = 0, nbjets = 0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	btagArray.Reset();         jptArray.Reset();        jetaArray.Reset();   npt15jets=0;
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
	  mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

	  if(doJetUnc) jet->pt *= (jetunc==kDown) ? (1-jet->unc) : (1+jet->unc);

          if(toolbox::deltaR(jet->eta,jet->phi,lep1.Eta(),lep1.Phi()) < 0.3) continue;
          if(toolbox::deltaR(jet->eta,jet->phi,lep2.Eta(),lep2.Phi()) < 0.3) continue;

          if(fabs(jet->eta) > 4.5) continue;

	  Int_t btagopt = 0;
	  if(isdata||isemb) btagopt = 1;
	  if(isfall11) btagopt = 2;

	  Bool_t btagged = isbtagged(jet,btagopt,btageff,mistag);

	  // look for b-jets
	  if((jet->pt > kBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
	    assert(npt15jets<kMaxPt15Jets);
	    btagArray.AddAt(jet->csv,npt15jets); npt15jets++;
	    if(btagged) {
	      nbjets++;
            if(!bjet || jet->pt > bjet->pt)
              bjet = jet; // leading b-jet
	    }
	  }

	  // look for jets
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

	/******** Candidate ********/        

	// calculate projection variables
	TVector3 m,e,metv;
	m.SetPtEtaPhi(mupt,0,muphi);
	e.SetPtEtaPhi(elept,0,elephi);
	metv.SetPtEtaPhi(rawMet,0,rawMetphi); // uncorrected met
	TVector3 bisector(m.Unit() + e.Unit());
	bisector = bisector.Unit();
	Double_t projVis  = (m+e).Dot(bisector);
	Double_t projMet  =  metv.Dot(bisector);
	rawprojvar  = 0.85*projVis - projMet;

	// recoil corrections
	double pU1      = 0;  //--
	double pU2      = 0;  //--

        if(corrector) corrector->CorrectAll(met,metphi,gen->vpt_a,gen->vphi_a,dilep.Pt(),dilep.Phi(), pU1, pU2, 0, 0, njets);
	metv.SetPtEtaPhi(met,0,metphi); // corrected met
	projMet  =  metv.Dot(bisector);

        // SVFit
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

        // get k-factor if necessary
        Double_t kf=1;
        if(reallyDoKf) kf = kfFHPValue(gen->vpt_a, hKFactors);

	// lepton ID corrections
	Double_t idscale = 1;
	if(doIdScale) idscale = muIDscale(mupt,mueta)*eleIDscale(elept,eleeta);

	// do vertex reweighting
	npuWgt = 1;
	if(!isdata && !isemb && doNpuRwgt) {
	  assert(puWeights);
	  Int_t npuxbin = puWeights->GetXaxis()->FindFixBin(TMath::Min(double(info->nPU), 34.499));
	  npuWgt = puWeights->GetBinContent(npuxbin);
	}

	// trigger scale factor for MC
	Double_t trigscale = 1;
	if(doTrigScale) trigscale=muTrigScale(mupt,mueta)*eleTrigScale(elept,eleeta);
	//if(isemb)       trigscale=muTrigEff(mupt,mueta)*eleTrigEff(elept,eleeta);

	// embedding weight for embedded sample
	Double_t embWgt = 1;
        Double_t pt1, eta1, phi1, pt2, eta2, phi2;
	if(isemb)    {
	  if(gen->pt_1_a > gen->pt_2_a) {
	    pt1 = gen->pt_1_a;
	    eta1 = gen->eta_1_a;
	    phi1 = gen->phi_1_a;
            pt2 = gen->pt_2_a;
            eta2 = gen->eta_2_a;
	    phi2 = gen->phi_2_a;
	  } else {
            pt2 = gen->pt_1_a;
            eta2 = gen->eta_1_a;
            phi2 = gen->phi_1_a;
            pt1 = gen->pt_2_a;
            eta1 = gen->eta_2_a;
            phi1 = gen->phi_2_a;
	  }
	  embWgt=info->embWeight*embUnfoldWgt(pt1,eta1,pt2,eta2);
	}

	// events passing selection in this file
	nsel[0]    += weight*kf*npuWgt*trigscale*idscale*embWgt;
	nselvar[0] += weight*weight*kf*kf*npuWgt*npuWgt*trigscale*trigscale*idscale*idscale*embWgt*embWgt;
	if(corrector && (gen->vmass_a < 50)) nlowmass += weight*kf*npuWgt*trigscale*idscale*embWgt;

	// passing events in whole sample 
        nSelEvents[0] += weight*kf*npuWgt*trigscale*idscale*embWgt;

        data.runNum   = info->runNum;
        data.evtNum   = info->evtNum;
        data.lumiSec  = info->lumiSec;
        data.nPV      = pvArr->GetEntriesFast();
        data.njets    = njets;
        data.nbjets   = nbjets;
	data.vpt      = (corrector) ? gen->vpt_a : 0;
        data.vphi     = (corrector) ? gen->vphi_a : 0;
	data.rawmet   = rawMet;
	data.rawmetphi= rawMetphi;
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
        data.genlpt1  = (isemb)? pt1 : 0;
        data.genleta1 = (isemb)? eta1: 0;
        data.genlphi1 = (isemb)? phi1: 0;
        data.genlpt2  = (isemb)? pt2 : 0;
        data.genleta2 = (isemb)? eta2: 0;
        data.genlphi2 = (isemb)? phi2: 0;
        data.weight   = (isdata) ? 1 : weight*kf*npuWgt*trigscale*idscale*embWgt/lumi;
        data.state    = finalState;  	   

	outtree.Fill();

      }

      printf("%8.2f +/- %-8.2f\n",nsel[0],sqrt(nselvar[0]));
      if(nlowmass > 0) printf("           ---> selected events with z mass < 50:  %10.3f (out of %15.3f)\n",nlowmass,nsel[0]);

      delete infile;
      if(corrector) {cout << "recoil corrections used" << endl; delete corrector;}
      delete fitter;
      infile=0, eventTree=0;    
    }
    outfile.Write();
    outfile.Close();

    if(samp->typev.size()>0 && samp->typev[0]==eMC)
      printf("    Yields for %1.2f/fb:",lumi/1000.);
    else
      printf("    Yields for data:    ");

    printf("%10.2f\n",nSelEvents[0]);
    cout << endl;
  }

  delete info;
  delete gen;
  delete muonArr;
  delete electronArr;
  delete jetArr;
  delete pvArr;
  delete svfitArr;


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  if(!doNpuRwgt) cout << endl << endl << "Not doing npv reweight!" << endl;
  cout << endl; 
  
  gBenchmark->Show("selectEmu");
}


//=== FUNCTION DEFINITIONS ======================================================================================
//----------------------------------------------------------------------------------------
TH1D* kfFHPInit(Int_t mH)
{
  TH1D *kfhist=0; 
  char kfilename[100];
  sprintf(kfilename, "$CMSSW_BASE/src/MitHtt/Utils/HiggsKFactors/weight_ptH_%d.root", mH);
  cout << "Getting k-factors from " << kfilename << endl;
  TFile *kfile = TFile::Open(kfilename); assert(kfile->IsOpen());
  TDirectory *kfdir = (TDirectory*)kfile->FindObjectAny("powheg_weight");
  char kfhistname[100];
  sprintf(kfhistname, "weight_hqt_fehipro_fit_%d", mH);
  cout << "kfactor histogram: " << kfhistname << endl;
  kfhist = (TH1D*)(kfdir->Get(kfhistname)); assert(kfhist);
  return kfhist;
} 
//--------------------------------------------------------------------------------------------------
Double_t kfFHPValue(Double_t pt, TH1D* hKF)
{ 
  return hKF->Interpolate(pt);
}
//----------------------------------------------------------------------------------------
Bool_t isbtagged(mithep::TJet *jet, Int_t isdata, UInt_t btageff, UInt_t mistag)
{

    // new scale factors
    // TCHEM    btag eff: 0.96 \pm 0.04         mistag rate: 0.0286 \pm 0.0003          mistag scale factor: 1.20 \pm 0.14
    // CSVM     btag eff: 0.97 \pm 0.04         mistag rate: 0.0152 \pm 0.0002          mistag scale factor: 1.10 \pm 0.11

  Bool_t btagged;
  Double_t demoteProb=0; // ~probability to demote from tagged 
  if(btageff==kNo)        demoteProb = 1-0.97; //1-0.93;  // SF = 0.97 => 0.03 = (prob to demote from tagged status)
  else if(btageff==kDown) demoteProb = 1-0.97+0.04; //1-0.93+0.07;
  else if(btageff==kUp)   demoteProb = 1-0.97-0.04; //1-0.93-0.07;
  Double_t promoteProb=0; // ~probability to promote to tagged
  if(mistag==kNo)         promoteProb = (1.10-1)*0.0152/(1-0.0152); //(1.21-1)*0.0145/(1-0.0145);  // (1-SF)*mistag = (prob. to promote to tagged status)*(1-mistag)
  else if(mistag==kDown)  promoteProb = (1.10-1+0.11)*0.0152/(1-0.0152);
  else if(mistag==kUp)    promoteProb = (1.10-1-0.11)*0.0152/(1-0.0152);

  UInt_t jetflavor = 0;
  if(isdata == 1) {
    if(jet->csv>0.679) btagged = kTRUE;
    else               btagged = kFALSE;
  } else { // MC
    if(isdata == 0)jetflavor = abs(jet->mcFlavor);
    else jetflavor = abs(jet->matchedId);
    if(jetflavor==5) {
      if(jet->csv>0.679) {
      if(randm.Uniform()>demoteProb) btagged = kTRUE;  // leave it tagged
      else                           btagged = kFALSE; // demote it
      } else                           btagged = kFALSE; // leave it untagged
    } else { // not bjet
      if(jet->csv>0.679)                   btagged = kTRUE;  // leave it tagged
      else if(randm.Uniform()<promoteProb) btagged = kTRUE;  // promote to tagged
      else                                 btagged = kFALSE; // leave it untagged
    }
  }

  return btagged;
}  
//----------------------------------------------------------------------------------------
Double_t embUnfoldWgt(Double_t pt1, Double_t eta1, Double_t pt2, Double_t eta2)
{
  TFile *unfFile1   = TFile::Open("data/unfold/v8/Unfold2D_1.root"); assert(unfFile1->IsOpen());
  TH2F  *unfWeight1 = (TH2F*) unfFile1->FindObjectAny("UnfoldDen1");
  TFile *unfFile2   = TFile::Open("data/unfold/v8/Unfold2D_2.root"); assert(unfFile2->IsOpen());
  TH2F  *unfWeight2 = (TH2F*) unfFile2->FindObjectAny("UnfoldDen2");
  double weight1 = unfWeight1->GetBinContent(unfWeight1->GetXaxis()->FindBin(eta1),unfWeight1->GetYaxis()->FindBin(pt1));
  double weight2 = unfWeight2->GetBinContent(unfWeight2->GetXaxis()->FindBin(eta2),unfWeight2->GetYaxis()->FindBin(pt2));
  unfFile1->Close();
  unfFile2->Close();
  return weight1*weight2;
}
//----------------------------------------------------------------------------------------
Double_t muIDscale(Double_t mupt, Double_t mueta)
{
  if((fabs(mueta) > 2.4) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  else if(mupt > 20) {
    if(fabs(mueta) < 1.479)     return 0.9930;
    else                          return 0.9981;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 1.479)     return 0.9455;
    else                          return 0.9604;
  }
  else {
    if(fabs(mueta) < 1.479)     return 0.9226;
    else                          return 0.9856;
  }
}
//----------------------------------------------------------------------------------------    
Double_t eleIDscale(Double_t elept, Double_t eleeta)
{
  if((fabs(eleeta) > 2.5) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  else if(elept > 20) {
    if(fabs(eleeta) < 1.479) return 0.9896;
    else                       return 1.0532;
  }
  else if(elept > 15) {
    if(fabs(eleeta) < 1.479) return 0.9783;
    else                       return 1.0623;
  }
  else {
    if(fabs(eleeta) < 1.479) return 1.1134;
    else                       return 1.1946;
  }
}
//----------------------------------------------------------------------------------------
Double_t eleTrigScale(Double_t elept, Double_t eleeta)
{
  if((fabs(eleeta) > 2.5) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  else if(elept > 30) {
    if(fabs(eleeta) < 1.479) return 1.0031;
    else                       return 1.0078;
  }
  else if(elept > 20) {
    if(fabs(eleeta) < 1.479) return 1.0012;
    else                       return 1.0040;
  }
  else if(elept > 15) {
    if(fabs(eleeta) < 1.479) return 1.0026;
    else                       return 1.0504;
  }
  else {
    if(fabs(eleeta) < 1.479) return 0.9769;
    else                       return 0.9696;
  }

}
//----------------------------------------------------------------------------------------
Double_t muTrigScale(Double_t mupt, Double_t mueta)
{
  if((fabs(mueta) > 2.4) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  else if(mupt > 30) {
    if(fabs(mueta) < 1.479) return 0.9922;
    else                      return 1.0550;
  }
  else if(mupt > 20) {
    if(fabs(mueta) < 1.479) return 0.9936;
    else                      return 1.0358;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 1.479) return 0.9918;
    else                      return 1.0712;
  }
  else {
    if(fabs(mueta) < 1.479)  return 1.0052;
    else                       return 1.0277;
  }
}
//----------------------------------------------------------------------------------------
Double_t eleTrigEff(Double_t elept, Double_t eleeta)
{
  if((fabs(eleeta) > 2.5) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  else if(elept > 30) {
    if(fabs(eleeta) < 1.479) return 0.9927;
    else                       return 0.9891;
  }
  else if(elept > 20) {
    if(fabs(eleeta) < 1.479) return 0.9862;
    else                       return 0.9931;
  }
  else if(elept > 15) {
    if(fabs(eleeta) < 1.479) return 0.9694;
    else                       return 0.9932;
  }
  else {
    if(fabs(eleeta) < 1.479) return 0.9358;
    else                       return 0.9211;
  }

}
//----------------------------------------------------------------------------------------
Double_t muTrigEff(Double_t mupt, Double_t mueta)
{
  if((fabs(mueta) > 2.4) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  else if(mupt > 30) {
    if(fabs(mueta) < 1.479) return 0.9660;
    else                      return 0.9419;
  }
  else if(mupt > 20) {
    if(fabs(mueta) < 1.479) return 0.9644;
    else                      return 0.9404;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 1.479) return 0.9625;
    else                      return 0.9474;
  }
  else {
    if(fabs(mueta) < 1.479)  return 0.9710;
    else                       return 0.8989;
  }
}
//----------------------------------------------------------------------------------------
Double_t unskimmedEntries(TString skimname)
{
  Double_t entries;
  
  skimname.ReplaceAll("_emu_skim.root","_ntuple.root");
  skimname.ReplaceAll("_emunod0_skim.root","_ntuple.root");
  TFile unskimmed(skimname);
  assert(unskimmed.IsOpen());
  TTree *tree = 0;
  unskimmed.GetObject("Events",tree);
  assert(tree);
  entries = (Double_t)tree->GetEntries();
  unskimmed.Close();

  return entries;
}
