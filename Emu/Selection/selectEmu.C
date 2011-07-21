//          mistag                         scale factor
// TCHEM  0.0175 \pm .0003 \pm .0038      1.21 \pm .02 \pm .17
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TH1.h>                    // histogram base class
#include <TNtuple.h>                  // class to access ntuples
#include <TTree.h>                  // class to access ntuples
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

#include "Common/MitStyleRemix.hh"  // style settings for drawing
#include "Common/CSample.hh"        // helper class for organizing input ntuple files
#include "Common/MyTools.hh"        // miscellaneous helper functions
#include "Common/CPlot.hh"          // helper class for plots

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh" 
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"   
#include "MitHtt/Ntupler/interface/TVertex.hh"   

#include "MitHtt/Utils/RecoilCorrector.hh"
#include "MitHtt/Emu/EScale/EScale.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh"

// define structure for output ntuple
#include "EmuData.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================
const Double_t pi = 3.14159265358979;
//----------------------------------------------------------------------------------------
// return the -(z-boost) of the boson, as approximated from the theta coordinates
// of the two leptons
Double_t v(Double_t t1, Double_t t2)
{
  return (-cos(t1) - cos(t2)) / (1+cos(t1+t2));
}
//----------------------------------------------------------------------------------------
// get eta starting from y
Double_t eta(Double_t pt, Double_t y, Double_t phi, Double_t m)
{
  Double_t a  = (1+exp(2*y))/(exp(2*y)-1); // intermediate term
  if(a*a<1) { cout << "a too small" << endl; assert(0); }
  Double_t E  = sqrt( a*a*(pt*pt+m*m)/(a*a-1) );
  Double_t pz = E*E - pt*pt - m*m;
  if(pz<0) { cout << "imag. pz" << endl; assert(0); }
  pz = sqrt(pz);
  if(y<0) pz *= -1;
  TLorentzVector v;
  v.SetPxPyPzE(pt*cos(phi),pt*sin(phi),pz,E);
  Double_t th = v.Theta();
  return -log(tan(th/2));
}
 
// Initialize k-factors (not implemented)
TH1F* kfInit(const TString kfdata);
// Get k-factor
Double_t kfValue(const Double_t pt, const TH1F* hKF);

// Initialize npu weights
vector<Double_t> generate_flat10_weights(TString datafname, TString mcfname);

// lepton id eff.
Double_t eleIDscale(const mithep::TElectron *ele);
Double_t muIDscale(const mithep::TMuon *mu);
Double_t eleIDscale42x(const mithep::TElectron *ele);
Double_t muIDscale42x(const mithep::TMuon *mu);

// trig. eff. numbers from kevin
Double_t eleTrigEff(const mithep::TElectron *ele);
Double_t muTrigEff(const mithep::TMuon *mu);

// get number of entries in unskimmed tree (hard-coded to look in /scratch)
Double_t unskimmedEntries(TString skimname);

//=== MAIN MACRO =================================================================================================

void selectEmu(const TString conf,         // input file
               const TString outputDir,    // output directory
	       const Double_t lumi,        // luminosity pb^-1
	       const UInt_t jetunc         // jet energy uncertainties config.
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
  enum { kMuMu, kEleEle, kEleMu, kMuEle };      // final state type
  enum { kNo, kDown, kUp };                     // jet energy uncertainties
  
  const Double_t kMuonPt1Min = 20;
  const Double_t kMuonPt2Min = 10;
  
  const Double_t kElePt1Min  = 20;
  const Double_t kElePt2Min  = 10;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 20;
  
  Bool_t doKFactors = kFALSE;
  TString kfdata("/home/ksung/releases/CMSSW_4_1_3/src/MitPhysics/data/HWW_KFactors_PowhegToNNLL_160_7TeV.dat");

  Bool_t doNpuRwgt = kTRUE;
  // write out png's to see how the npu distrib looks
  Bool_t checkNpuHists = kFALSE;
  // make hists for future reweights
  Bool_t makeNpuHists  = kFALSE;

  // Set up NNLO-NNLL k-factor reweighting (if necessary) [ not implemented ]
  TH1F *hKFactors = (doKFactors) ? kfInit(kfdata) : 0;

  // // set up trigger efficiency corrections (old method)
  // TriggerEfficiency TEff;

  // // set up energy scale/smearing
  // UInt_t escale = kCenter; // enums defined in EScale.hh
  // EScale scaler("data/data-EnergyScale.root","data/mc-EnergyScale.root");

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
    Double_t rawMet,rawprojvar,npuWgt;
    UInt_t npt20jets;
    const UInt_t kMaxPt20Jets=35;
    TArrayF btagArray; btagArray.Set(kMaxPt20Jets); // array to hold b-tag values for pt-20 jets
    outtree.Branch("Events",&data.runNum,
"runNum/i:evtNum:lumiSec:nPV:njets:nbjets:met/F:metphi:mass:dphi:mt:pt:phi:pmet:pvis:lpt1:leta1:lphi1:lpt2:leta2:lphi2:jpt1:jeta1:jphi1:jpt2:jeta2:jphi2:bjpt:bjeta:bjphi:mjj:weight:state/I");
    // extra branches
    outtree.Branch("npt20jets",&npt20jets);
    outtree.Branch("btagArray",&btagArray);
    outtree.Branch("rawMet",&rawMet);
    outtree.Branch("rawprojvar",&rawprojvar);
    outtree.Branch("npuWgt",&npuWgt);

    // Double_t counter[30]; for(Int_t i=0; i<30; i++) { counter[i] = 0; }
    
    //
    // loop through files
    //
    cout <<  "processing " << snamev[isam] << ":" << endl;
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      // TH1F hb("hb","hb",50,-1.1,1.1);
      printf("        %-55s",(samp->fnamev[ifile]+"...").Data()); fflush(stdout);
      infile = new TFile(samp->fnamev[ifile]); 
      assert(infile);

      //
      // which corrections to apply where
      //
      TString sfname    = samp->fnamev[ifile];
      Bool_t isdata     = !(samp->typev[ifile]==eMC);
      Bool_t chkboost   = kFALSE;
      Bool_t is42mc     = sfname.Contains("s11-");
      Bool_t is41mc     = sfname.Contains("p11-");
      Bool_t hasTrigs   = isdata || is42mc;
      Bool_t doRecoil   = sfname.Contains("ztt") || sfname.Contains("-zll") || sfname.Contains("zjets")
	                                         || sfname.Contains("_sm_") || sfname.Contains("_mssm_");
      Bool_t reallyDoKf = doKFactors && sfname.Contains("-gf-");
      Bool_t ismadz     = sfname.Contains("-zll") || sfname.Contains("-zjets"); // madgraph z samples
      Bool_t isvvj      = sfname.Contains("-vvj-");
      Bool_t doIdScale  = is41mc;
      Bool_t doTrigEff  = !isdata;
      Bool_t getGen     = doRecoil || reallyDoKf || isvvj || ismadz || chkboost;
      Bool_t doJetUnc   = (jetunc!=kNo) && (isdata || is42mc);

      TString basename = sfname(sfname.Last('/')+1,sfname.Last('.') - sfname.Last('/') - 1);
      TString dataNPVfname, mcNPVfname = "npu/"+basename+"-npu.root";
      if(mcNPVfname.Contains("emu_skim")) mcNPVfname.ReplaceAll("emu_skim","ntuple");

      dataNPVfname = "data/Pileup_2011_EPS_8_jul.root"; // officially produced predicted npu distribution

      UInt_t npubins = 55;
      TH1D *hpu=0, *hpuRwgt=0; // npu before/after reweighting
      if(doNpuRwgt || makeNpuHists) {
	if(makeNpuHists && sfname.Contains("skim")) { cout << "make npu distribs *before* selection" << endl; assert(0); }
	hpu     = new TH1D("hpu","hpu",npubins,-0.5,npubins-0.5); hpu->Sumw2(); // histogram of the npu distribution in MC
	hpuRwgt = new TH1D("hpuRwgt","hpuRwgt",npubins,-0.5,npubins-0.5); hpuRwgt->Sumw2();
      }
      vector<Double_t> puwgtv;
      if(!isdata && (doNpuRwgt || makeNpuHists)) puwgtv = generate_flat10_weights(dataNPVfname, mcNPVfname);

      // setup selecting with JSON file, if necessary
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if((samp->jsonv.size()>0) && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	ifstream jsonchk; jsonchk.open(samp->jsonv[ifile].Data()); assert(jsonchk.is_open()); jsonchk.close();
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }

      RecoilCorrector *corrector=0;
      if(doRecoil)
	corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfitZDat_800pb.root","$CMSSW_BASE/src/MitHtt/Utils/recoilfitZMC.root", 0xDEADBEEF);

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");      
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
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
      }
      samp->weightv.push_back(weight);

      // counters
      Double_t nsel[3], nselvar[3]; for(Int_t i=0; i<3; i++)  { nsel[i]    = nselvar[i] = 0; } // events in this file
      Double_t nlowmass=0; // low mass z events (below 50)

      // loop over events
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	// Bool_t  boocount[30]; for(Int_t i=0; i<30; i++) { boocount[i] = kFALSE; }

        infoBr->GetEntry(ientry);

	if(!isdata && (doNpuRwgt || makeNpuHists)) {
	  assert(hpu);
	  hpu->Fill(info->nPU); // do this before *any* selections
	  Double_t tmpwgt = (info->nPU >= puwgtv.size()) ? 0 : puwgtv[info->nPU];
	  assert(hpuRwgt);
	  hpuRwgt->Fill(info->nPU,tmpwgt);
	}

	if(getGen)  genBr->GetEntry(ientry);

	// skip non-ww events in vvj sample
	if(isvvj && (gen->id != EGenType::kWW)) continue;

	// skip non-tau events in madgraph sample
	if(ismadz && (fabs(gen->id_1)<3 || fabs(gen->id_1)>6)) continue;

	// boocount[0]=kTRUE;

        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

	// boocount[1]=kTRUE;

	ULong_t trigger =  kHLT_Mu8_Ele17_CaloIdL | kHLT_Mu17_Ele8_CaloIdL;
	if(hasTrigs) if(!(info->triggerBits & trigger)) continue;  // no trigger accept? Skip to next event...

	// boocount[2]=kTRUE;

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

	  ULong_t trigmatch = muon->hltMatchBits & (kHLT_Mu17_Ele8_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdL_MuObj);
	  if(info->runNum<167000) // trigger matching broken after this run
	    if(isdata && !trigmatch)                     continue;

          if(muon->pt < kMuonPt2Min)                     continue;
	  if(fabs(muon->eta) > 2.1)                      continue;

	  if(passMuonID(muon))  goodMuonsv.push_back(muon);
        }
	
        // loop through electrons 
        vector<mithep::TElectron*> goodElectronsv;   
        electronArr->Clear();
        electronBr->GetEntry(ientry);
        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
	  mithep::TElectron *electron = (mithep::TElectron*)((*electronArr)[i]);

          if(electron->pt < kElePt2Min)                     continue;
          if(fabs(electron->eta) > 2.5)   	            continue;

	  // if(!isdata) electron->pt = scaler.pt(electron->eta,electron->pt,escale); // not really proper: applies 42x corrections to 41x

	  ULong_t trigmatch = electron->hltMatchBits & (kHLT_Mu17_Ele8_CaloIdL_EGObj | kHLT_Mu8_Ele17_CaloIdL_EGObj);
	  if(info->runNum<167000)
	    if(isdata && !trigmatch)                          continue;

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

	  if(passEleID(electron))    goodElectronsv.push_back(electron);
        }

        TLorentzVector lep1, lep2, dilep;  // lepton 4-vectors
        Int_t finalState=-1;	           // final state type

	//----------------------------------------------------------------------------------------
	// if(goodMuonsv.size()>0)     boocount[3]=kTRUE;
	// if(goodElectronsv.size()>0) boocount[4]=kTRUE;

	if(goodMuonsv.size()<1 || goodElectronsv.size()<1) continue;

	const mithep::TMuon *mu	   = goodMuonsv[0];
	const mithep::TElectron *ele = goodElectronsv[0];

	if(mu->pt < kMuonPt2Min  || ele->pt < kElePt2Min) continue;
	if(mu->pt < kMuonPt1Min  && ele->pt < kElePt1Min) continue;

	// trigger requirements
	if(isdata) {
	  if(mu->pt  < kMuonPt1Min) {
	    if(!(info->triggerBits & kHLT_Mu8_Ele17_CaloIdL)) continue; // if failed trig1
	  }
	  else if(ele->pt < kElePt1Min) {
	    if(!(info->triggerBits & kHLT_Mu17_Ele8_CaloIdL)) continue; // if failed trig2
	  }
	}

	if(mu->q == ele->q) continue; // skip same-sign events

	// boocount[5]=kTRUE;;

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

        // loop through jets      
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets = 0, nbjets = 0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	btagArray.Reset(); npt20jets=0;
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
	  mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

	  if(doJetUnc) jet->pt *= (jetunc==kDown) ? (1-jet->unc) : (1+jet->unc);

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
	
	/******** Candidate ********/        

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

	// recoil corrections
	rawMet = info->pfMET;
	Double_t met=info->pfMET,metphi=info->pfMETphi;
	if(corrector) corrector->Correct(met,metphi,gen->vpt,gen->vphi,dilep.Pt(),dilep.Phi());
	metv.SetPtEtaPhi(met,0,metphi); // corrected met
	projMet  =  metv.Dot(bisector);

	// if(0.85*projVis-projMet < 25) boocount[6]=kTRUE;

	// //----------------------------------------------------------------------------------------

	// // boost into the rest frame
	// TLorentzVector nlep,plep,bos;
	// Double_t boseta = eta(gen->vpt,gen->vy,gen->vphi,gen->vmass);
	// bos.SetPtEtaPhiM(gen->vpt,boseta,gen->vphi,gen->vmass);
	// if(mu->q > 0) {
	//   plep.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
	//   nlep.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
	// } else {
	//   plep.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
	//   nlep.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
	// }
	// // plep.Boost(0,0,v(lep1.Theta(),2*pi-lep2.Theta()));
	// // nlep.Boost(0,0,v(lep1.Theta(),2*pi-lep2.Theta()));
	// plep.Boost(-bos.BoostVector());
	// nlep.Boost(-bos.BoostVector());

	// // bos.Boost(0,0,v(lep1.Theta(),2*pi-lep2.Theta()));
	// // Double_t v_me = v(lep1.Theta(),2*pi-lep2.Theta());
	// // Double_t v_re = (-bos.BoostVector()).Pz();
	// // cout << "       " << v(lep1.Theta(),2*pi-lep2.Theta()) << endl;
	// // cout << "       "; (-bos.BoostVector()).Print();
	// // cout << endl;
	// Double_t cths = bos.Vect()*plep.Vect()/(bos.P()*plep.P());
	// // Double_t e1e2 = plep.E()/nlep.E();
	// hb.Fill(cths);
	// if(ientry>200000) break;
	// //----------------------------------------------------------------------------------------

        // get k-factor if necessary
        Double_t kf=1;
        if(reallyDoKf) kf = kfValue(gen->vpt, hKFactors);

	// lepton ID corrections
	Double_t idscale = 1;
	if(doIdScale) idscale = is42mc ? muIDscale42x(mu)*eleIDscale42x(ele) : muIDscale(mu)*eleIDscale(ele);

	// do vertex reweighting
	npuWgt = 1;
	if(!isdata && doNpuRwgt) npuWgt = (info->nPU >= puwgtv.size()) ? 0 : puwgtv[info->nPU];

	// multiply by trigger effic. in MC
	Double_t trigeff = 1;
	if(doTrigEff) {
	  // old way: bad binning in efficiency graphs (eff=0 in 17<pt<14 bin)
	  // Double_t t1effold = TEff.trigEff(mu->pt,mu->eta,ele->pt,ele->eta,"Mu8","Ele17");
	  // Double_t t2effold = TEff.trigEff(mu->pt,mu->eta,ele->pt,ele->eta,"Mu15","Ele8");
	  Double_t t1eff,t2eff;
	  if(is41mc) {
	    t1eff = muTrigEff(mu)*eleTrigEff(ele);
	    t2eff = muTrigEff(mu)*eleTrigEff(ele);
	  }
	  else if(is42mc) { // this is a scale factor, not an efficiency, for 42x
	    t1eff = t2eff = 0.991*0.991;
	  }
	  else { cout << "Error: no trigger efficiency defined." << endl; assert(0); }

	  if(mu->pt < kMuonPt1Min)        trigeff = t1eff;
	  else if(ele->pt > kElePt1Min)   trigeff = t1eff + t2eff*(1-t1eff);
	  else                            trigeff = t2eff;
	}

	// events passing selection in this file
	nsel[0]    += weight*kf*npuWgt*trigeff*idscale;
	nselvar[0] += weight*weight*kf*kf*npuWgt*npuWgt*trigeff*trigeff*idscale*idscale;
	if(corrector && (gen->vmass < 50)) nlowmass += weight*kf*npuWgt*trigeff*idscale;

	// passing events in whole sample 
        nSelEvents[0] += weight*kf*npuWgt*trigeff*idscale;

	// for(Int_t i=0; i<30; i++) {
	//   if(boocount[i]) counter[i]+=weight*kf*npuWgt*trigeff*idscale;
	// }

        data.runNum  = info->runNum;
        data.evtNum  = info->evtNum;
        data.lumiSec = info->lumiSec;
        data.nPV     = pvArr->GetEntriesFast();
        data.njets   = njets;
        data.nbjets  = nbjets;
        data.met     = met;
	data.metphi  = metphi;
        data.mass    = dilep.M();
	data.dphi    = toolbox::deltaPhi(lep1.Phi(),lep2.Phi());
	data.mt      = sqrt( 2.0 * (dilep.Pt()) * met * (1.0-cos(toolbox::deltaPhi(dilep.Phi(),metphi))) );
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
        data.weight  = (isam==0) ? 1 : weight*kf*npuWgt*trigeff*idscale/lumi;
        data.state   = finalState;  	   

	outtree.Fill();

      }
      printf("%8.2f +/- %-8.2f\n",nsel[0],sqrt(nselvar[0]));
      if(nlowmass > 0) printf("           ---> selected events with z mass < 50:  %10.3f (out of %15.3f)\n",nlowmass,nsel[0]);

      if(!isdata && (doNpuRwgt || makeNpuHists)) {
	hpu    ->Scale(1./    hpu->Integral(0,    hpu->GetNbinsX()+1));
	hpuRwgt->Scale(1./hpuRwgt->Integral(0,hpuRwgt->GetNbinsX()+1));

	// write out root files of npu distributions for later use
	if(makeNpuHists) {
	  TFile puoutfile("npu/"+basename+"-npu.root","recreate");
	  hpu->Write();
	  puoutfile.Close();
	}

	// make png's of npu
	TCanvas c3("c3","c3");
	TH1D* data_npu = 0; TFile *foofile = TFile::Open(dataNPVfname); foofile->GetObject("pileup",data_npu); assert(data_npu);
	data_npu->SetDirectory(0); data_npu->Sumw2(); foofile->Close();
	data_npu->Scale(1./data_npu->Integral(0,data_npu->GetNbinsX()+1));
	data_npu->SetMarkerStyle(20);
	data_npu->SetMarkerSize(0.9);
	data_npu->Draw("EP");
	hpuRwgt->SetLineColor(kBlue);
	hpuRwgt->Draw("histsame");
	hpu->SetLineColor(kRed);
	hpu->Draw("histsame");
	if(makeNpuHists || checkNpuHists) c3.SaveAs("npu/"+basename+"-npu.png");
	
	delete hpu; delete hpuRwgt;
      }

      // TCanvas c9("c9","c9");
      // hb.Draw();
      // c9.SaveAs(basename+".png");
      
      delete infile;
      if(corrector) delete corrector;
      infile=0, eventTree=0;    
    }
    outfile.Write();
    outfile.Close();

    // // write out cutflow
    // FILE *fcut = fopen(outputDir+"/cutflow.txt","a");
    // cout << (outputDir+"/cutflow.txt").Data() << endl;
    // Double_t kssFakeWgt = 1;
    // if(snamev[isam].Contains("ss-fakes")) kssFakeWgt = 1.312;
    // fprintf(fcut,"%45s%20.2f%20.2f%20.2f%20.2f%20.2f%20.2f%20.2f\n",
    // 	    snamev[isam].Data(),kssFakeWgt*counter[0],kssFakeWgt*counter[1],
    // 	    kssFakeWgt*counter[2],kssFakeWgt*counter[3],kssFakeWgt*counter[4],kssFakeWgt*counter[5],kssFakeWgt*counter[6]);
    // fclose(fcut);

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

//--------------------------------------------------------------------------------------------------
TH1F* kfInit(const TString kfdata)
{
  cout << endl;
  cout << "Initializing k-factors from " << kfdata << "...";
  cout << endl;
  
  Int_t nbins;
  Double_t xlow, xhigh;
  string line;
  ifstream ifs;
  ifs.open(kfdata.Data());
  assert(ifs.is_open());
  
  // read in header
  getline(ifs,line); stringstream ssnbins(line); ssnbins >> nbins;
  getline(ifs,line); stringstream ssxlow(line);  ssxlow  >> xlow;
  getline(ifs,line); stringstream ssxhigh(line); ssxhigh >> xhigh;
  getline(ifs,line); 
  getline(ifs,line); 
  getline(ifs,line); 
  
  TH1F *h = new TH1F("hKFactors","",nbins,xlow,xhigh);
  while(getline(ifs,line)) {
    stringstream ss(line);
    Int_t ibin;
    Double_t scale;
    ss >> ibin >> scale;
    h->SetBinContent(ibin,scale);
  }
  ifs.close();
  
  return h;
}

//--------------------------------------------------------------------------------------------------
Double_t kfValue(const Double_t pt, const TH1F* hKF)
{
  if(pt < hKF->GetBinLowEdge(1)) {
    return hKF->GetBinContent(0);
  
  } else if(pt > hKF->GetBinLowEdge(hKF->GetNbinsX())) {
    return hKF->GetBinContent(hKF->GetNbinsX()+1);
  
  } else {
    for(Int_t ibin=1; ibin<=hKF->GetNbinsX(); ibin++) {
      if(pt >= hKF->GetBinLowEdge(ibin) && pt < hKF->GetBinLowEdge(ibin+1)) {
        return hKF->GetBinContent(ibin);
      }
    }
  }
  return 1;
}
//----------------------------------------------------------------------------------------    
vector<Double_t> generate_flat10_weights(TString datafname, TString mcfname){
  TH1D* data_npu = 0;
  TH1D* mc_npu = 0;

  TFile *infile = TFile::Open(datafname); assert(infile->IsOpen());
  infile->GetObject("pileup",data_npu);   assert(data_npu);
  data_npu->SetDirectory(0);
  infile->Close();
  data_npu->Scale(1./data_npu->Integral(0,data_npu->GetNbinsX()+1));

  infile = TFile::Open(mcfname);
  if(!infile) {
    infile = TFile::Open("npu/p11-vvj-v1g1-pu_ntuple.root-npu.root");
    cout << endl << "Warning: using NPU fro -vvj- sample. Run again to use npu from this file." << endl;
    assert(infile->IsOpen());
  }
  infile->GetObject("hpu",mc_npu); assert(mc_npu);
  mc_npu->SetDirectory(0);
  infile->Close();


  const UInt_t nbins = TMath::Min(data_npu->GetNbinsX(),mc_npu->GetNbinsX());
  
  vector<Double_t> result(nbins);
  Double_t sum = 0;
  for(UInt_t npu=0; npu<nbins; ++npu){
    Double_t data_wgt = data_npu->GetBinContent(data_npu->GetXaxis()->FindBin(npu));                              
    Double_t   mc_wgt =   mc_npu->GetBinContent(mc_npu->GetXaxis()->FindBin(npu));                              
    result[npu] = (mc_wgt==0) ? 0 : data_wgt / mc_wgt;
    sum += result[npu];
  }

  return result;
}
//----------------------------------------------------------------------------------------
Double_t muIDscale(const mithep::TMuon *mu)
{
  if((fabs(mu->eta) > 2.4) || (mu->pt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  else if(mu->pt > 20) {
    if(fabs(mu->eta) < 1.479)  return 0.9943;
    else                       return 0.9820;
  }
  else if(mu->pt > 15) {
    if(fabs(mu->eta) < 1.479)  return 0.9718;
    else                       return 0.9532;
  }
  else {
    if(fabs(mu->eta) < 1.479)  return 0.9344;
    else                       return 0.9712;
  }
}
//----------------------------------------------------------------------------------------
Double_t muIDscale42x(const mithep::TMuon *mu)
{
  if((fabs(mu->eta) > 2.4) || (mu->pt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  else if(mu->pt > 20) {
    if(fabs(mu->eta) < 1.479)  return 0.9958;
    else                       return 0.9980;
  }
  else if(mu->pt > 15) {
    if(fabs(mu->eta) < 1.479)  return 0.9705;
    else                       return 0.9894;
  }
  else {
    if(fabs(mu->eta) < 1.479)  return 0.9487;
    else                       return 1.0083;
  }
}
//----------------------------------------------------------------------------------------    
Double_t eleIDscale(const mithep::TElectron *ele)
{
  if((fabs(ele->eta) > 2.5) || (ele->pt < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  else if(ele->pt > 20) {
    if(fabs(ele->eta) < 1.479) return 0.9718;
    else                       return 0.9518;
  }
  else if(ele->pt > 15) {
    if(fabs(ele->eta) < 1.479) return 0.9434;
    else                       return 0.8843;
  }
  else {
    if(fabs(ele->eta) < 1.479) return 0.8582;
    else                       return 0.8004;
  }
}
//----------------------------------------------------------------------------------------    
Double_t eleIDscale42x(const mithep::TElectron *ele)
{
  if((fabs(ele->eta) > 2.5) || (ele->pt < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  else if(ele->pt > 20) {
    if(fabs(ele->eta) < 1.479) return 0.9865;
    else                       return 1.0107;
  }
  else if(ele->pt > 15) {
    if(fabs(ele->eta) < 1.479) return 0.9764;
    else                       return 1.0082;
  }
  else {
    if(fabs(ele->eta) < 1.479) return 0.9865;
    else                       return 1.0433;
  }
}
//----------------------------------------------------------------------------------------
// numbers from kevin
Double_t eleTrigEff(const mithep::TElectron *ele)
{
  if((fabs(ele->eta) > 2.5) || (ele->pt < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  else if(ele->pt > 20) {
    if(fabs(ele->eta) < 1.479) return 0.9970;
    else                       return 0.999;
  }
  else if(ele->pt > 15) {
    if(fabs(ele->eta) < 1.479) return 0.9947;
    else                       return 1;
  }
  else {
    if(fabs(ele->eta) < 1.479) return 0.978;
    else                       return 1;
  }
  
}
//----------------------------------------------------------------------------------------
// numbers from kevin
Double_t muTrigEff(const mithep::TMuon *mu)
{
  if((fabs(mu->eta) > 2.4) || (mu->pt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  else if(fabs(mu->eta) > 1.2) {
    if(mu->pt > 20)  return 0.9548;
    else             return 0.9478;
  }
  else if(fabs(mu->eta) > 0.8) {
    if(mu->pt > 20)  return 0.9488;
    else             return 0.9609;
  }
  else {
    if(mu->pt > 20)  return 0.9784;
    else             return 0.9674;
  }
}
//----------------------------------------------------------------------------------------
Double_t unskimmedEntries(TString skimname)
{
  Double_t entries;
  
  skimname.ReplaceAll("_emu_skim.root","_ntuple.root");
  skimname.ReplaceAll("/tmp/","/scratch/");
  TFile unskimmed(skimname);
  assert(unskimmed.IsOpen());
  TTree *tree = 0;
  unskimmed.GetObject("Events",tree);
  assert(tree);
  entries = (Double_t)tree->GetEntries();
  unskimmed.Close();

  return entries;
}
