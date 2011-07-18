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
#define BYTETOBINARYPATTERN "%d%d%d%d%d%d%d%d"
#define BYTETOBINARY(byte)  \
  (byte & 0x80 ? 1 : 0),    \
  (byte & 0x40 ? 1 : 0),    \
  (byte & 0x20 ? 1 : 0),    \
  (byte & 0x10 ? 1 : 0),    \
  (byte & 0x08 ? 1 : 0),    \
  (byte & 0x04 ? 1 : 0),    \
  (byte & 0x02 ? 1 : 0),    \
  (byte & 0x01 ? 1 : 0)
//printf("Leading text "BYTETOBINARYPATTERN"\n", BYTETOBINARY(byte));

#include "MitHtt/Utils/RecoilCorrector.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh"

// define structure for output ntuple
#include "MitHtt/Emu/Selection/EmuData.hh"
#endif

//=== FUNCTION DECLARATIONS ======================================================================================

// Initialize k-factors
TH1F* kfInit(const TString kfdata);

// Get k-factor
Double_t kfValue(const Double_t pt, const TH1F* hKF);

// Initialize vertex weights (two different ways)
TH1F* npvInit(TString fname);
vector<double> generate_flat10_weights(TString fname);

// Get weight for N vertices
Double_t npvWgtValue(Int_t npv, TH1F *hNpvWgts);

// print UInt_t in base-2 
void printtrig(UInt_t ktrig);

//=== MAIN MACRO =================================================================================================

void selectee(const TString conf,         // input file
               const TString outputDir,    // output directory
	       const Double_t lumi         // luminosity pb^-1
) {

  gBenchmark->Start("selectee");
  
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
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(0);
      samplev.back()->xsecv.push_back(xsec);
    }
  }
  ifs.close();

  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);
  const TString jsonDir = outputDir + TString("/json");
  gSystem->mkdir(jsonDir,kTRUE);


  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };  // dataset type  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };      // final state type
  
  const Double_t kMuonPt1Min = 20;
  const Double_t kMuonPt2Min = 10;
  
  const Double_t kElePt1Min  = 20;
  const Double_t kElePt2Min  = 10;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 20;
  
  Bool_t doKFactors = kFALSE;
  TString kfdata("/home/ksung/releases/CMSSW_4_1_3/src/MitPhysics/data/HWW_KFactors_PowhegToNNLL_160_7TeV.dat");

  enum { kFalse, kProper, kHack };
  UInt_t doNpvRwgt = kFalse;
  TString npvhackfname("../Selection/data/nvtxhists.root");
  TString npvproperfname("../Selection/data/estpileup-good-prv4-m10v1.root");
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
   
  Bool_t hasData = (samplev[0]->fnamev.size()>0);

  //
  // Set up NNLO-NNLL k-factor reweighting (if necessary)
  // 
  TH1F *hKFactors = (doKFactors) ? kfInit(kfdata) : 0;

  // get hist of N vtx weights
  TH1F *hNpvWgts = (doNpvRwgt==kHack)  ? npvInit(npvhackfname) : 0;
  vector<double> puwgtv = generate_flat10_weights(npvproperfname);

  TriggerEfficiency TEff;
  Double_t trigeff = 1;

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

  //
  // loop over samples
  //
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if(isam==0 && !hasData) continue;

  
    CSample* samp = samplev[isam];

    Double_t nSelEvents[3];
    for(Int_t i=0; i<3; i++)
      nSelEvents[i]=0;	
	  
    //
    // Set up output ntuple file for the sample
    //
    TString outfname = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    TFile outfile(outfname,"RECREATE");
    TTree outtree("Events","Events");

    EmuData data;
    Double_t rawMet,rawprojvar;
    UInt_t npt20jets;
    const UInt_t kMaxPt20Jets=35;
    TArrayF btagArray; btagArray.Set(kMaxPt20Jets); // array to hold b-tag values for pt-20 jets
    outtree.Branch("Events",&data.runNum,
"runNum/i:evtNum:lumiSec:nPV:njets:nbjets:met/F:metphi:mass:dphi:mt:pt:phi:pmet:pvis:lpt1:leta1:lphi1:lpt2:leta2:lphi2:jpt1:jeta1:jphi1:jpt2:jeta2:jphi2:bjpt:bjeta:bjphi:mjj:weight:state/I");

    // extra branches
    outtree.Branch("npt20jets",&npt20jets);
    outtree.Branch("btagArray",&btagArray);
    outtree.Branch("trigeff",&trigeff);
    outtree.Branch("rawMet",&rawMet);
    outtree.Branch("rawprojvar",&rawprojvar);

    //
    // loop through files
    //

    cout <<  "processing " << snamev[isam] << ":" << endl;
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      printf("        %-55s",(samp->fnamev[ifile]+"...").Data()); fflush(stdout);
      infile = new TFile(samp->fnamev[ifile]); 
      assert(infile);
      Bool_t isdata = !(samp->typev[ifile]==eMC);
      
      // setup selecting with JSON file, if necessary
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if((samp->jsonv.size()>0) && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	ifstream jsonchk; jsonchk.open(samp->jsonv[ifile].Data()); assert(jsonchk.is_open()); jsonchk.close();
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }

      RecoilCorrector *corrector=0;
      if( (samp->fnamev[ifile].Contains("ztt")) || (samp->fnamev[ifile].Contains("zll")))   corrector = new RecoilCorrector;

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");      
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
      Bool_t getGen =  corrector                               ||
	(doKFactors && samp->fnamev[ifile].Contains("-gf-"))   ||
	samp->fnamev[ifile].Contains("-vvj-")                  ||
	samp->fnamev[ifile].Contains("-zll50-");
      TBranch *genBr=0;
      if(getGen) {
        eventTree->SetBranchAddress("Gen", &gen);
        genBr = eventTree->GetBranch("Gen");
      }
      
      Double_t weight = 1; // (only initialized for each *file*)
      if(!isdata) {
        weight = lumi*(samp->xsecv[ifile])/(Double_t)eventTree->GetEntries(); // (assumes you've merged filesets)
      }
      samp->weightv.push_back(weight);
                  
      Double_t nsel[3], nselvar[3]; // (weighted) predicted number of events selected
      for(Int_t i=0; i<3; i++) { nsel[i] = nselvar[i] = 0; }

      Double_t counter[30];
      for(Int_t i=0; i<30; i++) { counter[i] = 0; }

      // loop over events
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	if(ientry>1000) break;
        infoBr->GetEntry(ientry);

        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

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
	  
	  // if(isdata  && !(muon->hltMatchBits & trigger)) continue;
          if(muon->pt < kMuonPt2Min)                     continue;
	  if(fabs(muon->eta) > 2.1)                      continue;

	  if(passMuonID(muon))  goodMuonsv.push_back(muon);
        }
	
        // loop through electrons 
        vector<const mithep::TElectron*> goodElectronsv;   
        electronArr->Clear();
        electronBr->GetEntry(ientry);
        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
          const mithep::TElectron *electron = (mithep::TElectron*)((*electronArr)[i]);

	  // if(isdata && !(electron->hltMatchBits & trigger)) continue;
          if(electron->pt < kElePt2Min)                     continue;
          if(fabs(electron->eta) > 2.5)   	            continue;
      
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
	if(goodElectronsv.size()<2) continue;

	const mithep::TElectron *ele1 = goodElectronsv[0];
	const mithep::TElectron *ele2 = goodElectronsv[1];

	if(ele1->pt < kElePt1Min  || ele2->pt < kElePt2Min) continue;
	
	if(ele1->q == ele2->q) continue; // skip same-sign events

	lep1.SetPtEtaPhiM(ele1->pt, ele1->eta,  ele1->phi, 0.000511);
	lep2.SetPtEtaPhiM(ele2->pt, ele2->eta,  ele2->phi, 0.000511);
	dilep = lep1+lep2;

	if(dilep.M() < 60 || dilep.M() > 120) continue;

	finalState=kEleEle;

        // loop through jets      
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets = 0, nbjets = 0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	btagArray.Reset(); npt20jets=0;
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

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
	
	/******** We have a candidate! Hide the fruit! ********/        

	if(getGen)  genBr->GetEntry(ientry);
   
	// calculate projection variables
	TVector3 e1,e2,metv;
	e1.SetPtEtaPhi(ele1->pt,0,ele1->phi);
	e2.SetPtEtaPhi(ele2->pt,0,ele2->phi);
	metv.SetPtEtaPhi(info->pfMET,0,info->pfMETphi); // uncorrected met
	TVector3 bisector(e1.Unit() + e2.Unit());
	bisector = bisector.Unit();
	Double_t projVis  = (e1+e2).Dot(bisector);
	Double_t projMet  =   metv.Dot(bisector);
	rawprojvar  = 0.85*projVis - projMet;

	// recoil corrections
	rawMet = info->pfMET;
	Double_t met=info->pfMET,metphi=info->pfMETphi;
	if(corrector) corrector->Correct(met,metphi,gen->vpt,gen->vphi,dilep.Pt(),dilep.Phi());
	metv.SetPtEtaPhi(met,0,metphi); // corrected met
	projMet  =   metv.Dot(bisector);

	// skip non-ww events in vvj sample
	if( (samp->fnamev[ifile].Contains("-vvj-")) && (gen->id != EGenType::kWW) )    continue;

	// skip non-tau events in madgraph sample
	if( (samp->fnamev[ifile].Contains("-zll50-"))  && (fabs(gen->id_1)<3 || fabs(gen->id_1)>6) )    continue;

        // get k-factor if necessary
        Double_t kf=1;
        if(doKFactors && samp->fnamev[ifile].Contains("-gf-"))    kf = kfValue(gen->vpt, hKFactors);      

	// do vertex reweighting
	Double_t npvWgt = 1;
	if(doNpvRwgt==kProper    && !isdata)    npvWgt = puwgtv[info->nPU];
	else if(doNpvRwgt==kHack && !isdata)    npvWgt = npvWgtValue(pvArr->GetEntriesFast(), hNpvWgts);

	nsel[0]    += weight*kf*npvWgt*trigeff; // events passing selection in this file
        nselvar[0] += weight*weight*kf*kf*npvWgt*npvWgt*trigeff*trigeff;

	// passing events in whole sample 
        nSelEvents[0] += weight*kf*npvWgt*trigeff;

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
        data.weight  = (isam==0) ? 1 : weight*kf*npvWgt*trigeff/lumi;
        data.state   = finalState;  	   

	outtree.Fill();

      }
      printf("%8.2f +/- %-8.2f\n",nsel[0],sqrt(nselvar[0]));

      delete infile;
      if(corrector) delete corrector;
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

  if(doNpvRwgt==kHack) delete hNpvWgts;
     
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;    
  cout << endl;

  cout << "NOTE: npu reweighting is not applied here (needs to be updated)." << endl << endl;
  
  gBenchmark->Show("selectee");
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
TH1F* npvInit(TString fname)
{
  TFile npvfile(fname);
  TH1F* hdata=0;
  TH1F* hmc=0;
  npvfile.GetObject("npv_data",hdata); assert(hdata); hdata->SetDirectory(0);
  hdata->Scale(1./hdata->Integral());
  npvfile.GetObject("npv_ztt",hmc); assert(hmc); hmc->SetDirectory(0);
  hmc->Scale(1./hmc->Integral());
  npvfile.Close();
  TH1F *hNpvWgts = new TH1F(*hmc);
  hNpvWgts->Reset();
  hNpvWgts->SetName("hNpvWgts");
  hNpvWgts->SetTitle("weights for mc");

  for(Int_t i=0;i<hmc->GetNbinsX()+2;i++) {
    Double_t wgt = hmc->GetBinContent(i)==0 ? 0 : hdata->GetBinContent(i)/hmc->GetBinContent(i);
    hNpvWgts->SetBinContent(i,wgt);
  }

  return hNpvWgts;
}
//----------------------------------------------------------------------------------------
Double_t npvWgtValue(Int_t npv, TH1F *hNpvWgts)
{
  Double_t wgt=1;
  Int_t bin = hNpvWgts->FindBin(npv);
  wgt = hNpvWgts->GetBinContent(bin);
  return wgt;
}  
//----------------------------------------------------------------------------------------    
vector<double> generate_flat10_weights(TString fname){

  TFile infile(fname);
  TH1D* data_npu_estimated = 0;
  infile.GetObject("pileup",data_npu_estimated);
  assert(data_npu_estimated);
  data_npu_estimated->SetDirectory(0);
  infile.Close();
  
  // see SimGeneral/MixingModule/python/mix_E7TeV_FlatDist10_2011EarlyData_inTimeOnly_cfi.py; copy and paste from there:
  const double npu_probs[25] = {0.0698146584, 0.0698146584, 0.0698146584,0.0698146584,0.0698146584,0.0698146584,
				0.0698146584,0.0698146584,0.0698146584,0.0698146584,0.0698146584 /* <-- 10*/,
				0.0630151648,0.0526654164,0.0402754482,0.0292988928,0.0194384503,
				0.0122016783,0.007207042,0.004003637,0.0020278322,
				0.0010739954,0.0004595759,0.0002229748,0.0001028162,4.58337152809607E-05 /* <-- 24 */};
  vector<double> result(25);
  double s = 0.0;
  for(int npu=0; npu<25; ++npu){
    double npu_estimated = data_npu_estimated->GetBinContent(data_npu_estimated->GetXaxis()->FindBin(npu));                              
    result[npu] = npu_estimated / npu_probs[npu];
    s += npu_estimated;
  }
  // normalize weights such that the total sum of weights over thw whole sample is 1.0, i.e., sum_i  result[i] * npu_probs[i] should be 1.0 (!)
  for(int npu=0; npu<25; ++npu){
    result[npu] /= s;
  }
  return result;
}
//----------------------------------------------------------------------------------------    
void printtrig(UInt_t ktrig)
{
  printf("  "BYTETOBINARYPATTERN""BYTETOBINARYPATTERN""BYTETOBINARYPATTERN"",
	 BYTETOBINARY(ktrig>>8),BYTETOBINARY(ktrig>>16),BYTETOBINARY(ktrig>>24));
}
