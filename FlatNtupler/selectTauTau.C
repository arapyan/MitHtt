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

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh"

// scale factros
#include "MitHtt/Utils/DataMC.hh"

#include "Output.hh"

#endif

const Double_t pi = 3.14159265358979;

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

  const Double_t kTauPtMin = 30;
  const Double_t kJetPtMin   = 30.0;
  const Double_t kBJetPtMin  = 20;
  
  Bool_t doKFactors = kTRUE;
  if(is2012)
    doKFactors = kFALSE;      // not needed in Summer12
  
  Bool_t doNpuRwgt = kTRUE;

  // Access samples and fill histograms
  TFile *infile=0;
  TTree *eventTree=0;  
  TTree *lTree=0;
 
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVfit");
  TClonesArray *tauArr      = new TClonesArray("mithep::TPFTau"); 
  
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
      Bool_t isdata     = !(samp->typev[ifile]==eMC);
      Bool_t isemb      = snamev[isam].Contains("emb");
      Bool_t reallyDoKf = doKFactors && sfname.Contains("-gf-");
      Bool_t ismadz     = sfname.Contains("ztt") || sfname.Contains("-zjets"); // madgraph z samples
      Bool_t ismadzmm   = snamev[isam].Contains("zmm") && (sfname.Contains("-zll") || sfname.Contains("-zjets")); // madgraph z samples 
      Bool_t ismssm     = sfname.Contains("-ggh-") || sfname.Contains("-bbh-");
      Bool_t doIdScale  = !isdata;
      Bool_t doTrigScale= !isdata;
      Int_t  doRecoil   = (sfname.Contains("ztt") || sfname.Contains("-zll") || sfname.Contains("zjets")) && !isemb;
      if(snamev[isam].Contains("wjets") || snamev[isam].Contains("w1jets") ||  snamev[isam].Contains("w2jets") || snamev[isam].Contains("w3jets") || snamev[isam].Contains("w4jets")) doRecoil = 2;
      if(snamev[isam].Contains("_sm_") || snamev[isam].Contains("_mssm_")) doRecoil = 3;
      Bool_t getGen     = doRecoil || ismadz ||isemb || ismssm || ismadzmm;
      if(sfname.Contains("vtth")) getGen=0;
   
      out->setupRecoil(doRecoil);
     
      // PU reweighting
      TString pileupReweightFile;
      if(!is2012) {
	cout << "Fall11 sample!" << endl;
	pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_Fall11toFull2011_PixelLumi_50bins.root";
      } else pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_S1253XTo2012_12ifb.root";
      TH1F *puWeights = 0;
      TFile *pufile = new TFile(pileupReweightFile.Data());
      puWeights = (TH1F*)pufile->Get("puWeights");

      // setup selecting with JSON file, if necessary
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if(isdata && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
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

      cout << eventTree->GetEntries() << " events" << endl;
      
      // loop over events
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	if(ientry%100000 == 0) cout << "processing " << float(ientry)/float(eventTree->GetEntriesFast()) << endl;
	if(getGen)  genBr->GetEntry(ientry);
        infoBr->GetEntry(ientry);
	//cout << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;
	
       	infoBr->GetEntry(ientry);
	
	// skip non-tau events in madgraph sample
	//if(ismadz && !ismadzmm && (fabs(gen->id_1_a)<15 || fabs(gen->id_1_a)>19)) continue;
	
        // skip non-mumu events in madgraph sample for zmm
        //if(ismadzmm && (fabs(gen->id_1_a)>14 && fabs(gen->id_1_a)<20)) continue;
	// certified run selection
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;

	// trigger
 	if(!isemb && isdata && !(info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30])) continue;
// 	if(!(info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30] || info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30])) continue;
	
	//if(!(info->triggerBits[kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20] || info->triggerBits[kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20] || info->triggerBits[kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20])) continue;

	//if(!(info->triggerBits[kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL] || info->triggerBits[kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL])) continue;

	//if(!(info->triggerBits[kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT] || info->triggerBits[kHLT_Ele27_WP80])) continue;
	//if(!(info->triggerBits[kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL])) continue;

        // good primary vertex
        if(!info->hasGoodPV) continue;
	pvArr->Clear();
	pvBr->GetEntry(ientry);	


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
	    //if(!(tauIdElectronMVA(tau,tau->antiEleID))) continue;
	    // Tau Kinematics
	    if(!(tau->pt > kTauPtMin)) continue;
	    //Tau HLT
	    Bool_t trigmatch = ((info->triggerBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj]) || (info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30] && tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj]));
	    if(isdata && !trigmatch)     continue;
	 
	    // Tau Isolation
	    if(!(tau->ringIso > 0.5)) continue; //0.921 tight
           
	    //if(!(tau->hcalOverP + tau->ecalOverP > 0.2 ||
	    //	 tau->nSignalPFChargedHadrCands > 1 ||
	    //	 tau->nSignalPFGammaCands > 0)) continue;
	    
	    if(!leadTau || tau->ringIso > leadTau->ringIso)
	      {
		subTau = leadTau;
		leadTau = tau;
	      } else if(!subTau || tau->ringIso > subTau->ringIso) {
	      subTau = tau;
	    }
          }	

 	if(!(leadTau && subTau)) continue;

	if(subTau->pt < leadTau->pt)
	  {
	    out->fillTau(leadTau,1,leadTau->ringIso > 0.884);
	    out->fillTau(subTau,0,subTau->ringIso > 0.884);
	  }
	else
	  {
	    out->fillTau(leadTau,0,leadTau->ringIso > 0.884);
	    out->fillTau(subTau,1,subTau->ringIso > 0.884);
	  }
       
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
        //if(cov_00==0 && cov_01==0 && cov_10==0 && cov_11==0) continue;

        // loop through jets      
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets = 0, nbjets = 0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	out->btagArray.Reset();	out->jptArray.Reset();	out->jetaArray.Reset();	UInt_t npt20jets=0;
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
	  mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	  
          if(toolbox::deltaR(jet->eta,jet->phi,leadTau->eta,leadTau->phi) < 0.5) continue;
          if(toolbox::deltaR(jet->eta,jet->phi,subTau->eta,subTau->phi) < 0.5) continue;

          if(fabs(jet->eta) > 4.7) continue;
	  if(!jet->id) continue;

	  // look for b-jets
	  Int_t btagopt = 0;
	  if(isdata||isemb) btagopt = 1;
	  else btagopt = 2;
	  Bool_t btagged = isbtagged(jet,btagopt,0,0);
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

	Int_t nCentralJets=0;
	if(njets>1) {
          for(Int_t i=2; i<jetArr->GetEntriesFast(); i++) {
            mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	    if(toolbox::deltaR(jet->eta,jet->phi,leadTau->eta,leadTau->phi) < 0.5) continue;
	    if(toolbox::deltaR(jet->eta,jet->phi,subTau->eta,subTau->phi) < 0.5) continue;
	    if(!(jet->pt > kJetPtMin && fabs(jet->eta)<4.7 && jet->id==1)) continue;
	    if(jet1->eta > jet2->eta && jet->eta > jet2->eta && jet->eta < jet1->eta) nCentralJets++;
	    else if(jet2->eta > jet1->eta && jet->eta > jet1->eta && jet->eta < jet2->eta) nCentralJets++;
	  }
        }
	
	out->fillJets(jet1,jet2,bjet,njets,nbjets,npt20jets,nCentralJets);

        // get k-factor if necessary
        Double_t kf=1;
        if(reallyDoKf) kf = kfFHPValue(gen->vpt_a, hKFactors);

	// do vertex reweighting
	Double_t npuWgt = 1;
	if(!isdata && !isemb && doNpuRwgt) {
	  assert(puWeights);
	  Int_t npuxbin = puWeights->GetXaxis()->FindFixBin(TMath::Min(double(info->nPU), 59.499));
	  npuWgt = puWeights->GetBinContent(npuxbin);
	}

	//W+Jets
	if(doRecoil == 2 && is2012 && gen->npartons == 1) kf  *= 0.02203*treeEntries/36541049.;
        if(doRecoil == 2 && is2012 && gen->npartons == 2) kf  *= 0.05819*treeEntries/36541049.;
        if(doRecoil == 2 && is2012 && gen->npartons == 3) kf  *= 0.03891*treeEntries/36541049.;
        if(doRecoil == 2 && is2012 && gen->npartons == 4) kf  *= 0.01934*treeEntries/36541049.;

	// lepton ID corrections
	Double_t idscale = 1;
	
	// trigger scale factor for MC
	Double_t trigscale = 1;
	if(doTrigScale)
	  trigscale =  eff2012IsoTau12fb(leadTau->pt,leadTau->eta)*eff2012IsoTau12fb(subTau->pt,subTau->eta);
	if(doTrigScale && jet1) trigscale *= eff2012Jet12fb(jet1->pt,jet1->eta);

	// embedding weight for embedded sample
	Double_t embWgt = 1;
      	if(!isdata) out->fillGen(gen);
	if(isemb)  embWgt=info->embWeight;

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

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  if(!doNpuRwgt) cout << endl << endl << "Not doing npv reweight!" << endl;
  cout << endl; 
  gBenchmark->Show("selectTauTau");
}

