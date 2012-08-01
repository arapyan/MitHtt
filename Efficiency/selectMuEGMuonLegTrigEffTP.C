//================================================================================================
//
// Select probes for single muon trigger efficiency with Tag&Probe method
//
//  * outputs ROOT file with a TTree of probes
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // class for 4-vector calculations
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "MitHtt/Common/MyTools.hh"        // miscellaneous helper functions

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TPhoton.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// helper functions for lepton ID selection
#include "MitHtt/Utils/LeptonIDCuts.hh"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

// structure for output ntuple
#include "EffData.hh" 
#endif


//=== MAIN MACRO ================================================================================================= 

void selectMuEGMuonLegTrigEffTP(const TString conf,              // input file
                                    const TString outputDir,         // output directory
                                    Int_t RunRange = 0,              // Run Range
                                    const Bool_t  matchGen = kFALSE  // match to generator muons
  ) {
  gBenchmark->Start("selectSingleMuEffTP");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  // mass region
  Double_t massLo;
  Double_t massHi;

  Double_t lumi;              // luminosity (pb^-1)
  
  vector<TString>  fnamev;    // sample files 
  vector<Int_t>    typev;     // dataset type 
  vector<Double_t> xsecv;     // per file cross section
  vector<TString>  jsonv;     // per file JSON file

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
    
    if(state==0) {  // general settings
      stringstream ss1(line); ss1 >> lumi;
      getline(ifs,line);
      stringstream ss2(line); ss2 >> massLo >> massHi; 
      
    } else if(state==1) {  // define data sample
      string fname;
      Int_t type;
      Double_t xsec;
      string json;
      stringstream ss(line);
      ss >> fname >> type >> xsec >> json;
      fnamev.push_back(fname);
      typev.push_back(type);
      xsecv.push_back(xsec);
      jsonv.push_back(json);        
    }
  }
  ifs.close();
  
  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };  // dataset type


  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  Double_t nProbes = 0;
  
  //
  // Set up output ntuple
  //
  gSystem->mkdir(outputDir,kTRUE);
  TFile *outFile = new TFile(outputDir+TString("/probes.root"),"RECREATE"); 
  TTree *outTree = new TTree("Events","Events");
  EffData data;
  outTree->Branch("Events",&data.mass,"mass/F:pt:eta:phi:weight:q/I:npv/i:npu:pass:runNum:lumiSec:evtNum:rho/F");

  TFile *infile=0;
  TTree *eventTree=0;
  
  // Data structures to store info from TTrees
  mithep::TEventInfo *info = new mithep::TEventInfo();
  mithep::TGenInfo   *gen  = new mithep::TGenInfo();
  TClonesArray *muonArr    = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr= new TClonesArray("mithep::TElectron");
  TClonesArray *pvArr      = new TClonesArray("mithep::TVertex");
  TClonesArray *jetArr     = new TClonesArray("mithep::TJet");
  
  // loop over files  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  

    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);

    Bool_t hasJSON = kFALSE;
    mithep::RunLumiRangeMap rlrm;
    if(jsonv[ifile].CompareTo("NONE")!=0) { 
      hasJSON = kTRUE;
      rlrm.AddJSONFile(jsonv[ifile].Data()); 
    }
  
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);  
    eventTree->SetBranchAddress("Info", &info);            TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Muon", &muonArr);         TBranch *muonBr     = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PV",   &pvArr);           TBranch *pvBr       = eventTree->GetBranch("PV");
    eventTree->SetBranchAddress("PFJet", &jetArr);         TBranch *jetBr      = eventTree->GetBranch("PFJet");
    TBranch *genBr = 0;
    if(matchGen) {
      eventTree->SetBranchAddress("Gen", &gen);
      genBr = eventTree->GetBranch("Gen");
    }
    
    const Double_t xsec = xsecv[ifile];
    Double_t weight = 1;
    if(lumi>0) { 
      if(xsec>0) { weight = lumi*xsec/(Double_t)eventTree->GetEntries(); }      
    }

    cout << "Total Events = " << eventTree->GetEntries() << endl;

    // loop over events
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      infoBr->GetEntry(ientry);
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
      //********************************************************************************
      // check for certified runs
      //********************************************************************************
      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

      //********************************************************************************
      //RunRange Requirement
      //********************************************************************************
      Bool_t passRunRange = kTRUE;
      if (RunRange == 1) {
        if (!(info->runNum <= 170053)) passRunRange = kFALSE;
      } else if (RunRange == 2) {
        if (!(info->runNum > 173198)) passRunRange = kFALSE;
      }
      if (!passRunRange) continue;

      //********************************************************************************
      // Tag & Probe Trigger Requirement
      //********************************************************************************
      if(!(info->triggerBits[kHLT_Ele27_WP80])) continue;
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;
      
      if(matchGen) genBr->GetEntry(ientry);
      
      muonArr->Clear();
      electronArr->Clear();
      jetArr->Clear();
      pvArr->Clear();
      muonBr->GetEntry(ientry);
      electronBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);
      pvBr->GetEntry(ientry);

      Int_t NElectrons = 0;
      Int_t NMuons = 0;
      const mithep::TElectron *selectedElectron = 0;
      const mithep::TMuon *selectedMuon = 0;
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const mithep::TElectron *tag = (mithep::TElectron*)((*electronArr)[i]);
        
	if(matchGen) {
	  Bool_t match1 = (toolbox::deltaR(tag->eta, tag->phi, gen->eta_1_a, gen->phi_1_a) < 0.5);
	  Bool_t match2 = (toolbox::deltaR(tag->eta, tag->phi, gen->eta_2_a, gen->phi_2_a) < 0.5);
	  if(!match1 && !match2)
	    continue;
          if(fnamev[ifile].Contains("zjets") && (fabs(gen->id_1_a)<15 || fabs(gen->id_1_a)>19)) continue;
	}
	
        //********************************************************************************
        // Tag Requirement
        //********************************************************************************
	if(tag->pt        < 10)  continue;
	if(fabs(tag->scEta) > 2.5) continue;
	if(!(pass2012EleMVAID(tag, kMedium) && passEleIsoPU(tag))) continue;
   
        if(!((info->triggerBits[kHLT_Ele27_WP80]) && (tag->hltMatchBits[kHLT_Ele27_WP80_EleObj])))
          //!((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            //(tag->hltMatchBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj))
          //&&
          //!((info->triggerBits & kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            //(tag->hltMatchBits & kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj))
          //&&
          //!((info->triggerBits & kHLT_Ele52_CaloIdVT_TrkIdT) && 
            //(tag->hltMatchBits & kHLT_Ele52_CaloIdVT_TrkIdT_EleObj))
          //&&
          //!((info->triggerBits & kHLT_Ele65_CaloIdVT_TrkIdT) && 
            //(tag->hltMatchBits & kHLT_Ele65_CaloIdVT_TrkIdT_EleObj))
          //&&
          //!((info->triggerBits & kHLT_Ele80_CaloIdVT_TrkIdT) && 
            //(tag->hltMatchBits & kHLT_Ele80_CaloIdVT_TrkIdT_EleObj))
          //&&
          //!((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && 
            //(tag->hltMatchBits & kHLT_Ele32_WP70_EleObj))
          //)
          continue;   

        NElectrons++;
        selectedElectron = tag;
      }
	
      for(Int_t j=0; j<muonArr->GetEntriesFast(); j++) {
	  
	const mithep::TMuon *probe = (mithep::TMuon*)((*muonArr)[j]);

	if(matchGen) {
	  Bool_t match1 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_1_a, gen->phi_1_a) < 0.5);
	  Bool_t match2 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_2_a, gen->phi_2_a) < 0.5);
	  if(!match1 && !match2)
	    continue;
	}
        if(!(passTightPFMuonID(probe) && passMuonIsoPU(probe)))     continue;

        NMuons++;
        selectedMuon = probe;
      }

      if(!(NElectrons==1 && NMuons==1)) continue;
      if(selectedMuon->pt < 10.0) continue;
      if(selectedElectron->q == selectedMuon->q) continue;
      if (toolbox::deltaR(selectedElectron->eta, selectedElectron->phi, selectedMuon->eta, selectedMuon->phi) < 0.5) continue;

      //ttbar selection : 2 jets, at least 1 b-tagged
      Int_t NJets = 0;
      Int_t NBTags = 0;
      for(Int_t jetIndex=0; jetIndex<jetArr->GetEntriesFast(); ++jetIndex) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[jetIndex]);

        if (toolbox::deltaR(jet->eta, jet->phi, selectedElectron->eta, selectedElectron->phi) < 0.5) continue;
        if (toolbox::deltaR(jet->eta, jet->phi, selectedMuon->eta, selectedMuon->phi) < 0.5) continue;
        if (jet->pt > 30.0) {
          NJets++;
          if (jet->tche > 3.3) { NBTags++; }
        }
       }
      //if (NJets != 2 && NBTags >= 1) continue;
	  
      nProbes ++;
	  
      Bool_t pass = kFALSE;

      pass = (//(info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && 
		selectedMuon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj]
                ||
                //(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && 
		selectedMuon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj]);

      /*if(typev[ifile]==eMC)
        pass = (((info->triggerBits[kHLT_Mu17_Ele8_CaloIdL]) && (selectedMuon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdL_MuObj]))
                ||
                ((info->triggerBits[kHLT_Mu8_Ele17_CaloIdL]) && (selectedMuon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdL_MuObj]))
                ||
                ((info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL]) && (selectedMuon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj]))
                ||
                ((info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL]) && (selectedMuon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj]))
                ||
                (info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && selectedMuon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj])
                ||
                (info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && selectedMuon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj])) ;*/

      if(selectedMuon->pt < 10) pass=kFALSE;

      if(selectedMuon->pt < 20) {
        //if(!(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL])) pass = kFALSE;
	if(!(selectedMuon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj])) pass = kFALSE;
      }

	  
      // Fill tree
      data.mass    = 90;
      data.pt      = selectedMuon->pt;
      data.eta     = selectedMuon->eta;
      data.phi     = selectedMuon->phi;
      data.weight  = weight;
      data.q       = selectedMuon->q;
      data.npv     = pvArr->GetEntriesFast();
      data.npu     = info->nPU;
      data.pass    = (pass) ? 1 : 0;
      data.runNum  = info->runNum;
      data.lumiSec = info->lumiSec;
      data.evtNum  = info->evtNum;
      data.rho     = info->rho;
      outTree->Fill();
    }
    delete infile;
    infile=0, eventTree=0;    
  }
  delete info;
  delete gen;
  delete muonArr;
    
  //--------------------------------------------------------------------------------------------------------------
  // Output
  //==============================================================================================================
   
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  if(lumi>0) {
    cout << " L_int = " << lumi << "/pb" << endl;
    cout << endl;
  }
  cout << " Number of probes selected: " << nProbes << endl;
  
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  cout << endl;
  cout << "  <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;  
      
  gBenchmark->Show("selectSingleMuEffTP"); 
}
