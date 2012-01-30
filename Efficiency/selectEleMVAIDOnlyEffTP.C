//================================================================================================
//
// Select probes for electron ID efficiency with Tag&Probe method
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
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// Helper functions for Electron ID selection
#include "MitHtt/Utils/LeptonIDCuts.hh"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

// structure for output ntuple
#include "EffData.hh" 
#endif


//=== MAIN MACRO ================================================================================================= 

void selectEleMVAIDOnlyEffTP(const TString conf,              // input file
                      const TString outputDir,         // output directory
		      const Bool_t  matchGen = kFALSE  // match to generator muons
) {
  gBenchmark->Start("selectEleMVAIDOnlyEffTP");

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  //*****************************************************************************************
  //Setup MVA
  //*****************************************************************************************
  mithep::ElectronIDMVA *electronIDMVANoIPInfo = new mithep::ElectronIDMVA();
  electronIDMVANoIPInfo->Initialize("BDTG method",
                              "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml", 
                              "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml", 
                              "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml", 
                              "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml", 
                              "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml", 
                              "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml",
                              mithep::ElectronIDMVA::kNoIPInfo );

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
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  
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
    eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PV",   &pvArr);           TBranch *pvBr   = eventTree->GetBranch("PV");
    cout << "NEvents = " << eventTree->GetEntries() << endl;

    TBranch *genBr = 0;
    if(matchGen) {
      eventTree->SetBranchAddress("Gen", &gen);
      genBr = eventTree->GetBranch("Gen");
    }
    
    // Determine maximum number of events to consider
    // *** CASES ***
    // <> lumi < 0 => use all events in the sample
    // <> xsec = 0 => for data (use all events)
    const Double_t xsec = xsecv[ifile];
    Double_t weight = 1;
    if(lumi>0) { 
      if(xsec>0) { weight = lumi*xsec/(Double_t)eventTree->GetEntries(); }      
    }

    // loop over events
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      if (ientry % 1000000 == 0) cout << "Processed Event " << ientry << endl;
      infoBr->GetEntry(ientry);
    
      // check for certified runs
      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  

      //use only odd numbered events to evaluate efficiency. even numbered events were used for training
      if (typev[ifile]!=eMC && info->evtNum % 2 == 0) continue;

      // trigger requirement               
      ULong_t  trigger = 0;
      if(typev[ifile]==eDiEl) {
        trigger = kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 |
		  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17 |
	          kHLT_Ele17_CaloIdL_CaloIsoVL;
        
      } else if(typev[ifile]==eEl) {
        if(info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30)  continue;
	if(info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17)                         continue;
	if(info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)                              continue;
        
	trigger = kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT;
      }
      if(typev[ifile]!=eMC && !(info->triggerBits & trigger)) continue;     
      
      // good vertex requirement
      if(!(info->hasGoodPV)) continue;

      pvArr->Clear();
      pvBr->GetEntry(ientry);
      
      if(matchGen) genBr->GetEntry(ientry);
      
      electronArr->Clear();
      electronBr->GetEntry(ientry);
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const mithep::TElectron *tag = (mithep::TElectron*)((*electronArr)[i]);
	
	if(matchGen) {
	  Bool_t match1 = (toolbox::deltaR(tag->eta, tag->phi, gen->eta_1, gen->phi_1) < 0.5);
	  Bool_t match2 = (toolbox::deltaR(tag->eta, tag->phi, gen->eta_2, gen->phi_2) < 0.5);
	  if(!match1 && !match2)
	    continue;
	}
	
	if(tag->pt          < 20)  continue;
	if(fabs(tag->scEta) > 2.5) continue;
	if(!(passEleIsoPU(tag) && passEleMVAID(tag,
                                          electronIDMVANoIPInfo->MVAValue(
                                          tag->pt,tag->scEta,
                                          tag->sigiEtaiEta,
                                          tag->deltaEtaIn,
                                          tag->deltaPhiIn,
                                          tag->HoverE,
                                          tag->d0,
                                          tag->dz,
                                          tag->fBrem,
                                          tag->EoverP,
                                          tag->ESeedClusterOverPOut,
                                          TMath::Sqrt(tag->sigiPhiiPhi),
                                          tag->nBrem,
                                          (1.0/(tag->scEt * TMath::CosH(tag->scEta)) - 1/tag->p),
                                          tag->ESeedClusterOverPIn,
                                          tag->ip3d,
                                          tag->ip3dSig )))) continue;


        if(typev[ifile]==eDiEl &&
           !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (tag->hltMatchBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj)) &&
           !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (tag->hltMatchBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj)) &&
           !((info->triggerBits & kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17) && (tag->hltMatchBits & kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj)) &&
           !((info->triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL) && (tag->hltMatchBits & kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj)) )
          continue;

        if(typev[ifile]==eEl &&
           !((info->triggerBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj)) &&
           !((info->triggerBits & kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT) && (tag->hltMatchBits & kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj)) &&
           !((info->triggerBits & kHLT_Ele52_CaloIdVT_TrkIdT) && (tag->hltMatchBits & kHLT_Ele52_CaloIdVT_TrkIdT_EleObj)) )
          continue;
        
	const Double_t m = 0.000511;
	TLorentzVector vtag;
	vtag.SetPtEtaPhiM(tag->pt, tag->eta, tag->phi, m);
        
	for(Int_t j=0; j<electronArr->GetEntriesFast(); j++) {
	  if(i==j) continue;
	  
	  const mithep::TElectron *probe = (mithep::TElectron*)((*electronArr)[j]);
	  if(probe->q == tag->q) continue;

          if(typev[ifile]==eDiEl &&
             !((info->triggerBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30) && (probe->hltMatchBits & kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj)) &&
             !((info->triggerBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17) && (probe->hltMatchBits & kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj)) &&
             !((info->triggerBits & kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17) && (probe->hltMatchBits & kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj)))
            continue;

          if(matchGen) {
            Bool_t match1 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_1, gen->phi_1) < 0.5);
            Bool_t match2 = (toolbox::deltaR(probe->eta, probe->phi, gen->eta_2, gen->phi_2) < 0.5);
            if(!match1 && !match2)
              continue;
          }

	  if(fabs(probe->scEta) > 2.5) continue;	  
	  TLorentzVector vprobe;
	  vprobe.SetPtEtaPhiM(probe->pt, probe->eta, probe->phi, m);
	  
	  TLorentzVector vdielectron = vtag + vprobe;
	  if((vdielectron.M()<massLo) || (vdielectron.M()>massHi)) continue;	  	  

	  if(!passEleIsoPU(probe)) continue;

	  nProbes++;
	  
 	  Bool_t pass = passEleMVAID(probe, electronIDMVANoIPInfo->MVAValue(
                                          probe->pt,probe->scEta,
                                          probe->sigiEtaiEta, 
                                          probe->deltaEtaIn,
                                          probe->deltaPhiIn, 
                                          probe->HoverE,
                                          probe->d0,
                                          probe->dz, 
                                          probe->fBrem,
                                          probe->EoverP,
                                          probe->ESeedClusterOverPOut,
                                          TMath::Sqrt(probe->sigiPhiiPhi),
                                          probe->nBrem,
                                          (1.0/(probe->scEt * TMath::CosH(probe->scEta)) - 1/probe->p), 
                                          probe->ESeedClusterOverPIn,
                                          probe->ip3d,
                                          probe->ip3dSig ));


          // Fill tree
	  data.mass   = vdielectron.M();
	  data.pt     = probe->pt;
          data.eta    = probe->scEta;
          data.phi    = probe->phi;
          data.weight = weight;
	  data.q      = probe->q;
	  data.npv    = pvArr->GetEntriesFast();
	  data.npu    = info->nPU;
          data.pass   = (pass) ? 1 : 0;
          data.rho    = info->rho;
	  outTree->Fill();	  
	}
      }
    }
    delete infile;
    infile=0, eventTree=0;    
  }
  delete info;
  delete gen;
  delete electronArr;
  
     
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
      
  gBenchmark->Show("selectEleMVAIDOnlyEffTP"); 
}
