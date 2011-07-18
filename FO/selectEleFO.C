//================================================================================================
//
//  Measure fake rate using muon triggered sample
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                   // access to gROOT, entry point to ROOT system
#include <TSystem.h>                 // interface to OS
#include <TFile.h>                   // file handle class
#include <TTree.h>                   // class to access ntuples
#include <TClonesArray.h>            // ROOT array class
#include <TLorentzVector.h>           // 4-vector class
#include <TBenchmark.h>              // class to track macro running statistics
#include <iostream>                  // standard I/O
#include <iomanip>                   // functions to format standard I/O
#include <fstream>                   // functions for file I/O
#include <sstream>                   // class for parsing strings

#include "Common/MyTools.hh"         // custom helper function

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh" 
#endif


//=== FUNCTION DECLARATIONS ======================================================================================

// check if event passes Z veto
Bool_t passZVeto(TClonesArray *muonArr, const Int_t fover);

// calculate transverse mass
Double_t calcMt(const Double_t met, const Double_t metphi, const mithep::TMuon *muon);


//=== MAIN MACRO ================================================================================================= 

void selectEleFO(const TString  foname,     // settings file
                       const TString  outputDir,  // output directory
		       const Double_t jetPtMin,   // "away" jet pT threshold
		       const Int_t    fover=1     // FO definition version
) {
  gBenchmark->Start("selectEleFO");
  

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   

  vector<TString>  fnamev;
  vector<TString>  jsonv;
  
  //
  // parse .fo file
  //
  ifstream ifs;
  ifs.open(foname.Data());
  assert(ifs.is_open());
  string line;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
     
    string fname;
    string json;
    stringstream ss(line);
    ss >> fname >> json;
    fnamev.push_back(fname);
    jsonv.push_back(json);    
  }
  ifs.close();
  
  const Double_t foPtMin = 10;
  const Double_t foPtMax = 7000;
  
  const Double_t metMax  = 20;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  UInt_t nEvents=0;
  UInt_t nDenom=0;
  UInt_t nNumer=0;
  
  //
  // Set up output ntuple file for the sample
  //
  gSystem->mkdir(outputDir,kTRUE);
  TString outName = outputDir + TString("/") + TString("fo.root");
  TFile *outFile = new TFile(outName,"RECREATE");
  TTree *outTree = new TTree("FO","FO");
  Float_t fopt, foeta, fophi, fopass, fojpt;
  Float_t jpt;
  Int_t npv;
  outTree->Branch("pt",   &fopt,  "pt/F");
  outTree->Branch("eta",  &foeta, "eta/F");
  outTree->Branch("phi",  &fophi, "phi/F");
  outTree->Branch("pass", &fopass,"pass/F");
  outTree->Branch("fojpt",&fojpt, "fojpt/F");
  outTree->Branch("jpt",  &jpt,   "jpt/F");
  outTree->Branch("npv",  &npv,   "npv/I");
  
  TFile *infile=0;
  TTree *eventTree=0;
    
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");

  // loop over samples
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

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",  &info);        TBranch *infoBr = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFJet", &jetArr);      TBranch *jetBr  = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("PV",    &pvArr);       TBranch *pvBr   = eventTree->GetBranch("PV");

    // loop over events
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {     
      infoBr->GetEntry(ientry);

      // check for certified runs
      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  
      
      nEvents++;
      
      // trigger requirement               
      ULong_t trigger = kHLT_Ele8_CaloIdL_CaloIsoVL | kHLT_Ele17_CaloIdL_CaloIsoVL | kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40;
      // ULong_t trigObj = kHLT_Mu8_MuObj | kHLT_Mu15_MuObj;  
      if(!(info->triggerBits & trigger)) continue;
        

      // Good vertex requirement    
      // if(!(info->hasGoodPV)) continue;
	
      // MET cut
      if(info->pfMET > metMax) continue;
	
      // Z veto
      electronArr->Clear();
      electronBr->GetEntry(ientry);
      if(electronArr->GetEntriesFast() > 1) continue;
			    
      jetArr->Clear();
      jetBr->GetEntry(ientry);
              
      pvArr->Clear();
      pvBr->GetEntry(ientry);

      // count denominator and numerator objects
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {	
        const mithep::TElectron* electron = (mithep::TElectron*)((*electronArr)[i]);

	// Bool_t trigmatch = electron->hltMatchBits &
	//   (kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj | kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj | kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj);
	// // if(!trigmatch)                      continue;
	// cout << electron->hltMatchBits << " " <<
	//   kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj << " " <<  kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj << " " <<
	//   kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj << endl;

	// if(ientry>300) return;

	if(electron->pt	   < foPtMin)       continue;
        if(electron->pt	   > foPtMax)       continue;
        if(fabs(electron->eta) > 2.5)       continue;
        if(!isEleFO(electron))              continue;
        
        // Bool_t hasAwayJet = kFALSE;
	const mithep::TJet *jet1=0, *fojet=0;
	Double_t highPt=-1;
        for(Int_t ijet=0; ijet<jetArr->GetEntriesFast(); ijet++) {
          const mithep::TJet* jet = (mithep::TJet*)((*jetArr)[ijet]);	  
	  
	  if(toolbox::deltaR(electron->eta,electron->phi,jet->eta,jet->phi)<0.5)
	    fojet = jet;
	  
	  if(jet->pt > highPt &&
	     (toolbox::deltaR(electron->eta,electron->phi,jet->eta,jet->phi)>1.0)
	     ) {
	    highPt = jet->pt;
	    jet1   = jet;
	  }
	  
	  if(jet->pt < jetPtMin) continue;
	  // if((jetPtMin==0) || (toolbox::deltaR(electron->eta,electron->phi,jet->eta,jet->phi)>1.0))
          //   hasAwayJet = kTRUE;
        }
        // if(!hasAwayJet) continue;
	if(highPt<jetPtMin) continue;
        
        
	nDenom++;
        
        Bool_t pass = passEleID(electron);        
        if(pass)
          nNumer++;
        
        // fill output TTree	  
        fopt   = electron->pt;
        foeta  = electron->eta;
        fophi  = electron->phi;
        fopass = pass ? 1 : 0;
	fojpt  = (fojet) ? fojet->pt : 0;
	jpt    = (jet1) ? jet1->pt : 0;
        npv    = pvArr->GetEntriesFast();
        outTree->Fill();
      }
    } 		          
    delete infile;
    infile=0, eventTree=0;
  }
  outFile->Write();
  delete outTree;
  outFile->Close();
  delete outFile;

  delete info;
  delete electronArr;
  delete jetArr; 
  delete pvArr;  
  
  
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++)
    cout << " ***  " << fnamev[ifile] << endl;
  cout << " >>>    No. of events processed: " << nEvents << endl;
  cout << " >>> No. of denominator objects: " << nDenom << endl;
  cout << " >>>   No. of numerator objects: " << nNumer << endl;     
  cout << endl;

  ofstream txtfile;
  char txtfname[100];
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << endl;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++)
    txtfile << " ***  " << fnamev[ifile] << endl;
  txtfile << " >>>    No. of events processed: " << nEvents << endl;
  txtfile << " >>> No. of denominator objects: " << nDenom << endl;
  txtfile << " >>>   No. of numerator objects: " << nNumer << endl;     
  txtfile << endl;  
  txtfile.close();
      
  cout << " <> Output saved in " << outputDir << "/" << endl;	   
  cout << endl;

  gBenchmark->Show("selectEleFO"); 
}

//=== FUNCTION DECLARATIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
Bool_t passZVeto(TClonesArray *muonArr, const Int_t fover)
{
  for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
    const mithep::TMuon* mu1 = (mithep::TMuon*)((*muonArr)[i]);
    if(mu1->pt        < 20)  continue;
    if(fabs(mu1->eta) > 2.4) continue;
    if(!isMuonFO(mu1,fover)) continue;
    for(Int_t j=i+1; j<muonArr->GetEntriesFast(); j++) {
      const mithep::TMuon* mu2 = (mithep::TMuon*)((*muonArr)[j]);
      if(mu1->q == mu2->q)     continue;
      if(mu2->pt	< 20)  continue;
      if(fabs(mu2->eta) > 2.4) continue;
      if(!isMuonFO(mu2,fover)) continue;
/*  	
      const Double_t m = 0.105659369;
      TLorentzVector vMu1; vMu1.SetPtEtaPhiM(mu1->pt,mu1->eta,mu1->phi,m);
      TLorentzVector vMu2; vMu2.SetPtEtaPhiM(mu2->pt,mu2->eta,mu2->phi,m);
      TLorentzVector vDiMu = vMu1+vMu2;
      if((vDiMu.M() < 60) || (vDiMu.M() > 120)) continue;
*/  	
      return kFALSE;  // Z candidate => fail Z veto
    }
  }
  
  return kTRUE;  // No Z candidate => pass Z veto
}

//--------------------------------------------------------------------------------------------------
Double_t calcMt(const Double_t met, const Double_t metphi, const mithep::TMuon *muon)
{
  const Double_t m = 0.105659369;
  TLorentzVector vMuon; vMuon.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, m);
  TLorentzVector vMet;  vMet.SetPtEtaPhiM(met, 0, metphi, 0);
  Double_t et = (vMuon.E())*(vMuon.Pt())/(vMuon.P());
  
  return sqrt( (et+vMet.Perp())*(et+vMet.Perp()) - (vMuon.Px()+vMet.Px())*(vMuon.Px()+vMet.Px()) - (vMuon.Py()+vMet.Py())*(vMuon.Py()+vMet.Py()) );
}
