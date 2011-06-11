//root -l -b -q MitHtt/Hww/FakeRate/ComputeElectronFakeRate_Data.C+\(\)
//================================================================================================
//
// HWW selection macro
//
//  * plots distributions associated with selected events
//  * prints list of selected events from data
//  * outputs ROOT files of events passing selection for each sample, 
//    which can be processed by plotSelect.C
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TBranch.h>                // class to access branches in TTree
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TGraphAsymmErrors.h>     
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include "MitPlots/Style/interface/MitStyle.h"

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/MitHttDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TPhoton.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHtt/Utils/interface/EfficiencyUtils.h"
#include "MitHtt/Utils/interface/PlotUtils.h"

#endif

using namespace std;
using namespace mithep;

//*************************************************************************************************
//Convert int to string
//*************************************************************************************************
string IntToString(int i) {
  char temp[100];
  sprintf(temp, "%d", i);
  string str = temp;
  return str;
}


//=== FUNCTION DECLARATIONS ======================================================================================



// print event dump
Bool_t passElectronNumeratorCuts(const mithep::TElectron *ele);
Bool_t passElectronDenominatorCuts(Int_t triggerBits, const mithep::TElectron *ele, Int_t DenominatorType, string SampleType);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);
void DoComputeElectronFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename, Int_t IsoChoice = 0);

//------------------------------------------------------------------------------
// getTreeFromFile
//------------------------------------------------------------------------------
TTree* getTreeFromFile(const char* infname, const char* tname)
{
  bool verbose = false;

  if (verbose) {
    cout << "--- Open file " << infname << endl;
  }
  
  TFile* inf = new TFile(infname,"read");
  assert(inf);

  TTree* t = (TTree*)inf->Get(tname);
  
  if (!t) {
    TDirectory *dir = (TDirectory*)inf->FindObjectAny("HttNtuplerMod");
    if (!dir) {
      cout << "Cannot get Directory HttNtuplerMod from file " << infname << endl;
      assert(dir);
    }
    t = (TTree*)dir->Get(tname);
  }

  if (!t) {
    cout << "Cannot get Tree with name " << tname << " from file " << infname << endl;
  }
  assert(t);


  if (verbose) {
    cout << "---\tRecovered tree " << t->GetName()
	 << " with "<< t->GetEntries() << " entries" << endl;
  }
  
  return t;
}


//=== MAIN MACRO =================================================================================================
void ComputeElectronFakeRate_Data() {


    DoComputeElectronFakeRate("/data/blue/vdutta/hist/FR_0525/r11a-del-pr-v_ntuple.root","ElectronFakeRate","ElectronFakeRate.root");

}



void DoComputeElectronFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename, Int_t IsoChoice)
{  
  gBenchmark->Start("WWTemplate");

  Double_t LUMI = 26.5;

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  //*****************************************************************************************
  //Define Pt bins
  //*****************************************************************************************
  vector<double> ptbins;
//   ptbins.push_back(10);  
//   ptbins.push_back(80);  


  ptbins.push_back(10);  
//   ptbins.push_back(11);  
//   ptbins.push_back(12);  
//   ptbins.push_back(13);  
//   ptbins.push_back(14);  
  ptbins.push_back(15);  
//   ptbins.push_back(16);  
//   ptbins.push_back(17);  
//   ptbins.push_back(18);  
//   ptbins.push_back(19);  
  ptbins.push_back(20);  
//   ptbins.push_back(22.5);  
  ptbins.push_back(25);  
//   ptbins.push_back(27.5);  
  ptbins.push_back(30);  
//   ptbins.push_back(32.5);  
  ptbins.push_back(35);  
//   ptbins.push_back(37.5);  
  ptbins.push_back(40);  
  ptbins.push_back(50);  
//   ptbins.push_back(60);  
//   ptbins.push_back(70);  
//   ptbins.push_back(80);  
  ptbins.push_back(100);  


  vector<double> etabins;
  etabins.push_back(-2.50);
  etabins.push_back(-2.00);
  etabins.push_back(-1.50);
  etabins.push_back(-1.00);
  etabins.push_back(-0.50);
  etabins.push_back(0.00);
  etabins.push_back(0.50);
  etabins.push_back(1.00);
  etabins.push_back(1.50);
  etabins.push_back(2.00);
  etabins.push_back(2.50);

  vector<double> phibins;
  phibins.push_back(-3.25);
  phibins.push_back(-2.75);
  phibins.push_back(-2.25);
  phibins.push_back(-1.75);
  phibins.push_back(-1.25);
  phibins.push_back(-0.75);
  phibins.push_back(-0.25);
  phibins.push_back(0.25);
  phibins.push_back(0.75);
  phibins.push_back(1.25);
  phibins.push_back(1.75);
  phibins.push_back(2.25);
  phibins.push_back(2.75);
  phibins.push_back(3.25);

  vector<double> ptbins2D;
  ptbins2D.push_back(10);  
  ptbins2D.push_back(12.5);  
  ptbins2D.push_back(15);  
  ptbins2D.push_back(17.5);  
  ptbins2D.push_back(20);  
  ptbins2D.push_back(25);  
  ptbins2D.push_back(30);  
  ptbins2D.push_back(40);  
  ptbins2D.push_back(80);  

  vector<double> etabins2D;
  etabins2D.push_back(-2.50);
  etabins2D.push_back(-2.00);
  etabins2D.push_back(-1.50);
  etabins2D.push_back(-1.00);
  etabins2D.push_back(-0.50);
  etabins2D.push_back(0.00);
  etabins2D.push_back(0.50);
  etabins2D.push_back(1.00);
  etabins2D.push_back(1.50);
  etabins2D.push_back(2.00);
  etabins2D.push_back(2.50);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *leadingJetPt = new TH1F("leadingJetPt" , "; p_{T} [GeV/c] ; Number of Events ",  200, 0 , 200);

  //3D array, indices give: [denominatorType][SampleType][ptThreshold]

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt;
  vector<vector<vector<TH1F*> > > DenominatorVector_Eta;
  vector<vector<vector<TH1F*> > > DenominatorVector_Phi;
  vector<vector<vector<TH2F*> > > DenominatorVector_PtEta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt;
  vector<vector<vector<TH1F*> > > NumeratorVector_Eta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Phi;
  vector<vector<vector<TH2F*> > > NumeratorVector_PtEta;
  vector<vector<vector<TH1F*> > > LeptonJetPt;
  vector<vector<vector<TH1F*> > > DenominatorIsolation;

  vector<Double_t> denominatorType;
  denominatorType.push_back(1);
  denominatorType.push_back(2);
  denominatorType.push_back(3);
  denominatorType.push_back(4);
//   denominatorType.push_back(5);  
//   denominatorType.push_back(6);
//    denominatorType.push_back(100);
  vector<string> sampleLabel;
  sampleLabel.push_back("Ele8Sample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLSample");
  sampleLabel.push_back("Ele17CaloIdLCaloIsoVLSample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLJet40Sample");
  sampleLabel.push_back("PhotonJetsSample");
  vector<Double_t> ptThreshold;
  ptThreshold.push_back(0);
  ptThreshold.push_back(15);
  ptThreshold.push_back(20);
  ptThreshold.push_back(25);
  ptThreshold.push_back(30);
  ptThreshold.push_back(35);
  ptThreshold.push_back(40);
  ptThreshold.push_back(45);
  ptThreshold.push_back(50);
  ptThreshold.push_back(70);
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    vector<vector<TH1F*> > tmpDenominatorVector_Pt;
    vector<vector<TH1F*> > tmpDenominatorVector_Eta;
    vector<vector<TH1F*> > tmpDenominatorVector_Phi;
    vector<vector<TH2F*> > tmpDenominatorVector_PtEta;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt;
    vector<vector<TH1F*> > tmpNumeratorVector_Eta;
    vector<vector<TH1F*> > tmpNumeratorVector_Phi;
    vector<vector<TH2F*> > tmpNumeratorVector_PtEta;
    vector<vector<TH1F*> > tmpLeptonJetPt;
    vector<vector<TH1F*> > tmpDenominatorIsolation;
    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      vector<TH1F*>  tmptmpDenominatorVector_Pt;
      vector<TH1F*>  tmptmpDenominatorVector_Eta;
      vector<TH1F*>  tmptmpDenominatorVector_Phi;
      vector<TH2F*>  tmptmpDenominatorVector_PtEta;
      vector<TH1F*>  tmptmpNumeratorVector_Pt;
      vector<TH1F*>  tmptmpNumeratorVector_Eta;
      vector<TH1F*>  tmptmpNumeratorVector_Phi;
      vector<TH2F*>  tmptmpNumeratorVector_PtEta;
      vector<TH1F*>  tmptmpLeptonJetPt;
      vector<TH1F*>  tmptmpDenominatorIsolation;
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
        TH1F *histDenominator_Pt = new TH1F(("histDenominator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str(), "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominator_Eta = new TH1F(("histDenominator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, -3.0 , 3.0);
        TH1F *histDenominator_Phi = new TH1F(("histDenominator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH2F *histDenominator_PtEta = new TH2F(("histDenominator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, -3.5, 3.5);
        TH1F *histNumerator_Pt = new TH1F(("histNumerator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histNumerator_Eta = new TH1F(("histNumerator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, -3.0 , 3.0);
        TH1F *histNumerator_Phi = new TH1F(("histNumerator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH2F *histNumerator_PtEta = new TH2F(("histNumerator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, -3.5, 3.5);        
        TH1F *histLeptonJetPt = new TH1F(("histLeptonJetPt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominatorIsolation = new TH1F(("histDenominatorIsolation_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; PF RelIso ; Number of Events ",  100, 0 , 1.0);

        tmptmpDenominatorVector_Pt.push_back(histDenominator_Pt);
        tmptmpDenominatorVector_Eta.push_back(histDenominator_Eta);
        tmptmpDenominatorVector_Phi.push_back(histDenominator_Phi);
        tmptmpDenominatorVector_PtEta.push_back(histDenominator_PtEta);
        tmptmpNumeratorVector_Pt.push_back(histNumerator_Pt);
        tmptmpNumeratorVector_Eta.push_back(histNumerator_Eta);
        tmptmpNumeratorVector_Phi.push_back(histNumerator_Phi);
        tmptmpNumeratorVector_PtEta.push_back(histNumerator_PtEta);
        tmptmpLeptonJetPt.push_back(histLeptonJetPt);
        tmptmpDenominatorIsolation.push_back(histDenominatorIsolation);

      }
        tmpDenominatorVector_Pt.push_back(tmptmpDenominatorVector_Pt);
        tmpDenominatorVector_Eta.push_back(tmptmpDenominatorVector_Eta);
        tmpDenominatorVector_Phi.push_back(tmptmpDenominatorVector_Phi);
        tmpDenominatorVector_PtEta.push_back(tmptmpDenominatorVector_PtEta);
        tmpNumeratorVector_Pt.push_back(tmptmpNumeratorVector_Pt);
        tmpNumeratorVector_Eta.push_back(tmptmpNumeratorVector_Eta);
        tmpNumeratorVector_Phi.push_back(tmptmpNumeratorVector_Phi);
        tmpNumeratorVector_PtEta.push_back(tmptmpNumeratorVector_PtEta);
        tmpLeptonJetPt.push_back(tmptmpLeptonJetPt);
        tmpDenominatorIsolation.push_back(tmptmpDenominatorIsolation);
    }
    DenominatorVector_Pt.push_back(tmpDenominatorVector_Pt);
    DenominatorVector_Eta.push_back(tmpDenominatorVector_Eta);
    DenominatorVector_Phi.push_back(tmpDenominatorVector_Phi);
    DenominatorVector_PtEta.push_back(tmpDenominatorVector_PtEta);
    NumeratorVector_Pt.push_back(tmpNumeratorVector_Pt);
    NumeratorVector_Eta.push_back(tmpNumeratorVector_Eta);
    NumeratorVector_Phi.push_back(tmpNumeratorVector_Phi);
    NumeratorVector_PtEta.push_back(tmpNumeratorVector_PtEta);
    LeptonJetPt.push_back(tmpLeptonJetPt);
    DenominatorIsolation.push_back(tmpDenominatorIsolation);
  }

  ofstream eventListFile("eventList.txt");

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================  
  
  //
  // Access samples and fill histograms
  TTree *eventTree=0;  
   
  // Data structures to store info from TTrees
  mithep::TEventInfo *info    = new mithep::TEventInfo();
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *muonArr = new TClonesArray("mithep::TMuon");
  TClonesArray *jetArr = new TClonesArray("mithep::TJet");
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
//   rlrm.AddJSONFile("Cert_TopOct22_Merged_135821-148058_allPVT.txt"); 
  rlrm.AddJSONFile("/home/vdutta/cms/root/json/Cert_160404-163869_7TeV_PromptReco_Collisions11_JSON.txt"); 
//   hasJSON = kFALSE;

  //********************************************************
  // Get Tree
  //********************************************************
  eventTree = getTreeFromFile(inputFilename.c_str(),"Events"); 
  TBranch *infoBr;
  TBranch *electronBr;
  TBranch *muonBr;
  TBranch *jetBr;
  TBranch *photonBr;


  //*****************************************************************************************
  //Loop over Data Tree
  //*****************************************************************************************
  // Set branch address to structures that will store the info  
  eventTree->SetBranchAddress("Info",       &info);      infoBr       = eventTree->GetBranch("Info");
  eventTree->SetBranchAddress("Electron", &electronArr); electronBr = eventTree->GetBranch("Electron");
  eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
  eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
  eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

  cout << "Total Events: " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
      //  for(UInt_t ientry=0; ientry<1000000; ientry++) {
    infoBr->GetEntry(ientry);
		
    if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

    //Double_t eventweight = info->eventweight;

    mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
    if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     
//     if (!(info->nPV0 <= 4)) continue;

    //********************************************************
    // Load the branches
    //********************************************************
    electronArr->Clear(); 
    muonArr->Clear(); 
    photonArr->Clear(); 
    jetArr->Clear(); 
    electronBr->GetEntry(ientry);
    muonBr->GetEntry(ientry);
    photonBr->GetEntry(ientry);
    jetBr->GetEntry(ientry);

    Double_t met = info->pfMET;

    Int_t NElectrons = electronArr->GetEntries();
      
    for(Int_t i=0; i<electronArr->GetEntries(); i++) {
      const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
      
      Double_t leadingJetPt = -1; // pt of jet that is not matched to electron
      //pass event selection     
      for(Int_t j=0; j<jetArr->GetEntries(); j++) {
        const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
        if (jet->pt > leadingJetPt &&
            mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) > 1.0) {
          leadingJetPt = jet->pt;          
        }
      }

      //********************************************************
      // Photons
      //********************************************************
      Int_t NPhotons = 0;      
      Double_t photonEt = -1;
      for(Int_t p=0; p<photonArr->GetEntries(); p++) {
        const mithep::TPhoton *photon = (mithep::TPhoton*)((*photonArr)[p]);
        
        Bool_t isEle = kFALSE;
        for(Int_t e=0; e<electronArr->GetEntries(); e++) {
          const mithep::TElectron *tmpEle = (mithep::TElectron*)((*electronArr)[e]);   
          if (mithep::MathUtils::DeltaR(photon->phi, photon->eta, tmpEle->phi, tmpEle->eta) < 0.3) isEle = kTRUE;
        }
 


        //photon ID
        if ( photon->pt > 20
             && !isEle
             && !photon->hasPixelSeed
             && photon->emIso04 < 2.0 + 0.006*photon->pt
             && photon->hadIso04 < 2.0+0.0025*photon->pt
             && photon->trkIso04 < 1.5 + 0.001*photon->pt
             && ( (fabs(photon->eta) < 1.5 && photon->sigiEtaiEta < 0.01) || (fabs(photon->eta) >= 1.5 && photon->sigiEtaiEta < 0.028) )
          ) {
          continue;
        }


        mithep::FourVectorM phFourVector;
        mithep::FourVectorM eleFourVector;
        eleFourVector.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
        phFourVector.SetCoordinates(photon->pt, photon->eta, photon->phi, 0.0 );
        mithep::FourVectorM dilepton = phFourVector+eleFourVector;
        
//         cout << "photon " << p << " : " << photon->et << " " << photon->eta << " " << photon->phi 
//              << "Ele : " << ele->pt << " " << ele->eta << " " << ele->phi << " : " 
//              << dilepton.M() << endl;
        
        if ( fabs(dilepton.M() - 91) > 20 ) { // if not in the Z mass peak
          photonEt = photon->pt;
          NPhotons++;
        }
      }
      
      
      for ( UInt_t denominatorTypeIndex = 0 ; denominatorTypeIndex < denominatorType.size() ; ++denominatorTypeIndex ) {
        for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
          for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
          
            //********************************************************
            // Event Selection Cuts
            //********************************************************

            if (NElectrons > 1) continue;

            if (sampleLabel[sampleTypeIndex] == "Ele8Sample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLSample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLJet40Sample") {
              if (met > 20) continue;

              Bool_t passJetSelection = kFALSE;
              if (ptThreshold[ptThresholdIndex] == 0) passJetSelection = kTRUE;
              if (leadingJetPt > ptThreshold[ptThresholdIndex]) {
                passJetSelection = kTRUE;               
              }

              if (!passJetSelection) continue;
            }
      
            if (sampleLabel[sampleTypeIndex] == "PhotonJetsSample") {
              if (NPhotons != 1) continue;
              if (!(photonEt > ptThreshold[ptThresholdIndex])) continue;  // if photon Et not high enough
            }

            if (passElectronDenominatorCuts(info->triggerBits, ele, denominatorType[denominatorTypeIndex], sampleLabel[sampleTypeIndex])) {
              

              eventListFile << info->runNum << " " << info->lumiSec << " " << info->evtNum << " : " << ele->pt << " " << ele->eta << " " << ele->phi << " : " << passElectronNumeratorCuts(ele) << endl;
              

              DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt);
              DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->eta);
              DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi);
              DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, ele->eta);


               if (
                (IsoChoice == 0 && passElectronNumeratorCuts(ele))
                ) {
                NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt);
                NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->eta);
                NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi);
                NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, ele->eta);   
              }
            }

          } //loop over denominator types
        } //loop over sample types
      } //loop over ptThresholds

    } //loop over electrons

  } //end loop over data     

  eventListFile.close();

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;
  
  
  //*****************************************************************************************
  //Make Efficiency Plots
  //*****************************************************************************************
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {

        
        Int_t ErrorType = 2; //Clopper Pearson errors
        TGraphAsymmErrors *efficiency_pt = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt", ptbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_eta = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Eta", etabins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_phi = mithep::EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Phi", phibins, ErrorType, -99, -99, 0, 1);
        mithep::TH2DAsymErr *efficiency_PtEta = 
          mithep::EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", 
                                                          ptbins2D, etabins2D, ErrorType);
        
        TFile *file = new TFile(outputFilename.c_str(), "UPDATE");
        file->cd();
        
        file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_PtEta, efficiency_PtEta->GetName(), "WriteDelete");
        file->WriteTObject(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        

        file->Close();
        

      }
    }
  }
    
  gBenchmark->Show("WWTemplate");       
} 


Bool_t passElectronNumeratorCuts(const mithep::TElectron *ele) {
  
  Bool_t pass = kTRUE;
//   if (ele->pt < 15) pass = kFALSE;
  if (fabs(ele->eta) >= 2.5) pass = kFALSE;

  //Barrel 
  if (fabs(ele->scEta) < 1.479) {
    if (! ( (0==0)
            && ele->sigiEtaiEta < 0.01 
            && fabs(ele->deltaEtaIn) < 0.004
            && fabs(ele->deltaPhiIn) < 0.06
            && ele->HoverE < 0.04
            //&& (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
	    && ele->pfIso04 / ele->pt < 0.13
            && ele->nExpHitsInner <= 0
            && !ele->isConv
            && fabs(ele->d0) < 0.02
            && fabs(ele->dz) < 0.1
          )
      ) {
      pass = kFALSE;
    }      
  }
  //Endcap
  else {
    if (! (  (0==0)
             && ele->sigiEtaiEta < 0.03
             && fabs(ele->deltaEtaIn) < 0.007
             && fabs(ele->deltaPhiIn) < 0.03
             //&& (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
             && ele->pfIso04 / ele->pt < 0.09
             && ele->nExpHitsInner <= 0
             && !ele->isConv
             && fabs(ele->d0) < 0.02
             && fabs(ele->dz) < 0.1
          )
      ) {
      pass = kFALSE;
    }
  } 
  if (ele->pt < 20) {
    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.03
              && ele->HoverE < 0.025
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && fabs(ele->deltaEtaIn) < 0.005
               && fabs(ele->deltaPhiIn) < 0.02
            )
        ) {
        pass = kFALSE;
      }
    } 

    if (ele->fBrem <= 0.15) {
      if (fabs(ele->eta) > 1.0) {
        pass = kFALSE;
      } else {
        if (!( ele->EoverP > 0.95 )) pass = kFALSE;
      }
    }
  }



  return pass;
}


Bool_t passElectronDenominatorCuts(Int_t triggerBits, const mithep::TElectron *ele, Int_t DenominatorType, string SampleType) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) < 2.5)) pass = kFALSE;

  //match to HLT
  if (SampleType == "Ele8Sample") {
    if (!(triggerBits & kHLT_Ele8 
          && ele->hltMatchBits & kHLTObject_Ele8)) pass = kFALSE;
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLSample") {
    if (!( (triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL 
            && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL
             )
//            || 
//            (triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL
//             && ele->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL
//              )           
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele17CaloIdLCaloIsoVLSample") {
    if (!( 
          (triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL
           && ele->hltMatchBits & kHLTObject_Ele17_CaloIdL_CaloIsoVL
            )           
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLJet40Sample") {
    if (!( (triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40)
          && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "PhotonJetsSample") {
    if (!(
          (triggerBits & kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL)
          && ele->hltMatchBits & kHLTObject_Ele8_CaloIdL_CaloIsoVL)
      ) pass = kFALSE;
  }

  if (DenominatorType == 1) {
    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.007
              && fabs(ele->deltaPhiIn) < 0.15
              && ele->HoverE < 0.12
              && ele->nExpHitsInner <= 0
              && !ele->isConv
              && fabs(ele->dz) < 0.1
//               && (ele->trkIso03 ) / ele->pt < 0.2
//               && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.2
//               && (ele->hadIso03) / ele->pt < 0.2
           )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.009
               && fabs(ele->deltaPhiIn) < 0.10
               && ele->nExpHitsInner <= 0
               && !ele->isConv
               && fabs(ele->dz) < 0.1
//                && (ele->trkIso03 ) / ele->pt < 0.2
//                && (ele->emIso03 ) / ele->pt < 0.2
//                && (ele->hadIso03) / ele->pt < 0.2
            )
        ) {
        pass = kFALSE;
      }
    } 
  }

  if (DenominatorType == 2) {
    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.007
              && fabs(ele->deltaPhiIn) < 0.15
              && ele->HoverE < 0.12
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
              && ele->nExpHitsInner <= 0
              && !ele->isConv
              && fabs(ele->dz) < 0.1
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.009
               && fabs(ele->deltaPhiIn) < 0.10
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
               && ele->nExpHitsInner <= 0
               && !ele->isConv
               && fabs(ele->dz) < 0.1
            )
        ) {
        pass = kFALSE;
      }
    } 
  }


  if (DenominatorType == 3) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.004
              && fabs(ele->deltaPhiIn) < 0.06
              && ele->HoverE < 0.04
              && ele->nExpHitsInner <= 0
              && !ele->isConv
              && fabs(ele->d0) < 0.02
              && fabs(ele->dz) < 0.1
//               && (ele->trkIso03 ) / ele->pt < 0.2
//               && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.2
//               && (ele->hadIso03) / ele->pt < 0.2
           )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.007
               && fabs(ele->deltaPhiIn) < 0.03
               && ele->nExpHitsInner <= 0
               && !ele->isConv
               && fabs(ele->d0) < 0.02
               && fabs(ele->dz) < 0.1
 //               && (ele->trkIso03 ) / ele->pt < 0.2
//                && (ele->emIso03 ) / ele->pt < 0.2
//                && (ele->hadIso03) / ele->pt < 0.2
            )
        ) {
        pass = kFALSE;
      }
    } 

  }

  if (DenominatorType == 4) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.007
               && fabs(ele->deltaPhiIn) < 0.15
//               && fabs(ele->deltaPhiIn) < 0.10
              && ele->HoverE < 0.12
              && ele->nExpHitsInner <= 0
              && !ele->isConv
              && fabs(ele->dz) < 0.1
              && (ele->trkIso03 ) / ele->pt < 0.2
              && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.20
              && (ele->hadIso03) / ele->pt < 0.20
              
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->sigiEtaiEta < 0.03
               && fabs(ele->deltaEtaIn) < 0.009
               && fabs(ele->deltaPhiIn) < 0.10
               && ele->nExpHitsInner <= 0
               && !ele->isConv
               && fabs(ele->dz) < 0.1
               && (ele->trkIso03) / ele->pt < 0.2
               && (ele->emIso03) / ele->pt < 0.20
               && (ele->hadIso03) / ele->pt < 0.20
           )
        ) {
        pass = kFALSE;
      }
    } 

  }



  if (DenominatorType == 5) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->nExpHitsInner <= 0
              && !ele->isConv
              && fabs(ele->dz) < 0.1
              && (ele->trkIso03 + TMath::Max(ele->emIso03 - 1.0, 0.0) + ele->hadIso03) / ele->pt < 0.1
              && (ele->trkIso03 ) / ele->pt < 0.2
              && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.2
              && (ele->hadIso03) / ele->pt < 0.2

            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else {
      if (! (  (0==0)
               && ele->nExpHitsInner <= 0
               && !ele->isConv
               && fabs(ele->dz) < 0.1
               && (ele->trkIso03 + ele->emIso03  + ele->hadIso03) / ele->pt < 0.1
               && (ele->trkIso03 ) / ele->pt < 0.2
               && (ele->emIso03 ) / ele->pt < 0.2
               && (ele->hadIso03) / ele->pt < 0.2
            )
        ) {
        pass = kFALSE;
      }
    } 

  }


  if (DenominatorType == 6) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->nExpHitsInner <= 0
              && !ele->isConv
              && fabs(ele->dz) < 0.1
              && (ele->trkIso03 ) / ele->pt < 0.2
              && (TMath::Max(ele->emIso03 - 1.0, 0.0)) / ele->pt < 0.2
              && (ele->hadIso03) / ele->pt < 0.2
            )
        ) {
        pass = kFALSE;
      }      
    }
    //Endcap
    else  {
      if (! (  (0==0)
               && ele->nExpHitsInner <= 0
               && !ele->isConv
               && fabs(ele->dz) < 0.1
               && (ele->trkIso03 ) / ele->pt < 0.2
               && (ele->emIso03 ) / ele->pt < 0.2
               && (ele->hadIso03) / ele->pt < 0.2
           )
        ) {
        pass = kFALSE;
      }
    } 

  }

  if (DenominatorType == 100) {
    pass = kTRUE;
  }


  return pass;
}


//--------------------------------------------------------------------------------------------------
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2)
{
  ofs <<   runNum << " " ;
  ofs <<  lumiSec << " ";
  ofs << evtNum<< " ";
  ofs << mass<< " ";

//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
//   ofs << "    pt    |    eta    |    phi    |   iso    |    d0      | ntk | npx | nseg | nval | chi^2/ndf | TM | HLT" << endl;
//   ofs << "----------+-----------+-----------+----------+------------+-----+-----+------+------+-----------+----+------" << endl;
  ofs << " " ;
  ofs << setw(9) << pt1 << " |";
  ofs << setw(10) << eta1 << " |";
  ofs << setw(10) << phi1 << " |";
  ofs << setw(10) << leptonCharge1 << " |";
  ofs << setw(9) << pt2 << " |";
  ofs << setw(10) << eta2 << " |";
  ofs << setw(10) << phi2 << " |";
  ofs << setw(10) << leptonCharge2 << " |";
  ofs << endl;
  
}
