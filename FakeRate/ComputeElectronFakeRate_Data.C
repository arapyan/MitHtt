//root -l -b -q MitHtt/FakeRate/ComputeElectronFakeRate_Data.C+\(\)
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
//#include <MitStyle.h>

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TPhoton.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Utils/LeptonIDCuts.hh"

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHtt/Utils/EfficiencyUtils.hh"
#include "MitHtt/Utils/PlotUtils.hh"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"

#endif

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
                               const string outputFilename,  
                               Int_t Option = 0);

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
      cout << "Cannot get Directory HwwNtuplerMod from file " << infname << endl;
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

  DoComputeElectronFakeRate("Full2011","ElectronFakeRate","ElectronFakeRate.root",21);

}



void DoComputeElectronFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename, 
                               Int_t Option)
{  
  gBenchmark->Start("WWTemplate");

  Double_t LUMI = 2163;

  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

  mithep::ElectronIDMVA *electronIDMVANoIPInfo = 0;

  if (Option >= 20 && Option <= 29) {
    electronIDMVANoIPInfo = new mithep::ElectronIDMVA();
    electronIDMVANoIPInfo->Initialize("BDTG method",
                                      "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet0LowPt_NoIPInfo_BDTG.weights.xml", 
                                      "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet1LowPt_NoIPInfo_BDTG.weights.xml", 
                                      "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet2LowPt_NoIPInfo_BDTG.weights.xml", 
                                      "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet0HighPt_NoIPInfo_BDTG.weights.xml", 
                                      "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet1HighPt_NoIPInfo_BDTG.weights.xml", 
                                      "/home/vdutta/cms/cmssw/023_2/CMSSW_4_2_4_patch1/src/MitPhysics/data/ElectronMVAWeights/Subdet2HighPt_NoIPInfo_BDTG.weights.xml",
                                      mithep::ElectronIDMVA::kNoIPInfo );
    
  }


  //*****************************************************************************************
  //Define Pt bins
  //*****************************************************************************************
  vector<double> ptbins;

  ptbins.push_back(10);  
  ptbins.push_back(15);  
  ptbins.push_back(20);  
  ptbins.push_back(25);  
  ptbins.push_back(30);  
  ptbins.push_back(35);  
  ptbins.push_back(40);  
  ptbins.push_back(50);  
  ptbins.push_back(100);  


  vector<double> etabins;
  etabins.push_back(0.0);
  etabins.push_back(0.5);
  etabins.push_back(1.0);
  etabins.push_back(1.5);
  etabins.push_back(2.0);
  etabins.push_back(2.5);
  etabins.push_back(3.0);


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

  vector<double> nvtxbins;
  nvtxbins.push_back(0);
  nvtxbins.push_back(2);
  nvtxbins.push_back(4);
  nvtxbins.push_back(6);
  nvtxbins.push_back(8);
  nvtxbins.push_back(10);
  nvtxbins.push_back(12);
  nvtxbins.push_back(14);
  nvtxbins.push_back(16);
  nvtxbins.push_back(18);
  nvtxbins.push_back(20);
  nvtxbins.push_back(25);


  vector<double> rhobins;
  rhobins.push_back(0);
  rhobins.push_back(2);
  rhobins.push_back(4);
  rhobins.push_back(6);
  rhobins.push_back(8);
  rhobins.push_back(10);
  rhobins.push_back(12);
  rhobins.push_back(14);
  rhobins.push_back(16);
  rhobins.push_back(18);



  vector<double> ptbins2D;
  ptbins2D.push_back(10);  
  ptbins2D.push_back(15);  
  ptbins2D.push_back(20);  
  ptbins2D.push_back(25);  
  ptbins2D.push_back(30);  
  ptbins2D.push_back(35);  
  vector<double> etabins2D;
  etabins2D.push_back(0.0);
  etabins2D.push_back(1.0);
  etabins2D.push_back(1.479);
  etabins2D.push_back(2.0);
  etabins2D.push_back(2.5);


  //--------------------------------------------------------------------------------------------------------------
  // Histograms
  //==============================================================================================================  
  TH1F *histLeadingJetPt = new TH1F("leadingJetPt" , "; p_{T} [GeV/c] ; Number of Events ",  200, 0 , 200);

  //3D array, indices give: [denominatorType][SampleType][ptThreshold]

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt;
  vector<vector<vector<TH1F*> > > DenominatorVector_Eta;
  vector<vector<vector<TH1F*> > > DenominatorVector_Phi;
  vector<vector<vector<TH1F*> > > DenominatorVector_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Rho;
  vector<vector<vector<TH2F*> > > DenominatorVector_PtEta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt;
  vector<vector<vector<TH1F*> > > NumeratorVector_Eta;
  vector<vector<vector<TH1F*> > > NumeratorVector_Phi;
  vector<vector<vector<TH1F*> > > NumeratorVector_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Rho;
  vector<vector<vector<TH2F*> > > NumeratorVector_PtEta;
  vector<vector<vector<TH1F*> > > LeptonJetPt;
  vector<vector<vector<TH1F*> > > DenominatorIsolation;

  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > DenominatorVector_Pt20ToInf_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt10To20_Endcap_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Barrel_Rho;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_NVtx;
  vector<vector<vector<TH1F*> > > NumeratorVector_Pt20ToInf_Endcap_Rho;


  vector<Double_t> denominatorType;
  denominatorType.push_back(4);
  vector<string> sampleLabel;
//   sampleLabel.push_back("Ele8Sample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLSample");
  sampleLabel.push_back("Ele17CaloIdLCaloIsoVLSample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLJet40Sample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLPtCombinedSample");
  sampleLabel.push_back("Ele8CaloIdLCaloIsoVLCombinedSample");
//    sampleLabel.push_back("PhotonJetsSample");
  vector<Double_t> ptThreshold;
  ptThreshold.push_back(15);
  ptThreshold.push_back(20);
  ptThreshold.push_back(25);
  ptThreshold.push_back(30);
  ptThreshold.push_back(35);
  ptThreshold.push_back(40);
  ptThreshold.push_back(45);
  ptThreshold.push_back(50);
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    vector<vector<TH1F*> > tmpDenominatorVector_Pt;
    vector<vector<TH1F*> > tmpDenominatorVector_Eta;
    vector<vector<TH1F*> > tmpDenominatorVector_Phi;
    vector<vector<TH1F*> > tmpDenominatorVector_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Rho;
    vector<vector<TH2F*> > tmpDenominatorVector_PtEta;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt;
    vector<vector<TH1F*> > tmpNumeratorVector_Eta;
    vector<vector<TH1F*> > tmpNumeratorVector_Phi;
    vector<vector<TH1F*> > tmpNumeratorVector_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Rho;
    vector<vector<TH2F*> > tmpNumeratorVector_PtEta;
    vector<vector<TH1F*> > tmpLeptonJetPt;
    vector<vector<TH1F*> > tmpDenominatorIsolation;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpDenominatorVector_Pt20ToInf_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt10To20_Endcap_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Barrel_Rho;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
    vector<vector<TH1F*> > tmpNumeratorVector_Pt20ToInf_Endcap_Rho;


    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      vector<TH1F*>  tmptmpDenominatorVector_Pt;
      vector<TH1F*>  tmptmpDenominatorVector_Eta;
      vector<TH1F*>  tmptmpDenominatorVector_Phi;
      vector<TH1F*>  tmptmpDenominatorVector_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Rho;
      vector<TH2F*>  tmptmpDenominatorVector_PtEta;
      vector<TH1F*>  tmptmpNumeratorVector_Pt;
      vector<TH1F*>  tmptmpNumeratorVector_Eta;
      vector<TH1F*>  tmptmpNumeratorVector_Phi;
      vector<TH1F*>  tmptmpNumeratorVector_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Rho;
      vector<TH2F*>  tmptmpNumeratorVector_PtEta;
      vector<TH1F*>  tmptmpLeptonJetPt;
      vector<TH1F*>  tmptmpDenominatorIsolation;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt10To20_Endcap_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx;
      vector<TH1F*>  tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho;
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {
        TH1F *histDenominator_Pt = new TH1F(("histDenominator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str(), "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominator_Eta = new TH1F(("histDenominator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0 , 3.0);
        TH1F *histDenominator_Phi = new TH1F(("histDenominator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histDenominator_NVtx = new TH1F(("histDenominator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Rho = new TH1F(("histDenominator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histDenominator_PtEta = new TH2F(("histDenominator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0, 3.0);
        TH1F *histNumerator_Pt = new TH1F(("histNumerator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histNumerator_Eta = new TH1F(("histNumerator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0.0 , 3.0);
        TH1F *histNumerator_Phi = new TH1F(("histNumerator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histNumerator_NVtx = new TH1F(("histNumerator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Rho = new TH1F(("histNumerator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histNumerator_PtEta = new TH2F(("histNumerator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0, 3.0);        
        TH1F *histLeptonJetPt = new TH1F(("histLeptonJetPt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histDenominatorIsolation = new TH1F(("histDenominatorIsolation_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; PF RelIso ; Number of Events ",  100, 0 , 1.0);

        TH1F *histDenominator_Pt10To20_Barrel_NVtx = new TH1F(("histDenominator_Pt10To20_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Barrel_Rho = new TH1F(("histDenominator_Pt10To20_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt10To20_Endcap_NVtx = new TH1F(("histDenominator_Pt10To20_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt10To20_Endcap_Rho = new TH1F(("histDenominator_Pt10To20_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Barrel_NVtx = new TH1F(("histDenominator_Pt20ToInf_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Barrel_Rho = new TH1F(("histDenominator_Pt20ToInf_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histDenominator_Pt20ToInf_Endcap_NVtx = new TH1F(("histDenominator_Pt20ToInf_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Pt20ToInf_Endcap_Rho = new TH1F(("histDenominator_Pt20ToInf_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Barrel_NVtx = new TH1F(("histNumerator_Pt10To20_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Barrel_Rho = new TH1F(("histNumerator_Pt10To20_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt10To20_Endcap_NVtx = new TH1F(("histNumerator_Pt10To20_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt10To20_Endcap_Rho = new TH1F(("histNumerator_Pt10To20_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Barrel_NVtx = new TH1F(("histNumerator_Pt20ToInf_Barrel_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Barrel_Rho = new TH1F(("histNumerator_Pt20ToInf_Barrel_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH1F *histNumerator_Pt20ToInf_Endcap_NVtx = new TH1F(("histNumerator_Pt20ToInf_Endcap_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Pt20ToInf_Endcap_Rho = new TH1F(("histNumerator_Pt20ToInf_Endcap_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);



        tmptmpDenominatorVector_Pt.push_back(histDenominator_Pt);
        tmptmpDenominatorVector_Eta.push_back(histDenominator_Eta);
        tmptmpDenominatorVector_Phi.push_back(histDenominator_Phi);
        tmptmpDenominatorVector_NVtx.push_back(histDenominator_NVtx);
        tmptmpDenominatorVector_Rho.push_back(histDenominator_Rho);
        tmptmpDenominatorVector_PtEta.push_back(histDenominator_PtEta);
        tmptmpNumeratorVector_Pt.push_back(histNumerator_Pt);
        tmptmpNumeratorVector_Eta.push_back(histNumerator_Eta);
        tmptmpNumeratorVector_Phi.push_back(histNumerator_Phi);
        tmptmpNumeratorVector_NVtx.push_back(histNumerator_NVtx);
        tmptmpNumeratorVector_Rho.push_back(histNumerator_Rho);
        tmptmpNumeratorVector_PtEta.push_back(histNumerator_PtEta);
        tmptmpLeptonJetPt.push_back(histLeptonJetPt);
        tmptmpDenominatorIsolation.push_back(histDenominatorIsolation);
        tmptmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(histDenominator_Pt10To20_Barrel_NVtx);
        tmptmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(histDenominator_Pt10To20_Barrel_Rho);
        tmptmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(histDenominator_Pt10To20_Endcap_NVtx);
        tmptmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(histDenominator_Pt10To20_Endcap_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(histDenominator_Pt20ToInf_Barrel_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(histDenominator_Pt20ToInf_Barrel_Rho);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(histDenominator_Pt20ToInf_Endcap_NVtx);
        tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(histDenominator_Pt20ToInf_Endcap_Rho);
        tmptmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(histNumerator_Pt10To20_Barrel_NVtx);
        tmptmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(histNumerator_Pt10To20_Barrel_Rho);
        tmptmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(histNumerator_Pt10To20_Endcap_NVtx);
        tmptmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(histNumerator_Pt10To20_Endcap_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(histNumerator_Pt20ToInf_Barrel_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(histNumerator_Pt20ToInf_Barrel_Rho);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(histNumerator_Pt20ToInf_Endcap_NVtx);
        tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(histNumerator_Pt20ToInf_Endcap_Rho);

      }
      tmpDenominatorVector_Pt.push_back(tmptmpDenominatorVector_Pt);
      tmpDenominatorVector_Eta.push_back(tmptmpDenominatorVector_Eta);
      tmpDenominatorVector_Phi.push_back(tmptmpDenominatorVector_Phi);
      tmpDenominatorVector_NVtx.push_back(tmptmpDenominatorVector_NVtx);
      tmpDenominatorVector_Rho.push_back(tmptmpDenominatorVector_Rho);
      tmpDenominatorVector_PtEta.push_back(tmptmpDenominatorVector_PtEta);
      tmpNumeratorVector_Pt.push_back(tmptmpNumeratorVector_Pt);
      tmpNumeratorVector_Eta.push_back(tmptmpNumeratorVector_Eta);
      tmpNumeratorVector_Phi.push_back(tmptmpNumeratorVector_Phi);
      tmpNumeratorVector_NVtx.push_back(tmptmpNumeratorVector_NVtx);
      tmpNumeratorVector_Rho.push_back(tmptmpNumeratorVector_Rho);
      tmpNumeratorVector_PtEta.push_back(tmptmpNumeratorVector_PtEta);
      tmpLeptonJetPt.push_back(tmptmpLeptonJetPt);
      tmpDenominatorIsolation.push_back(tmptmpDenominatorIsolation);
      tmpDenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_NVtx);
      tmpDenominatorVector_Pt10To20_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Barrel_Rho);
      tmpDenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_NVtx);
      tmpDenominatorVector_Pt10To20_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt10To20_Endcap_Rho);
      tmpDenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
      tmpDenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Barrel_Rho);
      tmpDenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
      tmpDenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpDenominatorVector_Pt20ToInf_Endcap_Rho);
      tmpNumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_NVtx);
      tmpNumeratorVector_Pt10To20_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Barrel_Rho);
      tmpNumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_NVtx);
      tmpNumeratorVector_Pt10To20_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt10To20_Endcap_Rho);
      tmpNumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
      tmpNumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Barrel_Rho);
      tmpNumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
      tmpNumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmptmpNumeratorVector_Pt20ToInf_Endcap_Rho);
    }
    DenominatorVector_Pt.push_back(tmpDenominatorVector_Pt);
    DenominatorVector_Eta.push_back(tmpDenominatorVector_Eta);
    DenominatorVector_Phi.push_back(tmpDenominatorVector_Phi);
    DenominatorVector_NVtx.push_back(tmpDenominatorVector_NVtx);
    DenominatorVector_Rho.push_back(tmpDenominatorVector_Rho);
    DenominatorVector_PtEta.push_back(tmpDenominatorVector_PtEta);
    NumeratorVector_Pt.push_back(tmpNumeratorVector_Pt);
    NumeratorVector_Eta.push_back(tmpNumeratorVector_Eta);
    NumeratorVector_Phi.push_back(tmpNumeratorVector_Phi);
    NumeratorVector_NVtx.push_back(tmpNumeratorVector_NVtx);
    NumeratorVector_Rho.push_back(tmpNumeratorVector_Rho);
    NumeratorVector_PtEta.push_back(tmpNumeratorVector_PtEta);
    LeptonJetPt.push_back(tmpLeptonJetPt);
    DenominatorIsolation.push_back(tmpDenominatorIsolation);
    DenominatorVector_Pt10To20_Barrel_NVtx.push_back(tmpDenominatorVector_Pt10To20_Barrel_NVtx);
    DenominatorVector_Pt10To20_Barrel_Rho.push_back(tmpDenominatorVector_Pt10To20_Barrel_Rho);
    DenominatorVector_Pt10To20_Endcap_NVtx.push_back(tmpDenominatorVector_Pt10To20_Endcap_NVtx);
    DenominatorVector_Pt10To20_Endcap_Rho.push_back(tmpDenominatorVector_Pt10To20_Endcap_Rho);
    DenominatorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_NVtx);
    DenominatorVector_Pt20ToInf_Barrel_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Barrel_Rho);
    DenominatorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_NVtx);
    DenominatorVector_Pt20ToInf_Endcap_Rho.push_back(tmpDenominatorVector_Pt20ToInf_Endcap_Rho);
    NumeratorVector_Pt10To20_Barrel_NVtx.push_back(tmpNumeratorVector_Pt10To20_Barrel_NVtx);
    NumeratorVector_Pt10To20_Barrel_Rho.push_back(tmpNumeratorVector_Pt10To20_Barrel_Rho);
    NumeratorVector_Pt10To20_Endcap_NVtx.push_back(tmpNumeratorVector_Pt10To20_Endcap_NVtx);
    NumeratorVector_Pt10To20_Endcap_Rho.push_back(tmpNumeratorVector_Pt10To20_Endcap_Rho);
    NumeratorVector_Pt20ToInf_Barrel_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_NVtx);
    NumeratorVector_Pt20ToInf_Barrel_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Barrel_Rho);
    NumeratorVector_Pt20ToInf_Endcap_NVtx.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_NVtx);
    NumeratorVector_Pt20ToInf_Endcap_Rho.push_back(tmpNumeratorVector_Pt20ToInf_Endcap_Rho);
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
  TClonesArray *pvArr = new TClonesArray("mithep::TVertex");
  TClonesArray *photonArr = new TClonesArray("mithep::TPhoton");
  
  //********************************************************
  // Good RunLumi Selection
  //********************************************************
  Bool_t hasJSON = kTRUE;
  mithep::RunLumiRangeMap rlrm;
  rlrm.AddJSONFile("/home/vdutta/cms/root/json/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v3.txt");
  rlrm.AddJSONFile("/home/vdutta/cms/root/json/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_JSON_v2.txt");
  rlrm.AddJSONFile("/home/vdutta/cms/root/json/Cert_160404-180252_7TeV_PromptReco_Collisions11_JSON.txt");

  Int_t NEvents = 0;

  vector<string> inputfiles;
  if (inputFilename == "Full2011") {
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-del-m10-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-del-pr-v4_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-del-a05-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-del-o03-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11b-del-pr-v1_fakerate_skim.root");
  } else if (inputFilename == "RUN2011A") {
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-del-m10-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-del-pr-v4_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-del-a05-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-del-o03-v1_fakerate_skim.root");
  }  else if (inputFilename == "RUN2011B") {
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11b-del-pr-v1_fakerate_skim.root");
  }
  else {
    inputfiles.push_back(inputFilename);
  }

  for (UInt_t f = 0; f < inputfiles.size(); ++f) {

    //********************************************************
    // Get Tree
    //********************************************************
    eventTree = getTreeFromFile(inputfiles[f].c_str(),"Events"); 
    TBranch *infoBr;
    TBranch *electronBr;
    TBranch *muonBr;
    TBranch *pvBr;
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
    eventTree->SetBranchAddress("PV", &pvArr);             pvBr = eventTree->GetBranch("PV");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry < eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;
      
      //For MVA, only use odd event numbers because even event numbers were used for training
      if (Option >= 20 && Option < 30) {
        if (info->evtNum % 2 == 0) continue;
      }

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...

      NEvents++;

      //********************************************************
      // Load the branches
      //********************************************************
      electronArr->Clear(); 
      muonArr->Clear(); 
      photonArr->Clear(); 
      pvArr->Clear();
      jetArr->Clear(); 
      electronBr->GetEntry(ientry);
      muonBr->GetEntry(ientry);
      photonBr->GetEntry(ientry);
      pvBr->GetEntry(ientry);
      jetBr->GetEntry(ientry);


      //********************************************************
      // Met
      //********************************************************
      Double_t met = info->pfMET;

      Int_t NElectrons = electronArr->GetEntries();
      
      for(Int_t i=0; i<electronArr->GetEntries(); i++) {
        const mithep::TElectron *ele = (mithep::TElectron*)((*electronArr)[i]);
      

        //likelihood value
        mithep::FourVectorM tmpEleSC;
        tmpEleSC.SetCoordinates(ele->scEt, ele->scEta, ele->scPhi, 0.51099892e-3 );
      

        Double_t leadingJetPt = -1;
        Double_t leptonJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) > 1.0) {
            leadingJetPt = jet->pt;          
          }
          if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, ele->phi, ele->eta) < 0.5) {
            leptonJetPt = jet->pt;
          }
        }
      
        //if there's no jet ( == isolated lepton?) then take pt of lepton
        if (leptonJetPt < 0) {
          leptonJetPt = ele->pt;

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
          if ( photon->scEt > 20
               && !isEle
               && !photon->hasPixelSeed
               && photon->emIso04 < 2.0 + 0.006*photon->scEt
               && photon->hadIso04 < 2.0+0.0025*photon->scEt
               && photon->trkIso04 < 1.5 + 0.001*photon->scEt
               && ( (fabs(photon->eta) < 1.5 && photon->sigiEtaiEta < 0.01) || (fabs(photon->eta) >= 1.5 && photon->sigiEtaiEta < 0.028) )
            ) {
            continue;
          }


          mithep::FourVectorM phFourVector;
          mithep::FourVectorM eleFourVector;
          eleFourVector.SetCoordinates(ele->pt, ele->eta, ele->phi, 0.51099892e-3 );
          phFourVector.SetCoordinates(photon->scEt, photon->eta, photon->phi, 0.0 );
          mithep::FourVectorM dilepton = phFourVector+eleFourVector;
        
          if ( fabs(dilepton.M() - 91) > 20 ) {
            photonEt = photon->scEt;
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

              if (sampleLabel[sampleTypeIndex] == "Ele8Sample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLSample" || sampleLabel[sampleTypeIndex] == "Ele17CaloIdLCaloIsoVLSample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLJet40Sample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLCombinedSample" || sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLPtCombinedSample") {
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
                if (!(photonEt > ptThreshold[ptThresholdIndex])) continue;  
              }

              if (passElectronDenominatorCuts(info->triggerBits, ele, denominatorType[denominatorTypeIndex], sampleLabel[sampleTypeIndex])) {
              

                DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt);
                DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta));
                DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi);
                DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, fabs(ele->eta));
              
                if (ele->pt < 20 && fabs(ele->eta) < 1.479) {
                  DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                } 
                if (ele->pt < 20 && fabs(ele->eta) >= 1.479) {
                  DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                } 
                if (ele->pt >= 20 && fabs(ele->eta) < 1.479) {
                  DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                } 
                if (ele->pt >= 20 && fabs(ele->eta) >= 1.479) {
                  DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                } 
              

		Double_t iso = eleIsoPU(ele)/ele->pt;
                LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(TMath::Min(TMath::Max(leptonJetPt, 0.01),99.9));
              
                DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(TMath::Min(TMath::Max(iso, 0.000001),0.9999999)); 
              
                if (
                  (Option == 0 && passElectronNumeratorCuts(ele))
                  ||
                  (Option == 21 && passEleIsoPU(ele) && passEleMVAID(ele, 
                                                   electronIDMVANoIPInfo->MVAValue(
                                                     ele->pt,ele->scEta,
                                                     ele->sigiEtaiEta, 
                                                     ele->deltaEtaIn,
                                                     ele->deltaPhiIn, 
                                                     ele->HoverE,
                                                     ele->d0,
                                                     ele->dz, 
                                                     ele->fBrem,
                                                     ele->EoverP,
                                                     ele->ESeedClusterOverPOut,
                                                     TMath::Sqrt(ele->sigiPhiiPhi),
                                                     ele->nBrem,
                                                     (1.0/(ele->scEt * TMath::CosH(ele->scEta)) - 1/ele->p), 
                                                     ele->ESeedClusterOverPIn,
                                                     ele->ip3d,
                                                     ele->ip3dSig )))
                  ) {
                  NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt);
                  NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(ele->eta));
                  NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->phi);
                  NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(ele->pt, fabs(ele->eta) );   
        
                  if (ele->pt < 20 && fabs(ele->eta) < 1.479) {
                    NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                    NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  }
                  if (ele->pt < 20 && fabs(ele->eta) >= 1.479) {
                    NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                    NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  } 
                  if (ele->pt >= 20 && fabs(ele->eta) < 1.479) {
                    NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                    NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  } 
                  if (ele->pt >= 20 && fabs(ele->eta) >= 1.479) {
                    NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                    NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  }
                }
              }

            } //loop over denominator types
          } //loop over sample types
        } //loop over ptThresholds

      } //loop over electrons

    } //end loop over data     
  } //end loop over files

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

        Bool_t printDebug = kFALSE;
        if (denominatorType[denominatorTypeIndex] == 4 
            // && sampleLabel[sampleTypeIndex] == "Ele8CaloIdLCaloIsoVLSample" && ptThreshold[ptThresholdIndex] == 30
          ) printDebug = kTRUE;
        cout << label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta" << endl;
        
        Int_t ErrorType = 2; //Clopper Pearson errors
        TGraphAsymmErrors *efficiency_pt = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt", ptbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_eta = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Eta", etabins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_phi = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Phi", phibins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_nvtx = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_rho = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        mithep::TH2DAsymErr *efficiency_PtEta = 
          EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", ptbins2D, etabins2D, ErrorType, printDebug);
        

        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_nvtx = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Barrel_rho = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Barrel_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_nvtx = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt10To20_Endcap_rho = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt10To20_Endcap_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_nvtx = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Barrel_rho = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Barrel_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_nvtx = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_Pt20ToInf_Endcap_rho = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt20ToInf_Endcap_Rho", rhobins, ErrorType, -99, -99, 0, 1);



        TFile *file = new TFile(outputFilename.c_str(), "UPDATE");
        file->cd();
        
        file->WriteTObject(efficiency_pt, efficiency_pt->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_eta, efficiency_eta->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_phi, efficiency_phi->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_nvtx, efficiency_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_PtEta, efficiency_PtEta->GetName(), "WriteDelete");
        file->WriteTObject(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        file->WriteTObject(DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->GetName(), "WriteDelete");
        
        file->WriteTObject(efficiency_Pt10To20_Barrel_nvtx, efficiency_Pt10To20_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Barrel_rho, efficiency_Pt10To20_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_nvtx, efficiency_Pt10To20_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt10To20_Endcap_rho, efficiency_Pt10To20_Endcap_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_nvtx, efficiency_Pt20ToInf_Barrel_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Barrel_rho, efficiency_Pt20ToInf_Barrel_rho->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_nvtx, efficiency_Pt20ToInf_Endcap_nvtx->GetName(), "WriteDelete");
        file->WriteTObject(efficiency_Pt20ToInf_Endcap_rho, efficiency_Pt20ToInf_Endcap_rho->GetName(), "WriteDelete");

        file->Close();
        
      }
    }
  }


  cout << "Total Events: " << NEvents << endl;

  gBenchmark->Show("WWTemplate");       
} 

Bool_t passElectronNumeratorCuts(const mithep::TElectron *ele) {

  if (fabs(ele->eta) >= 2.5) return kFALSE;

  if(fabs(ele->d0) > 0.02)   return kFALSE;
  if(fabs(ele->dz) > 0.1)    return kFALSE;

  // conversion rejection
  if(ele->nExpHitsInner > 0) return kFALSE;
  if(ele->isConv)            return kFALSE;

     
  // barrel/endcap dependent requirements      
  if(fabs(ele->scEta)<1.479) {
    //if(ele->pfIso04 > 0.13*(ele->pt)) return kFALSE;

    if(ele->pt>20) {
      if(ele->sigiEtaiEta      > 0.01)  return kFALSE;
      if(fabs(ele->deltaPhiIn) > 0.06)  return kFALSE;
      if(fabs(ele->deltaEtaIn) > 0.004) return kFALSE;
      if(ele->HoverE           > 0.04)  return kFALSE;

    } else {
      if(ele->sigiEtaiEta      > 0.01)  return kFALSE;
      if(fabs(ele->deltaPhiIn) > 0.03)  return kFALSE;
      if(fabs(ele->deltaEtaIn) > 0.004) return kFALSE;
      if(ele->HoverE           > 0.025) return kFALSE;

    }
  
  } else {
    //if(ele->pfIso04 > 0.09*(ele->pt)) return kFALSE;
     
    if(ele->pt>20) {
      if(ele->sigiEtaiEta      > 0.03)  return kFALSE;
      if(fabs(ele->deltaPhiIn) > 0.03)  return kFALSE;
      if(fabs(ele->deltaEtaIn) > 0.007) return kFALSE;
      if(ele->HoverE           > 0.10)  return kFALSE;

    } else {
      if(ele->sigiEtaiEta      > 0.03)  return kFALSE;
      if(fabs(ele->deltaPhiIn) > 0.02)  return kFALSE;
      if(fabs(ele->deltaEtaIn) > 0.005) return kFALSE;
      if(ele->HoverE           > 0.10)  return kFALSE;
     
    }
  }
  
  if(ele->pt < 20)
    return ((ele->fBrem>0.15) || (fabs(ele->eta)<1 && ele->EoverP>0.95));

  Double_t chargedIso = ele->pfIsoCharged;
  Double_t neutralIso = max(ele->pfIsoNeutral + ele->pfIsoGamma - 0.5 * ele->puIsoNoZ, 0.0);

  Double_t totalIso = chargedIso+neutralIso;

  //if(ele->pt > 20) {
    // barrel/endcap dependent requirments      
    if(fabs(ele->scEta)<1.479) {
      if(totalIso > 0.15*(ele->pt)) return kFALSE;
    } else {
      if(totalIso > 0.10*(ele->pt)) return kFALSE;
    }
  return kTRUE;
}

Bool_t passElectronDenominatorCuts(Int_t triggerBits, const mithep::TElectron *ele, Int_t DenominatorType, string SampleType) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(ele->pt > 10 && fabs(ele->eta) < 2.5)) pass = kFALSE;

  //match to HLT
  if (SampleType == "Ele8Sample") {
    if (!(triggerBits & kHLT_Ele8 
          && ele->hltMatchBits & kHLT_Ele8_EleObj)) pass = kFALSE;
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLSample") {
    if (!( (triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL 
            && ele->hltMatchBits & kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj
             )
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele17CaloIdLCaloIsoVLSample") {
    if (!( 
          (triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL
           && ele->hltMatchBits & kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj
            )
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLJet40Sample") {
    if (!( (triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40)
          && ele->hltMatchBits & kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLCombinedSample") {
    if (!(
          ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40)
           && ele->hltMatchBits & kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj
            )
          ||
          ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj
            )
          ||
          ((triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj
            )
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Ele8CaloIdLCaloIsoVLPtCombinedSample") {
    if (!(
          ((triggerBits & kHLT_Ele8_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj)
          ||
          ((triggerBits & kHLT_Ele17_CaloIdL_CaloIsoVL)
           && ele->hltMatchBits & kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj)
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "PhotonJetsSample") {
    if (!(
          (triggerBits & kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL)
          && ele->hltMatchBits & kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj)
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
              && ele->HoverE < 0.12
              && ele->nExpHitsInner <= 0
              && !ele->isConv
              && fabs(ele->dz) < 0.1
              && fabs(ele->d0) < 0.02
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
               && ele->HoverE < 0.10
               && ele->nExpHitsInner <= 0
               && !ele->isConv
               && fabs(ele->dz) < 0.1
               && fabs(ele->d0) < 0.02
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

//Boris Denominator
  if (DenominatorType == 20) {

    //Barrel 
    if (fabs(ele->scEta) < 1.479) {
      if (! ( (0==0)
              && ele->sigiEtaiEta < 0.01 
              && fabs(ele->deltaEtaIn) < 0.007
              && fabs(ele->deltaPhiIn) < 0.15
              && ele->HoverE < 0.12
              && (ele->trkIso03 ) / ele->pt < 0.2
              && (ele->emIso03 ) / ele->pt < 0.20
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
               && ele->HoverE < 0.10
               && (ele->trkIso03) / ele->pt < 0.2
               && (ele->emIso03) / ele->pt < 0.20
               && (ele->hadIso03) / ele->pt < 0.20
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
