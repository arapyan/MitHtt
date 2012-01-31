//root -l -b -q MitHtt/FakeRate/ComputeMuonFakeRate_Data.C+\(\)
//================================================================================================

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

// lumi section selection with JSON files
#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitHtt/Utils/EfficiencyUtils.hh"
#include "MitHtt/Utils/LeptonIDCuts.hh"

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
Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu);
Bool_t passMuonDenominatorCuts(Int_t triggerBits, const mithep::TMuon *mu, Int_t DenominatorType, string SampleType);
Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Int_t DenominatorType);
void eventDump(ofstream &ofs, const Int_t runNum, const Int_t lumiSec, const Int_t evtNum, Double_t mass,
               Double_t pt1, Double_t eta1, Double_t phi1, Int_t leptonCharge1, Double_t pt2, Double_t eta2, Double_t phi2, Int_t leptonCharge2);


//--------------------------------------------------------------------------------------------------
Bool_t passZVeto(TClonesArray *muonArr, Int_t DenominatorType)
{
  for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
    const mithep::TMuon* mu1 = (mithep::TMuon*)((*muonArr)[i]);
    if(mu1->pt        < 20)  continue;
    if(fabs(mu1->eta) > 2.4) continue;
    if(!passMuonDenominatorCuts(mu1, DenominatorType)) continue;
    for(Int_t j=i+1; j<muonArr->GetEntriesFast(); j++) {
      const mithep::TMuon* mu2 = (mithep::TMuon*)((*muonArr)[j]);
      if(mu1->q == mu2->q)     continue;
      if(mu2->pt	< 20)  continue;
      if(fabs(mu2->eta) > 2.4) continue;
      if(!passMuonDenominatorCuts(mu2, DenominatorType)) continue;

      return kFALSE;  // Z candidate => fail Z veto
    }
  }
  
  return kTRUE;  // No Z candidate => pass Z veto
}


//--------------------------------------------------------------------------------------------------
void DoComputeMuonFakeRate(const string inputFilename,
                           const string label, 
                           const string outputFilename);

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


//--------------------------------------------------------------------------------------------------
Double_t calcMt(const Double_t met, const Double_t metphi, const mithep::TMuon *muon)
{
  const Double_t m = 0.105659369;
  TLorentzVector vMuon; vMuon.SetPtEtaPhiM(muon->pt, muon->eta, muon->phi, m);
  TLorentzVector vMet;  vMet.SetPtEtaPhiM(met, 0, metphi, 0);
  Double_t et = (vMuon.E())*(vMuon.Pt())/(vMuon.P());
  
  return sqrt( (et+vMet.Perp())*(et+vMet.Perp()) - (vMuon.Px()+vMet.Px())*(vMuon.Px()+vMet.Px()) - (vMuon.Py()+vMet.Py())*(vMuon.Py()+vMet.Py()) );
}


//=== MAIN MACRO =================================================================================================
void ComputeMuonFakeRate_Data() {

  DoComputeMuonFakeRate("Full2011","MuonFakeRate","MuonFakeRate.root");

}



void DoComputeMuonFakeRate(const string inputFilename,
                               const string label, 
                               const string outputFilename)
{  
  gBenchmark->Start("WWTemplate");


  //*****************************************************************************************
  //Setup
  //*****************************************************************************************

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
  TH1F *leadingJetPt = new TH1F("leadingJetPt" , "; p_{T} [GeV/c] ; Number of Events ",  200, 0 , 200);

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
  //denominatorType.push_back(1);
  denominatorType.push_back(2);
  //denominatorType.push_back(3);
  vector<string> sampleLabel;
  //sampleLabel.push_back("Mu8Sample");
  //sampleLabel.push_back("Mu15Sample");
  sampleLabel.push_back("Mu8PtCombinedSample");
  //sampleLabel.push_back("Mu8Jet40Sample");
  //sampleLabel.push_back("PhotonJetsSample");
  vector<Double_t> ptThreshold;
  ptThreshold.push_back(0);
  ptThreshold.push_back(5);
  ptThreshold.push_back(10);
  ptThreshold.push_back(15);
  ptThreshold.push_back(20);
  ptThreshold.push_back(25);
  ptThreshold.push_back(30);
  
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
        TH1F *histDenominator_Eta = new TH1F(("histDenominator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0.0 , 3.0);
        TH1F *histDenominator_Phi = new TH1F(("histDenominator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histDenominator_NVtx = new TH1F(("histDenominator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histDenominator_Rho = new TH1F(("histDenominator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);

        TH2F *histDenominator_PtEta = new TH2F(("histDenominator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0.0, 3.0);
        TH1F *histNumerator_Pt = new TH1F(("histNumerator_Pt_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; Number of Events ",  100, 0 , 100);
        TH1F *histNumerator_Eta = new TH1F(("histNumerator_Eta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #eta ; Number of Events ",  100, 0.0 , 3.0);
        TH1F *histNumerator_Phi = new TH1F(("histNumerator_Phi_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #phi ; Number of Events ",  100, -3.2 , 3.2);
        TH1F *histNumerator_NVtx = new TH1F(("histNumerator_NVtx_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; Number of Primary Vertices ; Number of Events ",  30, -0.5 , 29.5);
        TH1F *histNumerator_Rho = new TH1F(("histNumerator_Rho_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; #rho (Energy Density) [GeV] ; Number of Events ",  90, 0 , 30);
        TH2F *histNumerator_PtEta = new TH2F(("histNumerator_PtEta_DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])).c_str() , "; p_{T} [GeV/c] ; #eta ; Number of Events ",  100, 0 , 100, 100, 0.0, 3.0);        
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
  TClonesArray *electronArr = new TClonesArray("mithep::TMuon");
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


  vector<string> inputfiles;
  if (inputFilename == "Full2011") {
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-smu-m10-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-smu-pr-v4_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-smu-a05-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-smu-o03-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11b-smu-pr-v1_fakerate_skim.root");
  } else if (inputFilename == "RUN2011A") {
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-smu-m10-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-smu-pr-v4_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-smu-a05-v1_fakerate_skim.root");
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11a-smu-o03-v1_fakerate_skim.root");
  }  else if (inputFilename == "RUN2011B") {
    inputfiles.push_back("/data/blue/vdutta/htt/2011/r11b-smu-pr-v1_fakerate_skim.root");
  } else {
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
    eventTree->SetBranchAddress("Muon", &electronArr); electronBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Muon", &muonArr);         muonBr = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Photon", &photonArr);     photonBr = eventTree->GetBranch("Photon");
    eventTree->SetBranchAddress("PV", &pvArr);             pvBr = eventTree->GetBranch("PV");
    eventTree->SetBranchAddress("PFJet", &jetArr);         jetBr = eventTree->GetBranch("PFJet");

    cout << "InputFile " << inputfiles[f] << " --- Total Events : " << eventTree->GetEntries() << endl;
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {       	
      infoBr->GetEntry(ientry);
		
      if (ientry % 100000 == 0) cout << "Event " << ientry << endl;

      //Double_t eventweight = info->eventweight;

      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  // not certified run? Skip to next event...
     

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
      pvArr->Clear();
      jetBr->GetEntry(ientry);
      pvBr->GetEntry(ientry);


      //********************************************************
      // Met
      //********************************************************
      Double_t met = info->pfMET;

      Int_t NMuons = muonArr->GetEntries();
      
      for(Int_t i=0; i<muonArr->GetEntries(); i++) {
        const mithep::TMuon *mu = (mithep::TMuon*)((*muonArr)[i]);
      
        Double_t leadingJetPt = -1;
        Double_t leptonJetPt = -1;
        //pass event selection     
        for(Int_t j=0; j<jetArr->GetEntries(); j++) {
          const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[j]);        
          if (jet->pt > leadingJetPt &&
              mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) > 1.0) {
            leadingJetPt = jet->pt;          
          }
          if (mithep::MathUtils::DeltaR(jet->phi, jet->eta, mu->phi, mu->eta) < 0.5) {
            leptonJetPt = jet->pt;
          }
        }
      
        //if there's no jet ( == isolated lepton?) then take pt of lepton
        if (leptonJetPt < 0) {
          leptonJetPt = mu->pt + mu->pfIso04; //mu->ChargedIso04 + mu->NeutralIso04_10Threshold;
        }


        //********************************************************
        // Photons
        //********************************************************
        Int_t NPhotons = 0;      
        Double_t photonEt = -1;
        for(Int_t p=0; p<photonArr->GetEntries(); p++) {
          const mithep::TPhoton *photon = (mithep::TPhoton*)((*photonArr)[p]);
        
          Bool_t isMu = kFALSE;
          for(Int_t e=0; e<electronArr->GetEntries(); e++) {
            const mithep::TMuon *tmpMu = (mithep::TMuon*)((*electronArr)[e]);   
            if (mithep::MathUtils::DeltaR(photon->phi, photon->eta, tmpMu->phi, tmpMu->eta) < 0.3) isMu = kTRUE;
          }
 


          //photon ID
          if ( photon->scEt > 20
               && !isMu
               && !photon->hasPixelSeed
               && photon->emIso04 < 2.0 + 0.006*photon->scEt
               && photon->hadIso04 < 2.0+0.0025*photon->scEt
               && photon->trkIso04 < 1.5 + 0.001*photon->scEt
               && ( (fabs(photon->eta) < 1.5 && photon->sigiEtaiEta < 0.01) || (fabs(photon->eta) >= 1.5 && photon->sigiEtaiEta < 0.028) )
            ) {
            continue;
          }


          mithep::FourVectorM phFourVector;
          mithep::FourVectorM muFourVector;
          muFourVector.SetCoordinates(mu->pt, mu->eta, mu->phi, 0.51099892e-3 );
          phFourVector.SetCoordinates(photon->scEt, photon->eta, photon->phi, 0.0 );
          mithep::FourVectorM dilepton = phFourVector+muFourVector;
        
        
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

              if (NMuons > 1) continue;

              if (sampleLabel[sampleTypeIndex] == "Mu8Sample" || sampleLabel[sampleTypeIndex] == "Mu15Sample" || sampleLabel[sampleTypeIndex] == "Mu8PtCombinedSample" || sampleLabel[sampleTypeIndex] == "Mu8Jet40Sample") {
                if (met > 20) continue;

                if (calcMt(info->pfMET, info->pfMETphi, mu) > 20) continue;

                if (!passZVeto(muonArr, denominatorType[denominatorTypeIndex])) continue;

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

              if (passMuonDenominatorCuts(info->triggerBits, mu, denominatorType[denominatorTypeIndex], sampleLabel[sampleTypeIndex])) {
              

                eventListFile << info->runNum << " " << info->lumiSec << " " << info->evtNum << " : " << mu->pt << " " << mu->eta << " " << mu->phi << " : " << passMuonNumeratorCuts(mu) << endl;
              

                DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt);
                DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(mu->eta));
                DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->phi);
                DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt, fabs(mu->eta));
                
                if (mu->pt < 20 && fabs(mu->eta) < 1.479) {
                  DenominatorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  DenominatorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                } 
                if (mu->pt < 20 && fabs(mu->eta) >= 1.479) {
                  DenominatorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  DenominatorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                } 
                if (mu->pt >= 20 && fabs(mu->eta) < 1.479) {
                  DenominatorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  DenominatorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                } 
                if (mu->pt >= 20 && fabs(mu->eta) >= 1.479) {
                  DenominatorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  DenominatorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                } 
              
                if (mu->pt < 30) {
                  Double_t iso = ( muonIsoPU(mu) ) / mu->pt;
                  LeptonJetPt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(TMath::Min(TMath::Max(leptonJetPt, 0.01),99.9));              
                  DenominatorIsolation[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(TMath::Min(TMath::Max(iso, 0.000001),0.9999999)); 
                }
              
                if (
                  (passMuonNumeratorCuts(mu))                     
                  ) {
                  NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt);
                  NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(fabs(mu->eta));
                  NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->phi);
                  NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                  NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(mu->pt, fabs(mu->eta) );   
     
                  if (mu->pt < 20 && fabs(mu->eta) < 1.479) {
                    NumeratorVector_Pt10To20_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                    NumeratorVector_Pt10To20_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  }
                  if (mu->pt < 20 && fabs(mu->eta) >= 1.479) {
                    NumeratorVector_Pt10To20_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                    NumeratorVector_Pt10To20_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  } 
                  if (mu->pt >= 20 && fabs(mu->eta) < 1.479) {
                    NumeratorVector_Pt20ToInf_Barrel_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                    NumeratorVector_Pt20ToInf_Barrel_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  } 
                  if (mu->pt >= 20 && fabs(mu->eta) >= 1.479) {
                    NumeratorVector_Pt20ToInf_Endcap_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(pvArr->GetEntriesFast());
                    NumeratorVector_Pt20ToInf_Endcap_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex]->Fill(info->rho);
                  }
                  
                }
              }

            } //loop over denominator types
          } //loop over sample types
        } //loop over ptThresholds

      } //loop over muons

    } //end loop over data     
  } //end loop over files
  eventListFile.close();

  delete info;
  delete electronArr;
  delete muonArr;
  delete jetArr;
  delete pvArr;
  
  
  //*****************************************************************************************
  //Make Efficiency Plots
  //*****************************************************************************************
  
  for (UInt_t denominatorTypeIndex = 0; denominatorTypeIndex < denominatorType.size(); ++denominatorTypeIndex) {
    for (UInt_t sampleTypeIndex = 0; sampleTypeIndex < sampleLabel.size(); ++sampleTypeIndex) {
      for (UInt_t ptThresholdIndex = 0; ptThresholdIndex < ptThreshold.size(); ++ptThresholdIndex) {

        
        Int_t ErrorType = 2; //Clopper Pearson errors
        TGraphAsymmErrors *efficiency_pt = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Pt[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Pt", ptbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_eta = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Eta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Eta", etabins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_phi = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Phi[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Phi", phibins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_nvtx = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_NVtx[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_NVtx", nvtxbins, ErrorType, -99, -99, 0, 1);
        TGraphAsymmErrors *efficiency_rho = EfficiencyUtils::createEfficiencyGraph(NumeratorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], DenominatorVector_Rho[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_Rho", rhobins, ErrorType, -99, -99, 0, 1);
        mithep::TH2DAsymErr *efficiency_PtEta = 
          EfficiencyUtils::createEfficiencyHist2D(NumeratorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], 
                                                          DenominatorVector_PtEta[denominatorTypeIndex][sampleTypeIndex][ptThresholdIndex], label+"DenominatorV"+IntToString(denominatorType[denominatorTypeIndex])+"_"+sampleLabel[sampleTypeIndex]+"_ptThreshold"+IntToString(ptThreshold[ptThresholdIndex])+"_PtEta", 
                                                          ptbins2D, etabins2D, ErrorType, kFALSE);
        

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
        file->WriteTObject(efficiency_rho, efficiency_rho->GetName(), "WriteDelete");
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
    
  gBenchmark->Show("WWTemplate");       
} 


Bool_t passMuonNumeratorCuts(const mithep::TMuon *mu) {

  if (mu->pt < 10) return kFALSE;
  if (fabs(mu->eta) > 2.4) return kFALSE;

  if(mu->nTkHits        < 11)    return kFALSE;
  if(mu->nPixHits       < 1)     return kFALSE;
  if(mu->muNchi2        > 10)    return kFALSE;
  if(mu->nMatch         < 2)     return kFALSE;
  if(mu->nValidHits     < 1)     return kFALSE;
  if(mu->ptErr/mu->pt > 0.1)   return kFALSE;
  if(fabs(mu->dz)       > 0.1)   return kFALSE;
  if(!(mu->typeBits & kGlobal))  return kFALSE;
  if(!(mu->typeBits & kTracker)) return kFALSE;
  
  if(fabs(mu->d0)>0.02)         return kFALSE;

  Double_t chargedIso = mu->pfIsoCharged;
  Double_t neutralIso = max(mu->pfIsoNeutral + mu->pfIsoGamma - 0.35 * mu->puIsoNoZ, 0.0);

  Double_t totalIso = chargedIso+neutralIso;

  if(fabs(mu->eta)<1.479) return (totalIso<0.15*(mu->pt));
  else                    return (totalIso<0.10*(mu->pt));

  return kTRUE;
}


Bool_t passMuonDenominatorCuts(Int_t triggerBits, const mithep::TMuon *mu, Int_t DenominatorType, string SampleType) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(mu->pt > 10 && fabs(mu->eta) < 2.4)) pass = kFALSE;

  //match to HLT
  if (SampleType == "Mu8Sample") {
    if (!( (triggerBits & kHLT_Mu8)
            && (mu->hltMatchBits & kHLT_Mu8_MuObj)
             
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Mu15Sample") {
    if (!( 
           (triggerBits & kHLT_Mu15)
            && (mu->hltMatchBits & kHLT_Mu15_MuObj)
                       
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Mu8PtCombinedSample") {
    if (!( 
          ((triggerBits & kHLT_Mu15)
            && (mu->hltMatchBits & kHLT_Mu15_MuObj)
            )
          ||
          ((triggerBits & kHLT_Mu8)
            && (mu->hltMatchBits & kHLT_Mu8_MuObj)
            )
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "Mu8Jet40Sample") {
    if (!( (triggerBits & kHLT_Mu8_Jet40)
           && (mu->hltMatchBits & kHLT_Mu8_MuObj)
          )
      ) {
      pass = kFALSE;
    }
  }
  if (SampleType == "PhotonJetsSample") {
    if (!(
          (triggerBits & kHLT_Mu8_Photon20_CaloIdVT_IsoT )
          && (mu->hltMatchBits & kHLT_Mu8_MuObj)
          )
      ) pass = kFALSE;
  }

  if (!passMuonDenominatorCuts(mu,DenominatorType)) pass = kFALSE;


  return pass;
}



Bool_t passMuonDenominatorCuts(const mithep::TMuon *mu, Int_t DenominatorType) {
  
  Bool_t pass = kTRUE;

  //eta , pt cut
  if (!(mu->pt > 10 && fabs(mu->eta) < 2.4)) pass = kFALSE;

  if (DenominatorType == 1) {

    if (! 
        ( (
            (Bool_t(mu->typeBits & kGlobal) 
             && mu->muNchi2 < 10.0
             && (mu->nValidHits > 0)
             && (mu->nMatch > 1 )
              )
            || 
            ( mu->typeBits & kTracker            
              && Bool_t(mu->qualityBits & kTMLastStationTight) 
              )
          ) 
          && mu->typeBits & kTracker
          && mu->nTkHits > 10                   
          && (mu->nPixHits > 0)
          && fabs(mu->d0) < 0.2
          && fabs(mu->dz) < 0.1
          && (mu->pfIso03) / mu->pt < 1.0
          && (mu->ptErr / mu->pt < 0.1)
          //&& (mu->TrkKink < 20)
          )
      ) pass = kFALSE;    
  }

  if (DenominatorType == 2) {
    if (! 
        ( //(
            //(Bool_t(mu->typeBits & kGlobal) 
	     mu->typeBits & kGlobal
	     && mu->muNchi2 < 10.0
             && (mu->nValidHits > 0)
             && (mu->nMatch > 1 )
              //)
            //|| 
            //( mu->typeBits & kTracker            
              //&& Bool_t(mu->qualityBits & kTMLastStationTight) 
              //)
          //)
          && mu->typeBits & kTracker
          && mu->nTkHits > 10
          && ( mu->nPixHits > 0)          
          && fabs(mu->d0) < 0.2
          && fabs(mu->dz) < 0.1
          && (mu->pfIso03) / mu->pt < 0.4
          && ( mu->ptErr / mu->pt < 0.1)
          //&& (mu->TrkKink < 20)
          )
      ) pass = kFALSE;    
  }

  if (DenominatorType == 3) {
    if (! 
        ( (
            (Bool_t(mu->typeBits & kGlobal) 
             && mu->muNchi2 < 10.0
             && (mu->nValidHits > 0)
             && (mu->nMatch > 1 )
              )
            || 
            ( mu->typeBits & kTracker            
              && Bool_t(mu->qualityBits & kTMLastStationTight) 
              )
          )
          && mu->typeBits & kTracker
          && mu->nTkHits > 10
          && ( mu->nPixHits > 0)
          && fabs(mu->d0) < 0.2
          && fabs(mu->dz) < 0.1
          && mu->trkIso03 / mu->pt < 0.3
          && mu->emIso03 / mu->pt < 0.3
          && mu->hadIso03 / mu->pt < 0.3
          && ( mu->ptErr / mu->pt < 0.1)
          //&& (mu->TrkKink < 20)
         )
      ) pass = kFALSE;    
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
