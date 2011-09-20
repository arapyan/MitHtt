//          mistag                         scale factor
// TCHEM  0.0175 \pm .0003 \pm .0038      1.21 \pm .02 \pm .17
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TH1.h>                    // histogram base class
#include <TNtuple.h>                // class to access ntuples
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TMath.h>                  // ROOT math functions
#include <TRegexp.h>                // ROOT regexp class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <map>                      // map tweedle dees to tweedle dums

#include "Common/MitStyleRemix.hh"  // style settings for drawing
#include "Common/CSample.hh"        // helper class for organizing input ntuple files
#include "Common/MyTools.hh"        // miscellaneous helper functions
#include "Common/CPlot.hh"          // helper class for plots

// define structures to read in ntuple
#include "MitHtt/NtupleDefs/interface/HiggsAnaDefs.hh"
#include "MitHtt/NtupleDefs/interface/TEventInfo.hh"
#include "MitHtt/NtupleDefs/interface/TGenInfo.hh"
#include "MitHtt/NtupleDefs/interface/TMuon.hh" 
#include "MitHtt/NtupleDefs/interface/TElectron.hh"
#include "MitHtt/NtupleDefs/interface/TJet.hh"   
#include "MitHtt/NtupleDefs/interface/TVertex.hh"   
#include "MitHtt/NtupleDefs/interface/TSVFit.hh"

#include "MitHtt/Utils/RecoilCorrector.hh"
#include "MitHtt/Emu/EScale/EScale.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh"

// define structure for output ntuple
#include "EmuData.hh"


enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };  // dataset type  
enum { kMuMu, kEleEle, kEleMu, kMuEle };      // final state type
enum { kNo, kDown, kUp };                     // jet energy uncertainties

const Double_t pi = 3.14159265358979;
  
class Selector
{
 public:
  Selector(TString conf, TString ouputDir, Double_t lumi);
  ~Selector();
  void     SampleLoop(); // Main processing loop over CSamples, files, and events

  UInt_t      fJetUnc;        // jet uncertainty configuration: vary energy up/down
  UInt_t      fBtagEff;       // vary btag efficiency up/down?
  UInt_t      fMistag;        // vary mistag rate up/down?
  Double_t    fMuonPt1Min;    // lead muon min pt
  Double_t    fEMuMuonPt2Min; // second muon min pt for emu selection
  Double_t    fMuMuMuonPt2Min;// second muon min pt for mumu selection
  Double_t    fElePt1Min;     // lead ele min pt  
  Double_t    fElePt2Min;     // second ele min pt
  Double_t    fJetPtMin;      // min jet pt for all jets
  Double_t    fBJetPtMin;     // min b-jet pt

  TH1F       *fFailHist;
  ULong64_t   fEvtFail;
    
 protected:
  void     ParseConfig		(TString conf);
  void     InitOutput		();
  void     InitNPU		(TString sfname,TH1D *&hpu,TH1D *&hpuRwgt,vector<Double_t> &puwgtv);
  Double_t GetWeight		(TString sfname, CSample *samp);
  Double_t UnskimmedEntries	(TString skimname);
  void     FillNPU		(TH1D *hpu, TH1D *hpuRwgt,vector<Double_t> puwgtv);
  void     WriteNPU             (TH1D *hpu, TH1D *hpuRwgt, TString basename);
  Bool_t   EvtFail              (ULong64_t failshift);
  void     MuonLoop		(Bool_t muMu, vector<const mithep::TMuon*> &goodMuonsv, vector<const mithep::TMuon*> &looseMuonsv);
  void     ElectronLoop		(vector<mithep::TElectron*> &goodElectronsv, vector<const mithep::TMuon*> looseMuonsv,
				 vector<const mithep::TMuon*> goodMuonsv);
  Bool_t   PassEmu              (vector<const mithep::TMuon*> goodMuonsv, vector<mithep::TElectron*> goodElectronsv,
				 TLorentzVector &lep1, TLorentzVector &lep2, TLorentzVector &dilep, Int_t &finalState,
				 const mithep::TMuon *&mu, mithep::TElectron *&ele);
  Bool_t   PassMuMu             (vector<const mithep::TMuon*> goodMuonsv, vector<mithep::TElectron*> goodElectronsv,
				 TLorentzVector &lep1, TLorentzVector &lep2, TLorentzVector &dilep, Int_t &finalState,
				 const mithep::TMuon *&mu, const mithep::TMuon *&mu_2);
  void     JetLoop		(TLorentzVector lep1, TLorentzVector lep2, UInt_t &njets, UInt_t &nbjets,
				 const mithep::TJet *&jet1, const mithep::TJet *&jet2, const mithep::TJet *&bjet);
  void     Projections		(TLorentzVector lep1, TLorentzVector lep2, TLorentzVector dilep,
				 Bool_t doRecoil, Double_t &projVis, Double_t &projMet, Double_t &met, Double_t &metphi);
  Double_t GetTrigEff		(const mithep::TMuon *mu, const mithep::TElectron *ele);
  
  Bool_t   IsBtagged		(mithep::TJet *jet);
  vector<Double_t> generate_flat10_weights(TString datafname, TString mcfname);
  Double_t eleTrigEff		(const mithep::TElectron *ele);
  Double_t muTrigEff		(const mithep::TMuon *mu);
  // k-factors not fully implemented...
  TH1D*    kfInit		(TString kfilename, Int_t mH);
  Double_t kfValue		(const Double_t pt, const TH1D* hKF);
  Int_t    higgsmass		(TString basename);
  void     RemoveSkimName       (TString &name);

  // per-CSample variables
  TTree            *fEventTree;		// pointer to input tree
  TFile            *fOutfile;		// output file
  TTree            *fOuttree;		// output ntuple
  UInt_t            fIsam;       	// current sample index

  // per-input-file variables
  UInt_t            fIfile;		// current file index
  Bool_t            fIsdata;		// is this file data?
  TString           fDataNPVfname;	// file with estimated data npu distribution
  TString           fmcNPVfname;	// file with weight histogram for npu reweighting for this mc sample

  // for corrections
  TRandom           fRandm;             // random generator for IsBtagged
  RecoilCorrector   fCorrector;         // class to apply recoil corrections

  // configuration values  
  vector<TString>   fSnamev;            // sample names (for output file)  
  vector<CSample*>  fSamplev;           // data/MC samples
  TString           fNtupDir;
  const Double_t    fLumi;
  Bool_t	    fDoNpuRwgt;         // reweight mc events to the "observed" npu distribution
  Bool_t	    fCheckNpuHists;     // write out png's to see how the npu distrib looks
  Bool_t	    fMakeNpuHists;      // make hists for future reweights
  Bool_t            fDebug;

  // structures to store info from input tree
  mithep::TEventInfo *fInfo;
  mithep::TGenInfo   *fGen;
  TClonesArray       *fMuonArr;
  TClonesArray       *fElectronArr;
  TClonesArray       *fJetArr;
  TClonesArray       *fPvArr;
  TClonesArray       *fSvfitArr;

  // structures for output ntuple
  EmuData             fData;
  Double_t            fRawMet,fRawprojvar,fNpuWgt;
  UInt_t              fNpt20jets;
  const UInt_t        kMaxPt20Jets;     // maximum number of b-jets in the array
  TArrayF             fBtagArray;       // array to hold b-tag discr. values for pt-20 jets

};

//----------------------------------------------------------------------------------------
// return the -(z-boost) of the boson, as approximated from the theta coordinates
// of the two leptons
Double_t v(Double_t t1, Double_t t2)
{
  return (-cos(t1) - cos(t2)) / (1+cos(t1+t2));
};
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
};
 

