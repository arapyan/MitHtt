#ifndef OUTPUT_HH
#define OUTPUT_HH

#include <TTree.h> 
#include <TBranch.h>
#include <TString.h>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access trees
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TMath.h>                  // ROOT math functions

#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TPFTau.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"   
#include "MitHtt/Ntupler/interface/TVertex.hh"   
#include "MitHtt/Ntupler/interface/TSVfit.h"

#include "MitHtt/Utils/RecoilCorrector.hh"
#include "MitHtt/Utils/HttMVA.hh" 

class Output {
public:

  Output(TString name);
  
  int fRun;      //Bookeeping
  int fLumi;     //Bookeeping
  int fEvt;      //Bookeeping

  //Event variables
  float fpvx;   // x-primary vertex
  float fpvy;   // y-primary vertex
  float fpvz;   // z-primary vertex

  int fNPV;      // NPV
  int fNPU;      // NPU
  float fRho;    // Rho

  //Event Weight
  float fMCWeight;  //MC Weight(xs/nevents*additional weight*kfactor(gg higgs)
  float fPUWeight;  //Pileup Weight
  float fEffWeight; //Efficiency Scale factor
  float fWeight;    //mcweight*puweight*effweight

  //Mass variables
  float fMass;      //gen mass if applicible 
  float fgenpt;     //gen pt of the boson
  float fgenphi;    //gen phi of the boson 
  float fMVis;      //visible mass
  
  //First Lepton: muon for muTau, electron for eTau, muon for emu, leading pT Tau for hadronic
  float fPt1;      //pt1
  float fPhi1;     //phi1
  float fEta1;     //Eta
  float fM1;       //Mass
  int   fq1;       //charge
  float fIso1;     //Isolation variable
  float fD01;      //d0 with respect to primary vertex
  float fDZ1;      //dZ with respect to primary vertex
  bool  fPassIso1; // passes default iso?
  float fMt1;      // mT of lepton wrt to pf met
  float fMVAMt1;   // mT of first leptron wrt to MVA met
  int   fngamma1;   //number of gamma candidates 
  int   fnprong1;  // number of prongs
  bool  fantiele1; //for tau: passes the antielectron discriminator mva

  ///Second lepton :  hadronic Tau for muTau had for eTau, electron for emu, Trailing (in pT)  Tau for TauTau
  float fPt2;      //pt2
  float fPhi2;     //phi2
  float fEta2;     //Eta
  float fM2;       //Mass
  int   fq2;       //charge
  float fIso2;     //Isolation variable
  float fD02;      //d0 with respect to primary vertex
  float fDZ2;      // dZ with respect to primary vertex
  bool  fPassIso2; // passes default iso?
  float fMt2;      // mT of lepton wrt to pf met
  float fMVAMt2;   // mT of first leptron wrt to MVA met
  int   fngamma2;    // number of gamma candidates
  int   fnprong2;    // number of prongs
  bool  fantiele2; //for tau: passes the antielectron discriminator mva

  float fdrll;     //dR between two leptons

  //Met related variables
  float fMet;          //met
  float fMetPhi;       //met phi
  float fMVAMet;       //mvamet
  float fMVAMetPhi;    //mvametphi
  float fPZetaVis;     //pZeta visible
  float fPZetaMiss;    //pZeta Missing
  float fPZetaMVAMiss; //pZdeta MVA Missing
  
  //MET covariance matrix
  float fMetCov00;
  float fMetCov01;
  float fMetCov10;
  float fMetCov11;
  //MVAMet covariance mmatrix
  float fMVACov00;
  float fMVACov01;
  float fMVACov10;
  float fMVACov11;
  
  //First Jet: leading jet after applying jet energy corrections
  float fJPt1;       //Jet Pt after corrections
  float fJEta1;      //Jet Eta
  float fJPhi1;      //Jet Phi
  float fJM1;        //Jet Mass
  float fJPtUnc1;    //Jet Unc (relative to Jet corrected pT)
  float fJMVA1;      //JetMVA id values
  float fJcsv1;      //CSV discriminator value 
  bool  fJPass1;     //Jet passes PU ID loose WP?

    
  //Second Jet: 2nd leading jet after applying jet energy corrections
  float fJPt2;      //Jet Pt after corrections
  float fJEta2;      //Jet Eta
  float fJPhi2;      //Jet Phi
  float fJM2;        //Jet Mass
  float fJPtUnc2;    //Jet Unc (relative to Jet corrected pT)
  float fJMVA2;      //JetMVA id values
  float fJcsv2;      //CSV discriminator value 
  bool  fJPass2;     //Jet passes PU ID loose WP?
  
  //Leading B Tagged Jet : leading btagged jet passing btag wp (pt > 20+cvs medium)
  float fBTagPt1;     //Corrected BTag Pt
  float fBTagEta1;    //Btag Eta
  float fBTagPhi1;    //Btag Phi
  float fBTagM1;      //Btag Mass
  float fbcsv1;       //B CSV

  //Second B Tagged Jet : 2nd leading btagged jet passing btag wp (pt > 20+cvs medium)
  float fBTagPt2;     //Corrected BTag Pt
  float fBTagEta2;    //Btag Eta
  float fBTagPhi2;    //Btag Phi
  float fBTagM2;      //Btag Mass
  float fbcsv2;       //B CSV
  
  //Di Jet kinematic variables (usefull for VBF selection)
  float fMJJ;        //Mass Di Jet system
  float fJDEta;      //|jet1-jet2|
  int fNJetInGap;    // # of Jets between the two leading jets
  float fMVA;        // VBF MVA value
  
  //Variables usefull for VBF MVA
  float fJDPhi;  //Delta Phi between two leading jets
  float fDiJetPt; //Pt of the di jet system
  float fDiJetPhi; //Phi of the di jet system
  float fHDJetPhi; //Phi of the di jet system
  float fVisJetEta; //TMath::Min(eta_vis - jeta,eta_vis,jeta2)
  float fPtVis;     //Pt Vis
  float fPtH;       //Pt of the system
  float fPtHMVA;       //Pt of the system mva
  
  //number of btags passing btag id (pt > 20)
  int fNBTag;   
  //number of jets passing jet id (pt > 30)
  int fNJets;

  //generator lepton pt, eta, phi for embedded
  float fGenPt1;  //pT leading
  float fGenPhi1; //Phi leading
  float fGenEta1; //Eta leading
  int   fGenId1;  //Pdg Id leading
  float fGenPt2;  //pT sub-leading
  float fGenPhi2; //Phi sub-leading
  float fGenEta2; //Eta sub-leading
  int   fGenId2;  //Pdg Id leading
  int   fGenMatch; //Matched to ll/l+Tau/Tau+Tau/l+Jet/Tau+Jet/None of above (1/2/3/4/5/0)
 
  int  doRecoil; //Recoil corrections  
  bool doEmu;    //use Emu version
  UInt_t npt20jets;
  TArrayF btagArray;
  TArrayF jptArray;
  TArrayF jetaArray;

  void fillMuon(const mithep::TMuon *muon, bool location, double iso, bool passiso);
  void fillElectron(const mithep::TElectron *ele, bool location, double iso, bool passiso, unsigned int scale=0);
  void fillTau(const mithep::TPFTau *tau, bool first,  bool passiso);
  void fillCov(mithep::TSVfit *svfit); 	
  void fillGen(mithep::TGenInfo *gen);
  void fillEvent(mithep::TEventInfo *event, HttMVA *vbfmva, int npv, double scalecorr=0);  //always to be called the last
  void fillJets(const mithep::TJet *jet1,const mithep::TJet *jet2,const mithep::TJet *bjet1, const mithep::TJet *bjet2, int njets, int bjets, int npt20, int nCentralJets);
  void setupRecoil(int doRec, bool is2012=true, bool isEmu=false);
  void save();
  void cd();
protected:
  /// output file
  TFile* fOutputFile;           
  /// output tree
  TTree* fEventTree;
  void setupOutput(TString name);
   // recoil corrections
  RecoilCorrector *corrector;	
};
#endif
