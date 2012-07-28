#ifndef HTTMVA_HH
#define HTTMVA_HH

#include "MitCommon/DataFormats/interface/Types.h"
#include <TError.h>
#include <TFile.h>
#include <TString.h>
#include <TMath.h>
#include <string>
#include <vector>
#include <utility>
#include <map>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace std;

using namespace TMVA;


class HttMVA {
public:
  HttMVA();
  ~HttMVA();

  enum MVAType {
    kNoSpinIP = 1,     // muon and electron d0, d0sig, ip3d, ip3dsig
    kNoSpinIPUB = 2,   // muon and electron d0, d0sig, ip3d, ip3dsig
    kNoSpin = 3,       // muon and electron d0sig, ip3dsig
    kWithSpin = 4,     // kNoSpin + polarization variables
    kMet = 5,          // topological variables
    kIPPlusMet = 6,    // comine ip + topological variables
    kIPMetNjets = 7,   // 6 + jet multiplicity
    kIPMetJets = 8,    // 7 + jet kinematics
    kVBF = 9,          // VBF
    kVBF2 = 10
  };

  void     Initialize(TString methodName, TString weightFile, HttMVA::MVAType type);

  Bool_t   IsInitialized() const { return fIsInitialized; }

  /*Double_t MVAValue(Double_t EleD0, Double_t EleD0Sig,
    Double_t EleIP3D, Double_t EleIP3DSig,
    Double_t MuD0, Double_t MuD0Sig,
    Double_t MuIP3D, Double_t MuIP3DSig);*/

  Double_t MVAValue(Double_t EleD0Sig , Double_t EleIP3DSig,
		    Double_t MuD0Sig, Double_t MuIP3DSig);

  //Double_t MVAValue(Double_t EleD0Sig , Double_t EleIP3DSig,
  //                  Double_t MuD0Sig, Double_t MuIP3DSig,
  //                  Double_t E1RF, Double_t E2RF);

  //Double_t MVAValue(Double_t ProjVis, Double_t ProjMet,
  //  		Double_t MtEle, Double_t MtMu,
  //		Double_t DilepMt, Double_t LepDPhi);

  Double_t MVAValue(Double_t Mjj, Double_t DEtajj,
		    Double_t DPhijj, Double_t Ptjj,
		    Double_t PtH, Double_t DPhiHjj,
		    Double_t DEtaVisJet, Double_t PtVis);
    
  Double_t MVAValue(Double_t EleD0, Double_t EleD0Sig,
		    Double_t EleIP3D, Double_t EleIP3DSig,
		    Double_t MuD0, Double_t MuD0Sig,
		    Double_t MuIP3D, Double_t MuIP3DSig,
		    Double_t ProjVis, Double_t ProjMet,
		    Double_t MtEle, Double_t MtMu,
		    Double_t DilepMt, Double_t LepDPhi);

  Double_t MVAValue(Double_t EleD0, Double_t EleD0Sig,
		    Double_t EleIP3D, Double_t EleIP3DSig,
		    Double_t MuD0, Double_t MuD0Sig,
		    Double_t MuIP3D, Double_t MuIP3DSig,
		    Double_t ProjVis, Double_t ProjMet,
		    Double_t MtEle, Double_t MtMu,
		    Double_t DilepMt, Double_t LepDPhi,
		    UInt_t NJets);

  Double_t MVAValue(Double_t EleD0, Double_t EleD0Sig,
		    Double_t EleIP3D, Double_t EleIP3DSig,
		    Double_t MuD0, Double_t MuD0Sig,
		    Double_t MuIP3D, Double_t MuIP3DSig,
		    Double_t ProjVis, Double_t ProjMet,
		    Double_t MtEle, Double_t MtMu,
		    Double_t DilepMt, Double_t LepDPhi,
		    UInt_t NJets, Double_t Jet1Pt,
		    Double_t Jet1Eta, Double_t Jet2Pt,
		    Double_t Jet2Eta);

protected:
  TMVA::Reader             *fTMVAReader;

  TString                   fMethodName;
  Bool_t                    fIsInitialized;

  Float_t                   fMVAVar_EleD0;
  Float_t                   fMVAVar_EleD0Sig;
  Float_t                   fMVAVar_EleIP3D;
  Float_t                   fMVAVar_EleIP3DSig;
  Float_t                   fMVAVar_MuD0;
  Float_t                   fMVAVar_MuD0Sig;
  Float_t                   fMVAVar_MuIP3D;
  Float_t                   fMVAVar_MuIP3DSig;
  Float_t                   fMVAVar_E1RF;
  Float_t                   fMVAVar_E2RF;
  Float_t                   fMVAVar_ProjVis;
  Float_t                   fMVAVar_ProjMet;
  Float_t                   fMVAVar_MtEle;
  Float_t                   fMVAVar_MtMu;
  Float_t                   fMVAVar_DilepMt;
  Float_t                   fMVAVar_LepDPhi;
  Float_t                   fMVAVar_NJets;
  Float_t                   fMVAVar_Mjj;
  Float_t                   fMVAVar_DEtajj;
  Float_t                   fMVAVar_DPhijj;
  Float_t                   fMVAVar_Ptjj;
  Float_t                   fMVAVar_PtH;
  Float_t                   fMVAVar_DPhiHjj;
  Float_t                   fMVAVar_DEtaVisJet;
  Float_t                   fMVAVar_PtVis;
  Float_t                   fMVAVar_Jet1Pt;
  Float_t                   fMVAVar_Jet1Eta;
  Float_t                   fMVAVar_Jet2Pt;
  Float_t                   fMVAVar_Jet2Eta;

};
 
#endif
