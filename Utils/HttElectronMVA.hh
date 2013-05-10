#ifndef HTTELECTRONMVA_HH
#define HTTELECTRONMVA_HH

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


class HttElectronMVA {
public:
  HttElectronMVA();
  ~HttElectronMVA();

  void     initialize(TString methodName, vector<TString> weightFiles);

  Bool_t   isInitialized() const { return fIsInitialized; }
  void     bindVariables();
  UInt_t   getMVABin(double eta,double pt ) const;

  Double_t mvaValue(Double_t fbrem, 
                      Double_t kfchi2,
                      Int_t    kfhits,
                      Double_t gsfchi2,
                      Double_t deta,
                      Double_t dphi,
                      Double_t detacalo,
                      Double_t see,
                      Double_t spp,
                      Double_t etawidth,
                      Double_t phiwidth,
                      Double_t e1x5e5x5,
                      Double_t R9,
                      Double_t HoE,
                      Double_t EoP,
                      Double_t IoEmIoP,
                      Double_t eleEoPout,
		      Double_t rho,
                      Double_t PreShowerOverRaw,
                      Double_t eta,
                      Double_t pt,
                      Bool_t printDebug = kFALSE );


protected:
  vector<TMVA::Reader*>      fTMVAReader;

  TString                    fMethodName;
  Bool_t                     fIsInitialized;
  UInt_t                     fNMVABins;

  Float_t                    fMVAVar_fbrem;
  Float_t                    fMVAVar_kfchi2;
  Float_t                    fMVAVar_kfhits;    //number of layers
  Float_t                    fMVAVar_gsfchi2;
  Float_t                    fMVAVar_deta;
  Float_t                    fMVAVar_dphi;
  Float_t                    fMVAVar_detacalo;
  Float_t                    fMVAVar_see;
  Float_t                    fMVAVar_spp;
  Float_t                    fMVAVar_etawidth;
  Float_t                    fMVAVar_phiwidth;
  Float_t                    fMVAVar_OneMinusE1x5E5x5;
  Float_t                    fMVAVar_R9;
  Float_t                    fMVAVar_HoE;
  Float_t                    fMVAVar_EoP;
  Float_t                    fMVAVar_IoEmIoP;
  Float_t                    fMVAVar_eleEoPout;
  Float_t                    fMVAVar_EoPout; 
  Float_t                    fMVAVar_PreShowerOverRaw;
  Float_t                    fMVAVar_eta;
  Float_t                    fMVAVar_pt;
  Float_t                    fMVAVar_rho;

};
 
#endif
