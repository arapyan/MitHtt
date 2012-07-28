#include "HttMVA.hh"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
  

using namespace mithep;
using namespace std;

//--------------------------------------------------------------------------------------------------
HttMVA::HttMVA() :
fMethodName("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
  fTMVAReader = 0;
}

//--------------------------------------------------------------------------------------------------
HttMVA::~HttMVA()
{
  if (fTMVAReader) delete fTMVAReader;
}

//--------------------------------------------------------------------------------------------------
void HttMVA::Initialize(TString methodName, TString weightFile, HttMVA::MVAType type) {

  fIsInitialized = kTRUE;

  fMethodName = methodName;

  if (fTMVAReader) delete fTMVAReader;

  fTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );
  fTMVAReader->SetVerbose(kTRUE);

  if(type==kIPPlusMet || type==kIPMetNjets || type==kIPMetJets) {
    fTMVAReader->AddVariable( "projvis",         &fMVAVar_ProjVis   );
    fTMVAReader->AddVariable( "projmet",         &fMVAVar_ProjMet   );
    fTMVAReader->AddVariable( "mtele",           &fMVAVar_MtEle     );
    fTMVAReader->AddVariable( "mtmu",            &fMVAVar_MtMu      );
    fTMVAReader->AddVariable( "eled0",           &fMVAVar_EleD0     );
    fTMVAReader->AddVariable( "eled0sig",        &fMVAVar_EleD0Sig  );
    fTMVAReader->AddVariable( "eleip3d",         &fMVAVar_EleIP3D   );
    fTMVAReader->AddVariable( "eleip3dsig",      &fMVAVar_EleIP3DSig);
    fTMVAReader->AddVariable( "mtll",            &fMVAVar_DilepMt   );
    fTMVAReader->AddVariable( "dphill",          &fMVAVar_LepDPhi   );
    fTMVAReader->AddVariable( "mud0",            &fMVAVar_MuD0      );
    fTMVAReader->AddVariable( "mud0sig",         &fMVAVar_MuD0Sig   );
    fTMVAReader->AddVariable( "muip3d",          &fMVAVar_MuIP3D    );
    fTMVAReader->AddVariable( "muip3dsig",       &fMVAVar_EleIP3DSig);
    if(type==kIPMetNjets || type==kIPMetJets) fTMVAReader->AddVariable( "njets",       &fMVAVar_NJets);
    if(type==kIPMetJets) {
      fTMVAReader->AddVariable( "ptj1",          &fMVAVar_Jet1Pt    );
      fTMVAReader->AddVariable( "etaj1",         &fMVAVar_Jet1Eta   );
      fTMVAReader->AddVariable( "ptj2",          &fMVAVar_Jet2Pt    );
      fTMVAReader->AddVariable( "etaj2",         &fMVAVar_Jet2Eta   );
    }
  } else if(type==kVBF) {
    fTMVAReader->AddVariable( "mjj",             &fMVAVar_Mjj       );
    fTMVAReader->AddVariable( "detajj",          &fMVAVar_DEtajj    );
    fTMVAReader->AddVariable( "dphijj",          &fMVAVar_DPhijj    );
    fTMVAReader->AddVariable( "ptjj",            &fMVAVar_Ptjj      );
    fTMVAReader->AddVariable( "hpt",             &fMVAVar_PtH       );
    fTMVAReader->AddVariable( "dphihjj",         &fMVAVar_DPhiHjj   );
    fTMVAReader->AddVariable( "detavisj",        &fMVAVar_DEtaVisJet);
    fTMVAReader->AddVariable( "ptvis",           &fMVAVar_PtVis     );
  } else if(type==kVBF2) {
    fTMVAReader->AddVariable( "mjj",             &fMVAVar_Mjj       );
    fTMVAReader->AddVariable( "dEta",            &fMVAVar_DEtajj    );
    fTMVAReader->AddVariable( "dPhi",            &fMVAVar_DPhijj    );
    fTMVAReader->AddVariable( "ditau_pt",        &fMVAVar_PtH       );
    fTMVAReader->AddVariable( "dijet_pt",        &fMVAVar_Ptjj      );
    fTMVAReader->AddVariable( "dPhi_hj",         &fMVAVar_DPhiHjj   );
    fTMVAReader->AddVariable( "C1",              &fMVAVar_DEtaVisJet);
    fTMVAReader->AddVariable( "C2",              &fMVAVar_PtVis     );
  } else if(type==kMet) {
    fTMVAReader->AddVariable( "projvis",         &fMVAVar_ProjVis   );
    fTMVAReader->AddVariable( "projmet",         &fMVAVar_ProjMet   );
    fTMVAReader->AddVariable( "mtele",           &fMVAVar_MtEle     );
    fTMVAReader->AddVariable( "mtmu",            &fMVAVar_MtMu      );
    fTMVAReader->AddVariable( "mtll",            &fMVAVar_DilepMt   );
    fTMVAReader->AddVariable( "dphill",          &fMVAVar_LepDPhi   );
  } else if(type==kNoSpinIPUB) {
    fTMVAReader->AddVariable( "eleubd0",         &fMVAVar_EleD0   );
    fTMVAReader->AddVariable( "eleubd0sig",         &fMVAVar_EleD0Sig   );
    fTMVAReader->AddVariable( "eleubip3d",       &fMVAVar_EleIP3D );
    fTMVAReader->AddVariable( "eleubip3dsig",       &fMVAVar_EleIP3DSig );
    fTMVAReader->AddVariable( "muubd0",          &fMVAVar_MuD0    );
    fTMVAReader->AddVariable( "muubd0sig",          &fMVAVar_MuD0Sig    );
    fTMVAReader->AddVariable( "muubip3d",        &fMVAVar_MuIP3D  );
    fTMVAReader->AddVariable( "muubip3dsig",        &fMVAVar_MuIP3DSig  );
  } else {
    if(type==kNoSpinIP) fTMVAReader->AddVariable( "eled0",         &fMVAVar_EleD0   );
    fTMVAReader->AddVariable( "eled0sig",         &fMVAVar_EleD0Sig   );
    if(type==kNoSpinIP) fTMVAReader->AddVariable( "eleip3d",       &fMVAVar_EleIP3D );
    fTMVAReader->AddVariable( "eleip3dsig",       &fMVAVar_EleIP3DSig );
    if(type==kNoSpinIP) fTMVAReader->AddVariable( "mud0",          &fMVAVar_MuD0    );
    fTMVAReader->AddVariable( "mud0sig",          &fMVAVar_MuD0Sig    );
    if(type==kNoSpinIP) fTMVAReader->AddVariable( "muip3d",        &fMVAVar_MuIP3D  );
    fTMVAReader->AddVariable( "muip3dsig",        &fMVAVar_MuIP3DSig  );
    if(type==kWithSpin) {
      fTMVAReader->AddVariable( "e1rf",             &fMVAVar_E1RF       );
      fTMVAReader->AddVariable( "e2rf",             &fMVAVar_E2RF       );
    }
  }


  fTMVAReader->BookMVA(fMethodName , weightFile);


  cout << "HttMVA Initialization\n";
  cout << "MethodName : " << fMethodName << endl;
  cout << "Load weights file : " << weightFile << endl;

}

//--------------------------------------------------------------------------------------------------
/*Double_t HttMVA::MVAValue(Double_t EleD0, Double_t EleD0Sig,
                        Double_t EleIP3D, Double_t EleIP3DSig,
                        Double_t MuD0, Double_t MuD0Sig,
                        Double_t MuIP3D, Double_t MuIP3DSig
  ) {

  if (!fIsInitialized) {
    std::cout << "Error: HttMVA not properly initialized.\n";
    return -9999;
  }
  
  //set all input variables
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleD0Sig = EleD0Sig;
  fMVAVar_EleIP3D = EleIP3D;
  fMVAVar_EleIP3DSig = EleIP3DSig;
  fMVAVar_MuD0 = MuD0;
  fMVAVar_MuD0Sig = MuD0Sig;
  fMVAVar_MuIP3D = MuIP3D;
  fMVAVar_MuIP3DSig = MuIP3DSig;
  
  Double_t mva = -9999;
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;
  
  mva = reader->EvaluateMVA( fMethodName );
  
  return mva;
}*/

//--------------------------------------------------------------------------------------------------
Double_t HttMVA::MVAValue(Double_t EleD0Sig , Double_t EleIP3DSig,
                        Double_t MuD0Sig, Double_t MuIP3DSig
  ) {

  if (!fIsInitialized) {
    std::cout << "Error: HttMVA not properly initialized.\n";
    return -9999;
  }

  //set all input variables
  fMVAVar_EleD0Sig = EleD0Sig;
  fMVAVar_EleIP3DSig = EleIP3DSig;
  fMVAVar_MuD0Sig = MuD0Sig;
  fMVAVar_MuIP3DSig = MuIP3DSig;

  Double_t mva = -9999;
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;

  mva = reader->EvaluateMVA( fMethodName );

  return mva;
}

/*//--------------------------------------------------------------------------------------------------
Double_t HttMVA::MVAValue(Double_t EleD0Sig , Double_t EleIP3DSig,
                        Double_t MuD0Sig, Double_t MuIP3DSig,
                        Double_t E1RF, Double_t E2RF
  ) {

  if (!fIsInitialized) {
    std::cout << "Error: HttMVA not properly initialized.\n";
    return -9999;
  }

  //set all input variables
  fMVAVar_EleD0Sig = EleD0Sig;
  fMVAVar_EleIP3DSig = EleIP3DSig;
  fMVAVar_MuD0Sig = MuD0Sig;
  fMVAVar_MuIP3DSig = MuIP3DSig;
  fMVAVar_E1RF = E1RF;
  fMVAVar_E2RF = E2RF;

  Double_t mva = -9999;
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;

  mva = reader->EvaluateMVA( fMethodName );

  return mva;
}*/

/*//--------------------------------------------------------------------------------------------------
Double_t HttMVA::MVAValue(Double_t ProjVis, Double_t ProjMet,
                        Double_t MtEle, Double_t MtMu,
                        Double_t DilepMt, Double_t LepDPhi
  ) {
    
  if (!fIsInitialized) {
    std::cout << "Error: HttMVA not properly initialized.\n";
    return -9999;
  }
  
  //set all input variables
  fMVAVar_ProjVis = ProjVis;
  fMVAVar_ProjMet = ProjMet;
  fMVAVar_MtEle = MtEle;
  fMVAVar_MtMu = MtMu;
  fMVAVar_DilepMt = DilepMt;
  fMVAVar_LepDPhi = LepDPhi;
  
  Double_t mva = -9999;
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;

  mva = reader->EvaluateMVA( fMethodName );

  return mva;
}*/

//--------------------------------------------------------------------------------------------------
Double_t HttMVA::MVAValue(Double_t Mjj, Double_t DEtajj,
                        Double_t DPhijj, Double_t Ptjj,
                        Double_t PtH, Double_t DPhiHjj,
			Double_t DEtaVisJet, Double_t PtVis
  ) {

  if (!fIsInitialized) {
    std::cout << "Error: HttMVA not properly initialized.\n";
    return -9999;
  }

  //set all input variables
  fMVAVar_Mjj = Mjj;
  fMVAVar_DEtajj = DEtajj;
  fMVAVar_DPhijj = DPhijj;
  fMVAVar_Ptjj = Ptjj;
  fMVAVar_PtH = PtH;
  fMVAVar_DPhiHjj = DPhiHjj;
  fMVAVar_DEtaVisJet = DEtaVisJet;
  fMVAVar_PtVis = PtVis;
  

  Double_t mva = -9999;
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;

  mva = reader->EvaluateMVA( fMethodName );

  return mva;
}
//--------------------------------------------------------------------------------------------------
Double_t HttMVA::MVAValue(Double_t EleD0, Double_t EleD0Sig,
                        Double_t EleIP3D, Double_t EleIP3DSig,
                        Double_t MuD0, Double_t MuD0Sig,
                        Double_t MuIP3D, Double_t MuIP3DSig,
			Double_t ProjVis, Double_t ProjMet,
                        Double_t MtEle, Double_t MtMu,
                        Double_t DilepMt, Double_t LepDPhi
  ) {

  if (!fIsInitialized) {
    std::cout << "Error: HttMVA not properly initialized.\n";
    return -9999;
  }

  //set all input variables
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleD0Sig = EleD0Sig;
  fMVAVar_EleIP3D = EleIP3D;
  fMVAVar_EleIP3DSig = EleIP3DSig;
  fMVAVar_MuD0 = MuD0;
  fMVAVar_MuD0Sig = MuD0Sig;
  fMVAVar_MuIP3D = MuIP3D;
  fMVAVar_MuIP3DSig = MuIP3DSig;
  fMVAVar_ProjVis = ProjVis;
  fMVAVar_ProjMet = ProjMet;
  fMVAVar_MtEle = MtEle;
  fMVAVar_MtMu = MtMu;
  fMVAVar_DilepMt = DilepMt;
  fMVAVar_LepDPhi = LepDPhi;

  Double_t mva = -9999;
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;

  mva = reader->EvaluateMVA( fMethodName );

  return mva;
}

//--------------------------------------------------------------------------------------------------
Double_t HttMVA::MVAValue(Double_t EleD0, Double_t EleD0Sig,
                        Double_t EleIP3D, Double_t EleIP3DSig,
                        Double_t MuD0, Double_t MuD0Sig,
                        Double_t MuIP3D, Double_t MuIP3DSig,
                        Double_t ProjVis, Double_t ProjMet,
                        Double_t MtEle, Double_t MtMu,
                        Double_t DilepMt, Double_t LepDPhi,
			UInt_t NJets
  ) {

  if (!fIsInitialized) {
    std::cout << "Error: HttMVA not properly initialized.\n";
    return -9999;
  }

  //set all input variables
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleD0Sig = EleD0Sig;
  fMVAVar_EleIP3D = EleIP3D;
  fMVAVar_EleIP3DSig = EleIP3DSig;
  fMVAVar_MuD0 = MuD0;
  fMVAVar_MuD0Sig = MuD0Sig;
  fMVAVar_MuIP3D = MuIP3D;
  fMVAVar_MuIP3DSig = MuIP3DSig;
  fMVAVar_ProjVis = ProjVis;
  fMVAVar_ProjMet = ProjMet;
  fMVAVar_MtEle = MtEle;
  fMVAVar_MtMu = MtMu;
  fMVAVar_DilepMt = DilepMt;
  fMVAVar_LepDPhi = LepDPhi;
  fMVAVar_NJets = NJets;

  Double_t mva = -9999;
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;

  mva = reader->EvaluateMVA( fMethodName );

  return mva;
}

//--------------------------------------------------------------------------------------------------
Double_t HttMVA::MVAValue(Double_t EleD0, Double_t EleD0Sig,
                        Double_t EleIP3D, Double_t EleIP3DSig,
                        Double_t MuD0, Double_t MuD0Sig,
                        Double_t MuIP3D, Double_t MuIP3DSig,
                        Double_t ProjVis, Double_t ProjMet,
                        Double_t MtEle, Double_t MtMu,
                        Double_t DilepMt, Double_t LepDPhi,
                        UInt_t NJets, Double_t Jet1Pt,
			Double_t Jet1Eta, Double_t Jet2Pt,
			Double_t Jet2Eta
  ) {

  if (!fIsInitialized) {
    std::cout << "Error: HttMVA not properly initialized.\n";
    return -9999;
  }

  //set all input variables
  fMVAVar_EleD0 = EleD0;
  fMVAVar_EleD0Sig = EleD0Sig;
  fMVAVar_EleIP3D = EleIP3D;
  fMVAVar_EleIP3DSig = EleIP3DSig;
  fMVAVar_MuD0 = MuD0;
  fMVAVar_MuD0Sig = MuD0Sig;
  fMVAVar_MuIP3D = MuIP3D;
  fMVAVar_MuIP3DSig = MuIP3DSig;
  fMVAVar_ProjVis = ProjVis;
  fMVAVar_ProjMet = ProjMet;
  fMVAVar_MtEle = MtEle;
  fMVAVar_MtMu = MtMu;
  fMVAVar_DilepMt = DilepMt;
  fMVAVar_LepDPhi = LepDPhi;
  fMVAVar_NJets = NJets;
  fMVAVar_Jet1Pt = Jet1Pt;
  fMVAVar_Jet1Eta = Jet1Eta;
  fMVAVar_Jet2Pt = Jet2Pt;
  fMVAVar_Jet2Eta = Jet2Eta;

  Double_t mva = -9999;
  TMVA::Reader *reader = 0;
  reader = fTMVAReader;

  mva = reader->EvaluateMVA( fMethodName );

  return mva;
}
