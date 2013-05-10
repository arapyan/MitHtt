#include "HttElectronMVA.hh"
#include <TFile.h>
#include <TRandom3.h>
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
  

using namespace std;

//--------------------------------------------------------------------------------------------------
HttElectronMVA::HttElectronMVA() :
fMethodName("BDTG method"),
fIsInitialized(kFALSE)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
HttElectronMVA::~HttElectronMVA()
{
  for(UInt_t i=0; i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void HttElectronMVA::initialize(TString methodName, vector<TString> weightFiles) {

  //clean up first
  for (uint i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
  fTMVAReader.clear();

  //initialize
  fIsInitialized = kTRUE;
  fMethodName = methodName;

  //Define expected number of bins
  fNMVABins = 6;

  //Check number of weight files given
  if (fNMVABins != weightFiles.size() ) {
    cout << "Error: Expected Number of bins = " << fNMVABins << " does not equal to weightFiles.size() = " 
              << weightFiles.size() << endl;
    assert(fNMVABins == weightFiles.size());
  }

  for(UInt_t i=0; i<fNMVABins; ++i) {
    TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    tmpTMVAReader->SetVerbose(kTRUE);

     
    // Pure tracking variables
    tmpTMVAReader->AddVariable("fbrem",           &fMVAVar_fbrem);
    tmpTMVAReader->AddVariable("kfchi2",          &fMVAVar_kfchi2);
    tmpTMVAReader->AddVariable("kfhits",          &fMVAVar_kfhits);
    tmpTMVAReader->AddVariable("gsfchi2",         &fMVAVar_gsfchi2);

    // Geometrical matchings
    tmpTMVAReader->AddVariable("deta",            &fMVAVar_deta);
    tmpTMVAReader->AddVariable("dphi",            &fMVAVar_dphi);
    tmpTMVAReader->AddVariable("detacalo",        &fMVAVar_detacalo);
         
    // Pure ECAL -> shower shapes
    tmpTMVAReader->AddVariable("see",             &fMVAVar_see);
    tmpTMVAReader->AddVariable("spp",             &fMVAVar_spp);
    tmpTMVAReader->AddVariable("etawidth",        &fMVAVar_etawidth);
    tmpTMVAReader->AddVariable("phiwidth",        &fMVAVar_phiwidth);
    tmpTMVAReader->AddVariable("e1x5e5x5",        &fMVAVar_OneMinusE1x5E5x5);
    tmpTMVAReader->AddVariable("R9",              &fMVAVar_R9);
     
    // Energy matching
    tmpTMVAReader->AddVariable("HoE",             &fMVAVar_HoE);
    tmpTMVAReader->AddVariable("EoP",             &fMVAVar_EoP); 
    tmpTMVAReader->AddVariable("IoEmIoP",         &fMVAVar_IoEmIoP);
    tmpTMVAReader->AddVariable("eleEoPout",       &fMVAVar_eleEoPout);
    tmpTMVAReader->AddVariable("rho",             &fMVAVar_rho);
     
    if(i == 2 || i == 5) 
      tmpTMVAReader->AddVariable("PreShowerOverRaw",&fMVAVar_PreShowerOverRaw);
      
    tmpTMVAReader->AddSpectator("eta",            &fMVAVar_eta);
    tmpTMVAReader->AddSpectator("pt",             &fMVAVar_pt);

    tmpTMVAReader->BookMVA(fMethodName , weightFiles[i] );
    cout << "MVABin " << i << " : MethodName = " << fMethodName 
              << "Load weight file : " << weightFiles[i] 
              << endl;
    fTMVAReader.push_back(tmpTMVAReader);
  }

}

//--------------------------------------------------------------------------------------------------
UInt_t HttElectronMVA::getMVABin(double eta, double pt) const {
  
  //Default is to return the first bin
  unsigned int bin = 0;

  if (pt < 20 && fabs(eta) < 0.8) bin = 0;
  if (pt < 20 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 1;
  if (pt < 20 && fabs(eta) >= 1.479) bin = 2;
  if (pt >= 20 && fabs(eta) < 0.8) bin = 3;
  if (pt >= 20 && fabs(eta) >= 0.8 && fabs(eta) < 1.479 ) bin = 4;
  if (pt >= 20 && fabs(eta) >= 1.479) bin = 5;

  return bin;

}

//--------------------------------------------------------------------------------------------------
Double_t HttElectronMVA::mvaValue(Double_t fbrem, 
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
					 Bool_t printDebug) {
  
  if (!fIsInitialized) { 
    cout << "Error: HttElectronMVA not properly initialized.\n"; 
    return -9999;
  }

  fMVAVar_fbrem           = fbrem; 
  fMVAVar_kfchi2          = kfchi2;
  fMVAVar_kfhits          = float(kfhits);   // BTD does not support int variables
  fMVAVar_gsfchi2         = gsfchi2;

  fMVAVar_deta            = deta;
  fMVAVar_dphi            = dphi;
  fMVAVar_detacalo        = detacalo;

  fMVAVar_see             = see;
  fMVAVar_spp             = spp;
  fMVAVar_etawidth        = etawidth;
  fMVAVar_phiwidth        = phiwidth;
  fMVAVar_OneMinusE1x5E5x5        = e1x5e5x5;
  fMVAVar_R9              = R9;

  fMVAVar_HoE             = HoE;
  fMVAVar_EoP             = EoP;
  fMVAVar_IoEmIoP         = IoEmIoP;
  fMVAVar_eleEoPout       = eleEoPout;
  fMVAVar_rho             = rho;
  fMVAVar_PreShowerOverRaw= PreShowerOverRaw;

  
  fMVAVar_eta             = eta;
  fMVAVar_pt              = pt;


  bindVariables();
  Double_t mva = -9999;  
  mva = fTMVAReader[getMVABin(fMVAVar_eta,fMVAVar_pt)]->EvaluateMVA(fMethodName);

  if(printDebug) {
    cout << " *** Inside the class fMethodName " << fMethodName << endl;
    cout << " fbrem " <<  fMVAVar_fbrem  
      	 << " kfchi2 " << fMVAVar_kfchi2  
	 << " mykfhits " << fMVAVar_kfhits  
	 << " gsfchi2 " << fMVAVar_gsfchi2  
	 << " deta " <<  fMVAVar_deta  
	 << " dphi " << fMVAVar_dphi  
      	 << " detacalo " << fMVAVar_detacalo  
	 << " see " << fMVAVar_see  
	 << " spp " << fMVAVar_spp  
	 << " etawidth " << fMVAVar_etawidth  
	 << " phiwidth " << fMVAVar_phiwidth  
	 << " e1x5e5x5 " << fMVAVar_OneMinusE1x5E5x5
	 << " R9 " << fMVAVar_R9  
	 << " HoE " << fMVAVar_HoE  
	 << " EoP " << fMVAVar_EoP  
	 << " IoEmIoP " << fMVAVar_IoEmIoP  
	 << " eleEoPout " << fMVAVar_eleEoPout 
	 << " rho " << fMVAVar_rho 
	 << " PreShowerOverRaw " << fMVAVar_PreShowerOverRaw  
       	 << " eta " << fMVAVar_eta  
	 << " pt " << fMVAVar_pt << endl;
    cout << " ### MVA " << mva << endl;
  }


  return mva;
}

void HttElectronMVA::bindVariables() {

  // this binding is needed for variables that sometime diverge. 


  if(fMVAVar_fbrem < -1.)
    fMVAVar_fbrem = -1.;	
  
  fMVAVar_deta = fabs(fMVAVar_deta);
  if(fMVAVar_deta > 0.06)
    fMVAVar_deta = 0.06;
  
  
  fMVAVar_dphi = fabs(fMVAVar_dphi);
  if(fMVAVar_dphi > 0.6)
    fMVAVar_dphi = 0.6;
  

  if(fMVAVar_EoP > 20.)
    fMVAVar_EoP = 20.;
  
  if(fMVAVar_eleEoPout > 20.)
    fMVAVar_eleEoPout = 20.;
  
  
  fMVAVar_detacalo = fabs(fMVAVar_detacalo);
  if(fMVAVar_detacalo > 0.2)
    fMVAVar_detacalo = 0.2;
  
  if(fMVAVar_OneMinusE1x5E5x5 < -1.)
    fMVAVar_OneMinusE1x5E5x5 = -1;
  
  if(fMVAVar_OneMinusE1x5E5x5 > 2.)
    fMVAVar_OneMinusE1x5E5x5 = 2.; 
  
  
  
  if(fMVAVar_R9 > 5)
    fMVAVar_R9 = 5;
  
  if(fMVAVar_gsfchi2 > 200.)
    fMVAVar_gsfchi2 = 200;
  
  
  if(fMVAVar_kfchi2 > 10.)
    fMVAVar_kfchi2 = 10.;
  
  
  // Needed for a bug in CMSSW_420, fixed in more recent CMSSW versions
  if(isnan(fMVAVar_spp))
    fMVAVar_spp = 0.;	
  
  
  return;
}
