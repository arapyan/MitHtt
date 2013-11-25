#ifndef DATAMC_HH
#define DATAMC_HH

#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1D.h"
#include "TH2D.h"
#include "BtagSF.hh"

enum { kNo, kDown, kUp };                     // systematic variation

// Get higgs mass point from sample name
Int_t higgsmass(TString basename)
{
  cout << basename << endl;
  stringstream ss(basename(TRegexp("[0-9][0-9]*"),3).Data());
  Int_t mass;
  ss >> mass;
  return mass;
}

TH2D* fEleTrigSF      = 0;
TH2D* fMuTrigSF      = 0;
TH2D* fTauTrigSF     = 0;
void setupTrigScale(Int_t is2012)
{
  TFile *lDFile = 0; 
  if(is2012)
    lDFile   = new TFile("$CMSSW_BASE/src/MitHtt/data/trigscale/weights.root");
  else
    lDFile   = new TFile("$CMSSW_BASE/src/MitHtt/data/trigscale/weights.root");
  assert(lDFile);
  fMuTrigSF  = (TH2D*) lDFile->FindObjectAny("IsoMuRatio");
  fMuTrigSF->SetDirectory(0);
  fEleTrigSF  = (TH2D*) lDFile->FindObjectAny("IsoMuRatio");
  fEleTrigSF->SetDirectory(0);
  fTauTrigSF = (TH2D*) lDFile->FindObjectAny("LooseIsoPFTauRatio");
  fTauTrigSF->SetDirectory(0);
  lDFile->Close();
} 

TH1D* kfFHPInit(Int_t mH);
Double_t kfFHPValue(Double_t pt, TH1D* hKF);

double efficiency(double m, double m0, double sigma, double alpha,double n, double norm);

Double_t unskimmedEntries(TString skimname)
{
  Double_t entries;

  skimname.ReplaceAll("_emu_skim.root","_ntuple.root");
  TFile unskimmed(skimname);
  assert(unskimmed.IsOpen());
  TTree *tree = 0;
  unskimmed.GetObject("hEvents",tree);
  assert(tree);
  entries = (Double_t)tree->GetEntries();
  unskimmed.Close();

  return entries;
}

TH1D* kfFHPInit(Int_t mH)
{
  TH1D *kfhist=0; 
  char kfilename[100];
  sprintf(kfilename, "$CMSSW_BASE/src/MitHtt/Utils/HiggsKFactors/weight_ptH_%d.root", mH);
  cout << "Getting k-factors from " << kfilename << endl;
  TFile *kfile = TFile::Open(kfilename); assert(kfile->IsOpen());
  TDirectory *kfdir = (TDirectory*)kfile->FindObjectAny("powheg_weight");
  char kfhistname[100];
  sprintf(kfhistname, "weight_hqt_fehipro_fit_%d", mH);
  cout << "kfactor histogram: " << kfhistname << endl;
  kfhist = (TH1D*)(kfdir->Get(kfhistname)); assert(kfhist);
  kfhist->SetDirectory(0);
  kfile->Close();
  return kfhist;
} 

//--------------------------------------------------------------------------------------------------
Double_t kfFHPValue(Double_t pt, TH1D* hKF)
{ 
  assert(hKF);
  return hKF->Interpolate(pt);
}

Bool_t passJetIDMVA(Double_t pt, Double_t eta, Double_t mvaVal)
{

  int lPtId = 0; 
  if(pt > 30) lPtId = 1;
  
  int lEtaId = 0;
  if(fabs(eta) > 2.5  && fabs(eta) < 2.75) lEtaId = 1; 
  if(fabs(eta) > 2.75 && fabs(eta) < 3.0 ) lEtaId = 2; 
  if(fabs(eta) > 3.0  && fabs(eta) < 5.0 ) lEtaId = 3; 

  double fMVACut[2][4];
  double Pt2030_Loose[4]    = {-0.63,-0.60,-0.55,-0.45};
  double Pt3050_Loose[4]    = {-0.63,-0.60,-0.55,-0.45};
  for(int i = 0; i < 4; i++) {
    fMVACut[0][i] = Pt2030_Loose[i];
    fMVACut[1][i] = Pt3050_Loose[i];
  }
  
  double lMVACut = fMVACut[lPtId][lEtaId];
  if(mvaVal < lMVACut) return false;
  return true;

}
//----------------------------------------------------------------
Double_t muIDscaleEmu(Double_t mupt, Double_t mueta, Int_t is2012)
{
  if((fabs(mueta) > 2.4) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  if(is2012) {
    if(mupt > 35) {
      if(fabs(mueta) < 0.8)   return 0.9841;
      else if(fabs(mueta) < 1.2)   return 0.9813;
      else if(fabs(mueta) < 1.6)   return 0.9919;
      else if(fabs(mueta) < 2.1)   return 0.9939;
      else                    return 1.0401;
    }
    else if(mupt > 30) {
      if(fabs(mueta) < 0.8)   return 0.9746;
      else if(fabs(mueta) < 1.2)   return 0.9797;
      else if(fabs(mueta) < 1.6)   return 0.9935;
      else if(fabs(mueta) < 2.1)   return 0.9987;
      else                    return 1.0785;
    }
    else if(mupt > 25) {
      if(fabs(mueta) < 0.8)   return 0.9691;
      else if(fabs(mueta) < 1.2)   return 0.9785;
      else if(fabs(mueta) < 1.6)   return 0.9909;
      else if(fabs(mueta) < 2.1)   return 0.9991;
      else                    return 1.0922;
    }
    else if(mupt > 20) {
      if(fabs(mueta) < 0.8)   return 0.9676;
      else if(fabs(mueta) < 1.2)   return 0.9785;
      else if(fabs(mueta) < 1.6)   return 0.9883;
      else if(fabs(mueta) < 2.1)   return 1.0031;
      else                    return 1.1362;
    }
    else if(mupt > 15) {
      if(fabs(mueta) < 0.8)   return 0.9556;
      else if(fabs(mueta) < 1.2)   return 0.9635;
      else if(fabs(mueta) < 1.6)   return 0.9806;
      else if(fabs(mueta) < 2.1)   return 1.0078;
      else                    return 1.1422;
    }
    else {
      if(fabs(mueta) < 0.8)   return 0.9811;
      else if(fabs(mueta) < 1.2)   return 0.9689;
      else if(fabs(mueta) < 1.6)   return 0.9757;
      else if(fabs(mueta) < 2.1)   return 1.0069;
      else                    return 1.1833;
    }
  } else {
    if(mupt > 20) {
      if(fabs(mueta) < 0.8)      return 0.9998;
      else if(fabs(mueta) < 1.2) return 1.0006;
      else                       return 1.0045;
    }
    else if(mupt > 15) {
      if(fabs(mueta) < 0.8)      return 1.0176;
      else if(fabs(mueta) < 1.2) return 1.0040;
      else                       return 1.0063;
    }
    else {
      if(fabs(mueta) < 0.8)      return 0.9303;
      else if(fabs(mueta) < 1.2) return 1.0125;
      else                       return 0.9994;
    }
  }
}

//----------------------------------------------------------------------------------------    
Double_t eleIDscaleEmu(Double_t elept, Double_t eleeta, Int_t is2012)
{
  if((fabs(eleeta) > 2.3) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  if(is2012) {
    if(elept > 35) {
      if(fabs(eleeta) < 0.8)        return 0.9533;                                                                                
      else if(fabs(eleeta) < 1.479) return 0.9496;
      else                          return 0.9389;                                                                                
    }     
    else if(elept > 30) {                                                                                                                
      if(fabs(eleeta) < 0.8)        return 0.9301;                                              
      else if(fabs(eleeta) < 1.479) return 0.9230;
      else                          return 0.8887;                                                
    }       
    else if(elept > 25) { 
      if(fabs(eleeta) < 0.8)        return 0.9069;                                                                                
      else if(fabs(eleeta) < 1.479) return 0.8896;                                         
      else                          return 1.0225;                                                                                
    }       
    else if(elept > 20) {
      if(fabs(eleeta) < 0.8)        return 0.8817;                                        
      else if(fabs(eleeta) < 1.479) return 0.8492;                                                                                
      else                          return 0.8057; 
    }         
    else if(elept > 15) {                      
      if(fabs(eleeta) < 0.8)        return 0.8437;                                                                                
      else if(fabs(eleeta) < 1.479) return 0.8447;                                              
      else                          return 0.7812;
    }
    else {
      if(fabs(eleeta) < 0.8)        return 0.7570;
      else if(fabs(eleeta) < 1.479) return 0.7807;
      else                          return 0.6276;
    } 
  } else {
    if(elept > 20) {
      if(fabs(eleeta) < 0.8)        return 0.9862;
      else if(fabs(eleeta) < 1.479) return 0.9786;
      else                          return 1.0136;
    }
    else if(elept > 15) {
      if(fabs(eleeta) < 0.8)        return 0.9612;
      else if(fabs(eleeta) < 1.479) return 0.9773;
      else                          return 1.0600;
    }
    else {
      if(fabs(eleeta) < 0.8)        return 1.0078;
      else if(fabs(eleeta) < 1.479) return 1.1236;
      else                          return 0.9336;
    }
  }
}

//----------------------------------------------------------------
Double_t muIDIsoscaleMuTau(Double_t mupt, Double_t mueta, Int_t is2012)
{
  if((fabs(mueta) > 2.1) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  if(mupt > 30) {
    if     (fabs(mueta) < 0.8)  return is2012 ? 0.9852*0.9883  : 0.9977*0.9895;
    else if(fabs(mueta) < 1.2)  return is2012 ? 0.9852*0.9937  : 0.9893*0.9936;
    else                        return is2012 ? 0.9884*0.9996  : 0.9829*0.9960;
  }
  else if(mupt > 20) {
    if     (fabs(mueta) < 0.8)  return is2012 ? 0.9818*0.9494  : 0.9962*1.0011;
    else if(fabs(mueta) < 1.2)  return is2012 ? 0.9829*0.9835  : 0.9904*0.9834;
    else                        return is2012 ? 0.9869*0.9923  : 0.9828*0.9975;
  }
  else {
    if     (fabs(mueta) < 0.8)  return 0.9963*0.9910;
    else if(fabs(mueta) < 1.2)  return 0.9846*0.9643;
    else                        return 0.9830*0.9504;
  }
}
//----------------------------------------------------------------
Double_t muIDscaleMuTau(Double_t mupt, Double_t mueta, Int_t is2012)
{
  if((fabs(mueta) > 2.1) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  if(mupt > 30) {
    if     (fabs(mueta) < 0.8)  return is2012 ? 0.9852  : 0.9977;
    else if(fabs(mueta) < 1.2)  return is2012 ? 0.9852 : 0.9893;
    else                        return is2012 ? 0.9884  : 0.9829;
  }
  else if(mupt > 20) {
    if     (fabs(mueta) < 0.8)  return is2012 ? 0.9818 : 0.9962;
    else if(fabs(mueta) < 1.2)  return is2012 ? 0.9829 : 0.9904;
    else                        return is2012 ? 0.9869 : 0.9828;
  }
  else {
    if     (fabs(mueta) < 0.8)  return 0.9963;
    else if(fabs(mueta) < 1.2)  return 0.9846;
    else                        return 0.9830;
  }
}


//----------------------------------------------------------------
Double_t eleIDIsoscaleETau(Double_t ept, Double_t eeta, Int_t is2012)
{
  if((fabs(eeta) > 2.1) || (ept < 10)) { cout << "electron kinematics out of range" << endl; assert(0); }
  if(ept > 30) {
    if(fabs(eeta) < 1.479)  return is2012 ? 0.9486*0.9804 : 0.9826*0.9845;
    else                    return is2012 ? 0.8866*0.99 : 0.9689*0.9971;
  }
  else {
    if(fabs(eeta) < 1.479) return is2012 ? 0.8999*0.9417 : 0.9590*0.9907;
    else                   return is2012 ? 0.7945*0.9471 : 0.9462*0.9875;
  }
}
//----------------------------------------------------------------
Double_t eleIDScaleETau(Double_t ept, Double_t eeta, Int_t is2012)
{
  if((fabs(eeta) > 2.1) || (ept < 10)) { cout << "electron kinematics out of range" << endl; assert(0); }
  if(ept > 30) {
    if(fabs(eeta) < 1.479)  return is2012 ? 0.9486 : 0.9826;
    else                    return is2012 ? 0.8866 : 0.9689;
  }
  else {
    if(fabs(eeta) < 1.479) return is2012 ? 0.8999 : 0.9590;
    else                   return is2012 ? 0.7945 : 0.9462;
  }
}
//----------------------------------------------------------------------------------------
Double_t eleTrigScaleEmu(Double_t elept, Double_t eleeta, Int_t is2012)
{
  if((fabs(eleeta) > 2.3) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  if(is2012) {
    if(elept > 35) {
      if(fabs(eleeta) < 0.8)        return 1.0069;
      else if(fabs(eleeta) < 1.479) return 1.0049;
      else                          return 0.9989;
    } 
    else if(elept > 30) {
      if(fabs(eleeta) < 0.8)        return 1.0084;
      else if(fabs(eleeta) < 1.479) return 0.9900;
      else                          return 0.9817;
    } 
    else if(elept > 25) {
      if(fabs(eleeta) < 0.8)        return 0.9772;
      else if(fabs(eleeta) < 1.479) return 0.9916;
      else                          return 0.9609;
    } 
    else if(elept > 20) {
      if(fabs(eleeta) < 0.8)        return 0.9716;
      else if(fabs(eleeta) < 1.479) return 0.9702;
      else                          return 0.9726;
    }
    else if(elept > 15) {
      if(fabs(eleeta) < 0.8)        return 0.9841;
      else if(fabs(eleeta) < 1.479) return 0.9699;
      else                          return 0.9286;
    }
    else {
      if(fabs(eleeta) < 0.8)        return 0.9529;
      else if(fabs(eleeta) < 1.479) return 0.8858;
      else                          return 0.9259;
    }
  } else {
    if(elept > 30) {
      if(fabs(eleeta) < 0.8)        return 0.99;
      else if(fabs(eleeta) < 1.479) return 1.00;
      else                          return 0.99;
    }
    else if(elept > 20) {
      if(fabs(eleeta) < 0.8)        return 0.98;
      else if(fabs(eleeta) < 1.479) return 0.99;
      else                          return 0.96;
    }
    else if(elept > 15) {
      if(fabs(eleeta) < 0.8)        return 1.10;
      else if(fabs(eleeta) < 1.479) return 1.08;
      else                          return 1.05;
    } 
    else {
      if(fabs(eleeta) < 0.8)        return 0.97;
      else if(fabs(eleeta) < 1.479) return 0.98;
      else                          return 0.81;
    }
  }
}

//----------------------------------------------------------------------------------------
Double_t muTrigScaleEmu(Double_t mupt, Double_t mueta, Int_t is2012)
{
  if((fabs(mueta) > 2.4) || (mupt < 10)) { cout << "muon kinematics out of range" << endl; assert(0); }
  if(is2012) {
    if(mupt > 35) {
      if(fabs(mueta) < 0.8)        return 0.9991;
      else if(fabs(mueta) < 1.2)   return 0.9626;
      else if(fabs(mueta) < 1.6)   return 0.9611;
      else if(fabs(mueta) < 2.1)   return 0.9314;
      else                         return 0.8541;
    }
    else if(mupt > 30) {
      if(fabs(mueta) < 0.8)        return 0.9930;
      else if(fabs(mueta) < 1.2)   return 0.9800;
      else if(fabs(mueta) < 1.6)   return 0.9958;
      else if(fabs(mueta) < 2.1)   return 0.9428;
      else                         return 0.7810;
    }
    else if(mupt > 25) {
      if(fabs(mueta) < 0.8)        return 0.9856;
      else if(fabs(mueta) < 1.2)   return 0.9818;
      else if(fabs(mueta) < 1.6)   return 0.9684;
      else if(fabs(mueta) < 2.1)   return 0.9642;
      else                         return 0.8370;
    }
    else if(mupt > 20) {
      if(fabs(mueta) < 0.8)        return 0.9937;
      else if(fabs(mueta) < 1.2)   return 0.9594;
      else if(fabs(mueta) < 1.6)   return 0.9692;
      else if(fabs(mueta) < 2.1)   return 0.9438;
      else                         return 0.7761;
    }
    else if(mupt > 15) {
      if(fabs(mueta) < 0.8)        return 0.9846;
      else if(fabs(mueta) < 1.2)   return 0.9834;
      else if(fabs(mueta) < 1.6)   return 0.9793;
      else if(fabs(mueta) < 2.1)   return 0.9257;
      else                         return 0.8199;
    }
    else {
      if(fabs(mueta) < 0.8)        return 0.9841;
      else if(fabs(mueta) < 1.2)   return 0.9742;
      else if(fabs(mueta) < 1.6)   return 0.9955;
      else if(fabs(mueta) < 2.1)   return 0.9151;
      else                         return 0.8067;
    }
  } else {
    if(mupt > 30) {
      if(fabs(mueta) < 0.8)        return 0.9783;
      else if(fabs(mueta) < 1.2)   return 0.9669;
      else                         return 0.9674;
    }
    else if(mupt > 25) {
      if(fabs(mueta) < 0.8)        return 0.9812;
      else if(fabs(mueta) < 1.2)   return 0.9783;
      else                         return 0.9700;
    }
    else if(mupt > 20) {
      if(fabs(mueta) < 0.8)        return 1.0035;
      else if(fabs(mueta) < 1.2)   return 0.9690;
      else                         return 0.9694;
    }
    else if(mupt > 15) {
      if(fabs(mueta) < 0.8)        return 0.9808;
      else if(fabs(mueta) < 1.2)   return 0.9791;
      else                         return 0.9948;
    }
    else {
      if(fabs(mueta) < 0.8)        return 0.9808;
      else if(fabs(mueta) < 1.2)   return 0.9623;
      else                         return 0.9602;
    }
  }
}

//----------------------------------------------------------------------------------------
Double_t eleTrigEffEmu(Double_t elept, Double_t eleeta, Int_t is2012)
{
  if((fabs(eleeta) > 2.3) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  if(is2012) {
    if(elept > 35) {
      if(fabs(eleeta) < 0.8)        return 0.9690;
      else if(fabs(eleeta) < 1.479) return 0.9811;
      else                          return 0.9777;
    }
    else if(elept > 30) {
      if(fabs(eleeta) < 0.8)        return 0.9528;
      else if(fabs(eleeta) < 1.479) return 0.9652;
      else                          return 0.9693;
    }
    else if(elept > 25) {
      if(fabs(eleeta) < 0.8)        return 0.9383;
      else if(fabs(eleeta) < 1.479) return 0.9681;
      else                          return 0.9383;
    }
    else if(elept > 20) {
      if(fabs(eleeta) < 0.8)        return 0.9145;
      else if(fabs(eleeta) < 1.479) return 0.9455;
      else                          return 0.9354;
    }
    else if(elept > 15) {
      if(fabs(eleeta) < 0.8)        return 0.8758;
      else if(fabs(eleeta) < 1.479) return 0.9061;
      else                          return 0.8475;
    }
    else {
      if(fabs(eleeta) < 0.8)        return 0.7217;
      else if(fabs(eleeta) < 1.479) return 0.7262;
      else                          return 0.7093;
    }
  } else {
    if(elept > 30) {
      if(fabs(eleeta) < 0.8)        return 0.9643;
      else if(fabs(eleeta) < 1.479) return 0.9778;
      else                          return 0.9737;
    }
    else if(elept > 25) {
      if(fabs(eleeta) < 0.8)        return 0.9394;
      else if(fabs(eleeta) < 1.479) return 0.9674;
      else                          return 0.9286;
    }
    else if(elept > 20) {
      if(fabs(eleeta) < 0.8)        return 0.9200;
      else if(fabs(eleeta) < 1.479) return 0.9515;
      else                          return 0.9323;
    }
    else if(elept > 15) {
      if(fabs(eleeta) < 0.8)        return 0.8874;
      else if(fabs(eleeta) < 1.479) return 0.9177;
      else                          return 0.8500;
    }
    else {
      if(fabs(eleeta) < 0.8)        return 0.7633;
      else if(fabs(eleeta) < 1.479) return 0.7356;
      else                          return 0.7010;
    }
  }
}

//----------------------------------------------------------------------------------------
Double_t muTrigEffEmu(Double_t mupt, Double_t mueta, Int_t is2012)
{
  if((fabs(mueta) > 2.4) || (mupt < 10)) { cout << "muon kinematics out of range" << endl; assert(0); }
  if(is2012) {
    if(mupt > 35) {
      if(fabs(mueta) < 0.8)        return 0.9683;
      else if(fabs(mueta) < 1.2)   return 0.9345;
      else if(fabs(mueta) < 1.6)   return 0.9098;
      else if(fabs(mueta) < 2.1)   return 0.9003;
      else                         return 0.7068;
    }
    else if(mupt > 30) {
      if(fabs(mueta) < 0.8)        return 0.9758;
      else if(fabs(mueta) < 1.2)   return 0.9314;
      else if(fabs(mueta) < 1.6)   return 0.9152;
      else if(fabs(mueta) < 2.1)   return 0.9019;
      else                         return 0.7278;
    }
    else if(mupt > 25) {
      if(fabs(mueta) < 0.8)        return 0.9711;
      else if(fabs(mueta) < 1.2)   return 0.9439;
      else if(fabs(mueta) < 1.6)   return 0.9222;
      else if(fabs(mueta) < 2.1)   return 0.8909;
      else                         return 0.7609;
    }
    else if(mupt > 20) {
      if(fabs(mueta) < 0.8)        return 0.9749;
      else if(fabs(mueta) < 1.2)   return 0.9423;
      else if(fabs(mueta) < 1.6)   return 0.9375;
      else if(fabs(mueta) < 2.1)   return 0.9114;
      else                         return 0.6909;
    }
    else if(mupt > 15) {
      if(fabs(mueta) < 0.8)        return 0.9706;
      else if(fabs(mueta) < 1.2)   return 0.9284;
      else if(fabs(mueta) < 1.6)   return 0.9306;
      else if(fabs(mueta) < 2.1)   return 0.8921;
      else                         return 0.7389;
    }
    else {
      if(fabs(mueta) < 0.8)        return 0.9713;
      else if(fabs(mueta) < 1.2)   return 0.9399;
      else if(fabs(mueta) < 1.6)   return 0.9299;
      else if(fabs(mueta) < 2.1)   return 0.8614;
      else                         return 0.6832;
    }
  } else {
    if(mupt > 30) {
      if(fabs(mueta) < 0.8)        return 0.9671;
      else if(fabs(mueta) < 1.2)   return 0.9502;
      else                         return 0.9383;
    }
    else if(mupt > 25) {
      if(fabs(mueta) < 0.8)        return 0.9680;
      else if(fabs(mueta) < 1.2)   return 0.9550;
      else                         return 0.9416;
    }
    else if(mupt > 20) {
      if(fabs(mueta) < 0.8)        return 0.9878;
      else if(fabs(mueta) < 1.2)   return 0.9495;
      else                         return 0.9379;
    }
    else if(mupt > 15) {
      if(fabs(mueta) < 0.8)        return 0.9668;
      else if(fabs(mueta) < 1.2)   return 0.9556;
      else                         return 0.9613;
    }
    else {
      if(fabs(mueta) < 0.8)        return 0.9660;
      else if(fabs(mueta) < 1.2)   return 0.9314;
      else                         return 0.9207;
    }
  }
}

//----------------------------------------------------------------------------------------
double efficiency(double m, double m0, double sigma, double alpha,
                  double n, double norm) {
  
  const double sqrtPiOver2 = 1.2533141373;
  const double sqrt2 = 1.4142135624;
  double sig = fabs((double) sigma);
  double t = (m - m0)/sig;
  if(alpha < 0)   t = -t;
  double absAlpha = fabs(alpha/sig);
  double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
  double b = absAlpha - n/absAlpha;
  double ApproxErf;
  double arg = absAlpha / sqrt2;
  if (arg > 5.) ApproxErf = 1;
  else if (arg < -5.) ApproxErf = -1;
  else ApproxErf = TMath::Erf(arg);
  double leftArea = (1 + ApproxErf) * sqrtPiOver2;
  double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
  double area = leftArea + rightArea;
  if( t <= absAlpha ){
    arg = t / sqrt2;
    if(arg > 5.) ApproxErf = 1;
    else if (arg < -5.) ApproxErf = -1;
    else ApproxErf = TMath::Erf(arg);
    return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
  }
  else{
    return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) -
                                   1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
  }
}
//2011 ltau triger efficiency
double eff2011TrigEle(double e_pt, double eta,bool turnon)
{
  double ele_trg_data=1.0;
  double ele_trg_mc=1.0;
  if(fabs(eta)<1.479)
    {
      double ele15_eb = efficiency(e_pt, 14.8474, 0.634431, 0.547211, 2.11286, 0.978835);
      double ele18_eb = efficiency(e_pt, 18.318, 0.811329, 0.938142, 1.85834, 0.927328);
      double ele20_eb = efficiency(e_pt, 20.244, 0.783044, 0.672118, 2.35819, 0.972123);
      ele_trg_mc = efficiency(e_pt, 16.9487, 0.162306, 0.0506008, 2.37608, 0.973829);
      ele_trg_data = (0.41 * ele15_eb) + (0.39 * ele18_eb) + (0.2 * ele20_eb);
    }
  else
    {
      double ele15_ee = efficiency(e_pt, 16.47, 1.20691, 1.54605, 1.68873, 0.985916);
      double ele18_ee = efficiency(e_pt, 16.9487, 0.162306, 0.0506008, 2.37608, 0.973829);
      double ele20_ee = efficiency(e_pt, 22.4559, 1.70818, 0.917021, 3.86463, 0.964978);
      ele_trg_mc = efficiency(e_pt, 18.05, 1.66609, 2.3836, 1.49316, 1.01229);
      ele_trg_data = (0.41 * ele15_ee) + (0.39 * ele18_ee) + (0.2 * ele20_ee);
    }
  if(turnon)
    return ele_trg_data;
  else
    return ele_trg_data/ele_trg_mc;
}
double eff2011TrigMu(double pt, double eta,bool turnon)
{
  double mu_trg_data=1.0;
  double mu_trg_mc=1.0;
  if(fabs(eta)<0.8)
    {
      double mu12 = 0.920;
      double mu15 = 0.917;
      double mu15_2p1 = efficiency(pt, 15.9877, 2.90938e-07, 2.63922e-11, 5.81194, 0.906943);
      mu_trg_mc = 0.923;
      mu_trg_data = (0.034 * mu12) + (0.368 * mu15) + (0.598 * mu15_2p1);
    }
  else if(fabs(eta)>0.8 && fabs(eta)<1.2)
    {
      double mu12 = 0.868;
      double mu15 = 0.871;
      double mu15_2p1 = efficiency(pt, 15.9995, 1.35931e-07, 7.88264e-11, 4.60253, 0.855461);
      mu_trg_data = (0.034 * mu12) + (0.368 * mu15) + (0.598 * mu15_2p1);
      mu_trg_mc = 0.879;
    }
  else
    {
      double mu12 = 0.845;
      double mu15 = 0.864;
      double mu15_2p1 = efficiency(pt, 15.9084, 2.27242e-12, 8.77174e-14, 1.00241, 12.9909);
      mu_trg_data = (0.034 * mu12) + (0.368 * mu15) + (0.598 * mu15_2p1);
      mu_trg_mc = 0.839;
    }
  if(turnon)
    return mu_trg_data;
  else
    return mu_trg_data/mu_trg_mc;
}
double eff2011TrigMuTau(double t_pt, double eta,bool turnon){
  double tau_trg_data = 1.0;
  double tau_trg_mc = 1.0;
  if(fabs(eta)<1.5)
    {
      double tau10l_eb = efficiency(t_pt, 13.6046,   1.66291,   1.71551,   141.929,   0.910686);
      double tau15l_eb = efficiency(t_pt, 13.9694,   0.084835,  0.057743,  1.50674,   0.984976);
      double tau20l_eb = efficiency(t_pt, 19.2102,   1.26519,   2.48994,   1.04699,  1.3492);
      tau_trg_data = (0.043 * tau10l_eb) + (0.359 * tau15l_eb) + (0.598 * tau20l_eb);
      tau_trg_mc = efficiency(t_pt, 14.4601, 0.0485272, 0.03849, 1.48324, 0.965257);
    }
  else
    {
      double tau10l_ee = efficiency(t_pt, -0.392211,   7.90467,   5.48228,   134.599,   0.925858);
      double tau15l_ee = efficiency(t_pt, 14.435,  1.34952,   2.43996,   1.03631,   1.79081);
      double tau20l_ee = efficiency(t_pt, 19.2438,   1.37298,   1.76448,   1.73935,   0.901291);
      tau_trg_data = (0.043 * tau10l_ee) + (0.359 * tau15l_ee) + (0.598 * tau20l_ee);
      tau_trg_mc = efficiency(t_pt, 14.4451, 0.0790573, 0.0732472, 1.47046, 0.942028);
    }
  if(turnon)
    return tau_trg_data;
  else
    return tau_trg_data/tau_trg_mc;
}
double eff2011TrigEleTau(double t_pt, double eta,bool turnon){
  double tau_trg_data = 1.0;
  double tau_trg_mc = 1.0;
  if(fabs(eta)<1.5)
    {
      double tau20l_eb = efficiency(t_pt, 19.3916,  0.996964,  1.70131,   1.38002,   0.903245);
      double tau20m_eb = efficiency(t_pt, 19.5667,  1.15203 ,  1.68126,   1.40025,   0.848033);
      double tau20t_eb = efficiency(t_pt, 19.6013,  0.987317,  1.08015,   1.88592,   0.776894);
      tau_trg_mc = efficiency(t_pt, 19.468, 0.0615381, 0.0349325, 1.59349, 0.860096);
      tau_trg_data = (0.25 * tau20l_eb) + (0.59 * tau20m_eb) + (0.16 * tau20t_eb);
    }
  else
    {
      double tau20l_ee = efficiency(t_pt, 18.8166,  0.526632,  0.20666,   6.80392,   0.903245);
      double tau20m_ee = efficiency(t_pt, 18.8476,  0.528963,  0.16717,   3.65814,   0.749759);
      double tau20t_ee = efficiency(t_pt, 18.8859,  0.271301,  0.128008,  1.50993,   0.825122);
      tau_trg_mc = efficiency(t_pt, 19.3862, 0.247148, 0.123187, 2.87108, 0.790894);
      tau_trg_data = (0.25 * tau20l_ee) + (0.59 * tau20m_ee) + (0.16 * tau20t_ee);
    }
  if(turnon)
    return tau_trg_data;
  else
    return tau_trg_data/tau_trg_mc;
}
double eff2012IsoTau12fb(double pt, double eta){
  return (808.411*(0.764166*0.5*(TMath::Erf((pt-33.2236)/2./0.97289/sqrt(pt))+1.))+
	  4428.0*(0.802387*0.5*(TMath::Erf((pt-38.0971)/2./0.82842/sqrt(pt))+1.))+
	  1783.003*(0.818051*0.5*(TMath::Erf((pt-37.3669)/2./0.74847/sqrt(pt))+1.))+
	  5109.155*(0.796086*0.5*(TMath::Erf((pt-37.3302)/2./0.757558/sqrt(pt))+1.))
	  )/(808.411+4428.0+1783.003+5109.155);
}
  
double eff2012Jet12fb(double pt, double eta){
  return (abs(eta)<=2.1)*
    ((808.411*(0.99212*0.5*(TMath::Erf((pt-31.3706)/2./1.22821/sqrt(pt))+1.))+
      4428.0*(0.99059*0.5*(TMath::Erf((pt-32.1104)/2./1.23292/sqrt(pt))+1.))+
      1783.003*(0.988256*0.5*(TMath::Erf((pt-31.3103)/2./1.18766/sqrt(pt))+1.))+
      5109.155*(0.988578*0.5*(TMath::Erf((pt-31.6391)/2./1.22826/sqrt(pt))+1.))
      )/(808.411+4428.0+1783.003+5109.155))+
    (abs(eta)>2.1)*
    ((808.411*(0.969591*0.5*(TMath::Erf((pt-36.8179)/2./0.904254/sqrt(pt))+1.))+
      4428.0*(0.975932*0.5*(TMath::Erf((pt-37.2121)/2./0.961693/sqrt(pt))+1.))+
      1783.003*(0.990305*0.5*(TMath::Erf((pt-36.3096)/2./0.979524/sqrt(pt))+1.))+
      5109.155*(0.971612*0.5*(TMath::Erf((pt-36.2294)/2./0.871726/sqrt(pt))+1.))
      )/(808.411+4428.0+1783.003+5109.155));
}

double eff2012Jet19fb(double pt, double eta){
return (abs(eta)<=2.1)*
((808.411*(0.99212*0.5*(TMath::Erf((pt-31.3706)/2./1.22821/sqrt(pt))+1.))+
4428.0*(0.99059*0.5*(TMath::Erf((pt-32.1104)/2./1.23292/sqrt(pt))+1.))+
1783.003*(0.988256*0.5*(TMath::Erf((pt-31.3103)/2./1.18766/sqrt(pt))+1.))+
5109.155*(0.988578*0.5*(TMath::Erf((pt-31.6391)/2./1.22826/sqrt(pt))+1.))+
4131*(0.989049*0.5*(TMath::Erf((pt-31.9836)/2./1.23871/sqrt(pt))+1.))+
3143*(0.988047*0.5*(TMath::Erf((pt-31.6975)/2./1.25372/sqrt(pt))+1.))
)/(808.411+4428.0+1783.003+5109.155+4131+3143))+
(abs(eta)>2.1)*
((808.411*(0.969591*0.5*(TMath::Erf((pt-36.8179)/2./0.904254/sqrt(pt))+1.))+
4428.0*(0.975932*0.5*(TMath::Erf((pt-37.2121)/2./0.961693/sqrt(pt))+1.))+
1783.003*(0.990305*0.5*(TMath::Erf((pt-36.3096)/2./0.979524/sqrt(pt))+1.))+
5109.155*(0.971612*0.5*(TMath::Erf((pt-36.2294)/2./0.871726/sqrt(pt))+1.))+
4131*(0.977958*0.5*(TMath::Erf((pt-37.131)/2./0.987523/sqrt(pt))+1.))+
3143*(0.968457*0.5*(TMath::Erf((pt-36.3159)/2./0.895031/sqrt(pt))+1.))
)/(808.411+4428.0+1783.003+5109.155+4131+3143));
}

double eff2012IsoTau19fb(double pt, double eta){
return (808.411*(0.764166*0.5*(TMath::Erf((pt-33.2236)/2./0.97289/sqrt(pt))+1.))+
4428.0*(0.802387*0.5*(TMath::Erf((pt-38.0971)/2./0.82842/sqrt(pt))+1.))+
1783.003*(0.818051*0.5*(TMath::Erf((pt-37.3669)/2./0.74847/sqrt(pt))+1.))+
5109.155*(0.796086*0.5*(TMath::Erf((pt-37.3302)/2./0.757558/sqrt(pt))+1.))+
4131*(0.828182*0.5*(TMath::Erf((pt-37.6596)/2./0.830682/sqrt(pt))+1.))+
3143*(0.833004*0.5*(TMath::Erf((pt-37.634)/2./0.777843/sqrt(pt))+1.))
)/(808.411+4428.0+1783.003+5109.155+4131+3143);
}


// Tau Parked with HLT_DoubleMediumIsoPFTau35_Trk*_eta2p1_v*
double eff2012IsoParkedTau19fb_Simone(double pt, double eta){
  
  // for real Taus mT<20
  if ( fabs(eta) < 1.4 )
    {
      return (  0.883869 * 0.5 * (TMath::Erf((pt-43.8723)/2./0.946593 /sqrt(pt))+1.) ) ; 
      
    }
  
    else
      {
	return (  0.798480 * 0.5 * (TMath::Erf((pt-43.1362)/2./1.04861  /sqrt(pt))+1.) ) ;
      }
  
}
 

double eff2012IsoParkedTau19fbMC(double pt, double eta)
  {
    return ( 0.813769 * 0.5 * (TMath::Erf((pt-39.9322)/2./0.819354  /sqrt(pt))+1.) ) ;
  }

double eff2012IsoParkedTau19fbData(double pt, double eta)
  {
    return ( 0.826969 * 0.5 * (TMath::Erf((pt-42.2274)/2./0.783258  /sqrt(pt))+1.) ) ;
  }

double eff2012IsoParkedTau19fbMCMSSM(double pt, double eta)
  {
    if(pt<140)
      return ( 0.813769 * 0.5 * (TMath::Erf((pt-39.9322)/2./0.819354  /sqrt(pt))+1.) ) ;
    else
      return 1.0;
  }

double eff2012IsoParkedTau19fbDataMSSM(double pt, double eta)
  {
    if(pt<140)
      return ( 0.826969 * 0.5 * (TMath::Erf((pt-42.2274)/2./0.783258  /sqrt(pt))+1.) ) ;
    else if(pt>400) return 2.03467;
    else if (pt>300) return 1.31593;
    else if (pt>250) return 1.25698;
    else if (pt>200) return 1.18941;
    else if (pt>180) return 1.17448;
    else if (pt>160) return 1.0964;
    else return 1.09279;
  }
  
double eff2012IsoTau19fbData(double pt, double eta){

    // for real Taus mT<20
    if ( fabs(eta) < 1.4 )
    {
      return (  808.411  * ( 0.764166 * 0.5 * (TMath::Erf((pt-33.2236)/2./0.97289 /sqrt(pt))+1.))   // 2012A by Bastian not split in eta
              + 4428.0   * ( 0.75721  * 0.5 * (TMath::Erf((pt-39.0836)/2./1.07753 /sqrt(pt))+1.))   // 2012B
              + 6892.158 * ( 0.791464 * 0.5 * (TMath::Erf((pt-38.4932)/2./1.01232 /sqrt(pt))+1.))   // 2012C measured in v2 only
              + 7274.    * ( 0.779446 * 0.5 * (TMath::Erf((pt-38.4603)/2./1.01071 /sqrt(pt))+1.)) ) // 2012D measured in one go
              /( 808.411 + 4428.0 + 6892.158 + 7274. );
    }
    
    else
    {
      return (  808.411  * ( 0.764166 * 0.5 * (TMath::Erf((pt-33.2236)/2./0.97289 /sqrt(pt))+1.))   // 2012A by Bastian not split in eta
              + 4428.0   * ( 0.693788 * 0.5 * (TMath::Erf((pt-37.7719)/2./1.09202 /sqrt(pt))+1.))   // 2012B
              + 6892.158 * ( 0.698909 * 0.5 * (TMath::Erf((pt-36.5533)/2./1.05743 /sqrt(pt))+1.))   // 2012C measured in v2 only
              + 7274.    * ( 0.703532 * 0.5 * (TMath::Erf((pt-38.8609)/2./1.05514 /sqrt(pt))+1.)) ) // 2012D measured in one go
              /( 808.411 + 4428.0 + 6892.158 + 7274. );
    }
    
  }

  double eff2012IsoTau19fbMC(double pt, double eta){

    // for real Taus using ggH120
    if ( fabs(eta) < 1.4 )
    {
      return ( 0.807425 * 0.5 * (TMath::Erf((pt-35.2214)/2./1.04214  /sqrt(pt))+1.) ) ;
    }
    
    else
    {
      return ( 0.713068 * 0.5 * (TMath::Erf((pt-33.4584)/2./0.994692 /sqrt(pt))+1.) ) ;
    }
    
  }




#endif
