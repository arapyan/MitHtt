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

Double_t unskimmedEntries(TString skimname);
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
Double_t muIDscaleMuTau(Double_t mupt, Double_t mueta, Int_t is2012)
{
  if((fabs(mueta) > 2.1) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  if(mupt > 30) {
    if     (fabs(mueta) < 0.8)  return is2012 ? 0.9872*0.9857  : 1.030*1.010;
    else if(fabs(mueta) < 1.2)  return is2012 ? 0.9924*0.9805  : 0.997*0.990;
    else                        return is2012 ? 1.0012*0.9900  : 0.997*0.990;
  }
  else if(mupt > 20) {
    if     (fabs(mueta) < 0.8)  return is2012 ? 0.9685*0.9853  : 0.977*0.995;
    else if(fabs(mueta) < 1.2)  return is2012 ? 0.9808*0.9818  : 0.984*0.986;
    else                        return is2012 ? 0.9972*0.9899  : 0.984*0.986;
  }
  else if(mupt > 17 && !is2012) {
    if(fabs(mueta) < 1.)   return 0.997*0.930;
    else                  return 0.986*0.929;
  }
  else {
    if(fabs(mueta) < 1.6)   return 0.945*0.989;
    else                 return 0.977*1.047;
  }
}

//----------------------------------------------------------------
Double_t eleIDscaleETau(Double_t ept, Double_t eeta, Int_t is2012)
{
  if((fabs(eeta) > 2.1) || (ept < 10)) { cout << "electron kinematics out of range" << endl; assert(0); }
  if(ept > 30) {
    if(fabs(eeta) < 1.479)  return is2012 ? 0.982*0.949 : 1.044*0.984;
    else                    return is2012 ? 0.995*0.926 : 0.989*0.977;
  }
  else {
    if(fabs(eeta) < 1.479) return is2012 ? 0.947*0.9100 : 0.955*0.980 ;
    else                   return is2012 ? 0.959*0.8244 : 0.967*0.938;
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

#endif
