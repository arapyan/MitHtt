#ifndef DATAMC_HH
#define DATAMC_HH

#include <sstream>
#include "TFile.h"
#include "TMath.h"
#include "TRegexp.h"
#include "TTree.h"
#include "TBranch.h"
#include "TRandom1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "BtagSF.hh"

TRandom1 randm(0xDEADBEEF);
enum { kNo, kDown, kUp };                     // systematic variation

// Get higgs mass point from sample name
Int_t higgsmass(TString basename)
{
  cout << basename << endl;
  stringstream ss(basename(TRegexp("[0-9][0-9]*"),3).Data());
  Int_t mass;
  ss >> mass;
  //if(basename.Contains("-gf-")) assert(mass>85 && mass<1200);
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

//Bool_t isbtagged(mithep::TJet *jet, Int_t isdata, UInt_t btageff, UInt_t mistag);

Double_t embUnfoldWgt(Double_t pt1, Double_t eta1, Double_t pt2, Double_t eta2);
Double_t unskimmedEntries(TString skimname);


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
//----------------------------------------------------------------------------------------
/*
Bool_t isbtagged(mithep::TJet *jet, Int_t isdata, UInt_t btageff, UInt_t mistag)
{

  //          mistag                         scale factor
  // TCHEM  0.0175 \pm .0003 \pm .0038      1.21 \pm .02 \pm .17
  //          btag eff.                      scale factor
  // TCHEM  0.63 \pm 0.01                   0.93 \pm 0.02 \pm 0.07

  // new scale factors
  // TCHEM	btag eff: 0.96 \pm 0.04		mistag rate: 0.0286 \pm 0.0003		mistag scale factor: 1.20 \pm 0.14
  // CSVM	btag eff: 0.97 \pm 0.04		mistag rate: 0.0152 \pm 0.0002		mistag scale factor: 1.10 \pm 0.11

  Bool_t btagged;
  Double_t demoteProb=0; // ~probability to demote from tagged
  if(btageff==kNo)        demoteProb = fabs(1-0.97); //1-0.93;  // SF = 0.93 -> 0.07 = (prob to demote from tagged status)
  else if(btageff==kDown) demoteProb = fabs(1-0.97+0.04); //1-0.93+0.07;
  else if(btageff==kUp)   demoteProb = fabs(1-0.97-0.04); //1-0.93-0.07;
  Double_t promoteProb=0; // ~probability to promote to tagged
  if(mistag==kNo)         promoteProb = fabs(1.10-1)*0.0152/(1-0.0152); //(1.21-1)*0.0145/(1-0.0145);  // (1-SF)*mistag = (prob. to promote to tagged status)*(1-mistag)
  else if(mistag==kDown)  promoteProb = fabs(1.10-1+0.11)*0.0152/(1-0.0152);
  else if(mistag==kUp)    promoteProb = fabs(1.10-1-0.11)*0.0152/(1-0.0152);

  UInt_t jetflavor = 0;
                   
  if(isdata == 1) {
    if(jet->csv>0.679) btagged = kTRUE;
    else               btagged = kFALSE;
  } else { // MC
    //if(isdata == 0)jetflavor = abs(jet->mcFlavor);
    jetflavor = abs(jet->matchedId);
    if(jetflavor==5) {
      if(jet->csv>0.679) {
      if(randm.Uniform()>demoteProb) btagged = kTRUE;  // leave it tagged
      else                           btagged = kFALSE; // demote it
      } else                           btagged = kFALSE; // leave it untagged
    } else { // not bjet
      if(jet->csv>0.679)                   btagged = kTRUE;  // leave it tagged
      else if(randm.Uniform()<promoteProb) btagged = kTRUE;  // promote to tagged
      else                                 btagged = kFALSE; // leave it untagged
    }
  }

  return btagged;
} 
*/ 
//----------------------------------------------------------------------------------------
Double_t embUnfoldWgt(Double_t pt1, Double_t eta1, Double_t pt2, Double_t eta2)
{
  TFile *unfFile1   = TFile::Open("mutau.root"); assert(unfFile1->IsOpen());
  TH2F  *unfWeight1 = (TH2F*) unfFile1->FindObjectAny("hadunf");
  TH2F  *unfWeight2 = (TH2F*) unfFile1->FindObjectAny("lepunf");
  double weight1 = unfWeight1->GetBinContent(unfWeight1->GetXaxis()->FindBin(pt1),unfWeight1->GetYaxis()->FindBin(eta1));
  double weight2 = unfWeight2->GetBinContent(unfWeight2->GetXaxis()->FindBin(pt2),unfWeight2->GetYaxis()->FindBin(eta2));
  unfFile1->Close();
  return weight1*weight2;
}

//----------------------------------------------------------------

Double_t muIDscaleEmu(Double_t mupt, Double_t mueta, Int_t is2012)
{
  if((fabs(mueta) > 2.1) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  if(mupt > 20) {
    if(fabs(mueta) < 1.5)  return is2012 ? 1.006 : 0.992;
    else   return is2012 ? 1.014 : 0.994;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 1.5)   return is2012 ? 1.017 : 0.995;
    else               return is2012 ? 1.025 : 1.002;
  }
  else {
    if(fabs(mueta) < 1.5)   return is2012 ? 0.990 : 0.991;
    else                 return is2012 ? 1.030 : 1.036;
  }
}
//----------------------------------------------------------------------------------------    
Double_t eleIDscaleEmu(Double_t elept, Double_t eleeta, Int_t is2012)
{
  if((fabs(eleeta) > 2.3) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  else if(elept > 20) {
    if(fabs(eleeta) < 0.8) return is2012 ? 0.959 : 0.985 ;
    else if(fabs(eleeta) < 1.479) return is2012 ?  0.954 : 0.985;
    else         return is2012 ? 0.968 : 1.012;
  }
  else if(elept > 15) {
    if(fabs(eleeta) < 0.8) return is2012  ? 0.926 : 0.962;
    else if(fabs(eleeta) < 1.479) return  is2012 ? 0.853 : 0.962;
    else                     return  is2012 ? 0.838 : 1.148;
  }
  else {
    if(fabs(eleeta) < 0.8) return is2012 ? 0.840 : 1.040;
    else if(fabs(eleeta) < 1.479) return is2012 ? 0.837 : 1.040;
    else                     return is2012 ? 0.722 : 0.976;
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
  else if(elept > 30) {
    if(fabs(eleeta) < 0.8)        return is2012 ? 1.00 : 1.003;
    else if(fabs(eleeta) < 1.479) return is2012 ? 0.99 : 1.003;
    else                          return is2012 ? 0.97 : 1.008;
  } 
  else if(elept > 25 && is2012) {
    if(fabs(eleeta) < 0.8)        return 0.97;
    else if(fabs(eleeta) < 1.479) return 0.98;
    else                          return 1.01;
  } 
  else if(elept > 20) {
    if(fabs(eleeta) < 0.8)        return is2012 ? 0.98 : 1.001;
    else if(fabs(eleeta) < 1.479) return is2012 ? 0.96 : 1.001;
    else                          return is2012 ? 1.01 : 1.00;
  } 
  else if(elept > 15) {
    if(fabs(eleeta) < 0.8)        return is2012 ? 0.99 : 1.00;
    else if(fabs(eleeta) < 1.479) return is2012 ? 1.00 : 1.00;
    else                          return is2012 ? 1.07 : 1.05;
  } 
  else {
    if(fabs(eleeta) < 0.8)        return is2012 ? 0.99 : 0.98;
    else if(fabs(eleeta) < 1.479) return is2012 ? 0.82 : 0.98;
    else                          return is2012 ? 0.96 : 0.97;
  }
}
//----------------------------------------------------------------------------------------
Double_t muTrigScaleEmu(Double_t mupt, Double_t mueta, Int_t is2012)
{
 if((fabs(mueta) > 2.1) || (mupt < 10)) { cout << "muon kinematics out of range" << endl; assert(0); }
  else if(mupt > 30) {
    if(fabs(mueta) < 0.8)        return is2012 ? 1.07 : 0.992;
    else if(fabs(mueta) < 1.2) return is2012 ? 1.12 : 0.992;
    else                       return is2012 ? 1.12 : 1.06;
  } 
  else if(mupt > 25 && is2012) {
    if(fabs(mueta) < 0.8)        return 1.00;
    else if(fabs(mueta) < 1.2) return 0.106;
    else                          return 1.02;
  } 
  else if(mupt > 20) {
    if(fabs(mueta) < 0.8)        return is2012 ? 1.01 : 0.99;
    else if(fabs(mueta) < 1.2) return is2012 ? 0.98 : 0.99;
    else                          return is2012 ? 0.96 : 1.04;
  } 
  else if(mupt > 15) {
    if(fabs(mueta) < 0.8)        return is2012 ? 1.00 : 0.99;
    else if(fabs(mueta) < 1.2) return is2012 ? 1.04 : 0.99;
    else                          return is2012 ? 1.02 : 1.07;
  } 
  else {
    if(fabs(mueta) < 0.8)        return is2012 ? 1.00 : 1.01;
    else if(fabs(mueta) < 1.2) return is2012 ? 0.99 : 1.01;
    else                          return is2012 ? 0.98 : 1.03;
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
