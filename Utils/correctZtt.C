#include <vector>
#include <sstream>
#include <iostream>
#include <string>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom1.h"
#include "/scratch/dkralph/Recoil/NYStyle.h"

TFile *fZMCFileNVtx = new TFile("/scratch/dkralph/Recoil/ztt_select_recoil_novtx.root");
TTree *fZMCTreeNVtx = (TTree*) fZMCFileNVtx->FindObjectAny("Events");

TFile *fZMCFile = new TFile("/scratch/dkralph/Recoil/ztt_select_recoil.root");
TTree *fZMCTree = (TTree*) fZMCFile->FindObjectAny("Events");

TFile *fWWMCFile = new TFile("/scratch/dkralph/Recoil/ww_select_recoil.root");
TTree *fWWMCTree = (TTree*) fWWMCFile->FindObjectAny("Events");

TFile *fQCDMCFile = new TFile("/scratch/dkralph/Recoil/qcd_select_recoil.root");
TTree *fQCDMCTree = (TTree*) fQCDMCFile->FindObjectAny("Events");

TFile *fTTMCFile = new TFile("/scratch/dkralph/Recoil/ttbar_select_recoil.root");
TTree *fTTMCTree = (TTree*) fTTMCFile->FindObjectAny("Events");

TFile *fDataFile = new TFile("/scratch/dkralph/Recoil/data_select_recoil.root");
TTree *fZDataTree = (TTree*) fDataFile->FindObjectAny("Events"); 

double calculate(int iMet,double iEPt,double iEPhi,double iWPhi,double iU1,double iU2) { 
  double lMX = -iEPt*cos(iEPhi) - iU1*cos(iWPhi) + iU2*sin(iWPhi);
  double lMY = -iEPt*sin(iEPhi) - iU1*sin(iWPhi) - iU2*cos(iWPhi);
  if(iMet == 0) return sqrt(lMX*lMX + lMY*lMY);
  if(iMet == 1) {if(lMX > 0) {return atan(lMY/lMX);} return (fabs(lMY)/lMY)*3.14159265 + atan(lMY/lMX); } 
  if(iMet == 2) return lMX;
  if(iMet == 3) return lMY;
  return lMY;
}
double getCorError2(double iVal,TF1 *iFit) { 
  double lE = sqrt(iFit->GetParError(0))  + iVal*sqrt(iFit->GetParError(2));
  if(fabs(iFit->GetParError(4)) > 0) lE += iVal*iVal*sqrt(iFit->GetParError(4));
  return lE*lE;
}
double getError2(double iVal,TF1 *iFit) { 
  double lE2 = iFit->GetParError(0) + iVal*iFit->GetParError(1) + iVal*iVal*iFit->GetParError(2);
  if(fabs(iFit->GetParError(3)) > 0) lE2 += iVal*iVal*iVal*     iFit->GetParError(3);
  if(fabs(iFit->GetParError(4)) > 0) lE2 += iVal*iVal*iVal*iVal*iFit->GetParError(4);
  if(fabs(iFit->GetParError(5)) > 0 && iFit->GetParameter(3) == 0) lE2 += iVal*iVal*               iFit->GetParError(5);
  if(fabs(iFit->GetParError(5)) > 0 && iFit->GetParameter(3) != 0) lE2 += iVal*iVal*iVal*iVal*iVal*iFit->GetParError(5);
  if(fabs(iFit->GetParError(6)) > 0) lE2 += iVal*iVal*iVal*iVal*iVal*iVal*iFit->GetParError(6);
  return lE2;
}
double getError(double iVal,TF1 *iWFit,TF1 *iZDatFit,TF1 *iZMCFit,bool iRescale=true) {
  double lEW2  = getError2(iVal,iWFit);
  if(!iRescale) return sqrt(lEW2);
  double lEZD2 = getError2(iVal,iZDatFit);
  double lEZM2 = getError2(iVal,iZMCFit);
  double lZDat = iZDatFit->Eval(iVal);
  double lZMC  = iZMCFit->Eval(iVal);
  double lWMC  = iWFit->Eval(iVal);
  double lR    = lZDat/lZMC;
  double lER   = lR*lR/lZDat/lZDat*lEZD2 + lR*lR/lZMC/lZMC*lEZM2;
  double lVal  = lR*lR*lEW2 + lWMC*lWMC*lER;
  return sqrt(lVal);//+(iZMCFit->Eval(iVal)-iWFit->Eval(iVal);
}
const int readRecoil(std::vector<double> &iSumEt,
		     std::vector<TF1*> &iU1Fit,std::vector<TF1*> &iU1MRMSFit,std::vector<TF1*> &iU1RMS1Fit,std::vector<TF1*> &iU1RMS2Fit,std::vector<TF1*> &iU1Sig3Fit,
		     std::vector<TF1*> &iU2Fit,std::vector<TF1*> &iU2MRMSFit,std::vector<TF1*> &iU2RMS1Fit,std::vector<TF1*> &iU2RMS2Fit,std::vector<TF1*> &iU2Sig3Fit,
		     std::string iFName = "/scratch/dkralph/Recoil/recoilfit.root") { 
  TFile *lFile  = new TFile(iFName.c_str());
  TGraph *lGraph = (TGraph *) lFile->FindObjectAny("sumet");
  const int lNBins = lGraph->GetN();
  for(int i0 = 0; i0 < lNBins; i0++) {
    iSumEt.push_back(lGraph->GetY()[i0]);
    std::stringstream pSS1,pSS2,pSS3,pSS4,pSS5,pSS6,pSS7,pSS8,pSS9,pSS10;
    pSS1  << "u1Mean_"    << i0;  iU1Fit.push_back    ( (TF1*) lFile->FindObjectAny(pSS1.str().c_str())); //iU1Fit[i0]->SetDirectory(0);
    pSS2  << "u1MeanRMS_" << i0;  iU1MRMSFit.push_back( (TF1*) lFile->FindObjectAny(pSS2.str().c_str())); //iU1RMSFit[i0]->SetDirectory(0);
    pSS3  << "u1RMS1_"    << i0;  iU1RMS1Fit.push_back( (TF1*) lFile->FindObjectAny(pSS3.str().c_str())); //iU1RMSFit[i0]->SetDirectory(0);
    pSS4  << "u1RMS2_"    << i0;  iU1RMS2Fit.push_back( (TF1*) lFile->FindObjectAny(pSS4.str().c_str())); //iU1RMSFit[i0]->SetDirectory(0);
    pSS5  << "u1Sig3_"    << i0;  iU1Sig3Fit.push_back( (TF1*) lFile->FindObjectAny(pSS5.str().c_str())); //iU2RMSFit[i0]->SetDirectory(0);
    pSS6  << "u2Mean_"    << i0;  iU2Fit    .push_back( (TF1*) lFile->FindObjectAny(pSS6.str().c_str())); //iU2Fit[i0]->SetDirectory(0);
    pSS7  << "u2MeanRMS_" << i0;  iU2MRMSFit.push_back( (TF1*) lFile->FindObjectAny(pSS7.str().c_str())); //iU2RMSFit[i0]->SetDirectory(0);
    pSS8  << "u2RMS1_"    << i0;  iU2RMS1Fit.push_back( (TF1*) lFile->FindObjectAny(pSS8.str().c_str())); //iU2RMSFit[i0]->SetDirectory(0);
    pSS9  << "u2RMS2_"    << i0;  iU2RMS2Fit.push_back( (TF1*) lFile->FindObjectAny(pSS9.str().c_str())); //iU2RMSFit[i0]->SetDirectory(0);
    pSS10 << "u2Sig3_"    << i0;  iU2Sig3Fit.push_back( (TF1*) lFile->FindObjectAny(pSS10.str().c_str())); //iU2RMSFit[i0]->SetDirectory(0);
  }
  lFile->Close();
  iSumEt.push_back(1000);  
  return lNBins;
}
void metDistribution(double &iMet,double &iMPhi,double iGenPt,double iGenPhi,
		     double iLepPt,double iLepPhi,TRandom1 *iRand,
		     TF1 *iU1RWFit   ,TF1 *iU1RZDatFit  ,TF1 *iU1RZMCFit,
		     TF1 *iU1MSWFit  ,TF1 *iU1MSZDatFit ,TF1 *iU1MSZMCFit, 
		     TF1 *iU1S1WFit  ,TF1 *iU1S1ZDatFit ,TF1 *iU1S1ZMCFit,
		     TF1 *iU1S2WFit  ,TF1 *iU1S2ZDatFit ,TF1 *iU1S2ZMCFit,
		     TF1 *iU13SWFit  ,TF1 *iU13SZDatFit ,TF1 *iU13SZMCFit,
		     TF1 *iU2MSWFit  ,TF1 *iU2MSZDatFit ,TF1 *iU2MSZMCFit, 		      
		     TF1 *iU2S1WFit  ,TF1 *iU2S1ZDatFit ,TF1 *iU2S1ZMCFit, 
		     TF1 *iU2S2WFit  ,TF1 *iU2S2ZDatFit ,TF1 *iU2S2ZMCFit, 
		     TF1 *iU23SWFit  ,TF1 *iU23SZDatFit ,TF1 *iU23SZMCFit, 
		     int iFluc=0,int iMetType=0) {
  //Important constants re-scaling of sigma on left and mean wpt of W resbos on right
  double lRescale = sqrt((TMath::Pi())/2.); double lPtMean = 16.3; //==> tuned for W bosons
  double pU1      = iU1RWFit->Eval(iGenPt);
  double pU2      = 0; //Right guys are for cumulants => code deleted
  double pSigma1_1 = iU1MSWFit->Eval(iGenPt)*lRescale*iU1S1WFit->Eval(iGenPt);
  double pSigma1_2 = iU1MSWFit->Eval(iGenPt)*lRescale*iU1S2WFit->Eval(iGenPt);
  double pFrac1    = iU1MSWFit->Eval(iGenPt)*lRescale;
  double pCorr1    = iU13SWFit->Eval(iGenPt)/iU13SWFit->Eval(lPtMean);
  double pSigma2_1 = iU2MSWFit->Eval(iGenPt)*lRescale*iU2S1WFit->Eval(iGenPt);
  double pSigma2_2 = iU2MSWFit->Eval(iGenPt)*lRescale*iU2S2WFit->Eval(iGenPt);
  double pFrac2    = iU2MSWFit->Eval(iGenPt)*lRescale;
  double pCorr2    = iU23SWFit->Eval(iGenPt)/iU23SWFit->Eval(lPtMean);
  
  //Left is Cumulant right is standard
  double pZMSigma1_1 = iU1S1ZMCFit->Eval(iGenPt)*iU1MSZMCFit->Eval(iGenPt);
  double pZMSigma1_2 = iU1S2ZMCFit->Eval(iGenPt)*iU1MSZMCFit->Eval(iGenPt);
  double pZMFrac1    = iU1MSZMCFit->Eval(iGenPt);
  double pZMSigma2_1 = iU2S1ZMCFit->Eval(iGenPt)*iU2MSZMCFit->Eval(iGenPt);
  double pZMSigma2_2 = iU2S2ZMCFit->Eval(iGenPt)*iU2MSZMCFit->Eval(iGenPt);
  double pZMFrac2    = iU2MSZMCFit->Eval(iGenPt);
  double pZMCorr1    = iU13SZMCFit->Eval(iGenPt)/iU13SZMCFit->Eval(lPtMean);
  double pZMCorr2    = iU23SZMCFit->Eval(iGenPt)/iU23SZMCFit->Eval(lPtMean);
  
  double pZDSigma1_1 = iU1S1ZDatFit->Eval(iGenPt)*iU1MSZDatFit->Eval(iGenPt);
  double pZDSigma1_2 = iU1S2ZDatFit->Eval(iGenPt)*iU1MSZDatFit->Eval(iGenPt);
  double pZDFrac1    = iU1MSZDatFit->Eval(iGenPt);
  double pZDSigma2_1 = iU2S1ZDatFit->Eval(iGenPt)*iU2MSZDatFit->Eval(iGenPt);
  double pZDSigma2_2 = iU2S2ZDatFit->Eval(iGenPt)*iU2MSZDatFit->Eval(iGenPt);
  double pZDFrac2    = iU2MSZDatFit->Eval(iGenPt);
  double pZDCorr1    = iU13SZDatFit->Eval(iGenPt)/iU13SZDatFit->Eval(lPtMean);
  double pZDCorr2    = iU23SZDatFit->Eval(iGenPt)/iU23SZDatFit->Eval(lPtMean);

  pU1       =  pU1*    (iU1RZDatFit->Eval(iGenPt)/iU1RZMCFit->Eval(iGenPt));
  pFrac1    =  pZDFrac1/pZMFrac1 * pFrac1;
  pFrac2    =  pZDFrac2/pZMFrac2 * pFrac2;
  pSigma1_1 =  pSigma1_1*pZDSigma1_1/pZMSigma1_1;
  pSigma1_2 =  pSigma1_2*pZDSigma1_2/pZMSigma1_2;
  pSigma2_1 =  pSigma2_1*pZDSigma2_1/pZMSigma2_1;
  pSigma2_2 =  pSigma2_2*pZDSigma2_2/pZMSigma2_2;
  pCorr1    =  pCorr1   *pZDCorr1/pZMCorr1;
  pCorr2    =  pCorr2   *pZDCorr2/pZMCorr2;
  
  //Uncertainty propagation
  if(iFluc != 0) { 
    double lEUR1    = getError(iGenPt,iU1RWFit ,iU1RZDatFit ,iU1RZMCFit ,true);
    double lEUS1_1  = getError(0.,iU1S1WFit,iU1S1ZDatFit,iU1S1ZMCFit    ,true);
    double lEUS1_2  = getError(0.,iU1S2WFit,iU1S2ZDatFit,iU1S2ZMCFit    ,true);
    double lEU1Frac = getError(iGenPt,iU1MSWFit,iU1MSZDatFit,iU1MSZMCFit,true);
    double lEU1Corr = getError(iGenPt,iU13SWFit,iU13SZDatFit,iU13SZMCFit,true);
    double lEUS2_1  = getError(0.,iU2S1WFit,iU2S1ZDatFit,iU2S1ZMCFit    ,true);
    double lEUS2_2  = getError(0.,iU2S2WFit,iU2S2ZDatFit,iU2S2ZMCFit    ,true);
    double lEU2Frac = getError(iGenPt,iU2MSWFit,iU2MSZDatFit,iU2MSZMCFit,true);
    double lEU2Corr = getError(iGenPt,iU23SWFit,iU23SZDatFit,iU23SZMCFit,true);
    
    //Modify all the different parameters the choice of signs makes it maximal
    pU1       = pU1       - iFluc*lEUR1;              //Recoil
    pFrac1    = pFrac1    + iFluc*(lEU1Frac);        //Mean RMS 
    pSigma1_1 = pSigma1_1 - iFluc*lEUS1_1*pFrac1;    //Sigma 1 smalles sigma
    pSigma1_2 = pSigma1_2 + iFluc*lEUS1_2*pFrac1;    //Sigma 2 (Maximal when oppsite sigma 1)
    pCorr1    = pCorr1    + iFluc*(lEU1Corr);        //Correction from > 3 <sigma> to be consistent must be same sign as sigma 2
    pFrac2    = pFrac2    + iFluc*(lEU2Frac);        //Mean RMS for U2
    pSigma2_1 = pSigma2_1 - iFluc*lEUS2_1*pFrac2;    //Sigma 1 U2
    pSigma2_2 = pSigma2_2 + iFluc*(lEUS2_2)*pFrac2;
    pCorr2    = pCorr2    + iFluc*(lEU2Corr);
  }
  //Now calcualte recoil
  double pMS1 = pFrac1;
  double pMS2 = pFrac2;
  //Caculat the proper fraction
  pFrac1 = (pFrac1-pSigma1_2)/(pSigma1_1-pSigma1_2);
  pFrac2 = (pFrac2-pSigma2_2)/(pSigma2_1-pSigma2_2);
  //Apply the correction to sigma 2
  pSigma1_2*=pCorr1;
  pSigma2_2*=pCorr2;   
  //Constant Sigma correction
  if(iMetType == 0) {
    if(iFluc != 0) pSigma1_1 *= (1 + iFluc*0.08);  pSigma2_1 *= 1 + iFluc*0.08 ;
    //Systematic on Simga variation 0.1 for muons 0.04 for electrons
    pFrac1 = (pMS1-pSigma1_2)/(pSigma1_1-pSigma1_2); //Constant Sigma 
    pFrac2 = (pMS2-pSigma2_2)/(pSigma2_1-pSigma2_2);
  }
  //Constat Fraction correction
  if(iMetType == 1) {
    if(iFluc != 0 ) {pFrac2 = pFrac2 - iFluc*0.04; pFrac1 = pFrac1 - iFluc*0.04;}
    //Systematic on fraction variation 0.1 for mu 0.04 for e-
    pSigma1_1 = (pMS1-(1-pFrac1)*pSigma1_2)/pFrac1;//(pFrac1-pSigma1_2)/(pSigma1_1-pSigma1_2); //V1
    pSigma2_1 = (pMS2-(1-pFrac2)*pSigma2_2)/pFrac2;//(pFrac2-pSigma2_2)/(pSigma2_1-pSigma2_2);
  }
  //Now sample for the MET distribution
  double lVal0 = iRand->Uniform(0,1);
  double lVal1 = iRand->Uniform(0,1);
  pU1   = (lVal0 < pFrac1)*iRand->Gaus(pU1,pSigma1_1)+(lVal0 > pFrac1)*iRand->Gaus(pU1,pSigma1_2);
  pU2   = (lVal1 < pFrac2)*iRand->Gaus(pU2,pSigma2_1)+(lVal1 > pFrac2)*iRand->Gaus(pU2,pSigma2_2);
  iMet  = calculate(0,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  iMPhi = calculate(1,iLepPt,iLepPhi,iGenPhi,pU1,pU2);
  return;
}
void format(TH1F *iH) {
  iH->SetLineWidth(1);
  iH->SetFillStyle(3001);
}

void correctZtt(std::string iNameZDat = "/scratch/dkralph/Recoil/recoilfitZDat.root",std::string iNameZMC = "/scratch/dkralph/Recoil/recoilfitZMC.root") { 
  Prep();
  double lLumi = 190;
  TRandom1 *lRandom = new TRandom1(0xDEADBEEF);
  vector<double> lWSumEt;  vector<double> lWWeights; 
  vector<double> lZDSumEt; vector<double> lZDWeights; 
  vector<double> lZMSumEt; vector<double> lZMWeights; 
  vector<TF1*>   lWU1Fit;  vector<TF1*> lWU1RMSSMFit;  vector<TF1*> lWU1RMS1Fit;  vector<TF1*> lWU1RMS2Fit;  vector<TF1*> lWU13SigFit; 
  vector<TF1*>   lWU2Fit;  vector<TF1*> lWU2RMSSMFit;  vector<TF1*> lWU2RMS1Fit;  vector<TF1*> lWU2RMS2Fit;  vector<TF1*> lWU23SigFit; 
  vector<TF1*>   lZDU1Fit; vector<TF1*> lZDU1RMSSMFit; vector<TF1*> lZDU1RMS1Fit; vector<TF1*> lZDU1RMS2Fit; vector<TF1*> lZDU13SigFit; 
  vector<TF1*>   lZDU2Fit; vector<TF1*> lZDU2RMSSMFit; vector<TF1*> lZDU2RMS1Fit; vector<TF1*> lZDU2RMS2Fit; vector<TF1*> lZDU23SigFit;  
  vector<TF1*>   lZMU1Fit; vector<TF1*> lZMU1RMSSMFit; vector<TF1*> lZMU1RMS1Fit; vector<TF1*> lZMU1RMS2Fit; vector<TF1*> lZMU13SigFit; 
  vector<TF1*>   lZMU2Fit; vector<TF1*> lZMU2RMSSMFit; vector<TF1*> lZMU2RMS1Fit; vector<TF1*> lZMU2RMS2Fit; vector<TF1*> lZMU23SigFit; 
 
  const int lWNBins  =readRecoil(lWSumEt ,lWU1Fit ,lWU1RMSSMFit ,lWU1RMS1Fit ,lWU1RMS2Fit ,lWU13SigFit ,lWU2Fit ,lWU2RMSSMFit ,lWU2RMS1Fit ,lWU2RMS2Fit ,lWU23SigFit);
  const int lZDNBins =readRecoil(lZDSumEt,lZDU1Fit,lZDU1RMSSMFit,lZDU1RMS1Fit,lZDU1RMS2Fit,lZDU13SigFit,lZDU2Fit,lZDU2RMSSMFit,lZDU2RMS1Fit,lZDU2RMS2Fit,lZDU23SigFit,iNameZDat);
  const int lZMNBins =readRecoil(lZMSumEt,lZMU1Fit,lZMU1RMSSMFit,lZMU1RMS1Fit,lZMU1RMS2Fit,lZMU13SigFit,lZMU2Fit,lZMU2RMSSMFit,lZMU2RMS1Fit,lZMU2RMS2Fit,lZMU23SigFit,iNameZMC);
  
  double lGenPt  = 0; double lGenPhi  = 0; 
  double lLepPt1 = 0; double lLepPhi1 = 0;  double lLepEta1 = 0;
  double lLepPt2 = 0; double lLepPhi2 = 0;  double lLepEta2 = 0;
  double lMet    = 0; double lMPhi    = 0;  double lWeight  = 0;
  fZMCTreeNVtx->SetBranchAddress("vpt"   ,&lGenPt);
  fZMCTreeNVtx->SetBranchAddress("vphi"  ,&lGenPhi);
  fZMCTreeNVtx->SetBranchAddress("lpt1"  ,&lLepPt1);
  fZMCTreeNVtx->SetBranchAddress("lphi1" ,&lLepPhi1);
  fZMCTreeNVtx->SetBranchAddress("leta1" ,&lLepEta1);
  fZMCTreeNVtx->SetBranchAddress("lpt2"  ,&lLepPt2);
  fZMCTreeNVtx->SetBranchAddress("lphi2" ,&lLepPhi2);
  fZMCTreeNVtx->SetBranchAddress("leta2" ,&lLepEta2);
  fZMCTreeNVtx->SetBranchAddress("met"   ,&lMet);
  fZMCTreeNVtx->SetBranchAddress("metphi",&lMPhi);
  fZMCTreeNVtx->SetBranchAddress("wgt"   ,&lWeight);
  int iNBins = 60; double iMin = -100; double iMax = 200;
  TH1F* lDMet  = new TH1F("dMET","dMET"  ,iNBins,iMin,iMax); lDMet->Sumw2();  lDMet->SetMarkerStyle(20);
  TH1F* lMMet0 = new TH1F("mMET0","mMET0",iNBins,iMin,iMax); lMMet0->Sumw2(); lMMet0->SetLineColor(kRed);
  TH1F* lMMet1 = new TH1F("mMET1","mMET1",iNBins,iMin,iMax); lMMet1->Sumw2(); lMMet1->SetLineColor(kBlue);
  TH1F* lMMet2 = new TH1F("mMET2","mMET2",iNBins,iMin,iMax); lMMet2->Sumw2(); lMMet2->SetLineColor(kGreen);
  TH1F* lTTMet = new TH1F("ttMET","ttMET",iNBins,iMin,iMax); lTTMet->Sumw2(); lTTMet->SetFillColor(kGreen);
  TH1F* lWWMet = new TH1F("wwMET","wwMET",iNBins,iMin,iMax); lWWMet->Sumw2(); lWWMet->SetFillColor(kViolet);
  TH1F* lQCMet = new TH1F("qcMET","qcMET",iNBins,iMin,iMax); lQCMet->Sumw2(); lQCMet->SetFillColor(kOrange);
  format(lTTMet); format(lWWMet); format(lMMet0); format(lMMet1); format(lMMet2); format(lDMet); format(lQCMet);
  double lCut = 0; double lTot = 0;
  for(int i0 = 0; i0 < fZMCTreeNVtx->GetEntries(); i0++) {
    //double pGenPt  = lRandom->Gaus(30,10);
    //double pGenPhi = lRandom->Uniform(-1.*TMath::Pi(),TMath::Pi());

    //double pLepPt  = lRandom->Gaus(50,10);
    //double pLepPhi = lRandom->Uniform(-1.*TMath::Pi(),TMath::Pi());

    //double pMX  = -pGenPt*cos(pGenPhi) - pLepPt*cos(pLepPhi);
    //double pMY  = -pGenPt*sin(pGenPhi) - pLepPt*sin(pLepPhi);
    //
    fZMCTreeNVtx->GetEntry(i0);
    
    TLorentzVector m,e,met;
    m.SetPtEtaPhiM(lLepPt1,lLepEta1,lLepPhi1,0.105658369);
    e.SetPtEtaPhiM(lLepPt2,lLepEta2,lLepPhi2,0.000511);
    met.SetPtEtaPhiM(lMet,0,lMPhi,0);
    TLorentzVector boson = m + e + met;
    Double_t phi1=m.Phi(), phi2=e.Phi();
    double lDPhi = fabs(m.Phi()-e.Phi()); if(lDPhi > 2*TMath::Pi()-lDPhi) lDPhi = 2*TMath::Pi()-lDPhi;
    double projPhi = 0;
    if(phi1>phi2)
      if(phi1-phi2 < TMath::Pi())
	projPhi = phi1 - lDPhi/2;
      else
	projPhi = phi1 + lDPhi/2;
    else
      if(phi2-phi1 < TMath::Pi())
	projPhi = phi2 - lDPhi/2;
      else
	projPhi = phi2 + lDPhi/2;
    
    Double_t projX  =  cos(projPhi);
    Double_t projY  =  sin(projPhi);
    //proj    = boson.Px()*projX  +   boson.Py()*projY;
    //projVis = dilep.Px()*projX  +   dilep.Py()*projY;

    double pMet  = lMet;
    double pMPhi = lMPhi; 
    double pMProj  = pMet*cos(pMPhi)*projX + pMet*sin(pMPhi)*projY;
    double pVProj  = (m+e).Px()*projX + (m+e).Py()*projY;

    //lMMet1->Fill(pMProj,lWeight*lLumi);
    printf("====> Met Before ==> %10.2f%10.2f\n",pMet,pMPhi);
    metDistribution(pMet,pMPhi,lGenPt,lGenPhi,(m+e).Pt(),(m+e).Phi(),lRandom,
		    lWU1Fit     [0],lZDU1Fit     [0],lZMU1Fit     [0],
		    lWU1RMSSMFit[0],lZDU1RMSSMFit[0],lZMU1RMSSMFit[0],
		    lWU1RMS1Fit [0],lZDU1RMS1Fit [0],lZMU1RMS1Fit [0],
		    lWU1RMS2Fit [0],lZDU1RMS2Fit [0],lZMU1RMS2Fit [0],
		    lWU13SigFit [0],lZDU13SigFit [0],lZMU13SigFit [0],
		    lWU2RMSSMFit[0],lZDU2RMSSMFit[0],lZMU2RMSSMFit[0],
		    lWU2RMS1Fit [0],lZDU2RMS1Fit [0],lZMU2RMS1Fit [0],
		    lWU2RMS2Fit [0],lZDU2RMS2Fit [0],lZMU2RMS2Fit [0],
		    lWU23SigFit [0],lZDU23SigFit [0],lZMU23SigFit [0]);
		    //,int iFluc=0,int iMetType=1) {
    printf("   => Met  After ==> %10.2f%10.2f\n",pMet,pMPhi);
    if(i0>1000) break;

    pMProj  = pMet*cos(pMPhi)*projX + pMet*sin(pMPhi)*projY;
    lMMet1->Fill(0.85*pVProj-pMProj,lWeight*lLumi);
    //lMMet1->Fill(pMet,lWeight*lLumi);
    //lMMet1->Fill(pMProj,lWeight*lLumi);
    lTot++; if(0.85*pVProj-pMProj < 25.) lCut++;
  }
  cout << "====> Yield ===> " << lCut/lTot << endl;
  //std::string iVar = "projMet"; std::stringstream iCut; iCut <<  "wgt*" << lLumi;
  std::string iVar = "0.85*projVis-projMet"; std::stringstream iCut; iCut <<  "wgt*" << lLumi;
  //std::string iVar = "met"; std::stringstream iCut; iCut <<  "wgt*" << lLumi;
  TCanvas *lC0 = new TCanvas("c0","c0",800,600); lC0->cd();
  fZDataTree  ->Draw((iVar+">>dMET").c_str());
  fZMCTree    ->Draw((iVar+">>mMET0").c_str(),iCut.str().c_str());
  fZMCTreeNVtx->Draw((iVar+">>mMET2").c_str(),iCut.str().c_str());
  fTTMCTree   ->Draw((iVar+">>ttMET").c_str(),iCut.str().c_str());
  fWWMCTree   ->Draw((iVar+">>wwMET").c_str(),iCut.str().c_str());
  fQCDMCTree  ->Draw((iVar+">>qcMET").c_str(),iCut.str().c_str());
  
  lWWMet->Add(lQCMet);
  lTTMet->Add(lWWMet);
  lMMet0->Add(lTTMet);
  lMMet1->Add(lTTMet);

  TLegend *lL = new TLegend(0.2,0.6,0.5,0.9); lL->SetFillColor(0); lL->SetBorderSize(0);
  //lL->AddEntry(lDMet,"data"         ,"lp");
  //lL->AddEntry(lMMet2,"Z#rightarrow#tau#tau " ,"lp");
  lL->AddEntry(lMMet0,"Z#rightarrow#tau#tau vertex weighted" ,"lp");
  lL->AddEntry(lMMet1,"Z#rightarrow#tau#tau recoil corrected","lp");
  lL->AddEntry(lTTMet,"t#bar{t}"    ,"f");
  lL->AddEntry(lWWMet,"EWK"         ,"f");
  lL->AddEntry(lQCMet,"QCD"         ,"f");
  
  lMMet0->SetLineStyle(kDashed);
  //lMMet0->Scale(lDMet->Integral()/lMMet0->Integral()); 
  //lMMet1->Scale(lDMet->Integral()/lMMet1->Integral()); 
  //lMMet2->Scale(lDMet->Integral()/lMMet2->Integral()); 
  //lDMet->GetYaxis()->SetRangeUser(0,0.05);
  std::stringstream lYaxis; lYaxis << "Events/" << ((lDMet->GetXaxis()->GetXmax()-lDMet->GetXaxis()->GetXmin())/lDMet->GetNbinsX()) << "GeV";


  lDMet->GetXaxis()->SetTitle("0.85*p_{#zeta}^{vis} - #slash{p}_{#zeta} [GeV]");//P_{vis}-0.85*P_{#zeta} (GeV)");
  //lDMet->GetXaxis()->SetTitle("#slash{E_{T}} [GeV]");//P_{vis}-0.85*P_{#zeta} (GeV)");
  //lDMet->GetXaxis()->SetTitle("#slash{p}_{#zeta} [GeV]");//P_{vis}-0.85*P_{#zeta} (GeV)");
  lDMet->GetYaxis()->SetTitle(lYaxis.str().c_str());
  lDMet->GetYaxis()->SetRangeUser(0,250);
  lDMet->Draw("EP");

  //lMMet0->GetXaxis()->SetTitle("#slash{p}_{#zeta} [GeV]");//P_{vis}-0.85*P_{#zeta} (GeV)");
  //lMMet0->GetYaxis()->SetTitle(lYaxis.str().c_str());
  //lMMet0->GetYaxis()->SetRangeUser(0,250);

  lTTMet->Draw("hist sames");
  lWWMet->Draw("hist sames");
  lQCMet->Draw("hist sames");
  lMMet0->Draw("hist sames");
  lDMet->Draw("EP sames");
  lMMet1->Draw("hist sames");
  //lMMet2->Draw("hist sames");
  lL->Draw();

  
}
/*
    double pLepPx = lLepPt1*cos(lLepPhi1) + lLepPt2*cos(lLepPhi2);
    double pLepPy = lLepPt1*sin(lLepPhi1) + lLepPt2*sin(lLepPhi2);
    double pLepPt = sqrt(pLepPx*pLepPx + pLepPy*pLepPy);
    double pLepPhi = 0; if(pLepPx > 0) {pLepPhi = atan(pLepPy/pLepPx);} else {pLepPhi = (fabs(pLepPy)/pLepPx)*3.14159265 + atan(pLepPy/pLepPx); }
    double pDPhi   = lLepPhi1-lLepPhi2; double pDCorrPhi = fabs(pDPhi); 
    if(fabs(pDPhi) > 2*TMath::Pi()-fabs(pDPhi)) pDCorrPhi = 2*TMath::Pi()-pDPhi;
    if(fabs(pDPhi) < TMath::Pi()) pDPhi = lLepPhi1 - pDCorrPhi/2.*pDPhi/fabs(pDPhi);
    if(fabs(pDPhi) > TMath::Pi()) pDPhi = lLepPhi1 + pDCorrPhi/2.*pDPhi/fabs(pDPhi);
*/
