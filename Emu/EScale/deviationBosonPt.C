#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <sstream>
#include "TLegend.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TRandom1.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooNumConvPdf.h"
#include "RooPolynomial.h"
#include "RooPlot.h"
#include "RooVoigtian.h"
#include "RooVoigtianShape.h"
#include "RooMCStudy.h"
#include "TNtuple.h"
#include "TRandom1.h"
#include "../Efficiency/Efficiency/interface/IdTools1.hh"
#include "../NYStyle/test/NYStyle.h"

using namespace RooFit;

//   1  mpar0       -4.02982e-01   4.65808e-02   2.33606e-02  -8.06838e-02
//   2  mpar1       -9.33897e-03   1.01796e-03   8.85938e-05  -1.86780e-03
//   3  mpar2       -1.83072e-01   1.55013e-01   7.85689e-02  -3.66226e-02
//   4  mpar3        9.24141e-03   9.33590e-04   7.25787e-04   1.84828e-03
//   5  mpar4       -2.31002e-02   1.45368e-01   7.89709e-02  -4.62006e-03

//   1  mpar0        6.11585e-02   3.61356e-02   2.35810e-02   1.22320e-02
//   2  mpar1       -9.36898e-03   9.20670e-04   4.63118e-04  -1.87380e-03
//   3  mpar2       -1.06297e-03   3.81337e-04   2.98644e-04  -2.12594e-04
//   4  mpar3        1.01272e-02   1.06699e-03   5.00514e-04   2.02544e-03
//   5  mpar4       -2.14348e-04   4.15047e-04   3.25035e-04  -4.28696e-05
//   6  mpar5        1.54704e-03   3.37067e-04   1.81255e-04   3.09408e-04
//   7  mpar6       -1.91883e-03   3.92682e-04   1.91277e-04  -3.83766e-04

//   1  mpar0        5.33434e-02   3.54348e-02   7.19299e-04   1.06689e-02
//   2  mpar1       -9.31989e-03   6.33885e-04   1.39045e-05  -1.86398e-03
//   3  mpar2       -8.43053e-04   4.13967e-04   9.30918e-06  -1.68611e-04
//   4  mpar3        1.00752e-02   6.45572e-04   1.40227e-05   2.01504e-03
//   5  mpar4       -1.13772e-04   4.27215e-04   9.30368e-06  -2.27545e-05
//   6  mpar5        1.40453e-03   2.51902e-04   5.65303e-06   2.80906e-04
//   7  mpar6       -1.98450e-03   2.66175e-04   5.64994e-06  -3.96900e-04
//   8  s1par0       1.30114e+00   1.62395e-02   6.59430e-04  -8.32732e-01
//   9  s1par1       5.09124e-01   1.87257e-02   4.27394e-04   1.02002e-01
//  10  s1par2       3.38297e-01   1.25673e-02   2.74296e-04   6.77111e-02
//  11  s1par3      -3.86663e-01   7.49179e-03   1.68251e-04  -7.74099e-02
//  12  s1par4      -1.86877e-01   4.01580e-03   9.77382e-05  -3.73842e-02


//   1  mpar0        5.30310e-02   3.88876e-02   2.58437e-02   1.06064e-02
//   2  mpar1       -9.24456e-03   1.33329e-03   5.30351e-04  -1.84891e-03
//   3  mpar2       -6.64266e-04   4.59133e-04   3.48520e-04  -1.32853e-04
//   4  mpar3        9.98830e-03   9.24758e-04   6.42868e-07   1.99766e-03
//   5  mpar4        1.43595e-04   4.39903e-04   4.24064e-05   2.87189e-05
//   6  mpar5        1.55671e-03   4.97718e-04   2.06109e-04   3.11342e-04
//   7  mpar6       -1.83657e-03   3.73729e-04   2.06383e-04  -3.67315e-04
//   8  s1par0       1.29895e+00   1.97207e-02   2.34508e-02  -8.33383e-01
//   9  s1par1       3.51482e-03   4.92370e-02   1.62733e-02   7.02964e-04
//  10  s1par2      -1.19386e-02   2.07310e-02   1.09138e-02  -2.38773e-03
//  11  s1par3      -1.43569e-02   1.94919e-02   6.37825e-03  -2.87138e-03
//  12  s1par4       3.67112e-02   6.05739e-03   3.69476e-03   7.34232e-03
//MC
//   1  mpar0       -2.79958e-01   1.42387e-02   2.25611e-02  -5.60209e-02
//   2  mpar1       -1.09243e-02   5.64497e-03   7.00850e-04  -2.18486e-03
//   3  mpar2       -1.36237e+00   4.53507e-01   6.57210e-02  -2.75963e-01
//   4  mpar3        9.36308e-03   5.39088e-03   7.14352e-04   1.87262e-03
//   5  mpar4       -1.27907e+00   5.93203e-01   7.51255e-02  -2.58690e-01

//   1  mpar0       -2.78531e-01   2.92677e-02   2.65133e-02  -5.57350e-02
//   2  mpar1        1.41613e-03   7.75707e-04   5.41385e-04   2.83226e-04
//   3  mpar2        1.24622e-04   3.43319e-04   3.60777e-04   2.49245e-05
//   4  mpar3       -1.91784e-03   8.15336e-04   3.57637e-06  -3.83567e-04
//   5  mpar4        4.15290e-04   3.36938e-04   3.63808e-04   8.30580e-05
//   6  mpar5       -1.89737e-03   2.94334e-04   2.18543e-04  -3.79475e-04
//   7  mpar6        2.03956e-03   3.03364e-04   2.83126e-07   4.07912e-04
//   8  s1par0       1.40042e+00   1.70463e-02   1.46760e-04  -8.03682e-01
//   9  s1par1       8.54439e-02   2.67302e-02   1.68671e-02   1.70896e-02
//  10  s1par2      -9.10097e-02   1.47824e-02   1.11748e-02  -1.82029e-02
//  11  s1par3      -4.72824e-02   1.06480e-02   6.87880e-03  -9.45662e-03
//  12  s1par4       5.66641e-02   4.22902e-03   3.90164e-03   1.13331e-02

double correct(double iM,double iPhi1,double iPhi2,double iEta1,double iEta2) { 
  return fabs(iM);
  //38
  //double   lM=(fabs(iM)+0.31)/sqrt((1+0.00473*cos(iPhi1+1.832))*(1+0.00504*cos(iPhi2-1.248))); ///-0.0168 -0.681  1.54 1.66;
  //double   lM=(fabs(iM)+0.31)/sqrt((1-0.01092*sin(iPhi1-1.36))*(1+0.00936*sin(iPhi2-1.28))); ///-0.0168 -0.681  1.54 1.66;
  //Data=>
  double lM=(fabs(iM)+0.37)/sqrt((1-0.00548*sin(iPhi1-0.5399))*(1+0.00937*sin(iPhi2+0.8056))); ///-0.0168 -0.681  1.54 1.66
  //Data=>
  //return lM/sqrt((1-0.009244*iEta1-0.0006644*iEta1*iEta1+0.00156*iEta1*iEta1*iEta1)*(1+0.009988*iEta2-0.000144*iEta2*iEta2-0.00184*iEta2*iEta2*iEta2)); 
  return lM;//MC=>sqrt((1+0.00141*iEta1-0.000124*iEta1*iEta1-0.001897*iEta1*iEta1*iEta1)*(1-0.001917*iEta2-0.0004153*iEta2*iEta2+0.00204*iEta2*iEta2*iEta2)); 
  //return lM;
  //Data double   
  //lM=(fabs(iM)+0.37)/sqrt((1-0.00484*sin(iPhi1-0.612))*(1-0.00987*sin(iPhi2+3.85))); ///-0.0168 -0.681  1.54 1.66
  //return lM/sqrt((1-0.00898*iEta1-0.00069*iEta1*iEta1+0.00134*iEta1*iEta1*iEta1)*(1+0.01009*iEta2+0.000139*iEta2*iEta2-0.00188*iEta2*iEta2*iEta2)); 
  //39
  //double   lM=(iM-0.2)/sqrt((1-0.00151*iPhi1-0.00153*iPhi1*iPhi1)*(1+0.001872*iPhi2-0.0003297*iPhi2*iPhi2));
  //return lM/sqrt((1+0.011370*iEta2+0.001549*iEta2*iEta2)*(1-0.009077*iEta1+0.00069194*iEta1*iEta1));//Quadratic //
  //return lM/sqrt((1+0.01*iEta2-0.0000001549*iEta2*iEta2)*(1-0.01255*iEta1-0.000000655*iEta1*iEta1));//Quadratic //
  //XXreturn lM/sqrt((1-0.00894*iEta1-0.000095344*iEta1*iEta1+0.007725*iEta1*iEta1*iEta1)*(1-0.0010271*iEta2-0.001348*iEta2*iEta2+0.001953*iEta2*iEta2*iEta2)); 
    //sqrt((1+0.0010271*iEta1+0.001348*iEta1*iEta1-0.001953*iEta1*iEta1*iEta1)*(1-0.0010271*iEta2-0.001348*iEta2*iEta2+0.001953*iEta2*iEta2*iEta2));
		 //(1-0.00894*iEta1-0.000095344*iEta1*iEta1+0.007725*iEta1*iEta1*iEta1)*
		 //(1-0.0010271*iEta2-0.001348*iEta2*iEta2+0.001953*iEta2*iEta2*iEta2)); //cubic
  //return (iM/sqrt((1-0.00207*iPhi1-0.00153*iPhi1*iPhi1)*(1+0.0038486*iPhi2-0.00221*iPhi2*iPhi2)));
  //double lM=(iM/sqrt((1-0.00207*iPhi1-0.00153*iPhi1*iPhi1)*(1+0.0038486*iPhi2-0.00221*iPhi2*iPhi2)));
  //double lM=(iM/sqrt((1-0.00179*iPhi1-0.00163*iPhi1*iPhi1)*(1+0.00399*iPhi2-0.00229*iPhi2*iPhi2)));
  //return lM/sqrt((1-0.00475*iEta1-0.00175186*iEta1*iEta1)*(1+0.00343*iEta2+0.0002176*iEta2*iEta2));
  //return (iM/sqrt((1+0.00449*iPhi1-0.0022072*iPhi1*iPhi1)*(1-0.00149349*iPhi2-0.0016897*iPhi2*iPhi2)));
}
void plotValue(std::string iName,RooMCStudy *iMC,RooRealVar &iVar) { 
  TCanvas *lC0 = new TCanvas((iName+" Value").c_str(),(iName+" Value").c_str(),600, 400); lC0->cd();
  RooPlot* lF0 = iMC->plotParam(iVar,Bins(40)) ;
  lF0->Draw();  

  TCanvas *lC1 = new TCanvas((iName+" Error").c_str(),(iName+" Error").c_str(),600, 400); lC1->cd();
  RooPlot* lF1 = iMC->plotError(iVar,Bins(40)) ;
  lF1->Draw();  

  TCanvas *lC2 = new TCanvas((iName+" Pull").c_str(),(iName+" Pull").c_str(),600, 400); lC2->cd();
  RooPlot* lF2 = iMC->plotPull(iVar,Bins(40),FitGauss(kTRUE)) ;
  lF2->Draw();  
}

void fitMass1(int iCharge,double iEtaMin,double iEtaMax,double iPhiMin,double iPhiMax,double &iRes,double &iShift,double &iResErr,double &iShiftErr) { 
  Prep();
  RooRealVar l1Sigma("sigma1","sigma1",1.7 ,0.,15.);  //l1Sigma.setConstant(kTRUE);
  RooRealVar l2Sigma("sigma2","sigma2",1.15,0.,15.);  l2Sigma.setConstant(kTRUE);
  RooRealVar l3Sigma("sigma3","sigma3",2.8,0.,35.);   l3Sigma.setConstant(kTRUE);
  RooRealVar lN     ("n"     ,"n"     ,1.5,-15,15.);  lN.setConstant(kTRUE);
  RooRealVar lExp   ("exp"   ,"exp"   ,-0.003,-15,15.); //lExp.setConstant(kTRUE);
  RooRealVar lR1Mean("mean","mean",90.8,60,150);   //lR1Mean.setConstant(kTRUE);
  RooRealVar lRXVar("XVar","mass(GeV/c^{2})",90,60,120); //lRXVar.setBins(500);
  RooRealVar lRYield("Yield","mass(GeV/c^{2})",10000,0,20000);
  RooVoigtianShape     lGAdd("Add","Add",lRXVar,lR1Mean,l1Sigma,l2Sigma,lN,l3Sigma,true);
  //RooExtendPdf lGEx("Ex","Ex",lGAdd,lRYield);
  //RooAddPdf      lGAdd("XAdd","XAdd",lGAdd1,lExpF,lCFrac);
  RooDataSet *lData = new RooDataSet("crap"  ,"crap",RooArgSet(lRXVar)); 
  RooDataSet *lXData = new RooDataSet("xcrap","xcrap",RooArgSet(lRXVar)); 
  RooDataSet *lYData = new RooDataSet("ycrap","ycrap",RooArgSet(lRXVar)); 
  TFile *lFile = new TFile("../Efficiency/Data/mTPNT8_v1.root");//G39TP.root");
  //TFile *lFile = new TFile("../Efficiency/Data/ZTP8_v2.root");//G39TP.root");
  TTree *lTree = (TTree*) lFile->FindObjectAny("WNtupleIdEffNT");
  float lMt   = 0;  lTree->SetBranchAddress("mt"    ,&lMt);
  int lCharge = 0;  lTree->SetBranchAddress("charge",&lCharge);
  float lEta  = 0;  lTree->SetBranchAddress("eta"   ,&lEta);
  float lPhi  = 0;  lTree->SetBranchAddress("phi"   ,&lPhi);
  float lPt   = 0;  lTree->SetBranchAddress("pt"    ,&lPt);
  float lOPt  = 0;  lTree->SetBranchAddress("jetpt" ,&lOPt);
  float lOPhi  = 0; lTree->SetBranchAddress("jetphi",&lOPhi);
  float lOEta  = 0; lTree->SetBranchAddress("jeteta",&lOEta);
  float lTrkIso   = 0; lTree->SetBranchAddress("trkiso",&lTrkIso);
  float lEcalIso  = 0; lTree->SetBranchAddress("ecaliso",&lEcalIso);
  float lHcalIso  = 0; lTree->SetBranchAddress("hcaliso",&lHcalIso);
  float lChi2     = 0; lTree->SetBranchAddress("chi2"   ,&lChi2);
  int   lNHit     = 0; lTree->SetBranchAddress("nhit"   ,&lNHit);
  Muon  lMuon;         lTree->SetBranchAddress("muon"   ,&lMuon);
  for(int i0 = 0; i0 < lTree->GetEntries();i0++) { 
    if(i0 > 200000) break;
    lTree->GetEntry(i0);
    if(fabs(lPt) < 25 || fabs(lOPt) < 25) continue; 
    if(lMt < 60) continue;
    if(lChi2 > 10) continue;
    if(lNHit < 10) continue;
    if(lMuon.NSeg < 2) continue;
    if(lMuon.NPixel == 0) continue;
    if(lMuon.NValid == 0) continue;
    if(lMuon.Type   != 3) continue;
    if((lTrkIso+lEcalIso+lHcalIso)/lPt > 0.15) continue;
    if(lCharge > 0 && iCharge < 0) continue;
    if(lCharge < 0 && iCharge > 0) continue;
    //if(lEta < iEtaMin || lEta > iEtaMax) continue;
    if(lPhi < iPhiMin || lPhi > iPhiMax) continue;
    //if(lEta < iPhiMin || lEta > iPhiMax) continue;
    if(lCharge > 0) lRXVar.setVal(correct(lMt,lPhi,lOPhi,lEta,lOEta));
    if(lCharge < 0) lRXVar.setVal(correct(lMt,lOPhi,lPhi,lOEta,lEta));
    //cout << "====> " << lMt << " ---> " << lRXVar.getVal() << " --- " << endl;
    //lXData->add(RooArgSet(lRXVar));
    //lRXVar.setVal(fabs(lMt));
    lYData->add(RooArgSet(lRXVar));
  }
  /*
    RooHistPdf  *lMYPdf  = new RooHistPdf ("MY","MY",RooArgList(lRXVar),*lYData->binnedClone()); 
    TFile *lQFile = new TFile("../Efficiency/Data/mTPNT8_v1.root");//G39TP.root");
    TTree *lQTree = (TTree*) lQFile->FindObjectAny("WNtupleIdEffNT");
    lQTree->SetBranchAddress("mt"    ,&lMt);
    lQTree->SetBranchAddress("charge",&lCharge);
    lQTree->SetBranchAddress("eta"   ,&lEta);
    lQTree->SetBranchAddress("phi"   ,&lPhi);
    lQTree->SetBranchAddress("pt"    ,&lPt);
    lQTree->SetBranchAddress("jetpt" ,&lOPt);
    lQTree->SetBranchAddress("jetphi",&lOPhi);
    lQTree->SetBranchAddress("jeteta",&lOEta);
    for(int i0 = 0; i0 < lQTree->GetEntries();i0++) { 
    lQTree->GetEntry(i0);
    //if(lOPt*lPt < 0) continue;
    if(lMt < 60) continue;
    // if(lCharge > 0 && iCharge < 0) continue;
    //if(lCharge < 0 && iCharge > 0) continue;
    //if(lEta < iEtaMin || lEta > iEtaMax) continue;
    //if(lPhi < iPhiMin || lPhi > iPhiMax) continue;
    if(lCharge > 0) lRXVar.setVal(correct(lMt,lPhi,lOPhi,lEta,lOEta));
    if(lCharge < 0) lRXVar.setVal(correct(lMt,lOPhi,lPhi,lOEta,lEta));
    lXData->add(RooArgSet(lRXVar));
    lRXVar.setVal(lMt);
    lData->add(RooArgSet(lRXVar));
  }
  RooHistPdf  *lMPdf   = new RooHistPdf ("MH","MH",RooArgList(lRXVar),*lXData->binnedClone()); 
  */
  lGAdd.fitTo(*lYData,Strategy(1),Minos());
  if(lR1Mean.getError() < 0.01) lGAdd.fitTo(*lYData,Strategy(2),Minos());
  lRXVar.setBins(60);
  iShift = 1./(lR1Mean.getVal()/91.2); iShiftErr = lR1Mean.getError()*91.2/lR1Mean.getVal()/lR1Mean.getVal();
  iRes   = l1Sigma.getVal(); iResErr   = l1Sigma.getError();
  return;
  
  TH1F *lH0 = new TH1F("A","A",1,-5,5); lH0->SetMarkerStyle(21); 
  TH1F *lH1 = new TH1F("B","B",1,-5,5); lH1->SetLineColor(kRed);
  TH1F *lH2 = new TH1F("C","C",1,-5,5); lH2->SetLineColor(kBlue);
  TLegend *lL = new TLegend(0.5,0.5,0.8,0.8); lL->SetFillColor(0); lL->SetBorderSize(0);
  lL->AddEntry(lH0,"data"         ,"lp");
  lL->AddEntry(lH1,"MC-corrected" ,"l");
  lL->AddEntry(lH2,"MC"           ,"l");

  RooPlot *lFrame1 = lRXVar.frame(RooFit::Title("XXX")) ;
  lYData->plotOn(lFrame1);
  lGAdd.plotOn(lFrame1);
  //lXData->plotOn(lFrame1,MarkerColor(kRed));
  //lMPdf->plotOn(lFrame1,LineColor(kRed));
  //lMYPdf->plotOn(lFrame1,LineColor(kBlue));
  TCanvas *iC =new TCanvas("A","A",800,600);
  iC->cd(); lFrame1->Draw();
  //lL->Draw();
  iC->SaveAs("Crap.png");
  //cin.get();
  
}

void fitMass2(int iCharge,double iEtaMin,double iEtaMax,double iPhiMin,double iPhiMax,double &iRes,double &iShift,double &iResErr,double &iShiftErr) { 
  Prep();
  RooRealVar    lXVar  ("XVar","mass(GeV/c^{2})",60,60,120); lXVar.setBins(1000);
  RooRealVar    lSPar  ("SPar","SPar", 1.,0., 2.);
  RooFormulaVar lXShift("uparshift","@0*@1",RooArgList(lXVar,lSPar));
  TFile *lMCFile = new TFile("../Efficiency/Data/ZTP8_v2.root");
  TTree *lMCTree = (TTree*) lMCFile->FindObjectAny("WNtupleIdEffNT"); 
  TH1F *lMass = new TH1F("M","M",100,60,120); 
  int   lCharge = 0;   lMCTree->SetBranchAddress("charge",&lCharge);
  float lEta    = 0;   lMCTree->SetBranchAddress("eta"   ,&lEta);
  float lPhi    = 0;   lMCTree->SetBranchAddress("phi"   ,&lPhi);
  float lMt     = 0;   lMCTree->SetBranchAddress("mt"    ,&lMt);
  float lPt     = 0;   lMCTree->SetBranchAddress("pt"    ,&lPt);
  float lOPt    = 0;   lMCTree->SetBranchAddress("jetpt" ,&lOPt);
  float lOEta   = 0;   lMCTree->SetBranchAddress("jeteta",&lOEta);
  float lOPhi   = 0;   lMCTree->SetBranchAddress("jetphi",&lOPhi);
  float lTrkIso   = 0; lMCTree->SetBranchAddress("trkiso",&lTrkIso);
  float lEcalIso  = 0; lMCTree->SetBranchAddress("ecaliso",&lEcalIso);
  float lHcalIso  = 0; lMCTree->SetBranchAddress("hcaliso",&lHcalIso);
  float lChi2     = 0; lMCTree->SetBranchAddress("chi2"   ,&lChi2);
  unsigned int   lNHit     = 0; lMCTree->SetBranchAddress("nhit"   ,&lNHit);
  Muon  lMuon;         lMCTree->SetBranchAddress("muon"   ,&lMuon);
  double lVPt   = 0; double lEt = 0; double lPx = 0; double lPy = 0;
  for(int i0 = 0; i0 < lMCTree->GetEntries(); i0++) {
    lMCTree->GetEntry(i0);
    if(lMt < 60)          continue;
    if(lChi2 > 10)        continue;
    if(lNHit < 10)        continue;
    if(lMuon.NSeg < 2)    continue;
    if(lMuon.NPixel == 0) continue;
    if(lMuon.NValid == 0) continue;
    if(lMuon.Type   != 3) continue;
    if((lTrkIso+lEcalIso+lHcalIso)/lPt > 0.15) continue;
    //if(lCharge > 0 && iCharge < 0) continue;
    //if(lCharge < 0 && iCharge > 0) continue;
    //if(lEta < iPhiMin || lEta > iPhiMax) continue;
    //if(fabs(lEta) < 1.2) continue;
    //if(lEta < iEtaMin || lEta > iEtaMax) continue;
    //if(lPhi*lOPhi > 0 || lEta*lOEta > 0) continue;
    if(lPhi < iPhiMin || lPhi > iPhiMax) continue;
    lEt = lOPt + lPt; lPx = fabs(lPt)*cos(lPhi) + fabs(lOPt)*cos(lOPhi); lPy = fabs(lPt)*sin(lPhi) + fabs(lOPt)*sin(lOPhi);
    lVPt = sqrt(lPx*lPx + lPy*lPy);
    //if(lVPt < iPhiMin || lVPt > iPhiMax) continue;
    lMass->Fill(fabs(lMt));     //lMass->Fill(lPt);
  }
  RooDataHist *lMHist = new RooDataHist("M" ,"M" ,RooArgSet(lXVar),lMass);
  RooHistPdf  *lMPdf  = new RooHistPdf ("MH","MH",lXShift,lXVar,*lMHist,5); 
  RooRealVar l1Sigma("sigma1","sigma1",0.2,0.,15.);  //l1Sigma.setConstant(kTRUE);
  RooRealVar lR0Mean("xmean","xmean",0,-10,10);    lR0Mean.setConstant(kTRUE);
  RooRealVar lExp   ("exp"   ,"exp"   ,-0.006,-15,15.); //lExp.setConstant(kTRUE);
  RooRealVar lFrac  ("frac","frac"    ,0.9,0.,1);
  RooGaussian   lGaus1("gaus1","gaus1",lXVar,lR0Mean,l1Sigma);
  RooExponential lExpF("Exp","Exp"  ,lXVar,lExp);
  RooFFTConvPdf  lConv("Conv","Conv",lXVar,*lMPdf,lGaus1); //lConv.setBufferStrategy(RooFFTConvPdf::Flat);
  RooAddPdf      lGAdd("Add","Add"  ,lConv,lExpF,lFrac);
  RooDataSet *lData = new RooDataSet("crap","crap",RooArgSet(lXVar)); 
  TFile *lFile = new TFile("../Efficiency/Data/mTPNT8_v1.root");
  TTree *lTree = (TTree*) lFile->FindObjectAny("WNtupleIdEffNT");
  lTree->SetBranchAddress("charge",&lCharge);
  lTree->SetBranchAddress("eta"   ,&lEta);
  lTree->SetBranchAddress("phi"   ,&lPhi);
  lTree->SetBranchAddress("mt"    ,&lMt);
  lTree->SetBranchAddress("pt"    ,&lPt);
  lTree->SetBranchAddress("jetpt" ,&lOPt);
  lTree->SetBranchAddress("jeteta",&lOEta);
  lTree->SetBranchAddress("jetphi",&lOPhi);
  lTree->SetBranchAddress("trkiso",&lTrkIso);
  lTree->SetBranchAddress("ecaliso",&lEcalIso);
  lTree->SetBranchAddress("hcaliso",&lHcalIso);
  lTree->SetBranchAddress("chi2"   ,&lChi2);
  lTree->SetBranchAddress("nhit"   ,&lNHit);
  lTree->SetBranchAddress("muon"   ,&lMuon);
  for(int i0 = 0; i0 < lTree->GetEntries();i0++) { 
    lTree->GetEntry(i0);
    if(lMt < 60)                   continue;
    if(lCharge > 0 && iCharge < 0) continue;
    if(lCharge < 0 && iCharge > 0) continue;
    if(lChi2 > 10) continue;
    if(lNHit < 10) continue;
    if(lMuon.NSeg < 2) continue;
    if(lMuon.NPixel == 0) continue;
    if(lMuon.NValid == 0) continue;
    if(lMuon.Type   != 3) continue;
    if((lTrkIso+lEcalIso+lHcalIso)/lPt > 0.15) continue;
    //if(fabs(lEta) < 1.2) continue;
    //if(lPhi*lOPhi > 0 || lEta*lOEta > 0) continue;
    //if(lEta < iPhiMin || lEta > iPhiMax) continue;
    if(lPhi < iPhiMin || lPhi > iPhiMax) continue;
    lEt = lOPt + lPt; lPx = fabs(lPt)*cos(lPhi) + fabs(lOPt)*cos(lOPhi); lPy = fabs(lPt)*sin(lPhi) + fabs(lOPt)*sin(lOPhi);
    lVPt = sqrt(lPx*lPx + lPy*lPy);
    //if(lVPt < iPhiMin || lVPt > iPhiMax) continue;
    lXVar.setVal(fabs(lMt));
    if(lCharge > 0) lXVar.setVal(correct(lMt,lPhi,lOPhi,lEta,lOEta));
    if(lCharge < 0) lXVar.setVal(correct(lMt,lOPhi,lPhi,lOEta,lEta));
    lData->add(RooArgSet(lXVar));
  }
  /*
  TFile *lQFile = new TFile("");
  TTree *lQTree = (TTree*) lQFile->FindObjectAny("WNtupleIdEffNT");
  lQTree->SetBranchAddress("mt"    ,&lMt);
  lQTree->SetBranchAddress("charge",&lCharge);
  lQTree->SetBranchAddress("eta"   ,&lEta);
  lQTree->SetBranchAddress("phi"   ,&lPhi);
  for(int i0 = 0; i0 < lQTree->GetEntries();i0++) { 
    lQTree->GetEntry(i0);
    if(lMt < 60) continue;
    if(lCharge > 0 && iCharge < 0) continue;
    if(lCharge < 0 && iCharge > 0) continue;
    if(lEta < iEtaMin || lEta > iEtaMax) continue;
    if(lPhi < iPhiMin || lPhi > iPhiMax) continue;
    lXVar.setVal(lMt);
    //lData->add(RooArgSet(lXVar));
  }
  */
  lConv.fitTo(*lData,Strategy(1));
  if(l1Sigma.getError() > 1) lConv.fitTo(*lData,Strategy(2));
  lXVar.setBins(60);
  RooPlot *lFrame1 = lXVar.frame(RooFit::Title("XXX")) ;
  lData->plotOn(lFrame1);
  lConv.plotOn(lFrame1);
  //lGAdd.plotOn(lFrame1,Components("Exp"),LineColor(kRed));
  TCanvas *iC =new TCanvas("A","A",800,600);
  iC->cd(); lFrame1->Draw();
  iC->SaveAs("Crap.png");
  //cin.get();
  iRes = l1Sigma.getVal(); iResErr   = l1Sigma.getError();
  iShift = lSPar.getVal(); iShiftErr = lSPar.getError();
}
void deviationBosonPt(int iNBin = 10) { 
  double dX = 2*3.14159265;//2*2.4
  double lInterval = dX/iNBin;
  double lPhi0 = -3.14159265;//-2.4;//3.14159265;
  double lX1Val[iNBin]; double lY1Val[iNBin]; double lX1EVal[iNBin]; double lY1EVal[iNBin];
  double lRY1Val[iNBin]; double lRY1EVal[iNBin];
  for(int i0 = 0; i0 < iNBin; i0++) { 
    lPhi0+=lInterval;
    lX1Val[i0]  = lPhi0-lInterval/2.;
    lX1EVal[i0] = lInterval/2.;
    //                                         resol.       shift      res. err.    shift err
    fitMass1(1,-2.5,2.5,lPhi0-lInterval,lPhi0,lRY1Val[i0],lY1Val[i0],lRY1EVal[i0],lY1EVal[i0]);
  }
  TGraphErrors *lGraph1 = new TGraphErrors(iNBin,lX1Val,lY1Val,lX1EVal,lY1EVal); lGraph1->SetLineColor(kRed); lGraph1->SetMarkerColor(kRed); lGraph1->SetMarkerStyle(20);
  TGraphErrors *lGraphR1 = new TGraphErrors(iNBin,lX1Val,lRY1Val,lX1EVal,lRY1EVal); lGraphR1->SetLineColor(kRed); lGraphR1->SetMarkerColor(kRed); lGraphR1->SetMarkerStyle(20);
  lGraph1->SetLineStyle(kDashed);     lGraphR1->SetLineStyle(kDashed);
  lGraph1->GetYaxis()->SetRangeUser(0.975,1.025); lGraphR1->GetYaxis()->SetRangeUser(0,3.0);
  TLegend *LL0 = new TLegend(0.35,0.2,0.8,0.5); LL0->SetBorderSize(0); LL0->SetFillColor(0);
  TLegend *LL1 = new TLegend(0.35,0.2,0.8,0.5); LL1->SetBorderSize(0); LL1->SetFillColor(0);
  LL0->AddEntry(lGraph1,"+ Muon","alp");
  LL1->AddEntry(lGraph1,"+ Muon","alp");

  lGraph1->GetXaxis()->SetTitle("#phi"); lGraph1->GetYaxis()->SetTitle("Energy Scale");
  TCanvas *lC0= new TCanvas("Res"  ,"Res",800,600);
  lGraph1->Draw("ap");
  LL0->Draw();

  lGraphR1->GetXaxis()->SetTitle("#phi"); lGraphR1->GetYaxis()->SetTitle("Additional Smearing (GeV/c^{2})");
  TCanvas *lC1= new TCanvas("Shift","Shift",800,600);
  lGraphR1->Draw("ap");
  lGraphR2->Draw("p");
  LL1->Draw();
  
  lGraph1 ->SetName("ScalePlus");  lGraph1->SetTitle("ScalePlus");
  lGraphR1->SetName("ResPlus");    lGraphR1->SetTitle("ResPlus");

  TFile *lFile = new TFile("EnergyScale.root","RECREATE");
  lGraph1->Write();
  lGraphR1->Write();
}    

void deviationBosonPtOld() { 
  TFile *lFile = new TFile("ZTPMC.root");
  TTree *lTree = (TTree*) lFile->FindObjectAny("WNtupleIdEffNT");
  //lTree->Print();
  //lTree->Draw("mt");
  std::string lZPt = "sqrt((abs(pt)*cos(phi)+abs(jetpt)*cos(jetphi))*(abs(pt)*cos(phi)+abs(jetpt)*cos(jetphi)) + (abs(pt)*sin(phi)+abs(jetpt)*sin(jetphi))*(abs(pt)*sin(phi)+abs(jetpt)*sin(jetphi)))";
  TProfile *lPProf = new TProfile("A","A",8,0,40); lPProf->SetMarkerColor(kBlue); lPProf->SetMarkerStyle(21); lPProf->SetLineColor(kBlue);
  TProfile *lMProf = new TProfile("B","B",8,0,40); lMProf->SetMarkerColor(kRed);  lMProf->SetMarkerStyle(21); lMProf->SetLineColor(kRed);
  lTree->Draw(std::string("mt:"+lZPt+">>A").c_str(),std::string("mt > 80 && charge > 0 ").c_str());
  lTree->Draw(std::string("mt:"+lZPt+">>B").c_str(),std::string("mt > 80 && charge < 0 ").c_str());
  //lTree->Draw(std::string("mt>>A").c_str(),"mt > 80 && charge > 0");
  //lTree->Draw(std::string("mt>>B").c_str(),"mt > 80 && charge < 0");
  
  TLegend *lL = new TLegend(0.6,0.6,0.9,0.9);
  lL->SetFillColor(0);
  lL->AddEntry(lPProf,"Plus","lp");
  lL->AddEntry(lMProf,"Minus","lp");
  lPProf->Draw();
  lMProf->Draw("same");
  lL->Draw();
}
//   1  s1par0       1.46010e+00   3.77879e-02   2.68182e-04  -7.86634e-01
//   2  s1par1       1.65505e-01   2.52491e-01   2.21493e-02   3.31070e-02
//   3  s1par2       3.70992e-02   5.43336e-02   1.42729e-04   7.41991e-03
//   4  s1par3      -2.90325e-01   1.36410e-01   8.18156e-05  -5.80976e-02
//   5  s1par4      -9.40729e-02   1.34435e-02   4.42113e-05  -1.88157e-02
//   6  s2par0       1.35149e+00   3.82580e-02   4.68436e-02  -8.17887e-01
//   7  s2par1      -9.97165e-02   1.59343e-01   1.81203e-04  -1.99446e-02
//   8  s2par2       2.46287e-02   2.49160e-02   2.06352e-02   4.92576e-03
//   9  s2par3       2.72721e-01   3.33808e-02   6.84386e-03   5.45713e-02
//  10  s2par4      -9.87116e-02   6.62868e-03   6.14061e-03  -1.97436e-02

void testConditional() { 
  RooRealVar lRXVar("XVar","XVar",90,60,120); lRXVar.setBins(2000);
  RooRealVar lSPar0("spar0","spar0",1.,0,10);  //lSPar0.setConstant(kTRUE);
  RooRealVar lS1Par0("s1par0","s1par0",1.0   ,0,10);  //lSPar0.setConstant(kTRUE);
  RooRealVar lS1Par1("s1par1","s1par1",0.02 ,-5,5); //lSPar1.setConstant(kTRUE);
  RooRealVar lS1Par2("s1par2","s1par2",-0.05,-5,5); //lSPar2.setConstant(kTRUE);
  RooRealVar lS1Par3("s1par3","s1par3", 0.0 ,-5,5);  //lSPar3.setConstant(kTRUE);
  RooRealVar lS1Par4("s1par4","s1par4",0.06,-5,5);   //lSPar4.setConstant(kTRUE);
  
  RooRealVar lS2Par0("s2par0","s2par0",1.26979,0,10);  //lSPar0.setConstant(kTRUE);
  RooRealVar lS2Par1("s2par1","s2par1",0.0127 ,-5,5); //lSPar1.setConstant(kTRUE);
  RooRealVar lS2Par2("s2par2","s2par2",-0.057,-5,5); //lSPar2.setConstant(kTRUE);
  RooRealVar lS2Par3("s2par3","s2par3",-0.19 ,-5,5);   //lSPar3.setConstant(kTRUE);
  RooRealVar lS2Par4("s2par4","s2par4",0.06  ,-5,5);    //SPar4.setConstant(kTRUE);
  
  RooRealVar lMPar0("mpar0","mpar0",0.,-5,5);
  RooRealVar lMPar1("mpar1","mpar1",0.,-5,5);
  RooRealVar lMPar2("mpar2","mpar2",0.,-5,5); //lMPar2.setConstant(kTRUE);
  RooRealVar lMPar3("mpar3","mpar3",0.,-5,5);
  RooRealVar lMPar4("mpar4","mpar4",0.,-5,5);// lMPar4.setConstant(kTRUE);
  RooRealVar lMPar5("mpar5","mpar5",0.,-5,5);
  RooRealVar lMPar6("mpar6","mpar6",0.,-5,5);
  RooRealVar lRPhi ("phi","phi"  ,0,-3.14,3.14);
  RooRealVar lRPhi2("phi2","phi2"  ,0,-3.14,3.14); 
  RooFormulaVar l1SigmaEta1("sigma1eta1","(@0+@1*@5+@2*@5*@5+@3*@5*@5*@5+@4*@5*@5*@5*@5)",RooArgList(lS1Par0,lS1Par1,lS1Par2,lS1Par3,lS1Par4,lRPhi));//,lSPar3,lSPar4,lRPhi2));    
  RooFormulaVar l1SigmaEta2("sigma1eta2","(@0+@1*@5+@2*@5*@5+@3*@5*@5*@5+@4*@5*@5*@5*@5)",RooArgList(lS1Par0,lS1Par1,lS1Par2,lS1Par3,lS1Par4,lRPhi2));//,lSPar3,lSPar4,lRPhi2));  
  //RooFormulaVar l1Sigma     ("sigma1"   ,"sqrt(@0*@0+@1*@1)"                             ,RooArgList(l1SigmaEta1,l1SigmaEta2));
  RooRealVar l1Sigma("sigma1","sigma1",1.47,0.,3.);    l1Sigma.setConstant(kTRUE);
  RooRealVar l2Sigma("sigma2","sigma2",1.1 ,0.,3.);    l2Sigma.setConstant(kTRUE);
  RooRealVar l3Sigma("sigma3","sigma3",2.8,0.,15.);    l3Sigma.setConstant(kTRUE);
  RooRealVar lN     ("n"     ,"n"     ,1.5,-15,15.);   lN.setConstant(kTRUE);
  //RooRealVar lR1Mean("mean","mean",91.2,60,250);   lR1Mean.setConstant(kTRUE);
  //RooFormulaVar lR1Mean("mean","91.2*sqrt((1+@1*sin(@3+@2))*(1+@1*sin(@6+@2)))+@0",RooArgList(lMPar0,lMPar1,lMPar2,lRPhi,lMPar3,lMPar4,lRPhi2));  //l1Sigma.setConstant(kTRUE);
  RooFormulaVar lR1Mean("mean","91.2*sqrt((1+@1*sin(@3+@2))*(1+@4*sin(@6+@5)))+@0",RooArgList(lMPar0,lMPar1,lMPar2,lRPhi,lMPar3,lMPar4,lRPhi2));  //l1Sigma.setConstant(kTRUE);
  //RooFormulaVar lR1Mean("mean","91.2*sqrt((1+@1*@3+@2*@3*@3+@7*@3*@3*@3)*(1+@4*@6+@5*@6*@6+@8*@6*@6*@6))+@0",RooArgList(lMPar0,lMPar1,lMPar2,lRPhi,lMPar3,lMPar4,lRPhi2,lMPar5,lMPar6));  //l1Sigma.setConstant(kTRUE);
  RooVoigtianShape     lConv("XXX","XXX",lRXVar,lR1Mean,l1Sigma,l2Sigma,lN,l3Sigma,false);
  RooDataSet     *lData = new RooDataSet("crap","crap",RooArgSet(lRXVar,lRPhi,lRPhi2)); 

  TFile *lFile = new TFile("../Efficiency/Data/mTPNT8_v1.root");
  //TFile *lFile = new TFile("../Efficiency/Data/ZTP8_v2.root");
  TTree *lTree = (TTree*) lFile->FindObjectAny("WNtupleIdEffNT");
  float lMt   = 0; lTree->SetBranchAddress("mt"    ,&lMt);
  int lCharge = 0; lTree->SetBranchAddress("charge",&lCharge);
  float lEta  = 0; lTree->SetBranchAddress("eta"   ,&lEta);
  float lPhi  = 0; lTree->SetBranchAddress("phi"   ,&lPhi);
  float lPt   = 0; lTree->SetBranchAddress("pt"    ,&lPt);
  float lOPt  = 0; lTree->SetBranchAddress("jetpt" ,&lOPt);
  float lOPhi = 0; lTree->SetBranchAddress("jetphi",&lOPhi);
  float lOEta = 0; lTree->SetBranchAddress("jeteta",&lOEta);
  float lTrkIso   = 0; lTree->SetBranchAddress("trkiso",&lTrkIso);
  float lEcalIso  = 0; lTree->SetBranchAddress("ecaliso",&lEcalIso);
  float lHcalIso  = 0; lTree->SetBranchAddress("hcaliso",&lHcalIso);
  float lChi2     = 0; lTree->SetBranchAddress("chi2"   ,&lChi2);
  int   lNHit     = 0; lTree->SetBranchAddress("nhit"   ,&lNHit);
  Muon  lMuon;         lTree->SetBranchAddress("muon"   ,&lMuon);
  for(int i0 = 0; i0 < lTree->GetEntries();i0++) { 
    if(i0 > 50000) continue;
    lTree->GetEntry(i0);
    if(lMt < 60) continue;
    if(lChi2 > 10) continue;
    if(lNHit < 10) continue;
    if(lMuon.NSeg < 2) continue;
    if(lMuon.NPixel == 0) continue;
    if(lMuon.NValid == 0) continue;
    if(lMuon.Type   != 3) continue;
    if((lTrkIso+lEcalIso+lHcalIso)/lPt > 0.15) continue;
    if(lMt < 60 || fabs(lEta) > 2.1 || fabs(lOEta) > 2.1 || fabs(lPt) < 25 || fabs(lOPt) < 25)  continue;
    lRXVar.setVal(lMt);

    //if(lCharge > 0) lRXVar.setVal(correct(lMt,lPhi,lOPhi,0,0));//lEta ,lOEta));
    //if(lCharge < 0) lRXVar.setVal(correct(lMt,lOPhi,lPhi,0,0));//lOEta,lEta));
    //lRPhi.setVal(lEta);  lRPhi2.setVal(lOEta);
    //if(lCharge < 0) {lRPhi.setVal(lOEta); lRPhi2.setVal(lEta);}
    lRPhi.setVal(lPhi);  lRPhi2.setVal(lOPhi);
    if(lCharge < 0) {lRPhi.setVal(lOPhi); lRPhi2.setVal(lPhi);}
    lData->add(RooArgSet(lRXVar,lRPhi,lRPhi2));//,lRPhi2));
  }
  lConv.fitTo(*lData,ConditionalObservables(RooArgSet(lRXVar,lRPhi,lRPhi2)),Strategy(1),Minos());
  //lConv.fitTo(*lData,Strategy(1));
  lRXVar.setBins(50);
  RooPlot *lFrame1 = lRXVar.frame(RooFit::Title("XXX")) ;
  lData->plotOn(lFrame1);
  lConv.plotOn(lFrame1);
  TCanvas *iC =new TCanvas("A","A",800,600);
  iC->cd(); lFrame1->Draw();
  iC->SaveAs("Crap.png");
}
/*
  RooRealVar lExp   ("exp"   ,"exp"   ,-0.006,-15,15.); //lExp.setConstant(kTRUE);
  RooRealVar lR0Mean("xmean","xmean",0,-10,10);    lR0Mean.setConstant(kTRUE);
  RooRealVar lFrac  ("frac","frac",0.46,0.,1);  lFrac.setConstant(kTRUE);
  RooRealVar lConst ("c","c"        ,1.,0.,2.); lConst.setConstant(kTRUE);
  RooRealVar lCFrac ("cfrac","cfrac",0.9,0.,1);
  RooPolynomial  lPoly("p","p",lRXVar,RooArgList(lConst),0) ;
  RooExponential lExpF("Exp","Exp",lRXVar,lExp);
  RooAddPdf      lAdd ("ExpAdd","ExpAdd",lExpF,lPoly,lCFrac);
  RooBreitWigner lBW("BW","BW",lRXVar,lR1Mean,l3Sigma);
  RooCBShape     lCB("CB","CB",lRXVar,lR0Mean,l1Sigma,l2Sigma,lN);
*/
//void makePdf() { 
//  TFile *lMCFile = new TFile("ZTPMC.root");
//  TTree *lMCTree = (TTree*) lMCFile->FindObjectAny("WNtupleIdEffNT"); 
//  TH1F *lMass = new TH1F("M","M",100,60,160); 
//  int   lCharge = 0; lMCTree->SetBranchAddress("charge",&lCharge);
//  float lEta    = 0; lMCTree->SetBranchAddress("eta"   ,&lEta);
//  float lPhi    = 0; lMCTree->SetBranchAddress("phi"   ,&lPhi);
//  float lMt     = 0; lMCTree->SetBranchAddress("mt"    ,&lMt);
//  float lPt     = 0; lMCTree->SetBranchAddress("pt"    ,&lPt);
//  float lOPt    = 0; lMCTree->SetBranchAddress("jetpt" ,&lOPt);
//  float lOEta   = 0; lMCTree->SetBranchAddress("jeteta",&lOEta);
//  float lOPhi   = 0; lMCTree->SetBranchAddress("jetphi",&lOPhi);
//  double lVPt   = 0; double lEt = 0; double lPx = 0; double lPy = 0;
//  for(int i0 = 0; i0 < lMCTree->GetEntries(); i0++) {
//    lMCTree->GetEntry(i0);
//    if(lMt < 60)                   continue;
//    if(lCharge > 0 && iCharge < 0) continue;
//    if(lCharge < 0 && iCharge > 0) continue;
//    if(lEta < iPhiMin || lEta > iPhiMax) continue;
//    //if(fabs(lEta) < 1.2) continue;
//    //if(lEta < iEtaMin || lEta > iEtaMax) continue;
//    //if(lPhi*lOPhi > 0 || lEta*lOEta > 0) continue;
//    //if(lPhi < iPhiMin || lPhi > iPhiMax) continue;
//    lEt = lOPt + lPt; lPx = fabs(lPt)*cos(lPhi) + fabs(lOPt)*cos(lOPhi); lPy = fabs(lPt)*sin(lPhi) + fabs(lOPt)*sin(lOPhi);
//    lVPt = sqrt(lPx*lPx + lPy*lPy);
//    //if(lVPt < iPhiMin || lVPt > iPhiMax) continue;
//    lMass->Fill(fabs(lMt));     //lMass->Fill(lPt);
//  }
//  RooDataHist *lMHist = new RooDataHist("M" ,"M" ,RooArgSet(lXVar),lMass);
//  RooHistPdf  *lMPdf  = new RooHistPdf ("MH","MH",lXShift,lXVar,*lMHist,0); 
//}
//void fitBins() { 
//  RooRealVar    lXVar  ("XVar","XVar",100,60,160); lXVar.setBins(1000);
//  RooRealVar    lSPar  ("SPar","SPar", 0.,-5,  5);
//  for(int i0 = 0; i0 < lNBins; i0++) { 
//    for(int i1 = 0; i0 <= i1; i1++) {
//      std::stringstream lSName;
//      
//      RooFormulaVar lXShift("uparshift","@0*@1",RooArgList(lXVar,lSPar));
//    }}
//
//  RooRealVar l1Sigma("sigma1","sigma1",0.2,0.,15.);  //l1Sigma.setConstant(kTRUE);
//  RooRealVar lR0Mean("xmean","xmean",0,-10,10);    lR0Mean.setConstant(kTRUE);
//  RooRealVar lExp   ("exp"   ,"exp"   ,-0.006,-15,15.); //lExp.setConstant(kTRUE);
//  RooRealVar lFrac  ("frac","frac"    ,0.9,0.,1);
//  RooGaussian   lGaus1("gaus1","gaus1",lXVar,lR0Mean,l1Sigma);
//  RooExponential lExpF("Exp","Exp"  ,lXVar,lExp);
//  RooFFTConvPdf  lConv("Conv","Conv",lXVar,*lMPdf,lGaus1,9); lConv.setBufferStrategy(RooFFTConvPdf::Flat);
//  RooAddPdf      lGAdd("Add","Add"  ,*lMPdf,lExpF,lFrac);
//  RooDataSet *lData = new RooDataSet("crap","crap",RooArgSet(lXVar)); 
//  TFile *lFile = new TFile("ZTP.root");
//  TTree *lTree = (TTree*) lFile->FindObjectAny("WNtupleIdEffNT");
//  lTree->SetBranchAddress("charge",&lCharge);
//  lTree->SetBranchAddress("eta"   ,&lEta);
//  lTree->SetBranchAddress("phi"   ,&lPhi);
//  lTree->SetBranchAddress("mt"    ,&lMt);
//  lTree->SetBranchAddress("pt"    ,&lPt);
//  lTree->SetBranchAddress("jetpt" ,&lOPt);
//  lTree->SetBranchAddress("jeteta",&lOEta);
//  lTree->SetBranchAddress("jetphi",&lOPhi);
//  for(int i0 = 0; i0 < lTree->GetEntries();i0++) { 
//    lTree->GetEntry(i0);
//    if(lMt < 60)                   continue;
//    if(lCharge > 0 && iCharge < 0) continue;
//    if(lCharge < 0 && iCharge > 0) continue;
//    if(lEta < iPhiMin || lEta > iPhiMax) continue;
//    //if(fabs(lEta) < 1.2) continue;
//    //if(lPhi*lOPhi > 0 || lEta*lOEta > 0) continue;
//    //if(lPhi < iPhiMin || lPhi > iPhiMax) continue;
//    lEt = lOPt + lPt; lPx = fabs(lPt)*cos(lPhi) + fabs(lOPt)*cos(lOPhi); lPy = fabs(lPt)*sin(lPhi) + fabs(lOPt)*sin(lOPhi);
//    lVPt = sqrt(lPx*lPx + lPy*lPy);
//    //if(lVPt < iPhiMin || lVPt > iPhiMax) continue;
//    lXVar.setVal(fabs(lMt));
//    lData->add(RooArgSet(lXVar));
//  }
//  lConv.fitTo(*lData,Strategy(2));
//  lXVar.setBins(20);
//  RooPlot *lFrame1 = lXVar.frame(RooFit::Title("XXX")) ;
//  lData->plotOn(lFrame1);
//  lConv.plotOn(lFrame1);
//  //lGAdd.plotOn(lFrame1,Components("Exp"),LineColor(kRed));
//  TCanvas *iC =new TCanvas("A","A",800,600);
//  iC->cd(); lFrame1->Draw();
//  iC->SaveAs("Crap.png");
//  //cin.get();
//  iRes = l1Sigma.getVal(); iResErr   = l1Sigma.getError();
//  iShift = lSPar.getVal(); iShiftErr = lSPar.getError();
//}
TH1F* getMass(int iCharge, double iPhiMin, double iPhiMax, double iEtaMin, double iEtaMax) { 
  TFile *lMCFile = new TFile("ZTPMC.root");
  TTree *lMCTree = (TTree*) lMCFile->FindObjectAny("WNtupleIdEffNT"); 
  TH1F *lMass = new TH1F("M","M",100,60,160); 
  int   lCharge = 0; lMCTree->SetBranchAddress("charge",&lCharge);
  float lEta    = 0; lMCTree->SetBranchAddress("eta"   ,&lEta);
  float lPhi    = 0; lMCTree->SetBranchAddress("phi"   ,&lPhi);
  float lMt     = 0; lMCTree->SetBranchAddress("mt"    ,&lMt);
  float lPt     = 0; lMCTree->SetBranchAddress("pt"    ,&lPt);
  float lOPt    = 0; lMCTree->SetBranchAddress("jetpt" ,&lOPt);
  float lOEta   = 0; lMCTree->SetBranchAddress("jeteta",&lOEta);
  float lOPhi   = 0; lMCTree->SetBranchAddress("jetphi",&lOPhi);
  double lVPt   = 0; double lEt = 0; double lPx = 0; double lPy = 0;
  for(int i0 = 0; i0 < lMCTree->GetEntries(); i0++) {
    lMCTree->GetEntry(i0);
    if(lMt < 60)                   continue;
    if(lCharge > 0 && iCharge < 0) continue;
    if(lCharge < 0 && iCharge > 0) continue;
    if(lEta < iPhiMin || lEta > iPhiMax) continue;
    if(lEta < iEtaMin || lEta > iEtaMax) continue;
    lMass->Fill(fabs(lMt));     //lMass->Fill(lPt);
  }
  return lMass;
}
void fillData(RooDataSet *iData,RooRealVar &lXVar,int iCharge,float iPhiMin,float iPhiMax,float iEtaMin,float iEtaMax) { 
  TFile *lFile = new TFile("ZTP.root");
  TTree *lTree = (TTree*) lFile->FindObjectAny("WNtupleIdEffNT");
  int   lCharge = 0; lTree->SetBranchAddress("charge",&lCharge);
  float lEta    = 0; lTree->SetBranchAddress("eta"   ,&lEta);
  float lPhi    = 0; lTree->SetBranchAddress("phi"   ,&lPhi);
  float lMt     = 0; lTree->SetBranchAddress("mt"    ,&lMt);
  float lPt     = 0; lTree->SetBranchAddress("pt"    ,&lPt);
  float lOPt    = 0; lTree->SetBranchAddress("jetpt" ,&lOPt);
  float lOEta   = 0; lTree->SetBranchAddress("jeteta",&lOEta);
  float lOPhi   = 0; lTree->SetBranchAddress("jetphi",&lOPhi);
  for(int i0 = 0; i0 < lTree->GetEntries();i0++) { 
    lTree->GetEntry(i0);
    if(lMt < 60)                   continue;
    if(lCharge > 0 && iCharge < 0) continue;
    if(lCharge < 0 && iCharge > 0) continue;
    if(lPhi < iPhiMin || lPhi > iPhiMax) continue;
    if(lEta < iEtaMin || lEta > iEtaMax) continue;
    lXVar.setVal(fabs(lMt));
    if(lCharge > 0) lXVar.setVal(correct(lMt,lPhi,lOPhi,lEta,lOEta));
    if(lCharge < 0) lXVar.setVal(correct(lMt,lOPhi,lPhi,lOEta,lEta));
    iData->add(RooArgSet(lXVar));
  }
}
void PEs(RooAbsPdf *iGen,RooAbsPdf *iFit,int iN,int iNEvents,RooRealVar &iVar,RooRealVar &iSig,RooRealVar &iMean,
	 RooRealVar &iScale,RooRealVar &iRes) { 
  double iM0 = iMean.getVal(); double iS0 = iSig.getVal();
  //iScale.setVal(iMeanScale); iRes.setVal(iSigScale);
  TRandom1 *lRand = new TRandom1(0xDEADBEEF);
  TNtuple * lDN= new TNtuple( "xxx","xxx","ntot:m_r:m:merr:sig_r:sig:sigerr");
  for(int i0=0;i0<iN;i0++){
    if(i0 % 10 == 0) cout << "+++++++++++++++++++++++++++ running ======> " << i0 << endl;
    int lN    = lRand->Poisson(iNEvents);
    RooDataSet * lSignal  = iGen->generate(iVar,lN);
    iMean.setVal(iM0); iSig.setVal(iS0);
    iFit->fitTo(*lSignal,Strategy(1));//,Save(kTRUE),PrintLevel(1));
    if(iMean.getError() < 0.05) iFit->fitTo(*lSignal,Strategy(2));
    Float_t values[]={
      (Float_t) lN,
      (Float_t) 90.78/iScale.getVal(),
      (Float_t) iMean.getVal(),
      (Float_t) iMean.getError(),
      (Float_t) iSig.getVal(),
      (Float_t) iSig.getVal(),
      (Float_t) iSig.getError()
    };
    lDN->Fill(values);
  }
  TFile *lF = new TFile("XXX.root","RECREATE");
  lDN->Write();
  lF->Close();
}
void Plot(RooAbsPdf *iGen,RooAbsPdf *iFit,int iN,int iNEvents,RooRealVar &iVar,RooRealVar &iSig,RooRealVar &iMean,
	  RooRealVar &iScale,RooRealVar &iRes,RooDataSet *iData=0) { 
  TRandom1 *lRand = new TRandom1(0xDEADBEEF);
  RooDataSet * lSignal  = iGen->generate(iVar,iNEvents);
  if(iData == 0) {
    iFit->fitTo(*lSignal,Strategy(1));
    if(iMean.getError() < 0.05) iFit->fitTo(*lSignal,Strategy(2));
  } else {
    iFit->fitTo(*iData,Strategy(1));
    if(iMean.getError() < 0.05) iFit->fitTo(*iData,Strategy(2));
  }    
  iVar.setBins(30);
  RooPlot *lFrame1 = iVar.frame(RooFit::Title("XXX")) ;
  if(iData == 0) lSignal->plotOn(lFrame1);
  if(iData != 0) iData->plotOn(lFrame1);
  iFit->plotOn(lFrame1);
  TCanvas *iC =new TCanvas("A","A",800,600);
  iC->cd(); lFrame1->Draw();
  iC->SaveAs("Crap.png");
  if(iData != 0) { 
    RooPlot *lFrame2 = iVar.frame(RooFit::Title("XXX")) ;
    iData->plotOn(lFrame2);
    iGen->plotOn(lFrame2);
    TCanvas *iC1 =new TCanvas("B","B",800,600);
    iC1->cd(); lFrame2->Draw();
    iC1->SaveAs("Crap.png");
  }
}
void testResolution() { 
  Prep();
  RooRealVar    lXVar  ("XVar","mass(GeV/c^{2})",100,60,150); lXVar.setBins(1000);
  RooRealVar l1Sigma("sigma1","sigma1",1.6 ,0.,15.);  //l1Sigma.setConstant(kTRUE);
  RooRealVar l2Sigma("sigma2","sigma2",1.6,0.,15.);   l2Sigma.setConstant(kTRUE);
  RooRealVar l3Sigma("sigma3","sigma3",2.90,0.,35.);  l3Sigma.setConstant(kTRUE);
  RooRealVar lN     ("n"     ,"n"     ,1.00,-15,15.); lN.setConstant(kTRUE);
  RooRealVar lExp   ("exp"   ,"exp"   ,-0.003,-15,15.); //lExp.setConstant(kTRUE);
  RooRealVar lR0Mean("xmean","xmean",0,-10,10);    lR0Mean.setConstant(kTRUE);
  RooRealVar lR1Mean("mean","mean",90.8,60,150);   //lR1Mean.setConstant(kTRUE);
  RooVoigtianShape     lGAdd("Add","Add",lXVar,lR1Mean,l1Sigma,l2Sigma,lN,l3Sigma,true);
  
  RooRealVar    lSPar  ("SPar","SPar", 1.,0., 2.);
  RooFormulaVar lXShift("uparshift","@0*@1",RooArgList(lXVar,lSPar));
  TH1F *lMass = getMass(0,-5,5,-1.5,1.5);
  RooDataHist *lMHist = new RooDataHist("M" ,"M" ,RooArgSet(lXVar),lMass);
  RooHistPdf  *lMPdf  = new RooHistPdf ("MH","MH",lXShift,lXVar,*lMHist,5); 
  RooRealVar lGSigma("gsigma","gsigma",1.6 ,0.,15.);  
  RooGaussian   lGaus1("gaus1","gaus1",lXVar,lR0Mean,lGSigma);
  RooFFTConvPdf  lConv("Conv","Conv",lXVar,*lMPdf,lGaus1);

  RooDataSet *lData = new RooDataSet("crap","crap",RooArgSet(lXVar)); 
  fillData(lData,lXVar,0,-5,5,-1.5,1.5);
  lConv.fitTo(*lData,Strategy(2));
  lGSigma.setVal(lGSigma.getVal());
  lSPar.setVal(lSPar.getVal()*1.01);
  //cout << "=====> Check " << l1Sigma.getVal() << " --- " << lR1Mean.getVal() << "----" << lSPar.getVal() << endl;
  PEs(&lConv,&lGAdd,2000,500,lXVar,l1Sigma,lR1Mean,lSPar,lGSigma);
  lData = 0;
  Plot(&lConv,&lGAdd,2000,50000,lXVar,l1Sigma,lR1Mean,lSPar,lGSigma,lData);
}
