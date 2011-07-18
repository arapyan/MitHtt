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

using namespace RooFit;

void plotVoigt()
{
  RooRealVar l1Sigma("sigma1","sigma1",1.523 ,0.,15.);  //l1Sigma.setConstant(kTRUE);     // sigma
  RooRealVar l2Sigma("sigma2","sigma2",0.750,0.,15.);      //l2Sigma.setConstant(kTRUE);      // alpha
  RooRealVar l3Sigma("sigma3","sigma3",2.4952,0.,35.);    l3Sigma.setConstant(kTRUE);    // width
  RooRealVar lN     ("n"     ,"n"     ,3.717,-15,15.);    lN.setConstant(kTRUE);         // n
  RooRealVar lR1Mean("mean","mean",90.836,70,110);          //lR1Mean.setConstant(kTRUE); // m0
  RooRealVar lRXVar("XVar","mass(GeV/c^{2})",90,70,110); //lRXVar.setBins(500);           // m
  //                                        m     m0     sigma   alpha  n   width
  RooVoigtianShape     lGAdd("Add","Add",lRXVar,lR1Mean,l1Sigma,l2Sigma,lN,l3Sigma,true);

  RooPlot *lFrame1 = lRXVar.frame(RooFit::Title("XXX")) ;
  
  TCanvas *iC =new TCanvas("A","A",800,600);
  lGAdd.plotOn(lFrame1);

  l2Sigma.setVal(1);
  lGAdd.plotOn(lFrame1,LineColor(kBlack),LineStyle(kDashed));

  l1Sigma.setVal(1.7);
  lGAdd.plotOn(lFrame1,LineColor(kRed));
  lFrame1->Draw();
}

