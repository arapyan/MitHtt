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
// #include "../Efficiency/Efficiency/interface/IdTools1.hh"
// #include "../NYStyle/test/NYStyle.h"

using namespace RooFit;

TString outputdir;
TString infile;

FILE *of=0;

void fitMass(double iXmin,double iXmax,double &iRes,double &iShift,double &iResErr,double &iShiftErr) { 
  // Prep();
  RooRealVar l1Sigma("sigma1","sigma1",1.523 ,0.,15.);  //l1Sigma.setConstant(kTRUE);     // sigma
  RooRealVar l2Sigma("sigma2","sigma2",0.750,0.,15.);      //l2Sigma.setConstant(kTRUE);      // alpha
  RooRealVar l3Sigma("sigma3","sigma3",2.4952,0.,35.);    l3Sigma.setConstant(kTRUE);    // width
  RooRealVar lN     ("n"     ,"n"     ,3.717,-15,15.);    lN.setConstant(kTRUE);         // n
  RooRealVar lR1Mean("mean","mean",90.836,70,110);          //lR1Mean.setConstant(kTRUE); // m0
  RooRealVar lRXVar("XVar","mass(GeV/c^{2})",90,70,110); //lRXVar.setBins(500);           // m
  //                                        m     m0     sigma   alpha  n   width
  RooVoigtianShape     lGAdd("Add","Add",lRXVar,lR1Mean,l1Sigma,l2Sigma,lN,l3Sigma,true);
  RooDataSet *lData = new RooDataSet("crap"  ,"crap",RooArgSet(lRXVar)); 
  RooDataSet *lXData = new RooDataSet("xcrap","xcrap",RooArgSet(lRXVar)); 
  RooDataSet *lYData = new RooDataSet("ycrap","ycrap",RooArgSet(lRXVar)); 
  TFile *lFile = new TFile(infile);
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  float lmass  = 0; lTree->SetBranchAddress("mass"    ,&lmass);
  float lpt1   = 0; lTree->SetBranchAddress("pt1"    ,&lpt1);
  float leta1  = 0; lTree->SetBranchAddress("eta1"   ,&leta1);
  float lphi1  = 0; lTree->SetBranchAddress("phi1"   ,&lphi1);
  float lpt2   = 0; lTree->SetBranchAddress("pt2" ,&lpt2);
  float leta2  = 0; lTree->SetBranchAddress("eta2",&leta2);
  float lphi2  = 0; lTree->SetBranchAddress("phi2",&lphi2);

  for(int i0 = 0; i0 < lTree->GetEntries();i0++) { 
    lTree->GetEntry(i0);

    if(lmass<70 || lmass>110) continue;

    if(leta1 < iXmin || leta1 > iXmax) continue;
    if(leta2 < iXmin || leta2 > iXmax) continue;

    lRXVar.setVal(lmass);
    lYData->add(RooArgSet(lRXVar));
  }

  RooDataHist *lYDataBinned = lYData->binnedClone("mass","mass");
  lGAdd.fitTo(*lYDataBinned,Strategy(1),Minos(),NumCPU(8));
  // if(lR1Mean.getError() < 0.01) lGAdd.fitTo(*lYData,Strategy(2),Minos());
  lRXVar.setBins(60);
  iShift = lR1Mean.getVal()/91.2; iShiftErr = lR1Mean.getError()/91.2;
  iRes   = l1Sigma.getVal(); iResErr   = l1Sigma.getError();

  fprintf(of,"%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f%10.4f\n",iXmin,
	  iShift,iShiftErr,iRes,iResErr,l2Sigma.getVal(),l2Sigma.getError(),lN.getVal(),lN.getError());

  
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

  TCanvas *iC =new TCanvas("A","A",800,600);
  iC->cd(); lFrame1->Draw();

  // // write out the fits for each eta bin
  // char canname[99];
  // sprintf(canname,"%.2f_Crap.png",iXmin);
  // iC->SaveAs(TString(canname).ReplaceAll("-","m"));

  return;
  
}

void escalefit(TString conf) { 

  int iNBin;
  double ldX;// = 2*2.5;
  // double ldX = 2*3.14159265;
  double lx0;// = -2.5;//3.14159265;

  
  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { state++; continue;}

    if(state==0) {
      stringstream ss(line);
      ss >> infile;
      continue;
    }
    if(state==1) {
      stringstream ss(line);
      ss >> outputdir;
      continue;
    }
    if(state==2) {
      stringstream ss(line);
      ss >> iNBin >> ldX >> lx0;
      continue;
    }
  }
  ifs.close();

  double lInterval = ldX/iNBin;
  double lXVal[iNBin];  double lYVal[iNBin]; double lXEVal[iNBin]; double lYEVal[iNBin];
  double lRYVal[iNBin]; double lRYEVal[iNBin];
  of = fopen("results.txt","w");
  for(int i0 = 0; i0 < iNBin; i0++) { 
    lx0 += lInterval;
    lXVal[i0]  = lx0-lInterval/2.;
    lXEVal[i0] = lInterval/2.;
    //                          resol.       shift      res. err.    shift err
    fitMass(lx0-lInterval,lx0,lRYVal[i0],lYVal[i0],lRYEVal[i0],lYEVal[i0]);
  }
  fclose(of);

  TGraphErrors *lGraph1  = new TGraphErrors(iNBin,lXVal,lYVal,lXEVal,lYEVal); lGraph1->SetLineColor(kRed); lGraph1->SetMarkerColor(kRed); lGraph1->SetMarkerStyle(20);
  TGraphErrors *lGraphR1 = new TGraphErrors(iNBin,lXVal,lRYVal,lXEVal,lRYEVal); lGraphR1->SetLineColor(kRed); lGraphR1->SetMarkerColor(kRed); lGraphR1->SetMarkerStyle(20);
  lGraph1->SetLineStyle(kDashed);     lGraphR1->SetLineStyle(kDashed);
  lGraph1->GetYaxis()->SetRangeUser(0.95,1.05); lGraphR1->GetYaxis()->SetRangeUser(0,5.0);

  lGraph1->GetXaxis()->SetTitle("#phi"); lGraph1->GetYaxis()->SetTitle("Energy Scale");
  TCanvas *lC0= new TCanvas("Res"  ,"Res",800,600);
  lGraph1->Draw("ap");

  lGraphR1->GetXaxis()->SetTitle("#phi"); lGraphR1->GetYaxis()->SetTitle("Additional Smearing (GeV/c^{2})");
  TCanvas *lC1= new TCanvas("Shift","Shift",800,600);
  lGraphR1->Draw("ap");
  
  lGraph1 ->SetName("ScalePlus");  lGraph1->SetTitle("ScalePlus");
  lGraphR1->SetName("ResPlus");    lGraphR1->SetTitle("ResPlus");

  TString prefix;
  if(infile.Contains("data-")) prefix = "data-";
  else if(infile.Contains("zee-")) prefix = "mc-";
  else assert(0);
  TFile *lFile = new TFile(outputdir+"/"+prefix+"EnergyScale.root","RECREATE");
  lGraph1->Write();
  lGraphR1->Write();
}    
