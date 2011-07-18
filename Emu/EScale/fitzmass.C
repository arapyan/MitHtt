#include <sstream>
#include <iostream>
#include "TMath.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TTree.h"
#include "TNtuple.h"
#include "RooRealVar.h"
#include "RooArgSet.h"
#include "RooArgList.h"
#include "RooDataSet.h"
#include "RooExponential.h"
#include "RooLandau.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooAddPdf.h"
#include "RooGaussian.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooHistFunc.h"
#include "RooMoment.h"
#include "RooFitResult.h"
#include "RooExtendPdf.h"
#include "RooGenericPdf.h"
#include "RooBreitWigner.h"
#include "RooBifurGauss.h"
#include "RooProdPdf.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"
#include "RooWorkspace.h"
#include "RooConstVar.h"


using namespace RooFit;

TH1D *puweights = 0;
float puweight(float npu) {
  if (npu<0) return 1.0;
  return puweights->GetBinContent(puweights->FindFixBin(npu));
}


void fitzmass() { 
  
  
  gROOT->Macro("MitStyle.C");
  gStyle->SetErrorX(0); 
  gStyle->SetOptStat(0);
  gROOT->ForceStyle();  
  
  
//  TFile *filepuest = new TFile("/home/bendavid/cms/root/pudists/puweightsDec22.root","READ");
  TFile *filepuest = new TFile("/home/bendavid/cms/root/pudists/PuWeightsJun12.root","READ");
  puweights= (TH1D*) filepuest->Get("pileup");

  
  
  gSystem->cd("./zmassfitsJun12");
//  gSystem->cd("./muresultsRegion2fixed");

    TFile *fdata = new TFile("/home/bendavid/cms/hist/hgg-v0/MergedData.root","READ");
    TDirectory *hdir = (TDirectory*) fdata->FindObjectAny("HGGModClassic");
    TTree *hdata = (TTree*)hdir->FindObjectAny("hHggNtuple");
 
    TFile *fmc = new TFile("/home/bendavid/cms/hist/hgg-v0/t2mit/filefi/merged/hgg-v0_s11-zeem20-v11-pu_noskim.root","READ");
    TDirectory *hdirmc = (TDirectory*) fmc->FindObjectAny("HGGModClassic");
    TTree *hmc = (TTree*)hdirmc->FindObjectAny("hHggNtuple");
    
    
  gStyle->SetOptStat(1110);
  
 
  
  //return;
  
  TCut loosesel = "pt1>40.0 && pt2>30.0 && nvtx>0";
  TCut vloosesel = "pt1>30.0 && pt2>25.0 && nvtx>0";

  TCut basesel = loosesel && "passelectronvetoconv1 && passelectronvetoconv2";
  TCut isb1 = "abs(sceta1)<1.5";
  TCut isb2 = "abs(sceta2)<1.5";
  
  TCut cat1l = vloosesel && (isb1&&isb2) && "(r91>0.94&&r92>0.94)";
  TCut cat2l = vloosesel && (isb1&&isb2) && "!(r91>0.94&&r92>0.94)";
  TCut cat3l = vloosesel && !(isb1&&isb2) && "(r91>0.94&&r92>0.94)";
  TCut cat4l = vloosesel && !(isb1&&isb2) && "!(r91>0.94&&r92>0.94)";
  
  TCut looseselebeb = loosesel && "iseb1&&iseb2";
  TCut looseselebee = loosesel && "((iseb1+iseb2)==1)";  
  TCut looseseleeee = loosesel && "isee1&&isee2";
  TCut looseseleeex = loosesel && "((isee1+isee2)>0)";
  TCut looseselebex = loosesel && "!(iseb1&&iseb2)";  

  TCut baseselebeb = basesel && "iseb1&&iseb2";
  TCut cat1 = basesel && (isb1&&isb2) && "(r91>0.94&&r92>0.94)";
  TCut cat2 = basesel && (isb1&&isb2) && "!(r91>0.94&&r92>0.94)";
  TCut cat3 = basesel && !(isb1&&isb2) && "(r91>0.94&&r92>0.94)";
  TCut cat4 = basesel && !(isb1&&isb2) && "!(r91>0.94&&r92>0.94)";
  
  TCut cat2a = basesel && "(iseb1&&iseb2)&&(r91<0.94&&r92<0.94)";    
  
  TCut higgscut = cat1 && "abs(genhz-vtxz)<0.5";
  TCut puweight = "puweight(ngenvtx-1)";
  //TCut higgscut = cat1;
  
   
  //define categories
  std::vector<TString> catnames;
  std::vector<TCut> catcuts;
  
 catnames.push_back("cat1");
 catnames.push_back("cat2");
 catnames.push_back("cat3");
catnames.push_back("cat4");
  
  catcuts.push_back(cat1l);
 catcuts.push_back(cat2l);
 catcuts.push_back(cat3l);
 catcuts.push_back(cat4l);   

  std::vector<double> cbsigmcs;
  std::vector<double> cbsigdatas;
  std::vector<double> cbsigmcerrs;
  std::vector<double> cbsigdataerrs;  
  
   //return;
  
  RooRealVar pt1("pt1","pT_1",0.0,10000.0,"GeV");
  RooRealVar pt2("pt2","pT_2",0.0,10000.0,"GeV");
  RooRealVar sceta1("sceta1","sceta1",-1000.0,1000.0);
  RooRealVar sceta2("sceta2","sceta2",-1000.0,1000.0);
  RooRealVar hmass("hmassvtx","m_{#gamma#gamma}",60.0,200.0,"GeV");
  RooRealVar nvtx("nvtx","nvtx",0.0,1000.0);
  RooRealVar passelectronvetoconv1("passelectronvetoconv1","passelectronvetoconv1",-2.0,2.0);
  RooRealVar passelectronvetoconv2("passelectronvetoconv2","passelectronvetoconv2",-2.0,2.0);
  RooRealVar iseb1("iseb1","iseb1",-2.0,2.0);
  RooRealVar iseb2("iseb2","iseb2",-2.0,2.0);
  RooRealVar isee1("isee1","isee1",-2.0,2.0);
  RooRealVar isee2("isee2","isee2",-2.0,2.0);  
  RooRealVar r91("r91","r91",0.0,10000.0);
  RooRealVar r92("r92","r92",0.0,10000.0);
  RooRealVar genhz("genhz","genhz",-10000,10000.0);
  RooRealVar vtxz("vtxz","vtxz",-10000,10000.0);
  RooRealVar ngenvtx("ngenvtx","ngenvtx",-10000,10000.0);
  
  
  hmass.setRange("zrange",60.0,120.0);
  hmass.setRange("higgsrange",100,130);
  //hmass.setRange("convRange",-20.0,20.0);  
  
  hmass.setBins(100);
  hmass.setBins(10000,"cache");

  //TFile *fhist = new TFile("datahist.root","RECREATE");
  TH1F *hmasshist = new TH1F("hmasshist","hmasshist",100,100.0,200.0);
  //hdata->Draw("hmass>>hmasshist",basesel);
  //hmasshist->Write();
  //fhist->Close();
  


  RooArgSet varlist;
  varlist.add(pt1);
  varlist.add(pt2);
  varlist.add(sceta1);
  varlist.add(sceta2);
  varlist.add(hmass);
  varlist.add(nvtx);
  varlist.add(passelectronvetoconv1);
  varlist.add(passelectronvetoconv2);
  varlist.add(iseb1);
  varlist.add(iseb2);
  varlist.add(isee1);
  varlist.add(isee2);
  varlist.add(r91);
  varlist.add(r92);
  varlist.add(genhz);
  varlist.add(vtxz);
  varlist.add(ngenvtx);
  
  RooFormulaVar puweightf("puweightv","puweightv","puweight(ngenvtx-1)",RooArgList(ngenvtx));

  

  
  
  RooRealVar mz("mz","mz",90,60,120);
  mz.removeRange();

  RooRealVar m0("m0","m0",0.0,-5.0,5.0);
  m0.removeRange();   
//   
  RooRealVar sigma("sigma","sigma",1.4,0.0,100.0);
  sigma.removeRange();  

  RooRealVar sigma2("sigma2","sigma2",1.0,0.0,10.0);
  sigma2.removeRange();  

  RooRealVar fgaus("fgaus","fgaus",0.9,0.5,1.0);    
  
  
  RooRealVar alpha("alpha","alpha",1.0,0.0,10.0);
  alpha.removeRange();    

  RooRealVar n("n","n",1.0,0.0,1000.0);
  n.removeRange();   
  

  
  
  
  const double widthzpdg = 2.4952;
  const double masszpdg = 91.1876;
  
  RooBreitWigner zbw("zbw","zbw",hmass,mz,RooConst(widthzpdg));
  RooCBShape zcb("zcb","zcb",hmass,RooConst(0),sigma,alpha,n);
  RooCBShape zcbnom("zcbnom","zcbnom",hmass,mz,sigma,alpha,n);
  
  RooFFTConvPdf zcbbw("zcbbw","zcbbw",hmass,zbw,zcb); 
  zcbbw.setBufferFraction(1.5); 
    
  for (UInt_t icat=0; icat<catcuts.size(); ++icat) {
    
     RooDataSet zdata(TString("zdata")+catnames.at(icat),"zdata",varlist,RooFit::Import(*hdata),RooFit::Cut(catcuts.at(icat))); 
    RooDataSet zmc(TString("zmc") +catnames.at(icat),"zmc",varlist,RooFit::Import(*hmc),RooFit::Cut(catcuts.at(icat))); 
    zmc.addColumn(puweightf);
    RooDataSet zmcw(TString("zmcw")+catnames.at(icat),"zmcw",*zmc.get(),RooFit::Import(zmc),RooFit::WeightVar("puweightv")); 
    
//     alpha.setVal(1.0);
//     n.setVal(1.0);
    
    if (icat>1) sigma.setVal(3.5);  
    
    zbw.fitTo(zmcw,Strategy(2),Minos(kFALSE),Range("zrange"),SumW2Error(kTRUE));
    zcbbw.fitTo(zmcw,Strategy(2),Minos(kFALSE),Range("zrange"),SumW2Error(kTRUE));
    
    
    TCanvas *czfit = new TCanvas;
    TString mcplotname = TString("zfitmc") + catnames.at(icat) + TString(".eps");
    RooPlot *zplot = hmass.frame(Bins(100),Range("zrange"));
    zmc.plotOn(zplot);
    zcbbw.plotOn(zplot,RooFit::LineColor(kBlue));  
    zplot->SetTitle(catnames.at(icat));      
    zcbbw.paramOn(zplot,Layout(0.55,0.93,0.7));                                   
    zplot->Draw();    
    TLegend *legmc = new TLegend(0.62,0.75,0.92,0.9);  
    legmc->AddEntry(zplot->getObject(0),"MC","LPE");
    legmc->AddEntry(zplot->getObject(1),"Fit","L");
    legmc->SetBorderSize(0);
    legmc->SetFillStyle(0);
    legmc->Draw();    
    czfit->SaveAs(mcplotname);
    
    double cbsigmc = sigma.getVal();
    cbsigmcs.push_back(sigma.getVal());
    cbsigmcerrs.push_back(sigma.getError());
    
    if (icat==1) {
      alpha.setVal(1.0);
      n.setVal(100.0);
      sigma.setVal(1.5);
    }

    
    if (icat==2) {
      alpha.setVal(1.0);
      n.setVal(1.0);    
      sigma.setVal(4.35);   
    }
    if (icat==3) sigma.setVal(4.5);  
    
    if (icat!=2) zbw.fitTo(zdata,Strategy(2),Minos(kTRUE),Range("zrange"));
    zcbbw.fitTo(zdata,Strategy(2),Minos(kTRUE),Range("zrange"));
    
    double cbsigdata = sigma.getVal();
    cbsigdatas.push_back(sigma.getVal());
    cbsigdataerrs.push_back(sigma.getError());
    
    TCanvas *czfitdata = new TCanvas;
    TString dataplotname = TString("zfitdata") + catnames.at(icat) + TString(".eps");    
    RooPlot *zplotdata = hmass.frame(Bins(100),Range("zrange"),NormRange("zrange"));
    zdata.plotOn(zplotdata);
    zcbbw.plotOn(zplotdata,RooFit::LineColor(kBlue));  
    zplotdata->SetTitle(catnames.at(icat));
    zcbbw.paramOn(zplotdata,Layout(0.55,0.93,0.7));                       
    zplotdata->Draw();
    TLegend *legdata = new TLegend(0.62,0.75,0.92,0.9);  
    legdata->AddEntry(zplotdata->getObject(0),"Data (239.1/pb)","LPE");
    legdata->AddEntry(zplotdata->getObject(1),"Fit","L");
    legdata->SetBorderSize(0);
    legdata->SetFillStyle(0);
    legdata->Draw();    
    czfitdata->SaveAs(dataplotname);
      
  }
    
  for (UInt_t icat=0; icat<catcuts.size(); ++icat) { 
    double cbsigmc = cbsigmcs.at(icat);
    double cbsigdata = cbsigdatas.at(icat);    
    double cbsigmcerr = cbsigmcerrs.at(icat);
    double cbsigdataerr = cbsigdataerrs.at(icat);
    
    double cbsmear = sqrt(cbsigdata*cbsigdata-cbsigmc*cbsigmc);
    double cbsmearerr = (1.0/cbsmear)*sqrt(cbsigdata*cbsigdata*cbsigdataerr*cbsigdataerr + cbsigmc*cbsigmc*cbsigmcerr*cbsigmcerr);
    
    printf("cat = %i, cbsigmc = %5f +- %5f, cbsigdata = %5f +- %5f, cbsmear = %5f +- %5f\n",icat,cbsigmc,cbsigmcerr,cbsigdata,cbsigdataerr, cbsmear, cbsmearerr);
  }
    
}

