#include <iostream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "axisTools.hh"
#include "plotTools.hh"
#include "../Common/MitStyleRemix.hh"
#include "plotHadTau.C"


void datacardHadTau(std::string iVar="m_vis",std::string iBaseCut="(pt_2 > 24)",int iTauId = 2,float iLumi=12000) { 
  SetStyle();
  loadfMap();
  std::stringstream lDCName; lDCName << "DataCard_" << iTauId << "_" << iVar << ".root";
  TFile *lFile = new TFile(lDCName.str().c_str(),"RECREATE");
  const int lN = 7;
  std::string lName = "etau2/ntuples/";

  fWId   = 2;
  fQCDId = 5;
  
  TTree **lTree = new TTree*[lN]; 
  TH1F**lH    = new TH1F*[lN]; 
  TH1F**lHSS  = new TH1F*[lN];
  TH1F**lHIso = new TH1F*[lN];
  TH1F**lHMT  = new TH1F*[lN]; 
  TH1F**lHNMT = new TH1F*[lN]; 
  TH1F**lHTIS = new TH1F*[lN]; 
  TH1F**lHScaleHigh = new TH1F*[lN]; 
  TH1F**lHScaleLow  = new TH1F*[lN]; 

  bool *lTauScale = new bool*[lN];
  fString = new std::string[lN]; fWeights = new std::string[lN]; fColor = new int[lN];  
  lTree[0]  = load(lName+"ztt-mad_select.root");        fString[0] = "ZTT";                  fColor[0] = 796;       lTauScale[0] = false; //Should be true
  lTree[1]  = load(lName+"ttbar-8TeV_select.root");     fString[1] = "TT";                   fColor[1] = 592;       lTauScale[1] = false;
  lTree[2]  = load(lName+"wjets_select.root");          fString[2] = "W";                    fColor[2] = 634;       lTauScale[2] = false;
  lTree[3]  = load(lName+"zmm_select.root");            fString[3] = "ZLJ";                  fColor[3] = kBlue;     lTauScale[3] = false; 
  lTree[4]  = load(lName+"zmm_select.root");            fString[4] = "ZLL";                  fColor[4] = kBlue;     lTauScale[4] = false; 
  lTree[5]  = load(lName+"data_select.root");           fString[5] = "QCD";                  fColor[5] = 606;       lTauScale[5] = false; 
  //lTree[5]  = load(lName+"f11-h120tt-gf.root");         fString[5] = "GGH120";             fColor[5] = kRed-4;    lTauScale[5] = false;
  lTree[lN-1]  = load(lName+"data_select.root");       fString[lN-1] = "data_obs";           fColor[lN-1] = kBlack; lTauScale[6] = false;
  
  std::stringstream lLumi; lLumi << iLumi;
  for(int i0 = 0; i0 < lN;   i0++) fWeights[i0]   = "(pt_1 > 20 && pt_2 > 20 && iso_1 < 0.1 && iso_2 > 0.785 )*"+iBaseCut;
  for(int i0 = 0; i0 < lN-1; i0++) if(i0 != fQCDId) fWeights[i0]  += "*puweight*mcweight*0.9*"+lLumi.str();//effweight*puweight*"+lLumi.str();
  //Z Scale Factors
  fWeights[0]  += "*1.01";
  fWeights[1]  += "*1.01*0.1";
  fWeights[3]  += "*1.01*(genmatch < 4)";
  fWeights[4]  += "*1.01*(genmatch > 3)";
  
  std::string lVar = iVar;
  TCanvas *lC0 = new TCanvas("A","A",400,400);
  drawSpec(lTree,lH,lHSS,lHMT,lHNMT,lN,lVar,"*( njets == 0 && pt_2 < 40 && mvamet > 30)",lFile,"SM0",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  drawSpec(lTree,lH,lHSS,lHMT,lHNMT,lN,lVar,"*( njets == 0 && pt_2 > 40 && mvamet > 30)",lFile,"SM1",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  //drawSpec(lTree,lH,lHSS,lHMT,lHNMT,lN,lVar,"*( njet == 1 && pt_2 < 40)",lFile,"SM2",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  //drawSpec(lTree,lH,lHSS,lHMT,lHNMT,lN,lVar,"*( njet == 1 && pt_2 > 40)",lFile,"SM3",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  //drawVBFSpec(lTree,lH,lHSS,lHMT,lHNMT,lN,lVar,"*( njet == 1 && pt_2 > 40)",lFile,"SM3",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
}
