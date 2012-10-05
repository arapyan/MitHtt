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

int fWId   = 2;
int fQCDId = 4;

void drawSpec(TTree **iTree,TH1F **iH,TH1F **iHSS,TH1F **iHMT,TH1F **iHNMT,int iN,std::string iVar,std::string iCut) { 
  for(int i0 = 0; i0 < iN; i0++) {
    iH   [i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mt_1 < 40)"," Main");
    iHSS [i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mt_1 < 40)"," Same Sign");
    iHMT [i0]   = draw(iVar,iTree[i0],i0,iCut+"*              (mt_1 > 70)"," m_{T} > 70 GeV");
    iHNMT[i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mt_1 >  0)"," No m_{T} Cut");
    cout << "====> " << fString[i0] << " -- " << iH[i0]->Integral() << endl;
  }
  //Comput MT Scale Factor
  TH1F *lMTMC = (TH1F*) iHMT[0]->Clone("mTTmp"); clear(lMTMC);
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 == fWId || i0 == fQCDId) continue;
    lMTMC->Add(iHMT[i0]);
  }
  double lDataInt = iHMT[iN-1]->Integral(-1,1000) - lMTMC->Integral(-1,1000); 
  double lWSF     = float(iHMT[fWId]->Integral(-1,1000))/lDataInt;
  cout << "===> W Boson Scale Factor : " << lWSF << " -W- " << endl;//iHMT[fWId]->Integral(-1,1000)/lDataInt;
  // << " -Data- "<< lDataInt << " - " << iHMT[iN-1]->Integral(-1,1000) << "  MC- " << lMTMC->Integral(-1,1000) << endl;
  iH   [fWId]->Scale(1./lWSF);
  iHSS [fWId]->Scale(1./lWSF);
  iHMT [fWId]->Scale(1./lWSF);
  iHNMT[fWId]->Scale(1./lWSF);
  //Compute QCD Shape
  TH1F *lSS = (TH1F*) iHSS[0]->Clone("SSTmp"); clear(lSS);
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 == fQCDId) continue;
    lSS->Add(iHSS[i0]);
  }
  iHSS[fQCDId]->Add(lSS,-1); 
  clear(iHMT[fQCDId]);
  for(int i0 = 0; i0 < iHSS[fQCDId]->GetNbinsX()+1; i0++) if(iHSS[fQCDId]->GetBinContent(i0) < 0) iHSS[fQCDId]->SetBinContent(i0,0);
  iHSS [fQCDId]->Scale(1.11);
  iH   [fQCDId] = iHSS[fQCDId];
  iHNMT[fQCDId] = iHSS[fQCDId];
  //Draw the plot
  draw(iH   ,iN,iVar+"A",iVar,5);
  draw(iHSS ,iN,iVar+"B",iVar,5);
  draw(iHMT ,iN,iVar+"C",iVar,5);
  draw(iHNMT,iN,iVar+"D",iVar,5);
}
void plotHadTau(std::string iVar="eta_2",std::string iCut="(pt_2 > 0)",int iTauId = 2,float iLumi=1200) { //"TMath::Min(abs(phi_1-metphi),2.*TMath::Pi()-abs(phi_1-metphi))",int iId = 0) { 
  SetStyle();
  loadfMap();
  std::stringstream lNameId; //lNameId << "Flat_" << lTauId << "_";
  const int lN = 6;
  std::string lName = "tmp/ntuples/"+lNameId.str();
  
  TTree **lTree = new TTree*[lN]; 
  TH1F**lH    = new TH1F*[lN]; 
  TH1F**lHSS  = new TH1F*[lN]; 
  TH1F**lHMT  = new TH1F*[lN]; 
  TH1F**lHNMT = new TH1F*[lN]; 
  fString = new std::string[lN]; fWeights = new std::string[lN]; fColor = new int[lN];
  lTree[0]  = load(lName+"Xztt-mad_select.root");        fString[0] = "Z#rightarrow#tau#tau ";          fColor[0] = 796;//kOrange-3;
  lTree[1]  = load(lName+"ttbar-8TeV_select.root");     fString[1] = "t#bar{t}";                       fColor[1] = 592;//kRed+4;
  lTree[2]  = load(lName+"wjets_select.root");          fString[2] = "W+Jets";                         fColor[2] = 634;//kBlue-5;
  lTree[3]  = load(lName+"zmm_select.root");            fString[3] = "Z#rightarrow#tau#tau fakes";     fColor[3] = kBlue;
  lTree[4]  = load(lName+"data_select.root");           fString[4] = "QCD";                            fColor[4] = 606;//kBlue+3;
  //lTree[5]  = load(lName+"f11-h120tt-gf.root");         fString[5] = "Higgs ";                         fColor[5] = kRed-4;
  lTree[lN-1]  = load(lName+"data_select.root");       fString[lN-1] = "Data"; fColor[lN-1] = kBlack;
  
  std::stringstream lLumi; lLumi << iLumi;
  for(int i0 = 0; i0 < lN;   i0++) fWeights[i0]   = "(pt_1 > 20 && pt_2 > 20 && iso_1 < 0.1 && iso_2 > 0.795)*"+iCut;
  for(int i0 = 0; i0 < lN-1; i0++) if(i0 != fQCDId) fWeights[i0]  += "*weight*"+lLumi.str();//effweight*puweight*"+lLumi.str();
  //Z Scale Factors
  fWeights[0]  += "*1000.*1.01";
  fWeights[3]  += "*1.01";
  
  //if(lTauId == 1) fWeights[0] += "*(abs(id1_l) + abs(id2_l) > 30)";
  //if(lTauId == 1) fWeights[3] += "*(abs(id1_l) + abs(id2_l) < 30)*(1.1*(abs(eta_2) < 1.5) + 1.02*(abs(eta_2) > 1.5))";
  //if(lTauId == 2) fWeights[0] += "*(abs(id1_l) + abs(id2_l) > 30)";//* 4612.25/1.84668e+08";//*29821.7/1.07845e+09*29821.9/34677.9";
  //if(lTauId == 2) fWeights[3] += "*(abs(id1_l) + abs(id2_l) < 30)";

  std::string lVar = iVar;
  TCanvas *lC0 = new TCanvas("A","A",400,400);
  drawSpec(lTree,lH,lHSS,lHMT,lHNMT,lN,lVar,"*( pt_2 > -10.105)");    lC0->cd();
}
