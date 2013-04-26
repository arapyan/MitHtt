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
//#include "axisToolsDC.hh"
//#include "plotToolsDC.hh"
//#include "/data/blue/pharris/CMSSW_5_2_5/src/MitHtt/Common/MitStyleRemix.hh"
#include "plotHadTau.C"


void datacardHadTau(std::string iVar="m_vis",std::string iBaseCut="(pt_2 > 20)",int iTauId = 2,float iLumi=19300) { 
  //  SetStyle();
  loadfMap();
  std::stringstream lDCName; lDCName << "DataCard_" << iTauId << "_" << iVar << ".root";
  if(iTauId == 2) fChanId = "mutau_8TeV";
  if(iTauId == 1) fChanId = "etau_8TeV";
  
  TFile *lFile = new TFile(lDCName.str().c_str(),"RECREATE");
  const int lN = 10;
  std::string lName = "newntuples/mutauisobits/ntuples/";//etau2/ntuples/";

  fWId   = 2;
  fQCDId = 5;
  
  TTree **lTree = new TTree*[lN]; 
  TH1F**lH    = new TH1F*[lN]; 
  TH1F**lHSS  = new TH1F*[lN];
  TH1F**lHSSMT= new TH1F*[lN];
  TH1F**lHIso = new TH1F*[lN];
  TH1F**lHMT  = new TH1F*[lN]; 
  TH1F**lHMTS = new TH1F*[lN]; 
  TH1F**lHMTSS = new TH1F*[lN]; 
  TH1F**lHNMT = new TH1F*[lN]; 
  TH1F**lHTIS = new TH1F*[lN]; 
  TH1F**lHScaleHigh = new TH1F*[lN]; 
  TH1F**lHScaleLow  = new TH1F*[lN]; 

  bool *lTauScale = new bool*[lN];
  fString = new std::string[lN]; fWeights = new std::string[lN+1]; fColor = new int[lN+1]; //Add One for special QCD  
  lTree[0]  = load(lName+"Scaleemb_select.root");                   fString[0] = "ZTT";                  fColor[0] = 796;       lTauScale[0] = true; //Should be true
  lTree[1]  = load(lName+"ttbar-8TeV_select.root");                 fString[1] = "TT";                   fColor[1] = 592;       lTauScale[1] = false;
  lTree[2]  = load(lName+"wjets_select.root");                      fString[2] = "W";                    fColor[2] = 634;       lTauScale[2] = false;
  lTree[3]  = load(lName+"ztt-mad_select.root");                    fString[3] = "ZJ";                   fColor[3] = kBlue;     lTauScale[3] = false; 
  lTree[4]  = load(lName+"ztt-mad_select.root");                    fString[4] = "ZL";                   fColor[4] = kBlue;     lTauScale[4] = false; 
  lTree[5]  = load(lName+"data_select.root");                       fString[5] = "QCD";                  fColor[5] = 606;       lTauScale[5] = false; 
  //ADD a for loop for all the masses 
  lTree[6]  = load(lName+"Scalehtt_gf_sm_125_select.root");         fString[6] = "ggH125";               fColor[6] = kRed-4;    lTauScale[6] = true;
  lTree[7]  = load(lName+"Scalehtt_vbf_et_sm_125_select.root");     fString[7] = "qqH125";               fColor[7] = kRed-4;    lTauScale[7] = true;
  lTree[8]  = load(lName+"Scalehtt_vtth_sm_125_select.root");       fString[8] = "VH125";                fColor[8] = kRed-4;    lTauScale[8] = true;
  lTree[lN-1]  = load(lName+"data_select.root");                    fString[lN-1] = "data_obs";          fColor[lN-1] = kBlack; lTauScale[lN-1] = false;
  
  
  std::stringstream lLumi; lLumi << iLumi;
  for(int i0 = 0; i0 < lN;   i0++) fWeights[i0]   = "(pt_1 > 20 && pt_2 > 20 && iso_1 < 0.1 && iso_2 > 0.785 )*"+iBaseCut;
  for(int i0 = 0; i0 < lN-1; i0++) if(i0 != fQCDId) fWeights[i0]  += "*weight*"+lLumi.str();//effweight*puweight*"+lLumi.str();
  //QCD Shape
  fWeights[lN]   = "(pt_1 > 20 && pt_2 > 20 && iso_1 > 0.2 && iso_1 < 0.5 && iso_2 > 0.785 )*"+iBaseCut;
  //Z Scale Factors
  //fWeights[0]  += "*1.01";
  //fWeights[1]  += "*1.01";
  fWeights[3]  += "*1.01";//*(genTaus == 0 && diLeptons == 0)";
  fWeights[4]  += "*1.01";//*(genTaus == 0 && diLeptons  > 0)";
  //VBF Loose selection
  for(int i0 = 0; i0 < lN;   i0++) fFreeWeights[i0]  += "(mjj > 300 && jdeta > 2.5 && iso_1 < 0.1 && iso_2 > 0.795 )";

  std::string lVar = iVar;
  TCanvas *lC0 = new TCanvas("A","A",400,400);
  drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lN,lVar,"*( njets == 0 && pt_2 < 40 )","A",lFile,"muTau_0jet_low",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lN,lVar,"*( njets == 0 && pt_2 > 40)" ,"B",lFile,"muTau_0jet_high",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lN,lVar,"*( njets >  0 && pt_2 < 40)" ,"C",lFile,"muTau_boost_low",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lN,lVar,"*( njets >  0 && pt_2 > 40)" ,"D",lFile,"muTau_boost_high",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  //Setup the new axis
  fAxis =  {0,20,40,60,80,100,120,140,160,180,200,250,300,350};
  for(int i0 = 0; i0 < fNBins.size(); i0++) fNBins[i0] = 13;
  //Now for the pain
  drawVBFSpec(lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lHMTSS,lN,lVar,"*( njets >  1 && pt_2 > 10)" ,true,lFile,"muTau_vbf",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
}
