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


void datacardHadTau(std::string iVar="m_sv",std::string iBaseCut="(pt_2 > 20)",int iTauId = 2,float iLumi=19300) { 
  //  SetStyle();
  loadfMap();
  std::stringstream lDCName; lDCName << "DataCard_" << iTauId << "_" << iVar << ".root";
  std::string lChan = "muTau";
  std::string lZ   = "zmm";
  if(iTauId == 2) fChanId = "mutau_8TeV";
  if(iTauId == 1) fChanId = "etau_8TeV";
  if(iTauId == 1) lChan  = "eleTau"; 
  if(iTauId == 1) lZ     = "zee";

  TFile *lFile = new TFile(lDCName.str().c_str(),"RECREATE");
  const int lN = 11;
  //std::string lName = "newntuples/mutauisobits/ntuples/";//etau2/ntuples/";
  std::string lName = "../sv/";//etau2/ntuples/";

  fWId   = 2;
  fQCDId = 5;
  
  TTree **lTree = new TTree*[lN]; 
  TH1F**lH     = new TH1F*[lN]; 
  TH1F**lHSS   = new TH1F*[lN];
  TH1F**lHSSMT = new TH1F*[lN];
  TH1F**lHIso  = new TH1F*[lN];
  TH1F**lHMT   = new TH1F*[lN]; 
  TH1F**lHMTS  = new TH1F*[lN]; 
  TH1F**lHMTSS = new TH1F*[lN]; 
  TH1F**lHNMT  = new TH1F*[lN]; 
  TH1F**lHTIS  = new TH1F*[lN]; 
  TH1F**lHTemp = new TH1F*[lN]; 
  TH1F**lHScaleHigh = new TH1F*[lN]; 
  TH1F**lHScaleLow  = new TH1F*[lN]; 

  bool *lTauScale = new bool*[lN];
  fString = new std::string[lN]; fWeights = new std::string[lN+1]; fColor = new int[lN+1]; fFreeWeights = new std::string[lN]; //Add One for special QCD  
  lTree[0]  = load(lName+"ScaleScaleemb_select_sv.root");                   fString[0] = "ZTT";                  fColor[0] = 796;       lTauScale[0] = true; //Should be true
  lTree[1]  = load(lName+"ttbar-8TeV_select_sv.root");                 fString[1] = "TT";                   fColor[1] = 592;       lTauScale[1] = false;
  lTree[2]  = load(lName+"wjets_skim_select.root");                    fString[2] = "W";                    fColor[2] = 634;       lTauScale[2] = false;
  lTree[3]  = load(lName+lZ+"_select_sv.root");                        fString[3] = "ZJ";                   fColor[3] = kBlue;     lTauScale[3] = false; 
  lTree[4]  = load(lName+lZ+"_select_sv.root");                        fString[4] = "ZL";                   fColor[4] = kBlue+2;   lTauScale[4] = false; 
  lTree[5]  = load(lName+"data_skim_select.root");                     fString[5] = "QCD";                  fColor[5] = 606;       lTauScale[5] = false; 
  lTree[6]  = load(lName+"ewk-8TeV_select_sv.root");                   fString[6] = "VV";                   fColor[6] = 634;       lTauScale[6] = false; 
  //ADD a for loop for all the masses 
  lTree[7]  = load(lName+"Scalehtt_gf_sm_125_select_sv.root");              fString[7] = "ggH125";               fColor[7] = kRed-4;    lTauScale[7] = true;
  lTree[8]  = load(lName+"Scalehtt_vbf_et_sm_125_select_sv.root");          fString[8] = "qqH125";               fColor[8] = kRed-4;    lTauScale[8] = true;
  lTree[9]  = load(lName+"Scalehtt_vtth_sm_125_select_sv.root");            fString[9] = "VH125";                fColor[9] = kRed-4;    lTauScale[9] = true;
  lTree[lN-1]  = load(lName+"data_skim_select.root");                       fString[lN-1] = "data_obs";          fColor[lN-1] = kBlack; lTauScale[lN-1] = false;
  //1-0.86 2-0.93
  std::stringstream lLumi; lLumi << iLumi;
  for(int i0 = 0; i0 < lN;   i0++) fWeights[i0]   = "(pt_1 > 20 && pt_2 > 20 && iso_1 < 0.1 && "+fIso+" && nbtag == 0 )*"+iBaseCut;
  for(int i0 = 0; i0 < lN-1; i0++) if(i0 != fQCDId) fWeights[i0]  += "*weight*"+lLumi.str();//effweight*puweight*"+lLumi.str();
  //QCD Shape
  fWeights[lN]   = "(pt_1 > 20 && pt_2 > 20 && iso_1 > 0.2 && iso_1 < 0.5 && "+fIso+" )*"+iBaseCut;
  //ETau
  if(iTauId == 1)  for(int i0 = 0; i0 < lN+1;   i0++) fWeights[i0]  += "*(pt_1 > 24)";
   //Z Scale Factors
  //fWeights[0]  += "*1.01";
  //fWeights[1]  += "*1.01";
  fWeights[3]  += "*(genmatch > 3)";
  fWeights[4]  += "*(genmatch < 3)";
  //VBF Loose selection
  for(int i0 = 0; i0 < lN;   i0++) fFreeWeights[i0]  += "(mjj > 300 && jdeta > 2.5 && iso_1 < 0.1 &&  "+fIso+" && mtMVA_1 < 20 )";
  fFreeWeights[fWId] += "*(weight < 0.0003)";
  std::string lMET = "";
  if(iTauId == 1) lMET = "*(mvamet > 30)";
  std::string lVar = iVar;
  TCanvas *lC0 = new TCanvas("A","A",400,400);
  drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lN,lVar,"*( njets == 0 && pt_2 < 40 )","A",lFile,lChan+"_0jet_low",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lN,lVar,"*( njets == 0 && pt_2 > 40)" ,"B",lFile,lChan+"_0jet_high",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lN,lVar,"*( njets >  0 && pt_2 < 40)"+lMET ,"C",lFile,lChan+"_boost_low",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
  drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lN,lVar,"*( njets >  0 && pt_2 > 40)"+lMET ,"D",lFile,lChan+"_boost_high",lTauScale,lHScaleHigh,lHScaleLow);    lC0->cd();
 
 //Setup the new axis  
  fAxis =  new double[14]; 
  for(int i0 = 0; i0 < 10; i0++) fAxis[i0]       = i0*20.;
  for(int i0 = 0; i0 < 4;  i0++) fAxis[i0+10]    = i0*50.+200.;
  for(int i0 = 0; i0 < fNBins.size(); i0++) fNBins[i0] = 13;
 
  //Setup the loose pre-selecxxtion
  for(int i0 = 0; i0 < lN;   i0++)                    fWeights[i0]   = "(pt_1 > 20 && pt_2 > 20 )*"+iBaseCut;
  for(int i0 = 0; i0 < lN-1; i0++) if(i0 != fQCDId)   fWeights[i0]  += "*weight*"+lLumi.str();//effweight*puweight*"+lLumi.str();
  if(iTauId == 1)  for(int i0 = 0; i0 < lN+1;   i0++) fWeights[i0]  += "*(pt_1 > 24)";
  //Now for the pain
  drawVBFSpec(lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHMTS,lHMTSS,lHTemp,lN,lVar,"*( njets >  1 && pt_2 > 10)" ,true,lFile,lChan+"_vbf",lTauScale,lHScaleHigh,lHScaleLow);    
}
