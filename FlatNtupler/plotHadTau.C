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
#include "axisToolsDC.hh"
#include "plotToolsDC.hh"
//#include "/data/blue/pharris/CMSSW_5_2_5/src/MitHtt/Common/MitStyleRemix.hh"

int fWId   ;//= 2;
int fQCDId ;//= 4;
std::string *fFreeWeights = 0;

void drawVBFSpec(TTree **iTree,TH1F **iH,TH1F **iHSS,TH1F **iHMT,TH1F **iHLIS,TH1F **iHTIS,TH1F **iHNMT,TH1F **iHMTSS,TH1F **iHTemp,int iN,std::string iVar,std::string iCut,bool is2012,TFile * iFile = 0,std::string iDirName="",bool *iUseScale,TH1F** iHHigh,TH1F** iHLow) { 
  std::string lVBF = "*(njetingap == 2 && mjj >  500 && jdeta > 3.5 )" 
  for(int i0 = 0; i0 < iN; i0++) {
      //iHNMT [i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 > -10)*(iso_1 < 0.1 && iso_2  >  0.795)               "+lVBF," nm_{T}");
      //iHMTSS[i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 > -10)*(iso_1 < 0.3 && iso_2  > -0.50 )               "+lVBF," nm_{T} SS");
    iH    [i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 < 20)*(iso_1 < 0.1 && iso_2  >  0.795 )               "+lVBF," Main");
    iHSS  [i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 < 20)*(iso_1 < 0.3 && iso_2  >  0.50  && iso_1 > 0.1 )"+lVBF," Same Sign");
    iHMT  [i0]   = draw(iVar,iTree[i0],i0,iCut+"*              (mtMVA_1 > 70)*(iso_1 < 0.3 && iso_2  >  0.795 )               "+lVBF," m_{T} > 70 GeV");
    iHNMT [i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 >  0)*(iso_1 < 0.3 && iso_2  >  0.795 )               "+lVBF," No m");
    iHMTSS[i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 > 70)*(iso_1 < 0.3 && iso_2  > -0.50  )               "+lVBF," XNomSS");
    iHTemp[i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*"+fFreeWeights[i0]                                                   ," XFree");
    iHLIS [i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 < 20)*(iso_1 < 0.3 && iso_2  >  0.5 && iso_1 > 0.1)"  ," X2 Jet Loose Iso SS Cut");
    iHTIS [i0]   = draw(iVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 < 20)*(iso_1 < 0.1 && iso_2  >  0.78 )"               ," Y2 Jet Tight Iso SS Cut");
    if(iFile == 0    ) continue;
    iHHigh[i0]   = 0;
    iHLow[i0]    = 0;
    if(!iUseScale[i0]) continue;
    iHHigh[i0]   = draw(iVar+"High",iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 < 20)"," Main High");
    iHLow [i0]   = draw(iVar+"Low" ,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 < 20)"," Main Low");
  }
  //Comput MT Scale Factor
  TH1F *lMTMC  = (TH1F*) iHMT  [0]->Clone("mTTmp");  clear(lMTMC);
  TH1F *lMTMCS = (TH1F*) iHMTSS[0]->Clone("mTTmpS"); clear(lMTMCS);
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 == fWId || i0 == fQCDId) continue;
    lMTMC ->Add(iHMT[i0]);
    lMTMCS->Add(iHMTSS[i0]);
  }
  double lDataInt   = iHMT[iN-1]->Integral(-1,1000) - lMTMC->Integral(-1,1000); 
  double lDataSSInt = iHMTSS[iN-1]->Integral(-1,1000) - lMTMCS->Integral(-1,1000); 
  double lWSF       = float(iHMT[fWId]->Integral(-1,1000))/lDataInt;
  double lWSFSS     = float(iHMTSS[fWId]->Integral(-1,1000))/lDataSSInt;
  cout << "===> W Boson Scale Factor : " << lWSF << " -W- " << lWSFSS << endl;//iHMT[fWId]->Integral(-1,1000)/lDataInt;
  // << " -Data- "<< lDataInt << " - " << iHMT[iN-1]->Integral(-1,1000) << "  MC- " << lMTMC->Integral(-1,1000) << endl;
  if(lWSF   == 0) lWSF   = 1.;
  if(lWSFSS == 0) lWSFSS = 1.;
  iH    [fWId]->Scale(1./lWSF);
  iHSS  [fWId]->Scale(1./lWSF);
  iHMT  [fWId]->Scale(1./lWSF);
  iHNMT [fWId]->Scale(1./lWSF);
  iHMTSS[fWId]->Scale(1./lWSF);
  //!!!!! Do we need to do this
  //iHLIS [fWId]->Scale(1./lWSF);
  //iHTIS [fWId]->Scale(1./lWSF);
  
  //Compute QCD Shape
  TH1F *lSS   = (TH1F*) iHSS  [0]->Clone("SSTmp0"); clear(lSS);
  TH1F *lSSL  = (TH1F*) iHLIS [0]->Clone("SSTmp1"); clear(lSSL);
  TH1F *lSST  = (TH1F*) iHTIS [0]->Clone("SSTmp2"); clear(lSST);
  TH1F *lSSMT = (TH1F*) iHMTSS[0]->Clone("SSTmp3"); clear(lSSMT);
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 == fQCDId) continue;
    lSS  ->Add(iHSS  [i0]);
    lSSL ->Add(iHLIS [i0]);
    lSST ->Add(iHTIS [i0]);
    lSSMT->Add(iHMTSS[i0]);
  }
  iHSS  [fQCDId]->Add(lSS   ,-1); 
  iHMTSS[fQCDId]->Add(lSSMT ,-1); 
  iHLIS [fQCDId]->Add(lSSL  ,-1); 
  iHTIS [fQCDId]->Add(lSST  ,-1); 
  clear(iHMT[fQCDId]);
  for(int i0 = 0; i0 < iHSS  [fQCDId]->GetNbinsX()+1; i0++) if(iHSS  [fQCDId]->GetBinContent(i0) < 0) iHSS  [fQCDId]->SetBinContent(i0,0);
  for(int i0 = 0; i0 < iHMTSS[fQCDId]->GetNbinsX()+1; i0++) if(iHMTSS[fQCDId]->GetBinContent(i0) < 0) iHMTSS[fQCDId]->SetBinContent(i0,0);
  for(int i0 = 0; i0 < iHLIS [fQCDId]->GetNbinsX()+1; i0++) if(iHLIS [fQCDId]->GetBinContent(i0) < 0) iHLIS [fQCDId]->SetBinContent(i0,0);
  for(int i0 = 0; i0 < iHTIS [fQCDId]->GetNbinsX()+1; i0++) if(iHTIS [fQCDId]->GetBinContent(i0) < 0) iHTIS [fQCDId]->SetBinContent(i0,0);
  iHTIS[fQCDId]->Scale(1.06);
  double lTightLooseRatio = iHTIS[fQCDId]->Integral()/iHLIS[fQCDId]->Integral(); 
  cout << "===> Tight/Loose QCD Ratio " << lTightLooseRatio << endl;
  iHSS   [fQCDId]->Scale(lTightLooseRatio);
  iHTIS  [fQCDId]->Scale(iHSS[fQCDId]->Integral()/iHTIS[fQCDId]->Integral());
  iH     [fQCDId] = iHSS  [fQCDId];
  //Generic Template Definition for the loose cuts
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 == fQCDId) continue;
    double pInt = iH[i0]->Integral();
    iH[i0] = iHTemp[i0];
    iH[i0]->Scale(pInt/iH[i0]->Integral());
  }
  //Make Data card
  if(iFile !=0) makeDataCard(iFile,iH,iHHigh,iHLow,iN,iDirName,fString);
  //Additional Crap
  //iHMTSS [fQCDId]->Scale(lTightLooseRatio);
  //iHNMT  [fQCDId] = iHMTSS[fQCDId];

  //Blind
  //for(int i0 = 0; i0 < iH[iN-1]->GetNbinsX()+1; i0++) if(iH[iN-1]->GetXaxis()->GetBinCenter(i0) > 60 && iH[iN-1]->GetXaxis()->GetBinCenter(i0) < 130) iH[iN-1]->SetBinContent(i0,0);
  //Draw the plot
  //draw(iH    ,iN,iVar+"VBFA",iVar);
  //draw(iHSS  ,iN,iVar+"VBFB",iVar);
  //draw(iHMT  ,iN,iVar+"VBFC",iVar);
  //draw(iHLIS ,iN,iVar+"VBFD",iVar);
  //draw(iHTIS ,iN,iVar+"VBFE",iVar);
  //draw(iHMTSS,iN,iVar+"VBFF",iVar);
  //draw(iHNMT ,iN,iVar+"VBFG",iVar);
}
void drawSpec(TTree **iTree,TH1F **iH,TH1F **iHSS,TH1F **iHMT,TH1F **iHNMT,TH1F **iHSSMT,TH1F **iHMTS,int iN,std::string iVar,std::string iCut,std::string iName,TFile * iFile = 0,std::string iDirName="",bool *iUseScale,TH1F** iHHigh,TH1F** iHLow) { 
  std::string lVar = iVar;
  TH1F *lQCDShape = 0;
  for(int i0 = 0; i0 < iN; i0++) {
    lVar = iVar;
    iH   [i0]   = draw(lVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 < 20)"," Main");
    iHSS [i0]   = draw(lVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 < 20)"," Same Sign");
    iHMT [i0]   = draw(lVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 > 70)"," OS m_{T} > 70 GeV");
    iHMTS[i0]   = draw(lVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 > 70)"," SS m_{T} > 70 GeV");
    cout << "====> " << fString[i0] << " -- " << iH[i0]->Integral() << " -- " << fWId << " -- " << fQCDId << endl;
    if(i0 == fQCDId && iName.find("A") == std::string::npos) lQCDShape =  draw(lVar,iTree[i0],iN,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 < 20)"," Same Sign Shape");
    
    //Additional Useless plots
    //iHNMT[i0]   = draw(lVar,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 >  0)"," No m_{T} Cut");
    //iHSSMT[i0]  = draw(lVar,iTree[i0],i0,iCut+"*(q_1*q_2 > 0)*(mtMVA_1 >  0)"," Same Sign");
    
    //Dealing with datacards
    if(iFile == 0    ) continue;
    iHHigh[i0]   = 0;
    iHLow[i0]    = 0;
    if(!iUseScale[i0]) continue;
    iHHigh[i0]   = draw(iVar+"High",iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 < 20)"," Main High");
    iHLow [i0]   = draw(iVar+"Low" ,iTree[i0],i0,iCut+"*(q_1*q_2 < 0)*(mtMVA_1 < 20)"," Main Low");
  }

  //Comput MT Scale Factor => separately for SS and OS
  TH1F *lMTMC  = (TH1F*) iHMT [0]->Clone("mTTmp");  clear(lMTMC);
  TH1F *lMTMCS = (TH1F*) iHMTS[0]->Clone("mTTmpS"); clear(lMTMCS);
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 == fWId || i0 == fQCDId) continue;
    lMTMC ->Add(iHMT [i0]);
    lMTMCS->Add(iHMTS[i0]);
  }
  double lDataInt   = iHMT [iN-1]->Integral(-1,1000) - lMTMC ->Integral(-1,1000); 
  double lWSF       = float(iHMT[fWId]->Integral(-1,1000))/lDataInt;
  double lDataSSInt = iHMTS[iN-1]->Integral(-1,1000) - lMTMCS->Integral(-1,1000); 
  double lWSFSS     = float(iHMTS[fWId]->Integral(-1,1000))/lDataSSInt;
  cout << "===> W Boson    Scale Factor : " << lWSF   << " -W- " << endl;//iHMT[fWId]->Integral(-1,1000)/lDataInt;
  cout << "===> W Boson SS Scale Factor : " << lWSFSS << " -W- " << endl;//iHMT[fWId]->Integral(-1,1000)/lDataInt;
  if(lWSF   == 0) lWSF   = 1.;
  if(lWSFSS == 0) lWSFSS = 1.;
  iH    [fWId]->Scale(1./lWSF);
  iHSS  [fWId]->Scale(1./lWSFSS);
  iHMT  [fWId]->Scale(1./lWSF);//No ss
  iHMTS [fWId]->Scale(1./lWSFSS);

  //Compute QCD Yield
  TH1F *lSS   = (TH1F*) iHSS[0]->Clone("SSTmp");   clear(lSS);
  for(int i0 = 0; i0 < iN-1; i0++) {if(i0 == fQCDId) continue;    lSS  ->Add(iHSS  [i0]);  }
  iHSS  [fQCDId]->Add(lSS  ,-1); 
  //Scale the Shape
  if(lQCDShape != 0) lQCDShape->Scale(iHSS[fQCDId]->Integral());
  if(lQCDShape != 0) iHSS  [fQCDId] = lQCDShape;
  
  for(int i0 = 0; i0 < iHSS  [fQCDId]->GetNbinsX()+1; i0++) if(iHSS  [fQCDId]->GetBinContent(i0) < 0) iHSS  [fQCDId]->SetBinContent(i0,0);
  TH1F *lXSS = (TH1F*) iHSS[fQCDId]->Clone("SSF");  lXSS->Scale(1.06);
  iH   [fQCDId] = lXSS;//iHSS[fQCDId];
  
  //Dealing with Addtional Plots
  //iHNMT [fWId]->Scale(1./lWSF); //No ss
  //iHSSMT[fWId]->Scale(1./lWSFSS);
  //TH1F *lSSMT = (TH1F*) iHSS[0]->Clone("SSTmpmT"); clear(lSSMT);
  //for(int i0 = 0; i0 < iN-1; i0++) {if(i0 == fQCDId) continue;    lSSMT->Add(iHSSMT[i0]);  }
  //iHSSMT[fQCDId]->Add(lSSMT,-1); 
  //clear(iHMT [fQCDId]);
  //clear(iHNMT[fQCDId]);
  //for(int i0 = 0; i0 < iHSSMT[fQCDId]->GetNbinsX()+1; i0++) if(iHSSMT[fQCDId]->GetBinContent(i0) < 0) iHSSMT[fQCDId]->SetBinContent(i0,0);
  //iHSSMT[fQCDId]->Scale(1.06);
  //iHNMT[fQCDId] = iHSSMT[fQCDId];
  
  //Make DataCards ==> If Asked
  if(iFile !=0) makeDataCard(iFile,iH,iHHigh,iHLow,iN,iDirName,fString);
  
  //Blind
  //for(int i0 = 0; i0 < iH[iN-1]->GetNbinsX()+1; i0++) if(iH[iN-1]->GetXaxis()->GetBinCenter(i0) > 100 && iH[iN-1]->GetXaxis()->GetBinCenter(i0) < 150) iH[iN-1]->SetBinContent(i0,0);
 
  //Draw the plot
  //draw(iH   ,iN,iVar+"A"+iName,iVar);
  //draw(iHSS ,iN,iVar+"B"+iName,iVar);
  //draw(iHMT ,iN,iVar+"C"+iName,iVar);
  //draw(iHNMT,iN,iVar+"D"+iName,iVar);
  //draw(iHMTS,iN,iVar+"E"+iName,iVar);
}
//m^2=E^2-p^2
//m^2=E^2-pz^2-px^2-py^2=(E1+E2+E3)^2
//eta = - log(tan(\theta/2))
//atan(2.*Exp(-eta/2)) = theta
std::string lPx    = "(pt_1*cos(phi_1)+pt_2*cos(phi_2)+mvamet*cos(mvametphi))"; 
std::string lPy    = "(pt_1*sin(phi_1)+pt_2*sin(phi_2)+mvamet*sin(mvametphi))"; 
std::string lPz    = "(pt_1*eta_1/abs(eta_1)*abs(cos(2*atan(exp(-eta_1))))+pt_2*eta_2/abs(eta_2)*abs(cos(2.*atan(exp(-eta_2)))))"; 
std::string lE     = "(mvamet + abs(pt_1/sin(2*atan(exp(-eta_1))))+abs(pt_2/sin(2.*atan(exp(-eta_2)))))"; 
std::string lMEff  = "sqrt("+lE+"*"+lE+" - pth*pth - "+lPz+"*"+lPz+")"; 

void plotHadTau(std::string iVar="TMath::Min(abs(jphi_1-jphi_2),6.28-abs(jphi_1-jphi_2))",std::string iCut="vbf",std::string iName="can",std::string iDir="mtau/2012/mtau/",bool is2012=true,int iTauId = 2,bool iEmbed=true) { 
  //void plotHadTau(std::string iVar="pt_2*(1.015 + 0.001 * TMath::Min(TMath::Max(pt_2-45,0),10))",std::string iCut="(nprong_2 == 1 && ngamma_2 > 0  && pt_1 > 20 && iso_1 < 0.1 && iso_2 > 0.795  && pt_1 > 20   && nbtag == 0)",std::string iName="can",std::string iDir="mtau/2012/mtau/",bool is2012=true,int iTauId = 2,bool iEmbed=true) { 
//void plotHadTau(std::string iVar="mvamet*((acos(cos(mvametphi-phi_1)) > TMath::Pi()/2.)+(acos(cos(mvametphi-phi_1)) < TMath::Pi()/2.)*cos(mvametphi-phi_1))",std::string iCut="(pt_1 > 24 && pt_2 > 20 && iso_1 < 0.1 && iso_2 > 0.785 )",std::string iName="can",std::string iDir="2012/etau",int iTauId = 2,bool is2012=true,bool iEmbed=false) { 
  float iLumi = 4800.;
  if(is2012) iLumi = 19100.;
  //SetStyle();
  loadfMap();
  std::stringstream lNameId; //lNameId << "Flat_" << lTauId << "_";
  const int lN = 6;
  std::string lName = iDir+"ntuples/";
  //std::string lName1 = "2011/mutau/ntuples/";
  //std::string lName = "svfit/"+iDir;
  //std::string lName = iDir;

  fWId   = 2;
  fQCDId = 4;
  
  TTree **lTree = new TTree*[lN]; 
  TH1F**lH     = new TH1F*[lN]; 
  TH1F**lHSS   = new TH1F*[lN];
  TH1F**lHIso  = new TH1F*[lN];
  TH1F**lHMT   = new TH1F*[lN]; 
  TH1F**lHNMT  = new TH1F*[lN]; 
  TH1F**lHSSMT = new TH1F*[lN]; 
  TH1F**lHTIS  = new TH1F*[lN]; 
  TH1F**lHLIS  = new TH1F*[lN]; 
  TH1F**lHTemp = new TH1F*[lN]; 
  fString = new std::string[lN]; fWeights = new std::string[lN]; fColor = new int[lN]; fFreeWeights = new std::string[lN];
  lTree[0]     = load(lName+"ztt-mad_select.root");        fString[0] = "Z#rightarrow#tau#tau ";          fColor[0] = 796;//kOrange-3;
  lTree[1]     = load(lName+"ttbar-8TeV_select.root");     fString[1] = "t#bar{t}";                       fColor[1] = 592;//kRed+4;
  lTree[2]     = load(lName+"wjets_select.root");          fString[2] = "W+Jets";                         fColor[2] = 634;//kBlue-5;
  lTree[3]     = load(lName+"zmm_select.root");            fString[3] = "Z#rightarrow#tau#tau fakes";     fColor[3] = kBlue;
  lTree[4]     = load(lName+"data_select.root");           fString[4] = "QCD";                            fColor[4] = 606;//kBlue+3;
  //lTree[5]     = load(lName+"higgs_select.root");          fString[5] = "Higgs ";                         fColor[5] = kBlack;
  lTree[lN-1]  = load(lName+"data_select.root");        fString[lN-1] = "Data"; fColor[lN-1] = kBlack;
  TTree *lEmbTree  = load(lName+"ztt-emb_select.root"); 

  std::stringstream lLumi; lLumi << iLumi;
  if(iCut != "vbf") for(int i0 = 0; i0 < lN;   i0++) fWeights[i0]   = "(pt_1 > 20 && pt_2 > 20 && abs(eta_1) < 2.1 && abs(eta_2) < 2.3 )*"+iCut;
  if(iCut == "vbf") for(int i0 = 0; i0 < lN;   i0++) fWeights[i0]  += "(pt_1 > 20 && pt_2 > 20 && abs(eta_1) < 2.1 && abs(eta_2) < 2.3 )";

  for(int i0 = 0; i0 < lN-1; i0++) if(i0 != fQCDId ) fWeights[i0]  += "*weight*"+lLumi.str();//effweight*puOBweight*"+lLumi.str();
  //Z Scale Factors
  if(!is2012) { 
    fWeights[0]  += "*0.95";
    //fWeights[1]  += "*0.1";
    fWeights[3]  += "*0.95";
  }
  if(is2012) fWeights[1]  += "*1.01";//*0.1";//0.9*0.1";
  if(is2012) fWeights[0]  += "*1.01";
  if(is2012) fWeights[3]  += "*1.01";
  //if(iCut == "vbf") fWeights[2] += "*(weight < 0.1)";
  if(iCut == "vbf") for(int i0 = 0; i0 < lN;   i0++) fFreeWeights[i0]  += "(mjj > 300 && jdeta > 2.5 && iso_1 < 0.1 && iso_2 > 0.795 )";

  if(iEmbed) {
    TH1F *lMC = new TH1F("pXX","pXX",21,0,1000);
    TH1F *lDa = new TH1F("pAA","pAA",21,0,1000);
    std::string lCut = "weight*(pt_1 > 20 && pt_2 > 20 && abs(eta_1) < 2.1 && abs(eta_2) < 2.3 && iso_1 < 0.1 && iso_2 > 0.795 && genmass > 70)*"+lLumi.str();
    lTree[0]->Draw("pt_2>>pXX" ,lCut.c_str());
    lEmbTree->Draw("pt_2>>pAA" ,lCut.c_str());
    cout << "===> " << lMC->Integral() << " - " << lDa->Integral() << endl;
    double lIncWeight = lMC->Integral()/lDa->Integral();
    std::stringstream pSS; pSS << lIncWeight << "*";
    fWeights[0]   = pSS.str()+fWeights[0];
  delete lTree[0];
  lTree[0] = lEmbTree;
  }
  //bool applyVBFVeto

  std::string lVar = iVar;
  TCanvas *lC0 = new TCanvas("A","A",400,400);
  if(iCut != "vbf") drawSpec   (lTree,lH,lHSS,lHMT,lHNMT,lHSSMT,lHTIS              ,lN,lVar,"*( pt_2 > -10.105)",iName );    lC0->cd();
  if(iCut == "vbf") drawVBFSpec(lTree,lH,lHSS,lHMT,lHLIS,lHTIS ,lHNMT,lHSSMT,lHTemp,lN,lVar,"*( pt_2 > -10.105)",is2012);    lC0->cd();
}
