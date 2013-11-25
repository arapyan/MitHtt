#include <iostream>
#include <string>
#include <sstream>
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom1.h"
#include "TGraphErrors.h"
#include "TLine.h"
#include "axisTools.hh"
#include "plotTools.hh"
#include "../Common/MitStyleRemix.hh"
//#include "HttStyles.h"
//#include "HttStyles.cc"

void drawSpec(TTree **iTree,TH1F **iHMM,TH1F **iHMMOS,TH1F **iHMMSST,TH1F **iHMMSSL,int iN,std::string iVar,int iCat,TFile * iFile = 0,std::string iDirName="",int iUseScale=0,int *ptrew=0,TH1F** iHHigh=0,TH1F** iHLow=0,TH1F** iHHVBF=0,TH1F** iHMMSSTVBF=0,TH1F** iHTheory=0,TH1F **iLTheory=0) { 
  cout << "=====>  -- " << fQCDId << endl;
  for(int i0 = 0; i0 < iN; i0++) {
    if(iCat==1)
      {
	iHMM   [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," Main");
	iHHVBF   [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(jdeta>2.0 && mjj>200 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," Main VBF");
	iHMMSSTVBF [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 > 0)*(!(byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," Same Sign Tight");
	iHMMSST [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(jdeta>2.0 && mjj>200 && njetingap==0)*(q_1*q_2 > 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," Same Sign Tight");
	iHMMOS [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(!(byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))*(jdeta>2.0 && mjj>200 && njetingap==0)"," Loose Iso");
	iHMMSSL[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 > 0)*(!(byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))*(jdeta>2.0 && mjj>200 && njetingap==0)"," Same Sign Loose");
      }
    else
      {
	iHMM   [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2<0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0)) "," Main");
	iHMMSST [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 > 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," Same Sign Tight");
	iHMMOS [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 <0)*(!(byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," Loose Iso");
	iHMMSSL[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 > 0)*(!(byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2 <1.0))"," Same Sign Loose");
      }
    cout << "====> " << fString[i0] << " -- " << iHMM[i0]->Integral() << endl;
    cout << "====> " << fString[i0] << " -- " << iHMMOS[i0]->Integral() << endl;
    if(iFile == 0    ) continue;
    iHHigh[i0]   = 0;
    iHLow[i0]    = 0;
    iHTheory[i0]   = 0;
    iLTheory[i0]    = 0;
    if(!iUseScale || i0==iN-1 || i0==fQCDId || i0==fWId || i0==1) continue;
    if(iCat==1)
      {   
	//iHHigh[i0]   = draw(iVar+"high",iCat,iTree[i0],i0,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," M",1);
	//iHLow [i0]   = draw(iVar+"low" ,iCat,iTree[i0],i0,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," K",2);
	iHHigh[i0]   = draw("m_vishigh",iCat,iTree[i0],i0,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," M",1);
	iHLow [i0]   = draw("m_vislow" ,iCat,iTree[i0],i0,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," K",2);
	if(ptrew[i0])
	  {
	    iHTheory[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," L",0,1);
	    iLTheory[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," N",0,2);
	  }
      }
    else
      {
	//iHHigh[i0]   = draw(iVar+"high",iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0)))"," M",1);
	//iHLow [i0]   = draw(iVar+"low" ,iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0)))"," K",2);
	iHHigh[i0]   = draw("m_vishigh",iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0)))"," M",1);
	iHLow [i0]   = draw("m_vislow" ,iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0)))"," K",2);
	if(ptrew[i0])
	  {
	    iHTheory[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0)))"," L",0,1);
	    iLTheory[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0)))"," N",0,2);
	  }
      }
  }
  
  //Compute QCD Shape
  TH1F *lHHVBF=0;
  TH1F *lMMOS = (TH1F*) iHMMOS[0]->Clone("OSTmp"); clear(lMMOS);
  TH1F *lMMSSL = (TH1F*) iHMMSSL[0]->Clone("SSLTmp"); clear(lMMSSL);
  TH1F *lMMSST = (TH1F*) iHMMSST[0]->Clone("SSTTmp"); clear(lMMSST);
  if(iCat==1)
    {
      lHHVBF = (TH1F*) iHMMSSTVBF[0]->Clone("SSLTmp"); clear(lHHVBF);
    }
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 == fQCDId) continue;
    if(i0 >5 && i0<iN-1 && i0!=30) continue;
    lMMOS->Add(iHMMOS[i0]);
    lMMSSL->Add(iHMMSSL[i0]);
    lMMSST->Add(iHMMSST[i0]);
    if(iCat==1)
      lHHVBF->Add(iHMMSSTVBF[i0]);
  }
  iHMMOS[fQCDId]->Add(lMMOS,-1); 
  iHMMSSL[fQCDId]->Add(lMMSSL,-1); 
  iHMMSST[fQCDId]->Add(lMMSST,-1); 
  if(iCat==1)
    iHMMSSTVBF[fQCDId]->Add(lHHVBF,-1);
  
  for(int i0 = -1; i0 < 1000; i0++) if(iHMMOS[fQCDId]->GetBinContent(i0) < 0) iHMMOS[fQCDId]->SetBinContent(i0,0);
  for(int i0 = -1; i0 < 1000; i0++) if(iHMMSST[fQCDId]->GetBinContent(i0) < 0) iHMMSST[fQCDId]->SetBinContent(i0,0);
  for(int i0 = -1; i0 < 1000; i0++) if(iHMMSSL[fQCDId]->GetBinContent(i0) < 0) iHMMSSL[fQCDId]->SetBinContent(i0,0);
  if(iCat==1)
    for(int i0 = -1; i0 < 1000; i0++) if(iHMMSSTVBF[fQCDId]->GetBinContent(i0) < 0) iHMMSSTVBF[fQCDId]->SetBinContent(i0,0);
  // //bin by bin
  for(int i0 = 0; i0 < iHMM[fQCDId]->GetNbinsX()+1; i0++) if(iHMM[fQCDId]->GetBinContent(i0) == 0) iHMM[fQCDId]->SetBinError(i0,1.0);
  for(int i0 = 0; i0 < iHMM[0]->GetNbinsX()+1; i0++) 
    if(iHMM[0]->GetBinContent(i0) == 0) 
      iHMM[0]->SetBinError(i0,1.0);

  double lQCDSF=1;
  double lQCDVBFSF=1;
  lQCDSF = iHMMSST[fQCDId]->Integral(-1,1000)/iHMMSSL[fQCDId]->Integral(-1,1000);
  if(iCat==1)
    lQCDVBFSF = iHMMSSTVBF[fQCDId]->Integral(-1,1000)/iHMMSSL[fQCDId]->Integral(-1,1000);

  TH1F *ltest = (TH1F*) iHMMOS[fQCDId]->Clone("OSTmp"); 
  //iHMM[fQCDId] =ltest;
  iHMM[fQCDId] =iHMMOS[fQCDId];
  iHMM[fQCDId]->Scale(lQCDSF*lQCDVBFSF);
 
  double lWSF =1;
  if(iCat==1)
    {
      TH1F *vbf  = draw(iVar,iCat,iTree[fWId],fWId,"*(jdeta>3.5 && mjj>500 && njetingap==0)*(q_1*q_2 < 0)*(!(byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," Same Sign Tight");
      cout << vbf->Integral() << endl;
      cout << iHMMSSTVBF[fWId]->Integral() << endl;
      cout << iHMMOS[fWId]->Integral() << endl;
      cout << iHMMSSL[fWId]->Integral() << endl;
      cout << iHMMSST[fWId]->Integral() << endl;
      cout << iHMM[fWId]->Integral() << endl;
      //lWSF = (vbf->Integral()+iHMMSSTVBF[fWId]->Integral())/(iHMMOS[fWId]->Integral()+iHMMSSL[fWId]->Integral());
      lWSF = (vbf->Integral())/(iHMMOS[fWId]->Integral());
      double lWSFIso = iHHVBF[fWId]->Integral()/iHMMOS[fWId]->Integral(); 
      TH1F *lw = (TH1F*) iHMMOS[fWId]->Clone("OSTmp");  
      iHMM[fWId] = lw;
      iHMM[fWId]->Scale(lWSF*lWSFIso);
    }
  else if(iCat==2 || iCat==3)
    {
      TH1F *llw = (TH1F*) iHMMOS[fWId]->Clone("OSTmp");
      double wwscale = iHMM[fWId]->Integral(-1,1000)/iHMMOS[fWId]->Integral(-1,1000);
      iHMM[fWId] = llw;
      iHMM[fWId]->Scale(wwscale);
    }

  if(iCat==1)
    {
    double ltt = iHMM[1]->Integral(-1,1000)/iHHVBF[1]->Integral(-1,1000);  
    TH1F *htt = (TH1F*) iHHVBF[1]->Clone("OSTmp");
    iHMM[1]=htt;
    iHMM[1]->Scale(ltt);
    }

    if(iCat==1)
    {
      //TH1F* high   = draw(iVar+"high",iCat,iTree[0],0,"*(jdeta>2.0 && mjj>200 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," M",1);
      //TH1F* low   = draw(iVar+"low" ,iCat,iTree[0],0,"*(jdeta>2.0 && mjj>200 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," K",2);
      TH1F* high   = draw("m_vishigh",iCat,iTree[0],0,"*(jdeta>2.0 && mjj>200 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," M",1);
      TH1F* low   = draw("m_vislow" ,iCat,iTree[0],0,"*(jdeta>2.0 && mjj>200 && njetingap==0)*(q_1*q_2 < 0)*((byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0  && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0))"," K",2);
      double ltt = iHMM[0]->Integral(-1,1000)/iHHVBF[0]->Integral(-1,1000);  
      double ltthigh = iHHigh[0]->Integral(-1,1000)/high->Integral(-1,1000);
      double lttlow = iHLow[0]->Integral(-1,1000)/low->Integral(-1,1000);
      TH1F *htt = (TH1F*) iHHVBF[0]->Clone("OSTmp");
      iHMM[0]=htt;
      iHMM[0]->Scale(ltt);
      iHHigh[0]=high;
      iHHigh[0]->Scale(ltthigh);
      iHLow[0]=low;
      iHLow[0]->Scale(lttlow);
    }
  
 
  //if(iCat==1 || iCat==4 || iCat==3)
  if(iCat==1)
    {
      double lZL = 1.0e-7;
      TH1F *hzl = (TH1F*) iHMMOS[3]->Clone("OSTmp");
      iHMM[3]=hzl;
      iHMM[3]->Scale(lZL/iHMM[3]->Integral());
      TH1F *hzj = (TH1F*) iHMMOS[30]->Clone("OSTmp");
      double lZJ = iHMM[30]->Integral()/iHMMOS[30]->Integral();
      iHMM[30]=hzj;
      iHMM[30]->Scale(lZL/iHMM[30]->Integral());
    }

  
  if(iCat==1)
    {
      TH1F *hvv = (TH1F*) iHHVBF[4]->Clone("OSTmp");
      double lvv = iHMM[4]->Integral(-1,1000)/iHHVBF[4]->Integral(-1,1000);
      iHMM[4]=hvv;
      iHMM[4]->Scale(lvv);
    }
  cout << "HERE2?" << endl << endl;
  cout << "QCD Normalization " << lQCDSF << endl;
  cout << "QCD normalization " << lQCDVBFSF << endl;
  cout << "WJets Normaliztion" << lWSF << endl;
  
  if(iUseScale) 
    {
      iHHigh[fQCDId]=iHMM[fQCDId];
      iHLow[fQCDId] = iHMM[fQCDId];
      iHHigh[fWId]= iHMM[fWId];
      iHLow[fWId] = iHMM[fWId];
    }
  Double_t error;
 
  cout << "====> " << fString[fQCDId] << " -- " << iHMM[fQCDId]->IntegralAndError(0,iHMM[fQCDId]->GetNbinsX(),error) << endl;
  cout << " " << error <<  endl; 
  //cout << "====> " << fString[fQCDId] << " -- " << iHMM[fQCDId]->Integral() << endl;
  //cout << "====> " << fString[fWId] << " -- " << iHMM[fWId]->Integral() << endl;  
  cout << "====> " << fString[fWId] << " -- " << iHMM[fWId]->IntegralAndError(0,iHMM[fWId]->GetNbinsX(),error) << endl;
  cout << " " << error <<  endl; 
  cout << "====> " << fString[0] << " -- " << iHMM[0]->IntegralAndError(0,iHMM[0]->GetNbinsX(),error) << endl;
  cout << " " << error <<  endl; 
  cout << "====> " << fString[3] << " -- " << iHMM[3]->IntegralAndError(0,iHMM[3]->GetNbinsX(),error) << endl;
  cout << " " << error <<  endl; 
  //cout << "====> " << fString[3] << " -- " << iHMM[3]->Integral() << endl; 
  cout << "====> " << fString[4] << " -- " << iHMM[4]->Integral() << endl; 
  cout << "====> " << fString[1] << " -- " << iHMM[1]->Integral() << endl; 
  cout << "====> " << fString[iN-1] << " -- " << iHMM[iN-1]->Integral() << endl; 
  if(iUseScale)
    {
      cout << "====> " << fString[6] << " -- " << iHMM[6]->Integral() << endl;
      cout << "====> " << fString[7] << " -- " << iHMM[7]->Integral() << endl;
    }
 
  if(iFile !=0) makeDataCard(iFile,iHMM,iHHigh,iHLow,iHTheory,iLTheory,iN,iDirName,fCard);
 
  //Draw the plot
  //draw(iHMM,iN,iVar+"A",iVar,6,iCat);
  if(!iUseScale)
    {
      if(iVar.compare("m_sv")==0 || iVar.compare("mvis")==0)
      	{
      	  for(int i0 = 0; i0 < iN; i0++) {
	    iHMM[i0]->Scale(1.0,"width");
	    //iHMMOS[i0]->Scale(1.0,"width");
	    //iHMMSST[i0]->Scale(1.0,"width");
	    //iHMMSSL[i0]->Scale(1.0,"width");
      	  }
	}
      //if(iCat<7)
      //	blind(iHMM[iN-1],100,150);
     
      //iHMM[0]->Scale(1.002555);
      //iHMM[1]->Scale(0.807180);
      //iHMM[2]->Scale(0.937000);
      //iHMM[3]->Scale(0.962700);
      //iHMM[4]->Scale(0.776722);
      //iHMM[5]->Scale(1.126000);
      //iHMM[6]->Scale(0.942424);
      //iHMM[7]->Scale(0.965586);
      //iHMM[8]->Scale(0.967813);
      draw(iHMM,iN,iVar+"A",iVar,6,iCat);
      //draw(iHMMOS ,iN,iVar+"B",iVar,6,iCat);
      //draw(iHMMSST,iN,iVar+"C",iVar,6,iCat);
      //draw(iHMMSSL,iN,iVar+"D",iVar,6,iCat);
    }
}

void plotTauTau(std::string iVar="m_vis",int iCat = 0,int makecard=0,TFile *outfile=0,std::string iCut="(pt_1 > 45 && pt_2 > 45 && iso_1 > 0.5 && iso_2 > 0.5 && jpt_1 > 50.0 && abs(jeta_1) < 3.0 ",float iLumi=12000) { 
  SetStyle();  
  gStyle->SetLineStyleString(11,"20 10");
  // defining the common canvas, axes pad styles
  //SetStyle(); 
  loadfMap();
  fQCDId = 5;
  fWId =2;
  TH1::SetDefaultSumw2(1);
  int lN = 10;
  int ptrew[55]={0};
  if(makecard) 
    {
      lN = 56;   //44
    }

  std::string lName  = "/data/blue/Bacon/029a/tautaunew/";
   
  TTree **lTree  = new TTree*[lN]; 
  TH1F**lHMM = new TH1F*[lN]; 
  TH1F**lHMMOS   = new TH1F*[lN];
  TH1F**lHMMSST  = new TH1F*[lN];
  TH1F**lHMMSSL  = new TH1F*[lN]; 
  TH1F**lHHH     = new TH1F*[lN]; 
  TH1F**lHHL     = new TH1F*[lN];
  TH1F**lHHVBF   = new TH1F*[lN];
  TH1F**lHHVBFST = new TH1F*[lN];
  TH1F**lHTheory = new TH1F*[lN];
  TH1F**lLTheory = new TH1F*[lN];

  fString = new std::string[lN]; fWeights = new std::string[lN]; fColor = new int[lN]; fCard = new std::string[lN];
  //lName = "/data/blue/Bacon/029a/notauescaletautau/";
  lTree[0]  = load(lName+"emb_select.root");        fString[0] = "Z#rightarrow#tau#tau ";   fCard[0]="ZTT";       fColor[0] = kOrange-4;//kOrange-3; 
  //lName  = "/data/blue/Bacon/029a/tautau/ntuples/";
  lTree[1]  = load(lName+"ttbar-8TeV_select.root");     fString[1] = "t#bar{t}";      fCard[1]="TT";                 fColor[1] = kBlue-8;//kRed+4;
  lTree[2]  = load(lName+"wjets_select.root");          fString[2] = "W+Jets";        fCard[2]="W";                 fColor[2] = kRed+2;//kBlue-5;
  lTree[3]  = load(lName+"zll-mad_select.root");  fString[3] = "Z#rightarrow#tau#tau fakes";     fColor[3] = 603;  fCard[3]="ZL";
  lTree[4]  = load(lName+"ewk-8TeV_select.root");            fString[4] = "electroweak";    fCard[4]="VV";  fColor[4] = kRed+2;
  lTree[5]  = load(lName+"data_select.root");           fString[5] = "QCD";                 fCard[5]="QCD";           fColor[5] = kMagenta-10;//kBlue+3;
  //lName = "/data/blue/Bacon/029a/notauescaletautau/";
  lTree[6]  = load(lName+"higgsrw/Scalehtt_gf_sm_125_select.root");         fString[6] = "H(125)->#tau#tau";     fCard[6]="ggH125";   fColor[6] = 0; ptrew[6]=1;
  lTree[7]  = load(lName+"htt_vbf_sm_125_select.root");         fString[7] = "Higgs ";     fCard[7]="qqH125";   fColor[7] = 0;
  lTree[8]  = load(lName+"htt_vtth_sm_125_select.root");         fString[8] = "Higgs ";     fCard[8]="VH125";    fColor[8] = 0;
  //lName  = "/data/blue/Bacon/029a/tautau/ntuples/";
  lTree[lN-1]  = load(lName+"data_select.root");       fString[lN-1] = "observed"; fColor[lN-1] = kBlack;  fCard[lN-1]="data_obs";
  if(makecard)
    {
      //lName = "/data/blue/Bacon/029a/notauescaletautau/";
      lTree[9]  = load(lName+"higgsrw/Scalehtt_gf_sm_110_select.root");         fString[9] = "Higgs ";     fCard[9]="ggH110";  ptrew[9]=1;
      lTree[10]  = load(lName+"htt_vbf_sm_110_select.root");         fString[10] = "Higgs ";     fCard[10]="qqH110";  
      lTree[11]  = load(lName+"htt_vtth_sm_110_select.root");         fString[11] = "Higgs ";     fCard[11]="VH110"; 
      lTree[12]  = load(lName+"higgsrw/Scalehtt_gf_sm_115_select.root");         fString[12] = "Higgs ";     fCard[12]="ggH115";   ptrew[12]=1;
      lTree[13]  = load(lName+"htt_vbf_sm_115_select.root");         fString[13] = "Higgs ";     fCard[13]="qqH115";  
      lTree[14]  = load(lName+"htt_vtth_sm_115_select.root");         fString[14] = "Higgs ";     fCard[14]="VH115";
      lTree[15]  = load(lName+"higgsrw/Scalehtt_gf_sm_120_select.root");         fString[15] = "Higgs ";     fCard[15]="ggH120";   ptrew[15]=1;
      lTree[16]  = load(lName+"htt_vbf_sm_120_select.root");         fString[16] = "Higgs ";     fCard[16]="qqH120";  
      lTree[17]  = load(lName+"htt_vtth_sm_120_select.root");         fString[17] = "Higgs ";     fCard[17]="VH120"; 
      lTree[18]  = load(lName+"higgsrw/Scalehtt_gf_sm_130_select.root");         fString[18] = "Higgs ";     fCard[18]="ggH130";   ptrew[18]=1;
      lTree[19]  = load(lName+"htt_vbf_sm_130_select.root");         fString[19] = "Higgs ";     fCard[19]="qqH130";  
      lTree[20]  = load(lName+"htt_vtth_sm_130_select.root");         fString[20] = "Higgs ";     fCard[20]="VH130";
      lTree[21]  = load(lName+"higgsrw/Scalehtt_gf_sm_135_select.root");         fString[21] = "Higgs ";     fCard[21]="ggH135";  ptrew[21]=1;
      lTree[22]  = load(lName+"htt_vbf_sm_135_select.root");         fString[22] = "Higgs ";     fCard[22]="qqH135";  
      lTree[23]  = load(lName+"htt_vtth_sm_135_select.root");         fString[23] = "Higgs ";     fCard[23]="VH135";
      lTree[24]  = load(lName+"higgsrw/Scalehtt_gf_sm_140_select.root");         fString[24] = "Higgs ";     fCard[24]="ggH140";   ptrew[24]=1;
      lTree[25]  = load(lName+"htt_vbf_sm_140_select.root");         fString[25] = "Higgs ";     fCard[25]="qqH140";  
      lTree[26]  = load(lName+"htt_vtth_sm_140_select.root");         fString[26] = "Higgs ";     fCard[26]="VH140"; 
      lTree[27]  = load(lName+"higgsrw/Scalehtt_gf_sm_145_select.root");         fString[27] = "Higgs ";     fCard[27]="ggH145";  ptrew[27]=1;
      lTree[28]  = load(lName+"htt_vbf_sm_145_select.root");         fString[28] = "Higgs ";     fCard[28]="qqH145";  
      lTree[29]  = load(lName+"htt_vtth_sm_145_select.root");         fString[29] = "Higgs ";     fCard[29]="VH145"; 
      //lName  = "/data/blue/Bacon/029a/tautau/ntuples/";
      lTree[30]  = load(lName+"zll-mad_select.root");  fString[30] = "Z#rightarrow#tau#tau fakes";     fColor[30] = 603;  fCard[30]="ZJ";
      //lName = "/data/blue/Bacon/029a/notauescaletautau/";
      lTree[31]  = load(lName+"higgsrw/Scalehtt_gf_sm_105_select.root");         fString[31] = "Higgs ";     fCard[31]="ggH105";  ptrew[31]=1;
      lTree[32]  = load(lName+"htt_vbf_sm_105_select.root");         fString[32] = "Higgs ";     fCard[32]="qqH105";  
      lTree[33]  = load(lName+"htt_vtth_sm_105_select.root");         fString[33] = "Higgs ";     fCard[33]="VH105";
      lTree[34]  = load(lName+"higgsrw/Scalehtt_gf_sm_100_select.root");         fString[34] = "Higgs ";     fCard[34]="ggH100";  ptrew[34]=1;
      lTree[35]  = load(lName+"htt_vbf_sm_100_select.root");         fString[35] = "Higgs ";     fCard[35]="qqH100";  
      lTree[36]  = load(lName+"htt_vtth_sm_100_select.root");         fString[36] = "Higgs ";     fCard[36]="VH100";
      lTree[37]  = load(lName+"higgsrw/Scalehtt_gf_sm_95_select.root");         fString[37] = "Higgs ";     fCard[37]="ggH95";  ptrew[37]=1;
      lTree[38]  = load(lName+"htt_vbf_sm_95_select.root");         fString[38] = "Higgs ";     fCard[38]="qqH95";  
      lTree[39]  = load(lName+"htt_vtth_sm_95_select.root");         fString[39] = "Higgs ";     fCard[39]="VH95";
      lTree[40]  = load(lName+"higgsrw/Scalehtt_gf_sm_90_select.root");         fString[40] = "Higgs ";     fCard[40]="ggH90";  ptrew[40]=1;
      lTree[41]  = load(lName+"htt_vbf_sm_90_select.root");         fString[41] = "Higgs ";     fCard[41]="qqH90";  
      lTree[42]  = load(lName+"htt_vtth_sm_90_select.root");         fString[42] = "Higgs ";     fCard[42]="VH90";
      lTree[43]  = load(lName+"higgsrw/Scalehtt_gf_sm_150_select.root");         fString[43] = "Higgs ";     fCard[43]="ggH150";  ptrew[43]=1;
      lTree[44]  = load(lName+"htt_vbf_sm_150_select.root");         fString[44] = "Higgs ";     fCard[44]="qqH150";  
      lTree[45]  = load(lName+"htt_vtth_sm_150_select.root");         fString[45] = "Higgs ";     fCard[45]="VH150";
      lTree[46]  = load(lName+"higgsrw/Scalehtt_gf_sm_155_select.root");         fString[46] = "Higgs ";     fCard[46]="ggH155";  ptrew[46]=1;
      lTree[47]  = load(lName+"htt_vbf_sm_155_select.root");         fString[47] = "Higgs ";     fCard[47]="qqH155";  
      lTree[48]  = load(lName+"htt_vtth_sm_155_select.root");         fString[48] = "Higgs ";     fCard[48]="VH155";
      lTree[49]  = load(lName+"higgsrw/Scalehtt_gf_sm_160_select.root");         fString[49] = "Higgs ";     fCard[49]="ggH160";  ptrew[49]=1;
      lTree[50]  = load(lName+"htt_vbf_sm_160_select.root");         fString[50] = "Higgs ";     fCard[50]="qqH160";  
      lTree[51]  = load(lName+"htt_vtth_sm_160_select.root");         fString[51] = "Higgs ";     fCard[51]="VH160";
      lTree[52]  = load(lName+"htt_gf_ww_sm_125_select.root");         fString[52] = "Higgs ";     fCard[52]="ggH_hww125";
      lTree[53]  = load(lName+"htt_vbf_ww_sm_125_select.root");         fString[53] = "Higgs ";     fCard[53]="qqH_hww125";
      lTree[54]  = load(lName+"htt_vtth_ww_sm_125_select.root");         fString[54] = "Higgs ";     fCard[54]="VH_hww125";
    }
  std::stringstream lLumi; lLumi << iLumi;

  std::string dir;
  if(iCat==0) 
    {
      //iCut=iCut+" && pthmva>100 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<2.5 && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<2.5)"; 
      iCut=iCut+" &&  byCombinedIsolationDeltaBetaCorrRaw3Hits_2<2.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<2.0)"; 
      dir="tauTau_inclusive";
    }
  else if(iCat==1) {
    iCut=iCut+"  && pthmva>100 && nbtag==0 && jpt_2 > 30 && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<10.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<10.0)";
    dir="tauTau_vbf";
  }
  else if(iCat==2)
    {
      iCut=iCut+"  && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<10.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<10.0 && njets>0 && nbtag==0 && pthmva>100 && pthmva<170 && (!(jpt_2 > 30 && mjj > 500 && jdeta > 3.5 && njetingap==0)))"; dir="tauTau_1jet_high_mediumhiggs";
    }
  else if(iCat==3)  
    {
      iCut=iCut+" && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<10.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<10.0  && nbtag==0 && njets>0 && pthmva>170 && (!(jpt_2 > 30 && mjj > 500 && jdeta > 3.5 && njetingap==0)))"; dir="tauTau_1jet_high_highhiggs";
    }
  else if(iCat==4) {
    iCut=iCut+" && nbtag==0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<10.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<10.0)";
    dir="tauTau_nobtag";
  }
  else if(iCat==5) {
    iCut=iCut+" && nbtag==1 && njets<2)";
    dir="tauTau_btag";
  }
  else
    {
      iCut=iCut+" && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<4.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<4.0)";
      dir="Inclusive";
    }

  TTree *ldy  = new TTree;
  ldy = load(lName+"ztt-mad_select.root");

  //fWeights[0]="(pt_1 > 45 && pt_2 > 45 &&  abs(eta_1) <2.1 && abs(eta_2)<2.1)*(weight*decaymodeweight)";
  //fWeights[3]="(pt_1 > 45 && pt_2 > 45 && abs(eta_1) <2.1  && abs(eta_2)<2.1)*weight*decaymodeweight";
  fWeights[0]="(pt_1 > 45 && pt_2 > 45 &&  abs(eta_1) <2.1 && abs(eta_2)<2.1)*(weight*decaymodeweight*(genmatch==3))";
  fWeights[3]="(pt_1 > 45 && pt_2 > 45 && abs(eta_1) <2.1  && abs(eta_2)<2.1)*(abs(genlid_1)>=15 && abs(genlid_2)>=15)*weight*decaymodeweight*(genmatch==3)";
  TH1F* lMC = new TH1F;
  TH1F* lEmb = new TH1F;
  lMC=draw("m_sv",iCat,ldy,3,"*(q_1*q_2 < 0)*(byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0 && passAntiEleNewWPMVA3_2>0.5) "," Main");
  lEmb=draw("m_sv",iCat,lTree[0],0,"*(q_1*q_2 < 0)*(byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0 && passAntiEleNewWPMVA3_2>0.5) "," Main");
  cout << lMC->Integral() << endl;
  cout << lEmb->Integral() << endl;
  double lfc = lMC->Integral()/lEmb->Integral();
  cout << lfc << endl;
  std::stringstream ifc; ifc << lfc;
  delete lMC;
  delete lEmb;
  delete ldy;
  
  for(int i0 = 0; i0 < lN;   i0++)  fWeights[i0]   = iCut;
  fWeights[0] += "*"+ifc.str();
  for(int i0 = 0; i0 < lN-1; i0++) if(i0 != fQCDId) fWeights[i0]  += "*(weight*decaymodeweight)*"+lLumi.str();
  
  fWeights[52] += "*(1.0/1.972e09)";
  fWeights[53] += "*(1.0/1.972e09)";
  fWeights[54] += "*(1.0/1.972e09)";
  //for(int i0 = 3; i0 < lN-1; i0++) 
  // if(i0 != fQCDId) fWeights[i0]  += "*1.0";
  //1.972e09
  //Z Scale Factors
  //fWeights[0]  += "*1.01*(embkinweight/genDiTauMassVsGenDiTauPt)*(genmatch==3)";
  fWeights[0]  += "*1.01*(genmatch==3)";
  //fWeights[0]  += "*1.01*(genmatch==3)";
  fWeights[1]  += "*1.11*0.96"; //0.605
  if(iCat==1)
    fWeights[2]  += "*1.03";
  else if(iCat==2 || iCat==0)
    fWeights[2]  += "*0.87";
  else if(iCat==3) 
    fWeights[2] += "*0.76";
  fWeights[5]  += "*1.0";
  if(makecard)
    fWeights[3]  += "*1.01*(genmatch==1 && !(abs(genlid_1)>=15 && abs(genlid_2)>=15))";
  else
    fWeights[3]  += "*1.01*(genmatch != 3)";
 
  if(!makecard)
    {
      fWeights[6]  += "*1.2337";
      fWeights[7]  += "*0.0997";
      fWeights[8]  += "*0.0772";
    }
  if(makecard) 
    {
      fWeights[30]  += "*(genmatch!=1 && genmatch!=3)";
      cout << "Hey Hey" << endl;
    }
     
  //TCanvas *lC0 = new TCanvas("A","A",700,700);
  drawSpec(lTree,lHMM,lHMMOS,lHMMSST,lHMMSSL,lN,iVar,iCat,outfile,dir,makecard,ptrew,lHHH,lHHL,lHHVBF,lHHVBFST,lHTheory,lLTheory); 
  //lC0->cd();
  delete [] lTree;
  delete [] lHMM;
  delete [] lHMMOS;
  delete [] lHMMSST;
  delete [] lHMMSSL;
  delete [] lHHH;
  delete [] lHHL;
  delete [] lHHVBF;
  delete [] lHHVBFST;
  delete [] lHTheory;
  delete [] lLTheory;
}
 
