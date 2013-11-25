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

void drawSpec(TTree **iTree,TH1F **iHMM,TH1F **iHMMSS,TH1F **iHMMMT,TH1F **iHMMSST,int iN,std::string iVar,int iCat,TFile * iFile = 0,std::string iDirName="",int iUseScale=0,int *ptrew=0,TH1F** iHHigh=0,TH1F** iHLow=0,TH1F** iHHVBF=0,TH1F** iHSST=0,TH1F** iHMMFB=0,TH1F **iHSSTFB=0,TH1F **iHMMSSFB=0) { 
  cout << "=====>  -- " << fQCDId << endl;
  for(int i0 = 0; i0 < iN; i0++) {
    iHMM   [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(mtMVA_1<30)*(iso_1<0.1)"," Main");
    iHSST  [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 > 0)*(mtMVA_1<30)*(iso_1<0.1)","Same Sign Iso");
    iHMMMT [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(mtMVA_1>70)*(iso_1<0.1)","mT sideband");	
    iHMMSST[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 > 0)*(mtMVA_1>70)*(iso_1<0.1)"," Same Sign mT Sideband");
    iHMMSS [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 > 0)*(mtMVA_1<30)*(iso_1>0.2 && iso_1<0.5)","Same Sign Loose Iso");
    //notes:  iHMM->W shape, iHMMMT/(iHSST+(iHMM))->W norm, iHSST->QCD norm, qcd shape-->iHSST, iHMMMST->W norm in ss
    //3 numbers-->inc..  QCD yield in OS, iHMMSS yield for tau iso/relax
    iHMMFB[i0] = 0;
    cout << endl << endl << "Here?" << endl;
    if(ptrew[i0] && iUseScale)
      {
	iHMMFB[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(mtMVA_1<30)*(iso_1<0.1)"," Main",0,0,1);
	iHSSTFB  [i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 > 0)*(mtMVA_1<30)*(iso_1<0.1)","Same Sign Iso",0,0,1);
	iHMMSSFB[i0]   = draw(iVar,iCat,iTree[i0],i0,"*(q_1*q_2 > 0)*(mtMVA_1<30)*(iso_1>0.2 && iso_1<0.5)","Same Sign Loose Iso",0,0,1);
      }
    
    cout << "====> " << fString[i0] << " -- " << iHMM[i0]->Integral() << endl;
    
    if(iFile == 0    ) continue;
    iHHigh[i0]   = 0;
    iHLow[i0]    = 0;
    //iHTheory[i0]   = 0;
    //iLTheory[i0]    = 0;
    if(!iUseScale || i0==iN-1 || i0==fQCDId || i0==fWId || i0==1) continue;
    if(i0==3)
      {
	iHHigh[i0]   = draw(iVar+"*1.02",iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(mtMVA_1<30)*(iso_1<0.1)"," Main");
	iHLow[i0]   = draw(iVar+"*0.98",iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(mtMVA_1<30)*(iso_1<0.1)"," Main");
      }
    else
      {
	iHHigh[i0]   = draw(iVar+"high",iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(mtMVA_1<30)*(iso_1<0.1)"," Main",1);
	iHLow[i0]   = draw(iVar+"low",iCat,iTree[i0],i0,"*(q_1*q_2 < 0)*(mtMVA_1<30)*(iso_1<0.1)"," Main",2);
      }
  }
  
  //MT scale factor
  TH1F* lMMMT = (TH1F *) iHMMMT[0]->Clone("OSMT"); clear(lMMMT);
  TH1F* lMMMST = (TH1F*) iHMMSST[0]->Clone("SSMT"); clear(lMMMST);
  for(int i =0; i<iN-1; i++)
    {
      if(i==fQCDId || i==fWId) continue;
      if(i>5 && i<iN-1 && i!=48) continue;
      lMMMT->Add(iHMMMT[i]);
      lMMMST->Add(iHMMSST[i]);      
    }
  iHMMMT[iN-1]->Add(lMMMT,-1);
  for(int i0 = -1; i0 < 1000; i0++) if(iHMMMT[iN-1]->GetBinContent(i0) < 0) iHMMMT[iN-1]->SetBinContent(i0,0);
  double lWSF=1.0;
  double lSWSF=1.0;
  iHMMSST[iN-1]->Add(lMMMST,-1);
  for(int i0 = -1; i0 < 1000; i0++) if(iHMMSST[iN-1]->GetBinContent(i0) < 0) iHMMSST[iN-1]->SetBinContent(i0,0);
  lSWSF = (iHSST[fWId]->Integral())/(iHMMSST[fWId]->Integral());
  iHSST[fWId]->Scale(lSWSF*iHMMSST[iN-1]->Integral()/iHSST[fWId]->Integral());
  
  lWSF = iHMM[fWId]->Integral()/iHMMMT[fWId]->Integral();
  iHMM[fWId]->Scale(lWSF*iHMMMT[iN-1]->Integral()/iHMM[fWId]->Integral());
  if(iUseScale)
    iHMMFB[fWId]->Scale(lWSF*iHMMMT[iN-1]->Integral()/iHMMFB[fWId]->Integral());
  
  //QCD norm and Shape
  TH1F *lSST = (TH1F*) iHSST[0]->Clone("OSTmp"); clear(lSST);
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 == fQCDId) continue;
    if(i0 >5 && i0<iN-1 && i0!=48) continue;
    lSST->Add(iHSST[i0]);
  }
  iHSST[fQCDId]->Add(lSST,-1); 
  for(int i0 = -1; i0 < 1000; i0++) if(iHSST[fQCDId]->GetBinContent(i0) < 0) iHSST[fQCDId]->SetBinContent(i0,0);

  //Compute QCD Shape
  if(iCat==8 || iCat==10)
    {
      TH1F *lSST1 = (TH1F*) iHSST[0]->Clone("SSTTmp"); clear(lSST1);
      lSST1->Add(iHSST[0]);
      lSST1->Add(iHSST[3]);
      lSST1->Add(iHSST[1]);
      lSST1->Add(iHSST[4]);
      if(iUseScale)
	lSST1->Add(iHSST[48]);
      lSST1->Add(iHSST[fWId]);
      iHSST[iN-1]->Add(lSST1,-1); 
      for(int i0 = -1; i0 < 1000; i0++) if(iHSST[iN-1]->GetBinContent(i0) < 0) iHSST[iN-1]->SetBinContent(i0,0);
      iHMM[fQCDId]=iHSST[iN-1];
      iHMM[fQCDId]->Scale(1.06*iHSST[fQCDId]->Integral()/iHMM[fQCDId]->Integral());
      if(iUseScale)
	{
	  TH1F *lSSTFB = (TH1F*) iHSSTFB[0]->Clone("SSTTmp"); clear(lSSTFB);
	  lSSTFB->Add(iHSSTFB[0]);
	  lSSTFB->Add(iHSSTFB[3]);
	  lSSTFB->Add(iHSSTFB[1]);
	  lSSTFB->Add(iHSSTFB[4]);
	  if(iUseScale)
	    lSSTFB->Add(iHSSTFB[48]);
	  lSSTFB->Add(iHSSTFB[fWId]);
	  iHSSTFB[iN-1]->Add(lSSTFB,-1); 
	  for(int i0 = -1; i0 < 1000; i0++) if(iHSSTFB[iN-1]->GetBinContent(i0) < 0) iHSSTFB[iN-1]->SetBinContent(i0,0);
	  iHMMFB[fQCDId]=iHSSTFB[iN-1];
	  iHMMFB[fQCDId]->Scale(1.06*iHSST[fQCDId]->Integral()/iHMMFB[fQCDId]->Integral());
	}
    }
  else
    {
      iHMM[fQCDId]=iHMMSS[iN-1];
      iHMM[fQCDId]->Scale(1.06*iHSST[fQCDId]->Integral()/iHMM[fQCDId]->Integral());
      if(iUseScale)
	{
	  iHMMFB[fQCDId]=iHMMSSFB[iN-1];
	  iHMMFB[fQCDId]->Scale(1.06*iHSST[fQCDId]->Integral()/iHMMFB[fQCDId]->Integral());
	}
    }

  // //bin by bin
  for(int i0 = 0; i0 < iHMM[fQCDId]->GetNbinsX()+1; i0++) if(iHMM[fQCDId]->GetBinContent(i0) == 0) iHMM[fQCDId]->SetBinError(i0,1.0);
  for(int i0 = 0; i0 < iHMM[0]->GetNbinsX()+1; i0++) 
    if(iHMM[0]->GetBinContent(i0) == 0) 
      iHMM[0]->SetBinError(i0,1.0);
   
  Double_t error;
 
  cout << "====> " << fString[fQCDId] << " -- " << iHMM[fQCDId]->IntegralAndError(0,iHMM[fQCDId]->GetNbinsX(),error) << endl;
  cout << " " << error <<  endl;
  if(iCat==10)
    {
      cout << "====> " << fString[fQCDId] << " -- " << iHMMSS [fQCDId]->IntegralAndError(0,iHMMSS[fQCDId]->GetNbinsX(),error) << endl;
      cout << " " << error <<  endl;
    } 
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
      if(iCat==8)
	{
	  for(int i0 = 0; i0 < iHMM[fQCDId]->GetNbinsX()+1; i0++) { 
	    if(iHMM[fQCDId]->GetXaxis()->GetBinCenter(i0) < 50) { 
	      iHMM[fQCDId]->SetBinContent(i0,1.1*iHMM[fQCDId]->GetBinContent(i0));
	    }
	  }
	  for(int i0 = 0; i0 < iHMMFB[fQCDId]->GetNbinsX()+1; i0++) { 
	    if(iHMMFB[fQCDId]->GetXaxis()->GetBinCenter(i0) <50) { 
	      iHMMFB[fQCDId]->SetBinContent(i0,1.1*iHMMFB[fQCDId]->GetBinContent(i0));
	    }
	  }
	  iHMM[1]->Scale((iHMM[1]->Integral()-0.015*iHMM[0]->Integral())/iHMM[1]->Integral());
	}
    }

  if(iFile !=0) makeDataCard(iFile,iHMM,iHHigh,iHLow,iHMMFB,iHMM,iN,iDirName,fCard);

  if(!iUseScale)
    {
      if(iVar.compare("m_sv")==0 || iVar.compare("m_vis")==0)
      	{
      	  for(int i0 = 0; i0 < iN; i0++) {
      	    iHMM[i0]->Scale(1.0,"width");
	    //iHMMMT[i0]->Scale(1.0,"width");
	    //iHMMSST[i0]->Scale(1.0,"width");
	    //iHMMSSL[i0]->Scale(1.0,"width");
      	  }
      	}
      // if(iCat<10)
// 	blind(iHMM[iN-1],100,150);
     
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
      //draw(iHMMMT ,iN,iVar+"B",iVar,6,iCat);
      //draw(iHMMOS ,iN,iVar+"B",iVar,6,iCat);
      //draw(iHMMSST,iN,iVar+"C",iVar,6,iCat);
      //draw(iHMMSSL,iN,iVar+"D",iVar,6,iCat);
    }
}

void plotLTauMSSM(std::string iVar="m_vis",int iCat = 0,int makecard=0,TFile *outfile=0,std::string iCut="(pt_1 > 45 && pt_2 > 45 && iso_1 > 0.5 && iso_2 > 0.5 && jpt_1 > 30.0 && abs(jeta_1) < 3.0 ",float iLumi=12000) { 
  SetStyle();  
  gStyle->SetLineStyleString(11,"20 10");
  // defining the common canvas, axes pad styles
  //SetStyle(); 
  loadfMap();
  fQCDId = 5;
  fWId =2;
  TH1::SetDefaultSumw2(1);
  int ptrew[52]={0};
  int lN = 10;
  if(makecard) 
    {
      lN = 53;   //44
    }

   std::string lName  = "/data/blue/Bacon/029a/mutau/ntuples/";
   
   TTree **lTree  = new TTree*[lN]; 
   TH1F**lHMM = new TH1F*[lN]; 
   TH1F**lHMMMT   = new TH1F*[lN];
   TH1F**lHMMMST  = new TH1F*[lN];
   TH1F**lHMMSSL  = new TH1F*[lN]; 
   TH1F**lHHH     = new TH1F*[lN]; 
   TH1F**lHHL     = new TH1F*[lN];
   TH1F**lHHVBF     = new TH1F*[lN];
   TH1F**lHHSST     = new TH1F*[lN];
   //TH1F**lHTheory = new TH1F*[lN];
   //TH1F**lLTheory = new TH1F*[lN];
   TH1F**lHMMOSFB   = new TH1F*[lN];
   TH1F**lHMMSSTFB  = new TH1F*[lN];
   TH1F**lHMMSSLFB  = new TH1F*[lN]; 
   
  fString = new std::string[lN]; fWeights = new std::string[lN]; fColor = new int[lN]; fCard = new std::string[lN];
  lTree[0]  = load(lName+"emb_select.root");        fString[0] = "Z#rightarrow#tau#tau ";   fCard[0]="ZTT";       fColor[0] = kOrange-4; ptrew[0]=1;
  lTree[1]  = load(lName+"ttbar-8TeV_select.root");     fString[1] = "t#bar{t}";      fCard[1]="TT";                 fColor[1] = kBlue-8;  ptrew[1]=1;//kRed+4;
  lTree[2]  = load(lName+"wjets_select.root");          fString[2] = "W+Jets";        fCard[2]="W";                 fColor[2] = kRed+2; ptrew[2]=1;//kBlue-5;
  lTree[3]  = load(lName+"zmm_select.root");  fString[3] = "Z#rightarrow#tau#tau fakes";     fColor[3] = 603;  fCard[3]="ZL"; ptrew[3]=1;
  lTree[4]  = load(lName+"ewk-8TeV_select.root");            fString[4] = "electroweak";    fCard[4]="VV";  fColor[4] = kRed+2; ptrew[4]=1;
  lTree[5]  = load(lName+"data_select.root");           fString[5] = "QCD";                 fCard[5]="QCD";           fColor[5] = kMagenta-10;  ptrew[5]=1;//kBlue+3; 
  lTree[6]  = load(lName+"htt_bbh_mssm_160_select.root");         fString[6] = "5X H(125)->#tau#tau";     fCard[6]="bbH160";   fColor[6] = 0;
  lTree[7]  = load(lName+"htt_ggh_mssm_160_select.root");         fString[7] = "Higgs ";     fCard[7]="ggH160";   fColor[7] = 0;
  lTree[lN-1]  = load(lName+"data_select.root");       fString[lN-1] = "observed"; fColor[lN-1] = kBlack;  fCard[lN-1]="data_obs";
  if(makecard)
    {
      lTree[8]  = load(lName+"htt_bbh_mssm_80_select.root");         fString[8] = "Higgs ";     fCard[8]="bbH80";  
      lTree[9]  = load(lName+"htt_ggh_mssm_80_select.root");         fString[9] = "Higgs ";     fCard[9]="ggH80";  
      lTree[10]  = load(lName+"htt_bbh_mssm_90_select.root");         fString[10] = "Higgs ";     fCard[10]="bbH90";  
      lTree[11]  = load(lName+"htt_ggh_mssm_90_select.root");         fString[11] = "Higgs ";     fCard[11]="ggH90";
      lTree[12]  = load(lName+"htt_bbh_mssm_100_select.root");         fString[12] = "Higgs ";     fCard[12]="bbH100";  
      lTree[13]  = load(lName+"htt_ggh_mssm_100_select.root");         fString[13] = "Higgs ";     fCard[13]="ggH100"; 
      lTree[14]  = load(lName+"htt_bbh_mssm_110_select.root");         fString[14] = "Higgs ";     fCard[14]="bbH110";  
      lTree[15]  = load(lName+"htt_ggh_mssm_110_select.root");         fString[15] = "Higgs ";     fCard[15]="ggH110"; 
      lTree[16]  = load(lName+"htt_bbh_mssm_120_select.root");         fString[16] = "Higgs ";     fCard[16]="bbH120";  
      lTree[17]  = load(lName+"htt_ggh_mssm_120_select.root");         fString[17] = "Higgs ";     fCard[17]="ggH120";
      lTree[18]  = load(lName+"htt_bbh_mssm_130_select.root");         fString[18] = "Higgs ";     fCard[18]="bbH130";  
      lTree[19]  = load(lName+"htt_ggh_mssm_130_select.root");         fString[19] = "Higgs ";     fCard[19]="ggH130";
      lTree[20]  = load(lName+"htt_bbh_mssm_140_select.root");         fString[20] = "Higgs ";     fCard[20]="bbH140";  
      lTree[21]  = load(lName+"htt_ggh_mssm_140_select.root");         fString[21] = "Higgs ";     fCard[21]="ggH140";
      lTree[22]  = load(lName+"htt_bbh_mssm_180_select.root");         fString[22] = "Higgs ";     fCard[22]="bbH180";  
      lTree[23]  = load(lName+"htt_ggh_mssm_180_select.root");         fString[23] = "Higgs ";     fCard[23]="ggH180";
      lTree[24]  = load(lName+"htt_bbh_mssm_200_select.root");         fString[24] = "Higgs ";     fCard[24]="bbH200";  
      lTree[25]  = load(lName+"htt_ggh_mssm_200_select.root");         fString[25] = "Higgs ";     fCard[25]="ggH200";
      lTree[26]  = load(lName+"htt_bbh_mssm_250_select.root");         fString[26] = "Higgs ";     fCard[26]="bbH250";  
      lTree[27]  = load(lName+"htt_ggh_mssm_250_select.root");         fString[27] = "Higgs ";     fCard[27]="ggH250";
      lTree[28]  = load(lName+"htt_bbh_mssm_300_select.root");         fString[28] = "Higgs ";     fCard[28]="bbH300";  
      lTree[29]  = load(lName+"htt_ggh_mssm_300_select.root");         fString[29] = "Higgs ";     fCard[29]="ggH300";
      lTree[30]  = load(lName+"htt_bbh_mssm_350_select.root");         fString[30] = "Higgs ";     fCard[30]="bbH350";  
      lTree[31]  = load(lName+"htt_ggh_mssm_350_select.root");         fString[31]= "Higgs ";     fCard[31]="ggH350";
      lTree[32]  = load(lName+"htt_bbh_mssm_400_select.root");         fString[32] = "Higgs ";     fCard[32]="bbH400";  
      lTree[33]  = load(lName+"htt_ggh_mssm_400_select.root");         fString[33] = "Higgs ";     fCard[33]="ggH400";
      lTree[34]  = load(lName+"htt_bbh_mssm_450_select.root");         fString[34] = "Higgs ";     fCard[34]="bbH450";  
      lTree[35]  = load(lName+"htt_ggh_mssm_450_select.root");         fString[35] = "Higgs ";     fCard[35]="ggH450";
      lTree[36]  = load(lName+"htt_bbh_mssm_500_select.root");         fString[36] = "Higgs ";     fCard[36]="bbH500";  
      lTree[37]  = load(lName+"htt_ggh_mssm_500_select.root");         fString[37] = "Higgs ";     fCard[37]="ggH500";
      lTree[38]  = load(lName+"htt_bbh_mssm_600_select.root");         fString[38] = "Higgs ";     fCard[38]="bbH600";  
      lTree[39]  = load(lName+"htt_ggh_mssm_600_select.root");         fString[39] = "Higgs ";     fCard[39]="ggH600";
      lTree[40]  = load(lName+"htt_bbh_mssm_700_select.root");         fString[40] = "Higgs ";     fCard[40]="bbH700";  
      lTree[41]  = load(lName+"htt_ggh_mssm_700_select.root");         fString[41] = "Higgs ";     fCard[41]="ggH700";
      lTree[42]  = load(lName+"htt_bbh_mssm_800_select.root");         fString[42] = "Higgs ";     fCard[42]="bbH800";  
      lTree[43]  = load(lName+"htt_ggh_mssm_800_select.root");         fString[43] = "Higgs ";     fCard[43]="ggH800";
      lTree[44]  = load(lName+"htt_bbh_mssm_900_select.root");         fString[44] = "Higgs ";     fCard[44]="bbH900";  
      lTree[45]  = load(lName+"htt_ggh_mssm_900_select.root");         fString[45] = "Higgs ";     fCard[45]="ggH900";
      lTree[46]  = load(lName+"htt_bbh_mssm_1000_select.root");         fString[46] = "Higgs ";     fCard[46]="bbH1000";  
      lTree[47]  = load(lName+"htt_ggh_mssm_1000_select.root");         fString[47] = "Higgs ";     fCard[47]="ggH1000";
      lTree[48]  = load(lName+"zll-mad_select.root");  fString[48] = "Z#rightarrow#tau#tau fakes";     fColor[48] = 603;  fCard[48]="ZJ"; ptrew[48]=1;
      lTree[49]  = load(lName+"htt_gf_sm_125_select.root");         fString[49] = "H(125)->#tau#tau";     fCard[49]="ggH_SM125";   fColor[49] = 0;
      lTree[50]  = load(lName+"htt_vbf_sm_125_select.root");         fString[50] = "Higgs ";     fCard[50]="qqH_SM125";   fColor[50] = 0;
      lTree[51]  = load(lName+"htt_vtth_sm_125_select.root");         fString[51] = "Higgs ";     fCard[51]="VH_SM125";    fColor[51] = 0;
    }

  std::stringstream lLumi; lLumi << iLumi;

  std::string dir;

  if(iCat==8) {
    iCut=iCut+" && nbtag==0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5)";
    dir="muTau_nobtag";
  }
  else if(iCat==9) {
     iCut=iCut+" && nbtag>=1 && jpt_2<30 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5)";
     dir="muTau_btag";
  }
  else {
    iCut=iCut+" && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5)";
    dir="muTau_inclusive";
  }

  TTree *ldy  = new TTree;
  ldy = load(lName+"ztt-mad_select.root");

  fWeights[0]="(pt_1 > 20 && pt_2 > 20 &&  abs(eta_1) <2.1 && abs(eta_2)<2.3 &&  abs(dZ_2)<0.2)*weight*decaymodeweight*emblid*(genmatch==3)";
  fWeights[3]="(pt_1 > 20 && pt_2 > 20 && abs(eta_1) <2.1  && abs(eta_2)<2.3 &&  abs(dZ_2)<0.2)*weight*decaymodeweight*(genmatch==3)";
  //fWeights[3]="(pt_1 > 20 && pt_2 > 20 && abs(eta_1) <2.1  && abs(eta_2)<2.3 &&  abs(dZ_2)<0.2)*weight*decaymodeweight*(genmatch==3)";
 
  TH1F* lMC = new TH1F;
  TH1F* lEmb = new TH1F;
  lMC=draw("m_sv",iCat,ldy,3,"*(q_1*q_2 < 0)*(iso_1<0.1 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5) "," Main");
  lEmb=draw("m_sv",iCat,lTree[0],0,"*(q_1*q_2 < 0)*(iso_1<0.1 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5) "," Main");
  double lfc = lMC->Integral()/lEmb->Integral();
  cout << lfc << endl;
  std::stringstream ifc; ifc << lfc;
  delete lMC;
  delete lEmb;
  
  for(int i0 = 0; i0 < lN;   i0++)  fWeights[i0]   = iCut;
  fWeights[0] +=  "*"+ifc.str(); //0.982 comes from having the non isolated trigger in the current sample for DY jets
  for(int i0 = 0; i0 < lN-1; i0++) if(i0 != fQCDId) fWeights[i0]  += "*(weight*decaymodeweight/1.972e09)*"+lLumi.str();
  
  //Z Scale Factors
  fWeights[0]  += "*1.01*emblid*(genmatch==3)";
  fWeights[1]  += "*1.11*0.96"; //0.988 comes from the non isolated trigger being used
  if(iCat==9)
    fWeights[1]  += "*1.0";
  else if(iCat==8)
    fWeights[1]  += "*1.0";
  fWeights[4]  += "*1.0"; //0.93 comes from the non isolated trigger
  if(makecard)
    fWeights[3]  += "*1.01*(abs(genlid_1)<15 && abs(genlid_2)<15 && genmatch==1)";
  else
    fWeights[3]  += "*1.01*(abs(genlid_1)<15 && abs(genlid_2)<15)";
 
  //if(!makecard)
  // {
  fWeights[49]  += "*1.2179";
  fWeights[50]  += "*0.0997";
  fWeights[51]  += "*0.0789"; //check
  // }
  if(makecard) 
    {
      fWeights[48]  += "*(genmatch!=3 && genmatch!=1)";
      cout << "Hey Hey" << endl;
    }
     
  //TCanvas *lC0 = new TCanvas("A","A",700,700);
  drawSpec(lTree,lHMM,lHMMSSL,lHMMMT,lHMMMST,lN,iVar,iCat,outfile,dir,makecard,ptrew,lHHH,lHHL,lHHVBF,lHHSST,lHMMOSFB,lHMMSSTFB,lHMMSSLFB);
 
  //lC0->cd();
  delete [] lTree;
  delete [] lHMM;
  delete [] lHMMMT;
  delete [] lHMMMST;
  delete [] lHMMSSL;
  delete [] lHHH;
  delete [] lHHL;
  delete [] lHHSST;
  //delete [] lHTheory;
  //delete [] lLTheory;
}
