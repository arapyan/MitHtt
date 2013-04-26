#include "TH2F.h"

std::string fChanId =  "";
std::string *fString  = 0;
std::string *fWeights = 0;
int *fColor           = 0;
int  fId              = 0;
int  fYId             = 0;

TTree * load(std::string iName) { 
  TFile *lFile = new TFile(iName.c_str());
  lFile->cd();
  TTree *lTree  = (TTree*) lFile->FindObjectAny("Events");
  if(lTree == 0) lTree = (TTree*) lFile->FindObjectAny("Flat"); 
  if(lTree == 0) lTree = (TTree*) lFile->FindObjectAny("TestTree"); 
  if(lTree == 0) lTree = (TTree*) lFile->FindObjectAny("TauCheck"); 
  return lTree;
}
void scale(TH1F** iH,double iInt,const int iN) {
  double lTotInt = iH[0]->Integral();
  cout <<  " Total " << lTotInt << " -- " << iInt/lTotInt << " - " << iN << endl;
  double lMax = 0;
  for(int i0 = 0; i0 < iN; i0++) {
    if(i0 != iN-1) iH[i0]->Scale(double(iInt/lTotInt));
    for(int i1 = 0; i1 < iH[i0]->GetNbinsX()+1; i1++) {
      double pMax = iH[i0]->GetBinContent(i1);
      if(lMax < pMax) lMax = pMax;
    }
    cout << " ==> " << i0 << " -- " << iH[i0]->Integral() << endl;
  }
  for(int i0 = 0; i0 < iN; i0++) iH[i0]->GetYaxis()->SetRangeUser(0.001,lMax*1.7);
}
void drawDifference(TH1* iH0,TH1 *iH1,TH1 *iHH=0,TH1 *iHL=0) {
  std::string lName = std::string(iH0->GetName());
  TH1F *lHDiff  = new TH1F((lName+"Diff").c_str(),(lName+"Diff").c_str(),iH0->GetNbinsX(),iH0->GetXaxis()->GetXmin(),iH0->GetXaxis()->GetXmax());
  TH1F *lHDiffH = new TH1F((lName+"DiffH").c_str(),(lName+"DiffH").c_str(),iH0->GetNbinsX(),iH0->GetXaxis()->GetXmin(),iH0->GetXaxis()->GetXmax()); lHDiffH->SetLineWidth(1); lHDiffH->SetLineColor(kRed);
  TH1F *lHDiffL = new TH1F((lName+"DiffL").c_str(),(lName+"DiffL").c_str(),iH0->GetNbinsX(),iH0->GetXaxis()->GetXmin(),iH0->GetXaxis()->GetXmax()); lHDiffL->SetLineWidth(1); lHDiffL->SetLineColor(kBlue);
  lHDiff->SetFillColor(kViolet); lHDiff->SetFillStyle(1001); lHDiff->SetLineWidth(1);
  TH1F *lXHDiff1 = new TH1F((lName+"XDiff1").c_str(),(lName+"XDiff1").c_str(),iH0->GetNbinsX(),iH0->GetXaxis()->GetXmin(),iH0->GetXaxis()->GetXmax());
  TH1F *lXHDiff2 = new TH1F((lName+"XDiff2").c_str(),(lName+"XDiff2").c_str(),iH0->GetNbinsX(),iH0->GetXaxis()->GetXmin(),iH0->GetXaxis()->GetXmax());
  int i1 = 0;
  lXHDiff1->SetLineWidth(2); lXHDiff1->SetLineColor(kRed);
  lXHDiff2->SetLineWidth(2); lXHDiff2->SetLineColor(kRed);

  lXHDiff1->GetYaxis()->SetTitle("Ratio");
  lXHDiff1->GetYaxis()->SetRangeUser(0.95,1.1);
  lXHDiff1->GetYaxis()->SetTitleOffset(0.4);
  lXHDiff1->GetYaxis()->SetTitleSize(0.2);
  lXHDiff1->GetYaxis()->SetLabelSize(0.11);
  for(int i0 = 0; i0 < lHDiff->GetNbinsX()+1; i0++) {
    double lXCenter = lHDiff->GetBinCenter(i0);
    double lXVal     = iH0   ->GetBinContent(i0);
    double lXValH    = iHH   ->GetBinContent(i0);
    double lXValL    = iHL   ->GetBinContent(i0);
    lXHDiff1->SetBinContent(i0, 1.0);
    lXHDiff2->SetBinContent(i0, 1.0);
    while(iH1->GetBinCenter(i1) < lXCenter) {i1++;}
    if(iH1->GetBinContent(i0) > 0) lHDiff->SetBinContent(i0,lXVal      /(iH1->GetBinContent(i0)));
    if(iH1->GetBinContent(i0) > 0) lHDiff->SetBinError  (i0,sqrt(lXVal)/(iH1->GetBinContent(i0)));
    if(iH1->GetBinContent(i0) > 0) lHDiffL->SetBinContent(i0,lXValL/(iH1->GetBinContent(i0)));
    if(iH1->GetBinContent(i0) > 0) lHDiffH->SetBinContent(i0,lXValH/(iH1->GetBinContent(i0)));
    //if(iH1->GetBinContent(i0) > 0)  cout << "unc" << lXVal << " -- " << sqrt(lXVal)/(iH1->GetBinContent(i0)) << endl;
  }
  lHDiff->SetMarkerStyle(kFullCircle);
  //lHDiff->Draw("EP");
  lXHDiff1->Draw("hist");
  lXHDiff2->Draw("hist sames");
  lHDiffH ->Draw("hist sames");
  lHDiffL ->Draw("hist sames");
}
void clear(TH1F *iH) { 
  for(int i0 = 0; i0 < iH->GetNbinsX()+1; i0++) iH->SetBinContent(i0,0);
}
TH1F* draw(std::string iVar,TTree *iTree,int iId,std::string iCut,std::string iSId="B",TH1F *iH=0) {
  std::stringstream lName;   lName   <<  "Met"  << iId << iSId;
  std::stringstream lVar;
  lVar   << iVar << ">>+" << lName.str();
  std::string lCut   = fWeights[iId];
  lCut+=iCut;
  cout << "var : " << iVar << "Cut : " << lCut << endl;
  TH1F *lMet = iH;
  if(iH == 0) lMet = new TH1F(lName.str().c_str(),lName.str().c_str(),getNBins(iVar),fAxis);
  lMet->SetLineColor(fColor[iId]);
  lMet->SetLineWidth(1);
  lMet->GetXaxis()->SetTitle(getXAxis(iVar));
  lMet->GetYaxis()->SetTitle(getYAxis(iVar));
  lMet->Sumw2();
  TString lDraw(lVar.str().c_str());
  iTree->Draw(lDraw,lCut.c_str());
  lMet->SetMarkerStyle(kFullCircle);
  lMet->SetLineWidth(2);
  return lMet;
}
void makeDataCard(TFile *iFile,TH1F** iH,TH1F** iHH,TH1F** iHL,const int iN,std::string iDirName,std::string* iHistName) { 
  iFile->cd();
  TDirectory* lTD = iFile->mkdir(iDirName.c_str());
  iFile->cd(lTD->GetPath());
  for(int i0 = 0; i0 < iN; i0++) { 
    iH[i0]->SetFillStyle(1001);
    iH[i0]->SetLineWidth(1); iH[i0]->SetLineColor(kBlack);
    std::string pName = iHistName[i0];
    iH[i0]->SetName (pName.c_str()); 
    iH[i0]->SetTitle(pName.c_str());
    iH[i0]->Write();
    if(iHH[i0] == 0 || iHL[i0] == 0) continue;
    iHH[i0]->SetName ((pName+std::string("_CMS_scale_t_"+fChanId+"Up")  ).c_str());
    iHH[i0]->SetTitle((pName+std::string("_CMS_scale_t_"+fChanId+"Up")  ).c_str());
    iHL[i0]->SetName ((pName+std::string("_CMS_scale_t_"+fChanId+"Down")).c_str());
    iHL[i0]->SetTitle((pName+std::string("_CMS_scale_t_"+fChanId+"Down")).c_str());
    iHL[i0]->Write();
    iHH[i0]->Write();
  }
  if(iDirName.find("boost_low") != std::string::npos) {
    TH1F* lHUp   = iH[fQCDId]->Clone(("QCD_CMS_htt_QCDShape_"+iDirName+"_"+fChanId+"Up").c_str());
    TH1F* lHDown = iH[fQCDId]->Clone(("QCD_CMS_htt_QCDShape_"+iDirName+"_"+fChanId+"Down").c_str());
    lHUp  ->SetName (("QCD_CMS_htt_QCDShape_"+iDirName+"_"+fChanId+"Up").c_str());
    lHDown->SetName (("QCD_CMS_htt_QCDShape_"+iDirName+"_"+fChanId+"Up").c_str());
    lHUp  ->SetTitle(("QCD_CMS_htt_QCDShape_"+iDirName+"_"+fChanId+"Up").c_str());
    lHDown->SetTitle(("QCD_CMS_htt_QCDShape_"+iDirName+"_"+fChanId+"Up").c_str());
    for(int i0 = 0; i0 < lHUp->GetNbinsX()+1; i0++) { 
      if(lHUp->GetXaxis()->GetBinCenter(i0) < 80) { 
	lHUp  ->SetBinContent(i0,1.2*lHUp  ->GetBinContent(i0));
	lHDown->SetBinContent(i0,0.8*lHDown->GetBinContent(i0));
      }
    }
    lHUp  ->Write();
    lHDown->Write();
  }
}
void fit(TH1F *iH0,TH1F *iH1, TH1F *iH2) { 
  RooRealVar lM("m","m",0,200);
  RooRealVar lN1("m1","m1",0,200000);
  RooRealVar lN2("m2","m2",0,200000);
  RooRealVar lScale("scale","scale",0,2.);
  iH1->Add(iH2,-1.);

  lN1.setVal(iH1->Integral()); lN1.setConstant(kTRUE);
  lN2.setVal(iH2->Integral()); lN2.setConstant(kTRUE);
  cout << "====> " << iH2->Integral() << endl;
  RooFormulaVar lNT("mT","mT","@0*@1",RooArgList(lN2,lScale));
  RooDataHist *pH0  =  new RooDataHist("Data","Data" ,RooArgList(lM),iH0);
  RooDataHist *pH1  =   new RooDataHist("MC0" ,"MC0" ,RooArgList(lM),iH1);
  RooDataHist *pH2  =   new RooDataHist("MC1" ,"MC1" ,RooArgList(lM),iH2);
  RooHistPdf lHP1("H1Pdf","H1Pdf",lM,*pH1,0);
  RooHistPdf lHP2("H2Pdf","H2Pdf",lM,*pH2,0);
  RooExtendPdf *pH1E = new RooExtendPdf ("EMC1","EMC1",lHP1,lN1); lN1.setConstant(kTRUE);
  RooExtendPdf *pH2E = new RooExtendPdf ("EMC2","EMC2",lHP2,lNT); 
  RooAddPdf    *pHA  = new RooAddPdf    ("Add" ,"Add",RooArgList(*pH1E,*pH2E));
  pHA->fitTo(*pH0,RooFit::Save(kTRUE),RooFit::SumW2Error(kFALSE),RooFit::Strategy(1),RooFit::Extended());
  TCanvas *lCrap = new TCanvas("A","A",800,600); lCrap->cd();
  RooPlot *lFrame1 = lM.frame(RooFit::Title("XXX")) ;                                                                                                                                               
  pH0 ->plotOn(lFrame1);//,RooFit::Cut("cat==cat::cat_0"));                                                                                                                                    
  pH2 ->plotOn(lFrame1);
  lFrame1->Draw();
  cin.get();
}
void draw(TH1F** iH,const int iN,std::string iName,std::string iFName,int iHiggs=-1) { 
  TCanvas *lC0 = new TCanvas(("XX"+iName).c_str(),("XX"+iName).c_str(),fId,fYId,400,400); if(fId + 400 > 1100) fYId = (fYId + 400) % 1000; fId = (fId + 400) % 1100; 
  TLegend *lL = new TLegend(0.6,0.65,0.9,0.9); lL->SetFillColor(0); lL->SetBorderSize(0);
  for(int i0 = 0; i0 < iN; i0++) { 
    iH[i0]->SetFillStyle(1001);
    if(i0 != iHiggs) iH[i0]->SetFillColor(fColor[i0]);
    for(int i1 = i0+1; i1 < iN-1; i1++)  {if(i1 != iHiggs) iH[i0]->Add(iH[i1]);} 
    if(i0 != iN-1 && iH[i0]->Integral() > 0.0001) lL->AddEntry(iH[i0],fString[i0].c_str(),"f");
    iH[i0]->SetLineWidth(1); iH[i0]->SetLineColor(kBlack);
  }
  if(iHiggs != -1) iH[iHiggs]->SetLineStyle(kDashed);
  if(iHiggs != -1) iH[iHiggs]->SetFillStyle(0);
  if(iHiggs != -1) iH[iHiggs]->SetLineWidth(1);
  if(iHiggs != -1 && iHiggs != 0) iH[iHiggs]->Add(iH[0]);
  if(iHiggs == 0) iH[iHiggs]->Add(iH[1]);
  //fit(iH[iN-1],iH[0],iH[iN-3]);
  
  double lMax = iH[iN-1]->GetMaximum();
  iH[iN-1]->SetMarkerSize(1.);
  //scale(iH,iH[iN-1]->Integral(),iN);
  iH[iN-1]->SetMarkerStyle(kFullCircle); iH[iN-1]->GetYaxis()->SetRangeUser(0.1,1.5*lMax);
  //lC0->Divide(1,2); lC0->cd();  lC0->cd(1)->SetPad(0,0.2,1.0,1.0); gPad->SetLeftMargin(0.2) ;
  iH[iN-1]->SetTitle("");
  iH[iN-1]->Draw("EP");
  double lVal = iH[0]->Integral();
  cout << "=====> Scale " << iH[iN-1]->Integral()/lVal << endl; 
  for(int i0 = 0; i0 < iN-1; i0++) {
    iH[i0]->Scale(iH[iN-1]->Integral()/lVal);
    //if(i0 == iHiggs) continue;
    cout << "================> " << fString[i0] << " -- " << iH[i0]->Integral() << " -- " << iH[i0]->GetMean() << " -- " << iH[i0]->GetRMS() << endl;
    iH[i0]->Draw("hist sames");
  }
  iH[iN-1]  ->Draw("EP sames");
  lL->Draw();
  //lC0->cd(2)->SetPad(0,0,1.0,0.2); gPad->SetLeftMargin(0.2) ;
  //drawDifference(iH[iN-1],iH[0]);//,iH[6],iH[7]);
  lC0->SaveAs((iName+".png").c_str());
}
TLegend * legend(int iN,TH1** iH) {
  TLegend *lL = new TLegend(0.20,0.7,0.35,0.9); lL->SetFillColor(0); lL->SetBorderSize(0);
  for(int i0 = 0; i0 < iN; i0++) {
    //if(i0 != iN-1 && i0 == 0 & i0 == 1) continue;
    if(iH[i0]->Integral() > -0.1 && i0  == 0) lL->AddEntry(iH[i0],fString[i0].c_str(),"l");
    if(iH[i0]->Integral() > -0.1 && i0  != 0) lL->AddEntry(iH[i0],fString[i0].c_str(),"l");
  }
  return lL;
}
TH2F* draw2D(std::string iVar,std::string iVar1,TTree *iTree,int iId,std::string iCut,TH2F *iH = 0) {
  std::stringstream lName;   lName   <<  "Met"  << iId;// << iVar;
  std::stringstream lVar;
  lVar   << iVar1 << ":" << iVar << ">>+" << lName.str();
  //std::string lCut   = "(pt_1 > 10 && pt_2 > 10 && abs(eta_1) < 2.5 && abs(eta_2) < 2.5 && q_1*q_2 < 0)*"+fWeights[iId]+iCut;
  std::string lCut   = fWeights[iId]+iCut;
  cout << "Cut : " << lCut << endl;
  TH2F *lMet = iH;
  if(iH == 0) lMet = new TH2F(lName.str().c_str(),lName.str().c_str(),getNBins(iVar),getXMin(iVar),getXMax(iVar),getNBins(iVar1),getXMin(iVar1),getXMax(iVar1));
  lMet->SetLineColor(fColor[iId]);
  lMet->SetLineWidth(1);
  lMet->GetXaxis()->SetTitle(getXAxis(iVar));
  lMet->GetYaxis()->SetTitle(getYAxis(iVar1));//"Events/5 GeV");     
  iTree->Draw(lVar.str()  .c_str()     ,(lCut    ).c_str());
  lMet->SetMarkerStyle(kFullCircle);
  lMet->SetLineWidth(2);
  return lMet;
}
TH1F* drawRaw(std::string iVar,TTree *iTree,int iId,std::string iCut,TH1F *iH = 0) {
  std::stringstream lName;   lName   <<  "Met"  << iId;
  std::stringstream lVar;
  lVar   << iVar << ">>+" << lName.str();
  std::string lCut   = iCut+"*"+fWeights[iId]+iCut;
  TH1F *lMet = iH;
  if(iH == 0) lMet = new TH1F(lName.str().c_str(),lName.str().c_str(),getNBins(iVar),getXMin(iVar),getXMax(iVar));
  lMet->SetLineColor(fColor[iId]);
  lMet->SetLineWidth(1);
  lMet->GetXaxis()->SetTitle(getXAxis(iVar));//"m_{T} (GeV/c^{2})");//m_{sv}[GeV/c^{2}]");                                                                                                                         
  lMet->GetYaxis()->SetTitle(getYAxis(iVar));//"Events/5 GeV");                                                                                                                                                    
  iTree->Draw(lVar.str()  .c_str()     ,(lCut    ).c_str());
  return lMet;
}
void drawBasic(TH1F** iH,const int iN,std::string iName) { 
  TCanvas *lC0 = new TCanvas(("XX"+iName).c_str(),("XX"+iName).c_str(),fId,fYId,600,600); if(fId + 500 > 1100) fYId = (fYId + 400) % 1000; fId = (fId + 500) % 1100; 
  TLegend *lL = new TLegend(0.6,0.65,0.9,0.9); lL->SetFillColor(0); lL->SetBorderSize(0);
  for(int i0 = 0; i0 < iN; i0++) { 
    std::stringstream pSS; pSS << i0;
    cout << "===> " << iH[i0] << " -- " << iH[i0]->Integral() << endl;
    if(iH[i0]->Integral() > 0.1) lL->AddEntry(iH[i0],fString[i0].c_str(),"f");
  }
  scale(iH,iH[iN-1]->Integral(),iN);
  //iH[iN-1]->SetMarkerStyle(kFullCircle); //iH[iN-1]->GetYaxis()->SetRangeUser(0.1,1500);
  //lC0->Divide(1,2); lC0->cd();  lC0->cd(1)->SetPad(0,0.2,1.0,1.0); gPad->SetLeftMargin(0.2) ;
  iH[iN-1]->Draw("hist");
  for(int i0 = 0; i0 < iN-1; i0++) {
    iH[i0]->Draw("hist sames");
    lL->Draw();
  }
  //iH[iN-1]->Draw("hist sames");
  //lC0->cd(2)->SetPad(0,0,1.0,0.2); gPad->SetLeftMargin(0.2) ;
  //drawDifference(iH[iN-1],iH[0]);//,iH[6],iH[7]);
  //lC0->SaveAs((iFName+".png").c_str());
}
