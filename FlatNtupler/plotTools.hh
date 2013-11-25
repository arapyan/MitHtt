#include "TH2F.h"

int fQCDId,fWId;
std::string *fString  = 0;
std::string *fCard    = 0;
std::string *fWeights = 0;
int *fColor           = 0;
int  fId              = 0;
int  fYId             = 0;

TTree * load(std::string iName) { 
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree  = (TTree*) lFile->FindObjectAny("Events");
  if(lTree == 0) lTree = (TTree*) lFile->FindObjectAny("Flat"); 
  if(lTree == 0) lTree = (TTree*) lFile->FindObjectAny("TestTree"); 
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
float maximum(TH1F* h, bool LOG=false){
  if(LOG){
    if(h->GetMaximum()>1000){ return 1000.*TMath::Nint(500*h->GetMaximum()/1000.); }
    if(h->GetMaximum()>  10){ return   10.*TMath::Nint( 50*h->GetMaximum()/  10.); }
    return 50*h->GetMaximum();
  }
  else{
    if(h->GetMaximum()>  12){ return 10.*TMath::Nint((1.3*h->GetMaximum()/10.)); }
    if(h->GetMaximum()> 1.2){ return TMath::Nint((1.6*h->GetMaximum())); }
    return 1.6*h->GetMaximum();
  }
}
void binbybin(TH1F* iH,double low,double high) {
  for(int i = 0; i < iH->GetNbinsX()+1; i++) {
    double mass = iH->GetBinCenter(i);
    if(mass <= high && mass >= low)
      {
	//iH->SetBinContent(i,0);
	iH->SetBinError(i,iH->GetBinContent(i)*0.1001);
      }
  }
}
void blind(TH1F* iH,double low,double high) {
  for(int i = 0; i < iH->GetNbinsX()+1; i++) {
    double mass = iH->GetBinCenter(i);
    if(mass <= high && mass >= low)
      {
	iH->SetBinContent(i,0);
	iH->SetBinError(i,0);
      }
  }
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

  lHDiff->GetYaxis()->SetTitle("Ratio");
  lHDiff->GetYaxis()->SetRangeUser(0.8,1.2);
  lHDiff->GetYaxis()->SetTitleOffset(0.4);
  lHDiff->GetYaxis()->SetTitleSize(0.2);
  lHDiff->GetYaxis()->SetLabelSize(0.11);
  for(int i0 = 0; i0 < lHDiff->GetNbinsX()+1; i0++) {
    double lXCenter = lHDiff->GetBinCenter(i0);
    double lXVal     = iH0   ->GetBinContent(i0);
    //double lXValH    = iHH   ->GetBinContent(i0);
    //double lXValL    = iHL   ->GetBinContent(i0);
    lXHDiff1->SetBinContent(i0, 1.0);
    lXHDiff2->SetBinContent(i0, 1.0);
    while(iH1->GetBinCenter(i1) < lXCenter) {i1++;}
    if(iH1->GetBinContent(i0) > 0) lHDiff->SetBinContent(i0,lXVal      /(iH1->GetBinContent(i0)));
    if(iH1->GetBinContent(i0) > 0) lHDiff->SetBinError  (i0,sqrt(lXVal)/(iH1->GetBinContent(i0)));
    //if(iH1->GetBinContent(i0) > 0) lHDiffL->SetBinContent(i0,lXValL/(iH1->GetBinContent(i0)));
    //if(iH1->GetBinContent(i0) > 0) lHDiffH->SetBinContent(i0,lXValH/(iH1->GetBinContent(i0)));
    //if(iH1->GetBinContent(i0) > 0)  cout << "unc" << lXVal << " -- " << sqrt(lXVal)/(iH1->GetBinContent(i0)) << endl;
  }
  lHDiff->SetMarkerStyle(kFullCircle);
  lHDiff->Draw("EP");
  lXHDiff1->Draw("hist sames");
  lXHDiff2->Draw("hist sames");
  //lHDiffH ->Draw("hist sames");
  //lHDiffL ->Draw("hist sames");
}
void clear(TH1F *iH) { 
  for(int i0 = 0; i0 < iH->GetNbinsX()+1; i0++) iH->SetBinContent(i0,0);
}
void makeDataCard(TFile *iFile,TH1F** iH,TH1F** iHH,TH1F** iHL,TH1F** iHTheory,TH1F** iLTheory,const int iN,std::string iDirName,std::string* iHistName,TH1F** TEffHigh=0,TH1F** TEffLow=0) { 
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
    cout << "===> " << i0 << endl;
    if(iHH[i0] == 0 || iHL[i0] == 0) continue;
    if(i0==3)
      {
	iHH[i0]->SetName ("ZL_CMS_htt_ZLScale_tautau_8TeVUp");
	iHH[i0]->SetTitle("ZL_CMS_htt_ZLScale_tautau_8TeVUp");
	iHL[i0]->SetName ("ZL_CMS_htt_ZLScale_tautau_8TeVDown");
	iHL[i0]->SetTitle("ZL_CMS_htt_ZLScale_tautau_8TeVDown");
      }
    else
      {
	iHH[i0]->SetName ((pName+std::string("_CMS_scale_t_tautau_8TeVUp")  ).c_str());
	iHH[i0]->SetTitle((pName+std::string("_CMS_scale_t_tautau_8TeVUp")  ).c_str());
	iHL[i0]->SetName ((pName+std::string("_CMS_scale_t_tautau_8TeVDown")).c_str());
	iHL[i0]->SetTitle((pName+std::string("_CMS_scale_t_tautau_8TeVDown")).c_str());
      }
    iHL[i0]->Write();
    iHH[i0]->Write();
    
//     if(TEffHigh[i0] != 0 && TEffLow[i0] !=0)
//       {
// 	TEffHigh[i0]->SetName ((pName+std::string("_CMS_eff_t_mssmHigh_tautau_8TeVUp")  ).c_str());
// 	TEffHigh[i0]->SetTitle((pName+std::string("_CMS_eff_t_mssmHigh_tautau_8TeVUp")  ).c_str());
// 	TEffLow[i0]->SetName ((pName+std::string("_CMS_eff_t_mssmHigh_tautau_8TeVDown")).c_str());
// 	TEffLow[i0]->SetTitle((pName+std::string("_CMS_eff_t_mssmHigh_tautau_8TeVDown")).c_str());
// 	TEffHigh[i0]->Write();
// 	TEffLow[i0]->Write();
//       }
   
   
    if(iHTheory[i0] == 0 || iLTheory[i0] == 0) continue;
    //if(iHTheory[i0] == 0) continue;
    //iHTheory[i0]->SetTitle((pName+"_fine_binning").c_str()); 
    //iHTheory[i0]->SetName ((pName+"_fine_binning").c_str()); 
    //iHTheory[i0]->Write();
    iHTheory[i0]->SetTitle((pName+"_QCDscale_ggH1inUp").c_str()); 
    iLTheory[i0]->SetTitle((pName+"_QCDscale_ggH1inDown").c_str()); 
    iHTheory[i0]->SetName ((pName+"_QCDscale_ggH1inUp").c_str()); 
    iLTheory[i0]->SetName ((pName+"_QCDscale_ggH1inDown").c_str()); 
    iHTheory[i0]->Write();
    iLTheory[i0]->Write();
  }
  //if(iDirName.find("muTau_1jet_medium") != std::string::npos || iDirName.find("muTau_1jet_high_lowhiggs") != std::string::npos || iDirName.find("muTau_vbf_loose") != std::string::npos) {
    if(iDirName.find("muTau_1jet_medium") != std::string::npos) {
    TH1F* lHUp   =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_1jet_medium_7TeVUp");
    TH1F* lHDown =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_1jet_medium_7TeVDown");
    lHUp  ->SetName ("QCD_CMS_htt_QCDShape_mutau_1jet_medium_7TeVUp");
    lHDown->SetName ("QCD_CMS_htt_QCDShape_mutau_1jet_medium_7TeVDown");
    lHUp  ->SetTitle("QCD_CMS_htt_QCDShape_mutau_1jet_medium_7TeVUp");
    lHDown->SetTitle("QCD_CMS_htt_QCDShape_mutau_1jet_medium_7TeVDown");
    for(int i0 = 0; i0 < lHUp->GetNbinsX()+1; i0++) { 
      if(lHUp->GetXaxis()->GetBinCenter(i0) < 50) { 
	lHUp  ->SetBinContent(i0,1.1*lHUp  ->GetBinContent(i0));
	lHDown->SetBinContent(i0,0.9*lHDown->GetBinContent(i0));
      }
    }
    lHUp  ->Write();
    lHDown->Write();
    }
    if(iDirName.find("muTau_1jet_high_lowhiggs") != std::string::npos) {
      TH1F* lHUp   =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_1jet_high_lowhiggs_7TeVUp");
    TH1F* lHDown =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_1jet_high_lowhiggs_7TeVDown");
    lHUp  ->SetName ("QCD_CMS_htt_QCDShape_mutau_1jet_high_lowhiggs_7TeVUp");
    lHDown->SetName ("QCD_CMS_htt_QCDShape_mutau_1jet_high_lowhiggs_7TeVDown");
    lHUp  ->SetTitle("QCD_CMS_htt_QCDShape_mutau_1jet_high_lowhiggs_7TeVUp");
    lHDown->SetTitle("QCD_CMS_htt_QCDShape_mutau_1jet_high_lowhiggs_7TeVDown");
    for(int i0 = 0; i0 < lHUp->GetNbinsX()+1; i0++) { 
      if(lHUp->GetXaxis()->GetBinCenter(i0) < 50) { 
	lHUp  ->SetBinContent(i0,1.1*lHUp  ->GetBinContent(i0));
	lHDown->SetBinContent(i0,0.9*lHDown->GetBinContent(i0));
      }
    }
    lHUp  ->Write();
    lHDown->Write();
    }
    if(iDirName.find("muTau_vbf_loose") != std::string::npos) {
      TH1F* lHUp   =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_etau_vbf_loose_8TeVUp");
      TH1F* lHDown =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_etau_vbf_loose_8TeVDown");
    lHUp  ->SetName ("QCD_CMS_htt_QCDShape_etau_vbf_loose_8TeVUp");
    lHDown->SetName ("QCD_CMS_htt_QCDShape_etau_vbf_loose_8TeVDown");
    lHUp  ->SetTitle("QCD_CMS_htt_QCDShape_etau_vbf_loose_8TeVUp");
    lHDown->SetTitle("QCD_CMS_htt_QCDShape_etau_vbf_loose_8TeVDown");
    for(int i0 = 0; i0 < lHUp->GetNbinsX()+1; i0++) { 
      if(lHUp->GetXaxis()->GetBinCenter(i0) < 40) { 
	lHUp  ->SetBinContent(i0,1.1*lHUp  ->GetBinContent(i0));
	lHDown->SetBinContent(i0,0.9*lHDown->GetBinContent(i0));
      }
    }
    lHUp  ->Write();
    lHDown->Write();
    }
    if(iDirName.find("muTau_vbf") != std::string::npos) {
      TH1F* lHUp   =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_vbf_7TeVUp");
      TH1F* lHDown =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_vbf_7TeVDown");
    lHUp  ->SetName ("QCD_CMS_htt_QCDShape_mutau_vbf_7TeVUp");
    lHDown->SetName ("QCD_CMS_htt_QCDShape_mutau_vbf_7TeVDown");
    lHUp  ->SetTitle("QCD_CMS_htt_QCDShape_mutau_vbf_7TeVUp");
    lHDown->SetTitle("QCD_CMS_htt_QCDShape_mutau_vbf_7TeVDown");
    for(int i0 = 0; i0 < lHUp->GetNbinsX()+1; i0++) { 
      if(lHUp->GetXaxis()->GetBinCenter(i0) < 40) { 
	lHUp  ->SetBinContent(i0,1.1*lHUp  ->GetBinContent(i0));
	lHDown->SetBinContent(i0,0.9*lHDown->GetBinContent(i0));
      }
    }
    lHUp  ->Write();
    lHDown->Write();
    }
  if(iDirName.find("btag") != std::string::npos) {
    TH1F* lHUp   =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_btag_8TeVUp");
    TH1F* lHDown =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_btag_8TeVDown");
    lHUp  ->SetName ("QCD_CMS_htt_QCDShape_mutau_btag_8TeVUp");
    lHDown->SetName ("QCD_CMS_htt_QCDShape_mutau_btag_8TeVDown");
    lHUp  ->SetTitle("QCD_CMS_htt_QCDShape_mutau_btag_8TeVUp");
    lHDown->SetTitle("QCD_CMS_htt_QCDShape_mutau_btag_8TeVDown");
    TH1F* lHUpfine   =  (TH1F *) iHTheory[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_btag_8TeVUp_fine_binning");
    TH1F* lHDownfine =  (TH1F *) iHTheory[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_btag_8TeVDown_finebinning");
    lHUpfine  ->SetName ("QCD_CMS_htt_QCDShape_mutau_btag_8TeVUp_fine_binning");
    lHDownfine->SetName ("QCD_CMS_htt_QCDShape_mutau_btag_8TeVDown_fine_binning");
    lHUpfine  ->SetTitle("QCD_CMS_htt_QCDShape_mutau_btag_8TeVUp_fine_binning");
    lHDownfine->SetTitle("QCD_CMS_htt_QCDShape_mutau_btag_8TeVDown_fine_binning");
    for(int i0 = 0; i0 < lHUpfine->GetNbinsX()+1; i0++) { 
      if(lHUpfine->GetXaxis()->GetBinCenter(i0) < 40) { 
	lHUpfine  ->SetBinContent(i0,1.1*lHUpfine  ->GetBinContent(i0));
	lHDownfine->SetBinContent(i0,0.9*lHDownfine->GetBinContent(i0));
      }
    }
    for(int i0 = 0; i0 < lHUp->GetNbinsX()+1; i0++) { 
      if(lHUp->GetXaxis()->GetBinCenter(i0) < 40) { 
	lHUp  ->SetBinContent(i0,1.1*lHUp  ->GetBinContent(i0));
	lHDown->SetBinContent(i0,0.9*lHDown->GetBinContent(i0));
      }
    }
    lHUp  ->Write();
    lHDown->Write();
    lHUpfine  ->Write();
    lHDownfine->Write();
  }
  if(iDirName.find("nobtag") != std::string::npos) {
    TH1F* lHUp   =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVUp");
    TH1F* lHDown =  (TH1F *) iH[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVDown");
    lHUp  ->SetName ("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVUp");
    lHDown->SetName ("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVDown");
    lHUp  ->SetTitle("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVUp");
    lHDown->SetTitle("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVDown");
    for(int i0 = 0; i0 < lHUp->GetNbinsX()+1; i0++) { 
      if(lHUp->GetXaxis()->GetBinCenter(i0) < 50) { 
	lHUp  ->SetBinContent(i0,1.1*lHUp  ->GetBinContent(i0));
	lHDown->SetBinContent(i0,0.9*lHDown->GetBinContent(i0));
      }
    }
    TH1F* lHUpfine   =  (TH1F *) iHTheory[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVUp_fine_binning");
    TH1F* lHDownfine =  (TH1F *) iHTheory[fQCDId]->Clone("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVDown_finebinning");
    lHUpfine  ->SetName ("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVUp_fine_binning");
    lHDownfine->SetName ("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVDown_fine_binning");
    lHUpfine  ->SetTitle("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVUp_fine_binning");
    lHDownfine->SetTitle("QCD_CMS_htt_QCDShape_mutau_nobtag_8TeVDown_fine_binning");
    for(int i0 = 0; i0 < lHUpfine->GetNbinsX()+1; i0++) { 
      if(lHUpfine->GetXaxis()->GetBinCenter(i0) < 50) { 
	lHUpfine  ->SetBinContent(i0,1.1*lHUpfine  ->GetBinContent(i0));
	lHDownfine->SetBinContent(i0,0.9*lHDownfine->GetBinContent(i0));
      }
    }
    lHUp  ->Write();
    lHDown->Write();
    lHUpfine  ->Write();
    lHDownfine->Write();
  }
}
TCanvas* MakeCanvasR(const char* name, const char *title, int dX, int dY)
{
  // Start with a canvas
  TCanvas *canvas = new TCanvas(name,title,0,0,dX,dY);
  // General overall stuff
  canvas->SetFillColor      (0);
  canvas->SetBorderMode     (0);
  canvas->SetBorderSize     (10);
  // Set margins to reasonable defaults
  canvas->SetLeftMargin     (0.18);
  canvas->SetRightMargin    (0.05);
  canvas->SetTopMargin      (0.08);
  canvas->SetBottomMargin   (0.15);
  // Setup a frame which makes sense
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);
  canvas->SetFrameFillStyle (0);
  canvas->SetFrameLineStyle (0);
  canvas->SetFrameBorderMode(0);
  canvas->SetFrameBorderSize(10);

  return canvas;
}
void cmsprelim(const char* channel, double lowX, double lowY)
{
  TPaveText* cmsprel  = new TPaveText(lowX, lowY, lowX+0.25, lowY+0.04, "NDC");
  cmsprel->SetBorderSize(   0 );
  cmsprel->SetFillStyle(    0 );
  cmsprel->SetTextAlign(   11 );
  cmsprel->SetTextSize ( 0.035);
  cmsprel->SetTextColor(    1 );
  cmsprel->SetTextFont (   62 );
  cmsprel->AddText("CMS Preliminary 2012, 19.4 fb^{-1}, #sqrt{s} = 8 TeV");
  cmsprel->Draw();

  // TPaveText* lumi     = new TPaveText(lowX+0.35, lowY, lowX+0.55, lowY+0.04, "NDC");
  // lumi->SetBorderSize(   0 );
  // lumi->SetFillStyle(    0 );
  // lumi->SetTextAlign(   11 );
  // lumi->SetTextSize ( 0.035);
  // lumi->SetTextColor(    1 );
  // lumi->SetTextFont (   42 );
  // lumi->AddText("4.6 fb^{-1}  #sqrt{s}=7 TeV");
  // lumi->Draw();

  // TPaveText* chan     = new TPaveText(lowX+0.68, lowY, lowX+0.73, lowY+0.04, "NDC");
  // chan->SetBorderSize(   0 );
  // chan->SetFillStyle(    0 );
  // chan->SetTextAlign(   11 );
  // chan->SetTextSize ( 0.035);
  // chan->SetTextColor(    1 );
  // chan->SetTextFont (   42 );
  // chan->AddText(channel);
  // chan->Draw();
}
void draw(TH1F** iH,const int iN,std::string iName,std::string iFName,int iHiggs=-1,int iCat=0) { 
  //TCanvas *lC0 = new TCanvas(("XX"+iName).c_str(),("XX"+iName).c_str(),fId,fYId,600,600); if(fId + 400 > 1100) fYId = (fYId + 400) % 1000; fId = (fId + 400) % 1100; 
  TCanvas *lC0 = MakeCanvasR("canv", "histograms", 600, 600);
  //TCanvas *lC0 = new TCanvas(("XX"+iName).c_str(),("XX"+iName).c_str(),fId,fYId,700,700); 
  //if(fId + 400 > 1100) fYId = (fYId + 400) % 1000; fId = (fId + 400) % 1100; 
  TLegend *lL = new TLegend(0.6,0.65,0.9,0.9); lL->SetFillColor(0); lL->SetBorderSize(0);
  //TLegend *lL = new TLegend(0.50, 0.65, 0.95, 0.88); lL->SetFillColor(0); lL->SetBorderSize(0);
  iH[iHiggs]->Add(iH[iHiggs+1]);
  //iH[iHiggs]->Add(iH[iHiggs+2]);
  double error;
  //cout << "Higgs Expected Signal  " << iH[iHiggs]->IntegralAndError(0,iH[iHiggs]->GetNbinsX(),error) << endl;
  cout << "Higgs Expected Signal  " << iH[iHiggs]->IntegralAndError(6,7,error) << endl;
  cout << "" << error << endl;
  //iH[iHiggs]->Scale(5);
  iH[4]->Add(iH[2]);
  //iH[0]->Add(iH[3]);
  for(int i0 = 0; i0 < iN-1; i0++) {
    //if(i0==2 || i0==3) continue;
    if(i0==2) continue;
    iH[i0]->SetLineColor(kBlack);
    iH[i0]->SetLineWidth(3);
    iH[i0]->SetFillColor(fColor[i0]);
    iH[i0]->SetFillStyle(1001);
    iH[i0]->SetMarkerStyle(20);
    iH[i0]->SetMarkerColor(fColor[i0]);
    iH[i0]->SetMarkerSize (0.6);
    //if(i0 != iHiggs) iH[i0]->SetFillColor(fColor[i0]);
    for(int i1 = i0+1; i1 < iHiggs; i1++)  {
      //if(i1==2 || i1==3) continue;
      if(i1==2) continue;
      iH[i0]->Add(iH[i1]);
    } 		    
    //if(i0 <= iHiggs && i0 !=2 && i0 !=3 && iH[i0]->Integral() > 0.0001) lL->AddEntry(iH[i0],fString[i0].c_str(),"f");
    //if(i0 < iHiggs && i0 !=2 && i0 !=3 && iH[i0]->Integral() > 0.0001) lL->AddEntry(iH[i0],fString[i0].c_str(),"f");
    //iH[i0]->SetLineWidth(1); iH[i0]->SetLineColor(kBlack);
  }
 
  iH[iHiggs]->Add(iH[0]);
  iH[iHiggs]->SetLineWidth(1); 
  iH[iHiggs]->SetLineColor(kBlue);
  iH[iN-1]->SetLineWidth(3);
  iH[iN-1]->SetMarkerSize(1.3); 
  iH[iN-1]->SetMarkerStyle(20);
  iH[iN-1]->SetMarkerColor(kBlack);
  lL->AddEntry(iH[iN-1],fString[iN-1].c_str(),"lp");
  lL->AddEntry(iH[iHiggs],"5X H(125)->#tau#tau","f");
  lL->AddEntry(iH[0],fString[0].c_str(),"f");
  lL->AddEntry(iH[3],fString[3].c_str(),"f");
  lL->AddEntry(iH[4],fString[4].c_str(),"f");
  lL->AddEntry(iH[1],fString[1].c_str(),"f");
  lL->AddEntry(iH[5],fString[5].c_str(),"f");
  //scale(iH,iH[iN-1]->Integral(),iN);
  //iH[iN-1]->GetYaxis()->SetRangeUser(0.1,1500);
  //lC0->Divide(1,2); lC0->cd();  lC0->cd(1)->SetPad(0,0.2,1.0,1.0); gPad->SetLeftMargin(0.2) ;
  TH1F* errorBand = (TH1F*)iH[0] ->Clone();
  errorBand  ->SetMarkerSize(0);
  errorBand  ->SetFillColor(1);
  errorBand  ->SetFillStyle(3013);
  errorBand  ->SetLineWidth(1);
  lL->AddEntry(errorBand,"bkg. uncertainty","f");
  //for(int idx=0; idx<errorBand->GetNbinsX(); ++idx){
   // if(errorBand->GetBinContent(idx)>0){
    //  errorBand->SetBinError(idx,0.084873*errorBand->GetBinContent(idx));      
   // }
   // else
    //  errorBand->SetBinError(idx,0);	
 // }
  
  iH[iN-1]->SetNdivisions(505);
  iH[iN-1]->SetMinimum(0.01);
  iH[iN-1]->SetMaximum(std::max(maximum(iH[iN-1], 0), maximum(iH[0], 0)));
  if(iFName.compare("eta_1") ==0 || iFName.compare("eta_2") ==0 || iFName.compare("jeta_1") ==0 || iFName.compare("jeta_2") ==0)
    {
      iH[iN-1]->SetMaximum(1.5*(std::max(maximum(iH[iN-1], 0), maximum(iH[0], 0))));
    }
  //iH[iN-1]->SetMaximum(3.0*(std::max(maximum(iH[iN-1], 0), maximum(iH[0], 0))));
  iH[iN-1]->Draw("EP");
  iH[iHiggs]->Draw("hist sames");
   
  for(int i0 = 0; i0 < iN-1; i0++) {
    if(i0 < iHiggs)
      //if(i0 !=2 && i0 != 3)
      if(i0 !=2)
	iH[i0]->Draw("hist sames");
  }
  iH[iN-1]  ->Draw("EP sames");
  errorBand->Draw("e2same"); 
  lL->Draw();

  cout << "Kolmogorov Test  " << iH[0]->KolmogorovTest(iH[iN-1]) << endl;
  cout << "Chi 2 test  " << iH[0]->Chi2Test(iH[iN-1]) << endl;
  //cmsprelim("#tau_{h}#tau_{h}", 0.18, 0.93);
  TPaveText* cmsprel  = new TPaveText(0.15, 0.93, 0.15+0.25, 0.93+0.04, "NDC");
  cmsprel->SetBorderSize(   0 );
  cmsprel->SetFillStyle(    0 );
  cmsprel->SetTextAlign(   12 );
  cmsprel->SetTextSize ( 0.04);
  cmsprel->SetTextColor(    1 );
  cmsprel->SetTextFont (   62 );
  cmsprel->AddText("CMS Preliminary, #sqrt{s} = 8 TeV, L = 19.8 fb^{-1}");
  cmsprel->Draw();
  
  TPaveText* chan     = new TPaveText(0.20, 0.74+0.061, 0.32, 0.74+0.161, "NDC");
  chan->SetBorderSize(   0 );
  chan->SetFillStyle(    0 );
  chan->SetTextAlign(   12 );
  chan->SetTextSize ( 0.05 );
  chan->SetTextColor(    1 );
  chan->SetTextFont (   62 );
  chan->AddText("#tau_{h}#tau_{h}");
  //chan->AddText("e#tau_{h}");
  chan->Draw();

  TPaveText* chan1     = new TPaveText(0.60, 0.50, 0.62, 0.50, "NDC");
  chan1->SetBorderSize(   0 );
  chan1->SetFillStyle(    0 );
  chan1->SetTextAlign(   12 );
  chan1->SetTextSize ( 0.04 );
  chan1->SetTextColor(    1 );
  chan1->SetTextFont (   62 );
  chan1->AddText("M_{#tau#tau}<100 GeV/c^{2}");
  //chan1->Draw();

  lC0->RedrawAxis();
  //lC0->GetYaxis()->SetRangeUser("")
  //lC0->SetLogy(1);

  //lC0->cd(2)->SetPad(0,0,1.0,0.2); gPad->SetLeftMargin(0.2) ;
  //drawDifference(iH[iN-1],iH[0]);//,iH[6],iH[7]);
  //lC0->Print((iFName+"_tt.png").c_str());
  if(iCat==0)
    lC0->SaveAs((iName+"_vbf_tight.png").c_str());
  else if(iCat==1) 
    lC0->SaveAs((iName+"_vbf_loose.png").c_str());
  else if(iCat==2)
    lC0->SaveAs((iName+"_0jet_low.png").c_str());
  else if(iCat==3)
    lC0->SaveAs((iName+"_0jet_medium.png").c_str());
  else if(iCat==4) 
    lC0->SaveAs((iName+"0jet_high.png").c_str());
  else if(iCat==5) 
    lC0->SaveAs((iName+"_1jet_medium.png").c_str());
  else if(iCat==6) 
    lC0->SaveAs((iName+"_1jet_mediumhiggs.png").c_str());
  else if(iCat==7) 
    lC0->SaveAs((iName+"_1jet_highhiggs.png").c_str());
  else
    lC0->SaveAs((iName+"_insclusive.png").c_str());
  // lC0->SaveAs("insclusive.png");
  //lC0->SaveAs("boosthigh.png");
  delete lC0;
  delete lL;
  delete cmsprel;
  delete chan;
}
TH1F* draw(std::string iVar,int iCat,TTree *iTree,int iId,std::string iCut,std::string iSId="B",int high=0,int theory=0,int mssm=0,TH1F *iH=0) {
  std::stringstream lName;   lName   <<  "Met"  << iId << iSId;
  std::stringstream lVar;
  lVar   << iVar << ">>+" << lName.str();
  std::string lCut   = fWeights[iId];
  lCut+=iCut;
  cout << "Cut : " << lCut << endl;
  std::string lHCutSH,lLCutSL;
  TString lNCut = lCut.c_str();
  if(high == 1)
    {
      lNCut.ReplaceAll("pt_1>45&&pt_2>45","pt_1high>45&&pt_2high>45");
      lHCutSH = lNCut.Data();
      cout << "Cut 1 : " << lHCutSH << endl;
    }
  else if (high==2)
    {
      lNCut.ReplaceAll("pt_1>45&&pt_2>45","pt_1low>45&&pt_2low>45");
      //lNCut.ReplaceAll("pt_2","pt_2low");       
      lLCutSL = lNCut.Data();
      cout << "Cut 2 : " << lLCutSL << endl;
    }
  
  //  if(high == 1)
//     {
//       lNCut.ReplaceAll("pt_1>17&&pt_2>30","pt_1high>17&&pt_2high>30");
//       lHCutSH = lNCut.Data();
//       cout << "Cut 1 : " << lHCutSH << endl;
//     }
//   else if (high==2)
//     {
//       lNCut.ReplaceAll("pt_1>17&&pt_2>30","pt_1low>17&&pt_2low>30");
//       //lNCut.ReplaceAll("pt_2","pt_2low");       
//       lLCutSL = lNCut.Data();
//       cout << "Cut 2 : " << lLCutSL << endl;
//     }

  if(theory == 1)
    {
      lNCut.ReplaceAll("*(weight","*(weighthigh");
      lHCutSH = lNCut.Data();
      cout << "Cut Theory high : " << lHCutSH << endl;
    }
  else if (theory==2)
    {
      lNCut.ReplaceAll("*(weight","*(weightlow");
      //lNCut.ReplaceAll("pt_2","pt_2low");       
      lLCutSL = lNCut.Data();
      cout << "Cut Thoery low : " << lLCutSL << endl;
    }

  TH1F *lMet = iH;
  double massLEdges[27]    = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,225.,250.,275.,300.,325.,350.};
  double massLEdgesVbf[14]    = {0.,20.,40.,60.,80.,100.,120.,140.,160.,180.,200.,250.,300.,350.};
  //double massLEdges[32]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500};
  //double massLEdges[32]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,500,700,1000,1500};
  //double massLEdges[56]    = {0.,5.,10.,15.,20.,25.,30.,35.,40.,45.,50.,55.,60.,65.,70.,75.,80.,85.,90.,95.,100.,105.,110.,115.,120.,125.,130.,135.,140.,145.,150.,155.,160.,165.,170.,175.,180.,185.,190.,195.,200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,300.,310.,320.,330.,340.,350.};
  //double massLEdges[48]={0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,180,190,200,225,250,275,300,325,350,400,450,500,550,600,650,700,750,800,850,900,950,1000,1100,1200,1300,1400,1500,1600,1800,2000};
  //double massLEdges[43]    = {0.,7.5,15.0,22.5,30.0,37.5,45.0,52.5,60.0,67.5,75.0,82.5,90.,97.5,105.,112.5,120.,127.5,135.,142.5,150.,157.5,165.,172.5,180.,187.5,195.,202.5,210.,220.,230.,240.,250.,260.,270.,280.,290.,300.,310.,320.,330.,340.,350.};
  //double massLEdgesVbf[19]={0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500,700,1000,1500};
  //double massLEdgesVbf[19]={0,20,40,60,80,100,120,140,160,180,200,250,300,350,400,500};
  //double massLEdgesVbf[27]    = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,225.,250.,275.,300.,325.,350.};
  if(iH == 0 && !mssm) 
    //if(iH == 0) 
    {
      if(iVar.compare("m_sv") ==0 || iVar.compare("m_svhigh")==0 || iVar.compare("m_svlow")==0 || iVar.compare("mvis")==0 || iVar.compare("m_vishigh")==0 || iVar.compare("m_vislow")==0 || iVar.compare("m_sv*0.98")==0 || iVar.compare("m_sv*1.02")==0) 
      	{
	  // if(iCat==1)
	  //   lMet = new TH1F(lName.str().c_str(),"",13,massLEdgesVbf);
	  //if(iCat==1 || iCat==9)
	  //if(iCat==1 || iCat==5)
	  if(iCat==1 || iCat==0)
	    lMet = new TH1F(lName.str().c_str(),"",13,massLEdgesVbf);
	  else
	    lMet = new TH1F(lName.str().c_str(),"",26,massLEdges);	  
	}
      else
	lMet = new TH1F(lName.str().c_str(),"",getNBins(iVar),getXMin(iVar),getXMax(iVar));
    }
  else
    if(iH == 0) lMet = new TH1F(lName.str().c_str(),"",getNBins(iVar),getXMin(iVar),getXMax(iVar));
  lMet->SetLineColor(fColor[iId]);
  lMet->SetLineWidth(1);
  lMet->GetXaxis()->SetTitle(getXAxis(iVar));
  lMet->GetYaxis()->SetTitle(getYAxis(iVar));
  TString lDraw(lVar.str().c_str());
  if(high==1 || theory==1)
    iTree->Draw(lDraw,lHCutSH.c_str());
  else if (high==2 || theory ==2)
    iTree->Draw(lDraw,lLCutSL.c_str());
  else
    iTree->Draw(lDraw,lCut.c_str());
  lMet->SetMarkerStyle(kFullCircle);
  lMet->SetLineWidth(2);
  return lMet;
}

TH2F* draw2D(std::string iVar,std::string iVar1,TTree *iTree,int iId,std::string iCut,TH2F *iH = 0) {
  std::stringstream lName;  lName   <<  "Met"  << iId;// << iVar;
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
TH1F * merge(std::string iName,double iMergePoint,TH1F *iH,TH1F *iFunc) { 
  TH1F *lH = (TH1F*) iH->Clone(iName.c_str());
  int lMergeBin = iH->GetXaxis()->FindBin(iMergePoint);
  double lVal  = iH->GetBinContent(lMergeBin);
  iFunc->Scale(lVal/iFunc->GetBinContent(lMergeBin));
  for(int i0 = 0;         i0 < lMergeBin;         i0++) lH->SetBinContent(i0,iH->GetBinContent(i0));
  for(int i0 = lMergeBin; i0 < iH->GetNbinsX()+1; i0++) lH->SetBinContent(i0,iFunc->GetBinContent(i0));
  return lH;
}
