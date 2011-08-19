#include <iostream>

#include <TH1F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <Rtypes.h>

#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TPaveText.h>

#include "HttStyles.h"

/**
   \class   mutauAfterFit_novbf mutauAfterFit_novbf.C "MitHtt/macros/mutauAfterFit_novbf.c"

   \brief   macro to create mvis plots with normalization after fit

   macro to create mvis plots for EPS2011 analysis after fit as used for the limit calculation. 
   A description how to retrieve the shifts as applied by the fit can be found on:

    + SWGuideHiggsAnalysisCombinedLimit#Maximum_likelihood_fits_and_diag

   The macro expects the following files to be present in your worling 
   directory:

    + HttStyles.h
    + HttStyles.cc
    + root/muTau_sm.root

   The root file is expected to contain the histograms that were the input to the limit cal-
   culation. To run the macro do the following:

   root -l 
   .L HttStyle.cc++
   .L mutauAfterFit_novbf.C++ 
   mutauAfterFit_b()

   There is a set of arguments that you can give to the function:

   + bool scaled : true  - scale to normalization after fit
                   false - leave unscaled 
   + bool log    : true  - plot in log 
                   false - plot in linear
*/

// re-fill histograms (this is only a little helper for the example histogram)
TH1F* refill(TH1F* hin)
{
  TH1F* hout = (TH1F*)hin->Clone(); hout->Clear();
  for(int i=0; i<hout->GetNbinsX(); ++i){
    hout->SetBinContent(i+1, hin->GetBinContent(i+1));
    hout->SetBinError  (i+1, 0.);
  }
  return hout;
}

// rescale histograms according to fit
void rescale(TH1F* hin, unsigned int idx)
{
  double lumi                  = 1.0013; // +0.03 * 1.045
  double CMS_eff_m             = 1.0223; // +2.23 * 1.01
  double CMS_eff_t             = 0.9598; // -0.67 * 1.06 
  double CMS_scale_j           = 0.9922; // +0.78 * 0.99
  double CMS_htt_zttNorm       = 0.9974; // -0.08 * 1.033
  double CMS_htt_ttbarNorm     = 1.0122; // +0.12 * 1.102
  double CMS_htt_DiBosonNorm   = 0.2600; // -0.74 * 2.00
  double CMS_htt_QCDNorm       = 1.0690; // +1.15 * 1.06
  double CMS_htt_QCDSyst       = 0.9310; // -1.21 * 1.019
  double CMS_htt_WNorm         = 1.0015; // +0.15 * 1.01
  double CMS_htt_WSyst         = 1.0583; // +0.94 * 1.062
  double CMS_htt_ZJFake        = 1.0250; // +0.21 * 1.119
  double CMS_htt_ZLFake        = 1.0696; // +0.27 * 1.258
  double pdf_gg                = 1.0000; //  0.00 * 1.03
  double pdf_qqbar             = 1.0000; //  0.00 * 1.03
  double QCDScale_qqH          = 1.0000; //  0.00 * 1.12
  double QCDScale_ggH          = 1.0000; //  0.00 * 1.035
  double ueps                  = 0.; //  0.00 * 0.96

  switch(idx){
  case 1: //ZTT 
    hin->Scale(CMS_eff_m*CMS_eff_t*CMS_htt_zttNorm*CMS_scale_j); break;
  case 2: // QCD
    hin->Scale(CMS_htt_QCDNorm*CMS_htt_QCDSyst); break;
  case 3: // W
    hin->Scale(CMS_htt_WNorm*CMS_htt_WSyst*CMS_scale_j); break;
  case 4: // ZJ
    hin->Scale(CMS_eff_m*CMS_htt_zttNorm*CMS_htt_ZJFake*CMS_scale_j); break;
  case 5: // ZL
    hin->Scale(CMS_eff_m*CMS_htt_zttNorm*CMS_htt_ZLFake*CMS_scale_j); break;
  case 6: // TT
    hin->Scale(CMS_eff_t*CMS_eff_m*CMS_htt_ttbarNorm*CMS_scale_j); break;
  case 7: // VV
    hin->Scale(CMS_eff_t*CMS_eff_m*CMS_htt_DiBosonNorm*CMS_scale_j); break;
  case 8: // ggH
    hin->Scale(lumi*CMS_eff_t*CMS_eff_m*CMS_scale_j*pdf_gg*QCDScale_ggH*ueps); break;
  case 9: // qqH
    hin->Scale(lumi*CMS_eff_t*CMS_eff_m*CMS_scale_j*pdf_qqbar*QCDScale_qqH*ueps); break;
  default :
    std::cout << "error histograms not known?!?" << std::endl;
  }
}

// examples macro
void 
mutauAfterFit_novbf(bool scaled = true, bool log = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/muTau_sm.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("muTau_SM0/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find hitogram " << "muTau_SM0/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("muTau_SM0/QCD"))              ; InitHist(Fakes, "", "", kMagenta-10, 1001);                   
  TH1F* EWK1  =  refill((TH1F*)exampleFile->Get("muTau_SM0/W"  ))              ; InitHist(EWK1 , "", "", kRed    + 2, 1001);
  TH1F* EWK2  =  refill((TH1F*)exampleFile->Get("muTau_SM0/ZJ" ))              ; InitHist(EWK2 , "", "", kRed    + 2, 1001);
  TH1F* EWK3  =  refill((TH1F*)exampleFile->Get("muTau_SM0/ZL" ))              ; InitHist(EWK3 , "", "", kRed    + 2, 1001);
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("muTau_SM0/VV" ))              ; InitHist(EWK  , "", "", kRed    + 2, 1001);
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("muTau_SM0/TT" ))              ; InitHist(ttbar, "", "", kBlue   - 8, 1001);
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("muTau_SM0/ZTT"))              ; InitHist(Ztt  , "", "", kOrange - 4, 1001);
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("muTau_SM0/SM120" ))           ; InitSignal(ggH); ggH ->Scale(10*16.63*0.071*29.671/ggH ->Integral());
  TH1F* qqH   =  refill((TH1F*)exampleFile->Get("muTau_SM0/VBF120"))           ; InitSignal(qqH); qqH ->Scale(10*1.269*0.071* 2.147/qqH ->Integral());
 
  if(scaled){
    rescale(Fakes, 2); 
    rescale(EWK1 , 3); 
    rescale(EWK2 , 4); 
    rescale(EWK3 , 5); 
    rescale(EWK  , 7); 
    rescale(ttbar, 6); 
    rescale(Ztt  , 1);
    rescale(ggH  , 8); 
    rescale(qqH  , 9);  
  }
  if(log){
    qqH  ->Add(ggH  );
    Fakes->Add(qqH  );
    EWK1 ->Add(Fakes);
    EWK2 ->Add(EWK1 );
    EWK3 ->Add(EWK2 );
    EWK  ->Add(EWK3 );
    ttbar->Add(EWK  );
    Ztt  ->Add(ttbar);
  }
  else{
    EWK1 ->Add(Fakes);
    EWK2 ->Add(EWK1 );
    EWK3 ->Add(EWK2 );
    EWK  ->Add(EWK3 );
    ttbar->Add(EWK  );
    Ztt  ->Add(ttbar);
    ggH  ->Add(Ztt  );
    qqH  ->Add(ggH  );
  }
  // define canvas
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);

  canv->cd();
  if(log){
    canv->SetLogy(1);
    data->SetMinimum(5.0);
    data->SetMaximum(10000000.);
  }
  else{
    data->SetMaximum(14000.);
  }
  data->SetNdivisions(505);
  data->Draw("e");

  if(log){
    Ztt->Draw("same");
    ttbar->Draw("same");
    EWK->Draw("same");
    Fakes->Draw("same");
    qqH->Draw("same");
  }
  else{
    qqH->Draw("same");
    Ztt->Draw("same");
    ttbar->Draw("same");
    EWK->Draw("same");
    Fakes->Draw("same");
  }
  data->Draw("esame");
  canv->RedrawAxis();

  CMSPrelim("#tau_{#mu}#tau_{h}", 0.45, 0.75);
  
  TLegend* leg = new TLegend(0.45, 0.45, 0.9, 0.75);
  SetLegendStyle(leg);
  leg->AddEntry(qqH  , "(10x) H#rightarrow#tau#tau" , "L" );
  leg->AddEntry(data , "Observed"                , "LP");
  leg->AddEntry(Ztt  , "Z#rightarrow#tau#tau"    , "F" );
  leg->AddEntry(ttbar, "t#bar{t}"                , "F" );
  leg->AddEntry(EWK  , "Electroweak"             , "F" );
  leg->AddEntry(Fakes, "QCD"                     , "F" );
  leg->Draw();

  TPaveText* mssm  = new TPaveText(0.78, 0.70, 0.90, 0.74, "NDC");
  mssm->SetBorderSize(   0 );
  mssm->SetFillStyle(    0 );
  mssm->SetTextAlign(   12 );
  mssm->SetTextSize ( 0.04 );
  mssm->SetTextColor(    1 );
  mssm->SetTextFont (   62 );
  mssm->AddText("m_{H}=120");
  mssm->Draw();

  if(log){
    if(scaled) canv->Print("mutau_rescaled_novbf_LOG.pdf"); else canv->Print("mutau_unscaled_novbf_LOG.pdf");
    if(scaled) canv->Print("mutau_rescaled_novbf_LOG.png"); else canv->Print("mutau_unscaled_novbf_LOG.png");

  }
  else{
    if(scaled) canv->Print("mutau_rescaled_novbf.pdf");     else canv->Print("mutau_unscaled_novbf.pdf");
    if(scaled) canv->Print("mutau_rescaled_novbf.png");     else canv->Print("mutau_unscaled_novbf.png");
  }
}
