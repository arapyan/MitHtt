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
   \class   emuAfterFit_novbf emuAfterFit_novbf.C "MitHtt/macros/emuAfterFit_novbf.c"

   \brief   macro to create mvis plots with normalization after fit

   macro to create mvis plots for EPS2011 analysis after fit as used for the limit calculation. 
   A description how to retrieve the shifts as applied by the fit can be found on:

    + SWGuideHiggsAnalysisCombinedLimit#Maximum_likelihood_fits_and_diag

   The macro expects the following files to be present in your worling 
   directory:

    + HttStyles.h
    + HttStyles.cc
    + root/eleMu.root

   The root file is expected to contain the histograms that were the input to the limit cal-
   culation. To run the macro do the following:

   root -l 
   .L MitStyle.cc++
   .L emuAfterFit_novbf.C++ 
   emuAfterFit_b()

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
  double lumi                = 1.0090; // +0.15 * 1.06
  double CMS_eff_e           = 1.0336; // +1.68 * 1.02
  double CMS_eff_m           = 1.0094; // +0.47 * 1.02 
  double CMS_scale_j         = 1.0000; // nan
  double CMS_htt_zttNorm     = 1.0318; // +1.06 * 1.03
  double CMS_htt_ttbarNorm   = 0.9350; // -0.65 * 1.10
  double CMS_htt_DiBosonNorm = 0.7840; // -0.72 * 1.30
  double CMS_hww_fakes_em    = 0.9370; // -0.21 * 1.30
  double pdf_gg              = 1.0561; // +1.87 * 1.03
  double pdf_qqbar           = 1.0561; // +1.87 * 1.03
  double QCDScale_qqHin      = 1.2244; // +1.87 * 1.12
  double QCDScale_ggH        = 1.0654; // +1.87 * 1.035
  double ueps                = 0.9252; // +1.87 * 1.04


  switch(idx){
  case 1: //Ztt 
    hin->Scale(CMS_eff_e*CMS_eff_m*CMS_scale_j*CMS_htt_zttNorm); break;
  case 2: // ttbar
    hin->Scale(CMS_eff_e*CMS_eff_m*CMS_scale_j*CMS_htt_ttbarNorm); break;
  case 3: // EWK
    hin->Scale(lumi*CMS_eff_e*CMS_eff_m*CMS_scale_j*CMS_htt_DiBosonNorm); break;
  case 4: // Fakes
    hin->Scale(CMS_eff_e*CMS_eff_m*CMS_scale_j*CMS_hww_fakes_em); break;
  case 5: // ggH
    hin->Scale(lumi*CMS_eff_e*CMS_eff_m*CMS_scale_j*pdf_gg*QCDScale_ggH*ueps); break;
  case 6: // qqH
    hin->Scale(lumi*CMS_eff_e*CMS_eff_m*CMS_scale_j*pdf_qqbar*QCDScale_qqHin*ueps); break;
  default :
    std::cout << "error histograms not known?!?" << std::endl;
  }
}

// examples macro
void 
emuAfterFit_novbf(bool scaled = true, bool log = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/eleMu.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("emu_novbf/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find hitogram " << "emu_novbf/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("emu_novbf/Fakes"))            ; InitHist(Fakes, "", "", kMagenta-10, 1001);                   
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("emu_novbf/EWK"  ))            ; InitHist(EWK  , "", "", kRed    + 2, 1001);
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("emu_novbf/ttbar"))            ; InitHist(ttbar, "", "", kBlue   - 8, 1001);
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("emu_novbf/Ztt"  ))            ; InitHist(Ztt  , "", "", kOrange - 4, 1001);
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("emu_novbf/Higgs_gf_sm_120"  )); InitSignal(ggH); ggH->Scale(10*16.63*0.071*10.2558/ggH->Integral());
  TH1F* qqH   =  refill((TH1F*)exampleFile->Get("emu_novbf/Higgs_vbf_sm_120" )); InitSignal(qqH); qqH->Scale(10*1.269*0.071* 0.8178/qqH->Integral());

  if(scaled){
    rescale(Fakes, 4); 
    rescale(EWK,   3); 
    rescale(ttbar, 2); 
    rescale(Ztt,   1); 
    rescale(ggH,   5);
    rescale(qqH,   6);
  }
  if(log){
    qqH  ->Add(ggH  );
    Fakes->Add(qqH  );
    EWK  ->Add(Fakes);
    ttbar->Add(EWK  );
    Ztt  ->Add(ttbar);
    //ggH  ->Add(Ztt  );
    //qqH  ->Add(ggH  );
  }
  else{
    EWK  ->Add(Fakes);
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
    data->SetMinimum(0.5);
    data->SetMaximum(5000.);
  }
  else{
    data->SetMaximum(1500.);
  }
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

  CMSPrelim("#tau_{e}#tau_{#mu}", 0.45, 0.75);
  
  TLegend* leg = new TLegend(0.45, 0.45, 0.9, 0.75);
  SetLegendStyle(leg);
  leg->AddEntry(qqH  , "(10x) H#rightarrow#tau#tau" , "L" );
  leg->AddEntry(data , "Observed"                , "LP");
  leg->AddEntry(Ztt  , "Z#rightarrow#tau#tau"    , "F" );
  leg->AddEntry(ttbar, "t#bar{t}"                , "F" );
  leg->AddEntry(EWK  , "Electroweak"             , "F" );
  leg->AddEntry(Fakes, "Fakes"                   , "F" );
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
    if(scaled) canv->Print("emu_rescaled_novbf_LOG.png"); else canv->Print("emu_unscaled_novbf_LOG.png");
    if(scaled) canv->Print("emu_rescaled_novbf_LOG.pdf"); else canv->Print("emu_unscaled_novbf_LOG.pdf");
  }
  else{
    if(scaled) canv->Print("emu_rescaled_novbf.png");     else canv->Print("emu_unscaled_novbf.png");
    if(scaled) canv->Print("emu_rescaled_novbf.pdf");     else canv->Print("emu_unscaled_novbf.pdf");
  }
}
