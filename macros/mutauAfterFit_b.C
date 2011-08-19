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
   \class   mutauAfterFit_b mutauAfterFit_b.C "MitHtt/macros/mutauAfterFit_b.c"

   \brief   macro to create mvis plots with normalization after fit

   macro to create mvis plots for EPS2011 analysis after fit as used for the limit calculation. 
   A description how to retrieve the shifts as applied by the fit can be found on:

    + SWGuideHiggsAnalysisCombinedLimit#Maximum_likelihood_fits_and_diag

   The macro expects the following files to be present in your worling 
   directory:

    + HttStyles.h
    + HttStyles.cc
    + root/muTau_mssm.root

   The root file is expected to contain the histograms that were the input to the limit cal-
   culation. To run the macro do the following:

   root -l 
   .L HttStyle.cc++
   .L mutauAfterFit_b.C++ 
   mutauAfterFit_b()

   There is a set of arguments that you can give to the function:

   + bool scaled : true  - scale to normalization after fit
                   false - leave unscaled 
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
  double lumi                  = 0.9946; // -0.12 * 1.045
  double CMS_eff_m             = 0.9901; // -0.99 * 1.01
  double CMS_eff_t             = 1.0090; // +0.15 * 1.06 
  double CMS_scale_j           = 0.8250; // -1.75 * 1.10
  double CMS_eff_b             = 1.0774; // +0.73 * 1.106
  double CMS_fake_b            = 1.0702; // +0.54 * 1.13
  double CMS_htt_zttNorm       = 1.0099; // +0.30 * 1.033
  double CMS_htt_ttbarNorm     = 1.0327; // +0.68 * 1.109
  double CMS_htt_DiBosonNorm   = 0.0000; // -1.31 * 2.00
  double CMS_htt_QCDNorm       = 1.0243; // +0.27 * 1.09
  double CMS_htt_QCDSyst       = 1.0061; // +0.32 * 1.019
  double CMS_htt_WNorm         = 1.0098; // +0.14 * 1.070
  double CMS_htt_ZJFake        = 1.0024; // +0.02 * 1.120
  double CMS_htt_ZLFake        = 0.9819; // -0.07 * 1.259

  switch(idx){
  case 1: //ZTT 
    hin->Scale(CMS_eff_m*CMS_eff_t*CMS_htt_zttNorm); break;
  case 2: // QCD
    hin->Scale(CMS_htt_QCDNorm*CMS_htt_QCDSyst); break;
  case 3: // W
    hin->Scale(CMS_htt_WNorm); break;
  case 4: // ZJ
    hin->Scale(CMS_eff_m*CMS_htt_zttNorm*CMS_htt_ZJFake); break;
  case 5: // ZL
    hin->Scale(CMS_eff_m*CMS_htt_zttNorm*CMS_htt_ZLFake); break;
  case 6: // TT
    hin->Scale(CMS_eff_t*CMS_eff_m*CMS_htt_ttbarNorm*CMS_eff_b); break;
  case 7: // VV
    hin->Scale(CMS_eff_t*CMS_eff_m*CMS_htt_DiBosonNorm); break;
  case 8: // ggH
    hin->Scale(lumi*CMS_eff_t*CMS_eff_m*CMS_scale_j*CMS_fake_b); break;
  case 9: // bbH
    hin->Scale(lumi*CMS_eff_t*CMS_eff_m*CMS_scale_j*CMS_eff_b); break;
  default :
    std::cout << "error histograms not known?!?" << std::endl;
  }
}

// examples macro
void 
mutauAfterFit_b(bool scaled = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/muTau_mssm.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("muTau_B/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find hitogram " << "muTau_B/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("muTau_B/QCD"))              ; InitHist(Fakes, "", "", kMagenta-10, 1001); 
  TH1F* EWK1  =  refill((TH1F*)exampleFile->Get("muTau_B/W"  ))              ; InitHist(EWK1 , "", "", kRed    + 2, 1001); 
  TH1F* EWK2  =  refill((TH1F*)exampleFile->Get("muTau_B/ZJ" ))              ; InitHist(EWK2 , "", "", kRed    + 2, 1001); 
  TH1F* EWK3  =  refill((TH1F*)exampleFile->Get("muTau_B/ZL" ))              ; InitHist(EWK3 , "", "", kRed    + 2, 1001); 
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("muTau_B/VV" ))              ; InitHist(EWK  , "", "", kRed    + 2, 1001); 
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("muTau_B/TT" ))              ; InitHist(ttbar, "", "", kBlue   - 8, 1001); 
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("muTau_B/ZTT"))              ; InitHist(Ztt  , "", "", kOrange - 4, 1001); 
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("muTau_B/GGH120" ))          ; InitSignal(ggH); ggH->Scale(43.1834*0.114107*0.083/ggH->Integral());
  TH1F* bbH   =  refill((TH1F*)exampleFile->Get("muTau_B/BBH120" ))          ; InitSignal(bbH); bbH->Scale(66.7725*0.114107*2.282/bbH->Integral()); 

  if(scaled){
    rescale(Fakes, 2); 
    rescale(EWK1 , 3); 
    rescale(EWK2 , 4); 
    rescale(EWK3 , 5); 
    rescale(EWK  , 7); 
    rescale(ttbar, 6); 
    rescale(Ztt  , 1);
    rescale(ggH  , 8);  
    rescale(bbH  , 9);
  }
  EWK1 ->Add(Fakes);
  EWK2 ->Add(EWK1 );
  EWK3 ->Add(EWK2 );
  EWK  ->Add(EWK3 );
  ttbar->Add(EWK  );
  Ztt  ->Add(ttbar);
  ggH  ->Add(Ztt  );
  bbH  ->Add(ggH  );


  // define canvas
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);

  canv->cd();
  data->SetMaximum(100.);
  data->SetNdivisions(505);
  data->Draw("e");

  bbH->Draw("same");
  Ztt->Draw("same");
  ttbar->Draw("same");
  EWK->Draw("same");
  Fakes->Draw("same");
  data->Draw("esame");
  canv->RedrawAxis();

  CMSPrelim("#tau_{#mu}#tau_{h}", 0.45, 0.75);
  
  TLegend* leg = new TLegend(0.45, 0.45, 0.9, 0.75);
  SetLegendStyle(leg);
  leg->AddEntry(bbH  , "#phi#rightarrow#tau#tau" , "L" );
  leg->AddEntry(data , "Observed"                , "LP");
  leg->AddEntry(Ztt  , "Z#rightarrow#tau#tau"    , "F" );
  leg->AddEntry(ttbar, "t#bar{t}"                , "F" );
  leg->AddEntry(EWK  , "Electroweak"             , "F" );
  leg->AddEntry(Fakes, "QCD"                     , "F" );
  leg->Draw();

  TPaveText* mssm  = new TPaveText(0.66, 0.70, 0.90, 0.74, "NDC");
  mssm->SetBorderSize(   0 );
  mssm->SetFillStyle(    0 );
  mssm->SetTextAlign(   12 );
  mssm->SetTextSize ( 0.03 );
  mssm->SetTextColor(    1 );
  mssm->SetTextFont (   62 );
  mssm->AddText("(m_{A}=120, tan#beta=20)");
  mssm->Draw();

  if(scaled) canv->Print("mutau_rescaled_b.pdf"); else canv->Print("mutau_unscaled_b.pdf");
  if(scaled) canv->Print("mutau_rescaled_b.png"); else canv->Print("mutau_unscaled_b.png");
}
