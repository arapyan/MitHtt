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
   \class   etauAfterFit_b etauAfterFit_b.C "MitHtt/macros/etauAfterFit_b.c"

   \brief   macro to create mvis plots with normalization after fit

   macro to create mvis plots for EPS2011 analysis after fit as used for the limit calculation. 
   A description how to retrieve the shifts as applied by the fit can be found on:

    + SWGuideHiggsAnalysisCombinedLimit#Maximum_likelihood_fits_and_diag

   The macro expects the following files to be present in your worling 
   directory:

    + HttStyles.h
    + HttStyles.cc
    + root/eleTau_mssm.root

   The root file is expected to contain the histograms that were the input to the limit cal-
   culation. To run the macro do the following:

   root -l 
   .L HttStyle.cc++
   .L etauAfterFit_b.C++ 
   etauAfterFit_b()

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
  double CMS_eff_e             = 1.0392; // +1.96 * 1.02
  double CMS_eff_t             = 1.0090; // +0.15 * 1.06 
  double CMS_scale_j           = 0.8775; // -1.75 * 1.07
  double CMS_eff_b             = 1.0774; // +0.73 * 1.106
  double CMS_fake_b            = 1.0756; // +0.54 * 1.14
  double CMS_htt_zttNorm       = 1.0099; // +0.30 * 1.033
  double CMS_htt_ttbarNorm     = 1.0816; // +0.68 * 1.12
  double CMS_htt_DiBosonNorm   = 0.0000; // -1.31 * 2.00
  double CMS_htt_QCDNorm       = 0.9130; // -0.87 * 1.10
  double CMS_htt_QCDSyst       = 1.0517; // +1.10 * 1.047
  double CMS_htt_WNorm         = 0.9894; // -0.15 * 1.071
  double CMS_htt_ZJFake        = 0.9898; // -0.09 * 1.113
  double CMS_htt_ZLFake        = 0.9789; // -0.27 * 1.078

  switch(idx){
  case 1: //ZTT 
    hin->Scale(CMS_eff_e*CMS_eff_t*CMS_htt_zttNorm); break;
  case 2: // QCD
    hin->Scale(CMS_htt_QCDNorm*CMS_htt_QCDSyst); break;
  case 3: // W
    hin->Scale(CMS_htt_WNorm); break;
  case 4: // ZJ
    hin->Scale(CMS_eff_e*CMS_htt_zttNorm*CMS_htt_ZJFake); break;
  case 5: // ZL
    hin->Scale(CMS_eff_e*CMS_htt_zttNorm*CMS_htt_ZLFake); break;
  case 6: // TT
    hin->Scale(CMS_eff_t*CMS_eff_e*CMS_htt_ttbarNorm*CMS_eff_b); break;
  case 7: // VV
    hin->Scale(CMS_eff_t*CMS_eff_e*CMS_htt_DiBosonNorm); break;
  case 8: // ggH
    hin->Scale(lumi*CMS_eff_t*CMS_eff_e*CMS_scale_j*CMS_fake_b); break;
  case 9: // bbH
    hin->Scale(lumi*CMS_eff_t*CMS_eff_e*CMS_scale_j*CMS_eff_b); break;
  default :
    std::cout << "error histograms not known?!?" << std::endl;
  }
}

// examples macro
void 
etauAfterFit_b(bool scaled = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/eleTau_mssm.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("eleTau_B/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find hitogram " << "eleTau_B/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("eleTau_B/QCD"))              ; InitHist(Fakes, "", "", kMagenta-10, 1001);                   
  TH1F* EWK1  =  refill((TH1F*)exampleFile->Get("eleTau_B/W"  ))              ; InitHist(EWK1 , "", "", kRed    + 2, 1001);
  TH1F* EWK2  =  refill((TH1F*)exampleFile->Get("eleTau_B/ZJ" ))              ; InitHist(EWK2 , "", "", kRed    + 2, 1001);
  TH1F* EWK3  =  refill((TH1F*)exampleFile->Get("eleTau_B/ZL" ))              ; InitHist(EWK3 , "", "", kRed    + 2, 1001);
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("eleTau_B/VV" ))              ; InitHist(EWK  , "", "", kRed    + 2, 1001);
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("eleTau_B/TT" ))              ; InitHist(ttbar, "", "", kBlue   - 8, 1001);
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("eleTau_B/ZTT"  ))            ; InitHist(Ztt  , "", "", kOrange - 4, 1001);
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("eleTau_B/GGH120"))           ; InitSignal(ggH); ggH ->Scale(43.1834*0.114107*0.042/ggH ->Integral());
  TH1F* bbH   =  refill((TH1F*)exampleFile->Get("eleTau_B/BBH120"))           ; InitSignal(bbH); bbH ->Scale(66.7725*0.114107*1.364/bbH ->Integral());

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
  data->SetMaximum(50.);
  data->SetNdivisions(505);
  data->Draw("e");

  bbH->Draw("same");
  Ztt->Draw("same");
  ttbar->Draw("same");
  EWK->Draw("same");
  Fakes->Draw("same");
  data->Draw("esame");
  canv->RedrawAxis();

  CMSPrelim("#tau_{e}#tau_{h}", 0.45, 0.75);
  
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

  if(scaled) canv->Print("etau_rescaled_b.pdf"); else canv->Print("etau_unscaled_b.pdf");
  if(scaled) canv->Print("etau_rescaled_b.png"); else canv->Print("etau_unscaled_b.png");
}
