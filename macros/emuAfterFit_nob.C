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
   \class   emuAfterFit_nob emuAfterFit_nob.C "MitHtt/macros/emuAfterFit_nob.c"

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
   .L emuAfterFit_nob.C++ 
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
  double lumi                  = 0.9934; // -0.11 * 1.06
  double CMS_eff_e             = 1.1006; // +5.03 * 1.02
  double CMS_eff_m             = 0.9986; // -0.07 * 1.02 
  double CMS_scale_j           = 1.0000; // nan
  double CMS_eff_b             = 1.0000; // nan
  double CMS_htt_zttNorm       = 1.0123; // +0.41 * 1.03
  double CMS_htt_ttbarNorm     = 0.9690; // -0.31 * 1.10
  double CMS_htt_DiBosonNorm   = 0.8440; // -0.52 * 1.30
  double CMS_hww_fakes_em      = 0.4480; // -1.84 * 1.30

  switch(idx){
  case 1: //Ztt 
    hin->Scale(CMS_eff_e*CMS_eff_m*CMS_eff_b*CMS_scale_j*CMS_htt_zttNorm); break;
  case 2: // ttbar
    hin->Scale(CMS_eff_e*CMS_eff_m*CMS_eff_b*CMS_scale_j*CMS_htt_ttbarNorm); break;
  case 3: // EWK
    hin->Scale(lumi*CMS_eff_e*CMS_eff_m*CMS_eff_b*CMS_scale_j*CMS_htt_DiBosonNorm); break;
  case 4: // Fakes
    hin->Scale(CMS_eff_e*CMS_eff_m*CMS_eff_b*CMS_scale_j*CMS_hww_fakes_em); break;
  case 5: // ggH
    hin->Scale(lumi*CMS_eff_e*CMS_eff_m*CMS_eff_b*CMS_scale_j); break;
  case 6: // bbH
    hin->Scale(lumi*CMS_eff_e*CMS_eff_m*CMS_eff_b*CMS_scale_j); break;
  default :
    std::cout << "error histograms not known?!?" << std::endl;
  }
}

// examples macro
void 
emuAfterFit_nob(bool scaled = true, bool log = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/eleMu.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("emu_nob/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find hitogram " << "emu_nob/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("emu_nob/Fakes"))            ; InitHist(Fakes, "", "", kMagenta-10, 1001);                   
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("emu_nob/EWK"  ))            ; InitHist(EWK  , "", "", kRed    + 2, 1001);
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("emu_nob/ttbar"))            ; InitHist(ttbar, "", "", kBlue   - 8, 1001);
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("emu_nob/Ztt"  ))            ; InitHist(Ztt  , "", "", kOrange - 4, 1001);
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("emu_nob/Higgs_gg_mssm_120")); InitSignal(ggH); ggH->Scale(43.1834*0.114107*2.9351/ggH->Integral());
  TH1F* bbH   =  refill((TH1F*)exampleFile->Get("emu_nob/Higgs_bb_mssm_120")); InitSignal(bbH); bbH->Scale(66.7725*0.114107*6.0323/bbH->Integral());

  if(scaled){
    rescale(Fakes, 4); 
    rescale(EWK,   3); 
    rescale(ttbar, 2); 
    rescale(Ztt,   1); 
    rescale(ggH,   5);
    rescale(bbH,   6);
  }
  if(log){
    bbH  ->Add(ggH  );
    Fakes->Add(bbH  );
    EWK  ->Add(Fakes);
    ttbar->Add(EWK  );
    Ztt  ->Add(ttbar);
    //ggH  ->Add(Ztt  );
    //bbH  ->Add(ggH  );
  }
  else{
    EWK  ->Add(Fakes);
    ttbar->Add(EWK  );
    Ztt  ->Add(ttbar);
    ggH  ->Add(Ztt  );
    bbH  ->Add(ggH  );
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
    bbH->Draw("same");
  }
  else{
    bbH->Draw("same");
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
  leg->AddEntry(bbH  , "#phi#rightarrow#tau#tau" , "L" );
  leg->AddEntry(data , "Observed"                , "LP");
  leg->AddEntry(Ztt  , "Z#rightarrow#tau#tau"    , "F" );
  leg->AddEntry(ttbar, "t#bar{t}"                , "F" );
  leg->AddEntry(EWK  , "Electroweak"             , "F" );
  leg->AddEntry(Fakes, "Fakes"                   , "F" );
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

  if(log){
    if(scaled) canv->Print("emu_rescaled_nob_LOG.png"); else canv->Print("emu_unscaled_nob_LOG.png");      
    if(scaled) canv->Print("emu_rescaled_nob_LOG.pdf"); else canv->Print("emu_unscaled_nob_LOG.pdf");      
  }
  else{
    if(scaled) canv->Print("emu_rescaled_nob.png");     else canv->Print("emu_unscaled_nob.png");
    if(scaled) canv->Print("emu_rescaled_nob.pdf");     else canv->Print("emu_unscaled_nob.pdf");
  }
}
