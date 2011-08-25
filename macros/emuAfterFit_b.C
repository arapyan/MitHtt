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
   \class   emuAfterFit_b emuAfterFit_b.C "MitHtt/macros/emuAfterFit_b.c"

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
   .L HttStyle.cc++
   .L emuAfterFit_b.C++ 
   emuAfterFit_b()

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
  gStyle->SetLineStyleString(11,"20 10");
  switch(idx){
  case 1: //Ztt 
std::cout<< "scaling by 1.055638"<<std::endl;hin->Scale(1.055638); break; 
  case 2: // ttbar
std::cout<< "scaling by 1.119990"<<std::endl;hin->Scale(1.119990); break; 
  case 3: // EWK
std::cout<< "scaling by 0.681870"<<std::endl;hin->Scale(0.681870); break; 
  case 4: // Fakes
std::cout<< "scaling by 1.155123"<<std::endl;hin->Scale(1.155123); break; 
  case 5: // ggH
std::cout<< "scaling by 1.055527"<<std::endl;hin->Scale(1.055527); break; 
  case 6: // bbH
std::cout<< "scaling by 1.055527"<<std::endl;hin->Scale(1.055527); break; 
  default :
    std::cout << "error histograms not known?!?" << std::endl;
  }
}

// examples macro
void 
emuAfterFit_b(bool scaled = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/eleMu.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("emu_b/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find hitogram " << "emu_b/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("emu_b/Fakes"))            ; InitHist(Fakes, "", "", kMagenta-10, 1001);                   
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("emu_b/EWK"  ))            ; InitHist(EWK  , "", "", kRed    + 2, 1001); 
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("emu_b/ttbar"))            ; InitHist(ttbar, "", "", kBlue   - 8, 1001); 
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("emu_b/Ztt"  ))            ; InitHist(Ztt  , "", "", kOrange - 4, 1001); 
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("emu_b/Higgs_gg_mssm_120")); InitSignal(ggH); ggH->Scale(43.1834*0.114107); 
  TH1F* bbH   =  refill((TH1F*)exampleFile->Get("emu_b/Higgs_bb_mssm_120")); InitSignal(bbH); bbH->Scale(66.7725*0.114107); 

  if(scaled){
    rescale(Fakes, 4); 
    rescale(EWK,   3); 
    rescale(ttbar, 2); 
    rescale(Ztt,   1); 
    rescale(ggH,   5);
    rescale(bbH,   6);
  }
  EWK  ->Add(Fakes);
  ttbar->Add(EWK  );
  Ztt  ->Add(ttbar);
  ggH  ->Add(Ztt  );
  bbH  ->Add(ggH  );

  // define canvas
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);

  canv->cd();
  data->SetMaximum(60.);
  data->Draw("e");

  bbH->Draw("same");
  Ztt->Draw("same");
  ttbar->Draw("same");
  EWK->Draw("same");
  Fakes->Draw("same");
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

  if(scaled) canv->Print("emu_rescaled_b.png"); else canv->Print("emu_unscaled_b.png");
  if(scaled) canv->Print("emu_rescaled_b.pdf"); else canv->Print("emu_unscaled_b.pdf");


 // TCanvas* canv1 = new TCanvas("Legend", "Legend", 600, 600);
 // canv1->cd();
 //
 // TLegend* leg1 = new TLegend(0.1, 0.1, 0.9, 0.9);
 // SetLegendStyle(leg1);
 // leg1->AddEntry(bbH  , "#phi#rightarrow#tau#tau" , "L" );
 // leg1->AddEntry(data , "Observed"                , "LP");
 // leg1->AddEntry(Ztt  , "Z#rightarrow#tau#tau"    , "F" );
 // leg1->AddEntry(ttbar, "t#bar{t}"                , "F" );
 // leg1->AddEntry(EWK  , "Electroweak"             , "F" );
 // leg1->AddEntry(Fakes, "Fakes"                   , "F" );
 // leg1->Draw();
 //
 // TPaveText* mssm1  = new TPaveText(0.55, 0.75, 0.90, 0.90, "NDC");
 // mssm1->SetBorderSize(   0 );
 // mssm1->SetFillStyle(    0 );
 // mssm1->SetTextAlign(   12 );
 // mssm1->SetTextSize ( 0.05 );
 // mssm1->SetTextColor(    1 );
 // mssm1->SetTextFont (   62 );
 // mssm1->AddText("(m_{A}=120, tan#beta=20)");
 // mssm1->Draw();
}
