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
   \class   emuAfterFit_vbf emuAfterFit_vbf.C "MitHtt/macros/emuAfterFit_vbf.c"

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
   .L emuAfterFit_vbf.C++ 
   emuAfterFit_vbf()

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
std::cout<< "scaling by 0.961886"<<std::endl;hin->Scale(0.961886); break; 
  case 2: // ttbar
std::cout<< "scaling by 0.951176"<<std::endl;hin->Scale(0.951176); break; 
  case 3: // EWK
std::cout<< "scaling by 0.703666"<<std::endl;hin->Scale(0.703666); break; 
  case 4: // Fakes
std::cout<< "scaling by 1.083042"<<std::endl;hin->Scale(1.083042); break; 
  case 5: // ggH
std::cout<< "scaling by 0.963926"<<std::endl;hin->Scale(0.963926); break; 
  case 6: // qqH
std::cout<< "scaling by 0.963926"<<std::endl;hin->Scale(0.963926); break; 
  default :
    std::cout << "error histograms not known?!?" << std::endl;
  }
}

// examples macro
void 
emuAfterFit_vbf(bool scaled = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/eleMu.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("emu_vbf/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find hitogram " << "emu_vbf/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("emu_vbf/Fakes"))            ; InitHist(Fakes, "", "", kMagenta-10, 1001);                   
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("emu_vbf/EWK"  ))            ; InitHist(EWK  , "", "", kRed    + 2, 1001);
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("emu_vbf/ttbar"))            ; InitHist(ttbar, "", "", kBlue   - 8, 1001);
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("emu_vbf/Ztt"  ))            ; InitHist(Ztt  , "", "", kOrange - 4, 1001);
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("emu_vbf/Higgs_gf_sm_120"  )); InitSignal(ggH); ggH->Scale(10*16.63*0.071);
  TH1F* qqH   =  refill((TH1F*)exampleFile->Get("emu_vbf/Higgs_vbf_sm_120" )); InitSignal(qqH); qqH->Scale(10*1.269*0.071);

  if(scaled){
    rescale(Fakes, 4); 
    rescale(EWK,   3); 
    rescale(ttbar, 2); 
    rescale(Ztt,   1); 
    rescale(ggH,   5);
    rescale(qqH,   6);
  }
  EWK  ->Add(Fakes);
  ttbar->Add(EWK  );
  Ztt  ->Add(ttbar);
  ggH  ->Add(Ztt  );
  qqH  ->Add(ggH  );

  // define canvas
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);

  canv->cd();
  data->SetMaximum(5.);
  data->Draw("e");

  qqH->Draw("same");
  Ztt->Draw("same");
  ttbar->Draw("same");
  EWK->Draw("same");
  Fakes->Draw("same");
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

  if(scaled) canv->Print("emu_rescaled_vbf.pdf"); else canv->Print("emu_unscaled_vbf.pdf");
  if(scaled) canv->Print("emu_rescaled_vbf.png"); else canv->Print("emu_unscaled_vbf.png");


 // TCanvas* canv1 = new TCanvas("Legend", "Legend", 600, 600);
 // canv1->cd();
 //
 // TLegend* leg1 = new TLegend(0.1, 0.1, 0.9, 0.9);
 // SetLegendStyle(leg1);
 // leg1->AddEntry(qqH  , "(10x) H#rightarrow#tau#tau  m_{H}=120" , "L" );
 // leg1->AddEntry(data , "Observed"                , "LP");
 // leg1->AddEntry(Ztt  , "Z#rightarrow#tau#tau"    , "F" );
 // leg1->AddEntry(ttbar, "t#bar{t}"                , "F" );
 // leg1->AddEntry(EWK  , "Electroweak"             , "F" );
 // leg1->AddEntry(Fakes, "Fakes"                   , "F" );
 // leg1->Draw();
 //
 // TPaveText* mssm1  = new TPaveText(0.67, 0.76, 0.90, 0.90, "NDC");
 // mssm1->SetBorderSize(   0 );
 // mssm1->SetFillStyle( 1001 );
 // mssm1->SetFillColor(kWhite);
 // mssm1->SetTextAlign(   12 );
 // mssm1->SetTextSize ( 0.05 );
 // mssm1->SetTextColor(    1 );
 // mssm1->SetTextFont (   62 );
 // mssm1->AddText("m_{H}=120");
 // mssm1->Draw();
}
