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
   \class   mutauAfterFit_vbf mutauAfterFit_vbf.C "MitHtt/macros/mutauAfterFit_vbf.c"

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
   .L mutauAfterFit_vbf.C++ 
   mutauAfterFit_vbf()

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
  case 1: //ZTT 
std::cout<<"ZTT scaling by 0.926084"<<std::endl;hin->Scale(0.926084);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 
  case 2: // QCD
std::cout<<"QCD scaling by 1.003610"<<std::endl;hin->Scale(1.003610);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 
  case 3: // W
std::cout<<"W scaling by 1.098560"<<std::endl;hin->Scale(1.098560);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 
  case 4: // ZJ
std::cout<<"ZJ scaling by 0.964887"<<std::endl;hin->Scale(0.964887);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 
  case 5: // ZL
std::cout<<"ZL scaling by 0.951843"<<std::endl;hin->Scale(0.951843);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 
  case 6: // TT
std::cout<<"TT scaling by 0.908412"<<std::endl;hin->Scale(0.908412);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 
  case 7: // VV
std::cout<<"VV scaling by 0.095106"<<std::endl;hin->Scale(0.095106);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 
  case 8: // ggH
std::cout<<"SM120 scaling by 0.945494"<<std::endl;hin->Scale(0.945494);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 
  case 9: // bbH
std::cout<<"VBF120 scaling by 0.945494"<<std::endl;hin->Scale(0.945494);std::cout<<hin->GetName()<<"Scaled Integral is "<<hin->Integral()<< std::endl; break; 

  default :
    std::cout << "error histograms not known?!?" << std::endl;
  }
}

// examples macro
void 
mutauAfterFit_vbf(bool scaled = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/muTau_sm.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("muTau_SM2/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find histogram " << "muTau_SM2/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("muTau_SM2/QCD"))              ; InitHist(Fakes, "", "", kMagenta-10, 1001);                   
  TH1F* EWK1  =  refill((TH1F*)exampleFile->Get("muTau_SM2/W"  ))              ; InitHist(EWK1 , "", "", kRed    + 2, 1001);
  TH1F* EWK2  =  refill((TH1F*)exampleFile->Get("muTau_SM2/ZJ" ))              ; InitHist(EWK2 , "", "", kRed    + 2, 1001);
  TH1F* EWK3  =  refill((TH1F*)exampleFile->Get("muTau_SM2/ZL" ))              ; InitHist(EWK3 , "", "", kRed    + 2, 1001);
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("muTau_SM2/VV" ))              ; InitHist(EWK  , "", "", kRed    + 2, 1001);
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("muTau_SM2/TT" ))              ; InitHist(ttbar, "", "", kBlue   - 8, 1001);
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("muTau_SM2/ZTT"))              ; InitHist(Ztt  , "", "", kOrange - 4, 1001);
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("muTau_SM2/SM120" ))           ; InitSignal(ggH); ggH ->Scale(10*1);
  TH1F* qqH   =  refill((TH1F*)exampleFile->Get("muTau_SM2/VBF120"))           ; InitSignal(qqH); qqH ->Scale(10*1);

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
  EWK1 ->Add(Fakes);
  EWK2 ->Add(EWK1 );
  EWK3 ->Add(EWK2 );
  EWK  ->Add(EWK3 );
  ttbar->Add(EWK  );
  Ztt  ->Add(ttbar);
  ggH  ->Add(Ztt  );
  qqH  ->Add(ggH  );

  // define canvas
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);

  canv->cd();
  data->SetMaximum(25.);
  data->SetNdivisions(505);
  data->Draw("e");

  qqH->Draw("same");
  Ztt->Draw("same");
  ttbar->Draw("same");
  EWK->Draw("same");
  Fakes->Draw("same");
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

  if(scaled) canv->Print("mutau_rescaled_vbf.pdf"); else canv->Print("mutau_unscaled_vbf.pdf");
  if(scaled) canv->Print("mutau_rescaled_vbf.png"); else canv->Print("mutau_unscaled_vbf.png");
}
