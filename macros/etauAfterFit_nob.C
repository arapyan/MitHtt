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
   \class   etauAfterFit_nob etauAfterFit_nob.C "MitHtt/macros/etauAfterFit_nob.c"

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
   .L etauAfterFit_nob.C++ 
   etauAfterFit_b()

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
  //std::cout << "What the hell? Print something out, please..." << hin->GetName() << hin->GetTitle() << idx << std::endl;
  gStyle->SetLineStyleString(11,"20 10");
  switch(idx){
  case 1: //ZTT 
std::cout<< "scaling by 1.042256"<<std::endl;hin->Scale(1.042256); break; 
  case 2: // QCD
std::cout<< "scaling by 1.015510"<<std::endl;hin->Scale(1.015510); break; 
  case 3: // W
std::cout<< "scaling by 1.005200"<<std::endl;hin->Scale(1.005200); break; 
  case 4: // ZJ
std::cout<< "scaling by 1.033018"<<std::endl;hin->Scale(1.033018); break; 
  case 5: // ZL
std::cout<< "scaling by 1.166692"<<std::endl;hin->Scale(1.166692); break; 
  case 6: // TT
std::cout<< "scaling by 1.019881"<<std::endl;hin->Scale(1.019881); break; 
  case 7: // VV
std::cout<< "scaling by 0.000000"<<std::endl;hin->Scale(0.000000); break; 
  case 8: // ggHnoJet
std::cout<< "scaling by 1.019218"<<std::endl;hin->Scale(1.019218); break; 
  case 9: // ggHJet
std::cout<< "scaling by 1.007803"<<std::endl;hin->Scale(1.007803); break; 
  case 10: // bbHNoJet
std::cout<< "scaling by 1.019218"<<std::endl;hin->Scale(1.019218); break; 
  case 11: // bbHJet
std::cout<< "scaling by 1.007803"<<std::endl;hin->Scale(1.007803); break; 
  case 12: // bbHBJet
std::cout<< "scaling by 0.945753"<<std::endl;hin->Scale(0.945753); break; 

  default :
    std::cout << "error histograms not known I AM HERE PLEASE?!?" << "IDX is: " << idx << " with name " << hin->GetName() << hin->GetTitle()<< std::endl;
  }
}

// examples macro
void 
etauAfterFit_nob(bool scaled = true, bool log = true)
{
  // defining the common canvas, axes pad styles
  SetStyle();

  // open example histogram file
  TFile* exampleFile = new TFile("root/eleTau_mssm.root");

  //load example histograms
  TH1F* data = (TH1F*)exampleFile->Get("eleTau_NoB/data_obs");
  if(data) {InitHist(data, "#bf{m_{vis} [GeV]}", "#bf{Events}"); InitData(data);} else{std::cout << "can't find hitogram " << "eleTau_NoB/data_obs" << std::endl;}

  TH1F* Fakes =  refill((TH1F*)exampleFile->Get("eleTau_NoB/QCD"))              ; InitHist(Fakes, "", "", kMagenta-10, 1001); 
  TH1F* EWK1  =  refill((TH1F*)exampleFile->Get("eleTau_NoB/W"  ))              ; InitHist(EWK1 , "", "", kRed    + 2, 1001); 
  TH1F* EWK2  =  refill((TH1F*)exampleFile->Get("eleTau_NoB/ZJ" ))              ; InitHist(EWK2 , "", "", kRed    + 2, 1001); 
  TH1F* EWK3  =  refill((TH1F*)exampleFile->Get("eleTau_NoB/ZL" ))              ; InitHist(EWK3 , "", "", kRed    + 2, 1001); 
  TH1F* EWK   =  refill((TH1F*)exampleFile->Get("eleTau_NoB/VV" ))              ; InitHist(EWK  , "", "", kRed    + 2, 1001); 
  TH1F* ttbar =  refill((TH1F*)exampleFile->Get("eleTau_NoB/TT" ))              ; InitHist(ttbar, "", "", kBlue   - 8, 1001); 
  TH1F* Ztt   =  refill((TH1F*)exampleFile->Get("eleTau_NoB/ZTT"  ))            ; InitHist(Ztt  , "", "", kOrange - 4, 1001); 
  TH1F* ggH1  =  refill((TH1F*)exampleFile->Get("eleTau_NoB/GGHNoJet120"))      ; InitSignal(ggH1); ggH1->Scale(43.1834*0.114107);
  TH1F* ggH   =  refill((TH1F*)exampleFile->Get("eleTau_NoB/GGHJet120"  ))      ; InitSignal(ggH ); ggH ->Scale(43.1834*0.114107); 
  TH1F* bbH1  =  refill((TH1F*)exampleFile->Get("eleTau_NoB/BBHNoJet120"))      ; InitSignal(bbH1); bbH1->Scale(66.7725*0.114107); 
  TH1F* bbH2  =  refill((TH1F*)exampleFile->Get("eleTau_NoB/BBHJet120"  ))      ; InitSignal(bbH2); bbH2->Scale(66.7725*0.114107); 
  TH1F* bbH   =  refill((TH1F*)exampleFile->Get("eleTau_NoB/BBHBJet120" ))      ; InitSignal(bbH ); bbH ->Scale(66.7725*0.114107); 

  if(scaled){
    rescale(Fakes, 2); 
    rescale(EWK1 , 3); 
    rescale(EWK2 , 4); 
    rescale(EWK3 , 5); 
    rescale(EWK  , 7); 
    rescale(ttbar, 6); 
    rescale(Ztt  , 1);
    rescale(ggH1 , 8); 
    rescale(ggH  , 9);  
    rescale(bbH1 ,10); 
    rescale(bbH2 ,11); 
    rescale(bbH  ,12);
  }
  if(log){
    ggH  ->Add(ggH1 );
    bbH1 ->Add(ggH  );
    bbH2 ->Add(bbH1 );
    bbH  ->Add(bbH2 );
    Fakes->Add(bbH  );
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
    ggH1 ->Add(Ztt  );
    ggH  ->Add(ggH1 );
    bbH1 ->Add(ggH  );
    bbH2 ->Add(bbH1 );
    bbH  ->Add(bbH2 );
  }

  // define canvas
  TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);

  canv->cd();
  if(log){
    canv->SetLogy(1);
    data->SetMinimum(5.0);
    data->SetMaximum(50000.);
  }
  else{
    data->SetMaximum(6000.);
  }
  data->SetNdivisions(505);
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

  if(log){
    if(scaled) canv->Print("etau_rescaled_nob_LOG.pdf"); else canv->Print("etau_unscaled_nob_LOG.pdf");
    if(scaled) canv->Print("etau_rescaled_nob_LOG.png"); else canv->Print("etau_unscaled_nob_LOG.png");
  }
  else{
    if(scaled) canv->Print("etau_rescaled_nob.pdf"); else canv->Print("etau_unscaled_nob.pdf");
    if(scaled) canv->Print("etau_rescaled_nob.png"); else canv->Print("etau_unscaled_nob.png");
  }
}
