#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <sstream>
#include "TLegend.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TRandom1.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "RooRealVar.h"
#include "RooFormulaVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooAddPdf.h"
#include "RooGenericPdf.h"
#include "RooExtendPdf.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "RooBreitWigner.h"
#include "RooCBShape.h"
#include "RooFFTConvPdf.h"
#include "RooNumConvPdf.h"
#include "RooPolynomial.h"
#include "RooPlot.h"
#include "RooVoigtian.h"
#include "RooVoigtianShape.h"
#include "RooMCStudy.h"
#include "TNtuple.h"
#include "Common/CPlot.hh"
#include "TRandom1.h"
#include "TLorentzVector.h"

#include "EScale.hh"

TString dir("865-41x");

EScale esc("$CMSSW_BASE/src/MitHtt/Emu/EScale/"+dir+"/data-EnergyScale.root",
	   "$CMSSW_BASE/src/MitHtt/Emu/EScale/"+dir+"/mc-EnergyScale.root");

void testCorrection() {
  TH1F *hmassraw = new TH1F("hmassraw","hmassraw",25,75,105);
  TH1F *hmass = new TH1F("hmass","hmass",25,75,105);
  TH1F *hmassScDown = new TH1F("hmassScDown","hmassScDown",25,75,105);
  TH1F *hmassScUp = new TH1F("hmassScUp","hmassScUp",25,75,105);
  TH1F *hmassResDown = new TH1F("hmassResDown","hmassResDown",25,75,105);
  TH1F *hmassResUp = new TH1F("hmassResUp","hmassResUp",25,75,105);
  TH1F *hmassdata = new TH1F("hmassdata","hmassdata",25,75,105);

  TFile *lFile = new TFile(dir+"/zee-flat.root");
  TTree *lTree = (TTree*) lFile->FindObjectAny("Events");
  float lmass  = 0; lTree->SetBranchAddress("mass"    ,&lmass);
  float lpt1   = 0; lTree->SetBranchAddress("pt1"    ,&lpt1);
  float leta1  = 0; lTree->SetBranchAddress("eta1"   ,&leta1);
  float lphi1  = 0; lTree->SetBranchAddress("phi1"   ,&lphi1);
  float lpt2   = 0; lTree->SetBranchAddress("pt2" ,&lpt2);
  float leta2  = 0; lTree->SetBranchAddress("eta2",&leta2);
  float lphi2  = 0; lTree->SetBranchAddress("phi2",&lphi2);

  for(int i0 = 0; i0 < lTree->GetEntries();i0++) { 
    lTree->GetEntry(i0);

    TLorentzVector l1; l1.SetPtEtaPhiM(lpt1,leta1,lphi1,0);
    TLorentzVector l2; l2.SetPtEtaPhiM(lpt2,leta2,lphi2,0);
    hmassraw->Fill((l1+l2).M());

    l1.SetPtEtaPhiM(esc.pt(leta1,lpt1,kCenter),leta1,lphi1,0);
    l2.SetPtEtaPhiM(esc.pt(leta2,lpt2,kCenter),leta2,lphi2,0);
    hmass->Fill((l1+l2).M());

    l1.SetPtEtaPhiM(esc.pt(leta1,lpt1,kScDown),leta1,lphi1,0);
    l2.SetPtEtaPhiM(esc.pt(leta2,lpt2,kScDown),leta2,lphi2,0);
    hmassScDown->Fill((l1+l2).M());

    l1.SetPtEtaPhiM(esc.pt(leta1,lpt1,kScUp),leta1,lphi1,0);
    l2.SetPtEtaPhiM(esc.pt(leta2,lpt2,kScUp),leta2,lphi2,0);
    hmassScUp->Fill((l1+l2).M());

    l1.SetPtEtaPhiM(esc.pt(leta1,lpt1,kResDown),leta1,lphi1,0);
    l2.SetPtEtaPhiM(esc.pt(leta2,lpt2,kResDown),leta2,lphi2,0);
    hmassResDown->Fill((l1+l2).M());

    l1.SetPtEtaPhiM(esc.pt(leta1,lpt1,kResUp),leta1,lphi1,0);
    l2.SetPtEtaPhiM(esc.pt(leta2,lpt2,kResUp),leta2,lphi2,0);
    hmassResUp->Fill((l1+l2).M());

  }

  lFile->Close();

  lFile = new TFile(dir+"/data-flat.root");
  lTree = (TTree*) lFile->FindObjectAny("Events");
  lTree->SetBranchAddress("mass"    ,&lmass);
  lTree->SetBranchAddress("pt1"    ,&lpt1);
  lTree->SetBranchAddress("eta1"   ,&leta1);
  lTree->SetBranchAddress("phi1"   ,&lphi1);
  lTree->SetBranchAddress("pt2" ,&lpt2);
  lTree->SetBranchAddress("eta2",&leta2);
  lTree->SetBranchAddress("phi2",&lphi2);

  for(int i0 = 0; i0 < lTree->GetEntries();i0++) { 
    lTree->GetEntry(i0);

    TLorentzVector l1; l1.SetPtEtaPhiM(lpt1,leta1,lphi1,0);
    TLorentzVector l2; l2.SetPtEtaPhiM(lpt2,leta2,lphi2,0);
    hmassdata->Fill((l1+l2).M());

  }
  lFile->Close();

  TCanvas *c = new TCanvas("c","c");

  CPlot plotmassSc("mass_sc","","m_{ee} [GeV]","Events");
  hmassScDown->Scale(hmassdata->Integral(0,101)/hmassScDown->Integral(0,101));  
  hmassScUp->Scale(hmassdata->Integral(0,101)/hmassScUp->Integral(0,101));  
  hmassResDown->Scale(hmassdata->Integral(0,101)/hmassResDown->Integral(0,101));  
  hmassResUp->Scale(hmassdata->Integral(0,101)/hmassResUp->Integral(0,101));  
  hmass->Scale(hmassdata->Integral(0,101)/hmass->Integral(0,101));  
  hmassraw->Scale(hmassdata->Integral(0,101)/hmassraw->Integral(0,101));  
  plotmassSc.AddHist1D(hmassScDown,"+/- scale","hist",kRed);
  plotmassSc.AddHist1D(hmassScUp,"",kRed);
  plotmassSc.AddHist1D(hmass,"corrected MC","hist",kBlue);
  plotmassSc.AddHist1D(hmassdata,"data","E");
  plotmassSc.TransLegend(0.08,0);
  plotmassSc.Draw(c,kTRUE,"png");

  CPlot plotmassRes("mass_res","","m_{ee} [GeV]","Events");
  plotmassRes.AddHist1D(hmassResDown,"+/- res.","hist",kRed);
  plotmassRes.AddHist1D(hmassResUp,"",kRed);
  plotmassRes.AddHist1D(hmass,"corrected MC","hist",kBlue);
  plotmassRes.AddHist1D(hmassdata,"data","E");
  plotmassRes.TransLegend(0.08,0);
  plotmassRes.Draw(c,kTRUE,"png");

  CPlot plotmasscorr("mass_corr","","m_{ee} [GeV]","Events");
  plotmasscorr.AddHist1D(hmass,"corrected MC","hist",kBlue);
  plotmasscorr.AddHist1D(hmassraw,"raw  MC","hist",kRed);
  plotmasscorr.AddHist1D(hmassdata,"data","E");
  plotmasscorr.TransLegend(0.08,0);
  plotmasscorr.Draw(c,kTRUE,"png");

  // print a table with the corrections and (statistical) uncertainties
  esc.print();
}
