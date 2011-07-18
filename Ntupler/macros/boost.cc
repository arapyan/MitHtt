#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TH1.h>                    // histogram base class
#include <TNtuple.h>                  // class to access ntuples
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TCanvas.h>
#include <TMath.h>                  // ROOT math functions
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "Common/MitStyleRemix.hh"  // style settings for drawing
#include "Common/CSample.hh"        // helper class for organizing input ntuple files
#include "Common/MyTools.hh"        // miscellaneous helper functions
#include "Common/CPlot.hh"          // helper class for plots

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh" 
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"   
#include "MitHtt/Ntupler/interface/TVertex.hh"   
#endif
using namespace std;

void boost()
{
  vector<TString> filev;

  filev.push_back("s11-zjetsm50-mg-v11-pu_0000_ntuple.root");
  filev.push_back("s11-h100tt-gf-v1g1-pu_0000_ntuple.root");

  vector<TH1F*> hCosThStv;
  vector<TH1F*> hE1v;

  char hname[100];
  for(UInt_t ifile=0; ifile<filev.size(); ifile++) {
    sprintf(hname,"hCosThSt_%i",ifile);    hCosThStv.push_back(new TH1F(hname,"",50,-1.1,1.1));   hCosThStv[ifile]->Sumw2();
    sprintf(hname,"hE1_%i",ifile);         hE1v.push_back(new TH1F(hname,"",50,0,1000));            hE1v[ifile]->Sumw2();
  }

  TLorentzVector *tau1 = new TLorentzVector;
  TLorentzVector *tau2 = new TLorentzVector;
  TLorentzVector ditau;
  for(UInt_t ifile=0; ifile<filev.size(); ifile++) {
    TFile infile(filev[ifile]);
    TTree *intree=0;
    infile.GetObject("Taus",intree); assert(intree);
    intree->SetBranchAddress("tau1",&tau1);
    intree->SetBranchAddress("tau2",&tau2);
    
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);

      Double_t en = (tau1->E() > tau2->E()) ? tau1->E() : tau2->E();

      ditau = *tau1 + *tau2;

      tau1->Boost(-ditau.BoostVector());
      tau2->Boost(-ditau.BoostVector());
      Float_t cths = tau1->Vect()*ditau.Vect()/(tau1->P()*ditau.P());

      hCosThStv[ifile]->Fill(cths);
      hE1v[ifile]->Fill(en);
    }
    
    infile.Close();
  }
  delete tau1; delete tau2;

  TCanvas c1("c1","c1");
  Double_t int0, int1;

  CPlot plotCosThSt("cosThSt","","cos(#theta^{*})","norm.");
  hCosThStv[0]->SetLineColor(kRed);
  int0 = hCosThStv[0]->Integral(0,hCosThStv[0]->GetNbinsX()+1);
  int1 = hCosThStv[1]->Integral(0,hCosThStv[1]->GetNbinsX()+1);
  hCosThStv[0]->Scale(int1/int0);
  plotCosThSt.AddHist1D(hCosThStv[0],filev[0](4,12),"hist");
  plotCosThSt.AddHist1D(hCosThStv[1],filev[1](4,12),"hist",kBlue);
  plotCosThSt.Draw(&c1,kTRUE,"png");

  CPlot plotE1("E1","","E1 [GeV]","norm.");
  hE1v[0]->SetLineColor(kRed);
  int0 = hE1v[0]->Integral(0,hE1v[0]->GetNbinsX()+1);
  int1 = hE1v[1]->Integral(0,hE1v[1]->GetNbinsX()+1);
  hE1v[0]->Scale(int1/int0);
  plotE1.AddHist1D(hE1v[0],filev[0](4,12),"hist");
  plotE1.AddHist1D(hE1v[1],filev[1](4,12),"hist",kBlue);
  plotE1.Draw(&c1,kTRUE,"png");
  plotE1.SetLogy();
  plotE1.SetName("E1-log");
  plotE1.Draw(&c1,kTRUE,"png");

  for(UInt_t ifile=0; ifile<filev.size(); ifile++) {
    delete hCosThStv[ifile];
    delete hE1v[ifile];
  }
  
}
