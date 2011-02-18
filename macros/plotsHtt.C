//--------------------------------------------------------------------------------------------------
// Perform a plot task using a specified set of samples. Nice and clean.
//
// Authors: C.Paus                                                                        (Aug 2010)
//--------------------------------------------------------------------------------------------------
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <iostream>
#include "TH1.h"
#include "TFile.h"
#include "TLatex.h"
#include "TPostScript.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "MitPlots/Style/interface/MitStyle.h"
#include "MitPlots/Input/interface/TaskSamples.h"
#include "MitPlots/Plot/interface/PlotTask.h"
#endif

using namespace std;
using namespace mithep;

void plot(const char *name, const char* title, int logy,
          double xmin, double xmax, double ymin, double ymax,
          int nRebin, double lumi, TString draw="", TString cut="", int nbins=100);

//==================================================================================================
void plotsHtt(double lumi = 3000.)
{
  // setup graphics stuff before starting
  MitStyle::Init();
  gROOT->Macro("$CMSSW_BASE/src/MitHtt/macros/plot.C+");


  TString cut = "pt1>40.0 && pt2>30.0";


  //plot from TTree named hHttNtuple
  //  plot("hHttNtuple","di-photon mass [GeV/c^{2}]",  0, 100.,  200., 0., -1.,1,lumi,"hmass",cut,100);

  //plot from TH1D
  //  plot("hDAllEvents","Events",                           0, 0.,  0., 0., 0.,1,lumi);

  plot("hTauPlay_NGenTaus",  "NGenTaus;N Gen Taus;",                             1, 0.,  0., 0., 0.,1,lumi);
  plot("hTauPlay_NPFTaus",   "NPFTaus;N PF Taus;",            	                 1, 0.,  0., 0., 0.,1,lumi);
  plot("hTauPlay_NHPSTaus",   "NHPSTaus;N HPS Taus;",         	                 1, 0.,  0., 0., 0.,1,lumi);
  plot("hTauPlay_PtGenTau",  "PtGenTau;Gen Tau Pt [GeV];",    	                 1, 0.,  0., 0., 0.,1,lumi);
  plot("hTauPlay_PtPFTau",   "PtPFTau;PF Tau Pt [GeV];",      	                 1, 0.,  0., 0., 0.,1,lumi);
  plot("hTauPlay_PtHPSTau",   "PtHPSTau;HPS Tau Pt [GeV];",   	                 1, 0.,  0., 0., 0.,1,lumi);
  plot("hTauPlay_PtDiff",    "PtDiff;#Delta Pt HPS-Gen Taus;",	                 1, 0.,  0., 0., 0.,1,lumi);
  plot("hTauPlay_Eff",       "Efficiency;NHPSTaus/NGenTaus;",                    1, 0.,  0., 0., 0.,1,lumi);    
  plot("hWWTauRepeat_PtLepton_pre",       "PtLepton_pre;Pt [GeV];",              1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_EtMet_pre",          "EtMet_pre;Et [GeV];",                 1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_PtHPSTau_pre",       "PtHPSTau_pre;Pt [GeV];",              1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_DeltaPhiLepTau_pre", "DeltaPhiLepTau_pre;#Delta #phi;",     1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_PtLepton_NMinusOne", "PtLepton_NMinusOne;Pt [GeV];",        1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_EtMet_NMinusOne",    "EtMet_NMinusOne;Et [GeV];",           1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_PtHPSTau_NMinusOne", "PtHPSTau_NMinusOne;Pt [GeV];",        1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_PtLepton_after",     "PtLepton_after;Pt [GeV];",            1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_EtMet_after",        "EtMet_after;Et [GeV];",               1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_PtHPSTau_after",     "PtHPSTau_after;Pt [GeV];lx",          1, 0.,  0., 0., 0.,2,lumi);
  plot("hWWTauRepeat_DeltaPhiLepTau_after", "DeltaPhiLepTau_after;#Delta #phi;", 1, 0.,  0., 0., 0.,2,lumi); 

  return;
}
