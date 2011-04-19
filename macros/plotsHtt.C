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

  plot("hTauPlay_TauPt","",                               1, 0.,  0., 0., 0.,10,lumi);
  plot("hTauPlay_TauEta","",                             0, 0.,  0., 0., 0.,10,lumi);

  return;
}
