//--------------------------------------------------------------------------------------------------
// $Id: EfficiencyUtils.h,v 1.12 2009/11/03 10:27:41 sixie Exp $
//
// EfficiencyUtils
//
// Math utility functions.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef EFFICIENCYUTILS_HH
#define EFFICIENCYUTILS_HH

#include "PlotUtils.hh"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TError.h>
#include <TMath.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitCommon/DataFormats/interface/TH3DAsymErr.h"
#include "PhysicsTools/RooStatsCms/interface/ClopperPearsonBinomialInterval.h"
#include "PhysicsTools/RooStatsCms/interface/FeldmanCousinsBinomialInterval.h"
#include <string>
#include <vector>
#include <utility>

using namespace mithep;
using namespace std;

  class EfficiencyUtils {
    public:
      static TGraphAsymmErrors* createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
                                               string histname, vector<Double_t> bins, 
                                               Int_t errorType = 2,
                                               Double_t xlow = -99, Double_t xhigh = -99, 
                                               Double_t ylow = -99, Double_t yhigh = -99);
      static TH2DAsymErr* createEfficiencyHist2D(TH2F* numerator, TH2F* denominator,
                                          string histname, 
                                          vector<Double_t> xbins, vector<Double_t> ybins, 
                                                 Int_t errorType = 0, Bool_t printDebug = kFALSE);
      static void createEfficiencyHist2D(TH2F* numerator, TH2F* denominator, 
                                         string histname, 
                                         vector<Double_t> xbins, vector<Double_t> ybins, 
                                         Int_t errorType, TFile *file );
       static TH3DAsymErr* createEfficiencyHist3D(TH3F* numerator, TH3F* denominator,
                                          string histname, 
                                          vector<Double_t> xbins, vector<Double_t> ybins, 
                                          vector<Double_t> zbins, 
                                          Int_t errorType = 0);
 
  };


#endif
