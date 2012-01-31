//--------------------------------------------------------------------------------------------------
// $Id $
//
// PlotUtils
//
// Utility functions for plotting
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef PLOTUTILS_HH
#define PLOTUTILS_HH

#include "MitCommon/DataFormats/interface/Types.h"
#include "MitAna/Utils/interface/SimpleTable.h"
#include <TError.h>
#include <TFile.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TF1.h>
#include <TF2.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <THStack.h>
#include <TPaveText.h>
#include <TGraphAsymmErrors.h>
#include "MitHtt/Common/MitStyleRemix.hh"
#include "MitAna/Utils/interface/SimpleTable.h"
#include "MitCommon/MathTools/interface/MathUtils.h"
#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include "MitCommon/DataFormats/interface/TH3DAsymErr.h"
#include <string>
#include <vector>
#include <utility>
#include <map>

using namespace mithep;
using namespace std;


  class PlotUtils {

    public:
      PlotUtils();
      PlotUtils(Double_t x);
      ~PlotUtils();

      void Init(Double_t x);
      void SetLuminosity(Double_t x)         { fIntegratedLuminosity = x ; }
      vector<Int_t> GetCOLORS()         { return fCOLORS;             }
      vector<Int_t> GetSYSCOLORS()      { return fSYSCOLORS;          }
      vector<Int_t> GetMARKERS()        { return fMARKERS;            }

      Double_t getWeight(string datasetFile, string datasetName);
      static TH1F *getHisto(string filename, string directoryname, string histoname);
      static TH2F* get2DHisto(string filename, string directoryname, string histoname);
      static TH1F* rebin(TH1F* hist, vector<Double_t> xlowedges);
      static TH1F* rebin(TH1F* hist, Int_t nbins);
      static TH2F* rebin(TH2F* hist, vector<Double_t> xlowedges, vector<Double_t> ylowedges);
      static TH3F* rebin(TH3F* hist, vector<Double_t> xlowedges, vector<Double_t> ylowedges, 
                  vector<Double_t> zlowedges);
      TH1F* addAllSamples(vector<string> datasetFiles, vector<string> datasetNames,
                    string dirName, string histName, vector<Double_t> bins = vector<Double_t>(0),
                    string xAxisTitle = "",string yAxisTitle = "");
      TH2F* addAllSamples2D(vector<string> datasetFiles, 
                            vector<string> datasetNames,
                            string dirName, string histName, 
                            vector<Double_t> xbins, 
                            vector<Double_t> ybins);
      THStack* addAllSamplesStacked(vector<string> datasetFiles, 
                                    vector<string> datasetNames, 
                                    string dirName, string histName, 
                                    vector<Double_t> bins = vector<Double_t>(0), 
                                    string xAxisTitle = "",
                                    string yAxisTitle = "");
      THStack* addAllSamplesStacked(vector<vector<string> > datasetFiles, 
                                    vector<vector<string> > datasetNames, 
                                    string dirName, string histName, 
                                    vector<Double_t> bins = vector<Double_t>(0), 
                                    string xAxisTitle = "",
                                    string yAxisTitle = "" );
      TH1F* createEfficiencyHist(vector<string> datasetFiles, 
                                 vector<string> datasetNames, 
                                 string dirName,string numeratorHistname, 
                                 string denominatorHistname, 
                                 string histname, vector<Double_t> bins, 
                                 Double_t xlow = -99, Double_t xhigh = -99, 
                                 Double_t ylow = -99, Double_t yhigh = -99);
      TGraphAsymmErrors* createEfficiencyGraph(vector<string> datasetFiles, 
                                             vector<string> datasetNames, 
                                             string dirName,string numeratorHistname, 
                                             string denominatorHistname, 
                                             string histname, vector<Double_t> bins, 
                                             Int_t errorType,
                                             Double_t xlow = -99, Double_t xhigh = -99, 
                                             Double_t ylow = -99, Double_t yhigh = -99);
      TGraphAsymmErrors* createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
                                               string histname, vector<Double_t> bins, 
                                               Int_t errorType = 0,
                                               Double_t xlow = -99, Double_t xhigh = -99, 
                                               Double_t ylow = -99, Double_t yhigh = -99);
     

      void drawStackedPlot(THStack *stackedHist , string plotname, 
                                  vector<string> legendNames,
                                  Bool_t logY = false, Double_t MaxY = -99, Double_t MinX = -99, Double_t MaxX = -99,
                                  Double_t legendX1 = -99, Double_t legendY1 = -99, 
                                  Double_t legendX2 = -99, Double_t legendY2 = -99,
                                  TH1F *hist = NULL, string histLegendLabel = "");
      TCanvas* DrawDataSignalBkgHistogram(TH1F* data, TH1F* sig, THStack *bkg , TLegend *legend, 
                                          string plotname, 
                                          Bool_t useLogY = false, Double_t MaxY = -99, 
                                          Double_t MinX = -99, Double_t MaxX = -99,
                                          Double_t legendX1 = -99, Double_t legendY1 = -99, 
                                          Double_t legendX2 = -99, Double_t legendY2 = -99);
      void makeDistributionComparisonPlot( vector<string> datasetfiles, 
                                           vector<string> datasetnames, 
                                           string dirName,
                                           vector<string> histNames, 
                                           vector<string> legendNames, 
                                           bool normalizeArea,
                                           Double_t MinX, Double_t MaxY,
                                           int nbins, string plotname );
      void makeXSliceDistributionComparisonPlot( vector<TH2F*> hists,                                       
                                                 vector<string> legendNames, 
                                                 bool normalizeArea,
                                                 string xAxisLabel = "",
                                                 string yAxisLabel = "",
                                                 Double_t xlow = -99, Double_t xhigh = -99, 
                                                 Double_t ylow = -99, Double_t yhigh = -99, 
                                                 Double_t legendX1 = -99, Double_t legendX2 = -99, 
                                                 Double_t legendY1 = -99, Double_t legendY2 = -99,
                                                 string plotname = "noNamePlot");
        void makeXSliceDistributionComparisonPlot( vector<string> datasetfiles, 
                                                   vector<string> datasetnames, 
                                                   string dirName,
                                                   vector<string> histNames, 
                                                   vector<string> legendNames, 
                                                   bool normalizeArea,
                                                   vector<Double_t> xbins, vector<Double_t> ybins,
                                                   string xAxisLabel = "",
                                                   string yAxisLabel = "",
                                                   Double_t xlow = -99, Double_t xhigh = -99, 
                                                   Double_t ylow = -99, Double_t yhigh = -99, 
                                                   Double_t legendX1 = -99, Double_t legendX2 = -99, 
                                                   Double_t legendY1 = -99, Double_t legendY2 = -99, 
                                                   string plotname = "noNamePlot");
      void makeYSliceDistributionComparisonPlot( vector<TH2F*> hists,     
                                                 vector<string> legendNames, 
                                                 bool normalizeArea,
                                                 string xAxisLabel = "",
                                                 string yAxisLabel = "",
                                                 Double_t xlow = -99, Double_t xhigh = -99, 
                                                 Double_t ylow = -99, Double_t yhigh = -99, 
                                                 Double_t legendX1 = -99, Double_t legendX2 = -99, 
                                                 Double_t legendY1 = -99, Double_t legendY2 = -99,
                                                 string plotname = "noNamePlot");
      void makeYSliceDistributionComparisonPlot( vector<string> datasetfiles, 
                                                 vector<string> datasetnames, 
                                                 string dirName,
                                                 vector<string> histNames, 
                                                 vector<string> legendNames, 
                                                 bool normalizeArea,
                                                 vector<Double_t> xbins, vector<Double_t> ybins,
                                                 string xAxisLabel = "",
                                                 string yAxisLabel = "",
                                                 Double_t xlow = -99, Double_t xhigh = -99, 
                                                 Double_t ylow = -99, Double_t yhigh = -99, 
                                                 Double_t legendX1 = -99, Double_t legendX2 = -99, 
                                                 Double_t legendY1 = -99, Double_t legendY2 = -99,
                                                 string plotname = "noNamePlot");
      void makeCrossDatasetComparisonPlot( vector<string> dataset1files, 
                                           vector<string> dataset1names, 
                                           string dataset1label,
                                           vector<string> dataset2files, 
                                           vector<string> dataset2names, 
                                           string dataset2label,
                                           string dirname1,
                                           vector<string> histNames1, 
                                           vector<string> legendNames1,
                                           string dirname2,
                                           vector<string> histNames2, 
                                           vector<string> legendNames2, 
                                           string plotname,
                                           bool normalizeArea = false, 
                                           string xAxisLabel = "",
                                           string yAxisLabel = "",
					   Double_t xlow = -99, Double_t xhigh = -99, 
                                           Double_t ylow = -99, Double_t yhigh = -99, 
                                           Double_t legendX1 = -99, Double_t legendX2 = -99, 
                                           Double_t legendY1 = -99, Double_t legendY2 = -99, 
                                           Bool_t useLogY= false, Int_t nbins = -1);
      void makeCrossDirComparisonPlot( vector<string> datasetfiles, 
                                       vector<string> datasetnames, 
                                       vector<string> dirnames, 
                                       vector<string> dirnamelabel, 
                                       vector<string> histNames,
                                       vector<string> legendNames, 
                                       vector<Double_t> bins, string plotname , 
                                       string xAxisLabel = "",
                                       string yAxisLabel = "",
                                       Double_t xlow = -99, Double_t xhigh = -99, 
                                       Double_t ylow = -99, Double_t yhigh = -99, 
                                       Double_t legendX1 = -99, Double_t legendX2 = -99, 
                                       Double_t legendY1 = -99, Double_t legendY2 = -99, 
                                       Bool_t useLogY = false);
      void makeComparisonPlotWithSystematics( vector<string> datasetfiles, 
                                              vector<string> datasetnames, 
                                              vector<string> dirnames, 
                                              vector<string> dirnamelabel, 
                                              vector<string> histNames, 
                                              vector<string> errorhistNames, 
                                              vector<string> legendNames, 
                                              vector<Double_t> bins, 
                                              string plotname,
                                              string xAxisLabel = "",
                                              string yAxisLabel = "",
                                              Double_t xlow = -99, Double_t xhigh = -99, 
                                              Double_t ylow = -99, Double_t yhigh = -99, 
                                              Double_t legendX1 = -99, Double_t legendX2 = -99, 
                                              Double_t legendY1 = -99, Double_t legendY2 = -99, 
                                              Bool_t useLogY = false);
   
      static void NormalizeHist(TH1F *hist);
      static TGraphAsymmErrors* MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, 
                                                        string name, 
                                                        Bool_t cutBelow = kTRUE );
      static TGraphAsymmErrors* MakeCurrentWPSigEffVsBkgEffGraph(Double_t signalEff, Double_t bkgEff, 
                                                                 string name );
      static TGraphAsymmErrors* MakeCurrentWPSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, 
                                                                 string name, Double_t myCutValue, 
                                                                 Bool_t cutBelow = kTRUE);
      static TGraphAsymmErrors* MakeSigEffVsCutValueGraph(TH1F* signalHist, string name , 
                                                          Bool_t cutBelow = kTRUE);
      static TGraphAsymmErrors* MakeCurrentWPSigEffVsCutValueGraph(TH1F* signalHist, string name, 
                                                                   Double_t myCutValue, 
                                                                   Bool_t cutBelow  = kTRUE);
      static Double_t FindCutValueAtFixedSignalEfficiency(TH1F* signalHist, Double_t targetSignalEff, 
                                                          Bool_t cutBelow = kTRUE );
      static Double_t FindBkgEffAtFixedSignalEfficiency(TH1F* signalHist, TH1F* bkgHist, 
                                                        Double_t targetSignalEff, 
                                                        Bool_t cutBelow = kTRUE );
      

    private:
      vector<Int_t> fCOLORS;
      vector<Int_t> fMARKERS; 
      vector<Int_t> fSYSCOLORS;
      Double_t fIntegratedLuminosity;
      SimpleTable *xstab;

  };

#endif
