// $Id: EfficiencyUtils.cc,v 1.11 2009/11/04 12:44:33 loizides Exp $

#include "EfficiencyUtils.hh"

//--------------------------------------------------------------------------------------------------
// Create Efficiency Histogram. 
//--------------------------------------------------------------------------------------------------
//--------------------------------------------------------------------------------------------------
// Create Efficiency Graph. Use RooStatsCms Errors. 
// errorType == 0 : ROOT BayesDivide
// errorType == 1 : Feldman Cousins Confidence Intervals
// errorType == 2 : Clopper Pearson Confidence Intervals
//--------------------------------------------------------------------------------------------------
TGraphAsymmErrors* EfficiencyUtils::createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
                                                    string histname, 
                                                    vector<Double_t> bins, Int_t errorType, 
                                                    Double_t xlow, Double_t xhigh, 
                                                    Double_t ylow, Double_t yhigh 
  ) {

  TH1F *n = PlotUtils::rebin(numerator,bins);
  TH1F *d = PlotUtils::rebin(denominator,bins);
  
  Int_t nbins = n->GetNbinsX();

  assert(nbins <= 200);
  Double_t x[200];
  Double_t y[200];
  Double_t xErr[200];
  Double_t yErrLow[200];
  Double_t yErrHigh[200];

  for (int i=0;i < 200; i++) {
    x[i] = 0;
    y[i] = 0;
    xErr[i] = 0;
    yErrLow[i] = 0;
    yErrHigh[i] = 0;
  }

  //don't take the overflow bins
  for (int b=0; b<nbins ; ++b) {

    x[b] = n->GetXaxis()->GetBinCenter(b+1);    
    xErr[b] = 0.0;

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     

    Double_t n1 = TMath::Nint(n->GetBinContent(b+1));
    Double_t n2 = TMath::Nint(d->GetBinContent(b+1));
    if (n1 > n2) n1 = n2;
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, errorType);

//     cerr << " done bin " << b << " " << x[b] << " : " << n1 << "(" << n->GetBinContent(b+1) << ")" << " / " << n2 << "(" << d->GetBinContent(b+1) << ")" << " = " << ratio << " " << errLow << " " << errHigh << endl;
    y[b] = ratio;
    yErrLow[b] = errLow;
    yErrHigh[b] = errHigh;
  }
  TGraphAsymmErrors *efficiency = new TGraphAsymmErrors(nbins, x, y, xErr, xErr, yErrLow,yErrHigh );
  efficiency->SetName(histname.c_str());
  efficiency->SetTitle(histname.c_str());
  efficiency->GetXaxis()->SetTitle(numerator->GetXaxis()->GetTitle());
  efficiency->GetYaxis()->SetTitle("Efficiency");

  if (yhigh != -99)
    efficiency->SetMaximum(yhigh);
  if (ylow != -99)
    efficiency->SetMinimum(ylow);
  if (xlow != -99 && xhigh != -99) 
    efficiency->GetXaxis()->SetRangeUser(xlow,xhigh);

  efficiency->SetMarkerSize(1);
  efficiency->SetLineWidth(2);

  return efficiency;
}


//--------------------------------------------------------------------------------------------------
// Create Efficiency Graph. Use RooStatsCms Errors. 
// errorType == 0 : ROOT BayesDivide
// errorType == 1 : Feldman Cousins Confidence Intervals
// errorType == 2 : Clopper Pearson Confidence Intervals
//--------------------------------------------------------------------------------------------------
TH2DAsymErr* EfficiencyUtils::createEfficiencyHist2D(TH2F* numerator, TH2F* denominator,
                                           string histname, 
                                           vector<Double_t> xbins, vector<Double_t> ybins, 
                                           Int_t errorType, Bool_t printDebug) {
  
  TH2F *n = numerator;
  TH2F *d = denominator;
  if (xbins.size() > 0 && ybins.size() > 0) {
    n = PlotUtils::rebin(numerator,xbins,ybins);
    d = PlotUtils::rebin(denominator,xbins,ybins);
  }
  
  TH2D *Efficiency = (TH2D*)n->Clone(histname.c_str());
  Efficiency->Divide(n, d, 1.0,1.0,"B");  
  mithep::TH2DAsymErr *result = new mithep::TH2DAsymErr(*Efficiency);

  for (int b=1; b<Efficiency->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<Efficiency->GetYaxis()->GetNbins()+1 ; ++c) {

      if (result->GetBinContent(b,c) > 1.0) result->SetBinContent(b,c, 1.0);
      Double_t num = TMath::Nint(n->GetBinContent(b,c));
      Double_t den = TMath::Nint(d->GetBinContent(b,c));
      Double_t ratio = 0.0;
      Double_t errLow = 0.0;
      Double_t errHigh = 0.0;
      mithep::MathUtils::CalcRatio(num , den, ratio, errLow, errHigh, errorType);      
      result->SetBinError(b,c,errLow,errHigh,0.0,0.0);
      result->SetBinError(b,c,errLow);

      if (printDebug) 
        cout << b << "," << c << " : " << num << " / " << den << " = " << ratio << " | " << result->GetBinStatErrorLow(b,c) << " " << result->GetBinStatErrorHigh(b,c) <<  " .. " << errLow << " " << errHigh << endl;
    }
  }
  return result;
}


//--------------------------------------------------------------------------------------------------
// Create Efficiency Graph. Use RooStatsCms Errors. 
// errorType == 0 : ROOT BayesDivide
// errorType == 1 : Feldman Cousins Confidence Intervals
// errorType == 2 : Clopper Pearson Confidence Intervals
//--------------------------------------------------------------------------------------------------
void EfficiencyUtils::createEfficiencyHist2D(TH2F* numerator, TH2F* denominator,
                                             string histname, 
                                             vector<Double_t> xbins, vector<Double_t> ybins, 
                                             Int_t errorType, TFile *file ) {
  
  TH2F *n = numerator;
  TH2F *d = denominator;
  if (xbins.size() > 0 && ybins.size() > 0) {
    n = PlotUtils::rebin(numerator,xbins,ybins);
    d = PlotUtils::rebin(denominator,xbins,ybins);
  }
  
  TH2F *eff = (TH2F*)n->Clone(histname.c_str());
  assert(eff);
  eff->Divide(n, d, 1.0,1.0,"B");  

  for (int b=1; b<eff->GetXaxis()->GetNbins()+1 ; ++b) {
    for (int c=1; c<eff->GetYaxis()->GetNbins()+1 ; ++c) {

      if (eff->GetBinContent(b,c) > 1.0) eff->SetBinContent(b,c, 1.0);
      Double_t num = TMath::Nint(n->GetBinContent(b,c));
      Double_t den = TMath::Nint(d->GetBinContent(b,c));
      Double_t ratio = 0.0;
      Double_t errLow = 0.0;
      Double_t errHigh = 0.0;
      mithep::MathUtils::CalcRatio(num , den, ratio, errLow, errHigh, errorType);      
      
//       eff->SetBinContent(b,c,ratio);
//       eff->SetBinError(b,c,(errLow+errHigh)/2);

      cout << b << "," << c << " : " << num << " / " << den << " = " << ratio << " | " << errLow << " " << errHigh << " " << (errLow+errHigh)/2 << " : " << eff->GetBinError(b,c) << endl;

    }
  }

  file->WriteTObject(eff, eff->GetName(), "WriteDelete");

  return ;
}


//--------------------------------------------------------------------------------------------------
// Create Efficiency Graph. Use RooStatsCms Errors. 
// errorType == 0 : ROOT BayesDivide
// errorType == 1 : Feldman Cousins Confidence Intervals
// errorType == 2 : Clopper Pearson Confidence Intervals
//--------------------------------------------------------------------------------------------------
TH3DAsymErr* EfficiencyUtils::createEfficiencyHist3D(TH3F* numerator, TH3F* denominator,
                                               string histname, 
                                               vector<Double_t> xbins, vector<Double_t> ybins, 
                                               vector<Double_t> zbins, 
                                               Int_t errorType) {
    
  TH3F *n = numerator;
  TH3F *d = denominator;
  if (xbins.size() > 0 && ybins.size() > 0 && zbins.size() > 0) {
    n = PlotUtils::rebin(numerator,xbins,ybins,zbins);
    d = PlotUtils::rebin(denominator,xbins,ybins,zbins);
  }
  
  TH3D *Efficiency = (TH3D*)n->Clone(histname.c_str());
  Efficiency->Divide(n, d, 1.0,1.0,"B");  
  mithep::TH3DAsymErr *result = new mithep::TH3DAsymErr(*Efficiency);

  for (int a=1; a<Efficiency->GetXaxis()->GetNbins()+1 ; ++a) {
    for (int b=1; b<Efficiency->GetYaxis()->GetNbins()+1 ; ++b) {
      for (int c=1; c<Efficiency->GetZaxis()->GetNbins()+1 ; ++c) {

        if (result->GetBinContent(a,b,c) > 1.0) result->SetBinContent(a,b,c,1.0);
        Double_t num = TMath::Nint(n->GetBinContent(a,b,c));
        Double_t den = TMath::Nint(d->GetBinContent(a,b,c));
        Double_t ratio = 0.0;
        Double_t errLow = 0.0;
        Double_t errHigh = 0.0;
        mithep::MathUtils::CalcRatio(num , den, ratio, errLow, errHigh, errorType);      
        result->SetBinError(a,b,c,errLow,errHigh,0.0,0.0);
        result->SetBinError(a,b,c,errLow);
        
//        cout << a << "," << b << "," << c << " : " << num << " / " << den << " = " << ratio << " | " << result->GetBinStatErrorLow(a,b,c) << " " << result->GetBinStatErrorHigh(a,b,c) <<  " .. " << errLow << " " << errHigh << endl;
      }
    }
  }
  return result;
}
