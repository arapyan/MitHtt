// $Id $

#include "PlotUtils.hh"
 

PlotUtils::PlotUtils() : xstab(0) {
  //Constructor. Fill with default color scheme

  Init(1.0);
}

PlotUtils::PlotUtils(Double_t x) : xstab(0) {
  //Constructor. Fill with specified integrated luminosity

  Init(x);
}

PlotUtils::~PlotUtils() {
  //Destructor

  delete xstab;
}              

void PlotUtils::Init( Double_t lumi ) {
  //Constructor. Fill with default color scheme

  //Load Cross-section Table
  xstab = new SimpleTable("$CMSSW_BASE/src/MitPhysics/data/xs.dat");

  //Set Default Color Scheme
  fCOLORS.push_back(kRed);   fCOLORS.push_back(kBlue);   fCOLORS.push_back(kMagenta);  
  fCOLORS.push_back(kGreen); fCOLORS.push_back(kOrange);   fCOLORS.push_back(kCyan);   
  fCOLORS.push_back(kGreen+3); fCOLORS.push_back(kRed+3);   fCOLORS.push_back(kBlue+3);  
  fCOLORS.push_back(kGray);    
  fMARKERS.push_back(20);  fMARKERS.push_back(21); fMARKERS.push_back(22); fMARKERS.push_back(23);  
  fMARKERS.push_back(24);  fMARKERS.push_back(25); fMARKERS.push_back(26); fMARKERS.push_back(27);  
  fMARKERS.push_back(28);  fMARKERS.push_back(30);
  fSYSCOLORS.push_back(kMagenta);  fSYSCOLORS.push_back(kRed);   fSYSCOLORS.push_back(kBlue);  
  fSYSCOLORS.push_back(kGreen); fSYSCOLORS.push_back(kOrange);   fSYSCOLORS.push_back(kCyan);  
  fSYSCOLORS.push_back(kGreen+3); fSYSCOLORS.push_back(kRed+3);   fSYSCOLORS.push_back(kBlue+3);  
  fSYSCOLORS.push_back(kGray);  

  fIntegratedLuminosity = lumi;
  return;
}


//--------------------------------------------------------------------------------------------------
// Get 1D Histogram function
//--------------------------------------------------------------------------------------------------
TH1F* PlotUtils::getHisto(string filename, string directoryname, string histoname) {
  // Get 1D Histogram from file

  TFile *file = new TFile(filename.c_str(),"READ");
  if (!file) {
    cout << "Could not open file " << filename << endl;
    return 0;
  }

  TDirectory *dir = (TDirectory*)file->FindObjectAny(directoryname.c_str());
  if (!dir) {
    cout << "Could not find directory " << directoryname 
         << " in file " << filename << endl;
    delete file;
    return 0;
  }

  TH1F *tmphist = (TH1F*)dir->Get(histoname.c_str());
  if (!tmphist) {
    cout << "Could not find histogram " <<  histoname << " in directory " << directoryname 
         << " in file " << filename << endl;
    delete dir;
    delete file;
    return 0;
  }
  TH1F *hist = (TH1F*)tmphist->Clone(histoname.c_str());
  hist->SetDirectory(0);

  delete tmphist;
  delete dir;
  delete file;
  return hist;

}

//--------------------------------------------------------------------------------------------------
// Get 2D Histogram function
//--------------------------------------------------------------------------------------------------
TH2F* PlotUtils::get2DHisto(string filename, string directoryname, string histoname) {
  // Get 2D Histogram from file

  TFile *file = new TFile(filename.c_str(),"READ");
  if (!file) {
    cout << "Could not open file " << filename << endl;
    return 0;
  }

  TDirectory *dir = (TDirectory*)file->FindObjectAny(directoryname.c_str());
  if (!dir) {
    cout << "Could not find directory " << directoryname 
         << " in file " << filename << endl;
    delete file;
    return 0;
  }

  TH2F *tmphist = (TH2F*)dir->Get(histoname.c_str());
  if (!tmphist) {
    cout << "Could not find histogram " <<  histoname << " in directory " << directoryname 
         << " in file " << filename << endl;
    delete dir;
    delete file;
    return 0;
  }

  TH2F *hist = (TH2F*)tmphist->Clone(histoname.c_str());
  hist->SetDirectory(0);

  delete tmphist;
  delete dir;
  delete file;
  return hist;

}

//--------------------------------------------------------------------------------------------------
// Rebin for 1D hists
//--------------------------------------------------------------------------------------------------
TH1F* PlotUtils::rebin(TH1F* hist, vector<Double_t> xlowedges) {

  Double_t *xbins;
  Double_t *BinContent;
  Double_t *BinError;
 
  xbins = new Double_t[xlowedges.size()];
  for (UInt_t i=0;i<xlowedges.size();i++) {
    xbins[i] = xlowedges[i];
//     cout << "bins: " << xbins[i] << endl;
  }

//   cout << "Nbins: " << xlowedges.size() << endl;

  BinContent = new Double_t[xlowedges.size()+2];
  BinError = new Double_t[xlowedges.size()+2];
  for (UInt_t binx=0; binx < xlowedges.size()+2 ;++binx) {
      BinContent[binx] = 0;
      BinError[binx] = 0;
  }

  TH1F *rebinHist = new TH1F(hist->GetName(), hist->GetTitle(), xlowedges.size()-1, xbins);

//   cout << "new nbins : " << rebinHist->GetXaxis()->GetNbins() << endl;

  //refill
  for (UInt_t i=0;int(i)<hist->GetXaxis()->GetNbins()+2;++i) {
    Double_t x = hist->GetXaxis()->GetBinCenter(i);
    UInt_t xbin = 0;

//     cout << "old bin " << i << " : " << x << endl;
    //Find which x rebinned bin we are in
    for (UInt_t binx=0; binx < xlowedges.size()+1 ; ++binx) {
      if (binx == 0) { //underflow 
        if (x <= xlowedges[0]) {
          xbin = 0;
          break;
        }
      } else if (binx < xlowedges.size()) {
        if (x > xlowedges[binx-1] && x <= xlowedges[binx]) {
          xbin = binx;
          break;
        }
      } else { //overflow
        if (x > xlowedges[binx]) {
          xbin = binx;
          break;
        }
      }
    }
    BinContent[xbin] += hist->GetBinContent(i);
    BinError[xbin] += hist->GetBinError(i)*hist->GetBinError(i);
  }
  
  for (UInt_t binx=0; binx < xlowedges.size()+2 ;++binx) {
    rebinHist->SetBinContent(binx,BinContent[binx]);
    rebinHist->SetBinError(binx,TMath::Sqrt(BinError[binx]));     
  } 
  
  delete [] xbins;
  delete [] BinContent;
  delete [] BinError;
  return rebinHist;
}

//--------------------------------------------------------------------------------------------------
// Rebin histogram into NBins
//--------------------------------------------------------------------------------------------------
TH1F* PlotUtils::rebin(TH1F* hist, Int_t nbins) {
  TH1F *result = hist;
  if (nbins > 0) {

//     cout << "nbins = " << nbins << endl;

    vector<Double_t> bins;
    for (int b=0; b<nbins+1; ++b) {
      Double_t binsize = (hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetNbins()) - 
                          hist->GetXaxis()->GetBinLowEdge(1)) / nbins;
      bins.push_back(hist->GetXaxis()->GetBinLowEdge(1) + binsize * b);
//       cout << "binsize = " << binsize << " : " << hist->GetXaxis()->GetBinLowEdge(1) + binsize * b << endl;
    }
    result = rebin(hist, bins);
  }
  return result;
}

//--------------------------------------------------------------------------------------------------
// Rebin for 2D hists
//--------------------------------------------------------------------------------------------------
TH2F* PlotUtils::rebin(TH2F* hist, vector<Double_t> xlowedges, vector<Double_t> ylowedges) {

  Double_t *xLow;
  Double_t *yLow;
  Double_t **BinContent;
  Double_t **BinError;

  xLow = new Double_t[xlowedges.size()];
  yLow = new Double_t[ylowedges.size()];
  BinContent = new Double_t*[xlowedges.size()+1];
  BinError = new Double_t*[xlowedges.size()+1];
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {    
    BinContent[binx] = new Double_t[ylowedges.size()+1];
    BinError[binx] = new Double_t[ylowedges.size()+1];
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
      BinContent[binx][biny] = 0;
      BinError[binx][biny] = 0;
    }
  }

  for (UInt_t i=0; i<xlowedges.size();++i)
    xLow[i] = xlowedges[i];
  for (UInt_t i=0; i<ylowedges.size();++i) {
    yLow[i] = ylowedges[i];
  }

  TH2F *rebinHist = new TH2F(hist->GetName(), hist->GetTitle(), xlowedges.size() - 1, xLow, ylowedges.size() - 1, yLow);

  //refill the histogram
  for (UInt_t i=0;Int_t(i)<=hist->GetXaxis()->GetNbins()+1;i++) {
    for (UInt_t j=0;Int_t(j)<=hist->GetYaxis()->GetNbins()+1;j++) {
      

      Double_t x = hist->GetXaxis()->GetBinCenter(i);
      Double_t y = hist->GetYaxis()->GetBinCenter(j);
      UInt_t xbin = 0;
      UInt_t ybin = 0;

      for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {
        if (binx == 0) { //underflow 
          if (x <= xlowedges[0]) {
            xbin = 0;
            break;
          }
        } else if (binx < xlowedges.size()) {
          if (x > xlowedges[binx-1] && x <= xlowedges[binx]) {
            xbin = binx;
            break;
          }
        } else { //overflow
          if (x > xlowedges[binx]) {
            xbin = binx;
            break;
          }
        }
      }
      
      for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
        if (biny == 0) { //underflow 
          if (y <= ylowedges[0]) {
            ybin = 0;
            break;
          }
        } else if (biny < ylowedges.size()) {
          if (y > ylowedges[biny-1] && y <= ylowedges[biny]) {
            ybin = biny;
            break;
          }
        } else { //overflow
          if (y > ylowedges[biny]) {
            ybin = biny;
            break;
          }
        }
      }
            
      BinContent[xbin][ybin] += hist->GetBinContent(i,j);
      BinError[xbin][ybin] += hist->GetBinError(i,j)*hist->GetBinError(i,j);     
    }
  }

  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
      rebinHist->SetBinContent(binx,biny,BinContent[binx][biny]);
      rebinHist->SetBinError(binx,biny,TMath::Sqrt(BinError[binx][biny])); 
    }
  }

  delete [] xLow;
  delete [] yLow;
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {    
    delete [] BinContent[binx];
    delete [] BinError[binx];
  }
  delete [] BinContent;
  delete [] BinError;

  return rebinHist;
}


//--------------------------------------------------------------------------------------------------
// Rebin for 3D hists
//--------------------------------------------------------------------------------------------------
TH3F* PlotUtils::rebin(TH3F* hist, vector<Double_t> xlowedges, vector<Double_t> ylowedges, vector<Double_t> zlowedges) {

  Double_t *xLow;
  Double_t *yLow;
  Double_t *zLow;
  Double_t ***BinContent;
  Double_t ***BinError;

  xLow = new Double_t[xlowedges.size()];
  yLow = new Double_t[ylowedges.size()];
  zLow = new Double_t[zlowedges.size()];
  BinContent = new Double_t**[xlowedges.size()+1];
  BinError = new Double_t**[xlowedges.size()+1];
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {    
    BinContent[binx] = new Double_t*[ylowedges.size()+1];
    BinError[binx] = new Double_t*[ylowedges.size()+1];
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
      BinContent[binx][biny] = new Double_t[zlowedges.size()+1];
      BinError[binx][biny] = new Double_t[zlowedges.size()+1];
      for (UInt_t binz=0; binz < zlowedges.size()+1 ;++binz) {
        BinContent[binx][biny][binz] = 0;
        BinError[binx][biny][binz] = 0;
      }
    }
  }

  for (UInt_t i=0; i<xlowedges.size();++i)
    xLow[i] = xlowedges[i];
  for (UInt_t i=0; i<ylowedges.size();++i) {
    yLow[i] = ylowedges[i];
  }
  for (UInt_t i=0; i<zlowedges.size();++i) {
    zLow[i] = zlowedges[i];
  }

  TH3F *rebinHist = new TH3F(hist->GetName(), hist->GetTitle(), xlowedges.size() - 1, xLow, ylowedges.size() - 1, yLow, zlowedges.size() - 1, zLow);

  //refill the histogram
  for (UInt_t i=0;Int_t(i)<=hist->GetXaxis()->GetNbins()+1;++i) {
    for (UInt_t j=0;Int_t(j)<=hist->GetYaxis()->GetNbins()+1;++j) {
      for (UInt_t k=0;Int_t(k)<=hist->GetZaxis()->GetNbins()+1;++k) {
      
        Double_t x = hist->GetXaxis()->GetBinCenter(i);
        Double_t y = hist->GetYaxis()->GetBinCenter(j);
        Double_t z = hist->GetZaxis()->GetBinCenter(k);
        UInt_t xbin = 0;
        UInt_t ybin = 0;
        UInt_t zbin = 0;
        
        for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {
          if (binx == 0) { //underflow 
            if (x <= xlowedges[0]) {
              xbin = 0;
              break;
            }
          } else if (binx < xlowedges.size()) {
            if (x > xlowedges[binx-1] && x <= xlowedges[binx]) {
              xbin = binx;
              break;
            }
          } else { //overflow
            if (x > xlowedges[binx]) {
              xbin = binx;
              break;
            }
          }
        }
        
        for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
          if (biny == 0) { //underflow 
            if (y <= ylowedges[0]) {
              ybin = 0;
              break;
            }
          } else if (biny < ylowedges.size()) {
            if (y > ylowedges[biny-1] && y <= ylowedges[biny]) {
              ybin = biny;
              break;
            }
          } else { //overflow
            if (y > ylowedges[biny]) {
              ybin = biny;
              break;
            }
          }
        }
        
        for (UInt_t binz=0; binz < zlowedges.size()+1 ;++binz) {
          if (binz == 0) { //underflow 
            if (z <= zlowedges[0]) {
              zbin = 0;
              break;
            }
          } else if (binz < zlowedges.size()) {
            if (z > zlowedges[binz-1] && z <= zlowedges[binz]) {
              zbin = binz;
              break;
            }
          } else { //overflow
            if (z > zlowedges[binz]) {
              zbin = binz;
              break;
            }
          }
        }
        
        BinContent[xbin][ybin][zbin] += hist->GetBinContent(i,j,k);
        BinError[xbin][ybin][zbin] += hist->GetBinError(i,j,k)*hist->GetBinError(i,j,k);
      }
    }
  }
  
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {
      for (UInt_t binz=0; binz < zlowedges.size()+1 ;++binz) {
        rebinHist->SetBinContent(binx,biny,binz,BinContent[binx][biny][binz]);
        rebinHist->SetBinError(binx,biny,binz,TMath::Sqrt(BinError[binx][biny][binz])); 
      }
    }
  }
  
  delete [] xLow;
  delete [] yLow;
  delete [] zLow;
  for (UInt_t binx=0; binx < xlowedges.size()+1 ;++binx) {    
    for (UInt_t biny=0; biny < ylowedges.size()+1 ;++biny) {    
      delete [] BinContent[binx][biny];
      delete [] BinError[binx][biny];
    }
    delete [] BinContent[binx];
    delete [] BinError[binx];
  }
  delete [] BinContent;
  delete [] BinError;
  
  return rebinHist;
}





Double_t PlotUtils::getWeight(string datasetFile, string datasetName) {

  Double_t weight = 1.0;
  double CrossSection = xstab->Get(datasetName.c_str());
  if (datasetName == "oct09-qcd15-hjl-sdj30") {
    weight = 1400.510000;
  } else if (datasetName == "oct09-qcd30-hjl-sdj30" || datasetName == "oct09-qcd30-hjl-sdj50" ) {
    weight = 115.014000;
  } else if (datasetName == "oct09-qcd80-hjl-sdj30" || datasetName == "oct09-qcd80-hjl-sdj50" ) {
    weight = 3.770120;
  }else if (datasetName == "minbias" ) {
    weight = 7.11332e+06 / 1.09635e+07;
  }
  else if (CrossSection > 0) {
    TH1F *tmpNEventsHist = getHisto(datasetFile,"AnaFwkMod", "hDEvents");
    double NEvents = tmpNEventsHist->Integral();
    
    //Oct Ex Skims have different input number of events
    if (datasetName == "oct09-qcdem20_30-7-hww-sde10eid") {
      NEvents = 33505929;
    }
    if (datasetName == "oct09-qcdem30_80-7-hww-sde10eid") {
      NEvents = 32168675;
    }
    if (datasetName == "oct09-qcdem80_170-7-hww-sde10eid") {
      NEvents = 5551386;
    }
    if (datasetName == "oct09-qcdbc20_30-7-hww-sde10eid") {
      NEvents = 2752942;
    }
    if (datasetName == "oct09-qcdbc30_80-7-hww-sde10eid") {
      NEvents = 2261916;
    }
    if (datasetName == "oct09-qcdbc80_170-7-hww-sde10eid") {
      NEvents = 1097829;
    }
    if (datasetName == "oct09-qcdem20_30-hww-sde10eid") {
      NEvents = 33638282;
    }
    if (datasetName == "oct09-qcdem30_80-hww-sde10eid") {
      NEvents = 38360886;
    }
    if (datasetName == "oct09-qcdem80_170-hww-sde10eid") {
      NEvents = 5729547;
    }
    if (datasetName == "oct09-qcdbce20_30-hww-sde10eid") {
      NEvents = 2383833;
    }
    if (datasetName == "oct09-qcdbce30_80-hww-sde10eid") {
      NEvents = 2035108;
    }
    if (datasetName == "oct09-qcdbce80_170-hww-sde10eid") {
      NEvents = 1038080;
    }    
    weight = CrossSection * fIntegratedLuminosity / NEvents;
  }
  return weight;
}


//--------------------------------------------------------------------------------------------------
// Add histograms form multiple files together weighted by cross-section
//--------------------------------------------------------------------------------------------------
TH1F* PlotUtils::addAllSamples(vector<string> datasetFiles, vector<string> datasetNames,
                    string dirName, string histName, vector<Double_t> bins, 
                    string xAxisTitle,string yAxisTitle) {

  assert(datasetFiles.size() > 0);
  TH1F *tmp = getHisto(datasetFiles[0], dirName, histName);
  if (bins.size() > 0) tmp = rebin(tmp, bins);
  assert(tmp);
  TH1F *finalHist = (TH1F*)tmp->Clone();
  finalHist->Sumw2();

  if (xAxisTitle != "") finalHist->SetXTitle(xAxisTitle.c_str());
  if (yAxisTitle != "") finalHist->SetYTitle(yAxisTitle.c_str());
  
  for (UInt_t i=0; i < datasetFiles.size() ;i++) {
    double weight = getWeight(datasetFiles[i],datasetNames[i]);   
    if (i==0) {
      for (int b=0;b<=finalHist->GetNbinsX()+1;b++) {
        finalHist->SetBinContent(b,finalHist->GetBinContent(b)*weight);
        finalHist->SetBinError(b,finalHist->GetBinError(b)*weight);
      }
    } else {
      TH1F *tmpHist = getHisto(datasetFiles[i], dirName, histName);
      if (bins.size() > 0) tmpHist = rebin(tmpHist, bins);
      for (int b=0;b<=finalHist->GetNbinsX();b++) {
        tmpHist->SetBinContent(b,tmpHist->GetBinContent(b)*weight);
        tmpHist->SetBinError(b,tmpHist->GetBinError(b)*weight);
      }    
      finalHist->Add(tmpHist);
      delete tmpHist;
    }
  }
  return finalHist;
}

//--------------------------------------------------------------------------------------------------
// Add histograms form multiple files together weighted by cross-section
//--------------------------------------------------------------------------------------------------
TH2F* PlotUtils::addAllSamples2D(vector<string> datasetFiles, vector<string> datasetNames,
                    string dirName, string histName, vector<Double_t> xbins, vector<Double_t> ybins) {

  assert(datasetFiles.size() > 0);
  TH2F *tmp = get2DHisto(datasetFiles[0], dirName, histName);
  TH2F *tmpRebinned = tmp;  
  if (xbins.size() > 0 && ybins.size() > 0)
    tmpRebinned = rebin(tmp,xbins,ybins);
  
  assert(tmp);
  TH2F *finalHist = (TH2F*)tmpRebinned->Clone();
  finalHist->Sumw2();

  for (UInt_t i=0; i < datasetFiles.size() ;i++) {
    double weight = getWeight(datasetFiles[i],datasetNames[i]);   
    if (i==0) {
      for (int b=0;b<=finalHist->GetNbinsX()+1;b++) {
        for (int c=0;c<=finalHist->GetNbinsY()+1;c++) {
          finalHist->SetBinContent(b,c,finalHist->GetBinContent(b,c)*weight);
          finalHist->SetBinError(b,c,finalHist->GetBinError(b,c)*weight);
        }
      }
    } else {
      TH2F *tmpHist = get2DHisto(datasetFiles[i], dirName, histName);
      TH2F *tmpRebinnedHist = rebin(tmpHist,xbins,ybins);
      for (int b=0;b<=tmpRebinnedHist->GetNbinsX()+1;b++) {
        for (int c=0;c<=tmpRebinnedHist->GetNbinsY()+1;c++) {
          tmpRebinnedHist->SetBinContent(b,c,tmpRebinnedHist->GetBinContent(b,c)*weight);
          tmpRebinnedHist->SetBinError(b,c,tmpRebinnedHist->GetBinError(b,c)*weight);
        }
      }
      finalHist->Add(tmpRebinnedHist);
      delete tmpHist;
      delete tmpRebinnedHist;
    }
  }
  return finalHist;
}

//--------------------------------------------------------------------------------------------------
// Add histograms in a stack from multiple files together weighted by cross-section 
//--------------------------------------------------------------------------------------------------
THStack* PlotUtils::addAllSamplesStacked(vector<string> datasetFiles, 
                                         vector<string> datasetNames, 
                                         string dirName, string histName, 
                                         vector<Double_t> bins, 
                                         string xAxisTitle,
                                         string yAxisTitle) {

  assert(datasetFiles.size() > 0);
  assert(datasetFiles.size() == datasetNames.size());

  vector<vector<string> > tmpDatasetFiles;
  vector<vector<string> > tmpDatasetNames;
  for (UInt_t i=0; i<datasetFiles.size() ; ++i) {
    vector<string> tmp1;
    vector<string> tmp2;
    tmp1.push_back(datasetFiles[i]);
    tmp2.push_back(datasetNames[i]);    
    tmpDatasetFiles.push_back(tmp1);
    tmpDatasetNames.push_back(tmp2);
  }
  return addAllSamplesStacked(tmpDatasetFiles, tmpDatasetNames, dirName, histName, bins, 
                              xAxisTitle, yAxisTitle);
}
//--------------------------------------------------------------------------------------------------
// Add histograms in a stack from multiple files together weighted by cross-section 
//--------------------------------------------------------------------------------------------------
THStack* PlotUtils::addAllSamplesStacked(vector<vector<string> > datasetFiles, 
                                         vector<vector<string> > datasetNames, 
                                         string dirName, string histName, 
                                         vector<Double_t> bins, 
                                         string xAxisTitle,
                                         string yAxisTitle) {

  assert(datasetFiles.size() > 0);
  assert(datasetFiles[0].size() > 0);
  TH1F *tmp = getHisto(datasetFiles[0][0], dirName, histName);
  assert(tmp);

  THStack *stackedHist = new THStack(tmp->GetName(),tmp->GetName());

  for (UInt_t i=0; i < datasetFiles.size() ;++i) {

    TH1F *tmpTotalHist = 0;
    for (UInt_t j=0; j < datasetFiles[i].size() ;++j) {
    
      double weight = getWeight(datasetFiles[i][j],datasetNames[i][j]);         
      TH1F *tmpHist = getHisto(datasetFiles[i][j], dirName, histName);
      TH1F *tmpRebinnedHist = tmpHist;
      if (bins.size() > 0 ) tmpRebinnedHist = rebin(tmpHist, bins);
      
      if (xAxisTitle != "") tmpRebinnedHist->SetXTitle(xAxisTitle.c_str());
      if (yAxisTitle != "") tmpRebinnedHist->SetYTitle(yAxisTitle.c_str());

      if (tmpRebinnedHist) {
        tmpRebinnedHist->SetTitle(datasetNames[i][j].c_str());
        
        //scale by weight
        for (int b=0;b<=tmpRebinnedHist->GetNbinsX();b++) {
          tmpRebinnedHist->SetBinContent(b,tmpRebinnedHist->GetBinContent(b)*weight);
          tmpRebinnedHist->SetBinError(b,tmpRebinnedHist->GetBinError(b)*weight);
        }    
        
        if (j==0) {
          tmpTotalHist = (TH1F*)tmpRebinnedHist->Clone();
          tmpTotalHist->Sumw2();
        } else {
          tmpTotalHist->Add(tmpRebinnedHist);
        }
      } else {
        cerr << "could not get histogram " << datasetNames[i][j] << "\n";
      }
      delete tmpHist;
      delete tmpRebinnedHist;
    }
    //add to stack
    tmpTotalHist->SetFillStyle(1001);
    tmpTotalHist->SetFillColor(fCOLORS[i]);
    tmpTotalHist->SetLineWidth(1);  
    stackedHist->Add(tmpTotalHist);    
  }
  return stackedHist;
}

//--------------------------------------------------------------------------------------------------
// Draw a Stack
// hist is drawn separately (used to distinguish signal from bkg stack)
//--------------------------------------------------------------------------------------------------
void PlotUtils::drawStackedPlot(THStack *stackedHist , string plotname, vector<string> legendNames,
                     Bool_t logY, Double_t MaxY, Double_t MinX, Double_t MaxX,
                     Double_t legendX1, Double_t legendY1, 
                     Double_t legendX2, Double_t legendY2,
                     TH1F *hist, string histLegendLabel) {

  if (stackedHist->GetStack()->GetEntries() != int(legendNames.size())) {
    cerr << "Number of entries in the stack is not equal to the number of legend labels given\n";
    assert(stackedHist->GetStack()->GetEntries() == int(legendNames.size()));
  }

  TCanvas *cv = MakeCanvas("cv", plotname.c_str(), 800, 600);

  TLegend *leg1=0;
  if (legendX1 == -99) {
    leg1 = new TLegend(0.65,0.70,0.95,0.95);   
  } else {
    leg1 = new TLegend(legendX1,legendY1,legendX2,legendY2);   
  }
  leg1->SetBorderSize(1);  
  leg1->SetTextSize(0.03);

  for (int i=0; i < stackedHist->GetStack()->GetEntries(); ++i) {
    leg1->AddEntry(stackedHist->GetStack()->At(i),legendNames[i].c_str(), "F"); 
  }

  stackedHist->Draw("hist");

  stackedHist->SetMinimum(0.001);
  if (MaxY > 0) stackedHist->SetMaximum(MaxY);
  if (MinX != -99 && MaxX != -99) stackedHist->GetXaxis()->SetRangeUser(MinX, MaxX);

  stackedHist->GetXaxis()->SetTitle(
    ((TH1F*)(stackedHist->GetHists()->At(0)))->GetXaxis()->GetTitle());
  stackedHist->GetYaxis()->SetTitle(
    ((TH1F*)(stackedHist->GetHists()->At(0)))->GetYaxis()->GetTitle());

  if (hist) {
    leg1->AddEntry(hist,histLegendLabel.c_str(), "F"); 
    hist->SetLineColor(kBlack);
    hist->SetMarkerColor(kBlack);
    hist->Draw("same");
  }

  leg1->Draw();
  if (logY) cv->SetLogy();

  cv->SaveAs((plotname+".gif").c_str());
  cv->SaveAs((plotname+".eps").c_str());
  cv->Delete();
}

//--------------------------------------------------------------------------------------------------
// Draw data, signal, and bkg stack
//--------------------------------------------------------------------------------------------------
TCanvas* PlotUtils::DrawDataSignalBkgHistogram(TH1F* data, TH1F* sig, THStack *bkg , TLegend *legend, 
                                           string plotname, 
                                           Bool_t useLogY, Double_t MaxY, 
                                           Double_t MinX, Double_t MaxX,
                                           Double_t legendX1, Double_t legendY1, 
                                           Double_t legendX2, Double_t legendY2) {

   string filename = plotname;
  if (useLogY) filename += "_logY";
//  filename += ".eps";

  TCanvas *cv = new TCanvas(plotname.c_str(), plotname.c_str(), 0,0,800,600);
  if (useLogY) cv->SetLogy();

  sig->SetLineWidth(1);
  sig->SetLineColor(kAzure+2);
  sig->SetFillStyle(1001);
  sig->SetFillColor(kAzure+6);

  legend->SetTextSize(0.03);
  legend->SetBorderSize(1);

//   legend->AddEntry(data, "Data", "P");  
//   legend->AddEntry(sig, "Z #rightarrow ee", "F");
  
  THStack *tmpBkg = (THStack*)bkg->Clone("bkgStack");
  tmpBkg->Add(sig);

   sig->Scale(data->Integral() / sig->Integral());



  //set best Y-range
  
  Double_t maxY = sig->GetMaximum();
  if (data->GetMaximum() > maxY )
    maxY = data->GetMaximum();

  tmpBkg->Draw("hist");
  if (tmpBkg->GetMaximum() > maxY)
    maxY = tmpBkg->GetMaximum();
  maxY = maxY*1.2;
  if (MaxY > -99) maxY = MaxY;

  tmpBkg->GetYaxis()->SetTitleOffset(1.1);
  sig->GetYaxis()->SetTitleOffset(1.1);
  data->GetYaxis()->SetTitleOffset(1.1);
  tmpBkg->GetXaxis()->SetTitleOffset(1.0);
  sig->GetXaxis()->SetTitleOffset(1.0);
  data->GetXaxis()->SetTitleOffset(1.0);

  tmpBkg->SetMinimum(0.01);
  tmpBkg->SetMaximum(maxY);  
  sig->SetMinimum(0.01);
  sig->SetMaximum(maxY);
  data->SetMinimum(0.01);
  data->SetMaximum(maxY);
  
//   if (minY != -999) {
//     tmpBkg->SetMinimum(minY);
//     sig->SetMinimum(minY);
//     data->SetMinimum(minY);    
//   }


  //CMS Preliminary label
  TPaveText *prelimLabel = new TPaveText(0.21,0.85,0.41,0.90,"NDC");
  prelimLabel->SetTextColor(kBlack);
  prelimLabel->SetFillColor(kWhite);
  prelimLabel->SetBorderSize(0);
  prelimLabel->SetTextAlign(12);
  prelimLabel->SetTextSize(0.03);
  prelimLabel->AddText("CMS Preliminary 2010 #sqrt{s} = 7 TeV");
  prelimLabel->Draw();

//   //Luminosity label
//   TPaveText *tb = new TPaveText(0.21,0.77,0.41,0.82,"NDC");
//   tb->SetTextColor(kBlack);
//   tb->SetFillColor(kWhite);
//   tb->SetBorderSize(0);
//   tb->SetTextAlign(12);
//   string lumi = DoubleToString(fIntegratedLuminosity);
//   tb->AddText((string("#int#font[12]{L}dt = ") + lumi + string(" nb^{ -1}")).c_str());
//   tb->Draw();


  sig->Draw("samehist");
  tmpBkg->Draw("samehist");
  data->Draw("sameE1");

  legend->Draw();
  cv->SaveAs((filename + ".gif").c_str());
  cv->SaveAs((filename + ".eps").c_str());
//   cv->SaveAs((filename + ".C").c_str());
  return cv;

}



//--------------------------------------------------------------------------------------------------
// 
//--------------------------------------------------------------------------------------------------
void PlotUtils::makeDistributionComparisonPlot( vector<string> datasetfiles, 
                                                vector<string> datasetnames, 
                                                string dirname,
                                                vector<string> histNames, 
                                                vector<string> legendNames, bool normalizeArea,
                                                Double_t MinY, Double_t MaxY,
                                                Int_t nbins, string plotname ) {
  assert(histNames.size() == legendNames.size());

  TCanvas *cv = MakeCanvas("cv", plotname.c_str(), 800, 600);
  TLegend *leg = new TLegend(0.25,0.80 - 0.025*histNames.size(),0.45,0.80);   
  leg->SetBorderSize(1);
  leg->SetTextSize(0.03);
  vector <TH1F*> hists;
  Double_t MAXY = 0.0;
  Int_t MaxIndex = -1;
  Double_t Normalization = 0;

  for(UInt_t i=0;i<histNames.size();i++) {

//     TH1F *hist = addAllSamples(datasetfiles, datasetnames, "FakeElectronAnalysisMod", 
//                                histNames[i] );
    TH1F *hist = addAllSamples( datasetfiles, datasetnames, dirname, 
                               histNames[i] );

    hist->SetMarkerColor(fCOLORS[i]);
    hist->SetMarkerSize(1.0);
    leg->AddEntry(hist, legendNames[i].c_str(), "LP"); 
    if (nbins > 0) {
      hist->Rebin(hist->GetXaxis()->GetNbins() / nbins);
    }

    hists.push_back(hist);
    if (hist->GetMaximum() > MAXY) {
      MAXY = hist->GetMaximum();
      MaxIndex = i;
      Normalization = hist->Integral();
    }
    
  }
    
  if (normalizeArea) {
    for(UInt_t i=0;i<hists.size();i++) {
      if (hists[i]->Integral() > 0) hists[i]->Scale(Normalization/hists[i]->Integral());      
    }
  
    MAXY = 0.0;
    for(UInt_t i=0;i<hists.size();i++) {
      if (hists[i]->GetMaximum() > MAXY) {
        MAXY = hists[i]->GetMaximum();
      }
    }
  }

  for(UInt_t i=0;i<hists.size();i++) {
    //do plots here
    if (i==0) {
      if (MaxY == -99) {
        hists[i]->SetMaximum(MAXY*1.1);
      } else {
        hists[i]->SetMaximum(MaxY);
      }
      if (MaxY != -99) {
        hists[i]->SetMinimum(MinY);        
      }

      hists[i]->Draw();      
    } else {
      hists[i]->Draw("same");        
    }


  }
  leg->Draw();

  string filename = plotname + ".gif";
  cv->SaveAs(filename.c_str());
  cv->Delete();

}

//--------------------------------------------------------------------------------------------------
// 
//--------------------------------------------------------------------------------------------------
void PlotUtils::makeXSliceDistributionComparisonPlot( vector<string> datasetfiles, 
                                                      vector<string> datasetnames, 
                                                      string dirname,
                                                      vector<string> histNames, 
                                                      vector<string> legendNames, bool normalizeArea,
                                                      vector<Double_t> xbins, vector<Double_t> ybins,
                                                      string xAxisLabel,
                                                      string yAxisLabel,
                                                      Double_t xlow, Double_t xhigh, 
                                                      Double_t ylow, Double_t yhigh, 
                                                      Double_t legendX1, Double_t legendX2, 
                                                      Double_t legendY1, Double_t legendY2, 
                                                      string plotname) {
  assert(histNames.size() == legendNames.size());
  assert(histNames.size() > 0);

  vector <TH2F*> hists;

  for(UInt_t i=0;i<histNames.size();i++) {    
    TH2F *hist = addAllSamples2D( datasetfiles, datasetnames, dirname, 
                                  histNames[i], xbins, ybins );    
    hists.push_back(hist);
  }
  
  makeXSliceDistributionComparisonPlot(hists,legendNames, normalizeArea,
                                       xAxisLabel,yAxisLabel,
                                       xlow,xhigh,ylow,yhigh,
                                       legendX1,legendX2,legendY1,legendY2,plotname);
  return;
}

//--------------------------------------------------------------------------------------------------
// 
//--------------------------------------------------------------------------------------------------
void PlotUtils::makeXSliceDistributionComparisonPlot( vector<TH2F*> hists,
                                                      vector<string> legendNames,
                                                      bool normalizeArea,
                                                      string xAxisLabel,
                                                      string yAxisLabel,
                                                      Double_t xlow, Double_t xhigh, 
                                                      Double_t ylow, Double_t yhigh, 
                                                      Double_t legendX1, Double_t legendX2, 
                                                      Double_t legendY1, Double_t legendY2, 
                                                      string plotname) {
  assert(hists.size() == legendNames.size());
  assert(hists.size() > 0);
  
  for(UInt_t i=1; int(i) <= hists[0]->GetXaxis()->GetNbins(); ++i) {
    char tmp[20]; 
    sprintf(tmp, " :  %.2f To %.2f",  hists[0]->GetXaxis()->GetBinLowEdge(i), 
            hists[0]->GetXaxis()->GetBinUpEdge(i));
    string sliceLabel = tmp;
    char tmp2[20]; 
    sprintf(tmp2, "%.d", i);
    string sliceLabelForFilename = tmp2;
    string plotLabel = plotname + " " + sliceLabel;
    string plotFilename = plotname + "_" + sliceLabelForFilename;

    vector<TH1F*> slices;
    Double_t Normalization = 0;
    Double_t MaxY = 0;
    for(UInt_t j=0; j < hists.size(); ++j) {
      char t[20]; 
      sprintf(t, "%d", j);
      string indexLabel = t;

      TH1F *slice = (TH1F*)hists[j]->ProjectionY((string("XSlice")+"_"+indexLabel).c_str(),i,i);
      slice->SetTitle(plotLabel.c_str());
      if (xAxisLabel != "") slice->GetXaxis()->SetTitle(xAxisLabel.c_str());
      slice->GetXaxis()->SetTitleOffset(1.0);
      if (yAxisLabel != "") slice->GetYaxis()->SetTitle(yAxisLabel.c_str());
      slice->GetYaxis()->SetTitleOffset(1.5);
      slices.push_back(slice);

      //normalize to first hist
      if (j==0) {
        Normalization = slice->Integral();
        MaxY = slice->GetMaximum();
      }
      if (j>0) {
        if (normalizeArea ) {
          if (slice->Integral() > 0) {
            slice->Scale(Normalization / slice->Integral());
          }
        }
      }
      if (slice->GetMaximum() > MaxY) MaxY = slice->GetMaximum();
    }

    if (yhigh != -99) MaxY = yhigh;
    
    TLegend *leg=0;
    if (legendX1 > -99) {
      leg = new TLegend(legendX1,legendY1,legendX2,legendY2);   
    } else {
      leg = new TLegend(0.25,0.75,0.55,0.9);   
    }
    leg->SetBorderSize(1);
    leg->SetTextSize(0.03);
    TCanvas *cv = MakeCanvas("cv", plotLabel.c_str(), 800, 900);

    for(UInt_t j=0; j < slices.size(); ++j) {
      leg->AddEntry(slices[j], legendNames[j].c_str(), "LP"); 
      slices[j]->SetMaximum(MaxY);
      if (ylow != -99) slices[j]->SetMinimum(ylow);
      slices[j]->SetLineColor(fCOLORS[j]);      
      slices[j]->SetMarkerColor(fCOLORS[j]);
      if (j==0) slices[j]->Draw("E1");
      else slices[j]->Draw("E1same");
    }
    leg->Draw();
    cv->SaveAs((plotFilename+".gif").c_str());
  }
}


//--------------------------------------------------------------------------------------------------
// 
//--------------------------------------------------------------------------------------------------
void PlotUtils::makeYSliceDistributionComparisonPlot( vector<string> datasetfiles, 
                                                      vector<string> datasetnames, 
                                                      string dirname,
                                                      vector<string> histNames, 
                                                      vector<string> legendNames, bool normalizeArea,
                                                      vector<Double_t> xbins, vector<Double_t> ybins,
                                                      string xAxisLabel,
                                                      string yAxisLabel,
                                                      Double_t xlow, Double_t xhigh, 
                                                      Double_t ylow, Double_t yhigh, 
                                                      Double_t legendX1, Double_t legendX2, 
                                                      Double_t legendY1, Double_t legendY2, 
                                                      string plotname) {
  assert(histNames.size() == legendNames.size());
  assert(histNames.size() > 0);

  vector <TH2F*> hists;

  for(UInt_t i=0;i<histNames.size();i++) {    
    TH2F *hist = addAllSamples2D( datasetfiles, datasetnames, dirname, 
                                  histNames[i], xbins, ybins );    
    hists.push_back(hist);
  }
  
  makeYSliceDistributionComparisonPlot(hists,legendNames, normalizeArea,
                                       xAxisLabel,yAxisLabel,
                                       xlow,xhigh,ylow,yhigh,
                                       legendX1,legendX2,legendY1,legendY2,plotname);
  return;
}

//--------------------------------------------------------------------------------------------------
// 
//--------------------------------------------------------------------------------------------------
void PlotUtils::makeYSliceDistributionComparisonPlot( vector<TH2F*> hists,
                                                      vector<string> legendNames,
                                                      bool normalizeArea,
                                                      string xAxisLabel,
                                                      string yAxisLabel,
                                                      Double_t xlow, Double_t xhigh, 
                                                      Double_t ylow, Double_t yhigh, 
                                                      Double_t legendX1, Double_t legendX2, 
                                                      Double_t legendY1, Double_t legendY2, 
                                                      string plotname) {
  assert(hists.size() == legendNames.size());
  assert(hists.size() > 0);
  
  for(UInt_t i=1; int(i) <= hists[0]->GetYaxis()->GetNbins(); ++i) {
    char tmp[20]; 
    sprintf(tmp, " :  %.2f To %.2f",  hists[0]->GetYaxis()->GetBinLowEdge(i), 
            hists[0]->GetYaxis()->GetBinUpEdge(i));
    string sliceLabel = tmp;
    char tmp2[20]; 
    sprintf(tmp2, "%.d", i);
    string sliceLabelForFilename = tmp2;
    string plotLabel = plotname + " " + sliceLabel;
    string plotFilename = plotname + "_" + sliceLabelForFilename;

    vector<TH1F*> slices;
    Double_t Normalization = 0;
    Double_t MaxY = 0;
    for(UInt_t j=0; j < hists.size(); ++j) {
      char t[20]; 
      sprintf(t, "%d", j);
      string indexLabel = t;

      TH1F *slice = (TH1F*)hists[j]->ProjectionX((string("YSlice")+"_"+indexLabel).c_str(),i,i);
      slice->SetTitle(plotLabel.c_str());
      if (xAxisLabel != "") slice->GetXaxis()->SetTitle(xAxisLabel.c_str());
      slice->GetXaxis()->SetTitleOffset(1.0);
      if (yAxisLabel != "") slice->GetYaxis()->SetTitle(yAxisLabel.c_str());
      slice->GetYaxis()->SetTitleOffset(1.5);
      slices.push_back(slice);

      //normalize to first hist
      if (j==0) {
        Normalization = slice->Integral();
        MaxY = slice->GetMaximum();
      }
      if (j>0) {
        if (normalizeArea ) {
          if (slice->Integral() > 0) {
            slice->Scale(Normalization / slice->Integral());
          }
        }
      }
      if (slice->GetMaximum() > MaxY) MaxY = slice->GetMaximum();
    }

    if (yhigh != -99) MaxY = yhigh;
    
    TLegend *leg=0;
    if (legendX1 > -99) {
      leg = new TLegend(legendX1,legendY1,legendX2,legendY2);   
    } else {
      leg = new TLegend(0.25,0.75,0.55,0.9);   
    }
    leg->SetBorderSize(1);
    leg->SetTextSize(0.03);
    TCanvas *cv = MakeCanvas("cv", plotLabel.c_str(), 800, 900);

    for(UInt_t j=0; j < slices.size(); ++j) {
      leg->AddEntry(slices[j], legendNames[j].c_str(), "LP"); 
      slices[j]->SetMaximum(MaxY);
      if (ylow != -99) slices[j]->SetMinimum(ylow);
      slices[j]->SetLineColor(fCOLORS[j]);      
      slices[j]->SetMarkerColor(fCOLORS[j]);
      if (j==0) slices[j]->Draw("E1");
      else slices[j]->Draw("E1same");
    }
    leg->Draw();
    cv->SaveAs((plotFilename+".gif").c_str());
  }
}

//--------------------------------------------------------------------------------------------------
// 
//--------------------------------------------------------------------------------------------------
void PlotUtils::makeCrossDatasetComparisonPlot( vector<string> dataset1files, 
                                                vector<string> dataset1names, string dataset1label,
                                                vector<string> dataset2files, 
                                                vector<string> dataset2names, string dataset2label,
                                                string dirname1 , vector<string> histNames1, vector<string> legendNames1, 
                                                string dirname2, vector<string> histNames2, vector<string> legendNames2, 
                                                string plotname ,
                                                bool normalizeArea, 
                                                string xAxisLabel,
                                                string yAxisLabel,
                                                Double_t xlow, Double_t xhigh, 
                                                Double_t ylow, Double_t yhigh, 
                                                Double_t legendX1, Double_t legendX2, 
                                                Double_t legendY1, Double_t legendY2, 
                                                Bool_t useLogY,
                                                Int_t nbins) {
  
  assert(histNames1.size() == legendNames1.size());
  assert(histNames2.size() == legendNames2.size());

  TCanvas *cv = MakeCanvas("cv", plotname.c_str(), 800, 600);

  TLegend *leg = 0;
  if(legendX1 == -99) {
    leg = new TLegend(0.65,0.6,0.9,0.8);   
  } else {
    leg = new TLegend(legendX1,legendY1,legendX2,legendY2);     
  }
  leg->SetBorderSize(1);
  leg->SetTextSize(0.03);
  
  int colorindex = 0;

  if (useLogY) cv->SetLogy();

  vector <TH1F*> hists1;
  vector <TH1F*> hists2;
  Double_t Normalization = 0;

  //determine Best Y-axis scale
  Double_t MAXY = 0;
  for(UInt_t i=0;i<histNames1.size();i++) {
    TH1F *hist1 = addAllSamples(dataset1files, dataset1names, dirname1, histNames1[i] );
    if (nbins > 0) {
      vector<Double_t> bins;
      for (int b=0; b<nbins; ++b) {
        Double_t binsize = (hist1->GetXaxis()->GetBinUpEdge(hist1->GetXaxis()->GetNbins()) - 
                            hist1->GetXaxis()->GetBinLowEdge(1)) / nbins;
        bins.push_back(hist1->GetXaxis()->GetBinLowEdge(1) + binsize * b);
      }
      hist1 = rebin(hist1, bins);
    }
    cerr << "hist1 Max: " << hist1->GetMaximum() << " MAXY=" << MAXY << endl;
    if (hist1->GetMaximum() > MAXY) {
      MAXY = hist1->GetMaximum();
      Normalization = hist1->Integral();
    }
    hists1.push_back(hist1);
  }

  for(UInt_t i=0;i<histNames2.size();i++) {
    TH1F *hist2 = addAllSamples(dataset2files, dataset2names, dirname2, histNames2[i] );    
    if (nbins > 0) {
      vector<Double_t> bins;
      for (int b=0; b<nbins; ++b) {
        Double_t binsize = (hist2->GetXaxis()->GetBinUpEdge(hist2->GetXaxis()->GetNbins()) - 
                            hist2->GetXaxis()->GetBinLowEdge(1)) / nbins;
        bins.push_back(hist2->GetXaxis()->GetBinLowEdge(1) + binsize * b);
      }
      hist2 = rebin(hist2, bins);
    }
    cerr << "hist2 Max: " << hist2->GetMaximum() << " MAXY=" << MAXY << endl;
    if (hist2->GetMaximum() > MAXY) {
      MAXY = hist2->GetMaximum();
      Normalization = hist2->Integral();
    }
    hists2.push_back(hist2);
  }
  
  //Normalize to histogram 1
  Normalization = hists1[0]->Integral();

  if (normalizeArea) {
    for(UInt_t i=0;i<hists1.size();i++) {
      if (hists1[i]->Integral() > 0) hists1[i]->Scale(Normalization/hists1[i]->Integral());      
    }
    for(UInt_t i=0;i<hists2.size();i++) {
      if (hists2[i]->Integral() > 0) hists2[i]->Scale(Normalization/hists2[i]->Integral());      
    }

    //recalculate max scale after normalization
    MAXY = 0;
    for(UInt_t i=0;i<hists1.size();i++) {
      if (hists1[i]->GetMaximum() > MAXY) {
        MAXY = hists1[i]->GetMaximum();
      }
    }
    for(UInt_t i=0;i<hists2.size();i++) {
      if (hists2[i]->GetMaximum() > MAXY) {
        MAXY = hists2[i]->GetMaximum();
      }
    }
  }  

  for(UInt_t i=0;i<hists1.size();i++) {
    hists1[i]->SetMarkerColor(fCOLORS[colorindex]);
    hists1[i]->SetMarkerStyle(fMARKERS[colorindex]);
    hists1[i]->SetMarkerSize(1.0);
    leg->AddEntry(hists1[i], (dataset1label+legendNames1[i]).c_str(), "LP"); 

    if (xlow != -99 && xhigh != -99) {
        hists1[i]->GetXaxis()->SetRangeUser(xlow,xhigh);
    } 
    if (ylow != -99) {
      hists1[i]->SetMinimum(ylow);
    }
    if (yhigh != -99) {
      hists1[i]->SetMaximum(yhigh);
    } else {
      hists1[i]->SetMaximum(MAXY*1.2);
    }

    if (xAxisLabel != "") {
      hists1[i]->SetXTitle(xAxisLabel.c_str());
    }
    if (yAxisLabel != "") {
      hists1[i]->SetYTitle(yAxisLabel.c_str());
    }

    //do plots here
    if (i==0) {
      hists1[i]->DrawCopy();
    } else {
      hists1[i]->DrawCopy("same");  
    }
    ++colorindex;
  }  
  for(UInt_t i=0;i<hists2.size();i++) {
    
    hists2[i]->SetMarkerColor(fCOLORS[colorindex]);
    hists2[i]->SetMarkerStyle(fMARKERS[colorindex]);
    hists2[i]->SetMarkerSize(1.0);
    leg->AddEntry(hists2[i], (dataset2label+legendNames2[i]).c_str(), "LP"); 

    if (xlow != -99 && xhigh != -99) {
        hists2[i]->GetXaxis()->SetRangeUser(xlow,xhigh);
    } 
    if (ylow != -99) {
      hists2[i]->SetMinimum(ylow);
    }
    if (yhigh != -99) {
      hists2[i]->SetMaximum(yhigh);
    } else {
      hists2[i]->SetMaximum(MAXY*1.2);
    }

    if (xAxisLabel != "") {
      hists2[i]->SetXTitle(xAxisLabel.c_str());      
    }
    if (yAxisLabel != "") {
      hists2[i]->SetYTitle(yAxisLabel.c_str());
    }

    hists2[i]->DrawCopy("same");
    ++colorindex;
  }

  leg->Draw();

  string filename = plotname + ".gif";
  cv->SaveAs(filename.c_str());
  cv->Delete();

}

//--------------------------------------------------------------------------------------------------
void PlotUtils::makeCrossDirComparisonPlot( vector<string> datasetfiles, 
                                            vector<string> datasetnames, 
                                            vector<string> dirnames, vector<string> dirnamelabel, 
                                            vector<string> histNames,vector<string> legendNames, 
                                            vector<Double_t> bins, string plotname , 
                                            string xAxisLabel,
                                            string yAxisLabel,
                                            Double_t xlow, Double_t xhigh, 
                                            Double_t ylow, Double_t yhigh, 
                                            Double_t legendX1, Double_t legendX2, 
                                            Double_t legendY1, Double_t legendY2, 
                                            Bool_t useLogY) {

  assert(histNames.size() > 0);
  assert(histNames.size() == legendNames.size());
  assert(datasetfiles.size() == datasetnames.size());
  assert(dirnames.size() == dirnamelabel.size());
   
  vector <TH1F*> hists;
  Double_t MAXY = 0.0;

  TCanvas *cv = MakeCanvas("cv", plotname.c_str(), 800, 600);

  TLegend *leg = 0;
  if(legendX1 == -99) {
    leg = new TLegend(0.55,0.75,0.9,0.9);     
  } else {
    leg = new TLegend(legendX1,legendY1,legendX2,legendY2);     
  }

  leg->SetBorderSize(1);
  leg->SetTextSize(0.03);

  int colorindex = 0;
  for(UInt_t i=0;i<histNames.size();++i) {
    for(UInt_t d=0;d<dirnames.size();++d) {
      
      TH1F *tmphist = addAllSamples(datasetfiles, datasetnames, dirnames[d], 
                                    histNames[i]);
      
      TH1F *hist = tmphist;
      if (bins.size() > 0) hist = rebin(tmphist,bins);
      hist->SetTitle(plotname.c_str());
      hist->SetMarkerColor(fCOLORS[colorindex]);
      hist->SetLineColor(fCOLORS[colorindex]);
      hist->SetMarkerStyle(fMARKERS[colorindex]);
      hist->SetMarkerSize(1.0);
      if (xlow != -99 && xhigh != -99) {
        hist->GetXaxis()->SetRangeUser(xlow,xhigh);
      } 
      if (xAxisLabel != "") hist->SetXTitle(xAxisLabel.c_str());
      if (yAxisLabel != "") hist->SetYTitle(yAxisLabel.c_str());


      leg->AddEntry(hist, (legendNames[i]+dirnamelabel[d]).c_str(), "LP");       
      hists.push_back(hist);
      if (hist->GetMaximum() > MAXY) MAXY = hist->GetMaximum();
      ++colorindex;

    }
  }

  if (useLogY) cv->SetLogy();

  for (UInt_t i=0; i < hists.size();++i) {

    cout << i << " " << hists[i] << endl;
    hists[i]->SetMaximum(1.2*MAXY);
    hists[i]->SetMinimum(0.0);

    if (ylow != -99) hists[i]->SetMinimum(ylow);
    if (yhigh != -99) hists[i]->SetMaximum(yhigh);

    if (i==0) 
      hists[i]->Draw("");
    else
      hists[i]->Draw("same");
  } 
  leg->Draw();
  
  cv->SaveAs((plotname+".gif").c_str());
  cv->SaveAs((plotname+".eps").c_str());
  cv->Delete();
}


//--------------------------------------------------------------------------------------------------
// Compare histograms from different directories with systematic errors added to the histogram
// in the first given directory.
//--------------------------------------------------------------------------------------------------
void PlotUtils::makeComparisonPlotWithSystematics( vector<string> datasetfiles, 
                                                   vector<string> datasetnames, 
                                                   vector<string> dirnames, 
                                                   vector<string> dirnamelabel, 
                                                   vector<string> histNames, 
                                                   vector<string> errorhistNames, 
                                                   vector<string> legendNames, 
                                                   vector<Double_t> bins, string plotname , 
                                                   string xAxisLabel,
                                                   string yAxisLabel,
                                                   Double_t xlow, Double_t xhigh, 
                                                   Double_t ylow, Double_t yhigh, 
                                                   Double_t legendX1, Double_t legendX2, 
                                                   Double_t legendY1, Double_t legendY2, 
                                                   Bool_t useLogY) {
  
  assert(histNames.size() == 1);
  assert(histNames.size() == legendNames.size());
  assert(datasetfiles.size() == datasetnames.size());
  assert(dirnames.size() == dirnamelabel.size());
    
  assert(errorhistNames.size() == histNames.size());
 
  vector <TH1F*> hists;
  Double_t MAXY = 0.0;

  TCanvas *cv = MakeCanvas("cv", plotname.c_str(), 800, 600);

  TLegend *leg = 0;
  if(legendX1 == -99) {
    leg = new TLegend(0.55,0.75,0.9,0.9);     
  } else {
    leg = new TLegend(legendX1,legendY1,legendX2,legendY2);     
  }

  leg->SetBorderSize(1);
  leg->SetTextSize(0.03);

  int colorindex = 0;
  for(UInt_t i=0;i<histNames.size();++i) {
    for(UInt_t d=0;d<dirnames.size();++d) {
      
      TH1F *tmphist = addAllSamples(datasetfiles, datasetnames, dirnames[d], 
                                    histNames[i]);
      TH1F *tmpErrorHist = addAllSamples(datasetfiles, datasetnames, dirnames[d], 
                                    errorhistNames[i]);
      
      if (tmphist->GetXaxis()->GetNbins() != tmpErrorHist->GetXaxis()->GetNbins()) {
        cerr << "Error: histNames[" << i << "] and errorhistNames[" << i << "] do not have the same bins.\n";
        assert(tmphist->GetXaxis()->GetNbins() == tmpErrorHist->GetXaxis()->GetNbins());
      }

      TH1F *hist = tmphist;
      TH1F *errorHist = tmpErrorHist;
      if (bins.size() > 0) { 
        hist = rebin(tmphist,bins);
        errorHist = rebin(tmpErrorHist,bins);
      }

      hist->SetTitle(plotname.c_str());
      hist->SetMarkerColor(fSYSCOLORS[colorindex]);
      hist->SetLineColor(fSYSCOLORS[colorindex]);
      hist->SetFillColor(fSYSCOLORS[colorindex]);
      hist->SetFillStyle(3544);
      hist->SetMarkerStyle(fMARKERS[colorindex]);
      hist->SetMarkerSize(1.0);
      if (xlow != -99 && xhigh != -99) {
        hist->GetXaxis()->SetRangeUser(xlow,xhigh);
      } 
      if (xAxisLabel != "") hist->SetXTitle(xAxisLabel.c_str());
      if (yAxisLabel != "") hist->SetYTitle(yAxisLabel.c_str());

      leg->AddEntry(hist, (legendNames[i]+dirnamelabel[d]).c_str(), "LP");       
      hists.push_back(hist);
      
      if (hist->GetMaximum() > MAXY) MAXY = hist->GetMaximum();

      //for the first directory, add systematic errors to the statistical error
      if (d==0) {
        for (int bin = 0; bin < hist->GetNbinsX()+1; ++bin) {
          Double_t totalError = TMath::Sqrt(hist->GetBinError(bin)*hist->GetBinError(bin) +
                                            errorHist->GetBinContent(bin)*errorHist->GetBinContent(bin));
          cout << "bin " << bin << " : " << totalError << ", " << hist->GetBinError(bin) << ", " << errorHist->GetBinContent(bin) << endl;
          if (TMath::IsNaN(errorHist->GetBinContent(bin))) cout << "NAN!!\n";
          hist->SetBinError(bin, totalError);
        }
      }      
      ++colorindex;
    }
  }

  if (useLogY) cv->SetLogy();

  for (UInt_t i=0; i < hists.size();++i) {

    cout << i << " " << hists[i] << endl;
    hists[i]->SetMaximum(1.2*MAXY);
    hists[i]->SetMinimum(0.0);

    if (ylow != -99) hists[i]->SetMinimum(ylow);
    if (yhigh != -99) hists[i]->SetMaximum(yhigh);

    if (i==0) 
      hists[i]->Draw("e5");
    else
      hists[i]->Draw("same");
  } 
  leg->Draw();
  
  cv->SaveAs((plotname+".gif").c_str());
  cv->SaveAs((plotname+".eps").c_str());
  cv->Delete();
}


//--------------------------------------------------------------------------------------------------
// Create Efficiency Histogram. Use binomial Errors
//--------------------------------------------------------------------------------------------------
TH1F* PlotUtils::createEfficiencyHist(vector<string> datasetFiles, vector<string> datasetNames, string dirName,
                                      string numeratorHistname, string denominatorHistname, string histname, 
                                      vector<Double_t> bins, Double_t xlow, Double_t xhigh, 
                                      Double_t ylow, Double_t yhigh 
  ) {

  //find the largest weight among all the samples
  Double_t maxWeight = 0;
  for (UInt_t s=0; s < datasetFiles.size() ; ++s ) {      
    Double_t CrossSection = xstab->Get(datasetNames[s].c_str());
    TH1F *tmpNEventsHist = getHisto(datasetFiles[s],"AnaFwkMod", "hDEvents");
    Double_t NEvents = tmpNEventsHist->Integral();
    Double_t weight = CrossSection * fIntegratedLuminosity / NEvents;
    if (weight > maxWeight || maxWeight == 0) maxWeight = weight;     
  }
  
  TH1F *denominator = addAllSamples(datasetFiles, datasetNames, dirName, denominatorHistname);
  TH1F *numerator = addAllSamples(datasetFiles, datasetNames, dirName, numeratorHistname);
  
  TH1F *n = rebin(numerator,bins);
  TH1F *d = rebin(denominator,bins);

  //create a fake numerator histogram where any bin that has 0 is filled by one event * maxWeight
  TH1F *n_AddedMaxWeightToZeroBins = (TH1F*)n->Clone(histname.c_str());
  for (int b=0;b<=n_AddedMaxWeightToZeroBins->GetNbinsX()+1;b++) {
    if (n_AddedMaxWeightToZeroBins->GetBinContent(b) == 0) {
      n_AddedMaxWeightToZeroBins->SetBinContent(b,maxWeight);
      n_AddedMaxWeightToZeroBins->SetBinError(b,maxWeight);
    }
  }

  TH1F *efficiency = (TH1F*)n->Clone(histname.c_str());
  efficiency->GetYaxis()->SetTitle("Efficiency");

  efficiency->Divide(n, d, 1.0,1.0,"B");  
  
  TH1F *efficiency_AddedMaxWeightToZeroBins = (TH1F*)n->Clone(histname.c_str());
  efficiency_AddedMaxWeightToZeroBins->GetYaxis()->SetTitle("Efficiency");
  efficiency_AddedMaxWeightToZeroBins->Divide(n_AddedMaxWeightToZeroBins, d, 1.0,1.0,"B");  
  
  //replace the errors in the bins with 0 by errors from the fake rate where I incremented
  //the bins with 0 numerator to maxWeight.
  for (int b=0;b<=efficiency->GetNbinsX()+1;b++) {
    if (efficiency->GetBinContent(b) == 0) {
      efficiency->SetBinError(b,efficiency_AddedMaxWeightToZeroBins->GetBinError(b));
    }
  } 
  
  if (yhigh != -99)
    efficiency->SetMaximum(yhigh);
  if (ylow != -99)
    efficiency->SetMinimum(ylow);
  if (xlow != -99 && xhigh != -99) 
    efficiency->GetXaxis()->SetRangeUser(xlow,xhigh);

  return efficiency;
}

//--------------------------------------------------------------------------------------------------
// Create Efficiency Graph. Use RooStatsCms Errors. 
// errorType == 0 : ROOT BayesDivide
// errorType == 1 : Feldman Cousins Confidence Intervals
// errorType == 2 : Clopper Pearson Confidence Intervals
//--------------------------------------------------------------------------------------------------
TGraphAsymmErrors* PlotUtils::createEfficiencyGraph(vector<string> datasetFiles, 
                                                    vector<string> datasetNames, string dirName,
                                                    string numeratorHistname, 
                                                    string denominatorHistname, string histname, 
                                                    vector<Double_t> bins, Int_t errorType, 
                                                    Double_t xlow, Double_t xhigh, 
                                                    Double_t ylow, Double_t yhigh 
  ) {

  TH1F *denominator = addAllSamples(datasetFiles, datasetNames, dirName, denominatorHistname);
  TH1F *numerator = addAllSamples(datasetFiles, datasetNames, dirName, numeratorHistname);
  
  return createEfficiencyGraph(numerator, denominator, histname, bins, errorType,
                               xlow , xhigh, ylow, yhigh );

}


//--------------------------------------------------------------------------------------------------
// Create Efficiency Graph. Use RooStatsCms Errors. 
// errorType == 0 : ROOT BayesDivide
// errorType == 1 : Feldman Cousins Confidence Intervals
// errorType == 2 : Clopper Pearson Confidence Intervals
//--------------------------------------------------------------------------------------------------
TGraphAsymmErrors* PlotUtils::createEfficiencyGraph(TH1F* numerator, TH1F* denominator,
                                                    string histname, 
                                                    vector<Double_t> bins, Int_t errorType, 
                                                    Double_t xlow, Double_t xhigh, 
                                                    Double_t ylow, Double_t yhigh 
  ) {

  TH1F *n = rebin(numerator,bins);
  TH1F *d = rebin(denominator,bins);
  
  Int_t nbins = n->GetNbinsX()+2;

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

  for (int b=0; b<nbins ; ++b) {

    x[b] = n->GetXaxis()->GetBinCenter(b);    
    xErr[b] = 0.0;

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;    

    Double_t n1 = TMath::Nint(n->GetBinContent(b));
    Double_t n2 = TMath::Nint(d->GetBinContent(b));
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, errorType);
     cerr << " done bin " << b << " " << n1 << "(" << n->GetBinContent(b) << ")" << " / " << n2 << "(" << d->GetBinContent(b) << ")" << " = " << ratio << " " << errLow << " " << errHigh << endl;
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


//*************************************************************************************************
//Normalize Hist
//*************************************************************************************************
void PlotUtils::NormalizeHist(TH1F *hist) {
  Double_t norm = 0;
  hist->SetTitle("");
  for (UInt_t b=0; int(b) < hist->GetXaxis()->GetNbins()+2; ++b) {
    norm += hist->GetBinContent(b);
  }
  for (UInt_t b=0; int(b) < hist->GetXaxis()->GetNbins()+2; ++b) {
    hist->SetBinContent(b,hist->GetBinContent(b) / norm);
    hist->SetBinError(b,hist->GetBinError(b) / norm);
  }

  return;
}


//--------------------------------------------------------------------------------------------------
// Create ROC curves
//--------------------------------------------------------------------------------------------------
TGraphAsymmErrors* PlotUtils::MakeSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, string name, Bool_t cutBelow ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  Double_t *SigEff = new Double_t[nPoints];
  Double_t *BkgEff = new Double_t[nPoints];
  Double_t *SigEffErrLow = new Double_t[nPoints];
  Double_t *SigEffErrHigh = new Double_t[nPoints];
  Double_t *BkgEffErrLow = new Double_t[nPoints];
  Double_t *BkgEffErrHigh = new Double_t[nPoints];
  Double_t NSigTotal = 0;
  Double_t NBkgTotal = 0;
  
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    Double_t nbkg = 0;

    if (cutBelow) {
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
        nbkg += bkgHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
        nbkg += bkgHist->GetBinContent(q);
      }
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;
//     SigEffErrLow[b] = errLow;
//     SigEffErrHigh[b] = errHigh;

    n1 = TMath::Nint(nbkg);
    n2 = TMath::Nint(NBkgTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    BkgEff[b] = ratio;
    BkgEffErrLow[b] = 0;
    BkgEffErrHigh[b] = 0;
//     BkgEffErrLow[b] = errLow;
//     BkgEffErrHigh[b] = errHigh;
  }

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (nPoints, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh, SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsBkgEff->SetMarkerSize(0.5);
  tmpSigEffVsBkgEff->SetMarkerStyle(20);

  delete [] SigEff;
  delete [] BkgEff;
  delete [] SigEffErrLow;
  delete [] SigEffErrHigh;
  delete [] BkgEffErrLow;
  delete [] BkgEffErrHigh;

  return tmpSigEffVsBkgEff;
}

//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* PlotUtils::MakeCurrentWPSigEffVsBkgEffGraph(Double_t signalEff, Double_t bkgEff, string name ) {
  //Make Met Plots
  Double_t SigEff[1];
  Double_t BkgEff[1];
  Double_t SigEffErrLow[1];
  Double_t SigEffErrHigh[1];
  Double_t BkgEffErrLow[1];
  Double_t BkgEffErrHigh[1];

  SigEff[0] = signalEff;
  SigEffErrLow[0] = 0;
  SigEffErrHigh[0] = 0;
  BkgEff[0] = bkgEff;
  BkgEffErrLow[0] = 0;
  BkgEffErrHigh[0] = 0;

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (1, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh , SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsBkgEff->SetMarkerColor(kBlack);
  tmpSigEffVsBkgEff->SetLineColor(kBlack);
  tmpSigEffVsBkgEff->SetMarkerSize(1.5);

  return tmpSigEffVsBkgEff;
}

//--------------------------------------------------------------------------------------------------
// Create ROC curves
//--------------------------------------------------------------------------------------------------
//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* PlotUtils::MakeCurrentWPSigEffVsBkgEffGraph(TH1F* signalHist, TH1F* bkgHist, 
                                                               string name, Double_t myCutValue, 
                                                               Bool_t cutBelow ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  Double_t SigEff[1] = {0};
  Double_t SigEffErrLow[1] = {0};
  Double_t SigEffErrHigh[1] = {0};
  Double_t BkgEff[1] = {0};
  Double_t BkgEffErrLow[1] = {0};
  Double_t BkgEffErrHigh[1] = {0};
  Double_t NSigTotal = 0;
  Double_t NBkgTotal = 0;
  Double_t cutValue = 0;

  Double_t effDiff = 9999;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NBkgTotal += bkgHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    Double_t nbkg = 0;

    if (cutBelow) {
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nbkg += bkgHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
      for (UInt_t q=0; q <= b; ++q) {
        nbkg += bkgHist->GetBinContent(q);
      }
    }

    Double_t sigEff;
    Double_t sigEffErrLow;
    Double_t sigEffErrHigh;     
//     mithep::MathUtils::CalcRatio(TMath::Nint(nsig) , TMath::Nint(NSigTotal), ratio, errLow, errHigh, 2);
    sigEff = nsig / NSigTotal;

    Double_t bkgEff;
    Double_t bkgEffErrLow;
    Double_t bkgEffErrHigh;     
//     mithep::MathUtils::CalcRatio(TMath::Nint(nbkg) , TMath::Nint(NBkgTotal), ratio, errLow, errHigh, 2);
    bkgEff = nbkg / NBkgTotal;
    
//       std::cout << myCutValue << " : " << signalHist->GetXaxis()->GetBinCenter(b) << " , " << cutValue[0] << endl;
    if (fabs(myCutValue - signalHist->GetXaxis()->GetBinCenter(b)) < fabs(myCutValue - cutValue)) {
      cutValue = signalHist->GetXaxis()->GetBinCenter(b);
      SigEff[0] = sigEff;
      SigEffErrLow[0] = 0;
      SigEffErrHigh[0] = 0;
//     SigEffErrLow[0] = sigEffErrLow;
//     SigEffErrHigh[0] = sigEffErrHigh;
      BkgEff[0] = bkgEff;
      BkgEffErrLow[0] = 0;
      BkgEffErrHigh[0] = 0;
//     BkgEffErrLow[0] = bkgEffErrLow;
//     BkgEffErrHigh[0] = bkgEffErrHigh;
    }
  }

//   std::cout << "Final: " << cutValue[0] << " , " << SigEff[0] << endl;

  TGraphAsymmErrors *tmpSigEffVsBkgEff = new TGraphAsymmErrors (1, BkgEff, SigEff, BkgEffErrLow, BkgEffErrHigh , SigEffErrLow, SigEffErrHigh );
  tmpSigEffVsBkgEff->SetName(name.c_str());
  tmpSigEffVsBkgEff->SetTitle("");
  tmpSigEffVsBkgEff->GetXaxis()->SetTitle("Bkg Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitle("Signal Eff");
  tmpSigEffVsBkgEff->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsBkgEff->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsBkgEff->SetMarkerColor(kBlack);
  tmpSigEffVsBkgEff->SetLineColor(kBlack);
  tmpSigEffVsBkgEff->SetMarkerSize(1.5);

  return tmpSigEffVsBkgEff;
}



//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* PlotUtils::MakeSigEffVsCutValueGraph(TH1F* signalHist, string name, 
                                                        Bool_t cutBelow  ) {

  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  Double_t *cutValue = new Double_t[nPoints];
  Double_t *cutValueErr = new Double_t[nPoints];
  Double_t *SigEff = new Double_t[nPoints];
  Double_t *SigEffErrLow = new Double_t[nPoints];
  Double_t *SigEffErrHigh = new Double_t[nPoints];
  Double_t NSigTotal = 0;
  
  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    cutValue[b] = signalHist->GetXaxis()->GetBinCenter(b);
    cutValueErr[b] = 0;
    Double_t nsig = 0;

    if (cutBelow) {
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    SigEff[b] = ratio;
    SigEffErrLow[b] = 0;
    SigEffErrHigh[b] = 0;
//     SigEffErrLow[b] = errLow;
//     SigEffErrHigh[b] = errHigh;

  }

  TGraphAsymmErrors *tmpSigEffVsCut = new TGraphAsymmErrors (nPoints, cutValue, SigEff, cutValueErr, cutValueErr, SigEffErrLow, SigEffErrHigh  );
  tmpSigEffVsCut->SetName(name.c_str());
  tmpSigEffVsCut->SetTitle("");
  tmpSigEffVsCut->GetXaxis()->SetTitle("Cut Value");
  tmpSigEffVsCut->GetYaxis()->SetTitle("Efficiency");
  tmpSigEffVsCut->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsCut->GetXaxis()->SetTitleOffset(1.05);

  delete [] cutValue;
  delete [] cutValueErr;
  delete [] SigEff;
  delete [] SigEffErrLow;
  delete [] SigEffErrHigh;

  return tmpSigEffVsCut;
}


//*************************************************************************************************
//
//*************************************************************************************************
TGraphAsymmErrors* PlotUtils::MakeCurrentWPSigEffVsCutValueGraph(TH1F* signalHist, string name, 
                                                                 Double_t myCutValue, 
                                                                 Bool_t cutBelow  ) {
  //Make Met Plots
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  Double_t cutValue[1] = {0};
  Double_t cutValueErr[1] = {0};
  Double_t SigEff[1] = {0};
  Double_t SigEffErrLow[1] = {0};
  Double_t SigEffErrHigh[1] = {0};
  Double_t NSigTotal = 0;
  
  Double_t effDiff = 9999;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }

  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;

    if (cutBelow) {
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
    }

    Double_t ratio;
    Double_t errLow;
    Double_t errHigh;     
    Double_t n1 = 0;
    Double_t n2 = 0;

    n1 = TMath::Nint(nsig);
    n2 = TMath::Nint(NSigTotal);
    mithep::MathUtils::CalcRatio(n1 , n2, ratio, errLow, errHigh, 2);
    
//       std::cout << myCutValue << " : " << signalHist->GetXaxis()->GetBinCenter(b) << " , " << cutValue[0] << endl;
    if (fabs(myCutValue - signalHist->GetXaxis()->GetBinCenter(b)) < fabs(myCutValue - cutValue[0])) {
      SigEff[0] = ratio;
      SigEffErrLow[0] = 0;
      SigEffErrHigh[0] = 0;
//     SigEffErrLow[0] = errLow;
//     SigEffErrHigh[0] = errHigh;
      cutValue[0] = signalHist->GetXaxis()->GetBinCenter(b);
      cutValueErr[0] = 0;
    }
  }

//   std::cout << "Final: " << cutValue[0] << " , " << SigEff[0] << endl;

  TGraphAsymmErrors *tmpSigEffVsCut = new TGraphAsymmErrors (1, cutValue, SigEff, cutValueErr, cutValueErr, SigEffErrLow, SigEffErrHigh  );
  tmpSigEffVsCut->SetName(name.c_str());
  tmpSigEffVsCut->SetTitle("");
  tmpSigEffVsCut->GetXaxis()->SetTitle("Cut Value");
  tmpSigEffVsCut->GetYaxis()->SetTitle("Efficiency");
  tmpSigEffVsCut->GetYaxis()->SetTitleOffset(1.1);
  tmpSigEffVsCut->GetXaxis()->SetTitleOffset(1.05);
  tmpSigEffVsCut->SetMarkerColor(kBlack);
  tmpSigEffVsCut->SetLineColor(kBlack);
  tmpSigEffVsCut->SetMarkerSize(1.5);

  return tmpSigEffVsCut;
}


//*************************************************************************************************
//
//*************************************************************************************************
Double_t PlotUtils::FindCutValueAtFixedSignalEfficiency(TH1F* signalHist, Double_t targetSignalEff, 
                                                        Bool_t cutBelow  ) {
  //Make Met Plots


  Double_t targetCutValue = -9999;
  Double_t bestCurrentSignalEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  Double_t NSigTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;

    if (cutBelow) {
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
    }

    Double_t ratio = nsig / NSigTotal;
//     std::cout << targetSignalEff << " : " << ratio << " , " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetSignalEff - ratio) < fabs(targetSignalEff - bestCurrentSignalEff)) {
      targetCutValue = signalHist->GetXaxis()->GetBinCenter(b);
      bestCurrentSignalEff = ratio;
    }
  }

  return targetCutValue;
}


//*************************************************************************************************
//
//*************************************************************************************************
Double_t PlotUtils::FindBkgEffAtFixedSignalEfficiency(TH1F* signalHist, TH1F* bkgHist, 
                                                      Double_t targetSignalEff, 
                                                      Bool_t cutBelow ) {
  //Make Met Plots


  Double_t targetBkgEff = 0;
  Double_t bestCurrentSignalEff = 0;
  const UInt_t nPoints = signalHist->GetXaxis()->GetNbins();
  Double_t NSigTotal = 0;
  Double_t NBkgTotal = 0;

  for (UInt_t q=0; q < nPoints+2; ++q) {
    NSigTotal += signalHist->GetBinContent(q);
    NBkgTotal += bkgHist->GetBinContent(q);
  }


  for(UInt_t b=0; b < nPoints; ++b) {
    Double_t nsig = 0;
    Double_t nbkg = 0;

    if (cutBelow) {
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
      for (UInt_t q=b; q < nPoints+2; ++q) {
        nbkg += bkgHist->GetBinContent(q);
      }
    } else {
      for (UInt_t q=0; q <= b; ++q) {
        nsig += signalHist->GetBinContent(q);
      }
      for (UInt_t q=0; q <= b; ++q) {
        nbkg += bkgHist->GetBinContent(q);
      }
    }

    Double_t ratio = nsig / NSigTotal;
    Double_t bkgEff = nbkg / NBkgTotal;
//     std::cout << targetSignalEff << " : " << ratio << " , " << bkgEff << " : " << signalHist->GetXaxis()->GetBinCenter(b) << endl;

    if (fabs(targetSignalEff - ratio) < fabs(targetSignalEff - bestCurrentSignalEff)) {
      bestCurrentSignalEff = ratio;
      targetBkgEff = bkgEff;
    }
  }

  return targetBkgEff;
}



