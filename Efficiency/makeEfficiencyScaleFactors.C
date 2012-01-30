
#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TGraphAsymmErrors.h>      // graphs
#include <TH2F.h>                   // 2D histograms
#include <TMath.h>                  // ROOT math library
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#endif

//=== MAIN MACRO =================================================================================================

void makeEfficiencyScaleFactors(string DataFilename = "data_eleWPEffTP/basic_76_106/eff.root", 
                   string MCFilename = "s11_zee_eleWPEffTP/basic_76_106/eff.root" , 
                   string outputDir = "data_eleWPEffTP/basic_76_106/", 
                   string histName = "h2_results_ele_selection",
                   string Label = "htt")
{
  gBenchmark->Start("printMuonWPEff");
    
  string label = Label;
  if (Label != "") label = "_" + Label;

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  TString outfname = (outputDir + "/eff_table.txt").c_str();
  TFile mcfile(MCFilename.c_str());
  TFile datafile(DataFilename.c_str());
  
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================   
  
  TH2F *hMCEff=0,    *hMCErrl=0,    *hMCErrh=0;
  TH2F *hDataEff1=0, *hDataErrl1=0, *hDataErrh1=0;
  TH2F *hDataEff2=0, *hDataErrl2=0, *hDataErrh2=0;
  
  hMCEff  = (TH2F*)mcfile.Get("hEffEtaPt");
  hMCErrl = (TH2F*)mcfile.Get("hErrlEtaPt");
  hMCErrh = (TH2F*)mcfile.Get("hErrhEtaPt");
  
  hDataEff  = (TH2F*)datafile.Get("hEffEtaPt");
  hDataErrl = (TH2F*)datafile.Get("hErrlEtaPt");
  hDataErrh = (TH2F*)datafile.Get("hErrhEtaPt");
  
  TGraphAsymmErrors* MCEffVsRho   = ( TGraphAsymmErrors*)mcfile.Get("grEffRho");
  TGraphAsymmErrors* DataEffVsRho = ( TGraphAsymmErrors*)datafile.Get("grEffRho");
 
  TGraphAsymmErrors* MCEffVsNVtx   = ( TGraphAsymmErrors*)mcfile.Get("grEffNPV");
  TGraphAsymmErrors* DataEffVsNVtx = ( TGraphAsymmErrors*)datafile.Get("grEffNPV");


  //--------------------------------------------------------------------------------------------------------------
  // Update root file histograms
  //==============================================================================================================   
  const Int_t nx = hMCEff->GetNbinsX();
  const Int_t ny = hMCEff->GetNbinsY();

  TFile *outputFile = new TFile(("efficiency_results"+label+".root").c_str(), "UPDATE");

  //Do Binning
  Double_t *ptbins = new Double_t[ny+1];
  Double_t *etabins = new Double_t[nx+1];
  for(Int_t iy=1; iy<=ny; iy++) {
      ptbins[iy-1] = hMCEff->GetYaxis()->GetBinLowEdge(iy);
  }
  for(Int_t ix=1; ix<=nx; ix++) {
    etabins[ix-1]= hMCEff->GetXaxis()->GetBinLowEdge(ix);
  }
  ptbins[ny] = 60;
  etabins[nx] = 2.5;
    
  TH2F *h2_results_selection = new TH2F(histName.c_str(),"",ny,ptbins,nx,etabins);
  for(Int_t iy=0; iy<=ny+2; iy++) {
    for(Int_t ix=0; ix<=nx+2; ix++) {
      h2_results_selection->SetCellContent(iy,ix, 1.0);
      h2_results_selection->SetCellError(iy,ix, 0.0);
    }
  }

  cout << "ptbins : ";
  for(Int_t i=0; i<ny+1;++i) {
    cout << ptbins[i] << " ";
  }
  cout << endl;
  cout << "etabins : ";
  for(Int_t i=0; i<nx+1;++i) {
    cout << etabins[i] << " ";
  }
  cout << endl;

  //--------------------------------------------------------------------------------------------------------------
  // Produce Text file table
  //==============================================================================================================   

  ofstream txtfile;
  txtfile.open(outfname.Data());
  assert(txtfile.is_open());
    
  
  txtfile << " pT        ";
  txtfile << " eta           ";
  txtfile << "    MC T&P           ";
  txtfile << "    Data T&P         ";  
  txtfile << "    Scale factor     ";
  txtfile << endl;
  txtfile << "----------------------------------------------------------------------------------------------------------------------------------------------------";
  txtfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {
      txtfile << "[" << setw(4) << hMCEff->GetYaxis()->GetBinLowEdge(iy) << "," << setw(4) << hMCEff->GetYaxis()->GetBinLowEdge(iy+1) << "]";
      txtfile << "[" << setw(6) << hMCEff->GetXaxis()->GetBinLowEdge(ix) << "," << setw(6) << hMCEff->GetXaxis()->GetBinLowEdge(ix+1) << "]";      
      ios_base::fmtflags flags = txtfile.flags();
      txtfile.precision(4);
      
      Double_t mceff  = hMCEff->GetCellContent(ix,iy);
      Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
      Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
      txtfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);

      Double_t dataeff  = hDataEff->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh->GetCellContent(ix,iy);
      txtfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      
      Double_t scale     = (hDataEff->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      txtfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);

      txtfile << endl;
      txtfile.flags(flags);

      cout << "Set " << iy << " " << ix << " " << scale << endl;
      h2_results_selection->SetCellContent(iy,ix, scale);
      h2_results_selection->SetCellError(iy,ix, (scaleErrl + scaleErrh)/2);
      //fill overflow bins with the same values as last bin
      if (ix == nx) {
      cout << "Set " << iy << " " << ix+1 << " " << scale << endl;
        h2_results_selection->SetCellContent(iy,nx+1, scale);
        h2_results_selection->SetCellError(iy,nx+1, (scaleErrl + scaleErrh)/2);
      }
    }
    if (iy == ny) {
      for(Int_t ix=1; ix<=nx; ix++) {
        Double_t mceff  = hMCEff->GetCellContent(ix,iy);
        Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
        Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
        Double_t dataeff  = hDataEff->GetCellContent(ix,iy);
        Double_t dataerrl = hDataErrl->GetCellContent(ix,iy);
        Double_t dataerrh = hDataErrh->GetCellContent(ix,iy);
        Double_t scale     = (hDataEff->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
        Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
        Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
        
        cout << "Set " << iy+1 << " " << ix << " " << scale << endl;
        h2_results_selection->SetCellContent(iy+1,ix, scale);
        h2_results_selection->SetCellError(iy+1,ix, (scaleErrl + scaleErrh)/2);
        //fill overflow bins with the same values as last bin
        if (ix == nx) {
          cout << "Set " << iy+1 << " " << ix+1 << " " << scale << endl;
          h2_results_selection->SetCellContent(iy+1,nx+1, scale);
          h2_results_selection->SetCellError(iy+1,nx+1, (scaleErrl + scaleErrh)/2);
        }
      }
    }

    txtfile << endl;
  }
  txtfile.close();
  
  cout << outfname << " created!" << endl;
  



  //--------------------------------------------------------------------------------------------------------------
  // Create TEX table
  //==============================================================================================================   

  ofstream texfile;
  texfile.open((outputDir + "/eff_table.tex").c_str());
  assert(texfile.is_open());

  texfile << " \\begin{table}[!ht]" << endl;
  texfile << " \\begin{center} " << endl;
  texfile << " \\begin{tabular}{|c|c|c|c|}" << endl;
  texfile << " \\hline\n";


  texfile << " $p_{T}$ / $\\eta$ bin    &  Monte Carlo Efficiency    &  Data Efficiency   &  MC to Data Scale Factor \\\\  ";
  texfile << " \\hline           ";
  texfile << endl;
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {

      string binLabel = Form("$%5.1f < p_{T} \\le %5.1f$ , $%5.1f  \\le |\\eta| < %5.1f$", 
                             hMCEff->GetYaxis()->GetBinLowEdge(iy), hMCEff->GetYaxis()->GetBinLowEdge(iy+1),
                             hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("$%5.1f < p_{T} $ , $%5.1f  \\le |\\eta| < %5.1f$", 
                        hMCEff->GetYaxis()->GetBinLowEdge(iy), 
                        hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      texfile << binLabel;
      texfile << "   &   ";

      ios_base::fmtflags flags = texfile.flags();
      texfile.precision(4);
      
      Double_t mceff  = hMCEff->GetCellContent(ix,iy);
      Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
      Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
      texfile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      texfile << "   &   ";

      Double_t dataeff  = hDataEff->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh->GetCellContent(ix,iy);
      texfile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      texfile << "   &   ";
     
      Double_t scale     = (hDataEff->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      texfile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      texfile << "   \\\\   ";
      texfile << endl;
      texfile << "\\hline";
      texfile << endl;

    }
  }
  texfile << "\\end{tabular}" << endl;
  texfile << "\\caption{CAPTION.}" << endl;
  texfile << "\\label{tab:eff_ele_offline}" << endl;
  texfile << "\\end{center}" << endl;
  texfile << "\\end{table}" << endl;
  texfile.close();
  cout << outputDir + "/eff_table.tex" << " created!" << endl;




  //--------------------------------------------------------------------------------------------------------------
  // Create TWIKI table
  //==============================================================================================================   

  ofstream twikifile;
  twikifile.open((outputDir + "/eff_table.twiki").c_str());
  assert(twikifile.is_open());
  

  twikifile << "| *pT bin* | *eta bin* | *MC Efficiency* | *Data Efficiency* | *ScaleFactor* | \n";
  
  for(Int_t iy=1; iy<=ny; iy++) {
    for(Int_t ix=1; ix<=nx; ix++) {

      string binLabel = Form("| *%4.1f < pT \\le %4.1f* | *%3.1f  <= eta < %3.1f* | ", 
                             hMCEff->GetYaxis()->GetBinLowEdge(iy), hMCEff->GetYaxis()->GetBinLowEdge(iy+1),
                             hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
      if (iy == ny) {
        binLabel = Form("| *%4.1f < pT* | *%3.1f  <= eta < %3.1f* | ", 
                        hMCEff->GetYaxis()->GetBinLowEdge(iy), 
                        hMCEff->GetXaxis()->GetBinLowEdge(ix), hMCEff->GetXaxis()->GetBinLowEdge(ix+1));
      }
      
      twikifile << binLabel ;

      ios_base::fmtflags flags = twikifile.flags();
      twikifile.precision(4);
      
      Double_t mceff  = hMCEff->GetCellContent(ix,iy);
      Double_t mcerrl = hMCErrl->GetCellContent(ix,iy);
      Double_t mcerrh = hMCErrh->GetCellContent(ix,iy);
      twikifile << " " << setw(9) << fixed << mceff << " +/- " << TMath::Max(mcerrl,mcerrh);
      twikifile << " | ";

      Double_t dataeff  = hDataEff->GetCellContent(ix,iy);
      Double_t dataerrl = hDataErrl->GetCellContent(ix,iy);
      Double_t dataerrh = hDataErrh->GetCellContent(ix,iy);
      twikifile << " " << setw(9) << fixed << dataeff << " +/- " << TMath::Max(dataerrl,dataerrh);
      twikifile << " | ";
     
      Double_t scale     = (hDataEff->GetCellContent(ix,iy))/(hMCEff->GetCellContent(ix,iy));
      Double_t scaleErrl = scale*sqrt(mcerrl*mcerrl/mceff/mceff + dataerrl*dataerrl/dataeff/dataeff);
      Double_t scaleErrh = scale*sqrt(mcerrh*mcerrh/mceff/mceff + dataerrh*dataerrh/dataeff/dataeff);
      twikifile << " " << setw(9) << fixed << scale << " +/- " << TMath::Max(scaleErrl,scaleErrh);
      twikifile << " |\n";

    }
  }
  twikifile.close();
  cout << outputDir + "/eff_table.twiki" << " created!" << endl;





  outputFile->WriteTObject(h2_results_selection, h2_results_selection->GetName(), "WriteDelete");
  outputFile->Close();

  gBenchmark->Show("printMuonWPEff"); 
}
