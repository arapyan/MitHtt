//================================================================================================
//
//________________________________________________________________________________________________


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                   // access to gROOT, entry point to ROOT system
#include <TSystem.h>                 // interface to OS
#include <TStyle.h>                  // class to handle ROOT plotting style
#include <TFile.h>                   // file handle class
#include <TTree.h>                   // class to access ntuples
#include <TCanvas.h>                 // class for drawing
#include <TH1F.h>                    // 1D histograms
#include <TH2D.h>                    // 2D histograms
#include <TGraphAsymmErrors.h>       // graph class with asymmetric errors
#include <TEfficiency.h>             // class to handle efficiency calculations
#include <TPaletteAxis.h>            // class to handle palette axis (color scale for contour plot)
#include <iostream>                  // standard I/O
#include <iomanip>                   // functions to format standard I/O
#include <fstream>                   // functions for file I/O
#include <sstream>                   // class for parsing strings

#include "Common/CPlot.hh"           // helper class for plots
#include "Common/MitStyleRemix.hh"   // style settings for drawing
#endif


//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
void makeHTML(const TString outDir);

TGraphAsymmErrors* computeFakeRate1D(const TH1F* hpass, const TH1F* htotal);
void computeFakeRate2D(const TH2D *hpass, const TH2D* htotal, 
                       TH2D *hresult, TH2D* herrl, TH2D* herrh);

TGraphAsymmErrors* combine1D(const vector<TGraphAsymmErrors*> &effv, const vector<Double_t> &weightv);
void combine2D(const vector<TH2D*> &heffv, const vector<TH2D*> &herrlv, const vector<TH2D*> &herrhv, const vector<Double_t> &weightv,
               TH2D* hresult, TH2D* herrorl, TH2D* herrorh); 

//=== MAIN MACRO ================================================================================================= 

void computeFakeRate(const TString input,    // input file
                     const TString format,   // plot format
		     const Bool_t  doAbsEta  // bin in |eta|?
) {

  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 
  
  // bin edges for kinematic variables
  vector<Double_t> ptBinEdgesv;
  vector<Double_t> etaBinEdgesv;
  vector<Double_t> phiBinEdgesv;
  
  TString outputDir;
  
  vector<TString>  fnamev;
  vector<Double_t> weightv;
  vector<Int_t>    colorv;
  vector<TString>  labelv;

  //
  // parse .fo file
  //
  ifstream ifs;
  ifs.open(input.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') {
      state++;
      continue;
    }
    
    if(state==0) {
      outputDir = line;
    
    } else if(state==1) {
      stringstream ss(line);
      string fname;
      Double_t weight;
      Int_t color;
      ss >> fname >> weight >> color;
      string label = line.substr(line.find('@')+1);
      fnamev.push_back(fname);
      weightv.push_back(weight);
      colorv.push_back(color);
      labelv.push_back(label);
    
    } else if(state>=2 && state<=4) {
      Double_t edge;
      stringstream ss(line);
      ss >> edge;
      if(state==2)      { etaBinEdgesv.push_back(edge); }
      else if(state==3) { ptBinEdgesv.push_back(edge);  }
      else if(state==4) { phiBinEdgesv.push_back(edge); }
    } 
  }
  ifs.close();
  
  CPlot::sOutDir = outputDir + TString("/plots");

  Double_t yhigh1 = 0.3;
  Double_t yhigh2 = 0.6;
  

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================
  
  //  
  // Set up histograms
  //
  TH1F *hNumerPt0=0, *hNumerEta0=0, *hNumerPhi0=0;
  TH1F *hDenomPt0=0, *hDenomEta0=0, *hDenomPhi0=0;
    
  vector<TH1F*> hNumerPtv;
  vector<TH1F*> hNumerEtav;
  vector<TH1F*> hNumerPhiv;
  vector<TH1F*> hNumerNPVv;

  vector<TH1F*> hDenomPtv; 
  vector<TH1F*> hDenomEtav;
  vector<TH1F*> hDenomPhiv;
  vector<TH1F*> hDenomNPVv;
  
  vector<TGraphAsymmErrors*> effPtv; 
  vector<TGraphAsymmErrors*> effEtav;
  vector<TGraphAsymmErrors*> effPhiv;
  vector<TGraphAsymmErrors*> effNPVv;
  
  vector<TH2D*> hNumerEtaPtv, hDenomEtaPtv, hEffEtaPtv, hErrlEtaPtv, hErrhEtaPtv;  

  Int_t nbins = ptBinEdgesv.size()-1;
  Double_t *ptbinning = new Double_t[ptBinEdgesv.size()];
  for(UInt_t i=0; i<ptBinEdgesv.size(); i++) { ptbinning[i] = ptBinEdgesv[i]; }
  hNumerPt0 = new TH1F("hNumerPt0","",nbins,ptbinning);
  hDenomPt0 = new TH1F("hDenomPt0","",nbins,ptbinning);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    char hname[100];
    sprintf(hname,"hNumerPt_%i",ifile);       hNumerPtv.push_back(new TH1F(hname,"",nbins,ptbinning));       hNumerPtv[ifile]->Sumw2();
    sprintf(hname,"hDenomPt_%i",ifile);       hDenomPtv.push_back(new TH1F(hname,"",nbins,ptbinning));       hDenomPtv[ifile]->Sumw2();
  }
    
  nbins = etaBinEdgesv.size()-1;
  Double_t *etabinning = new Double_t[etaBinEdgesv.size()];
  for(UInt_t i=0; i<etaBinEdgesv.size(); i++) { etabinning[i] = etaBinEdgesv[i]; } 
  hNumerEta0 = new TH1F("hNumerEta0","",nbins,etabinning);
  hDenomEta0 = new TH1F("hDenomEta0","",nbins,etabinning);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    char hname[100];
    sprintf(hname,"hNumerEta_%i",ifile); hNumerEtav.push_back(new TH1F(hname,"",nbins,etabinning)); hNumerEtav[ifile]->Sumw2();
    sprintf(hname,"hDenomEta_%i",ifile); hDenomEtav.push_back(new TH1F(hname,"",nbins,etabinning)); hDenomEtav[ifile]->Sumw2();
  }
        
  nbins = phiBinEdgesv.size()-1;
  Double_t *phibinning = new Double_t[phiBinEdgesv.size()]; 
  for(UInt_t i=0; i<phiBinEdgesv.size(); i++) { phibinning[i] = phiBinEdgesv[i]; }
  hNumerPhi0 = new TH1F("hNumerPhi0","",nbins,phibinning);
  hDenomPhi0 = new TH1F("hDenomPhi0","",nbins,phibinning);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    char hname[100];
    sprintf(hname,"hNumerPhi_%i",ifile); hNumerPhiv.push_back(new TH1F(hname,"",nbins,phibinning));       hNumerPhiv[ifile]->Sumw2();
    sprintf(hname,"hDenomPhi_%i",ifile); hDenomPhiv.push_back(new TH1F(hname,"",nbins,phibinning));       hDenomPhiv[ifile]->Sumw2();
  }   
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    char hname[100];
    sprintf(hname,"hNumerNPV_%i",ifile); hNumerNPVv.push_back(new TH1F(hname,"",15,-0.5,14.5)); hNumerNPVv[ifile]->Sumw2();
    sprintf(hname,"hDenomNPV_%i",ifile); hDenomNPVv.push_back(new TH1F(hname,"",15,-0.5,14.5)); hDenomNPVv[ifile]->Sumw2();
  }
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    char hname[100];
    sprintf(hname,"hNumerEtaPt_%i",ifile); 
    hNumerEtaPtv.push_back(new TH2D(hname,"",etaBinEdgesv.size()-1,etabinning,ptBinEdgesv.size()-1,ptbinning));
    hNumerEtaPtv[ifile]->Sumw2();
    
    sprintf(hname,"hDenomEtaPt_%i",ifile); 
    hDenomEtaPtv.push_back(new TH2D(hname,"",etaBinEdgesv.size()-1,etabinning,ptBinEdgesv.size()-1,ptbinning));
    hDenomEtaPtv[ifile]->Sumw2();
    
    sprintf(hname,"hEffEtaPt_%i",ifile);
    hEffEtaPtv.push_back(new TH2D(hname,"",etaBinEdgesv.size()-1,etabinning,ptBinEdgesv.size()-1,ptbinning));
    hEffEtaPtv[ifile]->Sumw2();

    sprintf(hname,"hErrlEtaPt_%i",ifile);
    hErrlEtaPtv.push_back(new TH2D(hname,"",etaBinEdgesv.size()-1,etabinning,ptBinEdgesv.size()-1,ptbinning));
    hErrlEtaPtv[ifile]->Sumw2();
    
    sprintf(hname,"hErrhEtaPt_%i",ifile);
    hErrhEtaPtv.push_back(new TH2D(hname,"",etaBinEdgesv.size()-1,etabinning,ptBinEdgesv.size()-1,ptbinning));
    hErrhEtaPtv[ifile]->Sumw2();
  }  
  
  delete [] ptbinning;
  delete [] etabinning;
  delete [] phibinning;
  
  vector<TH1F*> hJetPt1v;
  vector<TH1F*> hJetPt2v;
  vector<TH1F*> hFOJetPt1v;
  vector<TH1F*> hFOJetPt2v;
  
  vector<TH1F*> hNumerJetPtv, hNumerFOJetPtv;
  vector<TH1F*> hDenomJetPtv, hDenomFOJetPtv;
  
  vector<TGraphAsymmErrors*> effJetPtv, effFOJetPtv;

  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    char hname[100];
    sprintf(hname,"hJetPt1_%i",ifile);   hJetPt1v.push_back(new TH1F(hname,"",50,0,150));   hJetPt1v[ifile]->Sumw2();
    sprintf(hname,"hJetPt2_%i",ifile);   hJetPt2v.push_back(new TH1F(hname,"",50,0,150));   hJetPt2v[ifile]->Sumw2();
    sprintf(hname,"hFOJetPt1_%i",ifile); hFOJetPt1v.push_back(new TH1F(hname,"",50,0,150)); hFOJetPt1v[ifile]->Sumw2();    
    sprintf(hname,"hFOJetPt2_%i",ifile); hFOJetPt2v.push_back(new TH1F(hname,"",50,0,150)); hFOJetPt2v[ifile]->Sumw2();
    
    sprintf(hname,"hNumerJetPt_%i",ifile);   hNumerJetPtv.push_back(new TH1F(hname,"",30,15,90));   hNumerJetPtv[ifile]->Sumw2();
    sprintf(hname,"hNumerFOJetPt_%i",ifile); hNumerFOJetPtv.push_back(new TH1F(hname,"",30,15,90)); hNumerFOJetPtv[ifile]->Sumw2();
    
    sprintf(hname,"hDenomJetPt_%i",ifile);   hDenomJetPtv.push_back(new TH1F(hname,"",30,15,90));   hDenomJetPtv[ifile]->Sumw2();    
    sprintf(hname,"hDenomFOJetPt_%i",ifile); hDenomFOJetPtv.push_back(new TH1F(hname,"",30,15,90)); hDenomFOJetPtv[ifile]->Sumw2();
  }

  Float_t fopt, foeta, fophi, fopass, fojpt;
  Float_t jpt;
  Int_t npv;
  
  TFile *infile=0;
  TTree *intree=0;
  
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]);
    intree = (TTree*)infile->Get("FO");
    
    intree->SetBranchAddress("pt",   &fopt);
    intree->SetBranchAddress("eta",  &foeta);
    intree->SetBranchAddress("phi",  &fophi);
    intree->SetBranchAddress("pass", &fopass);
    intree->SetBranchAddress("fojpt",&fojpt);
    intree->SetBranchAddress("jpt",  &jpt);
    intree->SetBranchAddress("npv",  &npv);
    
    for(UInt_t ientry=0; ientry<intree->GetEntries(); ientry++) {
      intree->GetEntry(ientry);
      
      if(fopt < ptBinEdgesv.front()   || fopt > ptBinEdgesv.back())   continue;
      if(doAbsEta) {
        if(fabs(foeta) < etaBinEdgesv.front() || fabs(foeta) > etaBinEdgesv.back()) continue;
      } else {
        if(foeta < etaBinEdgesv.front() || foeta > etaBinEdgesv.back()) continue;
      }
      if(fophi < phiBinEdgesv.front() || fophi > phiBinEdgesv.back()) continue;
      
      Bool_t pass = (fopass==1);
				   
      if(jpt>0) {
        hJetPt1v[ifile]->Fill(jpt,weightv[ifile]);
	hJetPt2v[ifile]->Fill(jpt,weightv[ifile]);	
	
	hDenomJetPtv[ifile]->Fill(jpt);
	if(pass) hNumerJetPtv[ifile]->Fill(jpt);
      }
      
      if(fojpt>0) {
        hFOJetPt1v[ifile]->Fill(fojpt,weightv[ifile]);
	hFOJetPt2v[ifile]->Fill(fojpt,weightv[ifile]);
    
        hDenomFOJetPtv[ifile]->Fill(fojpt);
	if(pass) hNumerFOJetPtv[ifile]->Fill(fojpt);
      }
                  
      hDenomPt0 ->Fill(fopt,weightv[ifile]);
      hDenomEta0->Fill(doAbsEta ? fabs(foeta) : foeta,weightv[ifile]);
      hDenomPhi0->Fill(fophi,weightv[ifile]);
      
      hDenomPtv[ifile] ->Fill(fopt);
      hDenomEtav[ifile]->Fill(doAbsEta ? fabs(foeta) : foeta);
      hDenomPhiv[ifile]->Fill(fophi);
      hDenomNPVv[ifile]->Fill(npv);
      hDenomEtaPtv[ifile]->Fill(doAbsEta ? fabs(foeta) : foeta,fopt);      
      
      if(pass) {
        hNumerPt0 ->Fill(fopt,weightv[ifile]);
        hNumerEta0->Fill(doAbsEta ? fabs(foeta) : foeta,weightv[ifile]);
        hNumerPhi0->Fill(fophi,weightv[ifile]);
	
	hNumerPtv[ifile] ->Fill(fopt);
        hNumerEtav[ifile]->Fill(doAbsEta ? fabs(foeta) : foeta);
        hNumerPhiv[ifile]->Fill(fophi);
	hNumerNPVv[ifile]->Fill(npv);
	hNumerEtaPtv[ifile]->Fill(doAbsEta ? fabs(foeta) : foeta,fopt);	 
      }

    }
    
    effPtv.push_back(computeFakeRate1D(hNumerPtv[ifile],hDenomPtv[ifile]));
    effEtav.push_back(computeFakeRate1D(hNumerEtav[ifile],hDenomEtav[ifile]));
    effPhiv.push_back(computeFakeRate1D(hNumerPhiv[ifile],hDenomPhiv[ifile]));
    
    effNPVv.push_back(computeFakeRate1D(hNumerNPVv[ifile],hDenomNPVv[ifile]));
    effJetPtv.push_back(computeFakeRate1D(hNumerJetPtv[ifile],hDenomJetPtv[ifile]));
    effFOJetPtv.push_back(computeFakeRate1D(hNumerFOJetPtv[ifile],hDenomFOJetPtv[ifile]));
    
    computeFakeRate2D(hNumerEtaPtv[ifile],hDenomEtaPtv[ifile],hEffEtaPtv[ifile],hErrlEtaPtv[ifile],hErrhEtaPtv[ifile]);
    
    delete infile;
    infile=0, intree=0;
  }
  
  TGraphAsymmErrors *frPt        = combine1D(effPtv,       weightv); frPt->SetName("frPtv");
  TGraphAsymmErrors *frEta       = combine1D(effEtav,      weightv); frEta->SetName("frEta");
  TGraphAsymmErrors *frPhi       = combine1D(effPhiv,      weightv); frPhi->SetName("frPhi");
  TGraphAsymmErrors *frNPV       = combine1D(effNPVv,      weightv); frNPV->SetName("frNPV");
  TGraphAsymmErrors *frJetPt     = combine1D(effJetPtv,    weightv); frJetPt->SetName("frJetPt");
  TGraphAsymmErrors *frFOJetPt   = combine1D(effFOJetPtv,  weightv); frFOJetPt->SetName("frFOJetPt");
  
  TH2D *frEtaPt   = (TH2D*)hEffEtaPtv[0]->Clone("frEtaPt");
  TH2D *errlEtaPt = (TH2D*)hEffEtaPtv[0]->Clone("errlEtaPt");
  TH2D *errhEtaPt = (TH2D*)hEffEtaPtv[0]->Clone("errhEtaPt");
  combine2D(hEffEtaPtv,hErrlEtaPtv,hErrhEtaPtv,weightv,frEtaPt,errlEtaPt,errhEtaPt);

  //--------------------------------------------------------------------------------------------------------------
  // Make plots 
  //==============================================================================================================  
  TCanvas *c = MakeCanvas("c","c",800,600);
  char ylabel[100];
  TString etalabel = (doAbsEta) ? "|#eta|" : "#eta";

  //
  // Denominator distributions
  //
  sprintf(ylabel,"Events / %.1f GeV/c",hDenomPtv[0]->GetBinWidth(1));
  CPlot plotDenomPt("denompt","","denominator p_{T} [GeV/c]",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    hDenomPtv[ifile]->Scale(weightv[ifile]);
    plotDenomPt.AddToStack(hDenomPtv[ifile],labelv[ifile],colorv[ifile]);
  }
  plotDenomPt.TransLegend(0.1,0);
  plotDenomPt.Draw(c,kTRUE,format);
 
  plotDenomPt.SetName("denomptlog");
  plotDenomPt.SetLogy();
  plotDenomPt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f",hDenomEtav[0]->GetBinWidth(1));
  CPlot plotDenomEta("denometa","","denominator #eta",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    hDenomEtav[ifile]->Scale(weightv[ifile]);
    plotDenomEta.AddToStack(hDenomEtav[ifile],labelv[ifile],colorv[ifile]);
  }
  plotDenomEta.SetYRange(0,2.0*(plotDenomEta.GetStack()->GetMaximum()));
  plotDenomEta.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f",hDenomPhiv[0]->GetBinWidth(1));
  CPlot plotDenomPhi("denomphi","","denominator #phi",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    hDenomPhiv[ifile]->Scale(weightv[ifile]);
    plotDenomPhi.AddToStack(hDenomPhiv[ifile],labelv[ifile],colorv[ifile]);
  }
  plotDenomPhi.SetYRange(0,1.8*(plotDenomPhi.GetStack()->GetMaximum()));
  plotDenomPhi.Draw(c,kTRUE,format);

  //
  // Numerator distributions
  //
  sprintf(ylabel,"Events / %.1f GeV/c",hNumerPtv[0]->GetBinWidth(1));
  CPlot plotNumerPt("numerpt","","numerator p_{T} [GeV/c]",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    hNumerPtv[ifile]->Scale(weightv[ifile]);
    plotNumerPt.AddToStack(hNumerPtv[ifile],labelv[ifile],colorv[ifile]);
  }
  plotNumerPt.TransLegend(0.1,0);
  plotNumerPt.Draw(c,kTRUE,format);

  plotNumerPt.SetName("numerptlog");
  plotNumerPt.SetLogy();
  plotNumerPt.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f",hNumerEtav[0]->GetBinWidth(1));
  CPlot plotNumerEta("numereta","","numerator #eta",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    hNumerEtav[ifile]->Scale(weightv[ifile]);
    plotNumerEta.AddToStack(hNumerEtav[ifile],labelv[ifile],colorv[ifile]);
  }
  plotNumerEta.SetYRange(0,2.0*(plotNumerEta.GetStack()->GetMaximum()));
  plotNumerEta.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f",hNumerPhiv[0]->GetBinWidth(1));
  CPlot plotNumerPhi("numerphi","","numerator #phi",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    hNumerPhiv[ifile]->Scale(weightv[ifile]);
    plotNumerPhi.AddToStack(hNumerPhiv[ifile],labelv[ifile],colorv[ifile]);
  }
  plotNumerPhi.SetYRange(0,1.8*(plotNumerPhi.GetStack()->GetMaximum()));
  plotNumerPhi.Draw(c,kTRUE,format);


  //
  // Fakeable object kinematic distributions
  //
  sprintf(ylabel,"Events / %.1f GeV/c",hDenomPt0->GetBinWidth(1));
  CPlot plotPt("pt","","p_{T} [GeV/c]",ylabel);
  plotPt.AddHist1D(hDenomPt0,"loose","hist",kBlue,7);
  plotPt.AddHist1D(hNumerPt0,"tight","hist",kRed);
  plotPt.TransLegend(0.1,0);
  plotPt.Draw(c,kTRUE,format);
 
  plotPt.SetName("ptlog");
  plotPt.SetYRange(1e-4*(hDenomPt0->GetMaximum()),10*(hDenomPt0->GetMaximum()));
  plotPt.SetLogy();
  plotPt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f",hDenomEta0->GetBinWidth(1));
  CPlot plotEta("eta","","#eta",ylabel);
  plotEta.AddHist1D(hDenomEta0,"loose","hist",kBlue,7);
  plotEta.AddHist1D(hNumerEta0,"tight","hist",kRed);
  plotEta.SetYRange(0,2.0*(hDenomEta0->GetMaximum()));
  plotEta.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f",hDenomPhi0->GetBinWidth(1));
  CPlot plotPhi("phi","","#phi",ylabel);
  plotPhi.AddHist1D(hDenomPhi0,"loose","hist",kBlue,7);
  plotPhi.AddHist1D(hNumerPhi0,"tight","hist",kRed);
  plotPhi.SetYRange(0,1.8*(hDenomPhi0->GetMaximum()));
  plotPhi.Draw(c,kTRUE,format);
  
  //
  // Jet distributions
  //
  sprintf(ylabel,"Events / %.1f GeV/c",hJetPt1v[0]->GetBinWidth(1));
  CPlot plotJetPt1("jetpt","","leading jet p_{T} [GeV/c]",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++)
    plotJetPt1.AddToStack(hJetPt1v[ifile],labelv[ifile],colorv[ifile]);
  if(fnamev.size()>4)
    plotJetPt1.SetLegend(0.65,0.55,0.92,0.9);
  plotJetPt1.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hJetPt2v[0]->GetBinWidth(1));
  CPlot plotJetPt2("jetptlog","","leading jet p_{T} [GeV/c]",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++)
    plotJetPt2.AddToStack(hJetPt2v[ifile],labelv[ifile],colorv[ifile]);
  if(fnamev.size()>4)
    plotJetPt2.SetLegend(0.65,0.55,0.92,0.9);
  plotJetPt2.SetYRange(1e-5*(plotJetPt2.GetStack()->GetMaximum()),10*(plotJetPt2.GetStack()->GetMaximum()));  
  plotJetPt2.SetLogy();
  plotJetPt2.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f GeV/c",hFOJetPt1v[0]->GetBinWidth(1));
  CPlot plotFOJetPt1("fojetpt","","FO jet p_{T} [GeV/c]",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++)
    plotFOJetPt1.AddToStack(hFOJetPt1v[ifile],labelv[ifile],colorv[ifile]);
  if(fnamev.size()>4)
    plotFOJetPt1.SetLegend(0.65,0.55,0.92,0.9);
  plotFOJetPt1.Draw(c,kTRUE,format);
  
  sprintf(ylabel,"Events / %.1f GeV/c",hFOJetPt2v[0]->GetBinWidth(1));
  CPlot plotFOJetPt2("fojetptlog","","FO jet p_{T} [GeV/c]",ylabel);
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++)
    plotFOJetPt2.AddToStack(hFOJetPt2v[ifile],labelv[ifile],colorv[ifile]);
  if(fnamev.size()>4)
    plotFOJetPt2.SetLegend(0.65,0.55,0.92,0.9);
  plotFOJetPt2.SetYRange(1e-5*(plotFOJetPt2.GetStack()->GetMaximum()),10*(plotFOJetPt2.GetStack()->GetMaximum()));
  plotFOJetPt2.SetLogy();
  plotFOJetPt2.Draw(c,kTRUE,format);
        
  //
  // Fake rate plots
  //       
  CPlot plotFRPt("frpt","","p_{T} [GeV/c]","#varepsilon_{fake}");
  plotFRPt.AddGraph(frPt,"");
  plotFRPt.SetYRange(0,yhigh1);
  plotFRPt.Draw(c,kTRUE,format);
  
  CPlot plotFREta("freta","",etalabel,"#varepsilon_{fake}");
  plotFREta.AddGraph(frEta,"");
  plotFREta.SetYRange(0,yhigh1);
  plotFREta.Draw(c,kTRUE,format);
  
  CPlot plotFRPhi("frphi","","#phi","#varepsilon_{fake}");
  plotFRPhi.AddGraph(frPhi,"");
  plotFRPhi.SetYRange(0,yhigh1);
  plotFRPhi.Draw(c,kTRUE,format);

  CPlot plotFRNPV("frnpv","","N_{PV}","#varepsilon_{fake}");
  plotFRNPV.AddGraph(frNPV,"");
  plotFRNPV.SetYRange(0,yhigh1);
  plotFRNPV.Draw(c,kTRUE,format);

  CPlot plotFRJetPt("frjetpt","","leading jet p_{T} [GeV/c]","#varepsilon_{fake}");
  plotFRJetPt.AddGraph(frJetPt,"");
  plotFRJetPt.SetYRange(0,yhigh2);
  plotFRJetPt.Draw(c,kTRUE,format);
  
  CPlot plotFRFOJetPt("frfojetpt","","FO jet p_{T} [GeV/c]","#varepsilon_{fake}");
  plotFRFOJetPt.AddGraph(frFOJetPt,"");
  plotFRFOJetPt.SetYRange(0,yhigh2);
  plotFRFOJetPt.Draw(c,kTRUE,format);  
      
  gStyle->SetPalette(1);
  c->SetRightMargin(0.15);
  c->SetLeftMargin(0.15);  
  frEtaPt->SetTitleOffset(1.2,"Y");
//  TPaletteAxis *paxis = (TPaletteAxis*)frEtaPt->GetListOfFunctions()->FindObject("palette");
//  paxis->SetX1NDC(0.87);
//  paxis->SetX2NDC(0.92);
  CPlot plotFRPtEta("frpteta","",etalabel,"p_{T} [GeV/c]");
  plotFRPtEta.AddHist2D(frEtaPt,"COLZ");
  plotFRPtEta.Draw(c,kTRUE,format); 


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
    
  TFile outfile(outputDir + TString("/fr.root"),"RECREATE");
  frPt->Write();
  frEta->Write();
  frPhi->Write();
  frEtaPt->Write();
  errlEtaPt->Write();
  errhEtaPt->Write();
  frNPV->Write();
  frJetPt->Write();
  frFOJetPt->Write();
  outfile.Close();
  
  makeHTML(outputDir);
    
  cout << " <> Output saved in " << outputDir << "/" << endl;    
  cout << endl;
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
TGraphAsymmErrors* computeFakeRate1D(const TH1F* hpass, const TH1F* htotal)
{
  const Int_t nbins = htotal->GetNbinsX();
  Double_t xval[nbins], xerr[nbins], yval[nbins], yerrl[nbins], yerrh[nbins];
  for(Int_t ibin=1; ibin<=nbins; ibin++) {
    xval[ibin-1] = htotal->GetBinCenter(ibin);
    xerr[ibin-1] = 0.5*htotal->GetBinWidth(ibin);
  
    Int_t total  = htotal->GetBinContent(ibin);
    Int_t passed = hpass->GetBinContent(ibin);
    yval[ibin-1]   = (total>0) ? (Double_t)passed/(Double_t)total : 0;
    yerrl[ibin-1]  = (total>0) ? yval[ibin-1] - TEfficiency::ClopperPearson(total,passed,0.68269,kFALSE) : 0;
    yerrh[ibin-1]  = (total>0) ? TEfficiency::ClopperPearson(total,passed,0.68269,kTRUE) - yval[ibin-1]  : 0;
  }
  
  return new TGraphAsymmErrors(nbins,xval,yval,xerr,xerr,yerrl,yerrh);
}

//--------------------------------------------------------------------------------------------------
void computeFakeRate2D(const TH2D *hpass, const TH2D* htotal, 
                       TH2D *hresult, TH2D* herrl, TH2D* herrh)
{
  assert(hresult);
  assert(herrl);
  assert(herrh);
  
  const Int_t nbinsx = htotal->GetNbinsX();
  const Int_t nbinsy = htotal->GetNbinsY();
  for(Int_t ix=1; ix<=nbinsx; ix++) {
    for(Int_t iy=1; iy<=nbinsy; iy++) {
      Int_t total  = htotal->GetCellContent(ix,iy);
      Int_t passed = hpass->GetCellContent(ix,iy);
      Double_t eff  = (total>0) ? (Double_t)passed/(Double_t)total : 0;
      Double_t errl = (total>0) ? eff - TEfficiency::ClopperPearson(total,passed,0.68269,kFALSE) : 0;
      Double_t errh = (total>0) ? TEfficiency::ClopperPearson(total,passed,0.68269,kTRUE) - eff  : 0;
      hresult->SetCellContent(ix,iy,eff);
      herrl  ->SetCellContent(ix,iy,errl);
      herrh  ->SetCellContent(ix,iy,errh);
    }
  }
}

//--------------------------------------------------------------------------------------------------
TGraphAsymmErrors* combine1D(const vector<TGraphAsymmErrors*> &effv, const vector<Double_t> &weightv)
{
  const Int_t npts = effv[0]->GetN();
  Double_t yval[npts], yerrl[npts], yerrh[npts];
  for(Int_t i=0; i<npts; i++) {
    yval[i]=0;
    Double_t varl=0, varh=0;
    for(UInt_t ifile=0; ifile<effv.size(); ifile++) {
      yval[i] += weightv[ifile]*(effv[ifile]->GetY()[i]);
      varl += weightv[ifile]*weightv[ifile]*(effv[ifile]->GetErrorYlow(i))*(effv[ifile]->GetErrorYlow(i));
      varh += weightv[ifile]*weightv[ifile]*(effv[ifile]->GetErrorYhigh(i))*(effv[ifile]->GetErrorYhigh(i));
    }
    yerrl[i] = sqrt(varl);
    yerrh[i] = sqrt(varh);
  }
  
  return new TGraphAsymmErrors(npts,effv[0]->GetX(),yval,effv[0]->GetEXlow(),effv[0]->GetEXhigh(),yerrl,yerrh);
}

//--------------------------------------------------------------------------------------------------
void combine2D(const vector<TH2D*> &heffv, const vector<TH2D*> &herrlv, const vector<TH2D*> &herrhv, const vector<Double_t> &weightv,
               TH2D* hresult, TH2D* herrorl, TH2D* herrorh)
{
  assert(hresult);
  assert(herrorl);
  assert(herrorh);
  
  const Int_t nbinsx = heffv[0]->GetNbinsX();
  const Int_t nbinsy = heffv[0]->GetNbinsY();
  for(Int_t ix=1; ix<=nbinsx; ix++) {
    for(Int_t iy=1; iy<=nbinsy; iy++) {
      Double_t fr=0, varl=0, varh=0;
      for(UInt_t ifile=0; ifile<heffv.size(); ifile++) {
        fr += weightv[ifile]*heffv[ifile]->GetCellContent(ix,iy);
	varl += weightv[ifile]*weightv[ifile]*(herrlv[ifile]->GetCellContent(ix,iy))*(herrlv[ifile]->GetCellContent(ix,iy));
	varh += weightv[ifile]*weightv[ifile]*(herrhv[ifile]->GetCellContent(ix,iy))*(herrhv[ifile]->GetCellContent(ix,iy));
      }
      hresult->SetCellContent(ix,iy,fr);
      herrorl->SetCellContent(ix,iy,sqrt(varl));
      herrorh->SetCellContent(ix,iy,sqrt(varh));     
    }
  }  
}
	       		       
//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/fr.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;

  htmlfile << "<body bgcolor=\"EEEEEE\">" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/frpteta.png\"><img src=\"plots/frpteta.png\" alt=\"plots/frpteta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/frpt.png\"><img src=\"plots/frpt.png\" alt=\"plots/frpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/freta.png\"><img src=\"plots/freta.png\" alt=\"plots/freta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/frphi.png\"><img src=\"plots/frphi.png\" alt=\"plots/frphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/frnpv.png\"><img src=\"plots/frnpv.png\" alt=\"plots/frnpv.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/frptnpv.png\"><img src=\"plots/frptnpv.png\" alt=\"plots/frptnpv.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/fretanpv.png\"><img src=\"plots/fretanpv.png\" alt=\"plots/fretanpv.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/frphinpv.png\"><img src=\"plots/frphinpv.png\" alt=\"plots/frphinpv.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;    
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/denompt.png\"><img src=\"plots/denompt.png\" alt=\"plots/denompt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/denomptlog.png\"><img src=\"plots/denomptlog.png\" alt=\"plots/denomptlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/denometa.png\"><img src=\"plots/denometa.png\" alt=\"plots/denometa.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/denomphi.png\"><img src=\"plots/denomphi.png\" alt=\"plots/denomphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/numerpt.png\"><img src=\"plots/numerpt.png\" alt=\"plots/numerpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/numerptlog.png\"><img src=\"plots/numerptlog.png\" alt=\"plots/numerptlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/numereta.png\"><img src=\"plots/numereta.png\" alt=\"plots/numereta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/numerphi.png\"><img src=\"plots/numerphi.png\" alt=\"plots/numerphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt.png\"><img src=\"plots/pt.png\" alt=\"plots/pt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ptlog.png\"><img src=\"plots/ptlog.png\" alt=\"plots/ptlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta.png\"><img src=\"plots/eta.png\" alt=\"plots/eta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi.png\"><img src=\"plots/phi.png\" alt=\"plots/phi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl; 
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetpt.png\"><img src=\"plots/jetpt.png\" alt=\"plots/jetpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetptlog.png\"><img src=\"plots/jetptlog.png\" alt=\"plots/jetptlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/frjetpt.png\"><img src=\"plots/frjetpt.png\" alt=\"plots/frjetpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/fojetpt.png\"><img src=\"plots/fojetpt.png\" alt=\"plots/fojetpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/fojetptlog.png\"><img src=\"plots/fojetptlog.png\" alt=\"plots/fojetptlog.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/frfojetpt.png\"><img src=\"plots/frfojetpt.png\" alt=\"plots/frfojetpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"></td>" << endl;
  htmlfile << "</tr>" << endl;   
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
    
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close(); 
}  
