//================================================================================================
//
//
//________________________________________________________________________________________________

#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector3.h>               // 3D vector class
#include <TRegexp.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "Common/CPlot.hh"          // helper class for plots
#include "Common/MitStyleRemix.hh"  // style settings for drawing
#include "Common/MyTools.hh"        // miscellaneous helper functions
#include "Common/CSample.hh"        // helper class for organizing input ntuple files

// define structure for output ntuple
#include "EmuData.hh"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
void makeHTML(const TString outDir);


//=== MAIN MACRO =================================================================================================

void plotEmu(const TString  conf,         // input file
             const TString  ntupleDir,    // directory of input ntuples
	     const TString  outputDir,    // output directory
             const TString  format,       // plot file format
	     const Double_t lumi          // luminosity (pb^-1)
	     ) {  
  gBenchmark->Start("plotEmu");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  vector<TString>  snamev;    // sample name (for output file)  
  vector<CSample*> samplev;   // data/MC samples
    
  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    if(line[0]=='$') {
      samplev.push_back(new CSample());
      stringstream ss(line);
      string chr;
      string sname;
      Int_t color;
      ss >> chr >> sname >> color;
      string label = line.substr(line.find('@')+1);
      snamev.push_back(sname);
      samplev.back()->label = label;
      samplev.back()->color = color;
      continue;
    }
    
    if(state==0) {  // define data sample
      string fname;
      Double_t xsec;
      string json;
      Int_t type;
      stringstream ss(line);
      ss >> fname >> type >> xsec >> json;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(type);
      samplev.back()->xsecv.push_back(xsec);
      samplev.back()->jsonv.push_back(json);
    
    } else if(state==1) {  // define MC samples
      string fname;
      Double_t xsec;
      stringstream ss(line);
      ss >> fname >> xsec;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(0);
      samplev.back()->xsecv.push_back(xsec);
    }
  }
  ifs.close();

  CPlot::sOutDir = outputDir + TString("/plots");
  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };  // final state type

  //
  // scale factors
  //
  const Double_t kCat_nob  = 0.9993;
  const Double_t kCat_vbf  = 0.7949;
  const Double_t kCat_b    = 0.9958;
  
  //
  // VBF:
  //  * mjj cut
  //  * dEta(j,j) cut
  const Double_t mjjMin   = 450;
  const Double_t dEtaMin  = 3.5;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 
  
  const Double_t pi = 3.14159265358979;
  
  Bool_t hasData = (samplev[0]->fnamev.size()>0);
  
  //
  // Set up histograms
  //

  vector<TH2F*> hTrigEffv, hTrigEffEntriesv;
  vector<TH2F*> hTrigEffvMuPtv, hTrigEffvElePtv;

  // after lepton selection
  vector<TH1F*> hMassv, hMassLv, hMetv, hMetRawv;
  vector<TH2F*> hMassVpmissv;
  vector<TH1F*> hProjMetv, hProjVisv, hProjVarv,hRawProjVarv;
  vector<TH1F*> hNjetsv, hNbjetsv;
  vector<TH1F*> hBdiscrv;
  vector<TH1F*> hDPhiv, hMtv, hPtv;
  vector<TH1F*> hMetDPhiv;
  vector<TH1F*> hPt1v, hEta1v, hPhi1v;   	// leading lepton
  vector<TH1F*> hPt2v, hEta2v, hPhi2v;   	// trailing lepton
  vector<TH1F*> hPtMuv, hEtaMuv, hPhiMuv;  	// muon
  vector<TH1F*> hPtElev, hEtaElev, hPhiElev; 	// electron
  vector<TH1F*> hJetPt1v, hJetEta1v; 		// leading jet
  vector<TH1F*> hJetPt2v, hJetEta2v; 		// second pt jet
  vector<TH1F*> hJetDPhiv;           		// dphi of first two jets
  vector<TH1F*> hMjjv, hDEtav, hEtaProdv;       // kinematics of first two jets
  vector<TH1F*> hBJetPtv, hBJetEtav, hBJetPhiv; // leading b-jet (can be same as leading two jets)
  vector<TH1F*> hNPVv;                          // primary vertexes

  // inclusive
  vector<TH1F*> hMass_iv;
  vector<TH1F*> hMassL_iv;

  // Class 1: one or fewer pt 30 jets, no b-tag jet
  vector<TH1F*> hMass_nobv;
  vector<TH1F*> hMassL_nobv;

  // Class 2a: exactly two jets, vbf
  vector<TH1F*> hMass_vbfv;
  vector<TH1F*> hMassL_vbfv;

  // Class 2b: one or fewer pt 30 jets, at least on b-tag jet
  vector<TH1F*> hMass_bv;
  vector<TH1F*> hMassL_bv;

  vector<Double_t> nSelv;
  vector<Double_t> nSelVarv;
  vector<Double_t> nSel_iv;
  vector<Double_t> nSelVar_iv;
  vector<Double_t> nSel_nobv;
  vector<Double_t> nSelVar_nobv;
  vector<Double_t> nSel_vbfv;
  vector<Double_t> nSelVar_vbfv;
  vector<Double_t> nSel_bv;
  vector<Double_t> nSelVar_bv;
  
  vector<vector <Double_t> > countervv;

  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    // after lepton selection
    sprintf(hname,"hTrigEff_%i",isam);       hTrigEffv.push_back(new TH2F(hname,"",20,0,40,20,0,40));       hTrigEffv[isam]->Sumw2();
    sprintf(hname,"hTrigEffEntries_%i",isam);hTrigEffEntriesv.push_back(new TH2F(hname,"",20,0,40,20,0,40));hTrigEffEntriesv[isam]->Sumw2();
    sprintf(hname,"hTrigEffvMuPt_%i",isam);  hTrigEffvMuPtv.push_back(new TH2F(hname,"",20,10,40,20,0,1));  hTrigEffvMuPtv[isam]->Sumw2();
    sprintf(hname,"hTrigEffvElePt_%i",isam); hTrigEffvElePtv.push_back(new TH2F(hname,"",20,10,40,20,0,1)); hTrigEffvElePtv[isam]->Sumw2();

    sprintf(hname,"hMassVpmiss_%i",isam);    hMassVpmissv.push_back(new TH2F(hname,"",20,0,200,50,-100,100));     hMassVpmissv[isam]->Sumw2();

    sprintf(hname,"hMass_%i",isam);    hMassv.push_back(new TH1F(hname,"",20,0,200));     hMassv[isam]->Sumw2();
    sprintf(hname,"hMassL_%i",isam);   hMassLv.push_back(new TH1F(hname,"",50,0,500));     hMassLv[isam]->Sumw2();
    sprintf(hname,"hMet_%i",isam);     hMetv.push_back(new TH1F(hname,"",30,0,150));      hMetv[isam]->Sumw2();
    sprintf(hname,"hMetRaw_%i",isam);  hMetRawv.push_back(new TH1F(hname,"",30,0,150));   hMetRawv[isam]->Sumw2();
    sprintf(hname,"hProjMet_%i",isam);       hProjMetv.push_back(new TH1F(hname,"",50,-100,100));     hProjMetv[isam]->Sumw2();
    sprintf(hname,"hProjVis_%i",isam);       hProjVisv.push_back(new TH1F(hname,"",50,0,200));        hProjVisv[isam]->Sumw2();
    sprintf(hname,"hProjVar_%i",isam);       hProjVarv.push_back(new TH1F(hname,"",50,-100,200));     hProjVarv[isam]->Sumw2();
    sprintf(hname,"hRawProjVar_%i",isam);    hRawProjVarv.push_back(new TH1F(hname,"",50,-100,200));  hRawProjVarv[isam]->Sumw2();
    sprintf(hname,"hNjets_%i",isam);         hNjetsv.push_back(new TH1F(hname,"",5,-0.5,4.5));      hNjetsv[isam]->Sumw2();
    sprintf(hname,"hNbjets_%i",isam);        hNbjetsv.push_back(new TH1F(hname,"",5,-0.5,4.5));     hNbjetsv[isam]->Sumw2();
    sprintf(hname,"hBdiscr_%i",isam);        hBdiscrv.push_back(new TH1F(hname,"",50,0,15));          hBdiscrv[isam]->Sumw2();
    sprintf(hname,"hDPhi_%i",isam);    hDPhiv.push_back(new TH1F(hname,"",18,0,180));       hDPhiv[isam]->Sumw2();
    sprintf(hname,"hMt_%i",isam);      hMtv.push_back(new TH1F(hname,"",30,0,210));         hMtv[isam]->Sumw2();
    sprintf(hname,"hPt_%i",isam);      hPtv.push_back(new TH1F(hname,"",30,0,120));         hPtv[isam]->Sumw2();
    sprintf(hname,"hMetDPhi_%i",isam); hMetDPhiv.push_back(new TH1F(hname,"",30,0,180));    hMetDPhiv[isam]->Sumw2();
    sprintf(hname,"hPt1_%i",isam);     hPt1v.push_back(new TH1F(hname,"",30,0,150));        hPt1v[isam]->Sumw2();
    sprintf(hname,"hEta1_%i",isam);    hEta1v.push_back(new TH1F(hname,"",30,-3,3));        hEta1v[isam]->Sumw2();
    sprintf(hname,"hPhi1_%i",isam);    hPhi1v.push_back(new TH1F(hname,"",20,-3.2,3.2));    hPhi1v[isam]->Sumw2();  
    sprintf(hname,"hPt2_%i",isam);     hPt2v.push_back(new TH1F(hname,"",30,0,150));        hPt2v[isam]->Sumw2();
    sprintf(hname,"hEta2_%i",isam);    hEta2v.push_back(new TH1F(hname,"",30,-3,3));        hEta2v[isam]->Sumw2();
    sprintf(hname,"hPhi2_%i",isam);    hPhi2v.push_back(new TH1F(hname,"",20,-3.2,3.2));    hPhi2v[isam]->Sumw2();
    sprintf(hname,"hPtMu_%i",isam);    hPtMuv.push_back(new TH1F(hname,"",40,0,100));       hPtMuv[isam]->Sumw2();
    sprintf(hname,"hEtaMu_%i",isam);   hEtaMuv.push_back(new TH1F(hname,"",30,-3,3));       hEtaMuv[isam]->Sumw2();
    sprintf(hname,"hPhiMu_%i",isam);   hPhiMuv.push_back(new TH1F(hname,"",20,-3.2,3.2));   hPhiMuv[isam]->Sumw2();  
    sprintf(hname,"hPtEle_%i",isam);   hPtElev.push_back(new TH1F(hname,"",40,0,100));      hPtElev[isam]->Sumw2();
    sprintf(hname,"hEtaEle_%i",isam);  hEtaElev.push_back(new TH1F(hname,"",30,-3,3));      hEtaElev[isam]->Sumw2();
    sprintf(hname,"hPhiEle_%i",isam);  hPhiElev.push_back(new TH1F(hname,"",20,-3.2,3.2));  hPhiElev[isam]->Sumw2();
    sprintf(hname,"hJetPt1_%i",isam);  hJetPt1v.push_back(new TH1F(hname,"",30,0,300));     hJetPt1v[isam]->Sumw2();
    sprintf(hname,"hJetEta1_%i",isam); hJetEta1v.push_back(new TH1F(hname,"",20,-5,5));     hJetEta1v[isam]->Sumw2();
    sprintf(hname,"hJetPt2_%i",isam);  hJetPt2v.push_back(new TH1F(hname,"",30,0,300));     hJetPt2v[isam]->Sumw2();
    sprintf(hname,"hJetEta2_%i",isam); hJetEta2v.push_back(new TH1F(hname,"",20,-5,5));     hJetEta2v[isam]->Sumw2();
    sprintf(hname,"hJetDPhi_%i",isam); hJetDPhiv.push_back(new TH1F(hname,"",30,0,180));    hJetDPhiv[isam]->Sumw2();
    sprintf(hname,"hMjj_%i",isam);     hMjjv.push_back(new TH1F(hname,"",25,200,1200));     hMjjv[isam]->Sumw2();
    sprintf(hname,"hDEta_%i",isam);    hDEtav.push_back(new TH1F(hname,"",20,0,8));         hDEtav[isam]->Sumw2();
    sprintf(hname,"hEtaProd_%i",isam); hEtaProdv.push_back(new TH1F(hname,"",30,-7.5,7.5)); hEtaProdv[isam]->Sumw2();
    sprintf(hname,"hBJetPt_%i",isam);  hBJetPtv.push_back(new TH1F(hname,"",30,0,150));     hBJetPtv[isam]->Sumw2();
    sprintf(hname,"hBJetEta_%i",isam); hBJetEtav.push_back(new TH1F(hname,"",30,-3,3));     hBJetEtav[isam]->Sumw2();
    sprintf(hname,"hBJetPhi_%i",isam); hBJetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hBJetPhiv[isam]->Sumw2();
    sprintf(hname,"hNPV_%i",isam);     hNPVv.push_back(new TH1F(hname,"",20,-0.5,19.5));    hNPVv[isam]->Sumw2();

    // inclusive
    sprintf(hname,"hMass_i_%i",isam);    hMass_iv.push_back(new TH1F(hname,"",20,0,200));     hMass_iv[isam]->Sumw2();
    sprintf(hname,"hMassL_i_%i",isam);   hMassL_iv.push_back(new TH1F(hname,"",50,0,500));    hMassL_iv[isam]->Sumw2();

    // no b-tag
    sprintf(hname,"hMass_nob_%i",isam);    hMass_nobv.push_back(new TH1F(hname,"",20,0,200));     hMass_nobv[isam]->Sumw2();
    sprintf(hname,"hMassL_nob_%i",isam);   hMassL_nobv.push_back(new TH1F(hname,"",50,0,500));    hMassL_nobv[isam]->Sumw2();

    // vbf
    sprintf(hname,"hMass_vbf_%i",isam);    hMass_vbfv.push_back(new TH1F(hname,"",20,0,200));     hMass_vbfv[isam]->Sumw2();
    sprintf(hname,"hMassL_vbf_%i",isam);   hMassL_vbfv.push_back(new TH1F(hname,"",50,0,500));    hMassL_vbfv[isam]->Sumw2();

    // b-tag
    sprintf(hname,"hMass_b_%i",isam);    hMass_bv.push_back(new TH1F(hname,"",20,0,200));     hMass_bv[isam]->Sumw2();
    sprintf(hname,"hMassL_b_%i",isam);   hMassL_bv.push_back(new TH1F(hname,"",50,0,500));    hMassL_bv[isam]->Sumw2();

    nSelv.push_back(0);
    nSelVarv.push_back(0);
    nSel_iv.push_back(0);
    nSelVar_iv.push_back(0);
    nSel_nobv.push_back(0);
    nSelVar_nobv.push_back(0);
    nSel_vbfv.push_back(0);
    nSelVar_vbfv.push_back(0);
    nSel_bv.push_back(0);
    nSelVar_bv.push_back(0);

    vector<Double_t> *counter = new vector<Double_t>;
    countervv.push_back(*counter);
    for(int i=0;i<30;i++)
      countervv[isam].push_back(0);
    
  }

  EmuData data;
  Double_t trigeff,rawMet,rawprojvar;
  UInt_t npt20jets;
  const UInt_t kMaxPt20Jets=15;
  TArrayF *btagArray = new TArrayF;
  btagArray->Set(kMaxPt20Jets);
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0; 
  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if((isam==0) && !hasData) continue;
    
    const TString fname = ntupleDir + TString("/") + snamev[isam] + TString("_select.root");
    cout << "Processing " << fname << "..." << endl;   
    infile = new TFile(fname);
    assert(infile); 

    // Get the TTree and set branch address
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
    eventTree->SetBranchAddress("Events",&data);
    // extra branches
    eventTree->SetBranchAddress("npt20jets",&npt20jets);
    eventTree->SetBranchAddress("btagArray",&btagArray);
    eventTree->SetBranchAddress("trigeff",&trigeff);
    eventTree->SetBranchAddress("rawMet",&rawMet);
    eventTree->SetBranchAddress("rawprojvar",&rawprojvar);

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);

      Double_t wgt = 1;
      if(isam!=0)
	wgt = data.weight*lumi;

      // get muon, electron kinematics
      Float_t mupt,mueta,muphi,elept,eleta,elephi;
      if(data.state==kMuEle) {
	mupt  = data.lpt1; mueta = data.leta1; muphi = data.lphi1; elept = data.lpt2; eleta = data.leta2; elephi = data.lphi2;
      } else {
	mupt  = data.lpt2; mueta = data.leta2; muphi = data.lphi2; elept = data.lpt1; eleta = data.leta1; elephi = data.lphi1;
      }

      // fill trigger efficiency plots
      hTrigEffv[isam]        ->Fill(mupt,  elept,   trigeff);
      hTrigEffEntriesv[isam] ->Fill(mupt,  elept);    // keeps track of number of entries in htrigeffv
      hTrigEffvMuPtv[isam]   ->Fill(mupt,  trigeff);
      hTrigEffvElePtv[isam]  ->Fill(elept, trigeff);

      // fill plots after lepton selection
      hMassv[isam]       ->Fill(data.mass,   	wgt);
      hMassVpmissv[isam] ->Fill(data.mass, data.pmet, wgt);
      hMassLv[isam]      ->Fill(data.mass,   	wgt);      
      hMetv[isam]        ->Fill(data.met,    	wgt);
      hMetRawv[isam]     ->Fill(rawMet,      	wgt);
      hProjMetv[isam]    ->Fill(data.pmet,     	wgt);
      hProjVisv[isam]    ->Fill(data.pvis,     	wgt);
      hProjVarv[isam]    ->Fill(0.85*data.pvis - data.pmet, wgt);
      hRawProjVarv[isam] ->Fill(rawprojvar,     wgt);
      hNjetsv[isam]      ->Fill(data.njets,     wgt);
      hNbjetsv[isam]     ->Fill(data.nbjets,    wgt);

      assert(npt20jets<kMaxPt20Jets);
      for(UInt_t ib=0;ib<npt20jets;ib++)
	hBdiscrv[isam]   ->Fill((*btagArray)[ib], wgt);
      
      hDPhiv[isam]   ->Fill(data.dphi*180./pi,wgt);
      hMtv[isam]     ->Fill(data.mt,     wgt);
      hPtv[isam]     ->Fill(data.pt,     wgt);
      hMetDPhiv[isam]->Fill(toolbox::deltaPhi(data.phi,data.metphi)*180./pi,wgt);
      hPt1v[isam]    ->Fill(data.lpt1,   wgt);
      hEta1v[isam]   ->Fill(data.leta1,  wgt);
      hPhi1v[isam]   ->Fill(data.lphi1,  wgt);
      hPt2v[isam]    ->Fill(data.lpt2,   wgt);
      hEta2v[isam]   ->Fill(data.leta2,  wgt);
      hPhi2v[isam]   ->Fill(data.lphi2,  wgt);
      hPtMuv[isam]   ->Fill(mupt,        wgt);
      hEtaMuv[isam]  ->Fill(mueta,       wgt);
      hPhiMuv[isam]  ->Fill(muphi,       wgt);
      hPtElev[isam]  ->Fill(elept,       wgt);
      hEtaElev[isam] ->Fill(eleta,       wgt);
      hPhiElev[isam] ->Fill(elephi,      wgt);
      if(data.jpt1>0) {
	hJetPt1v[isam]	->Fill(data.jpt1,    wgt);
	hJetEta1v[isam]	->Fill(data.jeta1,   wgt);
      }
      if(data.jpt2>0) {
	hJetPt2v[isam]	->Fill(data.jpt2,    wgt);
	hJetEta2v[isam]	->Fill(data.jeta2,   wgt);
      }
      if(data.mjj>0) {
	hJetDPhiv[isam]	->Fill(toolbox::deltaPhi(data.jphi1,data.jphi2)*180./pi,   wgt);
	hMjjv[isam]	->Fill(data.mjj,     wgt);
	hDEtav[isam]	->Fill(fabs(data.jeta1 - data.jeta2),   wgt);
	hEtaProdv[isam]	->Fill(data.jeta1*data.jeta2,           wgt);
      }
      if(data.bjpt>0) {
	hBJetPtv[isam]	->Fill(data.bjpt,    wgt);
	hBJetEtav[isam]	->Fill(data.bjeta,   wgt);
	hBJetPhiv[isam]	->Fill(data.bjphi,   wgt);
      }
      hNPVv[isam]    ->Fill(data.nPV,    wgt);

      nSelv[isam]    += wgt;
      nSelVarv[isam] += wgt*wgt;

      if(0.85*data.pvis - data.pmet > 25) continue;

      // inclusive
      hMass_iv[isam]    ->Fill(data.mass,   wgt);
      hMassL_iv[isam]   ->Fill(data.mass,   wgt);

      nSel_iv[isam]     += wgt;
      nSelVar_iv[isam]  += wgt*wgt;

      // no b-tag
      if(data.njets<=1 && data.nbjets==0) {
	hMass_nobv[isam]   ->Fill(data.mass,   (isam==0) ? wgt : kCat_nob*wgt);      
	hMassL_nobv[isam]  ->Fill(data.mass,   (isam==0) ? wgt : kCat_nob*wgt);      

	nSel_nobv[isam]     += (isam==0) ? wgt : kCat_nob*wgt;
	nSelVar_nobv[isam]  += (isam==0) ? wgt : kCat_nob*wgt*kCat_nob*wgt;
      }

      // counters for vbf cuts
      countervv[isam][0]++;
      if(data.njets  ==    2        ) 			countervv[isam][1]++;
      if(data.nbjets ==    0        ) 		        countervv[isam][2]++;  
      if(data.mjj     >   mjjMin    ) 			countervv[isam][3]++;
      if(data.jeta1*data.jeta2       <  0         ) 	countervv[isam][4]++;
      if(fabs(data.jeta1-data.jeta2) >  dEtaMin   ) 	countervv[isam][5]++;

      if(data.njets == 0)  countervv[isam][6]++;
      if(data.njets == 1)  countervv[isam][7]++;
      if(data.njets == 2)  countervv[isam][8]++;
      if(data.njets >  2)  countervv[isam][9]++;
	 
      // vbf
      if(data.njets    ==    2        		&&
	 data.nbjets   ==    0        		&&
	 data.mjj       >   mjjMin    		&&
	 data.jeta1*data.jeta2       <  0       &&
	 fabs(data.jeta1-data.jeta2) >  dEtaMin     ) {

	hMass_vbfv[isam]   ->Fill(data.mass,   (isam==0) ? wgt : kCat_vbf*wgt);      
	hMassL_vbfv[isam]  ->Fill(data.mass,   (isam==0) ? wgt : kCat_vbf*wgt);      

	nSel_vbfv[isam]     += (isam==0) ? wgt : kCat_vbf*wgt;
	nSelVar_vbfv[isam]  += (isam==0) ? wgt : kCat_vbf*wgt*kCat_vbf*wgt;
      }

      // b-tag
      if(data.njets<=1 && data.nbjets>0) {
	hMass_bv[isam]   ->Fill(data.mass,   (isam==0) ? wgt : kCat_b*wgt);      
	hMassL_bv[isam]  ->Fill(data.mass,   (isam==0) ? wgt : kCat_b*wgt);      

	nSel_bv[isam]     += (isam==0) ? wgt : kCat_b*wgt;
	nSelVar_bv[isam]  += (isam==0) ? wgt : kCat_b*wgt*kCat_b*wgt;
      }

    }
    delete infile;
    infile=0, eventTree=0;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  TFile *fval=0;
  if(CPlot::sOutDir.Contains("limit"))
    fval = TFile::Open("data/HttSelectionPlotsFakePrediction-limit-bins.root");
  else
    fval = TFile::Open("data/HttSelectionPlotsFakePrediction.root");
  assert(fval);
  TH1D *fNVertices 		= (TH1D*)fval->Get("hNVertices"		);           fNVertices 	->SetDirectory(0);
  TH1D *fptm 			= (TH1D*)fval->Get("hEMUptm"		);           fptm 		->SetDirectory(0);
  TH1D *fpte 			= (TH1D*)fval->Get("hEMUpte"		);           fpte 		->SetDirectory(0);
  TH1D *fetam 			= (TH1D*)fval->Get("hEMUetam"		);	     fetam 		->SetDirectory(0);
  TH1D *fetae 			= (TH1D*)fval->Get("hEMUetae"		);	     fetae 		->SetDirectory(0);
  TH1D *fdphi 			= (TH1D*)fval->Get("hEMUdphi"		);	     fdphi 		->SetDirectory(0);
  TH1D *fprojVis 		= (TH1D*)fval->Get("hEMUprojVis"	);	     fprojVis 		->SetDirectory(0);
  TH1D *fprojMet 		= (TH1D*)fval->Get("hEMUprojMet"	);	     fprojMet 		->SetDirectory(0);
  TH1D *fpzetaVar 		= (TH1D*)fval->Get("hEMUpzetaVar"	);	     fpzetaVar 		->SetDirectory(0);
  TH1D *fnjetspt30 		= (TH1D*)fval->Get("hEMUnjetspt30"	);	     fnjetspt30 	->SetDirectory(0);
  TH1D *fnbtagjets 		= (TH1D*)fval->Get("hEMUnbtagjets"	);	     fnbtagjets 	->SetDirectory(0);
  TH1D *fbtag 			= (TH1D*)fval->Get("hEMUbtag"		);	     fbtag 		->SetDirectory(0);
  TH1D *fmet 			= (TH1D*)fval->Get("hEMUmet"		);           fmet 		->SetDirectory(0);
  TH1D *fvisMassCut0 		= (TH1D*)fval->Get("hHTTvisMassCut0"	);	     fvisMassCut0 	->SetDirectory(0);
  TH1D *fvisMassCut4 		= (TH1D*)fval->Get("hHTTvisMassCut4"	);	     fvisMassCut4 	->SetDirectory(0);
  TH1D *fvisMassCut5 		= (TH1D*)fval->Get("hHTTvisMassCut5"	);	     fvisMassCut5 	->SetDirectory(0);
  TH1D *fvisMassCut6 		= (TH1D*)fval->Get("hHTTvisMassCut6"	);	     fvisMassCut6 	->SetDirectory(0);
  TH1D *fvisMassCut7 		= (TH1D*)fval->Get("hHTTvisMassCut7"	);	     fvisMassCut7 	->SetDirectory(0);
  TH1D *fVBFDijetMass 		= (TH1D*)fval->Get("hVBFDijetMass"	);	     fVBFDijetMass 	->SetDirectory(0);
  TH1D *fVBFDijetDeltaEta 	= (TH1D*)fval->Get("hVBFDijetDeltaEta"	);	     fVBFDijetDeltaEta	->SetDirectory(0);
  TH1D *fnbtagjetsvbf 		= (TH1D*)fval->Get("hnbtagjetsvbf"	);	     fnbtagjetsvbf 	->SetDirectory(0);
  fval->Close(); delete fval;

  TCanvas *c = MakeCanvas("c","c",800,600);
  
  // string buffers
  char ylabel[100];   // y-axis label
  char lumitext[50];
  if(lumi>0) {
    if(lumi<100) { sprintf(lumitext,"#int#font[12]{L}dt = %.3g pb^{-1}",lumi); }
    else         { sprintf(lumitext,"#int#font[12]{L}dt = %.3g fb^{-1}",lumi/1000.); }
  }
  
  //
  // Overall plots
  //

  CPlot plotEffx("effx","","mu pt [GeV]","ele pt [GeV]");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    // calculate average efficiency
    for(Int_t xbin=1;xbin<hTrigEffv[isam]->GetNbinsX()+1;xbin++) {
      for(Int_t ybin=1;ybin<hTrigEffv[isam]->GetNbinsY()+1;ybin++) {
	Double_t tot      = hTrigEffv[isam]->GetBinContent(xbin,ybin);
	Double_t nentries = hTrigEffEntriesv[isam]->GetBinContent(xbin,ybin);
	hTrigEffv[isam]->SetBinContent(xbin, ybin, (nentries>0) ? tot/nentries : 0);
      }
    }
    plotEffx.AddHist2D(hTrigEffv[isam],"surf",kWhite,kBlue);
  }
  plotEffx.Draw(c,kTRUE,"C");

  CPlot plotEffvMuPtx("effvmupt","","mu pt [GeV]","eff.");
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotEffvMuPtx.AddHist2D(hTrigEffvMuPtv[isam],"surf",kWhite,kBlue);
  plotEffvMuPtx.Draw(c,kTRUE,"C");

  CPlot plotEffvElePtx("effvelept","","ele pt [GeV]","eff.");
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotEffvElePtx.AddHist2D(hTrigEffvElePtv[isam],"surf",kWhite,kBlue);
  plotEffvElePtx.Draw(c,kTRUE,"C");

  //
  // after lepton selection:
  //

  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMassv[0]->GetBinWidth(1));
  CPlot plotMass("vismass_lept","","m_{vis} [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); }
  plotMass.AddToStack((TH1F*)fvisMassCut0,"fakes",860); //cout << "fvisMassCut0: " << fvisMassCut0->GetNbinsX() << endl;
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color);
  plotMass.TransLegend(0.1,0);
  plotMass.SetXRange(0,200);
  if(lumi>0) plotMass.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotMass.AddTextBox("ele + mu",0.44,0.8,0.64,0.7,0);
  plotMass.Draw(c,kTRUE,format);

  CPlot plotMassVpmiss("massvpmiss","","m_{vis} [GeV]","#slash{p}_{#zeta} [GeV]");
//   for(UInt_t isam=1; isam<samplev.size(); isam++)
  assert(hMassVpmissv[0]);
  plotMassVpmiss.AddHist2D(hMassVpmissv[0],"surf",kWhite,kBlue);
  plotMassVpmiss.Draw(c,kTRUE,"C");

  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMetv[0]->GetBinWidth(1));
  CPlot plotMet("met","","met [GeV/c^{2}]",ylabel);
  if(hasData) { plotMet.AddHist1D(hMetv[0],samplev[0]->label,"E"); }
  plotMet.AddToStack((TH1F*)fmet,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMet.AddToStack(hMetv[isam],samplev[isam]->label,samplev[isam]->color);
  if(samplev.size()>5)
    plotMet.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMet.TransLegend(0.1,0);
  if(lumi>0) plotMet.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotMet.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMetRawv[0]->GetBinWidth(1));
  CPlot plotMetRaw("metraw","","raw met [GeV/c^{2}]",ylabel);
  if(hasData) { plotMetRaw.AddHist1D(hMetRawv[0],samplev[0]->label,"E"); }
  plotMetRaw.AddToStack((TH1F*)fmet,"fakes",860); //cout << "fmet: " << fmet->GetNbinsX() << " " << fmet->GetXaxis()->GetXmin() << " " << fmet->GetXaxis()->GetXmax() << endl;
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMetRaw.AddToStack(hMetRawv[isam],samplev[isam]->label,samplev[isam]->color);
  if(samplev.size()>5)
    plotMetRaw.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMetRaw.TransLegend(0.1,0);
  if(lumi>0) plotMetRaw.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotMetRaw.Draw(c,kTRUE,format);

  // projection variables
  sprintf(ylabel,"Events / %.2f",hProjMetv[0]->GetBinWidth(1));
  CPlot plotProjMet("pzetamiss","","#slash{p}_{#zeta} [GeV]",ylabel);
  if(hasData) { plotProjMet.AddHist1D(hProjMetv[0],samplev[0]->label,"E"); }
  plotProjMet.AddToStack((TH1F*)fprojMet,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotProjMet.AddToStack(hProjMetv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotProjMet.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotProjMet.TransLegend(0.15,0);
  if(samplev.size()>5)
    plotProjMet.SetLegend(0.55,0.55,0.88,0.9);
  plotProjMet.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hProjVisv[0]->GetBinWidth(1));
  CPlot plotProjVis("pzetavis","","p_{#zeta}^{vis} [GeV]",ylabel);
  if(hasData) { plotProjVis.AddHist1D(hProjVisv[0],samplev[0]->label,"E"); }
  plotProjVis.AddToStack((TH1F*)fprojVis,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotProjVis.AddToStack(hProjVisv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotProjVis.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotProjVis.SetLegend(0.55,0.55,0.88,0.9);
//   plotProjVis.SetLogy();
  plotProjVis.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hProjVarv[0]->GetBinWidth(1));
  CPlot plotProjVar("pzetavar","","0.85*p_{#zeta}^{vis} - #slash{p}_{#zeta} [GeV]",ylabel);
  if(hasData) { plotProjVar.AddHist1D(hProjVarv[0],samplev[0]->label,"E"); }
  plotProjVar.AddToStack((TH1F*)fpzetaVar,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotProjVar.AddToStack(hProjVarv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotProjVar.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotProjVar.SetLegend(0.55,0.55,0.88,0.9);
  plotProjVar.Draw(c,kTRUE,format);

  // plot corrected stack on top of the uncorrected stack:
  TH1F *hTmpStack = 0;
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==1) {
      hTmpStack = new TH1F(*hProjVarv[isam]);
      hTmpStack->Reset();
    }
    hTmpStack->Add(hProjVarv[isam]);
  }
  hTmpStack->Add(fpzetaVar);
  
  sprintf(ylabel,"Events / %.2f",hRawProjVarv[0]->GetBinWidth(1));
  CPlot plotRawProjVar("rawProjVar","","raw 0.85*p_{#zeta}^{vis} - #slash{p}_{#zeta} [GeV]",ylabel);
  if(hasData) { plotRawProjVar.AddHist1D(hRawProjVarv[0],samplev[0]->label,"E"); }
  plotRawProjVar.AddHist1D(hTmpStack,"corrected","hist",kBlue);
  plotRawProjVar.AddToStack((TH1F*)fpzetaVar,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotRawProjVar.AddToStack(hRawProjVarv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotRawProjVar.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotRawProjVar.SetLegend(0.55,0.55,0.88,0.9);
  plotRawProjVar.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hNjetsv[0]->GetBinWidth(1));
  CPlot plotNjets("njets","","number of jets",ylabel);
  if(hasData) { plotNjets.AddHist1D(hNjetsv[0],samplev[0]->label,"E"); }
  plotNjets.AddToStack((TH1F*)fnjetspt30,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotNjets.AddToStack(hNjetsv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotNjets.AddTextBox(lumitext,0.75,0.35,0.95,0.3,0);
  plotNjets.TransLegend(0.13,0);
  plotNjets.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hNbjetsv[0]->GetBinWidth(1));
  CPlot plotNbjets("nbjets","","number of b-jets",ylabel);
  if(hasData) { plotNbjets.AddHist1D(hNbjetsv[0],samplev[0]->label,"E"); }
  plotNbjets.AddToStack((TH1F*)fnbtagjets,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotNbjets.AddToStack(hNbjetsv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotNbjets.AddTextBox(lumitext,0.75,0.35,0.95,0.3,0);
  if(samplev.size()>5)
    plotNbjets.SetLegend(0.55,0.55,0.88,0.9);
  plotNbjets.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hBdiscrv[0]->GetBinWidth(1));
  CPlot plotBdiscr("btag","","b-tag discr.",ylabel);
  if(hasData) { plotBdiscr.AddHist1D(hBdiscrv[0],samplev[0]->label,"E"); }
  plotBdiscr.AddToStack((TH1F*)fbtag,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotBdiscr.AddToStack(hBdiscrv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotBdiscr.AddTextBox(lumitext,0.75,0.35,0.95,0.3,0);
  if(samplev.size()>5)
    plotBdiscr.SetLegend(0.55,0.55,0.88,0.9);
  plotBdiscr.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hDPhiv[0]->GetBinWidth(1));
  CPlot plotDPhi("dphi_emu","","#Delta^{}#phi_{emu} [deg]",ylabel);
  if(hasData) { plotDPhi.AddHist1D(hDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotDPhi.AddToStack(hDPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotDPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotDPhi.AddToStack((TH1F*)fdphi,"fakes",860);
  plotDPhi.TransLegend(-0.15,0);
  assert(plotDPhi.GetStack());
  plotDPhi.SetYRange(0,1.5*(plotDPhi.GetStack()->GetMaximum()));
  plotDPhi.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f GeV/c",hMtv[0]->GetBinWidth(1));
  CPlot plotMt("mt","","m_{T}(ll,#slash{E}_{T}) [GeV/c^{2}]",ylabel);
  if(hasData) { plotMt.AddHist1D(hMtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMt.AddToStack(hMtv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotMt.AddTextBox(lumitext,0.27,0.85,0.47,0.8,0);
  if(samplev.size()>5)
    plotMt.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMt.TransLegend(-.1,0);
  plotMt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f GeV/c",hPtv[0]->GetBinWidth(1));
  CPlot plotPt("pt","","p_{T}^{ll} [GeV/c]",ylabel);
  if(hasData) { plotPt.AddHist1D(hPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPt.AddToStack(hPtv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPt.AddTextBox(lumitext,0.37,0.85,0.57,0.8,0);
  plotPt.TransLegend(0.1,0);
  plotPt.Draw(c,kTRUE,format);

  // met 
  sprintf(ylabel,"Events / %.2f",hMetDPhiv[0]->GetBinWidth(1));
  CPlot plotMetDPhi("metdphi","","#Delta^{}#phi(ll,#slash{E}_{T}) [deg]",ylabel);
  if(hasData) { plotMetDPhi.AddHist1D(hMetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMetDPhi.AddToStack(hMetDPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotMetDPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMetDPhi.TransLegend(-0.1,0.01);
  plotMetDPhi.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.1f GeV/c",hPt1v[0]->GetBinWidth(1));
  CPlot plotPt1("pt1","","leading lepton p_{T} [GeV/c]",ylabel);
  if(hasData) { plotPt1.AddHist1D(hPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPt1.AddToStack(hPt1v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPt1.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  if(samplev.size()>5)
    plotPt1.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPt1.TransLegend(0.1,0);
  plotPt1.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.2f",hEta1v[0]->GetBinWidth(1));
  CPlot plotEta1("eta1","","leading lepton #eta",ylabel);
  if(hasData) { plotEta1.AddHist1D(hEta1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotEta1.AddToStack(hEta1v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotEta1.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotEta1.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotEta1.TransLegend(0.1,0);
//   plotEta1.SetYRange(0,1.5*(plotEta1.GetStack()->GetMaximum()));
  plotEta1.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.2f",hPhi1v[0]->GetBinWidth(1));
  CPlot plotPhi1("phi1","","leading lepton #phi",ylabel);
  if(hasData) { plotPhi1.AddHist1D(hPhi1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPhi1.AddToStack(hPhi1v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPhi1.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotPhi1.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPhi1.TransLegend(0.1,0);
  plotPhi1.SetYRange(0,2.0*(plotPhi1.GetStack()->GetMaximum()));
  plotPhi1.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f GeV/c",hPt2v[0]->GetBinWidth(1));
  CPlot plotPt2("pt2","","trailing lepton p_{T} [GeV/c]",ylabel);
  if(hasData) { plotPt2.AddHist1D(hPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPt2.AddToStack(hPt2v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPt2.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  if(samplev.size()>5)
    plotPt2.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPt2.TransLegend(0.1,0);
  plotPt2.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.2f",hEta2v[0]->GetBinWidth(1));
  CPlot plotEta2("eta2","","trailing lepton #eta",ylabel);
  if(hasData) { plotEta2.AddHist1D(hEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotEta2.AddToStack(hEta2v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotEta2.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotEta2.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotEta2.TransLegend(0.1,0);
//   plotEta2.SetYRange(0,1.5*(plotEta2.GetStack()->GetMaximum()));
  plotEta2.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.2f",hPhi2v[0]->GetBinWidth(1));
  CPlot plotPhi2("phi2","","trailing lepton #phi",ylabel);
  if(hasData) { plotPhi2.AddHist1D(hPhi2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPhi2.AddToStack(hPhi2v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPhi2.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotPhi2.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPhi2.TransLegend(0.1,0);
  plotPhi2.SetYRange(0,2.0*(plotPhi2.GetStack()->GetMaximum()));
  plotPhi2.Draw(c,kTRUE,format);

  // mu / electron kinematics
  sprintf(ylabel,"Events / %.1f GeV/c",hPtMuv[0]->GetBinWidth(1));
  CPlot plotPtMu("pt_mu","","mu p_{T} [GeV/c]",ylabel);
  if(hasData) { plotPtMu.AddHist1D(hPtMuv[0],samplev[0]->label,"E"); }
  plotPtMu.AddToStack((TH1F*)fptm,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPtMu.AddToStack(hPtMuv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPtMu.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  if(samplev.size()>5)
    plotPtMu.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPtMu.TransLegend(0.1,0);
  plotPtMu.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.2f",hEtaMuv[0]->GetBinWidth(1));
  CPlot plotEtaMu("eta_mu","","mu #eta",ylabel);
  if(hasData) { plotEtaMu.AddHist1D(hEtaMuv[0],samplev[0]->label,"E"); }
  plotEtaMu.AddToStack((TH1F*)fetam,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotEtaMu.AddToStack(hEtaMuv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotEtaMu.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotEtaMu.SetLegend(0.78,0.55,0.99,0.9);
//   plotEtaMu.SetYRange(0,1.5*(plotEtaMu.GetStack()->GetMaximum()));
  plotEtaMu.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.2f",hPhiMuv[0]->GetBinWidth(1));
  CPlot plotPhiMu("phimu","","mu #phi",ylabel);
  if(hasData) { plotPhiMu.AddHist1D(hPhiMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPhiMu.AddToStack(hPhiMuv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPhiMu.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotPhiMu.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPhiMu.TransLegend(0.1,0);
  plotPhiMu.SetYRange(0,2.0*(plotPhiMu.GetStack()->GetMaximum()));
  plotPhiMu.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.1f GeV/c",hPtElev[0]->GetBinWidth(1));
  CPlot plotPtEle("pt_e","","ele p_{T} [GeV/c]",ylabel);
  if(hasData) { plotPtEle.AddHist1D(hPtElev[0],samplev[0]->label,"E"); }
  plotPtEle.AddToStack((TH1F*)fpte,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPtEle.AddToStack(hPtElev[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPtEle.AddTextBox(lumitext,0.45,0.85,0.65,0.8,0);
  if(samplev.size()>5)
    plotPtEle.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPtEle.TransLegend(0.1,0);
  plotPtEle.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.2f",hEtaElev[0]->GetBinWidth(1));
  CPlot plotEtaEle("eta_e","","ele #eta",ylabel);
  if(hasData) { plotEtaEle.AddHist1D(hEtaElev[0],samplev[0]->label,"E"); }
  plotEtaEle.AddToStack((TH1F*)fetae,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotEtaEle.AddToStack(hEtaElev[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotEtaEle.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotEtaEle.SetLegend(0.78,0.57,0.99,0.91);
//   plotEtaEle.SetYRange(0,1.5*(plotEtaEle.GetStack()->GetMaximum()));
  plotEtaEle.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events / %.2f",hPhiElev[0]->GetBinWidth(1));
  CPlot plotPhiEle("phiele","","ele #phi",ylabel);
  if(hasData) { plotPhiEle.AddHist1D(hPhiElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotPhiEle.AddToStack(hPhiElev[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotPhiEle.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotPhiEle.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotPhiEle.TransLegend(0.1,0);
  plotPhiEle.SetYRange(0,2.0*(plotPhiEle.GetStack()->GetMaximum()));
  plotPhiEle.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hJetPt1v[0]->GetBinWidth(1));
  CPlot plotJetPt1("jetpt1","","leading jet pt [GeV]",ylabel);
  if(hasData) { plotJetPt1.AddHist1D(hJetPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotJetPt1.AddToStack(hJetPt1v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotJetPt1.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotJetPt1.TransLegend(0.1,0);
  plotJetPt1.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hJetPt2v[0]->GetBinWidth(1));
  CPlot plotJetPt2("jetpt2","","second jet pt [GeV]",ylabel);
  if(hasData) { plotJetPt2.AddHist1D(hJetPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotJetPt2.AddToStack(hJetPt2v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotJetPt2.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  if(samplev.size()>5)
    plotJetPt2.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotJetPt2.TransLegend(0.1,0);
  plotJetPt2.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hJetEta1v[0]->GetBinWidth(1));
  CPlot plotJetEta1("jeteta1","","leading jet #eta",ylabel);
  if(hasData) { plotJetEta1.AddHist1D(hJetEta1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotJetEta1.AddToStack(hJetEta1v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotJetEta1.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotJetEta1.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotJetEta1.TransLegend(0.1,0);
  plotJetEta1.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hJetEta2v[0]->GetBinWidth(1));
  CPlot plotJetEta2("jeteta2","","second jet #eta",ylabel);
  if(hasData) { plotJetEta2.AddHist1D(hJetEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotJetEta2.AddToStack(hJetEta2v[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotJetEta2.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotJetEta2.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotJetEta2.TransLegend(0.1,0);
  plotJetEta2.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hJetDPhiv[0]->GetBinWidth(1));
  CPlot plotJetDPhi("jetdphi","","#Delta#phi (jj)",ylabel);
  if(hasData) { plotJetDPhi.AddHist1D(hJetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotJetDPhi.AddToStack(hJetDPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotJetDPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotJetDPhi.TransLegend(-0.3,-0.14);
  plotJetDPhi.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hMjjv[0]->GetBinWidth(1));
  CPlot plotMjj("mjj","","dijet mass",ylabel);
  if(hasData) { plotMjj.AddHist1D(hMjjv[0],samplev[0]->label,"E"); }
  plotMjj.AddToStack((TH1F*)fVBFDijetMass,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMjj.AddToStack(hMjjv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotMjj.AddTextBox(lumitext,0.37,0.85,0.57,0.8,0);
  plotMjj.TransLegend(0.1,0);
  plotMjj.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hDEtav[0]->GetBinWidth(1));
  CPlot plotDEta("jetdeta","","#Delta#eta (jj)",ylabel);
  if(hasData) { plotDEta.AddHist1D(hDEtav[0],samplev[0]->label,"E"); }
  plotDEta.AddToStack((TH1F*)fVBFDijetDeltaEta,"fakes",860); cout << "deta: " << fVBFDijetDeltaEta->Integral() << endl;
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotDEta.AddToStack(hDEtav[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotDEta.AddTextBox(lumitext,0.33,0.85,0.57,0.8,0);
  plotDEta.TransLegend(0.1,0);
  plotDEta.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hEtaProdv[0]->GetBinWidth(1));
  CPlot plotEtaProd("jetEtaProd","","#eta1*#eta2(jj)",ylabel);
  if(hasData) { plotEtaProd.AddHist1D(hEtaProdv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotEtaProd.AddToStack(hEtaProdv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotEtaProd.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotEtaProd.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotEtaProd.TransLegend(0.1,0);
  plotEtaProd.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hBJetPtv[0]->GetBinWidth(1));
  CPlot plotBJetPt("bjetpt","","lead b-jet pt [GeV]",ylabel);
  if(hasData) { plotBJetPt.AddHist1D(hBJetPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotBJetPt.AddToStack(hBJetPtv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotBJetPt.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  if(samplev.size()>5)
    plotBJetPt.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotBJetPt.TransLegend(0.1,0);
  plotBJetPt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hBJetEtav[0]->GetBinWidth(1));
  CPlot plotBJetEta("bjeteta","","lead b-jet #eta",ylabel);
  if(hasData) { plotBJetEta.AddHist1D(hBJetEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotBJetEta.AddToStack(hBJetEtav[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotBJetEta.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  if(samplev.size()>5)
    plotBJetEta.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotBJetEta.TransLegend(0.1,0);
  plotBJetEta.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events / %.2f",hBJetPhiv[0]->GetBinWidth(1));
  CPlot plotBJetPhi("bjetphi","","lead b-jet #phi",ylabel);
  if(hasData) { plotBJetPhi.AddHist1D(hBJetPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotBJetPhi.AddToStack(hBJetPhiv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotBJetPhi.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotBJetPhi.TransLegend(-0.3,-0.14);
  plotBJetPhi.Draw(c,kTRUE,format);

  CPlot plotNPV("nvertices_reweighted","","N_{PV}","Events");
  if(hasData) { plotNPV.AddHist1D(hNPVv[0],samplev[0]->label,"E"); }
  plotNPV.AddToStack((TH1F*)fNVertices,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotNPV.AddToStack(hNPVv[isam],samplev[isam]->label,samplev[isam]->color);
  if(lumi>0) plotNPV.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotNPV.TransLegend(0.08,0);
  plotNPV.Draw(c,kTRUE,format);

// uncomment to make the hists for vertex reweighting
//   TFile fnvtx("data/nvtxhists.root","recreate");
//   for(UInt_t isam=1; isam<samplev.size(); isam++) {
//     if(snamev[isam] == "ztt") {
//       hNPVv[isam]->SetName("npv_ztt");
//       hNPVv[isam]->Write();
//     }
//   }
//   if(hasData) {
//     hNPVv[0]->SetName("npv_data");
//     hNPVv[0]->Write();
//   }
//   fnvtx.Close();

  //----------------------------------------------------------------------------------------
  // inclusive
  //----------------------------------------------------------------------------------------  
  
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMass_iv[0]->GetBinWidth(1));
  CPlot plotMass_i("vismass_incl","","m_{vis} [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass_i.AddHist1D(hMass_iv[0],samplev[0]->label,"E"); }
  plotMass_i.AddToStack((TH1F*)fvisMassCut5,"fakes",860); //cout << "fvisMassCut5: " <<  fvisMassCut5->GetNbinsX() << endl;
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMass_i.AddToStack(hMass_iv[isam],samplev[isam]->label,samplev[isam]->color);
  plotMass_i.SetLegend(0.75,0.55,0.98,0.9);
  if(lumi>0) plotMass_i.AddTextBox(lumitext,0.43,0.85,0.63,0.8,0);
  plotMass_i.AddTextBox("inclusive",0.44,0.8,0.64,0.7,0);
  plotMass_i.SetXRange(0,200);
  plotMass_i.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // no b-tag
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMass_nobv[0]->GetBinWidth(1));
  CPlot plotMass_nob("vismass_class1","","m_{vis} [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass_nob.AddHist1D(hMass_nobv[0],samplev[0]->label,"E"); }
  plotMass_nob.AddToStack((TH1F*)fvisMassCut6,"fakes",860); //cout << "fvisMassCut6: " <<  fvisMassCut6->GetNbinsX() << endl;
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMass_nob.AddToStack(hMass_nobv[isam],samplev[isam]->label,samplev[isam]->color);
  plotMass_nob.TransLegend(0.13,0);
  if(lumi>0) plotMass_nob.AddTextBox(lumitext,0.43,0.85,0.63,0.8,0);
  plotMass_nob.AddTextBox("no b-tag",0.44,0.8,0.64,0.7,0);
  plotMass_nob.SetXRange(0,200);
  plotMass_nob.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // vbf
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMass_vbfv[0]->GetBinWidth(1));
  CPlot plotMass_vbf("vismass_class2a","","m_{vis} [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass_vbf.AddHist1D(hMass_vbfv[0],samplev[0]->label,"E"); }
  plotMass_vbf.AddToStack((TH1F*)fvisMassCut4,"fakes",860); //cout << "fvisMassCut4: " <<  fvisMassCut4->GetNbinsX() << endl;
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMass_vbf.AddToStack(hMass_vbfv[isam],samplev[isam]->label,samplev[isam]->color);
  plotMass_vbf.SetLegend(0.65,0.55,0.85,0.9);
  if(lumi>0) plotMass_vbf.AddTextBox(lumitext,0.21,0.85,0.41,0.8,0);
  plotMass_vbf.AddTextBox("vbf",0.21,0.8,0.41,0.7,0);
  plotMass_vbf.SetXRange(0,200);
  plotMass_vbf.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // b-tag
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"Events / %.1f GeV/c^{2}",hMass_bv[0]->GetBinWidth(1));
  CPlot plotMass_b("vismass_class2b","","m_{vis} [GeV/c^{2}]",ylabel);
  if(hasData) { plotMass_b.AddHist1D(hMass_bv[0],samplev[0]->label,"E"); }
  plotMass_b.AddToStack((TH1F*)fvisMassCut7,"fakes",860);
  for(UInt_t isam=1; isam<samplev.size(); isam++)
    plotMass_b.AddToStack(hMass_bv[isam],samplev[isam]->label,samplev[isam]->color);
  if(samplev.size()>5)
    plotMass_b.SetLegend(0.75,0.55,0.98,0.9);
  else
    plotMass_b.TransLegend(0.1,0);
  if(lumi>0) plotMass_b.AddTextBox(lumitext,0.41,0.85,0.61,0.8,0);
  plotMass_b.AddTextBox("b-tagged",0.41,0.8,0.61,0.7,0);
  plotMass_b.SetXRange(0,200);
  plotMass_b.Draw(c,kTRUE,format);

  // write hists to limit file 
  TFile flimit(CPlot::sOutDir + "/limit-inputs.root","recreate");
  flimit.mkdir("emu_1");
  flimit.mkdir("emu_2a");
  flimit.mkdir("emu_2b");
  flimit.mkdir("emu_X");

  TString histname;
  for(UInt_t isam=0;isam<samplev.size();isam++) {
    CSample samp = (*samplev[isam]);
//     if     (snamev[isam].Contains("fakes", TString::kIgnoreCase))      	histname = "Fakes"; 
//     else
    if(snamev[isam].Contains("ewk",   TString::kIgnoreCase))   	        histname = "EWK"; 
    else if(snamev[isam].Contains("ttbar", TString::kIgnoreCase)) 	histname = "ttbar";
    else if(snamev[isam].Contains("Ztt",   TString::kIgnoreCase))   	histname = "Ztt"; 
    else if(snamev[isam].Contains("data",  TString::kIgnoreCase))  	histname = "data_obs";
    else if( (samp.fnamev[0].Contains("gg")) ||
	     (samp.fnamev[0].Contains("bb")) ||
	     (samp.fnamev[0].Contains("vbf"))  )                        histname = "Higgs_" + snamev[isam];
    

    flimit.cd("emu_1");
    hMassL_nobv[isam]->SetName(histname);
    hMassL_nobv[isam]->Write();
    flimit.cd("emu_2a");
    hMassL_vbfv[isam]->SetName(histname);
    hMassL_vbfv[isam]->Write();
    flimit.cd("emu_2b");
    hMassL_bv[isam]->SetName(histname);
    hMassL_bv[isam]->Write();
    flimit.cd("emu_X");
    hMassL_iv[isam]->SetName(histname);
    hMassL_iv[isam]->Write();
  }

  histname = "Fakes";
  cout << "ele + mu: " << fvisMassCut0->Integral() << endl;
  flimit.cd("emu_1");
  fvisMassCut6->SetName(histname); cout << "emu_1: " << fvisMassCut6->Integral() << endl;
  fvisMassCut6->Write();
  flimit.cd("emu_2a");
  fvisMassCut4->SetName(histname); cout << "emu_2a: " << fvisMassCut4->Integral() << endl;
  fvisMassCut4->Write();
  flimit.cd("emu_2b");
  fvisMassCut7->SetName(histname); cout << "emu_2b: " << fvisMassCut7->Integral() << endl;
  fvisMassCut7->Write();
  flimit.cd("emu_X");
  fvisMassCut5->SetName(histname); cout << "emu_X: " << fvisMassCut5->Integral() << endl;
  fvisMassCut5->Write();

  flimit.Close();
  
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  Double_t nTotal=0;
  Double_t nTotalVar=0;
  Double_t nTotal_i=0;
  Double_t nTotalVar_i=0;
  Double_t nTotal_nob=0;
  Double_t nTotalVar_nob=0;
  Double_t nTotal_vbf=0;
  Double_t nTotalVar_vbf=0;
  Double_t nTotal_b=0;
  Double_t nTotalVar_b=0;
  
  cout << setw(28) << "lepton sele." <<setw(20) << "inclusive" << setw(25) << "no btag" << setw(19)
       << "vbf" << setw(17) << "btag" << endl;


  cout << setw(10) << "fakes";
  cout << setw(10) << setprecision(3) << fixed << fvisMassCut0->Integral() << " +/- " << setw(6) << setprecision(3) << fixed << -1;
  nTotal +=fvisMassCut0->Integral();
  cout << setw(10) << setprecision(3) << fixed << fvisMassCut5->Integral() << " +/- " << setw(6) << setprecision(3) << fixed << -1;
  nTotal_i +=fvisMassCut5->Integral();
  cout << setw(10) << setprecision(3) << fixed << fvisMassCut6->Integral() << " +/- " << setw(6) << setprecision(3) << fixed << -1;
  nTotal_nob +=fvisMassCut6->Integral();
  cout << setw(10) << setprecision(3) << fixed << fvisMassCut4->Integral() << " +/- " << setw(6) << setprecision(3) << fixed << -1;
  nTotal_vbf +=fvisMassCut4->Integral();
  cout << setw(10) << setprecision(3) << fixed << fvisMassCut7->Integral() << " +/- " << setw(6) << setprecision(3) << fixed << -1;
  nTotal_b +=fvisMassCut7->Integral();
  cout << endl;

  if(samplev.size()>1) {
    for(UInt_t isam=1; isam<samplev.size(); isam++) {      
      cout << setw(10) << snamev[isam];
      cout << setw(10) << setprecision(3) << fixed << nSelv[isam]        << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[isam]);
      cout << setw(10) << setprecision(3) << fixed << nSel_iv[isam]      << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[isam]);
      cout << setw(10) << setprecision(3) << fixed << nSel_nobv[isam]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_nobv[isam]);
      cout << setw(10) << setprecision(3) << fixed << nSel_vbfv[isam]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[isam]);
      cout << setw(10) << setprecision(3) << fixed << nSel_bv[isam]      << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_bv[isam]);
      cout << endl;
      nTotal    	+= nSelv[isam];
      nTotalVar 	+= nSelVarv[isam];
      nTotal_i      	+= nSel_iv[isam];
      nTotalVar_i   	+= nSelVar_iv[isam];
      nTotal_nob      	+= nSel_nobv[isam];
      nTotalVar_nob   	+= nSelVar_nobv[isam];
      nTotal_vbf      	+= nSel_vbfv[isam];
      nTotalVar_vbf   	+= nSelVar_vbfv[isam];
      nTotal_b      	+= nSel_bv[isam];
      nTotalVar_b       += nSelVar_bv[isam];
    }
    cout << endl;
    cout << setw(10) << "MC";
    cout << setw(10) << setprecision(3) << fixed << nTotal     << " +/- " << sqrt(nTotalVar);
    cout << setw(10) << setprecision(3) << fixed << nTotal_i   << " +/- " << sqrt(nTotalVar_i);
    cout << setw(10) << setprecision(3) << fixed << nTotal_nob << " +/- " << sqrt(nTotalVar_nob);
    cout << setw(10) << setprecision(3) << fixed << nTotal_vbf << " +/- " << sqrt(nTotalVar_vbf);
    cout << setw(10) << setprecision(3) << fixed << nTotal_b   << " +/- " << sqrt(nTotalVar_b);
    cout << endl;
  }
  if(hasData) {
    cout << setw(10) << "Data";
    cout << setw(10) << setprecision(3) << fixed << nSelv[0]        << " +/- " << sqrt(nSelVarv[0]);
    cout << setw(10) << setprecision(3) << fixed << nSel_iv[0]      << " +/- " << sqrt(nSelVar_iv[0]);
    cout << setw(10) << setprecision(3) << fixed << nSel_nobv[0]    << " +/- " << sqrt(nSelVar_nobv[0]);
    cout << setw(10) << setprecision(3) << fixed << nSel_vbfv[0]    << " +/- " << sqrt(nSelVar_vbfv[0]);
    cout << setw(10) << setprecision(3) << fixed << nSel_bv[0]      << " +/- " << sqrt(nSelVar_bv[0]) << endl;
  }

  printf("\n\n%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s%10s\n","","tot","=2 jets","no bjets","mjj","n1*n2<0","deta",
	 "0 jets","1 jet","2 jets",">2 jets");  
  for(UInt_t isam=0;isam<samplev.size();isam++) {
    printf("%10s%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n",snamev[isam].Data(),countervv[isam][0],
	   countervv[isam][1]/countervv[isam][0],countervv[isam][2]/countervv[isam][0],
	   countervv[isam][3]/countervv[isam][0],countervv[isam][4]/countervv[isam][0],countervv[isam][5]/countervv[isam][0],
	   countervv[isam][6]/countervv[isam][0],countervv[isam][7]/countervv[isam][0],countervv[isam][8]/countervv[isam][0],
	   countervv[isam][9]/countervv[isam][0]);
  }

  makeHTML(outputDir);
  
  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  cout << endl;
  
  gBenchmark->Show("plotEmu");      
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/plots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>No Need For a Title</title></head>" << endl;
  htmlfile << "<body bgcolor=\"000000\">" << endl;
  htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">No Need For a Title</h3>" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt_mu.png\"><img src=\"plots/pt_mu.png\" alt=\"plots/pt_mu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt_e.png\"><img src=\"plots/pt_e.png\" alt=\"plots/pt_e.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta_mu.png\"><img src=\"plots/eta_mu.png\" alt=\"plots/eta_mu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta_e.png\"><img src=\"plots/eta_e.png\" alt=\"plots/eta_e.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phimu.png\"><img src=\"plots/phimu.png\" alt=\"plots/phimu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phiele.png\"><img src=\"plots/phiele.png\" alt=\"plots/phiele.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/dphi_emu.png\"><img src=\"plots/dphi_emu.png\" alt=\"plots/dphi_emu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt.png\"><img src=\"plots/pt.png\" alt=\"plots/pt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met.png\"><img src=\"plots/met.png\" alt=\"plots/met.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metraw.png\"><img src=\"plots/metraw.png\" alt=\"plots/metraw.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt.png\"><img src=\"plots/mt.png\" alt=\"plots/mt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metdphi.png\"><img src=\"plots/metdphi.png\" alt=\"plots/metdphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetamiss.png\"><img src=\"plots/pzetamiss.png\" alt=\"plots/pzetamiss.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetavis.png\"><img src=\"plots/pzetavis.png\" alt=\"plots/pzetavis.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetavar.png\"><img src=\"plots/pzetavar.png\" alt=\"plots/pzetavar.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/rawProjVar.png\"><img src=\"plots/rawProjVar.png\" alt=\"plots/rawProjVar.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetpt1.png\"><img src=\"plots/jetpt1.png\" alt=\"plots/jetpt1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetpt2.png\"><img src=\"plots/jetpt2.png\" alt=\"plots/jetpt2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jeteta1.png\"><img src=\"plots/jeteta1.png\" alt=\"plots/jeteta1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jeteta2.png\"><img src=\"plots/jeteta2.png\" alt=\"plots/jeteta2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/bjetpt.png\"><img src=\"plots/bjetpt.png\" alt=\"plots/bjetpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/bjeteta.png\"><img src=\"plots/bjeteta.png\" alt=\"plots/bjeteta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/bjetphi.png\"><img src=\"plots/bjetphi.png\" alt=\"plots/bjetphi.png\" width=\"100%\"></a></td>" << endl;

  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetdphi.png\"><img src=\"plots/jetdphi.png\" alt=\"plots/jetdphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mjj.png\"><img src=\"plots/mjj.png\" alt=\"plots/mjj.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetdeta.png\"><img src=\"plots/jetdeta.png\" alt=\"plots/jetdeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetEtaProd.png\"><img src=\"plots/jetEtaProd.png\" alt=\"plots/jetEtaProd.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/njets.png\"><img src=\"plots/njets.png\" alt=\"plots/njets.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nbjets.png\"><img src=\"plots/nbjets.png\" alt=\"plots/nbjets.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/btag.png\"><img src=\"plots/btag.png\" alt=\"plots/btag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nvertices_reweighted.png\"><img src=\"plots/nvertices_reweighted.png\" alt=\"plots/nvertices_reweighted.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;   
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_lept.png\"><img src=\"plots/vismass_lept.png\" alt=\"plots/vismass_lept.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_incl.png\"><img src=\"plots/vismass_incl.png\" alt=\"plots/vismass_incl.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class1.png\"><img src=\"plots/vismass_class1.png\" alt=\"plots/vismass_class1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class2a.png\"><img src=\"plots/vismass_class2a.png\" alt=\"plots/pt2_ee.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class2b.png\"><img src=\"plots/vismass_class2b.png\" alt=\"plots/vismass_class2b.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt1.png\"><img src=\"plots/pt1.png\" alt=\"plots/pt1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt2.png\"><img src=\"plots/pt2.png\" alt=\"plots/pt2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta1.png\"><img src=\"plots/eta1.png\" alt=\"plots/eta1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta2.png\"><img src=\"plots/eta2.png\" alt=\"plots/eta2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi1.png\"><img src=\"plots/phi1.png\" alt=\"plots/phi1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi2.png\"><img src=\"plots/phi2.png\" alt=\"plots/phi2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
            
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}
