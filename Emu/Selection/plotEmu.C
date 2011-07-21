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
      string json;
      Int_t type;
      stringstream ss(line);
      ss >> fname >> type >> json;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(type);
      samplev.back()->xsecv.push_back(0);
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

  // get indices of fakes and zmm, so we can add them together later
  UInt_t ifake=9999, izmm=9999;
  for(UInt_t isam=0; isam<snamev.size(); isam++) {
    if(snamev[isam].Contains("fakes")) ifake = isam;
    if(snamev[isam].Contains("zmm"))   izmm  = isam;
  }
  if(ifake>snamev.size() || izmm>snamev.size()) { cout << "error -- ifake: " << ifake << " izmm: " << izmm << endl << endl; return; }

  CPlot::sOutDir = outputDir + TString("/plots");
  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };  // final state type

  const Double_t kssFakeWgt = 1.312; // 595.8/454;    // ratio of fake-rate fakes to same-sign fakes
                                     // careful: this is hardcoded in selectEmu.C

  //
  // scale factors for each category
  //
  const Double_t kCat_novbf = 0.999;
  const Double_t kCat_vbf   = 1.401;
  const Double_t kCat_nob   = 0.999;
  const Double_t kCat_b     = 0.983;
  // no btag                      btag                       vbf                              no vbf
  // 0.999 +/- 0.008        0.983 +/- 0.064        1.401 +/- 0.411        0.999 /- 0.008

  const Double_t kExpB_b   = 0.97; // samples with expected b's, in the b-tag category
  const Double_t kExpB_nob = 1.26; // samples with expected b's, in the no b-tag category
  const Double_t kUnB_b    = 1.21; // samples with no expected b's, in the b-tag category
  const Double_t kUnB_nob  = 0.99; // samples with no expected b's, in the no b-tag category
                                   //   1a) [BTag w/ b's]  0.97
                                   //   1b) [non-BTag w/ b's]  1.26
                                   //   2a) [BTag w/o b's]  1.21
                                   //   2b) [non-BTag w/o b's]  0.99
  
  // vbf cut values
  const Double_t mjjMin   = 350;
  const Double_t dEtaMin  = 3.5;

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 
  
  const Double_t pi = 3.14159265358979;
  
  Bool_t hasData = (samplev[0]->fnamev.size()>0);
  
  //
  // Set up histograms
  //

  // after lepton selection
  vector<TH1F*> hMetv, hMetRawv;
  vector<TH2F*> hMassVpmissv;
  vector<TH1F*> hProjMetv, hProjVisv, hProjVis2v, hProjVarv,hRawProjVarv;
  vector<TH1F*> hNjetsv, hNbjetsv;
  vector<TH1F*> hBdiscrv;
  vector<TH1F*> hDPhiv, hMtv, hPtv, hLepDEtav;
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
  vector<TH1F*> hNPVv, hNPVrawv;                // primary vertexes

  vector<TH1F*> hMassv, hMassLv, hMassHighv;    // *L: plots for limits

  // inclusive
  vector<TH1F*> hMass_iv;
  vector<TH1F*> hMassL_iv;

  // Class novbf: one or fewer pt 30 jets, or exactly two pt 30 jets and !(vbf cuts)
  vector<TH1F*> hMass_novbfv;
  vector<TH1F*> hMassL_novbfv;

  // Class vbf: exactly two pt 30 jets and vbf cuts
  vector<TH1F*> hMass_vbfv;
  vector<TH1F*> hMassL_vbfv;

  // Class nob: one or fewer pt 30 jets, no b-tagged jet
  vector<TH1F*> hMass_nobv;
  vector<TH1F*> hMassL_nobv;

  // Class b: one or fewer pt 30 jets, at least one b-tag jet
  vector<TH1F*> hMass_bv;
  vector<TH1F*> hMassL_bv;

  vector<Double_t> nSelv,       nSelVarv;
  vector<Double_t> nSel_iv,     nSelVar_iv;
  vector<Double_t> nSel_novbfv, nSelVar_novbfv;
  vector<Double_t> nSel_vbfv,   nSelVar_vbfv;
  vector<Double_t> nSel_nobv,   nSelVar_nobv;
  vector<Double_t> nSel_bv,     nSelVar_bv;
  
  vector<Double_t> nUnskEventsv;

  vector<vector <Double_t> > countervv;
  vector<vector <Double_t> > vbfcountvv;

  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    // after lepton selection
    sprintf(hname,"hMet_%i",isam);     hMetv.push_back(new TH1F(hname,"",30,0,150));        hMetv[isam]->Sumw2();
    sprintf(hname,"hMetRaw_%i",isam);  hMetRawv.push_back(new TH1F(hname,"",30,0,150));     hMetRawv[isam]->Sumw2();
    sprintf(hname,"hProjMet_%i",isam);       hProjMetv.push_back(new TH1F(hname,"",30,-100,120));     hProjMetv[isam]->Sumw2();
    sprintf(hname,"hProjVis_%i",isam);       hProjVisv.push_back(new TH1F(hname,"",20,0,50));        hProjVisv[isam]->Sumw2();
    sprintf(hname,"hProjVis2_%i",isam);      hProjVis2v.push_back(new TH1F(hname,"",40,0,100));       hProjVis2v[isam]->Sumw2();
    sprintf(hname,"hProjVar_%i",isam);       hProjVarv.push_back(new TH1F(hname,"",30,-230,100));     hProjVarv[isam]->Sumw2();
    sprintf(hname,"hRawProjVar_%i",isam);    hRawProjVarv.push_back(new TH1F(hname,"",30,-230,100));  hRawProjVarv[isam]->Sumw2();
    sprintf(hname,"hNjets_%i",isam);         hNjetsv.push_back(new TH1F(hname,"",5,-0.5,4.5));      hNjetsv[isam]->Sumw2();
    sprintf(hname,"hNbjets_%i",isam);        hNbjetsv.push_back(new TH1F(hname,"",5,-0.5,4.5));     hNbjetsv[isam]->Sumw2();
    sprintf(hname,"hBdiscr_%i",isam);        hBdiscrv.push_back(new TH1F(hname,"",50,0,15));          hBdiscrv[isam]->Sumw2();
    sprintf(hname,"hDPhi_%i",isam);    hDPhiv.push_back(new TH1F(hname,"",36,0,180));       hDPhiv[isam]->Sumw2();
    sprintf(hname,"hMt_%i",isam);      hMtv.push_back(new TH1F(hname,"",30,0,210));         hMtv[isam]->Sumw2();
    sprintf(hname,"hPt_%i",isam);      hPtv.push_back(new TH1F(hname,"",30,0,120));         hPtv[isam]->Sumw2();
    sprintf(hname,"hLepDEta_%i",isam); hLepDEtav.push_back(new TH1F(hname,"",20,0,4));      hLepDEtav[isam]->Sumw2();
    sprintf(hname,"hMetDPhi_%i",isam); hMetDPhiv.push_back(new TH1F(hname,"",30,0,180));    hMetDPhiv[isam]->Sumw2();
    sprintf(hname,"hPt1_%i",isam);     hPt1v.push_back(new TH1F(hname,"",30,0,150));        hPt1v[isam]->Sumw2();
    sprintf(hname,"hEta1_%i",isam);    hEta1v.push_back(new TH1F(hname,"",20,-3,3));        hEta1v[isam]->Sumw2();
    sprintf(hname,"hPhi1_%i",isam);    hPhi1v.push_back(new TH1F(hname,"",20,-3.2,3.2));    hPhi1v[isam]->Sumw2();  
    sprintf(hname,"hPt2_%i",isam);     hPt2v.push_back(new TH1F(hname,"",60,0,150));        hPt2v[isam]->Sumw2();
    sprintf(hname,"hEta2_%i",isam);    hEta2v.push_back(new TH1F(hname,"",20,-3,3));        hEta2v[isam]->Sumw2();
    sprintf(hname,"hPhi2_%i",isam);    hPhi2v.push_back(new TH1F(hname,"",20,-3.2,3.2));    hPhi2v[isam]->Sumw2();
    sprintf(hname,"hPtMu_%i",isam);    hPtMuv.push_back(new TH1F(hname,"",50,0,100));       hPtMuv[isam]->Sumw2();
    sprintf(hname,"hEtaMu_%i",isam);   hEtaMuv.push_back(new TH1F(hname,"",20,-3,3));       hEtaMuv[isam]->Sumw2();
    sprintf(hname,"hPhiMu_%i",isam);   hPhiMuv.push_back(new TH1F(hname,"",20,-3.2,3.2));   hPhiMuv[isam]->Sumw2();  
    sprintf(hname,"hPtEle_%i",isam);   hPtElev.push_back(new TH1F(hname,"",30,0,100));      hPtElev[isam]->Sumw2();
    sprintf(hname,"hEtaEle_%i",isam);  hEtaElev.push_back(new TH1F(hname,"",20,-3,3));      hEtaElev[isam]->Sumw2();
    sprintf(hname,"hPhiEle_%i",isam);  hPhiElev.push_back(new TH1F(hname,"",20,-3.2,3.2));  hPhiElev[isam]->Sumw2();
    sprintf(hname,"hJetPt1_%i",isam);  hJetPt1v.push_back(new TH1F(hname,"",30,0,300));     hJetPt1v[isam]->Sumw2();
    sprintf(hname,"hJetEta1_%i",isam); hJetEta1v.push_back(new TH1F(hname,"",20,-5,5));     hJetEta1v[isam]->Sumw2();
    sprintf(hname,"hJetPt2_%i",isam);  hJetPt2v.push_back(new TH1F(hname,"",30,0,300));     hJetPt2v[isam]->Sumw2();
    sprintf(hname,"hJetEta2_%i",isam); hJetEta2v.push_back(new TH1F(hname,"",20,-5,5));     hJetEta2v[isam]->Sumw2();
    sprintf(hname,"hJetDPhi_%i",isam); hJetDPhiv.push_back(new TH1F(hname,"",30,0,180));    hJetDPhiv[isam]->Sumw2();
    sprintf(hname,"hMjj_%i",isam);     hMjjv.push_back(new TH1F(hname,"",50,0,1000));       hMjjv[isam]->Sumw2();
    sprintf(hname,"hDEta_%i",isam);    hDEtav.push_back(new TH1F(hname,"",20,0,8));         hDEtav[isam]->Sumw2();
    sprintf(hname,"hEtaProd_%i",isam); hEtaProdv.push_back(new TH1F(hname,"",30,-7.5,7.5)); hEtaProdv[isam]->Sumw2();
    sprintf(hname,"hBJetPt_%i",isam);  hBJetPtv.push_back(new TH1F(hname,"",30,0,150));     hBJetPtv[isam]->Sumw2();
    sprintf(hname,"hBJetEta_%i",isam); hBJetEtav.push_back(new TH1F(hname,"",20,-3,3));     hBJetEtav[isam]->Sumw2();
    sprintf(hname,"hBJetPhi_%i",isam); hBJetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hBJetPhiv[isam]->Sumw2();
    sprintf(hname,"hNPV_%i",isam);     hNPVv.push_back(new TH1F(hname,"",20,-0.5,19.5));    hNPVv[isam]->Sumw2();
    sprintf(hname,"hNPVraw_%i",isam);  hNPVrawv.push_back(new TH1F(hname,"",20,-0.5,19.5)); hNPVrawv[isam]->Sumw2();

    // lepton selection
    sprintf(hname,"hMass_%i",isam);    hMassv.push_back(new TH1F(hname,"",40,0,400));       hMassv[isam]->Sumw2();
    sprintf(hname,"hMassHigh_%i",isam);hMassHighv.push_back(new TH1F(hname,"",100,0,2000)); hMassHighv[isam]->Sumw2();
    sprintf(hname,"hMassL_%i",isam);   hMassLv.push_back(new TH1F(hname,"",50,0,500));      hMassLv[isam]->Sumw2();
    sprintf(hname,"hMassVpmiss_%i",isam);    hMassVpmissv.push_back(new TH2F(hname,"",50,0,500,50,-100,100));     hMassVpmissv[isam]->Sumw2();
    // inclusive
    sprintf(hname,"hMass_i_%i",isam);    hMass_iv.push_back(new TH1F(hname,"",40,0,400));     hMass_iv[isam]->Sumw2();
    sprintf(hname,"hMassL_i_%i",isam);   hMassL_iv.push_back(new TH1F(hname,"",50,0,500));    hMassL_iv[isam]->Sumw2();
    // no vbf
    sprintf(hname,"hMass_novbf_%i",isam);	  hMass_novbfv.push_back(new TH1F(hname,"",40,0,400));     hMass_novbfv[isam]->Sumw2();
    sprintf(hname,"hMassL_novbf_%i",isam);	  hMassL_novbfv.push_back(new TH1F(hname,"",50,0,500));    hMassL_novbfv[isam]->Sumw2();
    // vbf
    sprintf(hname,"hMass_vbf_%i",isam);		  hMass_vbfv.push_back(new TH1F(hname,"",20,0,400));       hMass_vbfv[isam]->Sumw2();
    sprintf(hname,"hMassL_vbf_%i",isam);	  hMassL_vbfv.push_back(new TH1F(hname,"",25,0,500));      hMassL_vbfv[isam]->Sumw2();
    // no b-tag
    sprintf(hname,"hMass_nob_%i",isam);		  hMass_nobv.push_back(new TH1F(hname,"",40,0,400));       hMass_nobv[isam]->Sumw2();
    sprintf(hname,"hMassL_nob_%i",isam);	  hMassL_nobv.push_back(new TH1F(hname,"",50,0,500));      hMassL_nobv[isam]->Sumw2();
    // b-tag
    sprintf(hname,"hMass_b_%i",isam);		  hMass_bv.push_back(new TH1F(hname,"",20,0,400));         hMass_bv[isam]->Sumw2();
    sprintf(hname,"hMassL_b_%i",isam);		  hMassL_bv.push_back(new TH1F(hname,"",25,0,500));        hMassL_bv[isam]->Sumw2();

    nSelv.push_back(0);       nSelVarv.push_back(0);
    nSel_iv.push_back(0);     nSelVar_iv.push_back(0);
    nSel_novbfv.push_back(0); nSelVar_novbfv.push_back(0);
    nSel_vbfv.push_back(0);   nSelVar_vbfv.push_back(0);
    nSel_nobv.push_back(0);   nSelVar_nobv.push_back(0);
    nSel_bv.push_back(0);     nSelVar_bv.push_back(0);

    nUnskEventsv.push_back(0);
    
    vector<Double_t> *counter = new vector<Double_t>;
    countervv.push_back(*counter);
    for(int i=0;i<30;i++)
      countervv[isam].push_back(0);
    vector<Double_t> *vbfcounter = new vector<Double_t>;
    vbfcountvv.push_back(*vbfcounter);
    for(int i=0;i<30;i++)
      vbfcountvv[isam].push_back(0);
    
  }

  EmuData data;
  Double_t rawMet=0,rawprojvar=0,npuWgt=1;
  UInt_t npt20jets;
  const UInt_t kMaxPt20Jets=35;
  TArrayF *btagArray = new TArrayF;
  btagArray->Set(kMaxPt20Jets);
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0; 
  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if((isam==0) && !hasData) continue;

    // expect b-jets in this sample?
    const Bool_t expB = snamev[isam].Contains("ttbar") || snamev[isam].Contains("bb_");
    
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
    eventTree->SetBranchAddress("rawMet",&rawMet);
    eventTree->SetBranchAddress("rawprojvar",&rawprojvar);
    if(isam!=0 && isam!=ifake)
      eventTree->SetBranchAddress("npuWgt",&npuWgt);

    assert(samplev[isam]->fnamev.size()>0);
    nUnskEventsv[isam] = lumi*samplev[isam]->xsecv[0];
    
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);

      Double_t wgt = 1;
      if(isam!=0 && isam!=ifake)
	wgt = data.weight*lumi;
      if(isam==ifake) {
	if(snamev[isam].Contains("ss-fakes")) { // same-sign fakes must have weight 1, but normalize them to fake-rate fakes
	  // if(lumi<500) {       // 190/pb
	  //   wgt = data.weight*(97.5/86.);    // pt10	      // wgt = data.weight*(63.6/52.); // pt15
	  // }
	  // else if(lumi<1000) { // 865/pb
	  //   wgt = data.weight*(447.8/343.);  // pt10     /// pt15: 	      wgt = data.weight*(283.8/192.);
	  // }
	  // else {              // 1084-1149/pb
	  wgt = data.weight*(kssFakeWgt);  // pt10
	  // }
	}
	else // fake-rate fakes have their weight stored in the ntuple
	  wgt = data.weight;
      }
	
      // get muon, electron kinematics
      Float_t mupt,mueta,muphi,elept,eleta,elephi;
      if(data.state==kMuEle) {
	mupt  = data.lpt1; mueta = data.leta1; muphi = data.lphi1; elept = data.lpt2; eleta = data.leta2; elephi = data.lphi2;
      } else {
	mupt  = data.lpt2; mueta = data.leta2; muphi = data.lphi2; elept = data.lpt1; eleta = data.leta1; elephi = data.lphi1;
      }

      // fill plots after lepton selection
      hMetv[isam]        ->Fill(data.met,    	wgt);
      hMetRawv[isam]     ->Fill(rawMet,      	wgt);
      hProjMetv[isam]    ->Fill(data.pmet,     	wgt);
      hProjVisv[isam]    ->Fill(data.pvis,     	wgt);
      hProjVis2v[isam]   ->Fill(data.pvis,     	wgt);
      hProjVarv[isam]    ->Fill(-0.85*data.pvis + data.pmet, wgt);
      hRawProjVarv[isam] ->Fill(-rawprojvar,    wgt);
      hNjetsv[isam]      ->Fill(data.njets,     wgt);
      hNbjetsv[isam]     ->Fill(data.nbjets,    wgt);

      assert(npt20jets<kMaxPt20Jets);
      for(UInt_t ib=0;ib<npt20jets;ib++)
	hBdiscrv[isam]   ->Fill((*btagArray)[ib], wgt);
      
      hDPhiv[isam]   ->Fill(data.dphi*180./pi,wgt);
      hMtv[isam]     ->Fill(data.mt,     wgt);
      hPtv[isam]     ->Fill(data.pt,     wgt);
      hLepDEtav[isam]->Fill(data.leta1-data.leta2,     wgt);
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
      // cout << npuWgt << endl;
      // if(ientry>10) break;
      if(npuWgt!=0)
	hNPVrawv[isam] ->Fill(data.nPV,    (isam==0 || isam==ifake) ? wgt : wgt/npuWgt);

      hMassv[isam]       ->Fill(data.mass,   	wgt);
      hMassLv[isam]      ->Fill(data.mass,   	wgt);      
      hMassHighv[isam]   ->Fill(data.mass,   	wgt);
      hMassVpmissv[isam] ->Fill(data.mass, data.pmet, wgt);

      nSelv[isam]    += wgt;
      nSelVarv[isam] += wgt*wgt;

      if(0.85*data.pvis - data.pmet > 25) continue;

      // inclusive
      hMass_iv[isam]    ->Fill(data.mass,   wgt);
      hMassL_iv[isam]   ->Fill(data.mass,   wgt);
      nSel_iv[isam]     += wgt;
      nSelVar_iv[isam]  += wgt*wgt;

      Bool_t vbfcuts = (data.mjj > mjjMin) && (data.jeta1*data.jeta2 < 0) && (fabs(data.jeta1-data.jeta2) > dEtaMin);

      // no vbf
      if((data.njets < 2) || (data.njets == 2 && !vbfcuts)) {
      	hMass_novbfv[isam]   ->Fill(data.mass,   (isam==0) ? wgt : kCat_novbf*wgt);      
      	hMassL_novbfv[isam]  ->Fill(data.mass,   (isam==0) ? wgt : kCat_novbf*wgt);      
      	nSel_novbfv[isam]     += (isam==0) ? wgt : kCat_novbf*wgt;
      	nSelVar_novbfv[isam]  += (isam==0) ? wgt : kCat_novbf*wgt*kCat_novbf*wgt;
      }
      // vbf
      if(data.njets==2 && vbfcuts) {
      	hMass_vbfv[isam]   ->Fill(data.mass,   (isam==0) ? wgt : kCat_vbf*wgt);      
      	hMassL_vbfv[isam]  ->Fill(data.mass,   (isam==0) ? wgt : kCat_vbf*wgt);      
      	nSel_vbfv[isam]     += (isam==0) ? wgt : kCat_vbf*wgt;
      	nSelVar_vbfv[isam]  += (isam==0) ? wgt : kCat_vbf*wgt*kCat_vbf*wgt;
      }
      Double_t bscale;
      // no b-tag
      if(data.njets<=1 && data.nbjets==0) {
	if(expB) bscale = kExpB_nob;
	else     bscale = kUnB_nob;
	hMass_nobv[isam]   ->Fill(data.mass,   (isam==0) ? wgt : bscale*kCat_nob*wgt);      
	hMassL_nobv[isam]  ->Fill(data.mass,   (isam==0) ? wgt : bscale*kCat_nob*wgt);      
	nSel_nobv[isam]     += (isam==0) ? wgt : bscale*kCat_nob*wgt;
	nSelVar_nobv[isam]  += (isam==0) ? wgt : bscale*kCat_nob*wgt*kCat_nob*wgt;
      }
      // b-tag
      if(data.njets<=1 && data.nbjets>0) {
	if(expB) bscale = kExpB_b;
	else     bscale = kUnB_b;
	hMass_bv[isam]   ->Fill(data.mass,   (isam==0) ? wgt : bscale*kCat_b*wgt);      
	hMassL_bv[isam]  ->Fill(data.mass,   (isam==0) ? wgt : bscale*kCat_b*wgt);      
	nSel_bv[isam]     += (isam==0) ? wgt : bscale*kCat_b*wgt;
	nSelVar_bv[isam]  += (isam==0) ? wgt : bscale*kCat_b*wgt*kCat_b*wgt;
      }
    }
    delete infile;
    infile=0, eventTree=0;
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  TCanvas *c = MakeCanvas("c","c",700,700);
  
  // string buffers
  char ylabel[100];   // y-axis label
  char lumitext[3][50];
  if(lumi>0) {
    sprintf(lumitext[0],"#bf{CMS Preliminary}");
    sprintf(lumitext[1],"#bf{%.1f fb^{-1}}",lumi/1000.);
    sprintf(lumitext[2],"#bf{#tau_{e}#tau_{mu}}");
  }
  
  //
  // Begin plots:
  //

  assert(hMetv[ifake] && hMetv[izmm]); // have to change some things if you want to plot without these samples

  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMetv[0]->GetBinWidth(1));
  CPlot plotMet("met","","#bf{met [GeV]}",ylabel);
  if(hMetv[ifake] && hMetv[izmm]) hMetv[ifake]->Add(hMetv[izmm]);
  if(hasData) { plotMet.AddHist1D(hMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotMet.AddToStack(hMetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotMet.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotMet.SetLegend(0.7,0.65,0.9,0.9);
  plotMet.Draw(c,kTRUE,format);

  // stack up the raw met hists
  TH1F *hTmpMetStack = 0;
  if(hMetRawv[ifake] && hMetRawv[izmm]) hMetRawv[ifake]->Add(hMetRawv[izmm]);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==1) {
      hTmpMetStack = new TH1F(*hMetv[isam]);
      hTmpMetStack->Reset();
    }
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    hTmpMetStack->Add(hMetRawv[isam]);
  }

  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMetv[0]->GetBinWidth(1));
  CPlot plotMetRaw("metraw","","#bf{raw met [GeV]}",ylabel);
  // if(hMetRawv[ifake] && hMetRawv[izmm]) hMetRawv[ifake]->Add(hMetRawv[izmm]);
  if(hasData) { plotMetRaw.AddHist1D(hMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotMetRaw.AddToStack(hMetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMetRaw.AddHist1D(hTmpMetStack,"#bf{un-corrected}","hist",kRed);
  if(lumi>0) plotMetRaw.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotMetRaw.SetLegend(0.7,0.65,0.9,0.9);
  plotMetRaw.Draw(c,kTRUE,format);

  // projection variables
  sprintf(ylabel,"#bf{Events / %.2f}",hProjMetv[0]->GetBinWidth(1));
  CPlot plotProjMet("pzetamiss","","#bf{#slash{p}_{#zeta} [GeV]}",ylabel);
  if(hProjMetv[ifake] && hProjMetv[izmm]) hProjMetv[ifake]->Add(hProjMetv[izmm]);
  if(hasData) { plotProjMet.AddHist1D(hProjMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotProjMet.AddToStack(hProjMetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotProjMet.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotProjMet.SetLegend(0.7,0.65,0.9,0.9);
  plotProjMet.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hProjVisv[0]->GetBinWidth(1));
  CPlot plotProjVis("pzetavis","","#bf{p_{#zeta}^{vis} [GeV]}",ylabel);
  if(hProjVisv[ifake] && hProjVisv[izmm]) hProjVisv[ifake]->Add(hProjVisv[izmm]);
  if(hasData) { plotProjVis.AddHist1D(hProjVisv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotProjVis.AddToStack(hProjVisv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotProjVis.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotProjVis.SetLegend(0.7,0.65,0.9,0.9);
  plotProjVis.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hProjVis2v[0]->GetBinWidth(1));
  CPlot plotProjVis2("pzetavis2","","#bf{p_{#zeta}^{vis} [GeV]}",ylabel);
  if(hProjVis2v[ifake] && hProjVis2v[izmm]) hProjVis2v[ifake]->Add(hProjVis2v[izmm]);
  if(hasData) { plotProjVis2.AddHist1D(hProjVis2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotProjVis2.AddToStack(hProjVis2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotProjVis2.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotProjVis2.SetLegend(0.7,0.65,0.9,0.9);
  plotProjVis2.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hProjVarv[0]->GetBinWidth(1));
  CPlot plotProjVar("pzetavar","","#bf{#slash{p}_{#zeta} - 0.85*p_{#zeta}^{vis} [GeV]}",ylabel);
  if(hProjVarv[ifake] && hProjVarv[izmm]) hProjVarv[ifake]->Add(hProjVarv[izmm]);
  if(hasData) { plotProjVar.AddHist1D(hProjVarv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotProjVar.AddToStack(hProjVarv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotProjVar.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotProjVar.SetLegend(0.21,0.45,0.41,0.7);
  plotProjVar.Draw(c,kTRUE,format);

  // stack up the raw projection variable hists
  TH1F *hTmpStack = 0;
  if(hRawProjVarv[ifake] && hRawProjVarv[izmm]) hRawProjVarv[ifake]->Add(hRawProjVarv[izmm]);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==1) {
      hTmpStack = new TH1F(*hProjVarv[isam]);
      hTmpStack->Reset();
    }
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    hTmpStack->Add(hRawProjVarv[isam]);
  }

  // then *re*plot the corrected variable, but with the raw one overlaid now
  sprintf(ylabel,"#bf{Events / %.2f}",hProjVarv[0]->GetBinWidth(1));
  CPlot plotRawProjVar("rawProjVar","","#bf{#slash{p}_{#zeta} - 0.85*p_{#zeta}^{vis} [GeV]}",ylabel);
  // if(hXXv[ifake] && hXXv[izmm]) hXXv[ifake]->Add(hXXv[izmm]); // don't add it in twice!
  if(hasData) { plotRawProjVar.AddHist1D(hProjVarv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==izmm) continue;
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    plotRawProjVar.AddToStack(hProjVarv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotRawProjVar.AddHist1D(hTmpStack,"#bf{un-corrected}","hist",kRed);
  if(lumi>0) plotRawProjVar.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotRawProjVar.SetLegend(0.21,0.45,0.41,0.7);
  plotRawProjVar.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hNjetsv[0]->GetBinWidth(1));
  CPlot plotNjets("njets","","#bf{number of jets}",ylabel);
  if(hNjetsv[ifake] && hNjetsv[izmm]) hNjetsv[ifake]->Add(hNjetsv[izmm]);
  if(hasData) { plotNjets.AddHist1D(hNjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotNjets.AddToStack(hNjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotNjets.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotNjets.SetYRange(0,1.5*(plotNjets.GetStack()->GetMaximum()));
  plotNjets.SetLegend(0.7,0.65,0.9,0.9);
  plotNjets.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hNbjetsv[0]->GetBinWidth(1));
  CPlot plotNbjets("nbjets","","#bf{number of b-jets}",ylabel);
  if(hNbjetsv[ifake] && hNbjetsv[izmm]) hNbjetsv[ifake]->Add(hNbjetsv[izmm]);
  if(hasData) { plotNbjets.AddHist1D(hNbjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotNbjets.AddToStack(hNbjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotNbjets.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotNbjets.SetLegend(0.7,0.65,0.9,0.9);
  plotNbjets.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hBdiscrv[0]->GetBinWidth(1));
  CPlot plotBdiscr("btag","","#bf{b-tag discr.}",ylabel);
  if(hBdiscrv[ifake] && hBdiscrv[izmm]) hBdiscrv[ifake]->Add(hBdiscrv[izmm]);
  if(hasData) { plotBdiscr.AddHist1D(hBdiscrv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotBdiscr.AddToStack(hBdiscrv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotBdiscr.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotBdiscr.SetLegend(0.7,0.65,0.9,0.9);
  plotBdiscr.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hDPhiv[0]->GetBinWidth(1));
  CPlot plotDPhi("dphi_emu","","#bf{#Delta^{}#phi_{emu} [deg]}",ylabel);
  if(hDPhiv[ifake] && hDPhiv[izmm]) hDPhiv[ifake]->Add(hDPhiv[izmm]);
  if(hasData) { plotDPhi.AddHist1D(hDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotDPhi.AddToStack(hDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotDPhi.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotDPhi.TransLegend(-0.15,0);
  assert(plotDPhi.GetStack());
  plotDPhi.SetLegend(0.7,0.65,0.9,0.9);
  plotDPhi.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMtv[0]->GetBinWidth(1));
  CPlot plotMt("mt","","#bf{m_{T}(ll,#slash{E}_{T}) [GeV]}",ylabel);
  if(hMtv[ifake] && hMtv[izmm]) hMtv[ifake]->Add(hMtv[izmm]);
  if(hasData) { plotMt.AddHist1D(hMtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotMt.AddToStack(hMtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotMt.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotMt.SetLegend(0.7,0.65,0.9,0.9);
  plotMt.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.1f GeV}",hPtv[0]->GetBinWidth(1));
  CPlot plotPt("pt","","#bf{p_{T}^{ll} [GeV]}",ylabel);
  if(hPtv[ifake] && hPtv[izmm]) hPtv[ifake]->Add(hPtv[izmm]);
  if(hasData) { plotPt.AddHist1D(hPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPt.AddToStack(hPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPt.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPt.SetLegend(0.7,0.65,0.9,0.9);
  plotPt.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.1f GeV}",hLepDEtav[0]->GetBinWidth(1));
  CPlot plotLepDEta("lepdeta","","#bf{#Delta^{}#eta(ll)}",ylabel);
  if(hLepDEtav[ifake] && hLepDEtav[izmm]) hLepDEtav[ifake]->Add(hLepDEtav[izmm]);
  if(hasData) { plotLepDEta.AddHist1D(hLepDEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotLepDEta.AddToStack(hLepDEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotLepDEta.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotLepDEta.SetLegend(0.7,0.65,0.9,0.9);
  plotLepDEta.Draw(c,kTRUE,format);

  // met 
  sprintf(ylabel,"#bf{Events / %.2f}",hMetDPhiv[0]->GetBinWidth(1));
  CPlot plotMetDPhi("metdphi","","#bf{#Delta^{}#phi(ll,#slash{E}_{T}) [deg]}",ylabel);
  if(hMetDPhiv[ifake] && hMetDPhiv[izmm]) hMetDPhiv[ifake]->Add(hMetDPhiv[izmm]);
  if(hasData) { plotMetDPhi.AddHist1D(hMetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotMetDPhi.AddToStack(hMetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotMetDPhi.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotMetDPhi.SetLegend(0.21,0.45,0.41,0.7);
  plotMetDPhi.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.1f GeV}",hPt1v[0]->GetBinWidth(1));
  CPlot plotPt1("pt1","","#bf{leading lepton p_{T} [GeV]}",ylabel);
  if(hPt1v[ifake] && hPt1v[izmm]) hPt1v[ifake]->Add(hPt1v[izmm]);
  if(hasData) { plotPt1.AddHist1D(hPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPt1.AddToStack(hPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPt1.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPt1.SetLegend(0.7,0.65,0.9,0.9);
  plotPt1.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.2f}",hEta1v[0]->GetBinWidth(1));
  CPlot plotEta1("eta1","","#bf{leading lepton #eta}",ylabel);
  if(hEta1v[ifake] && hEta1v[izmm]) hEta1v[ifake]->Add(hEta1v[izmm]);
  if(hasData) { plotEta1.AddHist1D(hEta1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotEta1.AddToStack(hEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotEta1.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotEta1.SetYRange(0,2.0*(plotEta1.GetStack()->GetMaximum()));
  plotEta1.SetLegend(0.7,0.65,0.9,0.9);
  plotEta1.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.2f}",hPhi1v[0]->GetBinWidth(1));
  CPlot plotPhi1("phi1","","#bf{leading lepton #phi}",ylabel);
  if(hPhi1v[ifake] && hPhi1v[izmm]) hPhi1v[ifake]->Add(hPhi1v[izmm]);
  if(hasData) { plotPhi1.AddHist1D(hPhi1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPhi1.AddToStack(hPhi1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPhi1.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPhi1.SetYRange(0,2.6*(plotPhi1.GetStack()->GetMaximum()));
  plotPhi1.SetLegend(0.7,0.65,0.9,0.9);
  plotPhi1.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.1f GeV}",hPt2v[0]->GetBinWidth(1));
  CPlot plotPt2("pt2","","#bf{trailing lepton p_{T} [GeV]}",ylabel);
  if(hPt2v[ifake] && hPt2v[izmm]) hPt2v[ifake]->Add(hPt2v[izmm]);
  if(hasData) { plotPt2.AddHist1D(hPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPt2.AddToStack(hPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPt2.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPt2.SetLegend(0.7,0.65,0.9,0.9);
  plotPt2.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.2f}",hEta2v[0]->GetBinWidth(1));
  CPlot plotEta2("eta2","","#bf{trailing lepton #eta}",ylabel);
  if(hEta2v[ifake] && hEta2v[izmm]) hEta2v[ifake]->Add(hEta2v[izmm]);
  if(hasData) { plotEta2.AddHist1D(hEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotEta2.AddToStack(hEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotEta2.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotEta2.SetYRange(0,2.0*(plotEta2.GetStack()->GetMaximum()));
  plotEta2.SetLegend(0.7,0.65,0.9,0.9);
  plotEta2.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.2f}",hPhi2v[0]->GetBinWidth(1));
  CPlot plotPhi2("phi2","","#bf{trailing lepton #phi}",ylabel);
  if(hPhi2v[ifake] && hPhi2v[izmm]) hPhi2v[ifake]->Add(hPhi2v[izmm]);
  if(hasData) { plotPhi2.AddHist1D(hPhi2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPhi2.AddToStack(hPhi2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPhi2.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPhi2.SetYRange(0,2.6*(plotPhi2.GetStack()->GetMaximum()));
  plotPhi2.SetLegend(0.7,0.65,0.9,0.9);
  plotPhi2.Draw(c,kTRUE,format);

  // mu / electron kinematics
  sprintf(ylabel,"#bf{Events / %.1f GeV}",hPtMuv[0]->GetBinWidth(1));
  CPlot plotPtMu("pt_mu","","#bf{mu p_{T} [GeV]}",ylabel);
  if(hPtMuv[ifake] && hPtMuv[izmm]) hPtMuv[ifake]->Add(hPtMuv[izmm]);
  if(hasData) { plotPtMu.AddHist1D(hPtMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPtMu.AddToStack(hPtMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPtMu.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPtMu.SetLegend(0.7,0.65,0.9,0.9);
  plotPtMu.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.2f}",hEtaMuv[0]->GetBinWidth(1));
  CPlot plotEtaMu("eta_mu","","#bf{mu #eta}",ylabel);
  if(hEtaMuv[ifake] && hEtaMuv[izmm]) hEtaMuv[ifake]->Add(hEtaMuv[izmm]);
  if(hasData) { plotEtaMu.AddHist1D(hEtaMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotEtaMu.AddToStack(hEtaMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotEtaMu.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotEtaMu.SetYRange(0,2.0*(plotEtaMu.GetStack()->GetMaximum()));
  plotEtaMu.SetLegend(0.7,0.65,0.9,0.9);
  plotEtaMu.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.2f}",hPhiMuv[0]->GetBinWidth(1));
  CPlot plotPhiMu("phimu","","#bf{mu #phi}",ylabel);
  if(hPhiMuv[ifake] && hPhiMuv[izmm]) hPhiMuv[ifake]->Add(hPhiMuv[izmm]);
  if(hasData) { plotPhiMu.AddHist1D(hPhiMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPhiMu.AddToStack(hPhiMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPhiMu.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPhiMu.SetYRange(0,2.6*(plotPhiMu.GetStack()->GetMaximum()));
  plotPhiMu.SetLegend(0.7,0.65,0.9,0.9);
  plotPhiMu.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.1f GeV/c}",hPtElev[0]->GetBinWidth(1));
  CPlot plotPtEle("pt_e","","#bf{ele p_{T} [GeV]}",ylabel);
  if(hPtElev[ifake] && hPtElev[izmm]) hPtElev[ifake]->Add(hPtElev[izmm]);
  if(hasData) { plotPtEle.AddHist1D(hPtElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPtEle.AddToStack(hPtElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPtEle.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPtEle.SetLegend(0.7,0.65,0.9,0.9);
  plotPtEle.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.2f}",hEtaElev[0]->GetBinWidth(1));
  CPlot plotEtaEle("eta_e","","#bf{ele #eta}",ylabel);
  if(hEtaElev[ifake] && hEtaElev[izmm]) hEtaElev[ifake]->Add(hEtaElev[izmm]);
  if(hasData) { plotEtaEle.AddHist1D(hEtaElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotEtaEle.AddToStack(hEtaElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotEtaEle.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotEtaEle.SetYRange(0,2.0*(plotEtaEle.GetStack()->GetMaximum()));
  plotEtaEle.SetLegend(0.7,0.65,0.9,0.9);
  plotEtaEle.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"#bf{Events / %.2f}",hPhiElev[0]->GetBinWidth(1));
  CPlot plotPhiEle("phiele","","#bf{ele #phi}",ylabel);
  if(hPhiElev[ifake] && hPhiElev[izmm]) hPhiElev[ifake]->Add(hPhiElev[izmm]);
  if(hasData) { plotPhiEle.AddHist1D(hPhiElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotPhiEle.AddToStack(hPhiElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotPhiEle.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotPhiEle.SetYRange(0,2.6*(plotPhiEle.GetStack()->GetMaximum()));
  plotPhiEle.SetLegend(0.7,0.65,0.9,0.9);
  plotPhiEle.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hJetPt1v[0]->GetBinWidth(1));
  CPlot plotJetPt1("jetpt1","","#bf{leading jet pt [GeV]}",ylabel);
  if(hJetPt1v[ifake] && hJetPt1v[izmm]) hJetPt1v[ifake]->Add(hJetPt1v[izmm]);
  if(hasData) { plotJetPt1.AddHist1D(hJetPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotJetPt1.AddToStack(hJetPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotJetPt1.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotJetPt1.SetLegend(0.7,0.65,0.9,0.9);
  plotJetPt1.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hJetPt2v[0]->GetBinWidth(1));
  CPlot plotJetPt2("jetpt2","","#bf{second jet pt [GeV]}",ylabel);
  if(hJetPt2v[ifake] && hJetPt2v[izmm]) hJetPt2v[ifake]->Add(hJetPt2v[izmm]);
  if(hasData) { plotJetPt2.AddHist1D(hJetPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotJetPt2.AddToStack(hJetPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotJetPt2.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotJetPt2.SetLegend(0.7,0.65,0.9,0.9);
  plotJetPt2.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hJetEta1v[0]->GetBinWidth(1));
  CPlot plotJetEta1("jeteta1","","#bf{leading jet #eta}",ylabel);
  if(hJetEta1v[ifake] && hJetEta1v[izmm]) hJetEta1v[ifake]->Add(hJetEta1v[izmm]);
  if(hasData) { plotJetEta1.AddHist1D(hJetEta1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotJetEta1.AddToStack(hJetEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotJetEta1.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotJetEta1.SetLegend(0.7,0.65,0.9,0.9);
  plotJetEta1.SetYRange(0,2.0*(plotJetEta1.GetStack()->GetMaximum()));
  plotJetEta1.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hJetEta2v[0]->GetBinWidth(1));
  CPlot plotJetEta2("jeteta2","","#bf{second jet #eta}",ylabel);
  if(hJetEta2v[ifake] && hJetEta2v[izmm]) hJetEta2v[ifake]->Add(hJetEta2v[izmm]);
  if(hasData) { plotJetEta2.AddHist1D(hJetEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotJetEta2.AddToStack(hJetEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotJetEta2.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotJetEta2.SetLegend(0.7,0.65,0.9,0.9);
  plotJetEta2.SetYRange(0,2.0*(plotJetEta2.GetStack()->GetMaximum()));
  plotJetEta2.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hJetDPhiv[0]->GetBinWidth(1));
  CPlot plotJetDPhi("jetdphi","","#bf{#Delta#phi (jj)}",ylabel);
  if(hJetDPhiv[ifake] && hJetDPhiv[izmm]) hJetDPhiv[ifake]->Add(hJetDPhiv[izmm]);
  if(hasData) { plotJetDPhi.AddHist1D(hJetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotJetDPhi.AddToStack(hJetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotJetDPhi.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotJetDPhi.SetLegend(0.21,0.45,0.41,0.7);
  plotJetDPhi.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hMjjv[0]->GetBinWidth(1));
  CPlot plotMjj("mjj","","#bf{dijet mass [GeV]}",ylabel);
  if(hMjjv[ifake] && hMjjv[izmm]) hMjjv[ifake]->Add(hMjjv[izmm]);
  if(hasData) { plotMjj.AddHist1D(hMjjv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotMjj.AddToStack(hMjjv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotMjj.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotMjj.SetLegend(0.7,0.65,0.9,0.9);
  plotMjj.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hDEtav[0]->GetBinWidth(1));
  CPlot plotDEta("jetdeta","","#bf{#Delta#eta (jj)}",ylabel);
  if(hDEtav[ifake] && hDEtav[izmm]) hDEtav[ifake]->Add(hDEtav[izmm]);
  if(hasData) { plotDEta.AddHist1D(hDEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotDEta.AddToStack(hDEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotDEta.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotDEta.SetLegend(0.7,0.65,0.9,0.9);
  plotDEta.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hEtaProdv[0]->GetBinWidth(1));
  CPlot plotEtaProd("jetEtaProd","","#bf{(#eta1*#eta2)(jj)}",ylabel);
  if(hEtaProdv[ifake] && hEtaProdv[izmm]) hEtaProdv[ifake]->Add(hEtaProdv[izmm]);
  if(hasData) { plotEtaProd.AddHist1D(hEtaProdv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotEtaProd.AddToStack(hEtaProdv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotEtaProd.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotEtaProd.SetLegend(0.7,0.65,0.9,0.9);
  plotEtaProd.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hBJetPtv[0]->GetBinWidth(1));
  CPlot plotBJetPt("bjetpt","","#bf{lead b-jet pt [GeV]}",ylabel);
  if(hBJetPtv[ifake] && hBJetPtv[izmm]) hBJetPtv[ifake]->Add(hBJetPtv[izmm]);
  if(hasData) { plotBJetPt.AddHist1D(hBJetPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotBJetPt.AddToStack(hBJetPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotBJetPt.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotBJetPt.SetYRange(0,2.0*(plotBJetPt.GetStack()->GetMaximum()));
  plotBJetPt.SetLegend(0.7,0.65,0.9,0.9);
  plotBJetPt.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hBJetEtav[0]->GetBinWidth(1));
  CPlot plotBJetEta("bjeteta","","#bf{lead b-jet #eta}",ylabel);
  if(hBJetEtav[ifake] && hBJetEtav[izmm]) hBJetEtav[ifake]->Add(hBJetEtav[izmm]);
  if(hasData) { plotBJetEta.AddHist1D(hBJetEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotBJetEta.AddToStack(hBJetEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotBJetEta.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotBJetEta.SetYRange(0,2.0*(plotBJetEta.GetStack()->GetMaximum()));
  plotBJetEta.SetLegend(0.7,0.65,0.9,0.9);
  plotBJetEta.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.2f}",hBJetPhiv[0]->GetBinWidth(1));
  CPlot plotBJetPhi("bjetphi","","#bf{lead b-jet #phi}",ylabel);
  if(hBJetPhiv[ifake] && hBJetPhiv[izmm]) hBJetPhiv[ifake]->Add(hBJetPhiv[izmm]);
  if(hasData) { plotBJetPhi.AddHist1D(hBJetPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotBJetPhi.AddToStack(hBJetPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotBJetPhi.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotBJetPhi.SetYRange(0,2.1*(plotBJetPhi.GetStack()->GetMaximum()));
  plotBJetPhi.SetLegend(0.7,0.65,0.9,0.9);
  plotBJetPhi.Draw(c,kTRUE,format);

  CPlot plotNPV("nvertices_reweighted","","#bf{N_{PV}}","#bf{Events}");
  if(hNPVv[ifake] && hNPVv[izmm]) hNPVv[ifake]->Add(hNPVv[izmm]);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotNPV.AddToStack(hNPVv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(hasData) { plotNPV.AddHist1D(hNPVv[0],samplev[0]->label,"E"); }
  if(lumi>0) plotNPV.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotNPV.SetLegend(0.7,0.65,0.9,0.9);
  plotNPV.Draw(c,kTRUE,format);

  CPlot plotNPVraw("nvertices_raw","","#bf{N_{PV}}","#bf{Events}");
  if(hNPVrawv[ifake] && hNPVrawv[izmm]) hNPVrawv[ifake]->Add(hNPVrawv[izmm]);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if(isam==izmm) continue;
    plotNPVraw.AddToStack(hNPVrawv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(hasData) { plotNPVraw.AddHist1D(hNPVrawv[0],samplev[0]->label,"E"); }
  if(lumi>0) plotNPVraw.AddTextBox(0.6,0.55,0.9,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotNPVraw.SetLegend(0.7,0.65,0.9,0.9);
  plotNPVraw.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMassv[0]->GetBinWidth(1));
  CPlot plotMass("vismass_lept","","#bf{m_{vis} [GeV]}",ylabel);
  if(hMassv[ifake] && hMassv[izmm]) hMassv[ifake]->Add(hMassv[izmm]);
  if(hasData) { plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==izmm) continue;
    if(snamev[isam].Contains("htt_gg")) plotMass.AddToStack(hMassv[isam],"",samplev[isam]->color,0,1,2);
    else if(snamev[isam].Contains("htt_bb")) plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(lumi>0) plotMass.AddTextBox(0.55,0.75,0.85,0.895,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]); //(0.6,0.754,0.93,0.9);
  // plotMass.AddTextBox("#bf{m_{A}=120, tan#beta=20}",0.737,0.537,0.897,0.637,0);
  plotMass.SetLegend(0.55,0.4,0.9,0.65);
  plotMass.Draw(c,kTRUE,format);

  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMassHighv[0]->GetBinWidth(1));
  CPlot plotMassHigh("vismass_lept-high","","m_{vis} [GeV]}",ylabel);
  if(hMassHighv[ifake] && hMassHighv[izmm]) hMassHighv[ifake]->Add(hMassHighv[izmm]);
  if(hasData) { plotMassHigh.AddHist1D(hMassHighv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==izmm) continue;
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    plotMassHigh.AddToStack(hMassHighv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassHigh.TransLegend(0.1,0);
  plotMassHigh.SetLogy();
  if(lumi>0) plotMassHigh.AddTextBox(0.71,0.45,0.91,0.4,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]);
  plotMass.SetLegend(0.7,0.65,0.9,0.9);
  plotMassHigh.Draw(c,kTRUE,format);

  CPlot plotMassVpmiss("massvpmiss","","m_{vis} [GeV]}","#slash{p}_{#zeta} [GeV]}");
  assert(hMassVpmissv[0]);
  plotMassVpmiss.AddHist2D((TH2D*)hMassVpmissv[0],"surf",kWhite,kBlue);

  //----------------------------------------------------------------------------------------
  // inclusive
  //----------------------------------------------------------------------------------------  
  
  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMass_iv[0]->GetBinWidth(1));
  CPlot plotMass_i("vismass_incl","","#bf{m_{vis} [GeV]}",ylabel);
  if(hMass_iv[ifake] && hMass_iv[izmm]) hMass_iv[ifake]->Add(hMass_iv[izmm]);
  if(hasData) { plotMass_i.AddHist1D(hMass_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==izmm) continue;
    if(snamev[isam].Contains("htt_gg")) plotMass_i.AddToStack(hMass_iv[isam],"",samplev[isam]->color,0,1,2);
    else if(snamev[isam].Contains("htt_bb")) plotMass_i.AddToStack(hMass_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_i.AddToStack(hMass_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_i.SetLegend(0.55,0.4,0.9,0.65);
  if(lumi>0) plotMass_i.AddTextBox(0.55,0.75,0.85,0.895,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]); //(0.6,0.754,0.93,0.9);
  // plotMass.AddTextBox("#bf{m_{A}=120, tan#beta=20}",0.737,0.537,0.897,0.637,0);
  plotMass.SetLegend(0.55,0.4,0.9,0.65);
  plotMass_i.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // no vbf
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMass_novbfv[0]->GetBinWidth(1));
  CPlot plotMass_novbf("vismass_class_novbf","","#bf{m_{vis} [GeV]}",ylabel);
  if(hMass_novbfv[ifake] && hMass_novbfv[izmm]) hMass_novbfv[ifake]->Add(hMass_novbfv[izmm]);
  if(hasData) { plotMass_novbf.AddHist1D(hMass_novbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==izmm) continue;
    if(snamev[isam].Contains("htt_gf")) {
      hMass_novbfv[isam]->Scale(10.);
      plotMass_novbf.AddToStack(hMass_novbfv[isam],"#bf{(10x)}"+samplev[isam]->label,samplev[isam]->color,1,1,3);
    }
    else if(snamev[isam].Contains("htt_vbf")) {
      hMass_novbfv[isam]->Scale(10.);
      plotMass_novbf.AddToStack(hMass_novbfv[isam],"",samplev[isam]->color,0,1,2);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_novbf.AddToStack(hMass_novbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_novbf.SetLegend(0.7,0.65,0.9,0.9);
  if(lumi>0) plotMass_novbf.AddTextBox(0.55,0.75,0.85,0.895,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]); //(0.6,0.754,0.93,0.9);
  // plotMass_novbf.AddTextBox("#bf{m_{H}=120}",0.737,0.537,0.837,0.637,0);
  plotMass_novbf.SetLegend(0.55,0.4,0.9,0.65);
  plotMass_novbf.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // vbf
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMass_vbfv[0]->GetBinWidth(1));
  CPlot plotMass_vbf("vismass_class_vbf","","#bf{m_{vis} [GeV]}",ylabel);
  if(hMass_vbfv[ifake] && hMass_vbfv[izmm]) hMass_vbfv[ifake]->Add(hMass_vbfv[izmm]);
  if(hasData) { plotMass_vbf.AddHist1D(hMass_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==izmm) continue;
    if(snamev[isam].Contains("htt_gf")) {
      hMass_vbfv[isam]->Scale(10.);
      plotMass_vbf.AddToStack(hMass_vbfv[isam],"#bf{(10x)}"+samplev[isam]->label,samplev[isam]->color,1,1,3);
    }
    else if(snamev[isam].Contains("htt_vbf")) {
      hMass_vbfv[isam]->Scale(10.);
      plotMass_vbf.AddToStack(hMass_vbfv[isam],"",samplev[isam]->color,0,1,2);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_vbf.AddToStack(hMass_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_vbf.SetYRange(0,2.1*(plotMass_vbf.GetStack()->GetMaximum()));
  plotMass_vbf.SetLegend(0.7,0.65,0.9,0.9);
  if(lumi>0) plotMass_vbf.AddTextBox(0.55,0.75,0.85,0.895,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]); //(0.6,0.754,0.93,0.9);
  // plotMass_vbf.AddTextBox("#bf{m_{H}=120}",0.737,0.537,0.837,0.637,0);
  plotMass_vbf.SetLegend(0.55,0.4,0.9,0.65);
  plotMass_vbf.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // no b-tag
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMass_nobv[0]->GetBinWidth(1));
  CPlot plotMass_nob("vismass_class_nob","","#bf{m_{vis} [GeV]}",ylabel);
  if(hMass_nobv[ifake] && hMass_nobv[izmm]) hMass_nobv[ifake]->Add(hMass_nobv[izmm]);
  if(hasData) { plotMass_nob.AddHist1D(hMass_nobv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==izmm) continue;
    if(snamev[isam].Contains("htt_gg")) plotMass_nob.AddToStack(hMass_nobv[isam],"",samplev[isam]->color,0,1,2);
    else if(snamev[isam].Contains("htt_bb")) plotMass_nob.AddToStack(hMass_nobv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_nob.AddToStack(hMass_nobv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_nob.SetLegend(0.7,0.65,0.9,0.9);
  if(lumi>0) plotMass_nob.AddTextBox(0.55,0.75,0.85,0.895,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]); //(0.6,0.754,0.93,0.9);
  // plotMass_nob.AddTextBox("#bf{m_{A}=120, tan#beta=20}",0.737,0.537,0.897,0.637,0);
  plotMass_nob.SetLegend(0.55,0.4,0.9,0.65);
  plotMass_nob.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // b-tag
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"#bf{Events / %.1f GeV}",hMass_bv[0]->GetBinWidth(1));
  CPlot plotMass_b("vismass_class_b","","#bf{m_{vis} [GeV]}",ylabel);
  if(hMass_bv[ifake] && hMass_bv[izmm]) hMass_bv[ifake]->Add(hMass_bv[izmm]);
  if(hasData) { plotMass_b.AddHist1D(hMass_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==izmm) continue;
    if(snamev[isam].Contains("htt_gg")) plotMass_b.AddToStack(hMass_bv[isam],"",samplev[isam]->color,0,1,2);
    else if(snamev[isam].Contains("htt_bb")) plotMass_b.AddToStack(hMass_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_b.AddToStack(hMass_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b.SetYRange(0,1.6*(plotMass_b.GetStack()->GetMaximum()));
  plotMass_b.SetLegend(0.7,0.65,0.9,0.9);
  if(lumi>0) plotMass_b.AddTextBox(0.21,0.75,0.51,0.9,0,1,-1,3,lumitext[0],lumitext[1],lumitext[2]); //(0.6,0.754,0.93,0.9);
  // plotMass_b.AddTextBox("#bf{m_{A}=120, tan#beta=20}",0.737,0.537,0.897,0.637,0);
  plotMass_b.SetLegend(0.55,0.4,0.9,0.65);
  plotMass_b.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // write hists to limit file
  //----------------------------------------------------------------------------------------
  TString flimitstr(CPlot::sOutDir + "/limit-inputs.root");
  TFile flimit(flimitstr,"recreate");
  flimit.mkdir("emu_X");
  flimit.mkdir("emu_novbf");
  flimit.mkdir("emu_vbf");
  flimit.mkdir("emu_nob");
  flimit.mkdir("emu_b");
    
  // add zmm into the fakes histogram
  if(hMass_iv[ifake] && hMass_iv[izmm]) { 
    hMassL_iv[ifake]     ->Add(hMassL_iv[izmm]);
    hMassL_novbfv[ifake] ->Add(hMassL_novbfv[izmm]);
    hMassL_vbfv[ifake]   ->Add(hMassL_vbfv[izmm]);
    hMassL_nobv[ifake]   ->Add(hMassL_nobv[izmm]);
    hMassL_bv[ifake]     ->Add(hMassL_bv[izmm]);	   
  }

  TString histname;
  for(UInt_t isam=0;isam<samplev.size();isam++) {

    if(isam==izmm) continue;

    if(snamev[isam].Contains("ewk",   TString::kIgnoreCase))   	        histname = "EWK"; 
    else if(snamev[isam].Contains("fakes", TString::kIgnoreCase)) 	histname = "Fakes";
    else if(snamev[isam].Contains("ttbar", TString::kIgnoreCase)) 	histname = "ttbar";
    else if(snamev[isam].Contains("Ztt",   TString::kIgnoreCase))   	histname = "Ztt"; 
    else if(snamev[isam].Contains("data",  TString::kIgnoreCase))  	histname = "data_obs";
    else if(snamev[isam].Contains("htt_")) continue;
    else if( (snamev[isam].Contains("gg_")) ||
	     (snamev[isam].Contains("bb_")) ||
	     (snamev[isam].Contains("gf_")) ||
	     (snamev[isam].Contains("vbf_"))  )                         histname = "Higgs_" + snamev[isam];
    else { cout << "error! name not found" << endl; assert(0); }
    
    // cout << "writing: " << snamev[isam] << " " << histname << endl;
    flimit.cd("emu_X");
    hMassL_iv[isam]->SetName(histname);
    hMassL_iv[isam]->Write();
    flimit.cd("emu_novbf");
    hMassL_novbfv[isam]->SetName(histname);
    hMassL_novbfv[isam]->Write();
    flimit.cd("emu_vbf");
    hMassL_vbfv[isam]->SetName(histname);
    hMassL_vbfv[isam]->Write();
    flimit.cd("emu_nob");
    hMassL_nobv[isam]->SetName(histname);
    hMassL_nobv[isam]->Write();
    flimit.cd("emu_b");
    hMassL_bv[isam]->SetName(histname);
    hMassL_bv[isam]->Write();

    // // write the energy scale up/down plots into the central-value file
    // if(outputDir.Contains("escale-down") || outputDir.Contains("escale-up") ||
    //    // outputDir.Contains("jet-down")    || outputDir.Contains("jet-up")    ||
    //    outputDir.Contains("eres-down")   || outputDir.Contains("eres-up")) {

    //   TString escale_suffix;
    //   if(outputDir.Contains("escale-down"))      escale_suffix = "_CMS_scale_eDown";
    //   else if(outputDir.Contains("escale-up"))   escale_suffix = "_CMS_scale_eUp";
    //   else if(outputDir.Contains("eres-down"))   escale_suffix = "_CMS_res_eDown";
    //   else if(outputDir.Contains("eres-up"))     escale_suffix = "_CMS_res_eUp";
    //   // else if(outputDir.Contains("jet-down"))    escale_suffix = "_CMS_scale_jetDown";
    //   // else if(outputDir.Contains("jet-up"))      escale_suffix = "_CMS_scale_jetUp";

    //   TFile fbase("/scratch/dkralph/htt/selections/1084/plots/limit-inputs.root","update");
    //   fbase.cd("emu_X");
    //   hMassL_iv[isam]->SetName(histname+escale_suffix);
    //   hMassL_iv[isam]->Write();
    //   fbase.cd("emu_novbf");
    //   hMassL_novbfv[isam]->SetName(histname+escale_suffix);
    //   hMassL_novbfv[isam]->Write();
    //   fbase.cd("emu_vbf");
    //   hMassL_vbfv[isam]->SetName(histname+escale_suffix);
    //   hMassL_vbfv[isam]->Write();
    //   fbase.cd("emu_nob");
    //   hMassL_nobv[isam]->SetName(histname+escale_suffix);
    //   hMassL_nobv[isam]->Write();
    //   fbase.cd("emu_b");
    //   hMassL_bv[isam]->SetName(histname+escale_suffix);
    //   hMassL_bv[isam]->Write();
    //   fbase.Close();
    // }

  }
  flimit.Close();
    

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  Double_t nTotal=0,       nTotalVar=0;    // background MC totals
  Double_t nTotal_i=0,     nTotalVar_i=0;
  Double_t nTotal_novbf=0, nTotalVar_novbf=0;
  Double_t nTotal_vbf=0,   nTotalVar_vbf=0;
  Double_t nTotal_nob=0,   nTotalVar_nob=0;
  Double_t nTotal_b=0,     nTotalVar_b=0;

  ofstream yieldfile;
  yieldfile.open(CPlot::sOutDir+"/yields.txt");
  
  yieldfile << setw(33) << "lepton sele." <<setw(20) << "inclusive" << setw(25) << "no vbf" << setw(19)
       << "vbf" << setw(17) << "no b-tag" << setw(17) << "b-tag" << endl;

  if(hMass_iv[ifake] && hMass_iv[izmm]) { // add zmm into the fakes
    nSelv[ifake]       += nSelv[izmm];       nSelVarv[ifake]       += nSelVarv[izmm];
    nSel_iv[ifake]     += nSel_iv[izmm];     nSelVar_iv[ifake]     += nSelVar_iv[izmm];
    nSel_novbfv[ifake] += nSel_novbfv[izmm]; nSelVar_novbfv[ifake] += nSelVar_novbfv[izmm];
    nSel_vbfv[ifake]   += nSel_vbfv[izmm];   nSelVar_vbfv[ifake]   += nSelVar_vbfv[izmm];
    nSel_nobv[ifake]   += nSel_nobv[izmm];   nSelVar_nobv[ifake]   += nSelVar_nobv[izmm];
    nSel_bv[ifake]     += nSel_bv[izmm];     nSelVar_bv[ifake]     += nSelVar_bv[izmm];
  }

  if(samplev.size()>1) {
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      if(isam==izmm) continue;

      yieldfile << setw(15) << snamev[isam];
      yieldfile << setw(10) << setprecision(3) << fixed << nSelv[isam]        << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[isam]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[isam]      << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[isam]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_novbfv[isam]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_novbfv[isam]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[isam]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[isam]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_nobv[isam]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_nobv[isam]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[isam]      << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_bv[isam]);
      yieldfile << endl;

      if(snamev[isam].Contains("_mssm_"))  continue; // don't add higgs samples to total
      if(snamev[isam].Contains("_sm_"))    continue;

      nTotal    	+= nSelv[isam];
      nTotalVar 	+= nSelVarv[isam];
      nTotal_i      	+= nSel_iv[isam];
      nTotalVar_i   	+= nSelVar_iv[isam];
      nTotal_novbf   	+= nSel_novbfv[isam];
      nTotalVar_novbf 	+= nSelVar_novbfv[isam];
      nTotal_vbf      	+= nSel_vbfv[isam];
      nTotalVar_vbf   	+= nSelVar_vbfv[isam];
      nTotal_nob      	+= nSel_nobv[isam];
      nTotalVar_nob   	+= nSelVar_nobv[isam];
      nTotal_b      	+= nSel_bv[isam];
      nTotalVar_b       += nSelVar_bv[isam];
    }
    yieldfile << endl;
    yieldfile << setw(15) << "bkg MC";
    yieldfile << setw(10) << setprecision(3) << fixed << nTotal       << " +/- " << sqrt(nTotalVar);
    yieldfile << setw(10) << setprecision(3) << fixed << nTotal_i     << " +/- " << sqrt(nTotalVar_i);
    yieldfile << setw(10) << setprecision(3) << fixed << nTotal_novbf << " +/- " << sqrt(nTotalVar_novbf);
    yieldfile << setw(10) << setprecision(3) << fixed << nTotal_vbf   << " +/- " << sqrt(nTotalVar_vbf);
    yieldfile << setw(10) << setprecision(3) << fixed << nTotal_nob   << " +/- " << sqrt(nTotalVar_nob);
    yieldfile << setw(10) << setprecision(3) << fixed << nTotal_b     << " +/- " << sqrt(nTotalVar_b);
    yieldfile << endl;
  }
  if(hasData) {
    yieldfile << setw(15) << "Data";
    yieldfile << setw(10) << setprecision(3) << fixed << nSelv[0]        << " +/- " << sqrt(nSelVarv[0]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[0]      << " +/- " << sqrt(nSelVar_iv[0]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_novbfv[0]  << " +/- " << sqrt(nSelVar_novbfv[0]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[0]    << " +/- " << sqrt(nSelVar_vbfv[0]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_nobv[0]    << " +/- " << sqrt(nSelVar_nobv[0]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[0]      << " +/- " << sqrt(nSelVar_bv[0]);
    yieldfile << endl;
  }

  // cat the yields to stdout
  TString cmd("cat "+CPlot::sOutDir+"/yields.txt");
  system(cmd.Data());

  // write out the acceptances
  yieldfile << endl << endl;
  yieldfile << setw(25) << " " << setw(15) << "initial events";
  yieldfile << setw(20) << "lepton sele." <<setw(20) << "inclusive" << setw(20) << "no vbf" << setw(20)
	    << "vbf" << setw(20) << "no b-tag" << setw(20) << "b-tag" << endl;
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if(!(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_"))) continue;
    if(snamev[isam].Contains("htt_")) continue;
    yieldfile << setw(25) << snamev[isam] << "  " <<
      setw(15) << nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSelv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_iv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_novbfv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_vbfv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_nobv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_bv[isam]/nUnskEventsv[isam] << endl;
  }
  yieldfile.close();

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
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lepdeta.png\"><img src=\"plots/lepdeta.png\" alt=\"plots/lepdeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
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
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetavis2.png\"><img src=\"plots/pzetavis2.png\" alt=\"plots/pzetavis2.png\" width=\"100%\"></a></td>" << endl;

  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetavar.png\"><img src=\"plots/pzetavar.png\" alt=\"plots/pzetavar.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lepdeta.png\"><img src=\"plots/lepdeta.png\" alt=\"plots/lepdeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/rawProjVar.png\"><img src=\"plots/rawProjVar.png\" alt=\"plots/rawProjVar.png\" width=\"100%\"></a></td>" << endl;
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
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nvertices_raw.png\"><img src=\"plots/nvertices_raw.png\" alt=\"plots/nvertices_raw.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;   
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_lept.png\"><img src=\"plots/vismass_lept.png\" alt=\"plots/vismass_lept.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_incl.png\"><img src=\"plots/vismass_incl.png\" alt=\"plots/vismass_incl.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_novbf.png\"><img src=\"plots/vismass_class_novbf.png\" alt=\"plots/vismass_class_novbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_vbf.png\"><img src=\"plots/vismass_class_vbf.png\" alt=\"plots/vismass_class_vbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_nob.png\"><img src=\"plots/vismass_class_nob.png\" alt=\"plots/vismass_class_nob.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_b.png\"><img src=\"plots/vismass_class_b.png\" alt=\"plots/vismass_class_b.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_lept-high.png\"><img src=\"plots/vismass_lept-high.png\" alt=\"plots/vismass_lept-high.png\" width=\"100%\"></a></td>" << endl;
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
