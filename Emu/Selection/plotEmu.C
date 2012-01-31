#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TRegexp.h>                // ROOT regexp class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "MitHtt/Common/CPlot.hh"          // helper class for plots
#include "MitHtt/Common/MitStyleRemix.hh"  // style settings for drawing
#include "MitHtt/Common/MyTools.hh"        // miscellaneous helper functions
#include "MitHtt/Common/CSample.hh"        // helper class for organizing input ntuple files

// define structure for output ntuple
#include "EmuData.hh"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

// rescale histograms with up/down shifts
Double_t rescale(TH1F* hist, TString catname, TString hname, TFile* file);

// generate web page
void makeHTML(const TString outDir);

enum { kCenter, kDown, kUp };

//=== MAIN MACRO =================================================================================================

// run first with ecorr=0 to get central histograms, then with ecorr=1,2 to get the up,down scale variations

void plotEmu(const TString  conf,         // input file
             const TString  ntupleDir,    // directory of input ntuples
	     const TString  outputDir,    // output directory
             const TString  format,       // plot file format
	     const Double_t lumi,         // luminosity (pb^-1)
             const Int_t ecorr=0          // energy scale shift: 0=center, 1=down, 2=up
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

  // bins for post-fit histograms
  const UInt_t nbinsPF_i = 41;
  const UInt_t nbinsPF_novbf = 29;
  const UInt_t nbinsPF_boost = 11;
  const UInt_t nbinsPF_vbf = 11;
  const UInt_t nbinsPF_nob = 30;
  const UInt_t nbinsPF_b = 14;

  Double_t massPFEdges_i[nbinsPF_i]         = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,300.,325.,350.,400.,450.,500.,550.,600.,650.,750.,1000.};
  Double_t massPFEdges_novbf[nbinsPF_novbf] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,275.,300.,350.};
  Double_t massPFEdges_boost[nbinsPF_boost] = {0.,25.,50.,75.,100.,125.,150.,175.,200.,250.,350.};
  Double_t massPFEdges_vbf[nbinsPF_vbf]     = {0.,25.,50.,75.,100.,125.,150.,175.,200.,250.,350.};
  Double_t massPFEdges_nob[nbinsPF_nob]     = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,220.,240.,260.,280.,300.,350.,400.,450.,500.};
  Double_t massPFEdges_b[nbinsPF_b]         = {0., 25., 50., 75., 100., 150., 175., 200., 225., 250., 300., 350., 400., 500.};

  // bins for limit input histograms
  const UInt_t massLHSNbins_SM = 30;
  const UInt_t massLHSNbins_MSSM = 35;
  const UInt_t massLLSNbins_SM = 16;
  const UInt_t massLLSNbins_MSSM = 21;

  Double_t massLEdges_i[massLHSNbins_SM]     = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,220.,240.,260.,280.,300.,350.,400.,450.,500.};
  Double_t massLEdges_novbf[massLHSNbins_SM] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,220.,240.,260.,280.,300.,350.,400.,450.,500.};
  Double_t massLEdges_nob[massLHSNbins_MSSM] = {0.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,220.,240.,260.,280.,300.,350.,400.,450.,500.,600.,700.,800.,900.,1000.};
  Double_t massLEdges_boost[massLLSNbins_SM] = {0.,20.,40.,60.,80.,100.,120.,140.,170.,200.,250.,300.,350.,400.,450.,500.};
  Double_t massLEdges_vbf[massLLSNbins_SM]   = {0.,20.,40.,60.,80.,100.,120.,140.,170.,200.,250.,300.,350.,400.,450.,500.};
  Double_t massLEdges_b[massLLSNbins_MSSM]   = {0.,20.,40.,60.,80.,100.,120.,140.,170.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.,900.,1000.};

  // get indices of samples
  UInt_t ifake=9999, izmm=9999, iztt=9999, ittbar=9999, imssm_gg=9999, imssm_bb=9999, ism_vbf=9999, ism_gf=9999, ism_vtth=9999, issfake=9999, iemb=9999, iewk=9999;
  for(UInt_t isam=0; isam<snamev.size(); isam++) {
    if(snamev[isam].Contains("fakes") && !snamev[isam].Contains("ss-fakes")) ifake = isam;
    if(snamev[isam].Contains("ss-fakes")) issfake = isam;
    if(snamev[isam].Contains("zmm"))   izmm  = isam;
    if(snamev[isam].Contains("ztt"))   iztt = isam;
    if(snamev[isam].Contains("emb"))   iemb = isam;
    if(snamev[isam].Contains("ttbar")) ittbar = isam;
    if(snamev[isam].Contains("ewk"))   iewk = isam;
    if(snamev[isam].Contains("htt_gg_mssm")) imssm_gg=isam;
    if(snamev[isam].Contains("htt_bb_mssm")) imssm_bb=isam;
    if(snamev[isam].Contains("htt_vbf_sm")) ism_vbf=isam;
    if(snamev[isam].Contains("htt_gf_sm")) ism_gf=isam;
    if(snamev[isam].Contains("htt_vtth_sm")) ism_vtth=isam;
  }

  CPlot::sOutDir = outputDir + TString("/plots");
  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };  // final state type

  const Double_t kssFakeWgt = 1.2;  // ratio of fake-rate fakes to same-sign fakes

  //
  // scale factors for each category (used only for madgraph ztt)
  //
  const Double_t kCat_novbf = 1.00;
  const Double_t kCat_boost = 0.92;
  const Double_t kCat_vbf   = 1.20;
  const Double_t kCat_nob   = 1.00;
  const Double_t kCat_b     = 0.99;

  // vbf cut values
  const Double_t mjjMin   = 400;
  const Double_t dEtaMin  = 4.0;


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
  vector<TH1F*> hProjMetv, hProjVisv, hProjVarv,hRawProjVarv;
  vector<TH1F*> hNjetsv, hNbjetsv;
  vector<TH1F*> hBdiscrv, hBdiscr_vbfv;
  vector<TH1F*> hDPhiv, hMtv, hPtv, hLepDEtav;
  vector<TH1F*> hMetDPhiv;
  vector<TH1F*> hPt1v, hEta1v, hPhi1v;   	// leading lepton
  vector<TH1F*> hPt2v, hEta2v, hPhi2v;   	// trailing lepton
  vector<TH1F*> hPtMuv, hEtaMuv, hPhiMuv;  	// muon
  vector<TH1F*> hPtElev, hEtaElev, hPhiElev; 	// electron
  vector<TH1F*> hMtElev, hMtMuv;                // mT
  vector<TH1F*> hJetPt1v, hJetEta1v; 		// leading jet
  vector<TH1F*> hJetPt2v, hJetEta2v; 		// second pt jet
  vector<TH1F*> hJetDPhiv;           		// dphi of first two jets
  vector<TH1F*> hMjjv, hDEtav, hEtaProdv;       // kinematics of first two jets
  vector<TH1F*> hBJetPtv, hBJetEtav, hBJetPhiv; // leading b-jet (can be same as leading two jets)
  vector<TH1F*> hNPVv, hNPVrawv;                // primary vertexes
  vector<TH1F*> hVBFJetPt1v, hVBFJetEta1v;      // leading jet
  vector<TH1F*> hVBFJetPt2v, hVBFJetEta2v;      // second pt jet
  vector<TH1F*> hBoostJetPtv, hBoostJetEtav;    // leading jet
  vector<TH1F*> hMassv, hMassVisv, hMassLv, hMassHighv, hMassVisHighv;    // *L: plots for limits

  // inclusive
  vector<TH1F*> hMass_iv;
  vector<TH1F*> hMassVis_iv;
  vector<TH1F*> hMassL_iv;
  vector<TH1F*> hMassPF_iv;

  // Class novbf: one or fewer pt 30 jets, or exactly two pt 30 jets and !(vbf cuts)
  vector<TH1F*> hMass_novbfv;
  vector<TH1F*> hMassVis_novbfv;
  vector<TH1F*> hMassL_novbfv;
  vector<TH1F*> hMassPF_novbfv;

  // Class boosted: one or fewer pt 30 jets + boost cuts
  vector<TH1F*> hMass_boostv;
  vector<TH1F*> hMassVis_boostv;
  vector<TH1F*> hMassL_boostv;
  vector<TH1F*> hMassPF_boostv;

  // Class vbf: exactly two pt 30 jets and vbf cuts
  vector<TH1F*> hMass_vbfv;
  vector<TH1F*> hMassVis_vbfv;
  vector<TH1F*> hMassL_vbfv;
  vector<TH1F*> hMassPF_vbfv;

  // Class nob: two or fewer pt 30 jets, no b-tagged jet
  vector<TH1F*> hMass_nobv;
  vector<TH1F*> hMassVis_nobv;
  vector<TH1F*> hMassL_nobv;
  vector<TH1F*> hMassPF_nobv;

  // Class b: two or fewer pt 30 jets, at least one b-tagged jet
  vector<TH1F*> hMass_bv;
  vector<TH1F*> hMassVis_bv;
  vector<TH1F*> hMassL_bv;
  vector<TH1F*> hMassPF_bv;

  vector<Double_t> nSelv,       nSelVarv;
  vector<Double_t> nSel_iv,     nSelVar_iv;
  vector<Double_t> nSel_novbfv, nSelVar_novbfv;
  vector<Double_t> nSel_boostv, nSelVar_boostv;
  vector<Double_t> nSel_vbfv,   nSelVar_vbfv;
  vector<Double_t> nSel_nobv,   nSelVar_nobv;
  vector<Double_t> nSel_bv,     nSelVar_bv;
  
  vector<Double_t> nUnskEventsv;

  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    // after lepton selection
    sprintf(hname,"hMet_%i",isam);     hMetv.push_back(new TH1F(hname,"",30,0,150));        hMetv[isam]->Sumw2();
    sprintf(hname,"hMetRaw_%i",isam);  hMetRawv.push_back(new TH1F(hname,"",30,0,150));     hMetRawv[isam]->Sumw2();
    sprintf(hname,"hProjMet_%i",isam);       hProjMetv.push_back(new TH1F(hname,"",30,-100,120));     hProjMetv[isam]->Sumw2();
    sprintf(hname,"hProjVis_%i",isam);       hProjVisv.push_back(new TH1F(hname,"",20,0,50));         hProjVisv[isam]->Sumw2();
    sprintf(hname,"hProjVar_%i",isam);       hProjVarv.push_back(new TH1F(hname,"",30,-230,100));     hProjVarv[isam]->Sumw2();
    sprintf(hname,"hRawProjVar_%i",isam);    hRawProjVarv.push_back(new TH1F(hname,"",30,-230,100));  hRawProjVarv[isam]->Sumw2();
    sprintf(hname,"hNjets_%i",isam);         hNjetsv.push_back(new TH1F(hname,"",11,-0.5,10.5));      hNjetsv[isam]->Sumw2();
    sprintf(hname,"hNbjets_%i",isam);        hNbjetsv.push_back(new TH1F(hname,"",6,-0.5,5.5));       hNbjetsv[isam]->Sumw2();
    sprintf(hname,"hBdiscr_%i",isam);        hBdiscrv.push_back(new TH1F(hname,"",50,0,15));          hBdiscrv[isam]->Sumw2();
    sprintf(hname,"hBdiscr_vbf_%i",isam);    hBdiscr_vbfv.push_back(new TH1F(hname,"",30,0,15));      hBdiscr_vbfv[isam]->Sumw2();
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
    sprintf(hname,"hPtMu_%i",isam);    hPtMuv.push_back(new TH1F(hname,"",20,0,100));       hPtMuv[isam]->Sumw2();
    sprintf(hname,"hEtaMu_%i",isam);   hEtaMuv.push_back(new TH1F(hname,"",20,-3,3));       hEtaMuv[isam]->Sumw2();
    sprintf(hname,"hPhiMu_%i",isam);   hPhiMuv.push_back(new TH1F(hname,"",20,-3.2,3.2));   hPhiMuv[isam]->Sumw2();  
    sprintf(hname,"hPtEle_%i",isam);   hPtElev.push_back(new TH1F(hname,"",20,0,100));      hPtElev[isam]->Sumw2();
    sprintf(hname,"hEtaEle_%i",isam);  hEtaElev.push_back(new TH1F(hname,"",20,-3,3));      hEtaElev[isam]->Sumw2();
    sprintf(hname,"hPhiEle_%i",isam);  hPhiElev.push_back(new TH1F(hname,"",20,-3.2,3.2));  hPhiElev[isam]->Sumw2();
    sprintf(hname,"hMtEle_%i",isam);   hMtElev.push_back(new TH1F(hname,"",30,0,150));      hMtElev[isam]->Sumw2();
    sprintf(hname,"hMtMu_%i",isam);    hMtMuv.push_back(new TH1F(hname,"",30,0,150));       hMtMuv[isam]->Sumw2();
    sprintf(hname,"hJetPt1_%i",isam);  hJetPt1v.push_back(new TH1F(hname,"",30,0,300));     hJetPt1v[isam]->Sumw2();
    sprintf(hname,"hJetEta1_%i",isam); hJetEta1v.push_back(new TH1F(hname,"",20,-5,5));     hJetEta1v[isam]->Sumw2();
    sprintf(hname,"hJetPt2_%i",isam);  hJetPt2v.push_back(new TH1F(hname,"",30,0,300));     hJetPt2v[isam]->Sumw2();
    sprintf(hname,"hJetEta2_%i",isam); hJetEta2v.push_back(new TH1F(hname,"",20,-5,5));     hJetEta2v[isam]->Sumw2();
    sprintf(hname,"hJetDPhi_%i",isam); hJetDPhiv.push_back(new TH1F(hname,"",30,0,180));    hJetDPhiv[isam]->Sumw2();
    sprintf(hname,"hMjj_%i",isam);     hMjjv.push_back(new TH1F(hname,"",20,0,1000));       hMjjv[isam]->Sumw2();
    sprintf(hname,"hDEta_%i",isam);    hDEtav.push_back(new TH1F(hname,"",20,0,8));         hDEtav[isam]->Sumw2();
    sprintf(hname,"hEtaProd_%i",isam); hEtaProdv.push_back(new TH1F(hname,"",30,-7.5,7.5)); hEtaProdv[isam]->Sumw2();
    sprintf(hname,"hBJetPt_%i",isam);  hBJetPtv.push_back(new TH1F(hname,"",30,0,150));     hBJetPtv[isam]->Sumw2();
    sprintf(hname,"hBJetEta_%i",isam); hBJetEtav.push_back(new TH1F(hname,"",20,-3,3));     hBJetEtav[isam]->Sumw2();
    sprintf(hname,"hBJetPhi_%i",isam); hBJetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hBJetPhiv[isam]->Sumw2();
    sprintf(hname,"hNPV_%i",isam);     hNPVv.push_back(new TH1F(hname,"",20,-0.5,19.5));    hNPVv[isam]->Sumw2();
    sprintf(hname,"hNPVraw_%i",isam);  hNPVrawv.push_back(new TH1F(hname,"",20,-0.5,19.5)); hNPVrawv[isam]->Sumw2();
    sprintf(hname,"hVBFJetPt1_%i",isam);     hVBFJetPt1v.push_back(new TH1F(hname,"",27,30,300));     hVBFJetPt1v[isam]->Sumw2();
    sprintf(hname,"hVBFJetEta1_%i",isam);    hVBFJetEta1v.push_back(new TH1F(hname,"",18,-4.5,4.5));  hVBFJetEta1v[isam]->Sumw2();
    sprintf(hname,"hVBFJetPt2_%i",isam);     hVBFJetPt2v.push_back(new TH1F(hname,"",27,30,300));     hVBFJetPt2v[isam]->Sumw2();
    sprintf(hname,"hVBFJetEta2_%i",isam);    hVBFJetEta2v.push_back(new TH1F(hname,"",18,-4.5,4.5));  hVBFJetEta2v[isam]->Sumw2();
    sprintf(hname,"hBoostJetPt_%i",isam);    hBoostJetPtv.push_back(new TH1F(hname,"",27,30,300));    hBoostJetPtv[isam]->Sumw2();
    sprintf(hname,"hBoostJetEta_%i",isam);   hBoostJetEtav.push_back(new TH1F(hname,"",18,-4.5,4.5)); hBoostJetEtav[isam]->Sumw2();

    // lepton selection
    sprintf(hname,"hMass_%i",isam);          hMassv.push_back(new TH1F(hname,"",40,0,400));                                  hMassv[isam]->Sumw2();
    sprintf(hname,"hMassVis_%i",isam);       hMassVisv.push_back(new TH1F(hname,"",40,0,400));                               hMassVisv[isam]->Sumw2();
    sprintf(hname,"hMassHigh_%i",isam);      hMassHighv.push_back(new TH1F(hname,"",100,0,2000));                            hMassHighv[isam]->Sumw2();
    sprintf(hname,"hMassVisHigh_%i",isam);   hMassVisHighv.push_back(new TH1F(hname,"",100,0,2000));                         hMassVisHighv[isam]->Sumw2();
    sprintf(hname,"hMassL_%i",isam);         hMassLv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_i));           hMassLv[isam]->Sumw2();
    sprintf(hname,"hMassVpmiss_%i",isam);    hMassVpmissv.push_back(new TH2F(hname,"",75,0,750,50,-100,100));                hMassVpmissv[isam]->Sumw2();
    // inclusive
    sprintf(hname,"hMass_i_%i",isam);        hMass_iv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_i));          hMass_iv[isam]->Sumw2();
    sprintf(hname,"hMassPF_i_%i",isam);      hMassPF_iv.push_back(new TH1F(hname,"",nbinsPF_i-1,massPFEdges_i));             hMassPF_iv[isam]->Sumw2();
    sprintf(hname,"hMassVis_i_%i",isam);     hMassVis_iv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_i));       hMassVis_iv[isam]->Sumw2();
    sprintf(hname,"hMassL_i_%i",isam);       hMassL_iv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_i));         hMassL_iv[isam]->Sumw2();
    // no vbf
    sprintf(hname,"hMass_novbf_%i",isam);    hMass_novbfv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_i));      hMass_novbfv[isam]->Sumw2();
    sprintf(hname,"hMassPF_novbf_%i",isam);  hMassPF_novbfv.push_back(new TH1F(hname,"",nbinsPF_novbf-1,massPFEdges_novbf)); hMassPF_novbfv[isam]->Sumw2();
    sprintf(hname,"hMassVis_novbf_%i",isam); hMassVis_novbfv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_i));   hMassVis_novbfv[isam]->Sumw2();
    sprintf(hname,"hMassL_novbf_%i",isam);   hMassL_novbfv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_novbf)); hMassL_novbfv[isam]->Sumw2();
    // boosted
    sprintf(hname,"hMass_boost_%i",isam);    hMass_boostv.push_back(new TH1F(hname,"",massLLSNbins_SM-1,massLEdges_vbf));    hMass_boostv[isam]->Sumw2();
    sprintf(hname,"hMassPF_boost_%i",isam);  hMassPF_boostv.push_back(new TH1F(hname,"",nbinsPF_boost-1,massPFEdges_boost)); hMassPF_boostv[isam]->Sumw2();
    sprintf(hname,"hMassVis_boost_%i",isam); hMassVis_boostv.push_back(new TH1F(hname,"",massLLSNbins_SM-1,massLEdges_vbf)); hMassVis_boostv[isam]->Sumw2();
    sprintf(hname,"hMassL_boost_%i",isam);   hMassL_boostv.push_back(new TH1F(hname,"",massLLSNbins_SM-1,massLEdges_boost)); hMassL_boostv[isam]->Sumw2();
    // vbf
    sprintf(hname,"hMass_vbf_%i",isam);      hMass_vbfv.push_back(new TH1F(hname,"",massLLSNbins_SM-1,massLEdges_vbf));      hMass_vbfv[isam]->Sumw2();
    sprintf(hname,"hMassPF_vbf_%i",isam);    hMassPF_vbfv.push_back(new TH1F(hname,"",nbinsPF_vbf-1,massPFEdges_vbf));       hMassPF_vbfv[isam]->Sumw2();
    sprintf(hname,"hMassVis_vbf_%i",isam);   hMassVis_vbfv.push_back(new TH1F(hname,"",massLLSNbins_SM-1,massLEdges_vbf));   hMassVis_vbfv[isam]->Sumw2();
    sprintf(hname,"hMassL_vbf_%i",isam);     hMassL_vbfv.push_back(new TH1F(hname,"",massLLSNbins_SM-1,massLEdges_vbf));     hMassL_vbfv[isam]->Sumw2();
    // no b-tag
    sprintf(hname,"hMass_nob_%i",isam);      hMass_nobv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_i));        hMass_nobv[isam]->Sumw2();
    sprintf(hname,"hMassPF_nob_%i",isam);    hMassPF_nobv.push_back(new TH1F(hname,"",nbinsPF_nob-1,massPFEdges_nob));       hMassPF_nobv[isam]->Sumw2();
    sprintf(hname,"hMassVis_nob_%i",isam);   hMassVis_nobv.push_back(new TH1F(hname,"",massLHSNbins_SM-1,massLEdges_i));     hMassVis_nobv[isam]->Sumw2();
    sprintf(hname,"hMassL_nob_%i",isam);     hMassL_nobv.push_back(new TH1F(hname,"",massLHSNbins_MSSM-1,massLEdges_nob));   hMassL_nobv[isam]->Sumw2();
    // b-tag
    sprintf(hname,"hMass_b_%i",isam);        hMass_bv.push_back(new TH1F(hname,"",massLLSNbins_SM-1,massLEdges_vbf));        hMass_bv[isam]->Sumw2();
    sprintf(hname,"hMassPF_b_%i",isam);      hMassPF_bv.push_back(new TH1F(hname,"",nbinsPF_b-1,massPFEdges_b));             hMassPF_bv[isam]->Sumw2();
    sprintf(hname,"hMassVis_b_%i",isam);     hMassVis_bv.push_back(new TH1F(hname,"",massLLSNbins_SM-1,massLEdges_vbf));     hMassVis_bv[isam]->Sumw2();
    sprintf(hname,"hMassL_b_%i",isam);	     hMassL_bv.push_back(new TH1F(hname,"",massLLSNbins_MSSM-1,massLEdges_b));       hMassL_bv[isam]->Sumw2();

    nSelv.push_back(0);       nSelVarv.push_back(0);
    nSel_iv.push_back(0);     nSelVar_iv.push_back(0);
    nSel_novbfv.push_back(0); nSelVar_novbfv.push_back(0);
    nSel_boostv.push_back(0); nSelVar_boostv.push_back(0);
    nSel_vbfv.push_back(0);   nSelVar_vbfv.push_back(0);
    nSel_nobv.push_back(0);   nSelVar_nobv.push_back(0);
    nSel_bv.push_back(0);     nSelVar_bv.push_back(0);

    nUnskEventsv.push_back(0);
    
  }

  EmuData data;
  Double_t rawMet=0,rawprojvar=0,npuWgt=1;
  UInt_t npt20jets;
  const UInt_t kMaxPt20Jets=35;
  TArrayF *btagArray = new TArrayF;
  TArrayF *jptArray = new TArrayF;
  TArrayF *jetaArray = new TArrayF;
  btagArray->Set(kMaxPt20Jets);
  jptArray->Set(kMaxPt20Jets);
  jetaArray->Set(kMaxPt20Jets);
  
  //
  // Access samples and fill histograms
  //  
  TFile *unfFile1   = TFile::Open("data/unfold/v8/Unfold2D_1.root"); assert(unfFile1->IsOpen());
  TH2F  *unfWeight1 = (TH2F*) unfFile1->FindObjectAny("UnfoldDen1");
  TFile *unfFile2   = TFile::Open("data/unfold/v8/Unfold2D_2.root"); assert(unfFile2->IsOpen());
  TH2F  *unfWeight2 = (TH2F*) unfFile2->FindObjectAny("UnfoldDen2");


  TFile *infile=0;
  TTree *eventTree=0;

  // Flat ntuple with information about selected events
  TFile   evtfile("emu-vbf.root","RECREATE");
  TTree   evttree("Events","Events");
  UInt_t  run;         evttree.Branch("run"       ,&run);
  UInt_t  lumisec;     evttree.Branch("lumisec"   ,&lumisec);
  UInt_t  event;       evttree.Branch("event"     ,&event);
  Float_t ptele;       evttree.Branch("ptele"     ,&ptele);
  Float_t etaele;      evttree.Branch("etaele"    ,&etaele);
  Float_t phiele;      evttree.Branch("phiele"    ,&phiele);
  Float_t relisoele;   evttree.Branch("relisoele"    ,&relisoele);
  Float_t ptmu;        evttree.Branch("ptmu"      ,&ptmu);
  Float_t etamu;       evttree.Branch("etamu"     ,&etamu);
  Float_t phimu;       evttree.Branch("phimu"     ,&phimu);
  Float_t relisomu;    evttree.Branch("relisomu"     ,&relisomu);
  Float_t pzeta;       evttree.Branch("pzeta"     ,&pzeta);
  Int_t   njets;       evttree.Branch("njets"     ,&njets);
  Int_t   npv;         evttree.Branch("npv"       ,&npv);
  Float_t svfmass;     evttree.Branch("svfmass"   ,&svfmass);
  Float_t vismass;     evttree.Branch("vismass"   ,&vismass);

  // event lists
  ofstream evtlistfile_novbf, evtlistfile_boost, evtlistfile_vbf, evtlistfile_nob, evtlistfile_b;
  evtlistfile_novbf.open("data/eventlist_novbf.txt");
  evtlistfile_boost.open("data/eventlist_boost.txt");
  evtlistfile_vbf.open("data/eventlist_vbf.txt");
  evtlistfile_nob.open("data/eventlist_nob.txt");
  evtlistfile_b.open("data/eventlist_b.txt");

  evtlistfile_novbf << "#run" << setw(13) << "lumi" << setw(15) << "event" << endl;
  evtlistfile_boost << "#run" << setw(13) << "lumi" << setw(15) << "event" << endl;
  evtlistfile_vbf << "#run" << setw(13) << "lumi" << setw(15) << "event" << endl;
  evtlistfile_nob << "#run" << setw(13) << "lumi" << setw(15) << "event" << endl;
  evtlistfile_b << "3run" << setw(13) << "lumi" << setw(15) << "event" << endl;
  
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
    eventTree->SetBranchAddress("jptArray",&jptArray);
    eventTree->SetBranchAddress("jetaArray",&jetaArray);
    eventTree->SetBranchAddress("rawMet",&rawMet);
    eventTree->SetBranchAddress("rawprojvar",&rawprojvar);
    if(isam!=0 && isam!=ifake && isam!=issfake)
      eventTree->SetBranchAddress("npuWgt",&npuWgt);

    assert(samplev[isam]->fnamev.size()>0);
    nUnskEventsv[isam] = lumi*samplev[isam]->xsecv[0];

    TH1D *helectron=0;
    char ehname[100];
    if (ecorr==kCenter) sprintf(ehname,"hecenter");
    else if (ecorr==kDown) sprintf(ehname,"hedown");
   
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);

      Double_t wgt = 1.0;
      if(isam > 0 && isam!=ifake && !snamev[isam].Contains("ss-fakes") && isam!=iemb)
	wgt = data.weight*lumi;
      if(snamev[isam].Contains("ss-fakes")) // same-sign fakes must have weight 1, but normalize them to fake-rate fakes
	wgt = data.weight*(kssFakeWgt);
      else if (isam==ifake) // fake-rate fakes have their weight stored in the ntuple
	wgt = data.weight;
      if(isam==iemb) {
	wgt=data.weight*lumi*0.04201335;        // normalized using madgraph
	Double_t pt1,eta1,pt2,eta2;
	pt1 = data.genlpt1; eta1 = data.genleta1; pt2 = data.genlpt2; eta2 = data.genleta2;

	double weight1 = unfWeight1->GetBinContent(unfWeight1->GetXaxis()->FindBin(eta1),unfWeight1->GetYaxis()->FindBin(pt1));
	double weight2 = unfWeight2->GetBinContent(unfWeight2->GetXaxis()->FindBin(eta2),unfWeight2->GetYaxis()->FindBin(pt2));
	double embUnfWgt=weight1*weight2;
	//wgt*=embUnfWgt;    // unfolding weight is applied in selectEmu
      }

      // get muon, electron kinematics
      Float_t mupt,mueta,muphi,elept,eleta,elephi;
      if(data.state==kMuEle) {
	mupt  = data.lpt1; mueta = data.leta1; muphi = data.lphi1; elept = data.lpt2; eleta = data.leta2; elephi = data.lphi2;
      } else {
	mupt  = data.lpt2; mueta = data.leta2; muphi = data.lphi2; elept = data.lpt1; eleta = data.leta1; elephi = data.lphi1;
      }

      if(fabs(eleta)>2.3) continue;                                 

      Float_t mtele, mtmu;

      mtele = sqrt( 2.0 * (elept*data.met*(1.0-cos(toolbox::deltaPhi(elephi,data.metphi))) )); 
      mtmu  = sqrt( 2.0 * (mupt *data.met*(1.0-cos(toolbox::deltaPhi(muphi,data.metphi))) ));               

      // fill plots after lepton selection
      if(isam!=ifake) {
	hMetv[isam]        ->Fill(data.met,    	 wgt);
	hMetRawv[isam]     ->Fill(rawMet,      	 wgt);
	hProjMetv[isam]    ->Fill(data.pmet,     	 wgt);
	hProjVisv[isam]    ->Fill(data.pvis,     	 wgt);
	hProjVarv[isam]    ->Fill(-0.85*data.pvis + data.pmet,  wgt);
	hRawProjVarv[isam] ->Fill(-rawprojvar,     wgt);
	hNjetsv[isam]      ->Fill(data.njets,      wgt);
	hNbjetsv[isam]     ->Fill(data.nbjets,     wgt);

	assert(npt20jets<kMaxPt20Jets);
	for(UInt_t ib=0;ib<npt20jets;ib++)
	  hBdiscrv[isam]   ->Fill((*btagArray)[ib],  wgt);

	hDPhiv[isam]   ->Fill(data.dphi*180./pi, wgt);
	hMtv[isam]     ->Fill(data.mt,      wgt);
	hPtv[isam]     ->Fill(data.pt,      wgt);
	hLepDEtav[isam]->Fill(data.leta1-data.leta2,      wgt);
	hMetDPhiv[isam]->Fill(toolbox::deltaPhi(data.phi,data.metphi)*180./pi, wgt);
	hPt1v[isam]    ->Fill(data.lpt1,    wgt);
	hEta1v[isam]   ->Fill(data.leta1,   wgt);
	hPhi1v[isam]   ->Fill(data.lphi1,   wgt);
	hPt2v[isam]    ->Fill(data.lpt2,    wgt);
	hEta2v[isam]   ->Fill(data.leta2,   wgt);
	hPhi2v[isam]   ->Fill(data.lphi2,   wgt);
	hPtMuv[isam]   ->Fill(mupt,         wgt);
	hEtaMuv[isam]  ->Fill(mueta,        wgt);
	hPhiMuv[isam]  ->Fill(muphi,        wgt);
	hPtElev[isam]  ->Fill(elept,        wgt);
	hEtaElev[isam] ->Fill(eleta,        wgt);
	hPhiElev[isam] ->Fill(elephi,       wgt);
	hMtElev[isam]  ->Fill(mtele,        wgt);
	hMtMuv[isam]   ->Fill(mtmu,         wgt);
	if(data.njets>0 && (0.85*data.pvis - data.pmet <= 25)) {
	  hJetPt1v[isam]	->Fill(data.jpt1,     wgt);
	  hJetEta1v[isam]	->Fill(data.jeta1,    wgt);
	  hBoostJetPtv[isam]	->Fill(data.jpt1,     wgt);
	  hBoostJetEtav[isam]	->Fill(data.jeta1,    wgt);
	}
	if(data.njets>1 && (0.85*data.pvis - data.pmet <= 25)) {
	  hJetPt2v[isam]	->Fill(data.jpt2,     wgt);
	  hJetEta2v[isam]	->Fill(data.jeta2,    wgt);
	  hVBFJetPt1v[isam]	->Fill(data.jpt1,     wgt);
	  hVBFJetEta1v[isam]	->Fill(data.jeta1,    wgt);
	  hVBFJetPt2v[isam]	->Fill(data.jpt2,     wgt);
	  hVBFJetEta2v[isam]	->Fill(data.jeta2,    wgt);
	}
	if(data.mjj>0 && (0.85*data.pvis - data.pmet <= 25)) {
	  hJetDPhiv[isam]	->Fill(toolbox::deltaPhi(data.jphi1,data.jphi2)*180./pi,    wgt);
	  hMjjv[isam]		->Fill(data.mjj,      wgt);
	  hDEtav[isam]		->Fill(fabs(data.jeta1 - data.jeta2),    wgt);
	  hEtaProdv[isam]	->Fill(data.jeta1*data.jeta2,            wgt);
	}
	if(data.bjpt>0) {
	  hBJetPtv[isam]	->Fill(data.bjpt,     wgt);
	  hBJetEtav[isam]	->Fill(data.bjeta,    wgt);
	  hBJetPhiv[isam]	->Fill(data.bjphi,    wgt);
	}
	hNPVv[isam]    ->Fill(data.nPV,     wgt);
	if(npuWgt!=0)
	  hNPVrawv[isam] ->Fill(data.nPV,    (isam<1 || snamev[isam].Contains("fake")) ? wgt : wgt/npuWgt);

	hMassv[isam]       ->Fill(data.svfmass,   	 wgt);
	hMassVisv[isam]    ->Fill(data.mass,             wgt);
	hMassLv[isam]      ->Fill(data.svfmass,   	 wgt);      
	hMassHighv[isam]   ->Fill(data.svfmass,   	 wgt);
	hMassVisHighv[isam]->Fill(data.mass,             wgt);
	hMassVpmissv[isam] ->Fill(data.svfmass, data.pmet, wgt);

	nSelv[isam]    += wgt;
	nSelVarv[isam] += wgt*wgt;
      }

      run        = data.runNum;
      lumisec    = data.lumiSec;
      event      = data.evtNum;
      ptele      = elept;    
      etaele     = eleta;   
      phiele     = elephi;   
      relisoele  = 1.0*data.eleiso/elept;   
      ptmu       = mupt;     
      etamu      = mueta;    
      phimu      = muphi;    
      relisomu   = 1.0*data.muiso/mupt;    
      pzeta      = data.pmet - 0.85*data.pvis;    
      njets      = data.njets;    
      npv        = data.nPV;      
      svfmass    = data.svfmass;  
      vismass    = data.mass;  

      Int_t nCentralJets=0;
      for(UInt_t ij=2;ij<data.njets;ij++){
        if(data.jeta1>data.jeta2 && (*jetaArray)[ij]>data.jeta2 && (*jetaArray)[ij]<data.jeta1) nCentralJets++;
        else if(data.jeta2>data.jeta1 && (*jetaArray)[ij]>data.jeta1 && (*jetaArray)[ij]<data.jeta2) nCentralJets++;
      }

      Bool_t vbfcuts = (data.mjj > mjjMin) && (fabs(data.jeta1-data.jeta2) > dEtaMin) && data.jeta1*data.jeta2<0 && nCentralJets==0;
      Bool_t boostcuts = data.jpt1>150;

      Bool_t passpzeta = 0.85*data.pvis - data.pmet <= 25;

      // inclusive
      if (passpzeta) {
	if(isam!=ifake) {
	  hMass_iv[isam]	->Fill(data.svfmass,   wgt);
	  hMassVis_iv[isam]	->Fill(data.mass,      wgt);
	  hMassL_iv[isam]	->Fill(data.svfmass,   wgt);
	  hMassPF_iv[isam]	->Fill(data.svfmass,   wgt);
	  nSel_iv[isam]     += wgt;
	  nSelVar_iv[isam]  += wgt*wgt;
	}
      }

      // 0-jet
      if(passpzeta && (data.njets==0 || (data.njets==1 && !boostcuts))) {
	if(isam!=ifake) {
	  hMass_novbfv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_novbf*wgt                : wgt);      
	  hMassVis_novbfv[isam]	->Fill(data.mass,    (isam==iztt) ? kCat_novbf*wgt                : wgt);
	  hMassL_novbfv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_novbf*wgt                : wgt);      
	  hMassPF_novbfv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_novbf*wgt                : wgt);
	  nSel_novbfv[isam]     +=                (isam==iztt) ? kCat_novbf*wgt                : wgt;
	  nSelVar_novbfv[isam]  +=                (isam==iztt) ? kCat_novbf*wgt*kCat_novbf*wgt : wgt*wgt;
	  if (isam==0) {
	    evtlistfile_novbf << data.runNum << setw(10) << data.lumiSec << setw(15) << data.evtNum << endl;
            //evttree.Fill();
          }
	}
      }

      // boosted
      if(passpzeta && boostcuts && data.nbjets==0 && !vbfcuts) {
	if(isam!=issfake) {
	  hMass_boostv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_boost*wgt                : wgt);
	  hMassVis_boostv[isam]	->Fill(data.mass,    (isam==iztt) ? kCat_boost*wgt                : wgt);
	  hMassL_boostv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_boost*wgt                : wgt);
	  hMassPF_boostv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_boost*wgt                : wgt);
	  nSel_boostv[isam]     +=                (isam==iztt) ? kCat_boost*wgt                : wgt;
	  nSelVar_boostv[isam]  +=                (isam==iztt) ? kCat_boost*wgt*kCat_boost*wgt : wgt*wgt;
	  if (isam==0) {
	    evtlistfile_boost << data.runNum << setw(10) << data.lumiSec << setw(15) << data.evtNum << endl;
	    //evttree.Fill();
	  }
	}
      }

      // vbf
      if(passpzeta && data.njets>=2 && vbfcuts) {
	if(isam!=issfake) {
	  hMass_vbfv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_vbf*wgt              : wgt);      
	  hMassVis_vbfv[isam]	->Fill(data.mass,    (isam==iztt) ? kCat_vbf*wgt              : wgt);
	  hMassL_vbfv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_vbf*wgt              : wgt);      
	  hMassPF_vbfv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_vbf*wgt              : wgt);
	  nSel_vbfv[isam]     +=                (isam==iztt) ? kCat_vbf*wgt              : wgt;
	  nSelVar_vbfv[isam]  +=                (isam==iztt) ? kCat_vbf*wgt*kCat_vbf*wgt : wgt*wgt;
	  for(UInt_t ib=0;ib<npt20jets;ib++)
	    hBdiscr_vbfv[isam]   ->Fill((*btagArray)[ib], (isam==iztt) ? kCat_vbf*wgt           : wgt);
	  if (isam==0) {
	    evtlistfile_vbf << data.runNum << setw(10) << data.lumiSec << setw(15) << data.evtNum << endl;
	    evttree.Fill();
	  }
	}
      }

      // no b-tag
      if(passpzeta && data.njets<=1 && data.nbjets==0) {
	if(isam!=ifake) {
	  hMass_nobv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_nob*wgt              : wgt);      
	  hMassVis_nobv[isam]	->Fill(data.mass,    (isam==iztt) ? kCat_nob*wgt              : wgt);
	  hMassL_nobv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_nob*wgt              : wgt);      
	  hMassPF_nobv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_nob*wgt              : wgt);
	  nSel_nobv[isam]     +=                (isam==iztt) ? kCat_nob*wgt              : wgt;
	  nSelVar_nobv[isam]  +=                (isam==iztt) ? kCat_nob*wgt*kCat_nob*wgt : wgt*wgt;
	  if (isam==0) {
	    evtlistfile_nob << data.runNum << setw(10) << data.lumiSec << setw(15) << data.evtNum << endl;
            //evttree.Fill();
	  }
	}
      }

      // b-tag
      if(passpzeta && data.njets<=1 && data.nbjets>=1) {
	if(isam!=issfake) {
	  hMass_bv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_b*wgt            : wgt);      
	  hMassVis_bv[isam]	->Fill(data.mass,    (isam==iztt) ? kCat_b*wgt            : wgt);
	  hMassL_bv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_b*wgt            : wgt);      
	  hMassPF_bv[isam]	->Fill(data.svfmass, (isam==iztt) ? kCat_b*wgt            : wgt);
	  nSel_bv[isam]     +=                (isam==iztt) ? kCat_b*wgt            : wgt;
	  nSelVar_bv[isam]  +=                (isam==iztt) ? kCat_b*wgt*kCat_b*wgt : wgt*wgt;
	  if (isam==0) {
	    evtlistfile_b << data.runNum << setw(10) << data.lumiSec << setw(15) << data.evtNum << endl;
	    //evttree.Fill();
	  }
	}
      }

    }

    delete infile;
    infile=0, eventTree=0;
    if(isam==0) {
      evtfile.Write();
      evtfile.Close();
    }
    evtlistfile_novbf.close(); evtlistfile_boost.close(); evtlistfile_vbf.close(); evtlistfile_nob.close(); evtlistfile_b.close();

  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  TCanvas *c = MakeCanvas("c","c",700,700);
  
  // string buffers
  char ylabel[100];   // y-axis label

  //----------------------------------------------------------------------------------------
  // write hists to input file for limit inputs and pre-fit plots
  //----------------------------------------------------------------------------------------

  TString flimitstr(CPlot::sOutDir + "/limit-inputs.root");

  TFile *flimit;
  if (ecorr == kCenter) {
    flimit = new TFile(flimitstr,"recreate");
    flimit->mkdir("emu_X");
    flimit->mkdir("emu_novbf");
    flimit->mkdir("emu_boost");
    flimit->mkdir("emu_vbf");
    flimit->mkdir("emu_nob");
    flimit->mkdir("emu_b");
  }
  else flimit = new TFile(flimitstr,"update");

  TString fprefitstr(CPlot::sOutDir + "/prefit-inputs.root");

  TFile *fprefit;
  if (ecorr == kCenter) {
    fprefit = new TFile(fprefitstr,"recreate");
    fprefit->mkdir("emu_X");
    fprefit->mkdir("emu_novbf");
    fprefit->mkdir("emu_boost");
    fprefit->mkdir("emu_vbf");
    fprefit->mkdir("emu_nob");
    fprefit->mkdir("emu_b");
  }
  else fprefit = new TFile(fprefitstr,"update");

  // add zmm into the fakes histogram
  if(hMassPF_iv[issfake] && hMassPF_iv[izmm]) {
    hMassPF_iv[issfake]		->Add(hMassPF_iv[izmm]);
    hMassPF_novbfv[issfake]	->Add(hMassPF_novbfv[izmm]);
    hMassPF_boostv[ifake]	->Add(hMassPF_boostv[izmm]);
    hMassPF_vbfv[ifake]		->Add(hMassPF_vbfv[izmm]);
    hMassPF_nobv[issfake]	->Add(hMassPF_nobv[izmm]);
    hMassPF_bv[ifake]		->Add(hMassPF_bv[izmm]);
    hMassL_iv[issfake]		->Add(hMassL_iv[izmm]);
    hMassL_novbfv[issfake]	->Add(hMassL_novbfv[izmm]);
    hMassL_boostv[ifake]	->Add(hMassL_boostv[izmm]);
    hMassL_vbfv[ifake]		->Add(hMassL_vbfv[izmm]);
    hMassL_nobv[issfake]	->Add(hMassL_nobv[izmm]);
    hMassL_bv[ifake]		->Add(hMassL_bv[izmm]);
  }

  TString histname;
  for(UInt_t isam=0;isam<samplev.size();isam++) {

    if(isam==izmm) continue;

    if(snamev[isam].Contains("ewk",   TString::kIgnoreCase))            histname = "EWK";
    else if(snamev[isam].Contains("fakes", TString::kIgnoreCase))       histname = "Fakes";
    else if(snamev[isam].Contains("ttbar", TString::kIgnoreCase))       histname = "ttbar";
    else if(snamev[isam].Contains("Zmm", TString::kIgnoreCase))         histname = "Zmm";
    else if(snamev[isam].Contains("Ztt",   TString::kIgnoreCase))       histname = "Ztt";
    else if(snamev[isam].Contains("emb",   TString::kIgnoreCase))       histname = "Ztt";
    else if((snamev[isam].Contains("data",  TString::kIgnoreCase)) && !(snamev[isam].Contains("emb",   TString::kIgnoreCase)))          histname = "data_obs";
    else if(snamev[isam].Contains("htt_")) continue;
    else if( (snamev[isam].Contains("gg_")) ||
             (snamev[isam].Contains("bb_")) ||
             (snamev[isam].Contains("gf_")) ||
             (snamev[isam].Contains("vtth_")) ||
             (snamev[isam].Contains("vbf_"))  )                         histname = "Higgs_" + snamev[isam];
    else { cout << "error! name not found" << endl; assert(0); }

    if (isam>0 && isam!=iewk && isam !=ittbar && isam!=ifake && isam!=issfake) {

      if (ecorr==kDown) {
	hMassL_iv[isam]->Scale(rescale(hMassL_iv[isam],"emu_X",histname,flimit));
        hMassL_novbfv[isam]->Scale(rescale(hMassL_novbfv[isam],"emu_novbf",histname,flimit));
        hMassL_boostv[isam]->Scale(rescale(hMassL_boostv[isam],"emu_boost",histname,flimit));
        hMassL_vbfv[isam]->Scale(rescale(hMassL_vbfv[isam],"emu_vbf",histname,flimit));
        hMassL_nobv[isam]->Scale(rescale(hMassL_nobv[isam],"emu_nob",histname,flimit));
        hMassL_bv[isam]->Scale(rescale(hMassL_bv[isam],"emu_b",histname,flimit));
        hMassPF_iv[isam]->Scale(rescale(hMassPF_iv[isam],"emu_X",histname,fprefit));
        hMassPF_novbfv[isam]->Scale(rescale(hMassPF_novbfv[isam],"emu_novbf",histname,fprefit));
        hMassPF_boostv[isam]->Scale(rescale(hMassPF_boostv[isam],"emu_boost",histname,fprefit));
        hMassPF_vbfv[isam]->Scale(rescale(hMassPF_vbfv[isam],"emu_vbf",histname,fprefit));
        hMassPF_nobv[isam]->Scale(rescale(hMassPF_nobv[isam],"emu_nob",histname,fprefit));
        hMassPF_bv[isam]->Scale(rescale(hMassPF_bv[isam],"emu_b",histname,fprefit));
        hMass_iv[isam]->Scale(rescale(hMass_iv[isam],"emu_X",histname,flimit));
        hMass_novbfv[isam]->Scale(rescale(hMass_novbfv[isam],"emu_novbf",histname,flimit));
        hMass_boostv[isam]->Scale(rescale(hMass_boostv[isam],"emu_boost",histname,flimit));
        hMass_vbfv[isam]->Scale(rescale(hMass_vbfv[isam],"emu_vbf",histname,flimit));
        hMass_nobv[isam]->Scale(rescale(hMass_nobv[isam],"emu_nob",histname,flimit));
        hMass_bv[isam]->Scale(rescale(hMass_bv[isam],"emu_b",histname,flimit));
        histname += "_CMS_res_eDown";
      }
      if (ecorr==kUp) {
        hMassL_iv[isam]->Scale(rescale(hMassL_iv[isam],"emu_X",histname,flimit));
        hMassL_novbfv[isam]->Scale(rescale(hMassL_novbfv[isam],"emu_novbf",histname,flimit));
        hMassL_boostv[isam]->Scale(rescale(hMassL_boostv[isam],"emu_boost",histname,flimit));
        hMassL_vbfv[isam]->Scale(rescale(hMassL_vbfv[isam],"emu_vbf",histname,flimit));
        hMassL_nobv[isam]->Scale(rescale(hMassL_nobv[isam],"emu_nob",histname,flimit));
        hMassL_bv[isam]->Scale(rescale(hMassL_bv[isam],"emu_b",histname,flimit));
        hMassPF_iv[isam]->Scale(rescale(hMassPF_iv[isam],"emu_X",histname,fprefit));
        hMassPF_novbfv[isam]->Scale(rescale(hMassPF_novbfv[isam],"emu_novbf",histname,fprefit));
        hMassPF_boostv[isam]->Scale(rescale(hMassPF_boostv[isam],"emu_boost",histname,fprefit));
        hMassPF_vbfv[isam]->Scale(rescale(hMassPF_vbfv[isam],"emu_vbf",histname,fprefit));
        hMassPF_nobv[isam]->Scale(rescale(hMassPF_nobv[isam],"emu_nob",histname,fprefit));
        hMassPF_bv[isam]->Scale(rescale(hMassPF_bv[isam],"emu_b",histname,fprefit));
        hMass_iv[isam]->Scale(rescale(hMass_iv[isam],"emu_X",histname,flimit));
        hMass_novbfv[isam]->Scale(rescale(hMass_novbfv[isam],"emu_novbf",histname,flimit));
        hMass_boostv[isam]->Scale(rescale(hMass_boostv[isam],"emu_boost",histname,flimit));
        hMass_vbfv[isam]->Scale(rescale(hMass_vbfv[isam],"emu_vbf",histname,flimit));
        hMass_nobv[isam]->Scale(rescale(hMass_nobv[isam],"emu_nob",histname,flimit));
        hMass_bv[isam]->Scale(rescale(hMass_bv[isam],"emu_b",histname,flimit));
        histname += "_CMS_res_eUp";
      }

    } else {if (ecorr!=kCenter) continue; }

    if (isam!=ifake) {
      flimit->cd("emu_X");
      hMassL_iv[isam]->SetName(histname);
      hMassL_iv[isam]->Write();
      flimit->cd("emu_novbf");
      hMassL_novbfv[isam]->SetName(histname);
      hMassL_novbfv[isam]->Write();
      flimit->cd("emu_nob");
      hMassL_nobv[isam]->SetName(histname);
      hMassL_nobv[isam]->Write();
    }
    if (isam!=issfake) {
      flimit->cd("emu_boost");
      hMassL_boostv[isam]->SetName(histname);
      hMassL_boostv[isam]->Write();
      flimit->cd("emu_vbf");
      hMassL_vbfv[isam]->SetName(histname);
      hMassL_vbfv[isam]->Write();
      flimit->cd("emu_b");
      hMassL_bv[isam]->SetName(histname);
      hMassL_bv[isam]->Write();
    }

    if (isam!=ifake) {
      fprefit->cd("emu_X");
      hMassPF_iv[isam]->SetName(histname);
      hMassPF_iv[isam]->Write();
      fprefit->cd("emu_novbf");
      hMassPF_novbfv[isam]->SetName(histname);
      hMassPF_novbfv[isam]->Write();
      fprefit->cd("emu_nob");
      hMassPF_nobv[isam]->SetName(histname);
      hMassPF_nobv[isam]->Write();
    }
    if (isam!=issfake) {
      fprefit->cd("emu_boost");
      hMassPF_boostv[isam]->SetName(histname);
      hMassPF_boostv[isam]->Write();
      fprefit->cd("emu_vbf");
      hMassPF_vbfv[isam]->SetName(histname);
      hMassPF_vbfv[isam]->Write();
      fprefit->cd("emu_b");
      hMassPF_bv[isam]->SetName(histname);
      hMassPF_bv[isam]->Write();
    }

  }

  flimit->Close();

  fprefit->Close();
  
  //
  // Begin plots:
  //

  sprintf(ylabel,"Events");
  CPlot plotMet("met","","#slash{E}_{T} [GeV]",ylabel);
  if(hMetv[issfake] && hMetv[izmm]) hMetv[issfake]->Add(hMetv[izmm]);
  if(hasData) { plotMet.AddHist1D(hMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotMet.AddToStack(hMetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMet.SetLegend(0.5,0.65,0.95,0.9);
  plotMet.Draw(c,kTRUE,format);

  // stack up the raw met hists
  TH1F *hTmpMetStack = 0;
  if(hMetRawv[issfake] && hMetRawv[izmm]) hMetRawv[issfake]->Add(hMetRawv[izmm]);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==1) {
      hTmpMetStack = new TH1F(*hMetv[isam]);
      hTmpMetStack->Reset();
    }
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    hTmpMetStack->Add(hMetRawv[isam]);
  }

  sprintf(ylabel,"Events");
  CPlot plotMetRaw("metraw","","#slash{E}_{T} [GeV]",ylabel);
  if(hMetv[ism_vbf] && hMetv[ism_gf]) {
    hMetv[ism_vbf]->Add(hMetv[ism_gf]);
  }
  if(hMetRawv[ism_vbf] && hMetRawv[ism_gf]) {
    hMetRawv[ism_vbf]->Add(hMetRawv[ism_gf]);
  }
  if(hasData) { plotMetRaw.AddHist1D(hMetv[ism_vbf],samplev[ism_vbf]->label+", recoil-corrected","hist",kBlack); }
  if(hasData) { plotMetRaw.AddHist1D(hMetRawv[ism_vbf],samplev[ism_vbf]->label+", uncorrected","hist",kRed+2); }
  plotMetRaw.SetLegend(0.45,0.6,0.95,0.9);
  plotMetRaw.DrawRatio(c,hMetv[ism_vbf],hMetRawv[ism_vbf],kTRUE,format);

  // projection variables
  sprintf(ylabel,"Events");
  CPlot plotProjMet("pzetamiss","","#slash{p}_{#zeta} [GeV]",ylabel);
  if(hProjMetv[issfake] && hProjMetv[izmm]) hProjMetv[issfake]->Add(hProjMetv[izmm]);
  if(hasData) { plotProjMet.AddHist1D(hProjMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotProjMet.AddToStack(hProjMetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotProjMet.SetLegend(0.2,0.65,0.65,0.9);
  plotProjMet.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotProjVis("pzetavis","","p_{#zeta}^{vis} [GeV]",ylabel);
  if(hProjVisv[issfake] && hProjVisv[izmm]) hProjVisv[issfake]->Add(hProjVisv[izmm]);
  if(hasData) { plotProjVis.AddHist1D(hProjVisv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotProjVis.AddToStack(hProjVisv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotProjVis.SetLegend(0.5,0.65,0.95,0.9);
  plotProjVis.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotProjVar("pzetavar","","#slash{p}_{#zeta} - 0.85 #times p_{#zeta}^{vis} [GeV]",ylabel);
  if(hProjVarv[issfake] && hProjVarv[izmm]) hProjVarv[issfake]->Add(hProjVarv[izmm]);
  if(hProjVarv[ism_vbf] && hProjVarv[ism_gf] && hProjVarv[ism_vtth]) {
    hProjVarv[ism_vbf]->Add(hProjVarv[ism_gf]);
    hProjVarv[ism_vbf]->Add(hProjVarv[ism_vtth]);
  }
  if(hasData) { plotProjVar.AddHist1D(hProjVarv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hProjVarv[isam]->Scale(5.);
      plotProjVar.AddToStack(hProjVarv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotProjVar.AddToStack(hProjVarv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotProjVar.SetLegend(0.2,0.6,0.65,0.9);
  plotProjVar.Draw(c,kTRUE,format);

  // stack up the raw projection variable hists
  TH1F *hTmpStack = 0;
  if(hRawProjVarv[issfake] && hRawProjVarv[izmm]) hRawProjVarv[issfake]->Add(hRawProjVarv[izmm]);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==1) {
      hTmpStack = new TH1F(*hProjVarv[isam]);
      hTmpStack->Reset();
    }
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    hTmpStack->Add(hRawProjVarv[isam]);
  }

  sprintf(ylabel,"Events");
  CPlot plotRawProjVar("rawProjVar","","#slash{p}_{#zeta} - 0.85*p_{#zeta}^{vis} [GeV]",ylabel);
  if(hRawProjVarv[ism_vbf] && hRawProjVarv[ism_gf]) {
    hRawProjVarv[ism_vbf]->Add(hRawProjVarv[ism_gf]);
  }
  if(hasData) { plotRawProjVar.AddHist1D(hProjVarv[ism_vbf],samplev[ism_vbf]->label+", recoil-corrected","hist",kBlack); }
  if(hasData) { plotRawProjVar.AddHist1D(hRawProjVarv[ism_vbf],samplev[ism_vbf]->label+", uncorrected","hist",kRed+2); }
  plotRawProjVar.SetLegend(0.45,0.6,0.95,0.9);
  plotRawProjVar.DrawRatio(c,hProjVarv[ism_vbf],hRawProjVarv[ism_vbf],kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotNjets("njets","","Number of Jets",ylabel);
  if(hNjetsv[issfake] && hNjetsv[izmm]) hNjetsv[issfake]->Add(hNjetsv[izmm]);
  if(hNjetsv[ism_vbf] && hNjetsv[ism_gf] && hNjetsv[ism_vtth]) {
    hNjetsv[ism_vbf]->Add(hNjetsv[ism_gf]);
    hNjetsv[ism_vbf]->Add(hNjetsv[ism_vtth]);
  }
  if(hasData) { plotNjets.AddHist1D(hNjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hNjetsv[isam]->Scale(5.);
      plotNjets.AddToStack(hNjetsv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotNjets.AddToStack(hNjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNjets.SetYRange(0,1.5*(plotNjets.GetStack()->GetMaximum()));
  plotNjets.SetLegend(0.5,0.6,0.95,0.9);
  plotNjets.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotNjets_log("njets_log","","Number of Jets",ylabel);
  if(hasData) { plotNjets_log.AddHist1D(hNjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotNjets_log.AddHist1D(hNjetsv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotNjets_log.AddToStack(hNjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNjets_log.SetYRange(0.2,30.0*(plotNjets_log.GetStack()->GetMaximum()));
  plotNjets_log.SetLogy();
  plotNjets_log.SetLegend(0.5,0.6,0.95,0.9);
  plotNjets_log.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotNbjets("nbjets","","Number of b-Tagged Jets",ylabel);
  if(hNbjetsv[issfake] && hNbjetsv[izmm]) hNbjetsv[issfake]->Add(hNbjetsv[izmm]);
  if(hNbjetsv[ism_vbf] && hNbjetsv[ism_gf] && hNbjetsv[ism_vtth]) {
    hNbjetsv[ism_vbf]->Add(hNbjetsv[ism_gf]);
    hNbjetsv[ism_vbf]->Add(hNbjetsv[ism_vtth]);
  }
  if(hasData) { plotNbjets.AddHist1D(hNbjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hNbjetsv[isam]->Scale(5.);
      plotNbjets.AddToStack(hNbjetsv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotNbjets.AddToStack(hNbjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNbjets.SetLegend(0.5,0.6,0.95,0.9);
  plotNbjets.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotNbjets_log("nbjets_log","","Number of b-Tagged Jets",ylabel);
  if(hasData) { plotNbjets_log.AddHist1D(hNbjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotNbjets_log.AddHist1D(hNbjetsv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotNbjets_log.AddToStack(hNbjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNbjets_log.SetYRange(0.2,100.0*(plotNbjets_log.GetStack()->GetMaximum()));
  plotNbjets_log.SetLogy();
  plotNbjets_log.SetLegend(0.5,0.6,0.95,0.9);
  plotNbjets_log.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotBdiscr("btag","","b-Tag Discriminator",ylabel);
  if(hBdiscrv[issfake] && hBdiscrv[izmm]) hBdiscrv[issfake]->Add(hBdiscrv[izmm]);
  if(hBdiscrv[ism_vbf] && hBdiscrv[ism_gf] && hBdiscrv[ism_vtth]) {
    hBdiscrv[ism_vbf]->Add(hBdiscrv[ism_gf]);
    hBdiscrv[ism_vbf]->Add(hBdiscrv[ism_vtth]);
  }
  if(hasData) { plotBdiscr.AddHist1D(hBdiscrv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hBdiscrv[isam]->Scale(5.);
      plotBdiscr.AddToStack(hBdiscrv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotBdiscr.AddToStack(hBdiscrv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBdiscr.SetLegend(0.5,0.6,0.95,0.9);
  plotBdiscr.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotBdiscr_vbf("btag_vbf","","b-tag discr.",ylabel);
  if(hBdiscr_vbfv[ifake] && hBdiscr_vbfv[izmm]) hBdiscr_vbfv[ifake]->Add(hBdiscr_vbfv[izmm]);
  if(hasData) { plotBdiscr_vbf.AddHist1D(hBdiscr_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==issfake)) continue;
    plotBdiscr_vbf.AddToStack(hBdiscr_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBdiscr_vbf.SetLegend(0.5,0.65,0.95,0.9);
  plotBdiscr_vbf.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotDPhi("dphi_emu","","#Delta^{}#phi_{e#mu} [deg]",ylabel);
  if(hDPhiv[issfake] && hDPhiv[izmm]) hDPhiv[issfake]->Add(hDPhiv[izmm]);
  if(hDPhiv[ism_vbf] && hDPhiv[ism_gf] && hDPhiv[ism_vtth]) {
    hDPhiv[ism_vbf]->Add(hDPhiv[ism_gf]);
    hDPhiv[ism_vbf]->Add(hDPhiv[ism_vtth]);
  }
  if(hasData) { plotDPhi.AddHist1D(hDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotDPhi.AddToStack(hDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotDPhi.TransLegend(-0.15,0);
  assert(plotDPhi.GetStack());
  plotDPhi.SetLegend(0.5,0.65,0.95,0.9);
  plotDPhi.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMt("mt","","m_{T}(ll,#slash{E}_{T}) [GeV]",ylabel);
  if(hMtv[issfake] && hMtv[izmm]) hMtv[issfake]->Add(hMtv[izmm]);
  if(hasData) { plotMt.AddHist1D(hMtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotMt.AddToStack(hMtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMt.SetLegend(0.5,0.65,0.95,0.9);
  plotMt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotPt("pt","","p_{T}^{ll} [GeV]",ylabel);
  if(hPtv[issfake] && hPtv[izmm]) hPtv[issfake]->Add(hPtv[izmm]);
  if(hasData) { plotPt.AddHist1D(hPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotPt.AddToStack(hPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPt.SetLegend(0.5,0.65,0.95,0.9);
  plotPt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotLepDEta("lepdeta","","#Delta^{}#eta(ll)",ylabel);
  if(hLepDEtav[issfake] && hLepDEtav[izmm]) hLepDEtav[issfake]->Add(hLepDEtav[izmm]);
  if(hLepDEtav[ism_vbf] && hLepDEtav[ism_gf] && hLepDEtav[ism_vtth]) {
    hLepDEtav[ism_vbf]->Add(hLepDEtav[ism_gf]);
    hLepDEtav[ism_vbf]->Add(hLepDEtav[ism_vtth]);
  }
  if(hasData) { plotLepDEta.AddHist1D(hLepDEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hLepDEtav[isam]->Scale(5.);
      plotLepDEta.AddToStack(hLepDEtav[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotLepDEta.AddToStack(hLepDEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotLepDEta.SetLegend(0.5,0.6,0.95,0.9);
  plotLepDEta.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMetDPhi("metdphi","","#Delta^{}#phi(ll,#slash{E}_{T}) [deg]",ylabel);
  if(hMetDPhiv[issfake] && hMetDPhiv[izmm]) hMetDPhiv[issfake]->Add(hMetDPhiv[izmm]);
  if(hasData) { plotMetDPhi.AddHist1D(hMetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotMetDPhi.AddToStack(hMetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMetDPhi.SetLegend(0.5,0.65,0.95,0.9);
  plotMetDPhi.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPt1("pt1","","leading lepton p_{T} [GeV]",ylabel);
  if(hPt1v[issfake] && hPt1v[izmm]) hPt1v[issfake]->Add(hPt1v[izmm]);
  if(hasData) { plotPt1.AddHist1D(hPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotPt1.AddToStack(hPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPt1.SetLegend(0.5,0.65,0.95,0.9);
  plotPt1.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEta1("eta1","","leading lepton #eta",ylabel);
  if(hEta1v[issfake] && hEta1v[izmm]) hEta1v[issfake]->Add(hEta1v[izmm]);
  if(hasData) { plotEta1.AddHist1D(hEta1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotEta1.AddToStack(hEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEta1.SetYRange(0,2.0*(plotEta1.GetStack()->GetMaximum()));
  plotEta1.SetLegend(0.5,0.65,0.95,0.9);
  plotEta1.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPhi1("phi1","","leading lepton #phi",ylabel);
  if(hPhi1v[issfake] && hPhi1v[izmm]) hPhi1v[issfake]->Add(hPhi1v[izmm]);
  if(hasData) { plotPhi1.AddHist1D(hPhi1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotPhi1.AddToStack(hPhi1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPhi1.SetYRange(0,2.6*(plotPhi1.GetStack()->GetMaximum()));
  plotPhi1.SetLegend(0.5,0.65,0.95,0.9);
  plotPhi1.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotPt2("pt2","","trailing lepton p_{T} [GeV]",ylabel);
  if(hPt2v[issfake] && hPt2v[izmm]) hPt2v[issfake]->Add(hPt2v[izmm]);
  if(hasData) { plotPt2.AddHist1D(hPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotPt2.AddToStack(hPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPt2.SetLegend(0.5,0.65,0.95,0.9);
  plotPt2.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEta2("eta2","","trailing lepton #eta",ylabel);
  if(hEta2v[issfake] && hEta2v[izmm]) hEta2v[issfake]->Add(hEta2v[izmm]);
  if(hasData) { plotEta2.AddHist1D(hEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotEta2.AddToStack(hEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEta2.SetYRange(0,2.0*(plotEta2.GetStack()->GetMaximum()));
  plotEta2.SetLegend(0.5,0.65,0.95,0.9);
  plotEta2.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPhi2("phi2","","trailing lepton #phi",ylabel);
  if(hPhi2v[issfake] && hPhi2v[izmm]) hPhi2v[issfake]->Add(hPhi2v[izmm]);
  if(hasData) { plotPhi2.AddHist1D(hPhi2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotPhi2.AddToStack(hPhi2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPhi2.SetYRange(0,2.6*(plotPhi2.GetStack()->GetMaximum()));
  plotPhi2.SetLegend(0.5,0.65,0.95,0.9);
  plotPhi2.Draw(c,kTRUE,format);

  // mu / electron kinematics
  sprintf(ylabel,"Events");
  CPlot plotPtMu("pt_mu","","#mu p_{T} [GeV]",ylabel);
  if(hPtMuv[issfake] && hPtMuv[izmm]) hPtMuv[issfake]->Add(hPtMuv[izmm]);
  if(hPtMuv[ism_vbf] && hPtMuv[ism_gf] && hPtMuv[ism_vtth]) {
    hPtMuv[ism_vbf]->Add(hPtMuv[ism_gf]);
    hPtMuv[ism_vbf]->Add(hPtMuv[ism_vtth]);
  }
  if(hasData) { plotPtMu.AddHist1D(hPtMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPtMuv[isam]->Scale(5.);
      plotPtMu.AddToStack(hPtMuv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotPtMu.AddToStack(hPtMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtMu.SetLegend(0.5,0.6,0.95,0.9);
  plotPtMu.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEtaMu("eta_mu","","#mu #eta",ylabel);
  if(hEtaMuv[issfake] && hEtaMuv[izmm]) hEtaMuv[issfake]->Add(hEtaMuv[izmm]);
  if(hEtaMuv[ism_vbf] && hEtaMuv[ism_gf] && hEtaMuv[ism_vtth]) {
    hEtaMuv[ism_vbf]->Add(hEtaMuv[ism_gf]);
    hEtaMuv[ism_vbf]->Add(hEtaMuv[ism_vtth]);
  }
  if(hasData) { plotEtaMu.AddHist1D(hEtaMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hEtaMuv[isam]->Scale(5.);
      plotEtaMu.AddToStack(hEtaMuv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotEtaMu.AddToStack(hEtaMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaMu.SetYRange(0,2.0*(plotEtaMu.GetStack()->GetMaximum()));
  plotEtaMu.SetLegend(0.5,0.6,0.95,0.9);
  plotEtaMu.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPhiMu("phimu","","#mu #phi",ylabel);
  if(hPhiMuv[issfake] && hPhiMuv[izmm]) hPhiMuv[issfake]->Add(hPhiMuv[izmm]);
  if(hPhiMuv[ism_vbf] && hPhiMuv[ism_gf] && hPhiMuv[ism_vtth]) {
    hPhiMuv[ism_vbf]->Add(hPhiMuv[ism_gf]);
    hPhiMuv[ism_vbf]->Add(hPhiMuv[ism_vtth]);
  }
  if(hasData) { plotPhiMu.AddHist1D(hPhiMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPhiMuv[isam]->Scale(5.);
      plotPhiMu.AddToStack(hPhiMuv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotPhiMu.AddToStack(hPhiMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPhiMu.SetYRange(0,2.6*(plotPhiMu.GetStack()->GetMaximum()));
  plotPhiMu.SetLegend(0.5,0.6,0.95,0.9);
  plotPhiMu.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotPtEle("pt_e","","e p_{T} [GeV]",ylabel);
  if(hPtElev[issfake] && hPtElev[izmm]) hPtElev[issfake]->Add(hPtElev[izmm]);
  if(hPtElev[ism_vbf] && hPtElev[ism_gf] && hPtElev[ism_vtth]) {
    hPtElev[ism_vbf]->Add(hPtElev[ism_gf]);
    hPtElev[ism_vbf]->Add(hPtElev[ism_vtth]);
  }
  if(hasData) { plotPtEle.AddHist1D(hPtElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPtElev[isam]->Scale(5.);
      plotPtEle.AddToStack(hPtElev[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotPtEle.AddToStack(hPtElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtEle.SetLegend(0.5,0.6,0.95,0.9);
  plotPtEle.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEtaEle("eta_e","","e #eta",ylabel);
  if(hEtaElev[issfake] && hEtaElev[izmm]) hEtaElev[issfake]->Add(hEtaElev[izmm]);
  if(hEtaElev[ism_vbf] && hEtaElev[ism_gf] && hEtaElev[ism_vtth]) {
    hEtaElev[ism_vbf]->Add(hEtaElev[ism_gf]);
    hEtaElev[ism_vbf]->Add(hEtaElev[ism_vtth]);
  }
  if(hasData) { plotEtaEle.AddHist1D(hEtaElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hEtaElev[isam]->Scale(5.);
      plotEtaEle.AddToStack(hEtaElev[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotEtaEle.AddToStack(hEtaElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaEle.SetYRange(0,2.0*(plotEtaEle.GetStack()->GetMaximum()));
  plotEtaEle.SetLegend(0.5,0.6,0.95,0.9);
  plotEtaEle.Draw(c,kTRUE,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPhiEle("phiele","","e #phi",ylabel);
  if(hPhiElev[issfake] && hPhiElev[izmm]) hPhiElev[issfake]->Add(hPhiElev[izmm]);
  if(hPhiElev[ism_vbf] && hPhiElev[ism_gf] && hPhiElev[ism_vtth]) {
    hPhiElev[ism_vbf]->Add(hPhiElev[ism_gf]);
    hPhiElev[ism_vbf]->Add(hPhiElev[ism_vtth]);
  }
  if(hasData) { plotPhiEle.AddHist1D(hPhiElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPhiElev[isam]->Scale(5.);
      plotPhiEle.AddToStack(hPhiElev[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotPhiEle.AddToStack(hPhiElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPhiEle.SetYRange(0,2.6*(plotPhiEle.GetStack()->GetMaximum()));
  plotPhiEle.SetLegend(0.5,0.6,0.95,0.9);
  plotPhiEle.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMtEle("mtele","","m_{T}_{e,MET}",ylabel);
  if(hMtElev[issfake] && hMtElev[izmm]) hMtElev[issfake]->Add(hMtElev[izmm]);
  if(hasData) { plotMtEle.AddHist1D(hMtElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotMtEle.AddToStack(hMtElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMtEle.SetYRange(0,2.6*(plotMtEle.GetStack()->GetMaximum()));
  plotMtEle.SetLegend(0.5,0.65,0.95,0.9);
  plotMtEle.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMtMu("mtmu","","m_{T}_{#mu,MET}",ylabel);
  if(hMtMuv[issfake] && hMtMuv[izmm]) hMtMuv[issfake]->Add(hMtMuv[izmm]);
  if(hasData) { plotMtMu.AddHist1D(hMtMuv[0],samplev[0]->label,"E"); } 
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotMtMu.AddToStack(hMtMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMtMu.SetYRange(0,2.6*(plotMtMu.GetStack()->GetMaximum()));
  plotMtMu.SetLegend(0.5,0.65,0.95,0.9);
  plotMtMu.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotJetPt1("jetpt1","","leading jet pt [GeV]",ylabel);
  if(hJetPt1v[issfake] && hJetPt1v[izmm]) hJetPt1v[issfake]->Add(hJetPt1v[izmm]);
  if(hasData) { plotJetPt1.AddHist1D(hJetPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotJetPt1.AddToStack(hJetPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetPt1.SetLegend(0.5,0.65,0.95,0.9);
  plotJetPt1.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotJetPt2("jetpt2","","second jet pt [GeV]",ylabel);
  if(hJetPt2v[issfake] && hJetPt2v[izmm]) hJetPt2v[issfake]->Add(hJetPt2v[izmm]);
  if(hasData) { plotJetPt2.AddHist1D(hJetPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotJetPt2.AddToStack(hJetPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetPt2.SetLegend(0.5,0.65,0.95,0.9);
  plotJetPt2.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotJetEta1("jeteta1","","leading jet #eta",ylabel);
  if(hJetEta1v[issfake] && hJetEta1v[izmm]) hJetEta1v[issfake]->Add(hJetEta1v[izmm]);
  if(hasData) { plotJetEta1.AddHist1D(hJetEta1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotJetEta1.AddToStack(hJetEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetEta1.SetLegend(0.5,0.65,0.95,0.9);
  plotJetEta1.SetYRange(0,2.0*(plotJetEta1.GetStack()->GetMaximum()));
  plotJetEta1.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotJetEta2("jeteta2","","second jet #eta",ylabel);
  if(hJetEta2v[issfake] && hJetEta2v[izmm]) hJetEta2v[issfake]->Add(hJetEta2v[izmm]);
  if(hasData) { plotJetEta2.AddHist1D(hJetEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotJetEta2.AddToStack(hJetEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetEta2.SetLegend(0.5,0.65,0.95,0.9);
  plotJetEta2.SetYRange(0,2.0*(plotJetEta2.GetStack()->GetMaximum()));
  plotJetEta2.Draw(c,kTRUE,format);

  CPlot plotJetDPhi("jetdphi","","#Delta#phi (jj)",ylabel);
  if(hJetDPhiv[issfake] && hJetDPhiv[izmm]) hJetDPhiv[issfake]->Add(hJetDPhiv[izmm]);
  if(hJetDPhiv[ism_vbf] && hJetDPhiv[ism_gf] && hJetDPhiv[ism_vtth]) {
    hJetDPhiv[ism_vbf]->Add(hJetDPhiv[ism_gf]);
    hJetDPhiv[ism_vbf]->Add(hJetDPhiv[ism_vtth]);
  }
  if(hasData) { plotJetDPhi.AddHist1D(hJetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hJetDPhiv[isam]->Scale(5.);
      plotJetDPhi.AddToStack(hJetDPhiv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotJetDPhi.AddToStack(hJetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetDPhi.SetLegend(0.2,0.6,0.65,0.9);
  plotJetDPhi.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMjj("mjj","","M(jj) [GeV]",ylabel);
  if(hMjjv[issfake] && hMjjv[izmm]) hMjjv[issfake]->Add(hMjjv[izmm]);
  if(hMjjv[ism_vbf] && hMjjv[ism_gf] && hMjjv[ism_vtth]) {
    hMjjv[ism_vbf]->Add(hMjjv[ism_gf]);
    hMjjv[ism_vbf]->Add(hMjjv[ism_vtth]);
  }
  if(hasData) { plotMjj.AddHist1D(hMjjv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMjjv[isam]->Scale(5.);
      plotMjj.AddToStack(hMjjv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMjj.AddToStack(hMjjv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMjj.SetLegend(0.5,0.6,0.95,0.9);
  plotMjj.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotDEta("jetdeta","","#Delta#eta(jj)",ylabel);
  if(hDEtav[issfake] && hDEtav[izmm]) hDEtav[issfake]->Add(hDEtav[izmm]);
  if(hDEtav[ism_vbf] && hDEtav[ism_gf] && hDEtav[ism_vtth]) {
    hDEtav[ism_vbf]->Add(hDEtav[ism_gf]);
    hDEtav[ism_vbf]->Add(hDEtav[ism_vtth]);
  }
  if(hasData) { plotDEta.AddHist1D(hDEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hDEtav[isam]->Scale(5.);
      plotDEta.AddToStack(hDEtav[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotDEta.AddToStack(hDEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotDEta.SetLegend(0.5,0.6,0.95,0.9);
  plotDEta.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotEtaProd("jetEtaProd","","(#eta1*#eta2)(jj)",ylabel);
  if(hEtaProdv[issfake] && hEtaProdv[izmm]) hEtaProdv[issfake]->Add(hEtaProdv[izmm]);
  if(hEtaProdv[ism_vbf] && hEtaProdv[ism_gf] && hEtaProdv[ism_vtth]) {
    hEtaProdv[ism_vbf]->Add(hEtaProdv[ism_gf]);
    hEtaProdv[ism_vbf]->Add(hEtaProdv[ism_vtth]);
  }
  if(hasData) { plotEtaProd.AddHist1D(hEtaProdv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hEtaProdv[isam]->Scale(5.);
      plotEtaProd.AddToStack(hEtaProdv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotEtaProd.AddToStack(hEtaProdv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaProd.SetLegend(0.2,0.6,0.65,0.9);
  plotEtaProd.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotBJetPt("bjetpt","","lead b-jet pt [GeV]",ylabel);
  if(hBJetPtv[issfake] && hBJetPtv[izmm]) hBJetPtv[issfake]->Add(hBJetPtv[izmm]);
  if(hasData) { plotBJetPt.AddHist1D(hBJetPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotBJetPt.AddToStack(hBJetPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBJetPt.SetYRange(0,2.0*(plotBJetPt.GetStack()->GetMaximum()));
  plotBJetPt.SetLegend(0.5,0.65,0.95,0.9);
  plotBJetPt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotBJetEta("bjeteta","","lead b-jet #eta",ylabel);
  if(hBJetEtav[issfake] && hBJetEtav[izmm]) hBJetEtav[issfake]->Add(hBJetEtav[izmm]);
  if(hasData) { plotBJetEta.AddHist1D(hBJetEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotBJetEta.AddToStack(hBJetEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBJetEta.SetYRange(0,2.0*(plotBJetEta.GetStack()->GetMaximum()));
  plotBJetEta.SetLegend(0.5,0.65,0.95,0.9);
  plotBJetEta.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotBJetPhi("bjetphi","","lead b-jet #phi",ylabel);
  if(hBJetPhiv[issfake] && hBJetPhiv[izmm]) hBJetPhiv[issfake]->Add(hBJetPhiv[izmm]);
  if(hasData) { plotBJetPhi.AddHist1D(hBJetPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotBJetPhi.AddToStack(hBJetPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBJetPhi.SetYRange(0,2.1*(plotBJetPhi.GetStack()->GetMaximum()));
  plotBJetPhi.SetLegend(0.5,0.65,0.95,0.9);
  plotBJetPhi.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotBoostJetPt("boostjetpt","","Leading jet p_{T} [GeV]",ylabel);
  if(hBoostJetPtv[issfake] && hBoostJetPtv[izmm]) hBoostJetPtv[issfake]->Add(hBoostJetPtv[izmm]);
  if(hBoostJetPtv[ism_vbf] && hBoostJetPtv[ism_gf] && hBoostJetPtv[ism_vtth]) {
    hBoostJetPtv[ism_vbf]->Add(hBoostJetPtv[ism_gf]);
    hBoostJetPtv[ism_vbf]->Add(hBoostJetPtv[ism_vtth]);
  }
  if(hasData) { plotBoostJetPt.AddHist1D(hBoostJetPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hBoostJetPtv[isam]->Scale(5.);
      plotBoostJetPt.AddToStack(hBoostJetPtv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotBoostJetPt.AddToStack(hBoostJetPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBoostJetPt.SetLegend(0.5,0.6,0.95,0.9);
  plotBoostJetPt.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotBoostJetEta("boostjeteta","","Leading jet #eta",ylabel);
  if(hBoostJetEtav[issfake] && hBoostJetEtav[izmm]) hBoostJetEtav[issfake]->Add(hBoostJetEtav[izmm]);
  if(hBoostJetEtav[ism_vbf] && hBoostJetEtav[ism_gf] && hBoostJetEtav[ism_vtth]) {
    hBoostJetEtav[ism_vbf]->Add(hBoostJetEtav[ism_gf]);
    hBoostJetEtav[ism_vbf]->Add(hBoostJetEtav[ism_vtth]);
  }
  if(hasData) { plotBoostJetEta.AddHist1D(hBoostJetEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hBoostJetEtav[isam]->Scale(5.);
      plotBoostJetEta.AddToStack(hBoostJetEtav[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotBoostJetEta.AddToStack(hBoostJetEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBoostJetEta.SetYRange(0,2.0*(plotBoostJetEta.GetStack()->GetMaximum()));
  plotBoostJetEta.SetLegend(0.5,0.6,0.95,0.9);
  plotBoostJetEta.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotVBFJetPt1("vbfjetpt1","","Leading jet p_{T} [GeV]",ylabel);
  if(hVBFJetPt1v[issfake] && hVBFJetPt1v[izmm]) hVBFJetPt1v[issfake]->Add(hVBFJetPt1v[izmm]);
  if(hVBFJetPt1v[ism_vbf] && hVBFJetPt1v[ism_gf] && hVBFJetPt1v[ism_vtth]) {
    hVBFJetPt1v[ism_vbf]->Add(hVBFJetPt1v[ism_gf]);
    hVBFJetPt1v[ism_vbf]->Add(hVBFJetPt1v[ism_vtth]);
  }
  if(hasData) { plotVBFJetPt1.AddHist1D(hVBFJetPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hVBFJetPt1v[isam]->Scale(5.);
      plotVBFJetPt1.AddToStack(hVBFJetPt1v[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotVBFJetPt1.AddToStack(hVBFJetPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetPt1.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetPt1.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotVBFJetPt2("vbfjetpt2","","Second jet p_{T} [GeV]",ylabel);
  if(hVBFJetPt2v[issfake] && hVBFJetPt2v[izmm]) hVBFJetPt2v[issfake]->Add(hVBFJetPt2v[izmm]);
  if(hVBFJetPt2v[ism_vbf] && hVBFJetPt2v[ism_gf] && hVBFJetPt2v[ism_vtth]) {
    hVBFJetPt2v[ism_vbf]->Add(hVBFJetPt2v[ism_gf]);
    hVBFJetPt2v[ism_vbf]->Add(hVBFJetPt2v[ism_vtth]);
  }
  if(hasData) { plotVBFJetPt2.AddHist1D(hVBFJetPt2v[0],samplev[0]->label,"E"); } 
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hVBFJetPt2v[isam]->Scale(5.);
      plotVBFJetPt2.AddToStack(hVBFJetPt2v[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotVBFJetPt2.AddToStack(hVBFJetPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetPt2.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetPt2.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotVBFJetEta1("vbfjeteta1","","Leading jet #eta",ylabel);
  if(hVBFJetEta1v[issfake] && hVBFJetEta1v[izmm]) hVBFJetEta1v[issfake]->Add(hVBFJetEta1v[izmm]);
  if(hVBFJetEta1v[ism_vbf] && hVBFJetEta1v[ism_gf] && hVBFJetEta1v[ism_vtth]) {
    hVBFJetEta1v[ism_vbf]->Add(hVBFJetEta1v[ism_gf]);
    hVBFJetEta1v[ism_vbf]->Add(hVBFJetEta1v[ism_vtth]);
  }
  if(hasData) { plotVBFJetEta1.AddHist1D(hVBFJetEta1v[0],samplev[0]->label,"E"); } 
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hVBFJetEta1v[isam]->Scale(5.);
      plotVBFJetEta1.AddToStack(hVBFJetEta1v[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotVBFJetEta1.AddToStack(hVBFJetEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetEta1.SetYRange(0,2.0*(plotVBFJetEta1.GetStack()->GetMaximum()));
  plotVBFJetEta1.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetEta1.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotVBFJetEta2("vbfjeteta2","","Second jet #eta",ylabel);
  if(hVBFJetEta2v[issfake] && hVBFJetEta2v[izmm]) hVBFJetEta2v[issfake]->Add(hVBFJetEta2v[izmm]);
  if(hVBFJetEta2v[ism_vbf] && hVBFJetEta2v[ism_gf] && hVBFJetEta2v[ism_vtth]) {
    hVBFJetEta2v[ism_vbf]->Add(hVBFJetEta2v[ism_gf]);
    hVBFJetEta2v[ism_vbf]->Add(hVBFJetEta2v[ism_vtth]);
  }
  if(hasData) { plotVBFJetEta2.AddHist1D(hVBFJetEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hVBFJetEta2v[isam]->Scale(5.);
      plotVBFJetEta2.AddToStack(hVBFJetEta2v[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotVBFJetEta2.AddToStack(hVBFJetEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetEta2.SetYRange(0,2.0*(plotVBFJetEta2.GetStack()->GetMaximum()));
  plotVBFJetEta2.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetEta2.Draw(c,kTRUE,format);

  CPlot plotNPV("nvertices_reweighted","","N_{PV}","Events");
  if(hNPVv[issfake] && hNPVv[izmm]) hNPVv[issfake]->Add(hNPVv[izmm]);
  if(hasData) { plotNPV.AddHist1D(hNPVv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotNPV.AddToStack(hNPVv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNPV.SetYRange(0,1.3*(hNPVv[0]->GetMaximum()));
  plotNPV.SetLegend(0.5,0.65,0.95,0.9);
  plotNPV.Draw(c,kTRUE,format);

  CPlot plotNPVraw("nvertices_raw","","N_{PV}","Events");
  if(hNPVrawv[issfake] && hNPVrawv[izmm]) hNPVrawv[issfake]->Add(hNPVrawv[izmm]);
  if(hasData) { plotNPVraw.AddHist1D(hNPVrawv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    if((isam==izmm) || (isam==ifake)) continue;
    plotNPVraw.AddToStack(hNPVrawv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNPVraw.SetYRange(0,1.3*(hNPVrawv[0]->GetMaximum()));
  plotNPVraw.SetLegend(0.5,0.65,0.95,0.9);
  plotNPVraw.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMass("svfitmass_lept","","m_{#tau#tau} [GeV]",ylabel);
  if(hMassv[issfake] && hMassv[izmm]) hMassv[issfake]->Add(hMassv[izmm]);
  if(hMassv[ism_vbf] && hMassv[ism_gf] && hMassv[ism_vtth]) {
    hMassv[ism_vbf]->Add(hMassv[ism_gf]);
    hMassv[ism_vbf]->Add(hMassv[ism_vtth]);
  }
  if(hasData) { plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassv[isam]->Scale(5.);
      plotMass.AddToStack(hMassv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass.SetLegend(0.5,0.6,0.95,0.9);
  plotMass.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis("vismass_lept","","m_{vis} [GeV]",ylabel);
  if(hMassVisv[issfake] && hMassVisv[izmm]) hMassVisv[issfake]->Add(hMassVisv[izmm]);
  if(hMassVisv[ism_vbf] && hMassVisv[ism_gf] && hMassVisv[ism_vtth]) {
    hMassVisv[ism_vbf]->Add(hMassVisv[ism_gf]);
    hMassVisv[ism_vbf]->Add(hMassVisv[ism_vtth]);
  }
  if(hasData) { plotMassVis.AddHist1D(hMassVisv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVisv[isam]->Scale(5.);
      plotMassVis.AddToStack(hMassVisv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis.AddToStack(hMassVisv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassHigh("svfitmass_lept-high","","m_{#tau#tau} [GeV]",ylabel);
  if(hMassHighv[issfake] && hMassHighv[izmm]) hMassHighv[issfake]->Add(hMassHighv[izmm]);
  if(hasData) { plotMassHigh.AddHist1D(hMassHighv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    plotMassHigh.AddToStack(hMassHighv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassHigh.SetLogy();
  plotMassHigh.SetLegend(0.5,0.65,0.95,0.9);
  plotMassHigh.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVisHigh("vismass_lept-high","","m_{vis} [GeV]",ylabel);
  if(hMassVisHighv[issfake] && hMassVisHighv[izmm]) hMassVisHighv[issfake]->Add(hMassVisHighv[izmm]);
  if(hasData) { plotMassVisHigh.AddHist1D(hMassVisHighv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    plotMassVisHigh.AddToStack(hMassVisHighv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  } 
  plotMassVisHigh.SetLogy();
  plotMassVisHigh.SetLegend(0.5,0.65,0.95,0.9);
  plotMassVisHigh.Draw(c,kTRUE,format);

  CPlot plotMassVpmiss("massvpmiss","","m_{#tau#tau} [GeV]}","#slash{p}_{#zeta} [GeV]}");
  assert(hMassVpmissv[0]);
  plotMassVpmiss.AddHist2D((TH2D*)hMassVpmissv[0],"surf",kWhite,kBlue);

  //----------------------------------------------------------------------------------------
  // inclusive
  //----------------------------------------------------------------------------------------  
  
  sprintf(ylabel,"Events");
  CPlot plotMass_i("svfitmass_incl","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_iv[issfake] && hMass_iv[izmm]) hMass_iv[issfake]->Add(hMass_iv[izmm]);
  if(hMass_iv[ism_vbf] && hMass_iv[ism_gf] && hMass_iv[ism_vtth]) {
    hMass_iv[ism_vbf]->Add(hMass_iv[ism_gf]);
    hMass_iv[ism_vbf]->Add(hMass_iv[ism_vtth]);
  }
  if(hasData) { plotMass_i.AddHist1D(hMass_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_iv[isam]->Scale(5.);
      plotMass_i.AddToStack(hMass_iv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_i.AddToStack(hMass_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_i.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_i.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_i("vismass_incl","","m_{vis} [GeV]",ylabel);
  if(hMassVis_iv[issfake] && hMassVis_iv[izmm]) hMassVis_iv[issfake]->Add(hMassVis_iv[izmm]);
  if(hMassVis_iv[ism_vbf] && hMassVis_iv[ism_gf] && hMassVis_iv[ism_vtth]) {
    hMassVis_iv[ism_vbf]->Add(hMassVis_iv[ism_gf]);
    hMassVis_iv[ism_vbf]->Add(hMassVis_iv[ism_vtth]);
  }
  if(hasData) { plotMassVis_i.AddHist1D(hMassVis_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_iv[isam]->Scale(5.);
      plotMassVis_i.AddToStack(hMassVis_iv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_i.AddToStack(hMassVis_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_i.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_i.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMass_i_log("svfitmass_incl_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_i_log.AddHist1D(hMass_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_i_log.AddHist1D(hMass_iv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_i_log.AddToStack(hMass_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_i_log.SetLogy();
  plotMass_i_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_i_log.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_i_log("vismass_incl_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_i_log.AddHist1D(hMassVis_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_i_log.AddHist1D(hMassVis_iv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_i_log.AddToStack(hMassVis_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_i_log.SetLogy();
  plotMassVis_i_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_i_log.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // no vbf
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"Events");
  CPlot plotMass_novbf("svfitmass_class_novbf","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_novbfv[issfake] && hMass_novbfv[izmm]) hMass_novbfv[issfake]->Add(hMass_novbfv[izmm]);
  if(hMass_novbfv[ism_vbf] && hMass_novbfv[ism_gf] && hMass_novbfv[ism_vtth]) {
    hMass_novbfv[ism_vbf]->Add(hMass_novbfv[ism_gf]);
    hMass_novbfv[ism_vbf]->Add(hMass_novbfv[ism_vtth]);
  }
  if(hasData) { plotMass_novbf.AddHist1D(hMass_novbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_novbfv[isam]->Scale(5.);
      plotMass_novbf.AddToStack(hMass_novbfv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_novbf.AddToStack(hMass_novbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_novbf.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_novbf.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_novbf("vismass_class_novbf","","m_{vis} [GeV]",ylabel);
  if(hMassVis_novbfv[issfake] && hMassVis_novbfv[izmm]) hMassVis_novbfv[issfake]->Add(hMassVis_novbfv[izmm]);
  if(hMassVis_novbfv[ism_vbf] && hMassVis_novbfv[ism_gf] && hMassVis_novbfv[ism_vtth]) {
    hMassVis_novbfv[ism_vbf]->Add(hMassVis_novbfv[ism_gf]);
    hMassVis_novbfv[ism_vbf]->Add(hMassVis_novbfv[ism_vtth]);
  }
  if(hasData) { plotMassVis_novbf.AddHist1D(hMassVis_novbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_novbfv[isam]->Scale(5.);
      plotMassVis_novbf.AddToStack(hMassVis_novbfv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_novbf.AddToStack(hMassVis_novbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_novbf.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_novbf.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMass_novbf_log("svfitmass_class_novbf_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_novbf_log.AddHist1D(hMass_novbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_novbf_log.AddHist1D(hMass_novbfv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_novbf_log.AddToStack(hMass_novbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_novbf_log.SetLogy();
  plotMass_novbf_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_novbf_log.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_novbf_log("vismass_class_novbf_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_novbf_log.AddHist1D(hMassVis_novbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_novbf_log.AddHist1D(hMassVis_novbfv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_novbf_log.AddToStack(hMassVis_novbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_novbf_log.SetLogy();
  plotMassVis_novbf_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_novbf_log.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // boosted
  //----------------------------------------------------------------------------------------

  sprintf(ylabel,"Events");
  CPlot plotMass_boost("svfitmass_class_boost","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_boostv[ifake] && hMass_boostv[izmm]) hMass_boostv[ifake]->Add(hMass_boostv[izmm]);
  if(hMass_boostv[ism_vbf] && hMass_boostv[ism_gf] && hMass_boostv[ism_vtth]) {
    hMass_boostv[ism_vbf]->Add(hMass_boostv[ism_gf]);
    hMass_boostv[ism_vbf]->Add(hMass_boostv[ism_vtth]);
  }
  if(hasData) { plotMass_boost.AddHist1D(hMass_boostv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_boostv[isam]->Scale(5.);
      plotMass_boost.AddToStack(hMass_boostv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_boost.AddToStack(hMass_boostv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_boost.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_boost.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_boost("vismass_class_boost","","m_{vis} [GeV]",ylabel);
  if(hMassVis_boostv[ifake] && hMassVis_boostv[izmm]) hMassVis_boostv[ifake]->Add(hMassVis_boostv[izmm]);
  if(hMassVis_boostv[ism_vbf] && hMassVis_boostv[ism_gf] && hMassVis_boostv[ism_vtth]) {
    hMassVis_boostv[ism_vbf]->Add(hMassVis_boostv[ism_gf]);
    hMassVis_boostv[ism_vbf]->Add(hMassVis_boostv[ism_vtth]);
  }
  if(hasData) { plotMassVis_boost.AddHist1D(hMassVis_boostv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_boostv[isam]->Scale(5.);
      plotMassVis_boost.AddToStack(hMassVis_boostv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_boost.AddToStack(hMassVis_boostv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_boost.SetYRange(0,1.3*(plotMassVis_boost.GetStack()->GetMaximum()));
  plotMassVis_boost.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_boost.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMass_boost_log("svfitmass_class_boost_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_boost_log.AddHist1D(hMass_boostv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_boost_log.AddHist1D(hMass_boostv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_boost_log.AddToStack(hMass_boostv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_boost_log.SetLogy();
  plotMass_boost_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_boost_log.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_boost_log("vismass_class_boost_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_boost_log.AddHist1D(hMassVis_boostv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_boost_log.AddHist1D(hMassVis_boostv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_boost_log.AddToStack(hMassVis_boostv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_boost_log.SetLogy();
  plotMassVis_boost_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_boost_log.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // vbf
  //----------------------------------------------------------------------------------------

  sprintf(ylabel,"Events");
  CPlot plotMass_vbf("svfitmass_class_vbf","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_vbfv[ifake] && hMass_vbfv[izmm]) hMass_vbfv[ifake]->Add(hMass_vbfv[izmm]);
  if(hMass_vbfv[ism_vbf] && hMass_vbfv[ism_gf] && hMass_vbfv[ism_vtth]) {
    hMass_vbfv[ism_vbf]->Add(hMass_vbfv[ism_gf]);
    hMass_vbfv[ism_vbf]->Add(hMass_vbfv[ism_vtth]);
  }
  if(hasData) { plotMass_vbf.AddHist1D(hMass_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_vbfv[isam]->Scale(5.);
      plotMass_vbf.AddToStack(hMass_vbfv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_vbf.AddToStack(hMass_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_vbf.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_vbf.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_vbf("vismass_class_vbf","","m_{vis} [GeV]",ylabel);
  if(hMassVis_vbfv[ifake] && hMassVis_vbfv[izmm]) hMassVis_vbfv[ifake]->Add(hMassVis_vbfv[izmm]);
  if(hMassVis_vbfv[ism_vbf] && hMassVis_vbfv[ism_gf] && hMassVis_vbfv[ism_vtth]) {
    hMassVis_vbfv[ism_vbf]->Add(hMassVis_vbfv[ism_gf]);
    hMassVis_vbfv[ism_vbf]->Add(hMassVis_vbfv[ism_vtth]);
  }
  if(hasData) { plotMassVis_vbf.AddHist1D(hMassVis_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_vbfv[isam]->Scale(5.);
      plotMassVis_vbf.AddToStack(hMassVis_vbfv[isam],"5#times "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_vbf.AddToStack(hMassVis_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_vbf.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_vbf.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMass_vbf_log("svfitmass_class_vbf_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_vbf_log.AddHist1D(hMass_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) plotMass_vbf_log.AddHist1D(hMass_vbfv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_vbf_log.AddToStack(hMass_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_vbf_log.SetLogy();
  plotMass_vbf_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_vbf_log.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_vbf_log("vismass_class_vbf_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_vbf_log.AddHist1D(hMassVis_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) plotMassVis_vbf_log.AddHist1D(hMassVis_vbfv[isam],"5#times "+samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_vbf_log.AddToStack(hMassVis_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_vbf_log.SetLogy();
  plotMassVis_vbf_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_vbf_log.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // no b-tag
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"Events");
  CPlot plotMass_nob("svfitmass_class_nob","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_nobv[issfake] && hMass_nobv[izmm]) hMass_nobv[issfake]->Add(hMass_nobv[izmm]);
  if(hMass_nobv[imssm_gg] && hMass_nobv[imssm_bb]) hMass_nobv[imssm_gg]->Add(hMass_nobv[imssm_bb]);
  if(hasData) { plotMass_nob.AddHist1D(hMass_nobv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_gg")) plotMass_nob.AddToStack(hMass_nobv[isam],samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_nob.AddToStack(hMass_nobv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_nob.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_nob.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_nob("vismass_class_nob","","m_{vis} [GeV]",ylabel);
  if(hMassVis_nobv[issfake] && hMassVis_nobv[izmm]) hMassVis_nobv[issfake]->Add(hMassVis_nobv[izmm]);
  if(hMassVis_nobv[imssm_gg] && hMassVis_nobv[imssm_bb]) hMassVis_nobv[imssm_gg]->Add(hMassVis_nobv[imssm_bb]);
  if(hasData) { plotMassVis_nob.AddHist1D(hMassVis_nobv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_gg")) plotMassVis_nob.AddToStack(hMassVis_nobv[isam],samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_nob.AddToStack(hMassVis_nobv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_nob.SetLegend(0.7,0.65,0.9,0.9);
  plotMassVis_nob.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_nob.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMass_nob_log("svfitmass_class_nob_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_nob_log.AddHist1D(hMass_nobv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_gg")) plotMass_nob_log.AddHist1D(hMass_nobv[isam],samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_nob_log.AddToStack(hMass_nobv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_nob_log.SetLogy();
  plotMass_nob_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_nob_log.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_nob_log("vismass_class_nob_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_nob_log.AddHist1D(hMassVis_nobv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_gg")) plotMassVis_nob_log.AddHist1D(hMassVis_nobv[isam],samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_nob_log.AddToStack(hMassVis_nobv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_nob_log.SetLogy();
  plotMassVis_nob_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_nob_log.Draw(c,kTRUE,format);

  //----------------------------------------------------------------------------------------
  // b-tag
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"Events");
  CPlot plotMass_b("svfitmass_class_b","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_bv[ifake] && hMass_bv[izmm]) hMass_bv[ifake]->Add(hMass_bv[izmm]);
  if(hMass_bv[imssm_gg] && hMass_bv[imssm_bb]) hMass_bv[imssm_gg]->Add(hMass_bv[imssm_bb]);
  if(hasData) { plotMass_b.AddHist1D(hMass_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_gg")) plotMass_b.AddToStack(hMass_bv[isam],samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_b.AddToStack(hMass_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_b.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_b("vismass_class_b","","m_{vis} [GeV]",ylabel);
  if(hMassVis_bv[ifake] && hMassVis_bv[izmm]) hMassVis_bv[ifake]->Add(hMassVis_bv[izmm]);
  if(hMassVis_bv[imssm_gg] && hMassVis_bv[imssm_bb]) hMassVis_bv[imssm_gg]->Add(hMassVis_bv[imssm_bb]);
  if(hasData) { plotMassVis_b.AddHist1D(hMassVis_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_gg")) plotMassVis_b.AddToStack(hMassVis_bv[isam],samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_b.AddToStack(hMassVis_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_b.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_b.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMass_b_log("svfitmass_class_b_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_b_log.AddHist1D(hMass_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_gg")) plotMass_b_log.AddHist1D(hMass_bv[isam],samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMass_b_log.AddToStack(hMass_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b_log.SetLogy();
  plotMass_b_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_b_log.Draw(c,kTRUE,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVis_b_log("vismass_class_b_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_b_log.AddHist1D(hMassVis_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_gg")) plotMassVis_b_log.AddHist1D(hMassVis_bv[isam],samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_")) continue;
    else plotMassVis_b_log.AddToStack(hMassVis_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_b_log.SetLogy();
  plotMassVis_b_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_b_log.Draw(c,kTRUE,format);

    
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
  Double_t nTotal_boost=0, nTotalVar_boost=0;
  Double_t nTotal_vbf=0,   nTotalVar_vbf=0;
  Double_t nTotal_nob=0,   nTotalVar_nob=0;
  Double_t nTotal_b=0,     nTotalVar_b=0;

  ofstream yieldfile;
  yieldfile.open(CPlot::sOutDir+"/yields.txt");
  
  yieldfile << setw(33) << "lepton sele." << setw(20) << "inclusive" << setw(20) << "no vbf" << setw(20) << "boosted" << setw(20) << "vbf" << setw(20) << "no b-tag" << setw(20) << "b-tag" << endl;

  if(hMass_iv[issfake] && hMass_iv[izmm]) { // add zmm into the fakes
    nSelv[issfake]       += nSelv[izmm];       nSelVarv[issfake]       += nSelVarv[izmm];
    nSel_iv[issfake]     += nSel_iv[izmm];     nSelVar_iv[issfake]     += nSelVar_iv[izmm];
    nSel_novbfv[issfake] += nSel_novbfv[izmm]; nSelVar_novbfv[issfake] += nSelVar_novbfv[izmm];
    nSel_boostv[ifake] += nSel_boostv[izmm]; nSelVar_boostv[ifake] += nSelVar_boostv[izmm];
    nSel_vbfv[ifake]   += nSel_vbfv[izmm];   nSelVar_vbfv[ifake]   += nSelVar_vbfv[izmm];
    nSel_nobv[issfake]   += nSel_nobv[izmm];   nSelVar_nobv[issfake]   += nSelVar_nobv[izmm];
    nSel_bv[ifake]     += nSel_bv[izmm];     nSelVar_bv[ifake]     += nSelVar_bv[izmm];
  }

  yieldfile << setw(15) << "fakes";
  yieldfile << setw(10) << setprecision(3) << fixed << nSelv[issfake]        << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[issfake]);
  yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[issfake]      << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[issfake]);
  yieldfile << setw(10) << setprecision(3) << fixed << nSel_novbfv[issfake]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_novbfv[issfake]);
  yieldfile << setw(10) << setprecision(3) << fixed << nSel_boostv[ifake]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boostv[ifake]);
  yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[ifake]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[ifake]);
  yieldfile << setw(10) << setprecision(3) << fixed << nSel_nobv[issfake]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_nobv[issfake]);
  yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[ifake]      << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_bv[ifake]);
  yieldfile << endl;

  nTotal            += nSelv[issfake];
  nTotalVar         += nSelVarv[issfake];
  nTotal_i          += nSel_iv[issfake];
  nTotalVar_i       += nSelVar_iv[issfake];
  nTotal_novbf      += nSel_novbfv[issfake];
  nTotalVar_novbf   += nSelVar_novbfv[issfake];
  nTotal_boost      += nSel_boostv[ifake];
  nTotalVar_boost   += nSelVar_boostv[ifake];
  nTotal_vbf        += nSel_vbfv[ifake];
  nTotalVar_vbf     += nSelVar_vbfv[ifake];
  nTotal_nob        += nSel_nobv[issfake];
  nTotalVar_nob     += nSelVar_nobv[issfake];
  nTotal_b          += nSel_bv[ifake];
  nTotalVar_b       += nSelVar_bv[ifake];


  if(samplev.size()>1) { 
    for(UInt_t isam=1; isam<samplev.size(); isam++) {
      if(isam==izmm) continue;
      else if(snamev[isam].Contains("htt_"))  continue;
      else if((isam==ifake) || (isam==issfake)) continue;
      else {
        yieldfile << setw(15) << snamev[isam];
        yieldfile << setw(10) << setprecision(3) << fixed << nSelv[isam]        << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[isam]);
        yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[isam]      << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[isam]);
        yieldfile << setw(10) << setprecision(3) << fixed << nSel_novbfv[isam]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_novbfv[isam]);
        yieldfile << setw(10) << setprecision(3) << fixed << nSel_boostv[isam]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boostv[isam]);
        yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[isam]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[isam]);
        yieldfile << setw(10) << setprecision(3) << fixed << nSel_nobv[isam]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_nobv[isam]);
        yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[isam]      << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_bv[isam]);
        yieldfile << endl;
      }

      if(snamev[isam].Contains("_mssm_"))  continue; // don't add higgs samples to total
      if(snamev[isam].Contains("_sm_"))    continue;

      nTotal    	+= nSelv[isam];
      nTotalVar 	+= nSelVarv[isam];
      nTotal_i      	+= nSel_iv[isam];
      nTotalVar_i   	+= nSelVar_iv[isam];
      nTotal_novbf   	+= nSel_novbfv[isam];
      nTotalVar_novbf 	+= nSelVar_novbfv[isam];
      nTotal_boost      += nSel_boostv[isam];
      nTotalVar_boost   += nSelVar_boostv[isam];
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
    yieldfile << setw(10) << setprecision(3) << fixed << nTotal_boost << " +/- " << sqrt(nTotalVar_boost);
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
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boostv[0]  << " +/- " << sqrt(nSelVar_boostv[0]);
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
  yieldfile << setw(20) << "lepton sele." <<setw(17) << "inclusive" << setw(19) << "no vbf" << setw(19) << "boosted" << setw(19) << "vbf" << setw(22) << "no b-tag" << setw(17) << "b-tag" << endl;

  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if(!(snamev[isam].Contains("_mssm_") || snamev[isam].Contains("_sm_"))) continue;
    if(snamev[isam].Contains("htt_")) continue;
    yieldfile << setw(25) << snamev[isam] << "  " <<
      setw(15) << nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSelv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_iv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_novbfv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_boostv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_vbfv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_nobv[isam]/nUnskEventsv[isam] <<
      setw(20) << setprecision(6) <<nSel_bv[isam]/nUnskEventsv[isam] << endl;
  }
  yieldfile.close();


  // need to clean up how this is done
  ofstream texfile;
  texfile.open(CPlot::sOutDir+"/yields.tex");


  texfile << " \\begin{table}[!ht]" << endl;
  texfile << " \\begin{center}" << endl;
  texfile << " \\begin{tabular}{|l|r@{$ \\,\\,\\pm\\,\\, $}l|r@{$\\,\\,\\pm\\,\\,$}l|r@{$\\,\\,\\pm\\,\\,$}l|r@{$\\,\\,\\pm\\,\\,$}l|r@{$\\,\\,\\pm\\,\\,$}l|}"  << endl;
  texfile << " \\cline{2-11}" << endl;

  texfile << " \\multicolumn{1}{c}{ } &  \\multicolumn{6}{|c|}{Standard Model}  & \\multicolumn{4}{|c|}{MSSM} \\\\" << endl;
  texfile << " \\hline" << endl;
  texfile << " Process & \\multicolumn{2}{|c|}{\\emph{0/1-Jet}} & \\multicolumn{2}{|c|}{\\emph{Boost}} & \\multicolumn{2}{|c|}{\\emph{VBF}} &  \\multicolumn{2}{|c|}{\\emph{Non B-Tag}} & \\multicolumn{2}{|c|}{\\emph{B-Tag}} \\\\" << endl;
  texfile << " \\hline" << endl;

  texfile << endl;
  texfile << " $Z\\rightarrow \\tau\\tau$       & " << setw(5) << setprecision(0) << fixed << nSel_novbfv[iemb] << " & " << setw(5) << setprecision(0) << fixed << 0.067*nSel_novbfv[iemb] << " & " << setw(5) << setprecision(0) << fixed << nSel_boostv[iemb] << " & " << setw(5) << setprecision(0) << fixed << 0.11*nSel_boostv[iemb] << " & " << setw(5) << setprecision(0) << fixed << nSel_vbfv[iemb] << " & " << setw(5) << setprecision(0) << fixed << 0.23*nSel_vbfv[iemb] << " & " << setw(5) << setprecision(0) << fixed << nSel_nobv[iemb] << " & "  << setw(5) << setprecision(0) << fixed << 0.068*nSel_nobv[iemb] << " & " << setw(5) << setprecision(0) << fixed << nSel_bv[iemb] << " & " << setw(5) << setprecision(0) << fixed << 0.099*nSel_bv[iemb] << "  \\\\";
  texfile << endl;
  texfile << " Fakes                         & " << setw(5) << setprecision(0) << fixed << nSel_novbfv[issfake] << " & " << setw(5) << setprecision(0) << fixed << 0.3*nSel_novbfv[issfake] << " & " << setw(5) << setprecision(0) << fixed << nSel_boostv[ifake] << " & " << setw(5) << setprecision(0) << fixed << 0.31*nSel_boostv[ifake] << " & " << setw(5) << setprecision(0) << fixed << nSel_vbfv[ifake] << " & " << setw(5) << setprecision(0) << fixed << 0.36*nSel_vbfv[ifake] << " & " << setw(5) << setprecision(0) << fixed << nSel_nobv[issfake] << " & "  << setw(5) << setprecision(0) << fixed << 0.31*nSel_nobv[issfake] << " & " << setw(5) << setprecision(0) << fixed << nSel_bv[ifake] << " & " << setw(5) << setprecision(0) << fixed << 0.31*nSel_bv[ifake] << "  \\\\";
  texfile << endl;
  texfile << " $t\\bar{t}$                    & " << setw(5) << setprecision(0) << fixed << nSel_novbfv[ittbar] << " & " << setw(5) << setprecision(0) << fixed << 0.095*nSel_novbfv[ittbar] << " & " << setw(5) << setprecision(0) << fixed << nSel_boostv[ittbar] << " & " << setw(5) << setprecision(0) << fixed << 0.12*nSel_boostv[ittbar] << " & " << setw(5) << setprecision(0) << fixed << nSel_vbfv[ittbar] << " & " << setw(5) << setprecision(0) << fixed << 0.22*nSel_vbfv[ittbar] << " & " << setw(5) << setprecision(0) << fixed << nSel_nobv[ittbar] << " & "  << setw(5) << setprecision(0) << fixed << 0.096*nSel_nobv[ittbar] << " & " << setw(5) << setprecision(0) << fixed << nSel_bv[ittbar] << " & " << setw(5) << setprecision(0) << fixed << 0.12*nSel_bv[ittbar] << "  \\\\";
  texfile << endl;
  texfile << " Di-Boson                      & " << setw(5) << setprecision(0) << fixed << nSel_novbfv[iewk] << " & " << setw(5) << setprecision(0) << fixed << 0.16*nSel_novbfv[iewk] << " & " << setw(5) << setprecision(0) << fixed << nSel_boostv[iewk] << " & " << setw(5) << setprecision(0) << fixed << 0.18*nSel_boostv[iewk] << " & " << setw(5) << setprecision(0) << fixed << nSel_vbfv[iewk] << " & " << setw(5) << setprecision(1) << fixed << 0.26*nSel_vbfv[iewk] << " & " << setw(5) << setprecision(0) << fixed << nSel_nobv[iewk] << " & "  << setw(5) << setprecision(0) << fixed << 0.16*nSel_nobv[iewk] << " & " << setw(5) << setprecision(0) << fixed << nSel_bv[iewk] << " & " << setw(5) << setprecision(0) << fixed << 0.18*nSel_bv[iewk] << "  \\\\";
  texfile << endl;


  texfile << " \\hline" << endl;
  texfile << " \\hline" << endl;

  Double_t nSystVar_novbf = sqrt(pow(0.16*nSel_novbfv[iewk],2)+pow(0.095*nSel_novbfv[ittbar],2)+pow(0.3*nSel_novbfv[issfake],2)+pow(0.067*nSel_novbfv[iemb],2));
  Double_t nSystVar_boost = sqrt(pow(0.18*nSel_boostv[iewk],2)+pow(0.12*nSel_boostv[ittbar],2)+pow(0.31*nSel_boostv[issfake],2)+pow(0.11*nSel_boostv[iemb],2));
  Double_t nSystVar_vbf = sqrt(pow(0.26*nSel_vbfv[iewk],2)+pow(0.22*nSel_vbfv[ittbar],2)+pow(0.36*nSel_vbfv[issfake],2)+pow(0.23*nSel_vbfv[iemb],2));
  Double_t nSystVar_nob = sqrt(pow(0.16*nSel_nobv[iewk],2)+pow(0.096*nSel_nobv[ittbar],2)+pow(0.31*nSel_nobv[issfake],2)+pow(0.068*nSel_nobv[iemb],2));
  Double_t nSystVar_b = sqrt(pow(0.18*nSel_bv[iewk],2)+pow(0.12*nSel_bv[ittbar],2)+pow(0.31*nSel_bv[issfake],2)+pow(0.099*nSel_bv[iemb],2));


  texfile << " Total Background              & " << setw(5) << setprecision(0) << fixed << nTotal_novbf << " & " << setw(5) << setprecision(0) << fixed << nSystVar_novbf << " & " << setw(5) << setprecision(0) << fixed << nTotal_boost << " & " << setw(5) << setprecision(0) << fixed << nSystVar_boost << " & " << setw(5) << setprecision(0) << fixed << nTotal_vbf << " & " << setw(5) << setprecision(0) << fixed << nSystVar_vbf << " & " << setw(5) << setprecision(0) << fixed << nTotal_nob << " & "  << setw(5) << setprecision(0) << fixed << nSystVar_nob << " & " << setw(5) << setprecision(0) << fixed << nTotal_b << " & " << setw(5) << setprecision(0) << fixed << nSystVar_b << "  \\\\";
  texfile << endl;
  texfile << " \\hline" << endl;

  Double_t nTotalSM120_novbf    = nSel_novbfv[ism_gf] + nSel_novbfv[ism_vbf];     
  Double_t nTotalSM120Var_novbf = 0.16*nSel_novbfv[ism_gf] + 0.11*nSel_novbfv[ism_vbf];
  Double_t nTotalSM120_boost    = nSel_boostv[ism_gf] + nSel_boostv[ism_vbf];     
  Double_t nTotalSM120Var_boost = 0.28*nSel_boostv[ism_gf] + 0.14*nSel_boostv[ism_vbf];
  Double_t nTotalSM120_vbf      = nSel_vbfv[ism_gf] + nSel_vbfv[ism_vbf];     
  Double_t nTotalSM120Var_vbf   = 0.26*nSel_vbfv[ism_gf] + 0.23*nSel_vbfv[ism_vbf];
  Double_t nTotalSM120_nob      = (nSel_nobv[imssm_gg] + nSel_nobv[imssm_bb]);
  Double_t nTotalSM120Var_nob   = 0.06*(nSel_nobv[imssm_gg] + nSel_nobv[imssm_bb]);
  Double_t nTotalSM120_b        = (nSel_bv[imssm_gg] + nSel_bv[imssm_bb]);     
  Double_t nTotalSM120Var_b     = 0.093*(nSel_bv[imssm_gg] + nSel_bv[imssm_bb]);


  texfile << " $H\\rightarrow \\tau\\tau$       & " << setw(5) << setprecision(0) << fixed << nTotalSM120_novbf << " & " << setw(5) << setprecision(0) << fixed << nTotalSM120Var_novbf << " & " << setw(5) << setprecision(0) << fixed << nTotalSM120_boost << " & " << setw(5) << setprecision(1) << fixed << nTotalSM120Var_boost << " & " << setw(5) << setprecision(0) << fixed << nTotalSM120_vbf << " & " << setw(5) << setprecision(1) << fixed << nTotalSM120Var_vbf << " & " << setw(5) << setprecision(0) << fixed << nTotalSM120_nob << " & "  << setw(5) << setprecision(0) << fixed << nTotalSM120Var_nob << " & " << setw(5) << setprecision(0) << fixed << nTotalSM120_b << " & " << setw(5) << setprecision(1) << fixed << nTotalSM120Var_b << "  \\\\";
  texfile << endl;


  texfile << " \\hline" << endl;

  texfile << " Data                          & \\multicolumn{2}{|c|}{" << setprecision(0) << fixed << nSel_novbfv[0] << "} & \\multicolumn{2}{|c|}{" << setprecision(0) << fixed << nSel_boostv[0] << "} & \\multicolumn{2}{|c|}{" << setprecision(0) << fixed << nSel_vbfv[0] << "} & \\multicolumn{2}{|c|}{" << setprecision(0) << fixed << nSel_nobv[0] << "} & \\multicolumn{2}{|c|}{" << setprecision(0) << fixed << nSel_bv[0] << "}  \\\\";
  texfile << endl;


  texfile << " \\hline" << endl;

  texfile << " \\multicolumn{6}{c}{ } \\\\" << endl;
  texfile << " \\multicolumn{3}{l}{Signal Efficiency } &  \\multicolumn{3}{c}{ } \\\\" << endl;

  texfile << " \\hline" << endl;

  texfile << " $gg\\rightarrow \\phi$           & \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_nobv[imssm_gg]/nUnskEventsv[imssm_gg] << "} & \\multicolumn{2}{|c|}{" <<  setprecision(2) << scientific << nSel_bv[imssm_gg]/nUnskEventsv[imssm_gg]  << "} \\\\";
  texfile << endl;
  texfile << " $bb\\rightarrow bb\\phi$         & \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_nobv[imssm_bb]/nUnskEventsv[imssm_bb] << "} & \\multicolumn{2}{|c|}{" <<  setprecision(2) << scientific << nSel_bv[imssm_bb]/nUnskEventsv[imssm_bb]  << "} \\\\";
  texfile << endl;
  texfile << " \\hline" << endl;
  texfile << " $gg\\rightarrow H$              & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_novbfv[ism_gf]/nUnskEventsv[ism_gf] << "} & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_boostv[ism_gf]/nUnskEventsv[ism_gf] << "} & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_vbfv[ism_gf]/nUnskEventsv[ism_gf] << "}& \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{-}  \\\\";
  texfile << endl;
  texfile << " $qq\\rightarrow H$              & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_novbfv[ism_vbf]/nUnskEventsv[ism_vbf] << "} & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_boostv[ism_vbf]/nUnskEventsv[ism_vbf] << "} & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_vbfv[ism_vbf]/nUnskEventsv[ism_vbf] << "}& \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{-}  \\\\";
  texfile << endl;
  texfile << " $qq\\rightarrow t\\bar{t}/VH$    & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_novbfv[ism_vtth]/nUnskEventsv[ism_vtth] << "} & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_boostv[ism_vtth]/nUnskEventsv[ism_vtth] << "} & \\multicolumn{2}{|c|}{" << setprecision(2) << scientific << nSel_vbfv[ism_vtth]/nUnskEventsv[ism_vtth] << "}& \\multicolumn{2}{|c|}{-} & \\multicolumn{2}{|c|}{-}  \\\\";
  texfile << endl;

  texfile << " \\hline" << endl;

  texfile << " \\end{tabular}" << endl;
  texfile << " \\caption{" << endl;
  texfile << " Number of expected and observed events in the event categories as described in the text and summarized in Tables~\\ref{tab:Event_Selection_SM} and~\\ref{tab:Event_Selection_MSSM} for the $e\\mu$-channel, for all considered backgrounds and for a SM Higgs boson with \\mbox{$m_{H}=120$~\\GeV}\\@. Also given are the signal acceptances for a MSSM Higgs boson with \\mbox{$m_{A}=120$~\\GeV} via gluon gluon fusion ($gg\\rightarrow \\phi$) and in association with $b$ quarks ($gg\\rightarrow \\phi bb$) and for a SM Higgs boson with \\mbox{$m_{H}=120$~\\GeV} via gluon gluon fusion ($gg\\rightarrow H$), \\emph{Vector Boson Fusion} ($qq\\rightarrow H$, VBF) and in association with a $t\\bar{t}$ pair or a vector boson ($qq\\rightarrow t\\bar{t}/VH$). All acceptances include the branching ratio into $\\tau\\tau$\\@." << endl;
  texfile << " }" << endl;
  texfile << " \\label{tab:cutflow_example}" << endl;
  texfile << " \\end{center}" << endl;
  texfile << " \\end{table}" << endl;

  texfile.close();


  makeHTML(outputDir);
  
  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  cout << endl;
  
  gBenchmark->Show("plotEmu");      
}


//=== FUNCTION DEFINITIONS ======================================================================================
//----------------------------------------------------------------------------------------
Double_t rescale(TH1F* hist, TString catname, TString histname, TFile* file){

  TH1F* hcenter = 0;
  
  TDirectory *dir = (TDirectory*)file->FindObjectAny(catname);
  hcenter = (TH1F*)(dir->Get(histname)); assert(hcenter);

  Double_t center = hcenter->Integral(0,hist->GetNbinsX()+1);
  Double_t integral = hist->Integral(0,hist->GetNbinsX()+1);

  Double_t scale = 1.0*center/integral;

  return scale;
}
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
  htmlfile << "<head><title>Htt->emu Plots</title></head>" << endl;
  htmlfile << "<body bgcolor=\"FFFFFF\">" << endl;
  htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">Htt->emu Plots</h3>" << endl;

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
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtele.png\"><img src=\"plots/mtele.png\" alt=\"plots/mtele.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtmu.png\"><img src=\"plots/mtmu.png\" alt=\"plots/mtmu.png\" width=\"100%\"></a></td>" << endl;
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
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/njets_log.png\"><img src=\"plots/njets_log.png\" alt=\"plots/njets_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nbjets_log.png\"><img src=\"plots/nbjets_log.png\" alt=\"plots/nbjets_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vbfjetpt1.png\"><img src=\"plots/vbfjetpt1.png\" alt=\"plots/vbfjetpt1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vbfjetpt2.png\"><img src=\"plots/vbfjetpt2.png\" alt=\"plots/vbfjetpt2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vbfjeteta1.png\"><img src=\"plots/vbfjeteta1.png\" alt=\"plots/vbfjeteta1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vbfjeteta2.png\"><img src=\"plots/vbfjeteta2.png\" alt=\"plots/vbfjeteta2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/boostjetpt.png\"><img src=\"plots/boostjetpt.png\" alt=\"plots/boostjetpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/boostjeteta.png\"><img src=\"plots/boostjeteta.png\" alt=\"plots/boostjeteta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/btag.png\"><img src=\"plots/btag.png\" alt=\"plots/btag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/btag_vbf.png\"><img src=\"plots/btag_vbf.png\" alt=\"plots/btag_vbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nvertices_reweighted.png\"><img src=\"plots/nvertices_reweighted.png\" alt=\"plots/nvertices_reweighted.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nvertices_raw.png\"><img src=\"plots/nvertices_raw.png\" alt=\"plots/nvertices_raw.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;   
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt1.png\"><img src=\"plots/pt1.png\" alt=\"plots/pt1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt2.png\"><img src=\"plots/pt2.png\" alt=\"plots/pt2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta1.png\"><img src=\"plots/eta1.png\" alt=\"plots/eta1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta2.png\"><img src=\"plots/eta2.png\" alt=\"plots/eta2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi1.png\"><img src=\"plots/phi1.png\" alt=\"plots/phi1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi2.png\"><img src=\"plots/phi2.png\" alt=\"plots/phi2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_lept.png\"><img src=\"plots/svfitmass_lept.png\" alt=\"plots/svfitmass_lept.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_lept-high.png\"><img src=\"plots/svfitmass_lept-high.png\" alt=\"plots/svfitmass_lept-high.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_incl.png\"><img src=\"plots/svfitmass_incl.png\" alt=\"plots/svfitmass_incl.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_incl_log.png\"><img src=\"plots/svfitmass_incl_log.png\" alt=\"plots/svfitmass_incl_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_novbf.png\"><img src=\"plots/svfitmass_class_novbf.png\" alt=\"plots/svfitmass_class_novbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_boost.png\"><img src=\"plots/svfitmass_class_boost.png\" alt=\"plots/svfitmass_class_boost.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_vbf.png\"><img src=\"plots/svfitmass_class_vbf.png\" alt=\"plots/svfitmass_class_vbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_nob.png\"><img src=\"plots/svfitmass_class_nob.png\" alt=\"plots/svfitmass_class_nob.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_b.png\"><img src=\"plots/svfitmass_class_b.png\" alt=\"plots/svfitmass_class_b.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_novbf_log.png\"><img src=\"plots/svfitmass_class_novbf_log.png\" alt=\"plots/svfitmass_class_novbf_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_boost_log.png\"><img src=\"plots/svfitmass_class_boost_log.png\" alt=\"plots/svfitmass_class_boost_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_vbf_log.png\"><img src=\"plots/svfitmass_class_vbf_log.png\" alt=\"plots/svfitmass_class_vbf_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_nob_log.png\"><img src=\"plots/svfitmass_class_nob_log.png\" alt=\"plots/svfitmass_class_nob_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_b_log.png\"><img src=\"plots/svfitmass_class_b_log.png\" alt=\"plots/svfitmass_class_b_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_lept.png\"><img src=\"plots/vismass_lept.png\" alt=\"plots/vismass_lept.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_lept-high.png\"><img src=\"plots/vismass_lept-high.png\" alt=\"plots/vismass_lept-high.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_incl.png\"><img src=\"plots/vismass_incl.png\" alt=\"plots/vismass_incl.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_incl_log.png\"><img src=\"plots/vismass_incl_log.png\" alt=\"plots/vismass_incl_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_novbf.png\"><img src=\"plots/vismass_class_novbf.png\" alt=\"plots/vismass_class_novbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_boost.png\"><img src=\"plots/vismass_class_boost.png\" alt=\"plots/vismass_class_boost.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_vbf.png\"><img src=\"plots/vismass_class_vbf.png\" alt=\"plots/vismass_class_vbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_nob.png\"><img src=\"plots/vismass_class_nob.png\" alt=\"plots/vismass_class_nob.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_b.png\"><img src=\"plots/vismass_class_b.png\" alt=\"plots/vismass_class_b.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_novbf_log.png\"><img src=\"plots/vismass_class_novbf_log.png\" alt=\"plots/vismass_class_novbf_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_boost_log.png\"><img src=\"plots/vismass_class_boost_log.png\" alt=\"plots/vismass_class_boost_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_vbf_log.png\"><img src=\"plots/vismass_class_vbf_log.png\" alt=\"plots/vismass_class_vbf_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_nob_log.png\"><img src=\"plots/vismass_class_nob_log.png\" alt=\"plots/vismass_class_nob_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_b_log.png\"><img src=\"plots/vismass_class_b_log.png\" alt=\"plots/vismass_class_b_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
            
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}
