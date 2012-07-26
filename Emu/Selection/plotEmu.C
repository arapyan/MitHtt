#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TRegexp.h>                // ROOT regexp class
#include <TEfficiency.h>            // efficiency class
#include <TMatrixD.h>               // matrix class
#include <TRandom3.h>               // ROOT math functions
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
#include "MitHtt/Utils/HttMVA.hh"          // MVA class
#include "MitHtt/Ntupler/interface/TSVfitter.h"
#include "MitHtt/Ntupler/interface/TSVSuperFitter.hh"
#include "CondFormats/EgammaObjects/interface/GBRForest.h"
#include "CondFormats/EgammaObjects/interface/GBRTree.h"

// define structure for output ntuple
#include "EmuData.hh"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

// rescale histograms with up/down shifts
Double_t rescale(TH1F* hist, TString catname, TString hname, TFile* file);

// generate web page
void makeHTML(const TString outDir);
TRandom3 randm(12345);

// Get higgs mass point from sample name
Int_t higgsmass(TString basename)
{ 
  stringstream ss(basename(TRegexp("[0-9][0-9]*"),3).Data());
  Int_t mass;
  ss >> mass;
  if(basename.Contains("gg-") || basename.Contains("bb-")) assert(mass>85 && mass<1200);
  return mass;
}

enum { kCenter, kDown, kUp };

//=== MAIN MACRO =================================================================================================

// run first with ecorr=0 to get central histograms, then with ecorr=1,2 to get the up,down scale variations

void plotEmu(const TString  conf,         // input file
             const TString  ntupleDir,    // directory of input ntuples
	     const TString  outputDir,    // output directory
             const TString  format,       // plot file format
	     const Double_t lumi,         // luminosity (pb^-1)
             const Int_t ecorr=0,         // energy scale shift: 0=center, 1=down, 2=up
             const Bool_t domssm=kFALSE
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


  // bins for limit input histograms
  const UInt_t massLNbins_SM = 13;
  const UInt_t massLNbins_MSSM = 19;

  Double_t massLEdges[massLNbins_SM]    = {0.,20.,40.,60.,80.,100.,120.,140.,170.,200.,250.,300.,350.};
  Double_t massLEdges_MSSM[massLNbins_MSSM]    = {0.,20.,40.,60.,80.,100.,120.,140.,170.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.};

  // get indices of samples
  UInt_t ifake=9999, izmm=9999, iztt=9999, ittbar=9999, imssm_gg=9999, imssm_bb=9999, ism_vbf=9999, ism_gf=9999, ism_vtth=9999, issfake=9999, iemb=9999, iewk=9999, iewk_7TeV=9999, ittbar_7TeV=9999;
  for(UInt_t isam=0; isam<snamev.size(); isam++) {
    if(snamev[isam].Contains("fakes") && !snamev[isam].Contains("ss-fakes")) ifake = isam;
    if(snamev[isam].Contains("ss-fakes")) issfake = isam;
    if(snamev[isam].Contains("zmm"))   izmm  = isam;
    if(snamev[isam].Contains("ztt"))   iztt = isam;
    if(snamev[isam].Contains("emb"))   iemb = isam;
    if(snamev[isam].Contains("ttbar-8TeV")) ittbar = isam;
    if(snamev[isam].Contains("ttbar-7TeV")) ittbar_7TeV = isam;
    if(snamev[isam].Contains("ewk-8TeV"))   iewk = isam;
    if(snamev[isam].Contains("ewk-7TeV"))   iewk_7TeV = isam;
    if(snamev[isam].Contains("htt_gg_mssm")) imssm_gg=isam;
    if(snamev[isam].Contains("htt_bb_mssm")) imssm_bb=isam;
    if(snamev[isam].Contains("htt_vbf_sm")) ism_vbf=isam;
    if(snamev[isam].Contains("htt_gf_sm")) ism_gf=isam;
    if(snamev[isam].Contains("htt_vtth_sm")) ism_vtth=isam;
  }

  //CPlot::sOutDir = outputDir + TString("/plots");
  const TString outDirName = outputDir + TString("/plots");
  CPlot::sOutDir = outDirName;
  gSystem->mkdir(outDirName,kTRUE);
  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };  // final state type

  const Double_t kssFakeWgt = 1.3; //1.11;  // ratio of fake-rate fakes to same-sign fakes

  //
  // scale factors for each category (used only for madgraph ztt)
  //
  const Double_t kCat_0jet = 1.00;
  const Double_t kCat_boost = 1.00;
  const Double_t kCat_b = 1.00;
  const Double_t kCat_2jet = 1.00;
  const Double_t kCat_vbf   = 1.00;

  // scale factors for normalizing 7TeV samples where used
  const Double_t ewk_0jet_lowpt_norm = 1.0037;
  const Double_t ewk_0jet_highpt_norm = 1.0971;
  const Double_t ewk_boost_highpt_norm = 1.5260;
  const Double_t ttbar_vbf_norm = 1.3369;

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
  //vector<TH2F*> hMassVpmissv, hProjVisVProjMetv;
  vector<TH1F*> hProjMetv, hProjVisv, hProjVarv,hRawProjVarv;
  vector<TH1F*> hNjetsv, hNbjetsv;
  vector<TH1F*> hBdiscrv, hBdiscr_vbfv;
  //vector<TH1F*> hDPhiv, hMtv, hPtv, hLepDEtav;
  //vector<TH1F*> hMetDPhiv;
  //vector<TH1F*> hPt1v, hEta1v, hPhi1v;   	// leading lepton
  //vector<TH1F*> hPt2v, hEta2v, hPhi2v;   	// trailing lepton
  vector<TH1F*> hPtMuv, hEtaMuv, hPhiMuv;  	// muon
  vector<TH1F*> hPtElev, hEtaElev, hPhiElev; 	// electron
  //vector<TH1F*> hD0Elev, hD0Muv;                // d0
  //vector<TH1F*> hMtElev, hMtMuv;                // mT
  //vector<TH1F*> hJetPt1v, hJetEta1v; 		// leading jet
  //vector<TH1F*> hJetPt2v, hJetEta2v; 		// second pt jet
  vector<TH1F*> hJetDPhiv, hJetDPhi_nobtagv;                    // dphi of first two jets
  vector<TH1F*> hMjjv, hDEtav, hEtaProdv;       // kinematics of first two jets
  vector<TH1F*> hPtjjv, hPtHv, hHDiJetDPhiv;    // inputs for VBF mva
  vector<TH1F*> hMjj_nobtagv, hDEta_nobtagv, hEtaProd_nobtagv;       // kinematics of first two jets
  vector<TH1F*> hPtjj_nobtagv, hPtH_nobtagv, hHDiJetDPhi_nobtagv;    // inputs for VBF mva
  vector<TH1F*> hBJetPtv, hBJetEtav, hBJetPhiv; // leading b-jet (can be same as leading two jets)
  vector<TH1F*> hNPVv, hNPVrawv;                // primary vertexes
  vector<TH1F*> hVBFJetPt1v, hVBFJetEta1v;      // leading jet
  vector<TH1F*> hVBFJetPt2v, hVBFJetEta2v;      // second pt jet
  vector<TH1F*> hBoostJetPtv, hBoostJetEtav;    // leading jet
  vector<TH1F*> hBDTv, hBDT_nobtagv;                          // VBF bdtg values
  vector<TH1F*> hMassv, hMassVisv, hMassLv, hMassHighv, hMassVisHighv;    // *L: plots for limits

  // inclusive
  vector<TH1F*> hMass_iv;
  vector<TH1F*> hMassVis_iv;
  vector<TH1F*> hMassL_iv;

  // Class 0-jet: no pt 30 jets
  vector<TH1F*> hMass_0jetv;
  vector<TH1F*> hMassVis_0jetv;
  vector<TH1F*> hMassL_0jetv;

  // Class 0-jet: no pt 30 jets
  vector<TH1F*> hMass_0jet_lowptv;
  vector<TH1F*> hMassVis_0jet_lowptv;
  vector<TH1F*> hMassL_0jet_lowptv;

  // Class 0-jet: no pt 30 jets
  vector<TH1F*> hMass_0jet_highptv;
  vector<TH1F*> hMassVis_0jet_highptv;
  vector<TH1F*> hMassL_0jet_highptv;

  // Class 1-jet, no btag: one pt 30 jet, no b-tagged jet
  vector<TH1F*> hMass_boostv;
  vector<TH1F*> hMassVis_boostv;
  vector<TH1F*> hMassL_boostv;

  // Class 1-jet, no btag: one pt 30 jet, no b-tagged jet
  vector<TH1F*> hMass_boost_lowptv;
  vector<TH1F*> hMassVis_boost_lowptv;
  vector<TH1F*> hMassL_boost_lowptv;

  // Class 1-jet, no btag: one pt 30 jet, no b-tagged jet
  vector<TH1F*> hMass_boost_highptv;
  vector<TH1F*> hMassVis_boost_highptv;
  vector<TH1F*> hMassL_boost_highptv;

  // Class 1-jet, btag: one pt 30 jet, one b-tagged jet
  vector<TH1F*> hMass_bv;
  vector<TH1F*> hMassVis_bv;
  vector<TH1F*> hMassL_bv;

  // Class 1-jet, btag: one pt 30 jet, one b-tagged jet
  vector<TH1F*> hMass_b_lowptv;
  vector<TH1F*> hMassVis_b_lowptv;
  vector<TH1F*> hMassL_b_lowptv;

  // Class 1-jet, btag: one pt 30 jet, one b-tagged jet
  vector<TH1F*> hMass_b_highptv;
  vector<TH1F*> hMassVis_b_highptv;
  vector<TH1F*> hMassL_b_highptv;

  // Class VH: two pt 30 jets, VH selection
  //vector<TH1F*> hMass_2jetv;
  //vector<TH1F*> hMassVis_2jetv;
  //vector<TH1F*> hMassL_2jetv;

  // Class vbf: two or more pt 30 jets and vbf cuts
  vector<TH1F*> hMass_vbfv;
  vector<TH1F*> hMassVis_vbfv;
  vector<TH1F*> hMassL_vbfv;

  vector<Double_t> nSelv,       nSelVarv;
  vector<Double_t> nSel_iv,     nSelVar_iv;
  vector<Double_t> nSel_0jetv,  nSelVar_0jetv, nSelpz_0jetv;
  vector<Double_t> nSel_0jet_lowptv,   nSelVar_0jet_lowptv;
  vector<Double_t> nSel_0jet_highptv,  nSelVar_0jet_highptv;
  vector<Double_t> nSel_boostv, nSelVar_boostv;
  vector<Double_t> nSel_boost_lowptv,  nSelVar_boost_lowptv;
  vector<Double_t> nSel_boost_highptv, nSelVar_boost_highptv;
  vector<Double_t> nSel_bv,     nSelVar_bv;
  vector<Double_t> nSel_b_lowptv,      nSelVar_b_lowptv;
  vector<Double_t> nSel_b_highptv,     nSelVar_b_highptv;
  vector<Double_t> nSel_2jetv,  nSelVar_2jetv;
  vector<Double_t> nSel_vbfv,   nSelVar_vbfv;
  
  vector<Double_t> nUnskEventsv;


  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    // after lepton selection
    sprintf(hname,"hMet_%i",isam);     hMetv.push_back(new TH1F(hname,"",30,0,150));        hMetv[isam]->Sumw2();
    sprintf(hname,"hMetRaw_%i",isam);  hMetRawv.push_back(new TH1F(hname,"",30,0,150));     hMetRawv[isam]->Sumw2();
    sprintf(hname,"hProjMet_%i",isam);       hProjMetv.push_back(new TH1F(hname,"",30,-100,120));     hProjMetv[isam]->Sumw2();
    sprintf(hname,"hProjVis_%i",isam);       hProjVisv.push_back(new TH1F(hname,"",20,0,50));         hProjVisv[isam]->Sumw2();
    //sprintf(hname,"hProjVisVProjMet_%i",isam);    hProjVisVProjMetv.push_back(new TH2F(hname,"",50,0,50,110,-100,120));                hProjVisVProjMetv[isam]->Sumw2();
    sprintf(hname,"hProjVar_%i",isam);       hProjVarv.push_back(new TH1F(hname,"",25,-200,100));     hProjVarv[isam]->Sumw2();
    sprintf(hname,"hRawProjVar_%i",isam);    hRawProjVarv.push_back(new TH1F(hname,"",25,-200,100));  hRawProjVarv[isam]->Sumw2();
    sprintf(hname,"hNjets_%i",isam);         hNjetsv.push_back(new TH1F(hname,"",11,-0.5,10.5));      hNjetsv[isam]->Sumw2();
    sprintf(hname,"hNbjets_%i",isam);        hNbjetsv.push_back(new TH1F(hname,"",6,-0.5,5.5));       hNbjetsv[isam]->Sumw2();
    sprintf(hname,"hBdiscr_%i",isam);        hBdiscrv.push_back(new TH1F(hname,"",20,0,1.0));          hBdiscrv[isam]->Sumw2();
    //sprintf(hname,"hBdiscr_vbf_%i",isam);    hBdiscr_vbfv.push_back(new TH1F(hname,"",10,0,1.0));      hBdiscr_vbfv[isam]->Sumw2();
    //sprintf(hname,"hDPhi_%i",isam);    hDPhiv.push_back(new TH1F(hname,"",36,0,180));       hDPhiv[isam]->Sumw2();
    //sprintf(hname,"hMt_%i",isam);      hMtv.push_back(new TH1F(hname,"",30,0,210));         hMtv[isam]->Sumw2();
    //sprintf(hname,"hPt_%i",isam);      hPtv.push_back(new TH1F(hname,"",30,0,120));         hPtv[isam]->Sumw2();
    //sprintf(hname,"hLepDEta_%i",isam); hLepDEtav.push_back(new TH1F(hname,"",20,0,4));      hLepDEtav[isam]->Sumw2();
    //sprintf(hname,"hMetDPhi_%i",isam); hMetDPhiv.push_back(new TH1F(hname,"",30,0,180));    hMetDPhiv[isam]->Sumw2();
    //sprintf(hname,"hPt1_%i",isam);     hPt1v.push_back(new TH1F(hname,"",30,0,150));        hPt1v[isam]->Sumw2();
    //sprintf(hname,"hEta1_%i",isam);    hEta1v.push_back(new TH1F(hname,"",20,-3,3));        hEta1v[isam]->Sumw2();
    //sprintf(hname,"hPhi1_%i",isam);    hPhi1v.push_back(new TH1F(hname,"",20,-3.2,3.2));    hPhi1v[isam]->Sumw2();  
    //sprintf(hname,"hPt2_%i",isam);     hPt2v.push_back(new TH1F(hname,"",60,0,150));        hPt2v[isam]->Sumw2();
    //sprintf(hname,"hEta2_%i",isam);    hEta2v.push_back(new TH1F(hname,"",20,-3,3));        hEta2v[isam]->Sumw2();
    //sprintf(hname,"hPhi2_%i",isam);    hPhi2v.push_back(new TH1F(hname,"",20,-3.2,3.2));    hPhi2v[isam]->Sumw2();
    sprintf(hname,"hPtMu_%i",isam);    hPtMuv.push_back(new TH1F(hname,"",20,0,100));       hPtMuv[isam]->Sumw2();
    sprintf(hname,"hEtaMu_%i",isam);   hEtaMuv.push_back(new TH1F(hname,"",20,-3,3));       hEtaMuv[isam]->Sumw2();
    sprintf(hname,"hPhiMu_%i",isam);   hPhiMuv.push_back(new TH1F(hname,"",20,-3.2,3.2));   hPhiMuv[isam]->Sumw2();  
    //sprintf(hname,"hD0Mu_%i",isam);    hD0Muv.push_back(new TH1F(hname,"",50,0,0.05));      hD0Muv[isam]->Sumw2();
    sprintf(hname,"hPtEle_%i",isam);   hPtElev.push_back(new TH1F(hname,"",20,0,100));      hPtElev[isam]->Sumw2();
    sprintf(hname,"hEtaEle_%i",isam);  hEtaElev.push_back(new TH1F(hname,"",20,-3,3));      hEtaElev[isam]->Sumw2();
    sprintf(hname,"hPhiEle_%i",isam);  hPhiElev.push_back(new TH1F(hname,"",20,-3.2,3.2));  hPhiElev[isam]->Sumw2();
    //sprintf(hname,"hD0Ele_%i",isam);   hD0Elev.push_back(new TH1F(hname,"",50,0,0.05));     hD0Elev[isam]->Sumw2();
    //sprintf(hname,"hMtEle_%i",isam);   hMtElev.push_back(new TH1F(hname,"",30,0,150));      hMtElev[isam]->Sumw2();
    //sprintf(hname,"hMtMu_%i",isam);    hMtMuv.push_back(new TH1F(hname,"",30,0,150));       hMtMuv[isam]->Sumw2();
    //sprintf(hname,"hJetPt1_%i",isam);  hJetPt1v.push_back(new TH1F(hname,"",30,0,300));     hJetPt1v[isam]->Sumw2();
    //sprintf(hname,"hJetEta1_%i",isam); hJetEta1v.push_back(new TH1F(hname,"",20,-5,5));     hJetEta1v[isam]->Sumw2();
    //sprintf(hname,"hJetPt2_%i",isam);  hJetPt2v.push_back(new TH1F(hname,"",30,0,300));     hJetPt2v[isam]->Sumw2();
    //sprintf(hname,"hJetEta2_%i",isam); hJetEta2v.push_back(new TH1F(hname,"",20,-5,5));     hJetEta2v[isam]->Sumw2();
    sprintf(hname,"hJetDPhi_%i",isam); hJetDPhiv.push_back(new TH1F(hname,"",15,0,180));    hJetDPhiv[isam]->Sumw2();
    sprintf(hname,"hJetDPhi_nobtag_%i",isam); hJetDPhi_nobtagv.push_back(new TH1F(hname,"",9,0,180));    hJetDPhi_nobtagv[isam]->Sumw2();
    sprintf(hname,"hMjj_%i",isam);     hMjjv.push_back(new TH1F(hname,"",6,0,1200));       hMjjv[isam]->Sumw2();
    sprintf(hname,"hPtjj_%i",isam);    hPtjjv.push_back(new TH1F(hname,"",20,0,500));       hPtjjv[isam]->Sumw2();
    sprintf(hname,"hPtH_%i",isam);     hPtHv.push_back(new TH1F(hname,"",20,0,500));        hPtHv[isam]->Sumw2();
    sprintf(hname,"hHDiJetDPhi_%i",isam); hHDiJetDPhiv.push_back(new TH1F(hname,"",15,0,180));    hHDiJetDPhiv[isam]->Sumw2();
    sprintf(hname,"hDEta_%i",isam);    hDEtav.push_back(new TH1F(hname,"",10,0,10));         hDEtav[isam]->Sumw2();
    sprintf(hname,"hEtaProd_%i",isam); hEtaProdv.push_back(new TH1F(hname,"",30,-7.5,7.5)); hEtaProdv[isam]->Sumw2();
    sprintf(hname,"hMjj_nobtag_%i",isam);     hMjj_nobtagv.push_back(new TH1F(hname,"",10,0,1000));       hMjj_nobtagv[isam]->Sumw2();
    sprintf(hname,"hPtjj_nobtag_%i",isam);    hPtjj_nobtagv.push_back(new TH1F(hname,"",10,0,500));       hPtjj_nobtagv[isam]->Sumw2();
    sprintf(hname,"hPtH_nobtag_%i",isam);     hPtH_nobtagv.push_back(new TH1F(hname,"",10,0,500));        hPtH_nobtagv[isam]->Sumw2();
    sprintf(hname,"hDEta_nobtag_%i",isam);    hDEta_nobtagv.push_back(new TH1F(hname,"",10,0,8));         hDEta_nobtagv[isam]->Sumw2();
    sprintf(hname,"hHDiJetDPhi_nobtag_%i",isam); hHDiJetDPhi_nobtagv.push_back(new TH1F(hname,"",9,0,180));    hHDiJetDPhi_nobtagv[isam]->Sumw2();
    sprintf(hname,"hEtaProd_nobtag_%i",isam); hEtaProd_nobtagv.push_back(new TH1F(hname,"",10,-7.5,7.5)); hEtaProd_nobtagv[isam]->Sumw2();
    sprintf(hname,"hBJetPt_%i",isam);  hBJetPtv.push_back(new TH1F(hname,"",30,0,150));     hBJetPtv[isam]->Sumw2();
    sprintf(hname,"hBJetEta_%i",isam); hBJetEtav.push_back(new TH1F(hname,"",20,-3,3));     hBJetEtav[isam]->Sumw2();
    sprintf(hname,"hBJetPhi_%i",isam); hBJetPhiv.push_back(new TH1F(hname,"",20,-3.2,3.2)); hBJetPhiv[isam]->Sumw2();
    sprintf(hname,"hNPV_%i",isam);     hNPVv.push_back(new TH1F(hname,"",35,-0.5,34.5));    hNPVv[isam]->Sumw2();
    sprintf(hname,"hNPVraw_%i",isam);  hNPVrawv.push_back(new TH1F(hname,"",35,-0.5,34.5)); hNPVrawv[isam]->Sumw2();
    sprintf(hname,"hVBFJetPt1_%i",isam);     hVBFJetPt1v.push_back(new TH1F(hname,"",27,30,300));     hVBFJetPt1v[isam]->Sumw2();
    sprintf(hname,"hVBFJetEta1_%i",isam);    hVBFJetEta1v.push_back(new TH1F(hname,"",18,-4.5,4.5));  hVBFJetEta1v[isam]->Sumw2();
    sprintf(hname,"hVBFJetPt2_%i",isam);     hVBFJetPt2v.push_back(new TH1F(hname,"",27,30,300));     hVBFJetPt2v[isam]->Sumw2();
    sprintf(hname,"hVBFJetEta2_%i",isam);    hVBFJetEta2v.push_back(new TH1F(hname,"",18,-4.5,4.5));  hVBFJetEta2v[isam]->Sumw2();
    sprintf(hname,"hBoostJetPt_%i",isam);    hBoostJetPtv.push_back(new TH1F(hname,"",27,30,300));    hBoostJetPtv[isam]->Sumw2();
    sprintf(hname,"hBoostJetEta_%i",isam);   hBoostJetEtav.push_back(new TH1F(hname,"",18,-4.5,4.5)); hBoostJetEtav[isam]->Sumw2();
    sprintf(hname,"hBDT_%i",isam);     hBDTv.push_back(new TH1F(hname,"",25,-1.0,1.0));    hBDTv[isam]->Sumw2();
    sprintf(hname,"hBDT_nobtag_%i",isam);     hBDT_nobtagv.push_back(new TH1F(hname,"",25,-1.0,1.0));    hBDT_nobtagv[isam]->Sumw2();

    // lepton selection
    sprintf(hname,"hMass_%i",isam);          hMassv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));                                  hMassv[isam]->Sumw2();
    sprintf(hname,"hMassVis_%i",isam);       hMassVisv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));                               hMassVisv[isam]->Sumw2();
    sprintf(hname,"hMassHigh_%i",isam);      hMassHighv.push_back(new TH1F(hname,"",100,0,2000));                            hMassHighv[isam]->Sumw2();
    sprintf(hname,"hMassVisHigh_%i",isam);   hMassVisHighv.push_back(new TH1F(hname,"",100,0,1000));                         hMassVisHighv[isam]->Sumw2();
    sprintf(hname,"hMassL_%i",isam);         hMassLv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));           hMassLv[isam]->Sumw2();
    //sprintf(hname,"hMassVpmiss_%i",isam);    hMassVpmissv.push_back(new TH2F(hname,"",75,0,750,50,-100,100));                hMassVpmissv[isam]->Sumw2();
    // inclusive
    sprintf(hname,"hMass_i_%i",isam);        hMass_iv.push_back(new TH1F(hname,"",30,0,300));          hMass_iv[isam]->Sumw2();
    sprintf(hname,"hMassVis_i_%i",isam);     hMassVis_iv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));       hMassVis_iv[isam]->Sumw2();
    sprintf(hname,"hMassL_i_%i",isam);       hMassL_iv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));         hMassL_iv[isam]->Sumw2();
    // 0-jet
    sprintf(hname,"hMass_0jet_%i",isam);    hMass_0jetv.push_back(new TH1F(hname,"",20,0,400));      hMass_0jetv[isam]->Sumw2();
    sprintf(hname,"hMassVis_0jet_%i",isam); hMassVis_0jetv.push_back(new TH1F(hname,"",20,0,400));   hMassVis_0jetv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_0jet_%i",isam);   hMassL_0jetv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_0jetv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_0jet_%i",isam);   hMassL_0jetv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_0jetv[isam]->Sumw2();}
    // 0-jet, low pT
    sprintf(hname,"hMass_0jet_lowpt_%i",isam);    hMass_0jet_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));      hMass_0jet_lowptv[isam]->Sumw2();
    sprintf(hname,"hMassVis_0jet_lowpt_%i",isam); hMassVis_0jet_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));   hMassVis_0jet_lowptv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_0jet_lowpt_%i",isam);   hMassL_0jet_lowptv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_0jet_lowptv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_0jet_lowpt_%i",isam);   hMassL_0jet_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_0jet_lowptv[isam]->Sumw2();}
    // 0-jet, high pT
    sprintf(hname,"hMass_0jet_highpt_%i",isam);    hMass_0jet_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));      hMass_0jet_highptv[isam]->Sumw2();
    sprintf(hname,"hMassVis_0jet_highpt_%i",isam); hMassVis_0jet_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));   hMassVis_0jet_highptv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_0jet_highpt_%i",isam);   hMassL_0jet_highptv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_0jet_highptv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_0jet_highpt_%i",isam);   hMassL_0jet_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_0jet_highptv[isam]->Sumw2();}
    // 1-jet
    sprintf(hname,"hMass_boost_%i",isam);    hMass_boostv.push_back(new TH1F(hname,"",20,0,400));      hMass_boostv[isam]->Sumw2();
    sprintf(hname,"hMassVis_boost_%i",isam); hMassVis_boostv.push_back(new TH1F(hname,"",20,0,400));   hMassVis_boostv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_boost_%i",isam);   hMassL_boostv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_boostv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_boost_%i",isam);   hMassL_boostv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_boostv[isam]->Sumw2();}
    // 1-jet, low pT
    sprintf(hname,"hMass_boost_lowpt_%i",isam);    hMass_boost_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));      hMass_boost_lowptv[isam]->Sumw2();
    sprintf(hname,"hMassVis_boost_lowpt_%i",isam); hMassVis_boost_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));   hMassVis_boost_lowptv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_boost_lowpt_%i",isam);   hMassL_boost_lowptv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_boost_lowptv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_boost_lowpt_%i",isam);   hMassL_boost_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_boost_lowptv[isam]->Sumw2();}
    // 1-jet, high pT
    sprintf(hname,"hMass_boost_highpt_%i",isam);    hMass_boost_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));      hMass_boost_highptv[isam]->Sumw2();
    sprintf(hname,"hMassVis_boost_highpt_%i",isam); hMassVis_boost_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));   hMassVis_boost_highptv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_boost_highpt_%i",isam);   hMassL_boost_highptv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_boost_highptv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_boost_highpt_%i",isam);   hMassL_boost_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_boost_highptv[isam]->Sumw2();}
    // 1-jet, btag
    sprintf(hname,"hMass_b_%i",isam);    hMass_bv.push_back(new TH1F(hname,"",20,0,400));      hMass_bv[isam]->Sumw2();
    sprintf(hname,"hMassVis_b_%i",isam); hMassVis_bv.push_back(new TH1F(hname,"",20,0,400));   hMassVis_bv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_b_%i",isam);   hMassL_bv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_bv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_b_%i",isam);   hMassL_bv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_bv[isam]->Sumw2();}
    // 1-jet, btag, low pT
    sprintf(hname,"hMass_b_lowpt_%i",isam);    hMass_b_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));      hMass_b_lowptv[isam]->Sumw2();
    sprintf(hname,"hMassVis_b_lowpt_%i",isam); hMassVis_b_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));   hMassVis_b_lowptv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_b_lowpt_%i",isam);   hMassL_b_lowptv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_b_lowptv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_b_lowpt_%i",isam);   hMassL_b_lowptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_b_lowptv[isam]->Sumw2();}
    // 1-jet, btag, high pT
    sprintf(hname,"hMass_b_highpt_%i",isam);    hMass_b_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));      hMass_b_highptv[isam]->Sumw2();
    sprintf(hname,"hMassVis_b_highpt_%i",isam); hMassVis_b_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));   hMassVis_b_highptv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_b_highpt_%i",isam);   hMassL_b_highptv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_b_highptv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_b_highpt_%i",isam);   hMassL_b_highptv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_b_highptv[isam]->Sumw2();}
    // VH
    //sprintf(hname,"hMass_2jet_%i",isam);    hMass_2jetv.push_back(new TH1F(hname,"",20,0,400));      hMass_2jetv[isam]->Sumw2();
    //sprintf(hname,"hMassVis_2jet_%i",isam); hMassVis_2jetv.push_back(new TH1F(hname,"",20,0,400));   hMassVis_2jetv[isam]->Sumw2();
    //sprintf(hname,"hMassL_2jet_%i",isam);   hMassL_2jetv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_2jetv[isam]->Sumw2();
    // VBF
    sprintf(hname,"hMass_vbf_%i",isam);    hMass_vbfv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));      hMass_vbfv[isam]->Sumw2();
    sprintf(hname,"hMassVis_vbf_%i",isam); hMassVis_vbfv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));   hMassVis_vbfv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_vbf_%i",isam);   hMassL_vbfv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_vbfv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_vbf_%i",isam);   hMassL_vbfv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_vbfv[isam]->Sumw2();}

    nSelv.push_back(0);       nSelVarv.push_back(0);
    nSel_iv.push_back(0);     nSelVar_iv.push_back(0);
    nSel_0jetv.push_back(0);  nSelVar_0jetv.push_back(0);
    nSel_0jet_lowptv.push_back(0);    nSelVar_0jet_lowptv.push_back(0);
    nSel_0jet_highptv.push_back(0);   nSelVar_0jet_highptv.push_back(0);
    nSel_boostv.push_back(0); nSelVar_boostv.push_back(0);
    nSel_boost_lowptv.push_back(0);   nSelVar_boost_lowptv.push_back(0);
    nSel_boost_highptv.push_back(0);  nSelVar_boost_highptv.push_back(0);
    nSel_bv.push_back(0);     nSelVar_bv.push_back(0);
    nSel_b_lowptv.push_back(0);       nSelVar_b_lowptv.push_back(0);
    nSel_b_highptv.push_back(0);      nSelVar_b_highptv.push_back(0);
    nSel_2jetv.push_back(0);  nSelVar_2jetv.push_back(0);
    nSel_vbfv.push_back(0);   nSelVar_vbfv.push_back(0);

    nUnskEventsv.push_back(0);

  }

  UInt_t npt20jets;
  const UInt_t kMaxPt20Jets=50;
  
  TArrayF *btagArray = new TArrayF;
  TArrayF *jptArray = new TArrayF;
  TArrayF *jetaArray = new TArrayF;
  btagArray->Set(kMaxPt20Jets);
  jptArray->Set(kMaxPt20Jets);
  jetaArray->Set(kMaxPt20Jets);
  
  // Access samples and fill histograms
  //HttMVA *vbfMVA = new HttMVA();  
  //vbfMVA->Initialize("BDTG method", "/home/vdutta/cms/cmssw/new/CMSSW_4_4_1/src/MitHtt/Emu/Selection/MVAVBF_v3/weights/TMVA_BDTG.weights.xml", HttMVA::kVBF2);   // vbf mva

  //mithep::TSVfitter *fitter = new mithep::TSVfitter();

  TFile *infile=0;
  TTree *eventTree=0;

  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if((isam==0) && !hasData) continue;
    if(!domssm && (snamev[isam].Contains("_mssm_") && !snamev[isam].Contains("htt_"))) continue;
    if(domssm && (snamev[isam].Contains("_sm_") && !snamev[isam].Contains("htt_"))) continue;
    if(ecorr !=0 && (isam==iewk || isam==ittbar || isam==iewk_7TeV || isam==ittbar_7TeV || isam==izmm || isam==ifake || isam==issfake)) continue;

    const TString fname = ntupleDir + TString("/") + snamev[isam] + TString("_select.root");
    cout << "Processing " << fname << "..." << endl;   
    infile = new TFile(fname);
    assert(infile);

    // Get the TTree and set branch address
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    int   lRun         = 0;    eventTree->SetBranchAddress("run"        ,&lRun           );//Run
    int   lLumi        = 0;    eventTree->SetBranchAddress("lumi"       ,&lLumi          );//Lumi
    int   lEvt         = 0;    eventTree->SetBranchAddress("evt"        ,&lEvt           );//Evt
    int   lNPV         = 0;    eventTree->SetBranchAddress("npv"        ,&lNPV           );//NPV
    int   lNPU         = 0;    eventTree->SetBranchAddress("npu"        ,&lNPU           );//NPU
    float lPUWeight    = 0;    eventTree->SetBranchAddress("puweight"   ,&lPUWeight      );//Pielup Weight
    float lWeight      = 0;    eventTree->SetBranchAddress("weight"     ,&lWeight        );//mcweight*puweight*effweight
    float lMGen        = 0;    eventTree->SetBranchAddress("mass"       ,&lMGen          );//gen mass for ggh, bbh
    float lMSV         = 0;    eventTree->SetBranchAddress("m_sv"       ,&lMSV           );//SV Fit using integration method
    float lMVis        = 0;    eventTree->SetBranchAddress("m_vis"      ,&lMVis          );//visible mass
    float lPt1         = 0;    eventTree->SetBranchAddress("pt_1"       ,&lPt1           ); //pT 
    float lPhi1        = 0;    eventTree->SetBranchAddress("phi_1"      ,&lPhi1          ); //Phi 
    float lEta1        = 0;    eventTree->SetBranchAddress("eta_1"      ,&lEta1          ); //Eta 
    float lM1          = 0;    eventTree->SetBranchAddress("m_1"        ,&lM1            ); //Mass 
    float lMt1         = 0;    eventTree->SetBranchAddress("mt_1"       ,&lMt1           );//mT of  first lepton wrt to MVA met
    float lPt2         = 0;    eventTree->SetBranchAddress("pt_2"       ,&lPt2           );//pT
    float lPhi2        = 0;    eventTree->SetBranchAddress("phi_2"      ,&lPhi2          );//Phi
    float lEta2        = 0;    eventTree->SetBranchAddress("eta_2"      ,&lEta2          );//Eta
    float lM2          = 0;    eventTree->SetBranchAddress("m_2"        ,&lM2            ); //Mass 
    float lMt2         = 0;    eventTree->SetBranchAddress("mt_2"       ,&lMt2           );//mT of 2nd lepton wrt to MVA met
    float lMet         = 0;    eventTree->SetBranchAddress("met"        ,&lMet           ); //pfmet
    float lMetPhi      = 0;    eventTree->SetBranchAddress("metphi"     ,&lMetPhi        ); //pfmet Phi
    float lMVAMet      = 0;    eventTree->SetBranchAddress("mvamet"     ,&lMVAMet        ); //mvamet
    float lMVAMetPhi   = 0;    eventTree->SetBranchAddress("mvametphi"  ,&lMVAMetPhi     ); //mvamet Phi
    float lPZetaVis    = 0;    eventTree->SetBranchAddress("pzetavis"   ,&lPZetaVis      ); //pZeta Visible
    float lPZetaMiss   = 0;    eventTree->SetBranchAddress("pzetamiss"  ,&lPZetaMiss     ); //pZeta Missing
    float lMetCov00    = 0;    eventTree->SetBranchAddress("metcov00"   ,&lMetCov00      ); //pf met covariance matrix 00 
    float lMetCov01    = 0;    eventTree->SetBranchAddress("metcov01"   ,&lMetCov01      ); //pf met covariance matrix 01 
    float lMetCov10    = 0;    eventTree->SetBranchAddress("metcov10"   ,&lMetCov10      ); //pf met covariance matrix 10 
    float lMetCov11    = 0;    eventTree->SetBranchAddress("metcov11"   ,&lMetCov11      ); //pf met covariance matrix 11 
    float lJPt1        = 0;    eventTree->SetBranchAddress("jpt_1"      ,&lJPt1          );//Jet Pt after corrections
    float lJEta1       = 0;    eventTree->SetBranchAddress("jeta_1"     ,&lJEta1         );//Jet Eta
    float lJPhi1       = 0;    eventTree->SetBranchAddress("jphi_1"     ,&lJPhi1         );//Jet Phi     
    float lJPt2        = 0;    eventTree->SetBranchAddress("jpt_2"      ,&lJPt2          );//Jet Pt after corrections
    float lJEta2       = 0;    eventTree->SetBranchAddress("jeta_2"     ,&lJEta2         );//Jet Eta
    float lJPhi2       = 0;    eventTree->SetBranchAddress("jphi_2"     ,&lJPhi2         );//Jet Phi
    float lBTagPt      = 0;    eventTree->SetBranchAddress("bpt"        ,&lBTagPt        );//Corrected BTag Pt
    float lBTagEta     = 0;    eventTree->SetBranchAddress("beta"       ,&lBTagEta       );//Btag Eta
    float lBTagPhi     = 0;    eventTree->SetBranchAddress("bphi"       ,&lBTagPhi       );//Btag Phi
    float lMJJ         = 0;    eventTree->SetBranchAddress("mjj"        ,&lMJJ           );//Mass Di Jet system  
    float lJDEta       = 0;    eventTree->SetBranchAddress("jdeta"      ,&lJDEta         );//|jeta_1-jeta_2| 
    float lVBFMVA      = 0;    eventTree->SetBranchAddress("mva"        ,&lVBFMVA        );//VBFMVA value
    float lJDPhi       = 0;    eventTree->SetBranchAddress("jdphi"      ,&lJDPhi         );//Delta Phi between two leading jets
    float lDiJetPt     = 0;    eventTree->SetBranchAddress("dijetpt"    ,&lDiJetPt       );//Pt of the di jet system
    float lDiJetPhi    = 0;    eventTree->SetBranchAddress("dijetphi"   ,&lDiJetPhi      );//Phi of the di jet system
    float lHDJetPhi    = 0;    eventTree->SetBranchAddress("hdijetphi"  ,&lHDJetPhi      );//Phi of the di jet system - Higgs system phi
    float lVisJetEta   = 0;    eventTree->SetBranchAddress("visjeteta"  ,&lVisJetEta     );
    float lPtVis       = 0;    eventTree->SetBranchAddress("ptvis"      ,&lPtVis         );//Pt Vis
    float lPtH         = 0;    eventTree->SetBranchAddress("pth"        ,&lPtH           );//Pt of the higgs system
    int   lNBTag       = 0;    eventTree->SetBranchAddress("nbtag"      ,&lNBTag         );
    int   lNJets       = 0;    eventTree->SetBranchAddress("njets"      ,&lNJets         );
    int   lNJetInGap   = 0;    eventTree->SetBranchAddress("njetingap"  ,&lNJetInGap     ); //number of gap jets
    float lGenPt1      = 0;    eventTree->SetBranchAddress("genlpt_1"   ,&lGenPt1        ); //pT 
    float lGenPt2      = 0;    eventTree->SetBranchAddress("genlpt_2"   ,&lGenPt2        ); //pT 

    // extra branches
    eventTree->SetBranchAddress("npt20jets",&npt20jets);
    eventTree->SetBranchAddress("btagArray",&btagArray);
    eventTree->SetBranchAddress("jptArray",&jptArray);
    eventTree->SetBranchAddress("jetaArray",&jetaArray);

    assert(samplev[isam]->fnamev.size()>0);
    nUnskEventsv[isam] = lumi*samplev[isam]->xsecv[0];

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);

      Double_t wgt = 1.0;
      if(isam > 0 && isam!=ifake && !snamev[isam].Contains("ss-fakes") && isam!=iemb)
	wgt = lWeight*lumi;
      if(snamev[isam].Contains("ss-fakes")) // same-sign fakes must have weight 1, but normalize them to fake-rate fakes
	wgt = kssFakeWgt;
      else if (isam==ifake) // fake-rate fakes have their weight stored in the ntuple
	wgt = lWeight;
      if(isam==iemb) {
	wgt=lWeight*lumi*0.081769112; // normalized using madgraph
      }
      if(isam==iztt) wgt *= 0.98;
      if(isam==ittbar && lNJets==0) wgt *= 0.97;
      if(isam==ittbar && lNJets==1) wgt *= 0.92;
      if(isam==ittbar && lNJets>=2) wgt *= 0.88;

      if(domssm && (snamev[isam].Contains("mssm") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH"))) {
	if (lMGen < 0.7*higgsmass(snamev[isam]) || lMGen > 1.3*higgsmass(snamev[isam])) continue;
      }

      Double_t vbfMvaValue = -1.1; 
      if(lMJJ>0) vbfMvaValue = lVBFMVA; //vbfMVA->MVAValue(lMJJ, lJDEta, lJDPhi, lDiJetPt, lPtH, lHDJetPhi, lVisJetEta, lPtVis);

      if(isam!=issfake) {
        if(lMJJ>0 && (0.85*lPZetaVis - lPZetaMiss <= 25)) {
	  if(isam>0 || (isam==0 && vbfMvaValue < 0.5)) {
	    hBDTv[isam]          ->Fill(vbfMvaValue,  wgt);
            if(lNBTag==0) hBDT_nobtagv[isam]          ->Fill(vbfMvaValue,  wgt);
	  }
	}
      }

      TVector3 m,e,metv;
      m.SetPtEtaPhi(lPt2,0,lPhi2);
      e.SetPtEtaPhi(lPt1,0,lPhi1);
      metv.SetPtEtaPhi(lMet,0,lMetPhi);
      TVector3 bisector(m.Unit() + e.Unit());
      bisector = bisector.Unit();
      Double_t pfProjMet  =  metv.Dot(bisector);

      // svfit
      TMatrixD covM(2,2);
      covM(0,0) = lMetCov00; //lMVACov00;
      covM(0,1) = lMetCov01; //lMVACov01;
      covM(1,0) = lMetCov10; //lMVACov10;
      covM(1,1) = lMetCov11; //lMVACov11;

      TLorentzVector     lv1;  lv1.SetPtEtaPhiM(lPt1,lEta1,lPhi1,lM1);
      TLorentzVector     lv2;  lv2.SetPtEtaPhiM(lPt2,lEta2,lPhi2,lM2);
      mithep::FourVectorM dau1, dau2;
      dau1.SetPxPyPzE(lv1.Px(),lv1.Py(),lv1.Pz(),lv1.E()); dau1.SetM(lv1.M());
      dau2.SetPxPyPzE(lv2.Px(),lv2.Py(),lv2.Pz(),lv2.E()); dau2.SetM(lv2.M());

      Double_t svfmass = 0;
      //TLorentzVector svfit = fitter->fit(covM, dau1, dau2, lMet, lMetPhi);
      //svfmass = svfit.M();

      // fill plots after lepton selection
      if(isam!=ifake) {
	hMetv[isam]        ->Fill(lMVAMet,    	 wgt);
	hMetRawv[isam]     ->Fill(lMet,      	 wgt);
	hProjMetv[isam]    ->Fill(lPZetaMiss,     	 wgt);
	hProjVisv[isam]    ->Fill(lPZetaVis,     	 wgt);
	//hProjVisVProjMetv[isam]->Fill(lPZetaVis,lPZetaMiss,wgt);
	hProjVarv[isam]    ->Fill(-0.85*lPZetaVis + lPZetaMiss,  wgt);
	//hRawProjVarv[isam] ->Fill(-0.85*lPZetaVis + pfProjMet,  wgt);
	hNjetsv[isam]      ->Fill(lNJets,      wgt);
	hNbjetsv[isam]     ->Fill(lNBTag,     wgt);

	assert(npt20jets<kMaxPt20Jets);
	for(UInt_t ib=0;ib<npt20jets;ib++)
	  hBdiscrv[isam]   ->Fill((*btagArray)[ib],  wgt);

	//hDPhiv[isam]   ->Fill(toolbox::deltaPhi(lPhi1,lPhi2)*180./pi, wgt);
	//hMtv[isam]     ->Fill(data.mt,      wgt);
	//hPtv[isam]     ->Fill(lPtVis,       wgt);
	//hLepDEtav[isam]->Fill(lEta1-lEta2,  wgt);
	//hMetDPhiv[isam]->Fill(toolbox::deltaPhi(data.phi,lMVAMetPhi)*180./pi, wgt);
	//hPt1v[isam]    ->Fill((lPt1 > lPt2 ? lPt1 : lPt2),    wgt);
	//hEta1v[isam]   ->Fill((lPt1 > lPt2 ? lEta1 : lEta2),   wgt);
	//hPhi1v[isam]   ->Fill((lPt1 > lPt2 ? lPhi1 : lPhi2),   wgt);
	//hPt2v[isam]    ->Fill((lPt1 > lPt2 ? lPt2 : lPt1),    wgt);
	//hEta2v[isam]   ->Fill((lPt1 > lPt2 ? lEta2 : lEta1),   wgt);
	//hPhi2v[isam]   ->Fill((lPt1 > lPt2 ? lPhi2 : lPhi1),   wgt);
	hPtMuv[isam]   ->Fill(lPt2,         wgt);
	hEtaMuv[isam]  ->Fill(lEta2,        wgt);
	hPhiMuv[isam]  ->Fill(lPhi2,        wgt);
        //hD0Muv[isam]   ->Fill(lD02,         wgt);
	hPtElev[isam]  ->Fill(lPt1,         wgt);
	hEtaElev[isam] ->Fill(lEta1,        wgt);
	hPhiElev[isam] ->Fill(lPhi1,        wgt);
        //hD0Elev[isam]  ->Fill(lD01,         wgt);
	//hMtElev[isam]  ->Fill(lMt1,         wgt);
	//hMtMuv[isam]   ->Fill(lMt2,         wgt);
	if(lNJets>0 && (0.85*lPZetaVis - lPZetaMiss <= 25)) {
	  //hJetPt1v[isam]	->Fill(lJPt1,     wgt);
	  //hJetEta1v[isam]	->Fill(lJEta1,    wgt);
	  hBoostJetPtv[isam]	->Fill(lJPt1,     wgt);
	  hBoostJetEtav[isam]	->Fill(lJEta1,    wgt);
	}
	if(lNJets>1 && (0.85*lPZetaVis - lPZetaMiss <= 25)) {
	  //hJetPt2v[isam]	->Fill(lJPt2,     wgt);
	  //hJetEta2v[isam]	->Fill(lJEta2,    wgt);
	  hVBFJetPt1v[isam]	->Fill(lJPt1,     wgt);
	  hVBFJetEta1v[isam]	->Fill(lJEta1,    wgt);
	  hVBFJetPt2v[isam]	->Fill(lJPt2,     wgt);
	  hVBFJetEta2v[isam]	->Fill(lJEta2,    wgt);
	}
      }
      if(isam!=issfake) {
	if(lMJJ>0 && (0.85*lPZetaVis - lPZetaMiss <= 25)) {
          if(lNJetInGap==0 && vbfMvaValue>0.5 && lNBTag==0) {
	    hJetDPhiv[isam]	->Fill(toolbox::deltaPhi(lJPhi1,lJPhi2)*180./pi,    wgt);
	    hMjjv[isam]		->Fill(lMJJ,      wgt);
	    hPtjjv[isam]		->Fill(lDiJetPt,  wgt);
	    hPtHv[isam]		->Fill(lPtH,      wgt);
	    hHDiJetDPhiv[isam]	->Fill(lHDJetPhi*180./pi, wgt);
	    hDEtav[isam]		->Fill(fabs(lJEta1 - lJEta2),    wgt);
	    hEtaProdv[isam]	->Fill(lJEta1*lJEta2,            wgt);
	  }
          if(lNBTag==0) {
            hJetDPhi_nobtagv[isam]      ->Fill(lJDPhi*180./pi,    wgt);
            hMjj_nobtagv[isam]          ->Fill(lMJJ,      wgt);
            hPtjj_nobtagv[isam]         ->Fill(lDiJetPt,  wgt);
            hPtH_nobtagv[isam]          ->Fill(lPtH,      wgt);
            hHDiJetDPhi_nobtagv[isam]   ->Fill(lHDJetPhi*180./pi, wgt);
            hDEta_nobtagv[isam]         ->Fill(lJDEta,    wgt);
            hEtaProd_nobtagv[isam]      ->Fill(lJEta1*lJEta2,            wgt);
	  }
	}
      }
      if(isam!=ifake) {
	if(lBTagPt>0) {
	  hBJetPtv[isam]	->Fill(lBTagPt,     wgt);
	  hBJetEtav[isam]	->Fill(lBTagEta,    wgt);
	  hBJetPhiv[isam]	->Fill(lBTagPhi,    wgt);
	}
	hNPVv[isam]    ->Fill(lNPV,     wgt);
	if(lPUWeight!=0)
	  hNPVrawv[isam] ->Fill(lNPV,    (isam<1 || snamev[isam].Contains("fake")) ? wgt : wgt/lPUWeight);

	hMassv[isam]       ->Fill(lMSV,   	 wgt);
	hMassVisv[isam]    ->Fill(lMVis,     wgt);
	hMassLv[isam]      ->Fill(lMSV,   	 wgt);      
	hMassHighv[isam]   ->Fill(lMSV,   	 wgt);
	hMassVisHighv[isam]->Fill(lMVis,     wgt);
	//hMassVpmissv[isam] ->Fill(lMSV, lPZetaMiss, wgt);

	nSelv[isam]    += wgt;
	nSelVarv[isam] += wgt*wgt;
      }

      Bool_t bjets = lNBTag>0;
      Bool_t vbfcuts = lNJets>=2 && (lMJJ>mjjMin) && (fabs(lJEta1-lJEta2)>dEtaMin) && lNJetInGap==0;
      Bool_t vbfmvacuts = lNJets>=2 && lNJetInGap==0 && vbfMvaValue>0.5;
      Bool_t vhcuts = lNJets==2 && lMJJ>70 && lMJJ<120 && lDiJetPt>150;
      Bool_t boostcuts = lNJets>=1 && lJPt1>30;
      Bool_t passpzeta = 0.85*lPZetaVis - lPZetaMiss <= 25;

      // inclusive
      if (passpzeta) {
	if(isam!=ifake) {
	  hMass_iv[isam]	->Fill(lMSV,   wgt);
	  hMassL_iv[isam]	->Fill(lMSV,   wgt);
          hMassVis_iv[isam]   ->Fill(lMVis,      wgt);
	  nSel_iv[isam]     += wgt;
	  nSelVar_iv[isam]  += wgt*wgt;
	}
      }

      // 0-jet
      if(passpzeta && !bjets && !boostcuts && !vbfcuts) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_0jetv[isam]    ->Fill(lMSV, (isam==iztt) ? kCat_0jet*wgt  : wgt);      
	  }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_0jetv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_0jet*wgt       : wgt);
	  }
          hMassL_0jetv[isam]   ->Fill(lMSV, (isam==iztt) ? kCat_0jet*wgt  : wgt);
          nSel_0jetv[isam]     +=                (isam==iztt) ? kCat_0jet*wgt    : wgt;
          nSelVar_0jetv[isam]  +=                (isam==iztt) ? kCat_0jet*wgt*kCat_0jet*wgt : wgt*wgt;
        }
      }

      // 0-jet, low pT
      if(passpzeta && !bjets && !boostcuts && !vbfcuts && lPt2 <=35) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_0jet_lowptv[isam]    ->Fill(lMSV, (isam==iewk_7TeV) ? ewk_0jet_lowpt_norm*wgt  : wgt);
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_0jet_lowptv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_0jet*wgt       : wgt);
	  }
          hMassL_0jet_lowptv[isam]   ->Fill(lMSV, (isam==iewk_7TeV) ? ewk_0jet_lowpt_norm*wgt  : wgt);
          nSel_0jet_lowptv[isam]     +=                (isam==iewk_7TeV) ? ewk_0jet_lowpt_norm*wgt    : wgt;
          nSelVar_0jet_lowptv[isam]  +=                (isam==iewk_7TeV) ? ewk_0jet_lowpt_norm*wgt*ewk_0jet_lowpt_norm*wgt : wgt*wgt;
        }
      }

      // 0-jet, high pT
      if(passpzeta && !bjets && !boostcuts && !vbfcuts && lPt2 > 35) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_0jet_highptv[isam]    ->Fill(lMSV, (isam==iewk_7TeV) ? ewk_0jet_highpt_norm*wgt  : wgt);
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_0jet_highptv[isam] ->Fill(lMVis,    (isam==iewk_7TeV) ? ewk_0jet_highpt_norm*wgt       : wgt);
	  }
          hMassL_0jet_highptv[isam]   ->Fill(lMSV, (isam==iewk_7TeV) ? ewk_0jet_highpt_norm*wgt  : wgt);
          nSel_0jet_highptv[isam]     +=                (isam==iewk_7TeV) ? ewk_0jet_highpt_norm*wgt    : wgt;
          nSelVar_0jet_highptv[isam]  +=                (isam==iewk_7TeV) ? ewk_0jet_highpt_norm*wgt*ewk_0jet_highpt_norm*wgt : wgt*wgt;
        }
      }

      // 1-jet, no btag
      if(passpzeta && !bjets && boostcuts && !vbfcuts) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_boostv[isam]    ->Fill(lMSV, (isam==iztt) ? kCat_boost*wgt  : wgt);      
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_boostv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_boost*wgt       : wgt);
	  }
          hMassL_boostv[isam]   ->Fill(lMSV, (isam==iztt) ? kCat_boost*wgt  : wgt);
          nSel_boostv[isam]     +=                (isam==iztt) ? kCat_boost*wgt    : wgt;
          nSelVar_boostv[isam]  +=                (isam==iztt) ? kCat_boost*wgt*kCat_boost*wgt : wgt*wgt;
        }
      }

      // 1-jet, no btag, low pT
      if(passpzeta && !bjets && boostcuts && !vbfcuts && lPt2 <=35) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_boost_lowptv[isam]    ->Fill(lMSV, (isam==iztt) ? kCat_boost*wgt  : wgt);
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_boost_lowptv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_boost*wgt       : wgt);
	  }
          hMassL_boost_lowptv[isam]   ->Fill(lMSV, (isam==iztt) ? kCat_boost*wgt  : wgt);
          nSel_boost_lowptv[isam]     +=                (isam==iztt) ? kCat_boost*wgt    : wgt;
          nSelVar_boost_lowptv[isam]  +=                (isam==iztt) ? kCat_boost*wgt*kCat_boost*wgt : wgt*wgt;
        }
      }

      // 1-jet, no btag, high pT
      if(passpzeta && !bjets && boostcuts && !vbfcuts && lPt2 > 35) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_boost_highptv[isam]    ->Fill(lMSV, (isam==iewk_7TeV) ? ewk_boost_highpt_norm*wgt  : wgt);
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_boost_highptv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_boost*wgt       : wgt);
	  }
	  hMassL_boost_highptv[isam]   ->Fill(lMSV, (isam==iewk_7TeV) ? ewk_boost_highpt_norm*wgt  : wgt);
          nSel_boost_highptv[isam]     +=                (isam==iewk_7TeV) ? ewk_boost_highpt_norm*wgt    : wgt;
          nSelVar_boost_highptv[isam]  +=                (isam==iewk_7TeV) ? ewk_boost_highpt_norm*wgt*ewk_boost_highpt_norm*wgt : wgt*wgt;
        }
      }

      // btag
      if(passpzeta && bjets && lNJets<2) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_bv[isam]    ->Fill(lMSV, (isam==iztt) ? kCat_b*wgt  : wgt);      
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_bv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_b*wgt       : wgt);
	  }
          hMassL_bv[isam]   ->Fill(lMSV, (isam==iztt) ? kCat_b*wgt  : wgt);
          nSel_bv[isam]     +=                (isam==iztt) ? kCat_b*wgt    : wgt;
          nSelVar_bv[isam]  +=                (isam==iztt) ? kCat_b*wgt*kCat_b*wgt : wgt*wgt;
        }
      }

      // btag, low pT
      if(passpzeta && bjets && lNJets<2 && lPt2 <= 35) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_b_lowptv[isam]    ->Fill(lMSV, (isam==iztt) ? kCat_b*wgt  : wgt);
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_b_lowptv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_b*wgt       : wgt);
          }
          hMassL_b_lowptv[isam]   ->Fill(lMSV, (isam==iztt) ? kCat_b*wgt  : wgt);
          nSel_b_lowptv[isam]     +=                (isam==iztt) ? kCat_b*wgt    : wgt;
          nSelVar_b_lowptv[isam]  +=                (isam==iztt) ? kCat_b*wgt*kCat_b*wgt : wgt*wgt;
        }
      }

      // btag, high pT
      if(passpzeta && bjets && lNJets<2 && lPt2 > 35) {
        if(isam!=ifake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_b_highptv[isam]    ->Fill(lMSV, (isam==iztt) ? kCat_b*wgt  : wgt);
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_b_highptv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_b*wgt       : wgt);
          }
          hMassL_b_highptv[isam]   ->Fill(lMSV, (isam==iztt) ? kCat_b*wgt  : wgt);
          nSel_b_highptv[isam]     +=                (isam==iztt) ? kCat_b*wgt    : wgt;
          nSelVar_b_highptv[isam]  +=                (isam==iztt) ? kCat_b*wgt*kCat_b*wgt : wgt*wgt;
        }
      }

      // VH
      /*if(passpzeta && !bjets && vhcuts && !vbfcuts) {
        if(isam!=issfake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_2jetv[isam]    ->Fill(lMSV, (isam==iztt) ? kCat_2jet*wgt  : wgt);      
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_2jetv[isam] ->Fill(lMVis,    (isam==iztt) ? kCat_2jet*wgt       : wgt);
	  }
          hMassL_2jetv[isam]   ->Fill(lMVis, (isam==iztt) ? kCat_2jet*wgt  : wgt);
          nSel_2jetv[isam]     +=                (isam==iztt) ? kCat_2jet*wgt    : wgt;
          nSelVar_2jetv[isam]  +=                (isam==iztt) ? kCat_2jet*wgt*kCat_2jet*wgt : wgt*wgt;
        }
      }*/

      // VBF
      if(passpzeta && vbfcuts) {
        if(isam!=issfake) {
          if(isam > 0 || (isam==0 && lMSV < 100.0)) {
            hMass_vbfv[isam]    ->Fill(lMSV, (isam==ittbar_7TeV) ? ttbar_vbf_norm*wgt  : wgt);      
          }
          if(isam > 0 || (isam==0 && lMVis < 60.0)) {
            hMassVis_vbfv[isam] ->Fill(lMVis,    (isam==ittbar_7TeV) ? ttbar_vbf_norm*wgt       : wgt);
	  }
          hMassL_vbfv[isam]   ->Fill(lMSV, (isam==ittbar_7TeV) ? ttbar_vbf_norm*wgt  : wgt);
          nSel_vbfv[isam]     +=                (isam==ittbar_7TeV) ? ttbar_vbf_norm*wgt    : wgt;
          nSelVar_vbfv[isam]  +=                (isam==ittbar_7TeV) ? ttbar_vbf_norm*wgt*ttbar_vbf_norm*wgt : wgt*wgt;
        }
      }
    
    }

    delete infile;
    infile=0, eventTree=0;

  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  TCanvas *c = MakeCanvas("c","c",700,700);
  TCanvas *cr = MakeCanvas("cr","cr",700,700);
  
  // string buffers
  char ylabel[100];   // y-axis label

  //----------------------------------------------------------------------------------------
  // write hists to input file for limit inputs and pre-fit plots
  //----------------------------------------------------------------------------------------

  TString flimitstr(CPlot::sOutDir + "/limit-inputs.root");

  TFile *flimit;
  if (ecorr == kCenter) {
    flimit = new TFile(flimitstr,"recreate");
    flimit->mkdir("emu_0jet");
    flimit->mkdir("emu_0jet_low");
    flimit->mkdir("emu_0jet_high");
    flimit->mkdir("emu_boost");
    flimit->mkdir("emu_boost_low");
    flimit->mkdir("emu_boost_high");
    flimit->mkdir("emu_btag");
    flimit->mkdir("emu_btag_low");
    flimit->mkdir("emu_btag_high");
    flimit->mkdir("emu_vbf");
  }
  else flimit = new TFile(flimitstr,"update");

  // add zmm into the fakes histogram
  if(hMassL_iv[iewk] && hMassL_iv[izmm]) {
    hMassL_0jetv[iewk_7TeV]   ->Add(hMassL_0jetv[izmm]);
    hMassL_0jet_lowptv[iewk_7TeV]   ->Add(hMassL_0jet_lowptv[izmm]);
    hMassL_0jet_highptv[iewk_7TeV]  ->Add(hMassL_0jet_highptv[izmm]);
    hMassL_boostv[iewk]  ->Add(hMassL_boostv[izmm]);
    hMassL_boost_lowptv[iewk]  ->Add(hMassL_boost_lowptv[izmm]);
    hMassL_boost_highptv[iewk_7TeV] ->Add(hMassL_boost_highptv[izmm]);
    hMassL_bv[iewk]      ->Add(hMassL_bv[izmm]);
    hMassL_b_lowptv[iewk]      ->Add(hMassL_b_lowptv[izmm]);
    hMassL_b_highptv[iewk]     ->Add(hMassL_b_highptv[izmm]);
    hMassL_vbfv[iewk]		  ->Add(hMassL_vbfv[izmm]);
  }

  double fakebin = hMassL_vbfv[ifake]->GetBinContent(5);
  double fakebinerr = hMassL_vbfv[ifake]->GetBinError(5);
  hMassL_vbfv[ifake]->SetBinContent(5,0.3*fakebin);
  hMassL_vbfv[ifake]->SetBinError(5,0.3*fakebinerr);
  hMassL_vbfv[ifake]->SetBinContent(3,hMassL_vbfv[ifake]->GetBinContent(3)+0.15*fakebin);
  hMassL_vbfv[ifake]->SetBinError(3,hMassL_vbfv[ifake]->GetBinError(3)+0.15*fakebinerr);
  hMassL_vbfv[ifake]->SetBinContent(4,hMassL_vbfv[ifake]->GetBinContent(4)+0.2*fakebin);
  hMassL_vbfv[ifake]->SetBinError(4,hMassL_vbfv[ifake]->GetBinError(4)+0.2*fakebinerr);
  hMassL_vbfv[ifake]->SetBinContent(6,hMassL_vbfv[ifake]->GetBinContent(6)+0.25*fakebin);
  hMassL_vbfv[ifake]->SetBinError(6,hMassL_vbfv[ifake]->GetBinError(6)+0.25*fakebinerr);
  hMassL_vbfv[ifake]->SetBinContent(7,hMassL_vbfv[ifake]->GetBinContent(7)+0.1*fakebin);
  hMassL_vbfv[ifake]->SetBinError(7,hMassL_vbfv[ifake]->GetBinError(7)+0.1*fakebinerr);
  fakebin = hMass_vbfv[ifake]->GetBinContent(5);
  fakebinerr = hMass_vbfv[ifake]->GetBinError(5);
  hMass_vbfv[ifake]->SetBinContent(5,0.3*fakebin);
  hMass_vbfv[ifake]->SetBinError(5,0.3*fakebinerr);
  hMass_vbfv[ifake]->SetBinContent(3,hMass_vbfv[ifake]->GetBinContent(3)+0.15*fakebin);
  hMass_vbfv[ifake]->SetBinError(3,hMass_vbfv[ifake]->GetBinError(3)+0.15*fakebinerr);
  hMass_vbfv[ifake]->SetBinContent(4,hMass_vbfv[ifake]->GetBinContent(4)+0.2*fakebin);
  hMass_vbfv[ifake]->SetBinError(4,hMass_vbfv[ifake]->GetBinError(4)+0.2*fakebinerr);
  hMass_vbfv[ifake]->SetBinContent(6,hMass_vbfv[ifake]->GetBinContent(6)+0.25*fakebin);
  hMass_vbfv[ifake]->SetBinError(6,hMass_vbfv[ifake]->GetBinError(6)+0.25*fakebinerr);
  hMass_vbfv[ifake]->SetBinContent(7,hMass_vbfv[ifake]->GetBinContent(7)+0.1*fakebin);
  hMass_vbfv[ifake]->SetBinError(7,hMass_vbfv[ifake]->GetBinError(7)+0.1*fakebinerr);

  TString histname;
  for(UInt_t isam=0;isam<samplev.size();isam++) {

    if(isam==izmm) continue;
    if(!domssm && snamev[isam].Contains("mssm")) continue;
    else if(domssm && snamev[isam].Contains("_sm")) continue;

    if(snamev[isam].Contains("ewk",   TString::kIgnoreCase))            histname = "EWK";
    else if(snamev[isam].Contains("fakes", TString::kIgnoreCase))       histname = "Fakes";
    else if(snamev[isam].Contains("ttbar", TString::kIgnoreCase))       histname = "ttbar";
    else if(snamev[isam].Contains("Zmm", TString::kIgnoreCase))         histname = "Zmm";
    else if(snamev[isam].Contains("Ztt",   TString::kIgnoreCase))       histname = "Ztt";
    else if(snamev[isam].Contains("emb",   TString::kIgnoreCase))       histname = "Ztt";
    else if((snamev[isam].Contains("data",  TString::kIgnoreCase)) && !(snamev[isam].Contains("emb",   TString::kIgnoreCase)))          {histname = "data_obs";}
    else if(snamev[isam].Contains("htt_")) continue;
    else if(snamev[isam].Contains("gf_sm")) histname = snamev[isam].ReplaceAll("gf_sm_","ggH");
    else if(snamev[isam].Contains("vbf_sm")) histname = snamev[isam].ReplaceAll("vbf_sm_","qqH");
    else if(snamev[isam].Contains("vtth_sm")) histname = snamev[isam].ReplaceAll("vtth_sm_","VH");
    else if(snamev[isam].Contains("gg_mssm")) histname = snamev[isam].ReplaceAll("gg_mssm_","ggH");
    else if(snamev[isam].Contains("bb_mssm")) histname = snamev[isam].ReplaceAll("bb_mssm_","bbH");
    else { cout << "error! name not found" << endl; assert(0); }

    if(hMassL_vbfv[isam]->Integral(0,hMassL_vbfv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_vbfv[isam]->SetBinContent(ibin,0.00001);
      hMassL_vbfv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_bv[isam]->Integral(0,hMassL_bv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_bv[isam]->SetBinContent(ibin,0.00001);
      hMassL_bv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_b_lowptv[isam]->Integral(0,hMassL_b_lowptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_b_lowptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_b_lowptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_b_highptv[isam]->Integral(0,hMassL_b_highptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_b_highptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_b_highptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_boostv[isam]->Integral(0,hMassL_boostv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_boostv[isam]->SetBinContent(ibin,0.00001);
      hMassL_boostv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_boost_lowptv[isam]->Integral(0,hMassL_boost_lowptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_boost_lowptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_boost_lowptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_boost_highptv[isam]->Integral(0,hMassL_boost_highptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_boost_highptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_boost_highptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_0jetv[isam]->Integral(0,hMassL_0jetv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_0jetv[isam]->SetBinContent(ibin,0.00001);
      hMassL_0jetv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_0jet_lowptv[isam]->Integral(0,hMassL_0jet_lowptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_0jet_lowptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_0jet_lowptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_0jet_highptv[isam]->Integral(0,hMassL_0jet_highptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_0jet_highptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_0jet_highptv[isam]->SetBinError(ibin,0.00001);
    }

    if (isam>0 && isam!=iewk && isam!=iewk_7TeV && isam !=ittbar && isam!=ittbar_7TeV && isam!=ifake && isam!=issfake) {

      if (ecorr==kDown) {
        hMassL_0jetv[isam]->Scale(rescale(hMassL_0jetv[isam],"emu_0jet_low",histname,flimit));
        hMassL_0jet_lowptv[isam]->Scale(rescale(hMassL_0jet_lowptv[isam],"emu_0jet_low",histname,flimit));
        hMassL_0jet_highptv[isam]->Scale(rescale(hMassL_0jet_highptv[isam],"emu_0jet_high",histname,flimit));
        hMassL_boostv[isam]->Scale(rescale(hMassL_boostv[isam],"emu_boost_low",histname,flimit));
        hMassL_boost_lowptv[isam]->Scale(rescale(hMassL_boost_lowptv[isam],"emu_boost_low",histname,flimit));
        hMassL_boost_highptv[isam]->Scale(rescale(hMassL_boost_highptv[isam],"emu_boost_high",histname,flimit));
        hMassL_bv[isam]->Scale(rescale(hMassL_bv[isam],"emu_btag_low",histname,flimit));
        hMassL_b_lowptv[isam]->Scale(rescale(hMassL_b_lowptv[isam],"emu_btag_low",histname,flimit));
        hMassL_b_highptv[isam]->Scale(rescale(hMassL_b_highptv[isam],"emu_btag_high",histname,flimit));
        hMassL_vbfv[isam]->Scale(rescale(hMassL_vbfv[isam],"emu_vbf",histname,flimit));

        histname += "_CMS_scale_e_8TeVDown";
      }
      if (ecorr==kUp) {
        hMassL_0jetv[isam]->Scale(rescale(hMassL_0jetv[isam],"emu_0jet_low",histname,flimit));
        hMassL_0jet_lowptv[isam]->Scale(rescale(hMassL_0jet_lowptv[isam],"emu_0jet_low",histname,flimit));
        hMassL_0jet_highptv[isam]->Scale(rescale(hMassL_0jet_highptv[isam],"emu_0jet_high",histname,flimit));
        hMassL_boostv[isam]->Scale(rescale(hMassL_boostv[isam],"emu_boost_low",histname,flimit));
        hMassL_boost_lowptv[isam]->Scale(rescale(hMassL_boost_lowptv[isam],"emu_boost_low",histname,flimit));
        hMassL_boost_highptv[isam]->Scale(rescale(hMassL_boost_highptv[isam],"emu_boost_high",histname,flimit));
        hMassL_bv[isam]->Scale(rescale(hMassL_bv[isam],"emu_btag_low",histname,flimit));
        hMassL_b_lowptv[isam]->Scale(rescale(hMassL_b_lowptv[isam],"emu_btag_low",histname,flimit));
        hMassL_b_highptv[isam]->Scale(rescale(hMassL_b_highptv[isam],"emu_btag_high",histname,flimit));
        hMassL_vbfv[isam]->Scale(rescale(hMassL_vbfv[isam],"emu_vbf",histname,flimit));

        histname += "_CMS_scale_e_8TeVUp";
      }

    } else {if (ecorr!=kCenter) continue; }

    if (isam!=ifake && isam!=ittbar_7TeV) {
      if(isam!=iewk) {
	flimit->cd("emu_0jet_low");
	hMassL_0jet_lowptv[isam]->SetName(histname);
	hMassL_0jet_lowptv[isam]->Write();
	flimit->cd("emu_0jet_high");
	hMassL_0jet_highptv[isam]->SetName(histname);
	hMassL_0jet_highptv[isam]->Write();
      }
      if(isam!=iewk_7TeV) {
	flimit->cd("emu_0jet");
	hMassL_0jetv[isam]->SetName(histname);
	hMassL_0jetv[isam]->Write();
	flimit->cd("emu_boost");
	hMassL_boostv[isam]->SetName(histname);
	hMassL_boostv[isam]->Write();
	flimit->cd("emu_boost_low");
	hMassL_boost_lowptv[isam]->SetName(histname);
	hMassL_boost_lowptv[isam]->Write();
      }
      if(isam!=iewk) {
	flimit->cd("emu_boost_high");
	hMassL_boost_highptv[isam]->SetName(histname);
	hMassL_boost_highptv[isam]->Write();
      }
      if(isam!=iewk_7TeV) {
	flimit->cd("emu_btag");
	hMassL_bv[isam]->SetName(histname);
	hMassL_bv[isam]->Write();
	flimit->cd("emu_btag_low");
	hMassL_b_lowptv[isam]->SetName(histname);
	hMassL_b_lowptv[isam]->Write();
	flimit->cd("emu_btag_high");
	hMassL_b_highptv[isam]->SetName(histname);
	hMassL_b_highptv[isam]->Write();
      }
    }
    if (isam!=issfake && isam!=ittbar && isam!=iewk_7TeV) {
      flimit->cd("emu_vbf");
      hMassL_vbfv[isam]->SetName(histname);
      hMassL_vbfv[isam]->Write();
    }

  }

  flimit->Close();

  //
  // Begin plots:
  //

  sprintf(ylabel,"Events");
  CPlot plotMet("met","","#slash{E}_{T} [GeV]",ylabel);
  if(hMetv[iewk] && hMetv[izmm]) hMetv[iewk]->Add(hMetv[izmm]);
  if(hMetv[ism_vbf] && hMetv[ism_gf] && hMetv[ism_vtth]) {
    hMetv[ism_vbf]->Add(hMetv[ism_gf]);
    hMetv[ism_vbf]->Add(hMetv[ism_vtth]);
  }
  if(hasData) { plotMet.AddHist1D(hMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMetv[isam]->Scale(5.);
      plotMet.AddToStack(hMetv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMet.AddToStack(hMetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMet.SetLegend(0.5,0.65,0.95,0.9);
  plotMet.Draw(c,ecorr==kCenter,format);

  if(iztt == 9999) iztt = iemb; 
  sprintf(ylabel,"Events");
  CPlot plotMetRatio("met_ratio","","#slash{E}_{T} [GeV]",ylabel);
  TH1F *hMet_MC = new TH1F("hMet_MC","",hMetv[0]->GetNbinsX(),hMetv[0]->GetXaxis()->GetXmin(),hMetv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hMetv[0]->GetNbinsX()+1; i++) {
    hMet_MC->SetBinContent(i,(hMetv[iztt]->GetBinContent(i)+hMetv[ittbar]->GetBinContent(i)+hMetv[iewk]->GetBinContent(i)+hMetv[issfake]->GetBinContent(i)));
    hMet_MC->SetBinError(i,(0.043*hMetv[iztt]->GetBinContent(i)+0.08*hMetv[ittbar]->GetBinContent(i)+0.153*hMetv[iewk]->GetBinContent(i)+0.301*hMetv[issfake]->GetBinContent(i)));
  }
  plotMetRatio.AddHist1D(hMetv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMetRatio.AddToStack(hMetv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMetRatio.AddToStack(hMetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMetRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotMetRatio.DrawRatioStack(cr,hMetv[0],hMet_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotMetRaw("metraw","","#slash{E}_{T} [GeV]",ylabel);
  if(hMetRawv[iewk] && hMetRawv[izmm]) hMetRawv[iewk]->Add(hMetRawv[izmm]);
  if(hasData) { plotMetRaw.AddHist1D(hMetRawv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotMetRaw.AddToStack(hMetRawv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMetRaw.SetLegend(0.5,0.65,0.95,0.9);
  plotMetRaw.Draw(c,ecorr==kCenter,format);

  // stack up the raw met hists
  /*TH1F *hTmpMetStack = 0;
  if(hMetRawv[iewk] && hMetRawv[izmm]) hMetRawv[iewk]->Add(hMetRawv[izmm]);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==1) {
      hTmpMetStack = new TH1F(*hMetv[isam]);
      hTmpMetStack->Reset();
    }
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
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
  plotMetRaw.DrawRatio(c,hMetv[ism_vbf],hMetRawv[ism_vbf],ecorr==kCenter,format);*/

  // projection variables
  sprintf(ylabel,"Events");
  CPlot plotProjMet("pzetamiss","","#slash{p}_{#zeta} [GeV]",ylabel);
  if(hProjMetv[iewk] && hProjMetv[izmm]) hProjMetv[iewk]->Add(hProjMetv[izmm]);
  if(hasData) { plotProjMet.AddHist1D(hProjMetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotProjMet.AddToStack(hProjMetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotProjMet.SetLegend(0.2,0.65,0.65,0.9);
  plotProjMet.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotProjVis("pzetavis","","p_{#zeta}^{vis} [GeV]",ylabel);
  if(hProjVisv[iewk] && hProjVisv[izmm]) hProjVisv[iewk]->Add(hProjVisv[izmm]);
  if(hasData) { plotProjVis.AddHist1D(hProjVisv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotProjVis.AddToStack(hProjVisv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotProjVis.SetLegend(0.5,0.65,0.95,0.9);
  plotProjVis.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotProjVar("pzetavar","","#slash{p}_{#zeta} - 0.85 #times p_{#zeta}^{vis} [GeV]",ylabel);
  if(hProjVarv[iewk] && hProjVarv[izmm]) hProjVarv[iewk]->Add(hProjVarv[izmm]);
  if(hProjVarv[ism_vbf] && hProjVarv[ism_gf] && hProjVarv[ism_vtth]) {
    hProjVarv[ism_vbf]->Add(hProjVarv[ism_gf]);
    hProjVarv[ism_vbf]->Add(hProjVarv[ism_vtth]);
  }
  if(hasData) { plotProjVar.AddHist1D(hProjVarv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hProjVarv[isam]->Scale(5.);
      plotProjVar.AddToStack(hProjVarv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotProjVar.AddToStack(hProjVarv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotProjVar.SetLegend(0.2,0.6,0.65,0.9);
  plotProjVar.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotProjVarRatio("pzetavar_ratio","","#slash{p}_{#zeta} - 0.85 #times p_{#zeta}^{vis} [GeV]",ylabel);
  TH1F *hProjVar_MC = new TH1F("hProjVar_MC","",hProjVarv[0]->GetNbinsX(),hProjVarv[0]->GetXaxis()->GetXmin(),hProjVarv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hProjVarv[0]->GetNbinsX()+1; i++) {
    hProjVar_MC->SetBinContent(i,(hProjVarv[iztt]->GetBinContent(i)+hProjVarv[ittbar]->GetBinContent(i)+hProjVarv[iewk]->GetBinContent(i)+hProjVarv[issfake]->GetBinContent(i)));
    hProjVar_MC->SetBinError(i,(0.043*hProjVarv[iztt]->GetBinContent(i)+0.08*hProjVarv[ittbar]->GetBinContent(i)+0.153*hProjVarv[iewk]->GetBinContent(i)+0.301*hProjVarv[issfake]->GetBinContent(i)));
  }
  plotProjVarRatio.AddHist1D(hProjVarv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotProjVarRatio.AddToStack(hProjVarv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotProjVarRatio.AddToStack(hProjVarv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotProjVarRatio.SetLegend(0.2,0.6,0.65,0.9);
  plotProjVarRatio.DrawRatioStack(cr,hProjVarv[0],hProjVar_MC,ecorr==kCenter,format);
    
  /*sprintf(ylabel,"Events");
  CPlot plotRawProjVar("pzetavar_pf","","#slash{p}_{#zeta} - 0.85 #times p_{#zeta}^{vis} [GeV]",ylabel);
  if(hRawProjVarv[iewk] && hRawProjVarv[izmm]) hRawProjVarv[iewk]->Add(hRawProjVarv[izmm]);
  if(hRawProjVarv[ism_vbf] && hRawProjVarv[ism_gf] && hRawProjVarv[ism_vtth]) {
    hRawProjVarv[ism_vbf]->Add(hRawProjVarv[ism_gf]);
    hRawProjVarv[ism_vbf]->Add(hRawProjVarv[ism_vtth]);
  }
  if(hasData) { plotRawProjVar.AddHist1D(hRawProjVarv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hRawProjVarv[isam]->Scale(5.);
      plotRawProjVar.AddToStack(hRawProjVarv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotRawProjVar.AddToStack(hRawProjVarv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotRawProjVar.SetLegend(0.2,0.6,0.65,0.9);
  plotRawProjVar.Draw(c,ecorr==kCenter,format);*/

  // stack up the raw projection variable hists
  /*TH1F *hTmpStack = 0;
  if(hRawProjVarv[iewk] && hRawProjVarv[izmm]) hRawProjVarv[iewk]->Add(hRawProjVarv[izmm]);
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(isam==1) {
      hTmpStack = new TH1F(*hProjVarv[isam]);
      hTmpStack->Reset();
    }
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
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
  plotRawProjVar.DrawRatio(c,hProjVarv[ism_vbf],hRawProjVarv[ism_vbf],ecorr==kCenter,format);*/

  sprintf(ylabel,"Events");
  CPlot plotNjets("njets","","Number of Jets",ylabel);
  if(hNjetsv[iewk] && hNjetsv[izmm]) hNjetsv[iewk]->Add(hNjetsv[izmm]);
  if(hNjetsv[ism_vbf] && hNjetsv[ism_gf] && hNjetsv[ism_vtth]) {
    hNjetsv[ism_vbf]->Add(hNjetsv[ism_gf]);
    hNjetsv[ism_vbf]->Add(hNjetsv[ism_vtth]);
  }
  if(hasData) { plotNjets.AddHist1D(hNjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hNjetsv[isam]->Scale(5.);
      plotNjets.AddToStack(hNjetsv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotNjets.AddToStack(hNjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNjets.SetYRange(0,1.5*(plotNjets.GetStack()->GetMaximum()));
  plotNjets.SetLegend(0.5,0.6,0.95,0.9);
  plotNjets.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotNjetsRatio("njets_ratio","","Number of Jets",ylabel);
  TH1F *hNjets_MC = new TH1F("hNjets_MC","",hNjetsv[0]->GetNbinsX(),hNjetsv[0]->GetXaxis()->GetXmin(),hNjetsv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hNjetsv[0]->GetNbinsX()+1; i++) {
    hNjets_MC->SetBinContent(i,(hNjetsv[iztt]->GetBinContent(i)+hNjetsv[ittbar]->GetBinContent(i)+hNjetsv[iewk]->GetBinContent(i)+hNjetsv[issfake]->GetBinContent(i)));
    hNjets_MC->SetBinError(i,(0.043*hNjetsv[iztt]->GetBinContent(i)+0.08*hNjetsv[ittbar]->GetBinContent(i)+0.153*hNjetsv[iewk]->GetBinContent(i)+0.301*hNjetsv[issfake]->GetBinContent(i)));
  }
  plotNjetsRatio.AddHist1D(hNjetsv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotNjetsRatio.AddToStack(hNjetsv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotNjetsRatio.AddToStack(hNjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNjetsRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotNjetsRatio.DrawRatioStack(cr,hNjetsv[0],hNjets_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotNjets_log("njets_log","","Number of Jets",ylabel);
  if(hasData) { plotNjets_log.AddHist1D(hNjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotNjets_log.AddHist1D(hNjetsv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotNjets_log.AddToStack(hNjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNjets_log.SetYRange(0.2,30.0*(plotNjets_log.GetStack()->GetMaximum()));
  plotNjets_log.SetLogy();
  plotNjets_log.SetLegend(0.5,0.6,0.95,0.9);
  plotNjets_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotNjets_log_Ratio("njets_log_ratio","","Number of Jets",ylabel);
  plotNjets_log_Ratio.AddHist1D(hNjetsv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotNjets_log_Ratio.AddHist1D(hNjetsv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotNjets_log_Ratio.AddToStack(hNjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNjets_log_Ratio.SetYRange(0.2,30.0*(plotNjets_log_Ratio.GetStack()->GetMaximum()));
  plotNjets_log_Ratio.SetLogy();
  plotNjets_log_Ratio.SetLegend(0.5,0.6,0.95,0.9);
  plotNjets_log_Ratio.DrawRatioStack(cr,hNjetsv[0],hNjets_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotNbjets("nbjets","","Number of b-Tagged Jets",ylabel);
  if(hNbjetsv[iewk] && hNbjetsv[izmm]) hNbjetsv[iewk]->Add(hNbjetsv[izmm]);
  if(hNbjetsv[ism_vbf] && hNbjetsv[ism_gf] && hNbjetsv[ism_vtth]) {
    hNbjetsv[ism_vbf]->Add(hNbjetsv[ism_gf]);
    hNbjetsv[ism_vbf]->Add(hNbjetsv[ism_vtth]);
  }
  if(hasData) { plotNbjets.AddHist1D(hNbjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hNbjetsv[isam]->Scale(5.);
      plotNbjets.AddToStack(hNbjetsv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotNbjets.AddToStack(hNbjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNbjets.SetLegend(0.5,0.6,0.95,0.9);
  plotNbjets.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotNbjetsRatio("nbjets_ratio","","Number of b-Tagged Jets",ylabel);
  TH1F *hNbjets_MC = new TH1F("hNbjets_MC","",hNbjetsv[0]->GetNbinsX(),hNbjetsv[0]->GetXaxis()->GetXmin(),hNbjetsv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hNbjetsv[0]->GetNbinsX()+1; i++) {
    hNbjets_MC->SetBinContent(i,(hNbjetsv[iztt]->GetBinContent(i)+hNbjetsv[ittbar]->GetBinContent(i)+hNbjetsv[iewk]->GetBinContent(i)+hNbjetsv[issfake]->GetBinContent(i)));
    hNbjets_MC->SetBinError(i,(0.043*hNbjetsv[iztt]->GetBinContent(i)+0.08*hNbjetsv[ittbar]->GetBinContent(i)+0.153*hNbjetsv[iewk]->GetBinContent(i)+0.301*hNbjetsv[issfake]->GetBinContent(i)));
  }
  plotNbjetsRatio.AddHist1D(hNbjetsv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotNbjetsRatio.AddToStack(hNbjetsv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotNbjetsRatio.AddToStack(hNbjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNbjetsRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotNbjetsRatio.DrawRatioStack(cr,hNbjetsv[0],hNbjets_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotNbjets_log("nbjets_log","","Number of b-Tagged Jets",ylabel);
  if(hasData) { plotNbjets_log.AddHist1D(hNbjetsv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotNbjets_log.AddHist1D(hNbjetsv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotNbjets_log.AddToStack(hNbjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNbjets_log.SetYRange(0.2,100.0*(plotNbjets_log.GetStack()->GetMaximum()));
  plotNbjets_log.SetLogy();
  plotNbjets_log.SetLegend(0.5,0.6,0.95,0.9);
  plotNbjets_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotNbjets_log_Ratio("nbjets_log_ratio","","Number of b-Tagged Jets",ylabel);
  plotNbjets_log_Ratio.AddHist1D(hNbjetsv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotNbjets_log_Ratio.AddHist1D(hNbjetsv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotNbjets_log_Ratio.AddToStack(hNbjetsv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNbjets_log_Ratio.SetYRange(0.2,100.0*(plotNbjets_log_Ratio.GetStack()->GetMaximum()));
  plotNbjets_log_Ratio.SetLogy();
  plotNbjets_log_Ratio.SetLegend(0.5,0.6,0.95,0.9);
  plotNbjets_log_Ratio.DrawRatioStack(cr,hNbjetsv[0],hNbjets_MC,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBdiscr("btag","","b-Tag Discriminator",ylabel);
  if(hBdiscrv[iewk] && hBdiscrv[izmm]) hBdiscrv[iewk]->Add(hBdiscrv[izmm]);
  if(hBdiscrv[ism_vbf] && hBdiscrv[ism_gf] && hBdiscrv[ism_vtth]) {
    hBdiscrv[ism_vbf]->Add(hBdiscrv[ism_gf]);
    hBdiscrv[ism_vbf]->Add(hBdiscrv[ism_vtth]);
  }
  if(hasData) { plotBdiscr.AddHist1D(hBdiscrv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hBdiscrv[isam]->Scale(5.);
      plotBdiscr.AddToStack(hBdiscrv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBdiscr.AddToStack(hBdiscrv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBdiscr.SetLegend(0.5,0.6,0.95,0.9);
  plotBdiscr.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBdiscrRatio("btag_ratio","","b-Tag Discriminator",ylabel);
  TH1F *hBdiscr_MC = new TH1F("hBdiscr_MC","",hBdiscrv[0]->GetNbinsX(),hBdiscrv[0]->GetXaxis()->GetXmin(),hBdiscrv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hBdiscrv[0]->GetNbinsX()+1; i++) {
    hBdiscr_MC->SetBinContent(i,(hBdiscrv[iztt]->GetBinContent(i)+hBdiscrv[ittbar]->GetBinContent(i)+hBdiscrv[iewk]->GetBinContent(i)+hBdiscrv[issfake]->GetBinContent(i)));
    hBdiscr_MC->SetBinError(i,(0.043*hBdiscrv[iztt]->GetBinContent(i)+0.08*hBdiscrv[ittbar]->GetBinContent(i)+0.153*hBdiscrv[iewk]->GetBinContent(i)+0.301*hBdiscrv[issfake]->GetBinContent(i)));
  }
  plotBdiscrRatio.AddHist1D(hBdiscrv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotBdiscrRatio.AddToStack(hBdiscrv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBdiscrRatio.AddToStack(hBdiscrv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBdiscrRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotBdiscrRatio.DrawRatioStack(cr,hBdiscrv[0],hBdiscr_MC,ecorr==kCenter,format);
    
  /*sprintf(ylabel,"Events");
  CPlot plotBdiscr_vbf("btag_vbf","","b-tag discr.",ylabel);
  if(hBdiscr_vbfv[iewk] && hBdiscr_vbfv[izmm]) hBdiscr_vbfv[iewk]->Add(hBdiscr_vbfv[izmm]);
  if(hasData) { plotBdiscr_vbf.AddHist1D(hBdiscr_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    plotBdiscr_vbf.AddToStack(hBdiscr_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBdiscr_vbf.SetLegend(0.5,0.65,0.95,0.9);
  plotBdiscr_vbf.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotDPhi("dphi_emu","","#Delta^{}#phi_{e#mu} [deg]",ylabel);
  if(hDPhiv[iewk] && hDPhiv[izmm]) hDPhiv[iewk]->Add(hDPhiv[izmm]);
  if(hDPhiv[ism_vbf] && hDPhiv[ism_gf] && hDPhiv[ism_vtth]) {
    hDPhiv[ism_vbf]->Add(hDPhiv[ism_gf]);
    hDPhiv[ism_vbf]->Add(hDPhiv[ism_vtth]);
  }
  if(hasData) { plotDPhi.AddHist1D(hDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotDPhi.AddToStack(hDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotDPhi.TransLegend(-0.15,0);
  assert(plotDPhi.GetStack());
  plotDPhi.SetLegend(0.5,0.65,0.95,0.9);
  plotDPhi.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotMt("mt","","m_{T}(ll,#slash{E}_{T}) [GeV]",ylabel);
  if(hMtv[iewk] && hMtv[izmm]) hMtv[iewk]->Add(hMtv[izmm]);
  if(hasData) { plotMt.AddHist1D(hMtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotMt.AddToStack(hMtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMt.SetLegend(0.5,0.65,0.95,0.9);
  plotMt.Draw(c,ecorr==kCenter,format);*/

  /*sprintf(ylabel,"Events");
  CPlot plotPt("pt","","p_{T}^{ll} [GeV]",ylabel);
  if(hPtv[iewk] && hPtv[izmm]) hPtv[iewk]->Add(hPtv[izmm]);
  if(hasData) { plotPt.AddHist1D(hPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotPt.AddToStack(hPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPt.SetLegend(0.5,0.65,0.95,0.9);
  plotPt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotLepDEta("lepdeta","","#Delta^{}#eta(ll)",ylabel);
  if(hLepDEtav[iewk] && hLepDEtav[izmm]) hLepDEtav[iewk]->Add(hLepDEtav[izmm]);
  if(hLepDEtav[ism_vbf] && hLepDEtav[ism_gf] && hLepDEtav[ism_vtth]) {
    hLepDEtav[ism_vbf]->Add(hLepDEtav[ism_gf]);
    hLepDEtav[ism_vbf]->Add(hLepDEtav[ism_vtth]);
  }
  if(hasData) { plotLepDEta.AddHist1D(hLepDEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hLepDEtav[isam]->Scale(5.);
      plotLepDEta.AddToStack(hLepDEtav[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotLepDEta.AddToStack(hLepDEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotLepDEta.SetLegend(0.5,0.6,0.95,0.9);
  plotLepDEta.Draw(c,ecorr==kCenter,format);*/

  /*sprintf(ylabel,"Events");
  CPlot plotMetDPhi("metdphi","","#Delta^{}#phi(ll,#slash{E}_{T}) [deg]",ylabel);
  if(hMetDPhiv[iewk] && hMetDPhiv[izmm]) hMetDPhiv[iewk]->Add(hMetDPhiv[izmm]);
  if(hasData) { plotMetDPhi.AddHist1D(hMetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotMetDPhi.AddToStack(hMetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMetDPhi.SetLegend(0.5,0.65,0.95,0.9);
  plotMetDPhi.Draw(c,ecorr==kCenter,format);*/
    
  /*sprintf(ylabel,"Events");
  CPlot plotPt1("pt1","","leading lepton p_{T} [GeV]",ylabel);
  if(hPt1v[iewk] && hPt1v[izmm]) hPt1v[iewk]->Add(hPt1v[izmm]);
  if(hasData) { plotPt1.AddHist1D(hPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotPt1.AddToStack(hPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPt1.SetLegend(0.5,0.65,0.95,0.9);
  plotPt1.Draw(c,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEta1("eta1","","leading lepton #eta",ylabel);
  if(hEta1v[iewk] && hEta1v[izmm]) hEta1v[iewk]->Add(hEta1v[izmm]);
  if(hasData) { plotEta1.AddHist1D(hEta1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotEta1.AddToStack(hEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEta1.SetYRange(0,2.0*(plotEta1.GetStack()->GetMaximum()));
  plotEta1.SetLegend(0.5,0.65,0.95,0.9);
  plotEta1.Draw(c,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPhi1("phi1","","leading lepton #phi",ylabel);
  if(hPhi1v[iewk] && hPhi1v[izmm]) hPhi1v[iewk]->Add(hPhi1v[izmm]);
  if(hasData) { plotPhi1.AddHist1D(hPhi1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotPhi1.AddToStack(hPhi1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPhi1.SetYRange(0,2.6*(plotPhi1.GetStack()->GetMaximum()));
  plotPhi1.SetLegend(0.5,0.65,0.95,0.9);
  plotPhi1.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotPt2("pt2","","trailing lepton p_{T} [GeV]",ylabel);
  if(hPt2v[iewk] && hPt2v[izmm]) hPt2v[iewk]->Add(hPt2v[izmm]);
  if(hasData) { plotPt2.AddHist1D(hPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotPt2.AddToStack(hPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPt2.SetLegend(0.5,0.65,0.95,0.9);
  plotPt2.Draw(c,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEta2("eta2","","trailing lepton #eta",ylabel);
  if(hEta2v[iewk] && hEta2v[izmm]) hEta2v[iewk]->Add(hEta2v[izmm]);
  if(hasData) { plotEta2.AddHist1D(hEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotEta2.AddToStack(hEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEta2.SetYRange(0,2.0*(plotEta2.GetStack()->GetMaximum()));
  plotEta2.SetLegend(0.5,0.65,0.95,0.9);
  plotEta2.Draw(c,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPhi2("phi2","","trailing lepton #phi",ylabel);
  if(hPhi2v[iewk] && hPhi2v[izmm]) hPhi2v[iewk]->Add(hPhi2v[izmm]);
  if(hasData) { plotPhi2.AddHist1D(hPhi2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotPhi2.AddToStack(hPhi2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPhi2.SetYRange(0,2.6*(plotPhi2.GetStack()->GetMaximum()));
  plotPhi2.SetLegend(0.5,0.65,0.95,0.9);
  plotPhi2.Draw(c,ecorr==kCenter,format);*/

  // mu / electron kinematics
  sprintf(ylabel,"Events");
  CPlot plotPtMu("pt_mu","","#mu p_{T} [GeV]",ylabel);
  if(hPtMuv[iewk] && hPtMuv[izmm]) hPtMuv[iewk]->Add(hPtMuv[izmm]);
  if(hPtMuv[ism_vbf] && hPtMuv[ism_gf] && hPtMuv[ism_vtth]) {
    hPtMuv[ism_vbf]->Add(hPtMuv[ism_gf]);
    hPtMuv[ism_vbf]->Add(hPtMuv[ism_vtth]);
  }
  if(hasData) { plotPtMu.AddHist1D(hPtMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPtMuv[isam]->Scale(5.);
      plotPtMu.AddToStack(hPtMuv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtMu.AddToStack(hPtMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtMu.SetLegend(0.5,0.6,0.95,0.9);
  plotPtMu.Draw(c,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPtMuRatio("pt_mu_ratio","","#mu p_{T} [GeV]",ylabel);
  TH1F *hPtMu_MC = new TH1F("hPtMu_MC","",hPtMuv[0]->GetNbinsX(),hPtMuv[0]->GetXaxis()->GetXmin(),hPtMuv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hPtMuv[0]->GetNbinsX()+1; i++) {
    hPtMu_MC->SetBinContent(i,(hPtMuv[iztt]->GetBinContent(i)+hPtMuv[ittbar]->GetBinContent(i)+hPtMuv[iewk]->GetBinContent(i)+hPtMuv[issfake]->GetBinContent(i)));
    hPtMu_MC->SetBinError(i,(0.043*hPtMuv[iztt]->GetBinContent(i)+0.08*hPtMuv[ittbar]->GetBinContent(i)+0.153*hPtMuv[iewk]->GetBinContent(i)+0.301*hPtMuv[issfake]->GetBinContent(i)));
  }
  plotPtMuRatio.AddHist1D(hPtMuv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotPtMuRatio.AddToStack(hPtMuv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtMuRatio.AddToStack(hPtMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtMuRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotPtMuRatio.DrawRatioStack(cr,hPtMuv[0],hPtMu_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEtaMu("eta_mu","","#mu #eta",ylabel);
  if(hEtaMuv[iewk] && hEtaMuv[izmm]) hEtaMuv[iewk]->Add(hEtaMuv[izmm]);
  if(hEtaMuv[ism_vbf] && hEtaMuv[ism_gf] && hEtaMuv[ism_vtth]) {
    hEtaMuv[ism_vbf]->Add(hEtaMuv[ism_gf]);
    hEtaMuv[ism_vbf]->Add(hEtaMuv[ism_vtth]);
  }
  if(hasData) { plotEtaMu.AddHist1D(hEtaMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hEtaMuv[isam]->Scale(5.);
      plotEtaMu.AddToStack(hEtaMuv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotEtaMu.AddToStack(hEtaMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaMu.SetYRange(0,2.0*(plotEtaMu.GetStack()->GetMaximum()));
  plotEtaMu.SetLegend(0.5,0.6,0.95,0.9);
  plotEtaMu.Draw(c,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEtaMuRatio("eta_mu_ratio","","#mu #eta",ylabel);
  TH1F *hEtaMu_MC = new TH1F("hEtaMu_MC","",hEtaMuv[0]->GetNbinsX(),hEtaMuv[0]->GetXaxis()->GetXmin(),hEtaMuv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hEtaMuv[0]->GetNbinsX()+1; i++) {
    hEtaMu_MC->SetBinContent(i,(hEtaMuv[iztt]->GetBinContent(i)+hEtaMuv[ittbar]->GetBinContent(i)+hEtaMuv[iewk]->GetBinContent(i)+hEtaMuv[issfake]->GetBinContent(i)));
    hEtaMu_MC->SetBinError(i,(0.043*hEtaMuv[iztt]->GetBinContent(i)+0.08*hEtaMuv[ittbar]->GetBinContent(i)+0.153*hEtaMuv[iewk]->GetBinContent(i)+0.301*hEtaMuv[issfake]->GetBinContent(i)));
  }
  plotEtaMuRatio.AddHist1D(hEtaMuv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotEtaMuRatio.AddToStack(hEtaMuv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotEtaMuRatio.AddToStack(hEtaMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaMuRatio.SetYRange(0,2.0*(plotEtaMuRatio.GetStack()->GetMaximum()));
  plotEtaMuRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotEtaMuRatio.DrawRatioStack(cr,hEtaMuv[0],hEtaMu_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPhiMu("phimu","","#mu #phi",ylabel);
  if(hPhiMuv[iewk] && hPhiMuv[izmm]) hPhiMuv[iewk]->Add(hPhiMuv[izmm]);
  if(hPhiMuv[ism_vbf] && hPhiMuv[ism_gf] && hPhiMuv[ism_vtth]) {
    hPhiMuv[ism_vbf]->Add(hPhiMuv[ism_gf]);
    hPhiMuv[ism_vbf]->Add(hPhiMuv[ism_vtth]);
  }
  if(hasData) { plotPhiMu.AddHist1D(hPhiMuv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPhiMuv[isam]->Scale(5.);
      plotPhiMu.AddToStack(hPhiMuv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPhiMu.AddToStack(hPhiMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPhiMu.SetYRange(0,2.6*(plotPhiMu.GetStack()->GetMaximum()));
  plotPhiMu.SetLegend(0.5,0.6,0.95,0.9);
  plotPhiMu.Draw(c,ecorr==kCenter,format);

  /*sprintf(ylabel,"Events");
  CPlot plotD0Mu("d0mu","","#mu d0",ylabel);
  if(hD0Muv[iewk] && hD0Muv[izmm]) hD0Muv[iewk]->Add(hD0Muv[izmm]);
  if(hD0Muv[ism_vbf] && hD0Muv[ism_gf] && hD0Muv[ism_vtth]) {
    hD0Muv[ism_vbf]->Add(hD0Muv[ism_gf]);
    hD0Muv[ism_vbf]->Add(hD0Muv[ism_vtth]);
  }
  if(hasData) { plotD0Mu.AddHist1D(hD0Muv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hD0Muv[isam]->Scale(5.);
      plotD0Mu.AddHist1D(hD0Muv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
      //plotD0Mu.AddToStack(hD0Muv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotD0Mu.AddToStack(hD0Muv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotD0Mu.SetLogy();
  plotD0Mu.SetLegend(0.5,0.6,0.95,0.9);
  plotD0Mu.Draw(c,ecorr==kCenter,format);*/

  sprintf(ylabel,"Events");
  CPlot plotPtEle("pt_e","","e p_{T} [GeV]",ylabel);
  if(hPtElev[iewk] && hPtElev[izmm]) hPtElev[iewk]->Add(hPtElev[izmm]);
  if(hPtElev[ism_vbf] && hPtElev[ism_gf] && hPtElev[ism_vtth]) {
    hPtElev[ism_vbf]->Add(hPtElev[ism_gf]);
    hPtElev[ism_vbf]->Add(hPtElev[ism_vtth]);
  }
  if(hasData) { plotPtEle.AddHist1D(hPtElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPtElev[isam]->Scale(5.);
      plotPtEle.AddToStack(hPtElev[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtEle.AddToStack(hPtElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtEle.SetLegend(0.5,0.6,0.95,0.9);
  plotPtEle.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotPtEleRatio("pt_e_ratio","","e p_{T} [GeV]",ylabel);
  TH1F *hPtEle_MC = new TH1F("hPtEle_MC","",hPtElev[0]->GetNbinsX(),hPtElev[0]->GetXaxis()->GetXmin(),hPtElev[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hPtElev[0]->GetNbinsX()+1; i++) {
    hPtEle_MC->SetBinContent(i,(hPtElev[iztt]->GetBinContent(i)+hPtElev[ittbar]->GetBinContent(i)+hPtElev[iewk]->GetBinContent(i)+hPtElev[issfake]->GetBinContent(i)));
    hPtEle_MC->SetBinError(i,(0.043*hPtElev[iztt]->GetBinContent(i)+0.08*hPtElev[ittbar]->GetBinContent(i)+0.153*hPtElev[iewk]->GetBinContent(i)+0.301*hPtElev[issfake]->GetBinContent(i)));
  }
  plotPtEleRatio.AddHist1D(hPtElev[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotPtEleRatio.AddToStack(hPtElev[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtEleRatio.AddToStack(hPtElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtEleRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotPtEleRatio.DrawRatioStack(cr,hPtElev[0],hPtEle_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEtaEle("eta_e","","e #eta",ylabel);
  if(hEtaElev[iewk] && hEtaElev[izmm]) hEtaElev[iewk]->Add(hEtaElev[izmm]);
  if(hEtaElev[ism_vbf] && hEtaElev[ism_gf] && hEtaElev[ism_vtth]) {
    hEtaElev[ism_vbf]->Add(hEtaElev[ism_gf]);
    hEtaElev[ism_vbf]->Add(hEtaElev[ism_vtth]);
  }
  if(hasData) { plotEtaEle.AddHist1D(hEtaElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hEtaElev[isam]->Scale(5.);
      plotEtaEle.AddToStack(hEtaElev[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotEtaEle.AddToStack(hEtaElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaEle.SetYRange(0,2.0*(plotEtaEle.GetStack()->GetMaximum()));
  plotEtaEle.SetLegend(0.5,0.6,0.95,0.9);
  plotEtaEle.Draw(c,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotEtaEleRatio("eta_e_ratio","","e #eta",ylabel);
  TH1F *hEtaEle_MC = new TH1F("hEtaEle_MC","",hEtaElev[0]->GetNbinsX(),hEtaElev[0]->GetXaxis()->GetXmin(),hEtaElev[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hEtaElev[0]->GetNbinsX()+1; i++) {
    hEtaEle_MC->SetBinContent(i,(hEtaElev[iztt]->GetBinContent(i)+hEtaElev[ittbar]->GetBinContent(i)+hEtaElev[iewk]->GetBinContent(i)+hEtaElev[issfake]->GetBinContent(i)));
    hEtaEle_MC->SetBinError(i,(0.043*hEtaElev[iztt]->GetBinContent(i)+0.08*hEtaElev[ittbar]->GetBinContent(i)+0.153*hEtaElev[iewk]->GetBinContent(i)+0.301*hEtaElev[issfake]->GetBinContent(i)));
  }
  plotEtaEleRatio.AddHist1D(hEtaElev[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotEtaEleRatio.AddToStack(hEtaElev[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotEtaEleRatio.AddToStack(hEtaElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaEleRatio.SetYRange(0,2.0*(plotEtaEleRatio.GetStack()->GetMaximum()));
  plotEtaEleRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotEtaEleRatio.DrawRatioStack(cr,hEtaElev[0],hEtaEle_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPhiEle("phiele","","e #phi",ylabel);
  if(hPhiElev[iewk] && hPhiElev[izmm]) hPhiElev[iewk]->Add(hPhiElev[izmm]);
  if(hPhiElev[ism_vbf] && hPhiElev[ism_gf] && hPhiElev[ism_vtth]) {
    hPhiElev[ism_vbf]->Add(hPhiElev[ism_gf]);
    hPhiElev[ism_vbf]->Add(hPhiElev[ism_vtth]);
  }
  if(hasData) { plotPhiEle.AddHist1D(hPhiElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPhiElev[isam]->Scale(5.);
      plotPhiEle.AddToStack(hPhiElev[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPhiEle.AddToStack(hPhiElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPhiEle.SetYRange(0,2.6*(plotPhiEle.GetStack()->GetMaximum()));
  plotPhiEle.SetLegend(0.5,0.6,0.95,0.9);
  plotPhiEle.Draw(c,ecorr==kCenter,format);

  /*sprintf(ylabel,"Events");
  CPlot plotD0Ele("d0ele","","e d0",ylabel);
  if(hD0Elev[iewk] && hD0Elev[izmm]) hD0Elev[iewk]->Add(hD0Elev[izmm]);
  if(hD0Elev[ism_vbf] && hD0Elev[ism_gf] && hD0Elev[ism_vtth]) {
    hD0Elev[ism_vbf]->Add(hD0Elev[ism_gf]);
    hD0Elev[ism_vbf]->Add(hD0Elev[ism_vtth]);
  }
  if(hasData) { plotD0Ele.AddHist1D(hD0Elev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hD0Elev[isam]->Scale(5.);
      plotD0Ele.AddHist1D(hD0Elev[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
      //plotD0Ele.AddToStack(hD0Elev[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotD0Ele.AddToStack(hD0Elev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotD0Ele.SetLogy();
  plotD0Ele.SetLegend(0.5,0.6,0.95,0.9);
  plotD0Ele.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotMtEle("mtele","","m_{T}_{e,MET}",ylabel);
  if(hMtElev[iewk] && hMtElev[izmm]) hMtElev[iewk]->Add(hMtElev[izmm]);
  if(hasData) { plotMtEle.AddHist1D(hMtElev[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotMtEle.AddToStack(hMtElev[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMtEle.SetYRange(0,2.6*(plotMtEle.GetStack()->GetMaximum()));
  plotMtEle.SetLegend(0.5,0.65,0.95,0.9);
  plotMtEle.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotMtMu("mtmu","","m_{T}_{#mu,MET}",ylabel);
  if(hMtMuv[iewk] && hMtMuv[izmm]) hMtMuv[iewk]->Add(hMtMuv[izmm]);
  if(hasData) { plotMtMu.AddHist1D(hMtMuv[0],samplev[0]->label,"E"); } 
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotMtMu.AddToStack(hMtMuv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMtMu.SetYRange(0,2.6*(plotMtMu.GetStack()->GetMaximum()));
  plotMtMu.SetLegend(0.5,0.65,0.95,0.9);
  plotMtMu.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotJetPt1("jetpt1","","leading jet pt [GeV]",ylabel);
  if(hJetPt1v[iewk] && hJetPt1v[izmm]) hJetPt1v[iewk]->Add(hJetPt1v[izmm]);
  if(hasData) { plotJetPt1.AddHist1D(hJetPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotJetPt1.AddToStack(hJetPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetPt1.SetLegend(0.5,0.65,0.95,0.9);
  plotJetPt1.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotJetPt2("jetpt2","","second jet pt [GeV]",ylabel);
  if(hJetPt2v[iewk] && hJetPt2v[izmm]) hJetPt2v[iewk]->Add(hJetPt2v[izmm]);
  if(hasData) { plotJetPt2.AddHist1D(hJetPt2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotJetPt2.AddToStack(hJetPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetPt2.SetLegend(0.5,0.65,0.95,0.9);
  plotJetPt2.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotJetEta1("jeteta1","","leading jet #eta",ylabel);
  if(hJetEta1v[iewk] && hJetEta1v[izmm]) hJetEta1v[iewk]->Add(hJetEta1v[izmm]);
  if(hasData) { plotJetEta1.AddHist1D(hJetEta1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotJetEta1.AddToStack(hJetEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetEta1.SetLegend(0.5,0.65,0.95,0.9);
  plotJetEta1.SetYRange(0,2.0*(plotJetEta1.GetStack()->GetMaximum()));
  plotJetEta1.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotJetEta2("jeteta2","","second jet #eta",ylabel);
  if(hJetEta2v[iewk] && hJetEta2v[izmm]) hJetEta2v[iewk]->Add(hJetEta2v[izmm]);
  if(hasData) { plotJetEta2.AddHist1D(hJetEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotJetEta2.AddToStack(hJetEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetEta2.SetLegend(0.5,0.65,0.95,0.9);
  plotJetEta2.SetYRange(0,2.0*(plotJetEta2.GetStack()->GetMaximum()));
  plotJetEta2.Draw(c,ecorr==kCenter,format);*/

  CPlot plotJetDPhi("jetdphi","","#Delta#phi (jj)",ylabel);
  if(hJetDPhiv[iewk] && hJetDPhiv[izmm]) hJetDPhiv[iewk]->Add(hJetDPhiv[izmm]);
  if(hJetDPhiv[ism_vbf] && hJetDPhiv[ism_gf] && hJetDPhiv[ism_vtth]) {
    hJetDPhiv[ism_vbf]->Add(hJetDPhiv[ism_gf]);
    hJetDPhiv[ism_vbf]->Add(hJetDPhiv[ism_vtth]);
  }
  if(hasData) { plotJetDPhi.AddHist1D(hJetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hJetDPhiv[isam]->Scale(5.);
      plotJetDPhi.AddToStack(hJetDPhiv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotJetDPhi.AddToStack(hJetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetDPhi.SetLegend(0.2,0.6,0.65,0.9);
  plotJetDPhi.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotJetDPhiRatio("jetdphi_ratio","","#Delta#phi (jj)",ylabel);
  TH1F *hJetDPhi_MC = new TH1F("hJetDPhi_MC","",hJetDPhiv[0]->GetNbinsX(),hJetDPhiv[0]->GetXaxis()->GetXmin(),hJetDPhiv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hJetDPhiv[0]->GetNbinsX()+1; i++) {
    hJetDPhi_MC->SetBinContent(i,(hJetDPhiv[iztt]->GetBinContent(i)+hJetDPhiv[ittbar]->GetBinContent(i)+hJetDPhiv[iewk]->GetBinContent(i)+hJetDPhiv[issfake]->GetBinContent(i)));
    hJetDPhi_MC->SetBinError(i,(0.205*hJetDPhiv[iztt]->GetBinContent(i)+0.215*hJetDPhiv[ittbar]->GetBinContent(i)+0.252*hJetDPhiv[iewk]->GetBinContent(i)+0.362*hJetDPhiv[issfake]->GetBinContent(i)));
  }
  plotJetDPhiRatio.AddHist1D(hJetDPhiv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotJetDPhiRatio.AddToStack(hJetDPhiv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotJetDPhiRatio.AddToStack(hJetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetDPhiRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotJetDPhiRatio.DrawRatioStack(cr,hJetDPhiv[0],hJetDPhi_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotJetDPhi_nobtag("jetdphi_nobtag","","#Delta#phi (jj)",ylabel);
  if(hJetDPhi_nobtagv[iewk] && hJetDPhi_nobtagv[izmm]) hJetDPhi_nobtagv[iewk]->Add(hJetDPhi_nobtagv[izmm]);
  if(hJetDPhi_nobtagv[ism_vbf] && hJetDPhi_nobtagv[ism_gf] && hJetDPhi_nobtagv[ism_vtth]) {
    hJetDPhi_nobtagv[ism_vbf]->Add(hJetDPhi_nobtagv[ism_gf]);
    hJetDPhi_nobtagv[ism_vbf]->Add(hJetDPhi_nobtagv[ism_vtth]);
  }
  if(hasData) { plotJetDPhi_nobtag.AddHist1D(hJetDPhi_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hJetDPhi_nobtagv[isam]->Scale(5.);
      plotJetDPhi_nobtag.AddToStack(hJetDPhi_nobtagv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotJetDPhi_nobtag.AddToStack(hJetDPhi_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotJetDPhi_nobtag.SetLegend(0.2,0.6,0.65,0.9);
  plotJetDPhi_nobtag.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotMjj("mjj","","M(jj) [GeV]",ylabel);
  if(hMjjv[iewk] && hMjjv[izmm]) hMjjv[iewk]->Add(hMjjv[izmm]);
  if(hMjjv[ism_vbf] && hMjjv[ism_gf] && hMjjv[ism_vtth]) {
    hMjjv[ism_vbf]->Add(hMjjv[ism_gf]);
    hMjjv[ism_vbf]->Add(hMjjv[ism_vtth]);
  }
  if(hasData) { plotMjj.AddHist1D(hMjjv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMjjv[isam]->Scale(5.);
      plotMjj.AddToStack(hMjjv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMjj.AddToStack(hMjjv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMjj.SetYRange(0.0, 1.8*(plotMjj.GetStack()->GetMaximum()));
  plotMjj.SetLegend(0.5,0.6,0.95,0.9);
  plotMjj.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotMjjRatio("mjj_ratio","","M(jj) [GeV]",ylabel);
  TH1F *hMjj_MC = new TH1F("hMjj_MC","",hMjjv[0]->GetNbinsX(),hMjjv[0]->GetXaxis()->GetXmin(),hMjjv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hMjjv[0]->GetNbinsX()+1; i++) {
    hMjj_MC->SetBinContent(i,(hMjjv[iztt]->GetBinContent(i)+hMjjv[ittbar]->GetBinContent(i)+hMjjv[iewk]->GetBinContent(i)+hMjjv[issfake]->GetBinContent(i)));
    hMjj_MC->SetBinError(i,(0.205*hMjjv[iztt]->GetBinContent(i)+0.215*hMjjv[ittbar]->GetBinContent(i)+0.252*hMjjv[iewk]->GetBinContent(i)+0.362*hMjjv[issfake]->GetBinContent(i)));
  }
  plotMjjRatio.AddHist1D(hMjjv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMjjRatio.AddToStack(hMjjv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMjjRatio.AddToStack(hMjjv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMjjRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotMjjRatio.DrawRatioStack(cr,hMjjv[0],hMjj_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotMjj_nobtag("mjj_nobtag","","M(jj) [GeV]",ylabel);
  if(hMjj_nobtagv[iewk] && hMjj_nobtagv[izmm]) hMjj_nobtagv[iewk]->Add(hMjj_nobtagv[izmm]);
  if(hMjj_nobtagv[ism_vbf] && hMjj_nobtagv[ism_gf] && hMjj_nobtagv[ism_vtth]) {
    hMjj_nobtagv[ism_vbf]->Add(hMjj_nobtagv[ism_gf]);
    hMjj_nobtagv[ism_vbf]->Add(hMjj_nobtagv[ism_vtth]);
  }
  if(hasData) { plotMjj_nobtag.AddHist1D(hMjj_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMjj_nobtagv[isam]->Scale(5.);
      plotMjj_nobtag.AddToStack(hMjj_nobtagv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMjj_nobtag.AddToStack(hMjj_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMjj_nobtag.SetLegend(0.5,0.6,0.95,0.9);
  plotMjj_nobtag.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotPtjj("ptjj","","p_{T}(jj) [GeV]",ylabel);
  if(hPtjjv[iewk] && hPtjjv[izmm]) hPtjjv[iewk]->Add(hPtjjv[izmm]);
  if(hPtjjv[ism_vbf] && hPtjjv[ism_gf] && hPtjjv[ism_vtth]) {
    hPtjjv[ism_vbf]->Add(hPtjjv[ism_gf]);
    hPtjjv[ism_vbf]->Add(hPtjjv[ism_vtth]);
  }
  if(hasData) { plotPtjj.AddHist1D(hPtjjv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPtjjv[isam]->Scale(5.);
      plotPtjj.AddToStack(hPtjjv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtjj.AddToStack(hPtjjv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtjj.SetLegend(0.5,0.6,0.95,0.9);
  plotPtjj.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotPtjjRatio("ptjj_ratio","","p_{T}(jj) [GeV]",ylabel);
  TH1F *hPtjj_MC = new TH1F("hPtjj_MC","",hPtjjv[0]->GetNbinsX(),hPtjjv[0]->GetXaxis()->GetXmin(),hPtjjv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hPtjjv[0]->GetNbinsX()+1; i++) {
    hPtjj_MC->SetBinContent(i,(hPtjjv[iztt]->GetBinContent(i)+hPtjjv[ittbar]->GetBinContent(i)+hPtjjv[iewk]->GetBinContent(i)+hPtjjv[issfake]->GetBinContent(i)));
    hPtjj_MC->SetBinError(i,(0.205*hPtjjv[iztt]->GetBinContent(i)+0.215*hPtjjv[ittbar]->GetBinContent(i)+0.252*hPtjjv[iewk]->GetBinContent(i)+0.362*hPtjjv[issfake]->GetBinContent(i)));
  }
  plotPtjjRatio.AddHist1D(hPtjjv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotPtjjRatio.AddToStack(hPtjjv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtjjRatio.AddToStack(hPtjjv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtjjRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotPtjjRatio.DrawRatioStack(cr,hPtjjv[0],hPtjj_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPtjj_nobtag("ptjj_nobtag","","p_{T}(jj) [GeV]",ylabel);
  if(hPtjj_nobtagv[iewk] && hPtjj_nobtagv[izmm]) hPtjj_nobtagv[iewk]->Add(hPtjj_nobtagv[izmm]);
  if(hPtjj_nobtagv[ism_vbf] && hPtjj_nobtagv[ism_gf] && hPtjj_nobtagv[ism_vtth]) {
    hPtjj_nobtagv[ism_vbf]->Add(hPtjj_nobtagv[ism_gf]);
    hPtjj_nobtagv[ism_vbf]->Add(hPtjj_nobtagv[ism_vtth]);
  }
  if(hasData) { plotPtjj_nobtag.AddHist1D(hPtjj_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPtjj_nobtagv[isam]->Scale(5.);
      plotPtjj_nobtag.AddToStack(hPtjj_nobtagv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtjj_nobtag.AddToStack(hPtjj_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtjj_nobtag.SetLegend(0.5,0.6,0.95,0.9);
  plotPtjj_nobtag.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotPtH("ditaupt","","p_{T}(#tau#tau) [GeV]",ylabel);
  if(hPtHv[iewk] && hPtHv[izmm]) hPtHv[iewk]->Add(hPtHv[izmm]);
  if(hPtHv[ism_vbf] && hPtHv[ism_gf] && hPtHv[ism_vtth]) {
    hPtHv[ism_vbf]->Add(hPtHv[ism_gf]);
    hPtHv[ism_vbf]->Add(hPtHv[ism_vtth]);
  }
  if(hasData) { plotPtH.AddHist1D(hPtHv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPtHv[isam]->Scale(5.);
      plotPtH.AddToStack(hPtHv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtH.AddToStack(hPtHv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtH.SetLegend(0.5,0.6,0.95,0.9);
  plotPtH.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotPtHRatio("ditaupt_ratio","","p_{T}(#tau#tau) [GeV]",ylabel);
  TH1F *hPtH_MC = new TH1F("hPtH_MC","",hPtHv[0]->GetNbinsX(),hPtHv[0]->GetXaxis()->GetXmin(),hPtHv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hPtHv[0]->GetNbinsX()+1; i++) {
    hPtH_MC->SetBinContent(i,(hPtHv[iztt]->GetBinContent(i)+hPtHv[ittbar]->GetBinContent(i)+hPtHv[iewk]->GetBinContent(i)+hPtHv[issfake]->GetBinContent(i)));
    hPtH_MC->SetBinError(i,(0.205*hPtHv[iztt]->GetBinContent(i)+0.215*hPtHv[ittbar]->GetBinContent(i)+0.252*hPtHv[iewk]->GetBinContent(i)+0.362*hPtHv[issfake]->GetBinContent(i)));
  }
  plotPtHRatio.AddHist1D(hPtHv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotPtHRatio.AddToStack(hPtHv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtHRatio.AddToStack(hPtHv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtHRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotPtHRatio.DrawRatioStack(cr,hPtHv[0],hPtH_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotPtH_nobtag("ditaupt_nobtag","","p_{T}(#tau#tau) [GeV]",ylabel);
  if(hPtH_nobtagv[iewk] && hPtH_nobtagv[izmm]) hPtH_nobtagv[iewk]->Add(hPtH_nobtagv[izmm]);
  if(hPtH_nobtagv[ism_vbf] && hPtH_nobtagv[ism_gf] && hPtH_nobtagv[ism_vtth]) {
    hPtH_nobtagv[ism_vbf]->Add(hPtH_nobtagv[ism_gf]);
    hPtH_nobtagv[ism_vbf]->Add(hPtH_nobtagv[ism_vtth]);
  }
  if(hasData) { plotPtH_nobtag.AddHist1D(hPtH_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hPtH_nobtagv[isam]->Scale(5.);
      plotPtH_nobtag.AddToStack(hPtH_nobtagv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotPtH_nobtag.AddToStack(hPtH_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotPtH_nobtag.SetLegend(0.5,0.6,0.95,0.9);
  plotPtH_nobtag.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotHDiJetDPhi("ditaujjdphi","","#Delta#phi (#tau#tau,jj)",ylabel);
  if(hHDiJetDPhiv[iewk] && hHDiJetDPhiv[izmm]) hHDiJetDPhiv[iewk]->Add(hHDiJetDPhiv[izmm]);
  if(hHDiJetDPhiv[ism_vbf] && hHDiJetDPhiv[ism_gf] && hHDiJetDPhiv[ism_vtth]) {
    hHDiJetDPhiv[ism_vbf]->Add(hHDiJetDPhiv[ism_gf]);
    hHDiJetDPhiv[ism_vbf]->Add(hHDiJetDPhiv[ism_vtth]);
  }
  if(hasData) { plotHDiJetDPhi.AddHist1D(hHDiJetDPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hHDiJetDPhiv[isam]->Scale(5.);
      plotHDiJetDPhi.AddToStack(hHDiJetDPhiv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotHDiJetDPhi.AddToStack(hHDiJetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotHDiJetDPhi.SetLegend(0.5,0.6,0.95,0.9);
  plotHDiJetDPhi.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotHDiJetDPhiRatio("ditaujjdphi_ratio","","#Delta#phi (#tau#tau,jj)",ylabel);
  TH1F *hHDiJetDPhi_MC = new TH1F("hHDiJetDPhi_MC","",hHDiJetDPhiv[0]->GetNbinsX(),hHDiJetDPhiv[0]->GetXaxis()->GetXmin(),hHDiJetDPhiv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hHDiJetDPhiv[0]->GetNbinsX()+1; i++) {
    hHDiJetDPhi_MC->SetBinContent(i,(hHDiJetDPhiv[iztt]->GetBinContent(i)+hHDiJetDPhiv[ittbar]->GetBinContent(i)+hHDiJetDPhiv[iewk]->GetBinContent(i)+hHDiJetDPhiv[issfake]->GetBinContent(i)));
    hHDiJetDPhi_MC->SetBinError(i,(0.205*hHDiJetDPhiv[iztt]->GetBinContent(i)+0.215*hHDiJetDPhiv[ittbar]->GetBinContent(i)+0.252*hHDiJetDPhiv[iewk]->GetBinContent(i)+0.362*hHDiJetDPhiv[issfake]->GetBinContent(i)));
  }
  plotHDiJetDPhiRatio.AddHist1D(hHDiJetDPhiv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotHDiJetDPhiRatio.AddToStack(hHDiJetDPhiv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotHDiJetDPhiRatio.AddToStack(hHDiJetDPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotHDiJetDPhiRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotHDiJetDPhiRatio.DrawRatioStack(cr,hHDiJetDPhiv[0],hHDiJetDPhi_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotHDiJetDPhi_nobtag("ditaujjdphi_nobtag","","#Delta#phi (#tau#tau,jj)",ylabel);
  if(hHDiJetDPhi_nobtagv[iewk] && hHDiJetDPhi_nobtagv[izmm]) hHDiJetDPhi_nobtagv[iewk]->Add(hHDiJetDPhi_nobtagv[izmm]);
  if(hHDiJetDPhi_nobtagv[ism_vbf] && hHDiJetDPhi_nobtagv[ism_gf] && hHDiJetDPhi_nobtagv[ism_vtth]) {
    hHDiJetDPhi_nobtagv[ism_vbf]->Add(hHDiJetDPhi_nobtagv[ism_gf]);
    hHDiJetDPhi_nobtagv[ism_vbf]->Add(hHDiJetDPhi_nobtagv[ism_vtth]);
  }
  if(hasData) { plotHDiJetDPhi_nobtag.AddHist1D(hHDiJetDPhi_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hHDiJetDPhi_nobtagv[isam]->Scale(5.);
      plotHDiJetDPhi_nobtag.AddToStack(hHDiJetDPhi_nobtagv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotHDiJetDPhi_nobtag.AddToStack(hHDiJetDPhi_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotHDiJetDPhi_nobtag.SetLegend(0.5,0.6,0.95,0.9);
  plotHDiJetDPhi_nobtag.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotDEta("jetdeta","","#Delta#eta(jj)",ylabel);
  if(hDEtav[iewk] && hDEtav[izmm]) hDEtav[iewk]->Add(hDEtav[izmm]);
  if(hDEtav[ism_vbf] && hDEtav[ism_gf] && hDEtav[ism_vtth]) {
    hDEtav[ism_vbf]->Add(hDEtav[ism_gf]);
    hDEtav[ism_vbf]->Add(hDEtav[ism_vtth]);
  }
  if(hasData) { plotDEta.AddHist1D(hDEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hDEtav[isam]->Scale(5.);
      plotDEta.AddToStack(hDEtav[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotDEta.AddToStack(hDEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotDEta.SetYRange(0.0, 1.8*(plotDEta.GetStack()->GetMaximum()));
  plotDEta.SetLegend(0.5,0.6,0.95,0.9);
  plotDEta.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotDEtaRatio("jetdeta_ratio","","#Delta#eta(jj)",ylabel);
  TH1F *hDEta_MC = new TH1F("hDEta_MC","",hDEtav[0]->GetNbinsX(),hDEtav[0]->GetXaxis()->GetXmin(),hDEtav[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hDEtav[0]->GetNbinsX()+1; i++) {
    hDEta_MC->SetBinContent(i,(hDEtav[iztt]->GetBinContent(i)+hDEtav[ittbar]->GetBinContent(i)+hDEtav[iewk]->GetBinContent(i)+hDEtav[issfake]->GetBinContent(i)));
    hDEta_MC->SetBinError(i,(0.205*hDEtav[iztt]->GetBinContent(i)+0.215*hDEtav[ittbar]->GetBinContent(i)+0.252*hDEtav[iewk]->GetBinContent(i)+0.362*hDEtav[issfake]->GetBinContent(i)));
  }
  plotDEtaRatio.AddHist1D(hDEtav[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotDEtaRatio.AddToStack(hDEtav[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotDEtaRatio.AddToStack(hDEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotDEtaRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotDEtaRatio.DrawRatioStack(cr,hDEtav[0],hDEta_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotDEta_nobtag("jetdeta_nobtag","","#Delta#eta(jj)",ylabel);
  if(hDEta_nobtagv[iewk] && hDEta_nobtagv[izmm]) hDEta_nobtagv[iewk]->Add(hDEta_nobtagv[izmm]);
  if(hDEta_nobtagv[ism_vbf] && hDEta_nobtagv[ism_gf] && hDEta_nobtagv[ism_vtth]) {
    hDEta_nobtagv[ism_vbf]->Add(hDEta_nobtagv[ism_gf]);
    hDEta_nobtagv[ism_vbf]->Add(hDEta_nobtagv[ism_vtth]);
  }
  if(hasData) { plotDEta_nobtag.AddHist1D(hDEta_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hDEta_nobtagv[isam]->Scale(5.);
      plotDEta_nobtag.AddToStack(hDEta_nobtagv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotDEta_nobtag.AddToStack(hDEta_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotDEta_nobtag.SetLegend(0.5,0.6,0.95,0.9);
  plotDEta_nobtag.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotEtaProd("jetEtaProd","","(#eta1*#eta2)(jj)",ylabel);
  if(hEtaProdv[iewk] && hEtaProdv[izmm]) hEtaProdv[iewk]->Add(hEtaProdv[izmm]);
  if(hEtaProdv[ism_vbf] && hEtaProdv[ism_gf] && hEtaProdv[ism_vtth]) {
    hEtaProdv[ism_vbf]->Add(hEtaProdv[ism_gf]);
    hEtaProdv[ism_vbf]->Add(hEtaProdv[ism_vtth]);
  }
  if(hasData) { plotEtaProd.AddHist1D(hEtaProdv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hEtaProdv[isam]->Scale(5.);
      plotEtaProd.AddToStack(hEtaProdv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotEtaProd.AddToStack(hEtaProdv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaProd.SetLegend(0.2,0.6,0.65,0.9);
  plotEtaProd.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotEtaProd_nobtag("jetEtaProd_nobtag","","(#eta1*#eta2)(jj)",ylabel);
  if(hEtaProd_nobtagv[iewk] && hEtaProd_nobtagv[izmm]) hEtaProd_nobtagv[iewk]->Add(hEtaProd_nobtagv[izmm]);
  if(hEtaProd_nobtagv[ism_vbf] && hEtaProd_nobtagv[ism_gf] && hEtaProd_nobtagv[ism_vtth]) {
    hEtaProd_nobtagv[ism_vbf]->Add(hEtaProd_nobtagv[ism_gf]);
    hEtaProd_nobtagv[ism_vbf]->Add(hEtaProd_nobtagv[ism_vtth]);
  }
  if(hasData) { plotEtaProd_nobtag.AddHist1D(hEtaProd_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hEtaProd_nobtagv[isam]->Scale(5.);
      plotEtaProd_nobtag.AddToStack(hEtaProd_nobtagv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotEtaProd_nobtag.AddToStack(hEtaProd_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotEtaProd_nobtag.SetLegend(0.2,0.6,0.65,0.9);
  plotEtaProd_nobtag.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBJetPt("bjetpt","","lead b-jet pt [GeV]",ylabel);
  if(hBJetPtv[iewk] && hBJetPtv[izmm]) hBJetPtv[iewk]->Add(hBJetPtv[izmm]);
  if(hasData) { plotBJetPt.AddHist1D(hBJetPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotBJetPt.AddToStack(hBJetPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBJetPt.SetYRange(0,2.0*(plotBJetPt.GetStack()->GetMaximum()));
  plotBJetPt.SetLegend(0.5,0.65,0.95,0.9);
  plotBJetPt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBJetEta("bjeteta","","lead b-jet #eta",ylabel);
  if(hBJetEtav[iewk] && hBJetEtav[izmm]) hBJetEtav[iewk]->Add(hBJetEtav[izmm]);
  if(hasData) { plotBJetEta.AddHist1D(hBJetEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotBJetEta.AddToStack(hBJetEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBJetEta.SetYRange(0,2.0*(plotBJetEta.GetStack()->GetMaximum()));
  plotBJetEta.SetLegend(0.5,0.65,0.95,0.9);
  plotBJetEta.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBJetPhi("bjetphi","","lead b-jet #phi",ylabel);
  if(hBJetPhiv[iewk] && hBJetPhiv[izmm]) hBJetPhiv[iewk]->Add(hBJetPhiv[izmm]);
  if(hasData) { plotBJetPhi.AddHist1D(hBJetPhiv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotBJetPhi.AddToStack(hBJetPhiv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBJetPhi.SetYRange(0,2.1*(plotBJetPhi.GetStack()->GetMaximum()));
  plotBJetPhi.SetLegend(0.5,0.65,0.95,0.9);
  plotBJetPhi.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBoostJetPt("boostjetpt","","Leading jet p_{T} [GeV]",ylabel);
  if(hBoostJetPtv[iewk] && hBoostJetPtv[izmm]) hBoostJetPtv[iewk]->Add(hBoostJetPtv[izmm]);
  if(hBoostJetPtv[ism_vbf] && hBoostJetPtv[ism_gf] && hBoostJetPtv[ism_vtth]) {
    hBoostJetPtv[ism_vbf]->Add(hBoostJetPtv[ism_gf]);
    hBoostJetPtv[ism_vbf]->Add(hBoostJetPtv[ism_vtth]);
  }
  if(hasData) { plotBoostJetPt.AddHist1D(hBoostJetPtv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hBoostJetPtv[isam]->Scale(5.);
      plotBoostJetPt.AddToStack(hBoostJetPtv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBoostJetPt.AddToStack(hBoostJetPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBoostJetPt.SetLegend(0.5,0.6,0.95,0.9);
  plotBoostJetPt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBoostJetPtRatio("boostjetpt_ratio","","Leading jet p_{T} [GeV]",ylabel);
  TH1F *hBoostJetPt_MC = new TH1F("hBoostJetPt_MC","",hBoostJetPtv[0]->GetNbinsX(),hBoostJetPtv[0]->GetXaxis()->GetXmin(),hBoostJetPtv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hBoostJetPtv[0]->GetNbinsX()+1; i++) {
    hBoostJetPt_MC->SetBinContent(i,(hBoostJetPtv[iztt]->GetBinContent(i)+hBoostJetPtv[ittbar]->GetBinContent(i)+hBoostJetPtv[iewk]->GetBinContent(i)+hBoostJetPtv[issfake]->GetBinContent(i)));
    hBoostJetPt_MC->SetBinError(i,(0.062*hBoostJetPtv[iztt]->GetBinContent(i)+0.092*hBoostJetPtv[ittbar]->GetBinContent(i)+0.159*hBoostJetPtv[iewk]->GetBinContent(i)+0.305*hBoostJetPtv[issfake]->GetBinContent(i)));
  }
  plotBoostJetPtRatio.AddHist1D(hBoostJetPtv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotBoostJetPtRatio.AddToStack(hBoostJetPtv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBoostJetPtRatio.AddToStack(hBoostJetPtv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBoostJetPtRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotBoostJetPtRatio.DrawRatioStack(cr,hBoostJetPtv[0],hBoostJetPt_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotBoostJetEta("boostjeteta","","Leading jet #eta",ylabel);
  if(hBoostJetEtav[iewk] && hBoostJetEtav[izmm]) hBoostJetEtav[iewk]->Add(hBoostJetEtav[izmm]);
  if(hBoostJetEtav[ism_vbf] && hBoostJetEtav[ism_gf] && hBoostJetEtav[ism_vtth]) {
    hBoostJetEtav[ism_vbf]->Add(hBoostJetEtav[ism_gf]);
    hBoostJetEtav[ism_vbf]->Add(hBoostJetEtav[ism_vtth]);
  }
  if(hasData) { plotBoostJetEta.AddHist1D(hBoostJetEtav[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hBoostJetEtav[isam]->Scale(5.);
      plotBoostJetEta.AddToStack(hBoostJetEtav[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBoostJetEta.AddToStack(hBoostJetEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBoostJetEta.SetYRange(0,2.0*(plotBoostJetEta.GetStack()->GetMaximum()));
  plotBoostJetEta.SetLegend(0.5,0.6,0.95,0.9);
  plotBoostJetEta.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBoostJetEtaRatio("boostjeteta_ratio","","Leading jet #eta",ylabel);
  TH1F *hBoostJetEta_MC = new TH1F("hBoostJetEta_MC","",hBoostJetEtav[0]->GetNbinsX(),hBoostJetEtav[0]->GetXaxis()->GetXmin(),hBoostJetEtav[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hBoostJetEtav[0]->GetNbinsX()+1; i++) {
    hBoostJetEta_MC->SetBinContent(i,(hBoostJetEtav[iztt]->GetBinContent(i)+hBoostJetEtav[ittbar]->GetBinContent(i)+hBoostJetEtav[iewk]->GetBinContent(i)+hBoostJetEtav[issfake]->GetBinContent(i)));
    hBoostJetEta_MC->SetBinError(i,(0.062*hBoostJetEtav[iztt]->GetBinContent(i)+0.092*hBoostJetEtav[ittbar]->GetBinContent(i)+0.159*hBoostJetEtav[iewk]->GetBinContent(i)+0.305*hBoostJetEtav[issfake]->GetBinContent(i)));
  }
  plotBoostJetEtaRatio.AddHist1D(hBoostJetEtav[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotBoostJetEtaRatio.AddToStack(hBoostJetEtav[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBoostJetEtaRatio.AddToStack(hBoostJetEtav[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBoostJetEtaRatio.SetYRange(0,2.0*(plotBoostJetEtaRatio.GetStack()->GetMaximum()));
  plotBoostJetEtaRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotBoostJetEtaRatio.DrawRatioStack(cr,hBoostJetEtav[0],hBoostJetEta_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotVBFJetPt1("vbfjetpt1","","Leading jet p_{T} [GeV]",ylabel);
  if(hVBFJetPt1v[iewk] && hVBFJetPt1v[izmm]) hVBFJetPt1v[iewk]->Add(hVBFJetPt1v[izmm]);
  if(hVBFJetPt1v[ism_vbf] && hVBFJetPt1v[ism_gf] && hVBFJetPt1v[ism_vtth]) {
    hVBFJetPt1v[ism_vbf]->Add(hVBFJetPt1v[ism_gf]);
    hVBFJetPt1v[ism_vbf]->Add(hVBFJetPt1v[ism_vtth]);
  }
  if(hasData) { plotVBFJetPt1.AddHist1D(hVBFJetPt1v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hVBFJetPt1v[isam]->Scale(5.);
      plotVBFJetPt1.AddToStack(hVBFJetPt1v[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotVBFJetPt1.AddToStack(hVBFJetPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetPt1.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetPt1.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotVBFJetPt1Ratio("vbfjetpt1_ratio","","Leading jet p_{T} [GeV]",ylabel);
  TH1F *hVBFJetPt1_MC = new TH1F("hVBFJetPt1_MC","",hVBFJetPt1v[0]->GetNbinsX(),hVBFJetPt1v[0]->GetXaxis()->GetXmin(),hVBFJetPt1v[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hVBFJetPt1v[0]->GetNbinsX()+1; i++) {
    hVBFJetPt1_MC->SetBinContent(i,(hVBFJetPt1v[iztt]->GetBinContent(i)+hVBFJetPt1v[ittbar]->GetBinContent(i)+hVBFJetPt1v[iewk]->GetBinContent(i)+hVBFJetPt1v[issfake]->GetBinContent(i)));
    hVBFJetPt1_MC->SetBinError(i,(0.205*hVBFJetPt1v[iztt]->GetBinContent(i)+0.215*hVBFJetPt1v[ittbar]->GetBinContent(i)+0.252*hVBFJetPt1v[iewk]->GetBinContent(i)+0.362*hVBFJetPt1v[issfake]->GetBinContent(i)));
  }
  plotVBFJetPt1Ratio.AddHist1D(hVBFJetPt1v[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotVBFJetPt1Ratio.AddToStack(hVBFJetPt1v[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotVBFJetPt1Ratio.AddToStack(hVBFJetPt1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetPt1Ratio.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetPt1Ratio.DrawRatioStack(cr,hVBFJetPt1v[0],hVBFJetPt1_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotVBFJetPt2("vbfjetpt2","","Second jet p_{T} [GeV]",ylabel);
  if(hVBFJetPt2v[iewk] && hVBFJetPt2v[izmm]) hVBFJetPt2v[iewk]->Add(hVBFJetPt2v[izmm]);
  if(hVBFJetPt2v[ism_vbf] && hVBFJetPt2v[ism_gf] && hVBFJetPt2v[ism_vtth]) {
    hVBFJetPt2v[ism_vbf]->Add(hVBFJetPt2v[ism_gf]);
    hVBFJetPt2v[ism_vbf]->Add(hVBFJetPt2v[ism_vtth]);
  }
  if(hasData) { plotVBFJetPt2.AddHist1D(hVBFJetPt2v[0],samplev[0]->label,"E"); } 
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hVBFJetPt2v[isam]->Scale(5.);
      plotVBFJetPt2.AddToStack(hVBFJetPt2v[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotVBFJetPt2.AddToStack(hVBFJetPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetPt2.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetPt2.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotVBFJetPt2Ratio("vbfjetpt2_ratio","","Second jet p_{T} [GeV]",ylabel);
  TH1F *hVBFJetPt2_MC = new TH1F("hVBFJetPt2_MC","",hVBFJetPt2v[0]->GetNbinsX(),hVBFJetPt2v[0]->GetXaxis()->GetXmin(),hVBFJetPt2v[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hVBFJetPt2v[0]->GetNbinsX()+1; i++) {
    hVBFJetPt2_MC->SetBinContent(i,(hVBFJetPt2v[iztt]->GetBinContent(i)+hVBFJetPt2v[ittbar]->GetBinContent(i)+hVBFJetPt2v[iewk]->GetBinContent(i)+hVBFJetPt2v[issfake]->GetBinContent(i)));
    hVBFJetPt2_MC->SetBinError(i,(0.205*hVBFJetPt2v[iztt]->GetBinContent(i)+0.215*hVBFJetPt2v[ittbar]->GetBinContent(i)+0.252*hVBFJetPt2v[iewk]->GetBinContent(i)+0.362*hVBFJetPt2v[issfake]->GetBinContent(i)));
  }
  plotVBFJetPt2Ratio.AddHist1D(hVBFJetPt2v[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotVBFJetPt2Ratio.AddToStack(hVBFJetPt2v[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotVBFJetPt2Ratio.AddToStack(hVBFJetPt2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetPt2Ratio.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetPt2Ratio.DrawRatioStack(cr,hVBFJetPt2v[0],hVBFJetPt2_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotVBFJetEta1("vbfjeteta1","","Leading jet #eta",ylabel);
  if(hVBFJetEta1v[iewk] && hVBFJetEta1v[izmm]) hVBFJetEta1v[iewk]->Add(hVBFJetEta1v[izmm]);
  if(hVBFJetEta1v[ism_vbf] && hVBFJetEta1v[ism_gf] && hVBFJetEta1v[ism_vtth]) {
    hVBFJetEta1v[ism_vbf]->Add(hVBFJetEta1v[ism_gf]);
    hVBFJetEta1v[ism_vbf]->Add(hVBFJetEta1v[ism_vtth]);
  }
  if(hasData) { plotVBFJetEta1.AddHist1D(hVBFJetEta1v[0],samplev[0]->label,"E"); } 
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hVBFJetEta1v[isam]->Scale(5.);
      plotVBFJetEta1.AddToStack(hVBFJetEta1v[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotVBFJetEta1.AddToStack(hVBFJetEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetEta1.SetYRange(0,2.0*(plotVBFJetEta1.GetStack()->GetMaximum()));
  plotVBFJetEta1.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetEta1.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotVBFJetEta1Ratio("vbfjeteta1_ratio","","Leading jet #eta",ylabel);
  TH1F *hVBFJetEta1_MC = new TH1F("hVBFJetEta1_MC","",hVBFJetEta1v[0]->GetNbinsX(),hVBFJetEta1v[0]->GetXaxis()->GetXmin(),hVBFJetEta1v[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hVBFJetEta1v[0]->GetNbinsX()+1; i++) {
    hVBFJetEta1_MC->SetBinContent(i,(hVBFJetEta1v[iztt]->GetBinContent(i)+hVBFJetEta1v[ittbar]->GetBinContent(i)+hVBFJetEta1v[iewk]->GetBinContent(i)+hVBFJetEta1v[issfake]->GetBinContent(i)));
    hVBFJetEta1_MC->SetBinError(i,(0.205*hVBFJetEta1v[iztt]->GetBinContent(i)+0.215*hVBFJetEta1v[ittbar]->GetBinContent(i)+0.252*hVBFJetEta1v[iewk]->GetBinContent(i)+0.362*hVBFJetEta1v[issfake]->GetBinContent(i)));
  }
  plotVBFJetEta1Ratio.AddHist1D(hVBFJetEta1v[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotVBFJetEta1Ratio.AddToStack(hVBFJetEta1v[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotVBFJetEta1Ratio.AddToStack(hVBFJetEta1v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetEta1Ratio.SetYRange(0,2.0*(plotVBFJetEta1Ratio.GetStack()->GetMaximum()));
  plotVBFJetEta1Ratio.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetEta1Ratio.DrawRatioStack(cr,hVBFJetEta1v[0],hVBFJetEta1_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotVBFJetEta2("vbfjeteta2","","Second jet #eta",ylabel);
  if(hVBFJetEta2v[iewk] && hVBFJetEta2v[izmm]) hVBFJetEta2v[iewk]->Add(hVBFJetEta2v[izmm]);
  if(hVBFJetEta2v[ism_vbf] && hVBFJetEta2v[ism_gf] && hVBFJetEta2v[ism_vtth]) {
    hVBFJetEta2v[ism_vbf]->Add(hVBFJetEta2v[ism_gf]);
    hVBFJetEta2v[ism_vbf]->Add(hVBFJetEta2v[ism_vtth]);
  }
  if(hasData) { plotVBFJetEta2.AddHist1D(hVBFJetEta2v[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hVBFJetEta2v[isam]->Scale(5.);
      plotVBFJetEta2.AddToStack(hVBFJetEta2v[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotVBFJetEta2.AddToStack(hVBFJetEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetEta2.SetYRange(0,2.0*(plotVBFJetEta2.GetStack()->GetMaximum()));
  plotVBFJetEta2.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetEta2.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotVBFJetEta2Ratio("vbfjeteta2_ratio","","Second jet #eta",ylabel);
  TH1F *hVBFJetEta2_MC = new TH1F("hVBFJetEta2_MC","",hVBFJetEta2v[0]->GetNbinsX(),hVBFJetEta2v[0]->GetXaxis()->GetXmin(),hVBFJetEta2v[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hVBFJetEta2v[0]->GetNbinsX()+1; i++) {
    hVBFJetEta2_MC->SetBinContent(i,(hVBFJetEta2v[iztt]->GetBinContent(i)+hVBFJetEta2v[ittbar]->GetBinContent(i)+hVBFJetEta2v[iewk]->GetBinContent(i)+hVBFJetEta2v[issfake]->GetBinContent(i)));
    hVBFJetEta2_MC->SetBinError(i,(0.205*hVBFJetEta2v[iztt]->GetBinContent(i)+0.215*hVBFJetEta2v[ittbar]->GetBinContent(i)+0.252*hVBFJetEta2v[iewk]->GetBinContent(i)+0.362*hVBFJetEta2v[issfake]->GetBinContent(i)));
  }
  plotVBFJetEta2Ratio.AddHist1D(hVBFJetEta2v[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotVBFJetEta2Ratio.AddToStack(hVBFJetEta2v[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotVBFJetEta2Ratio.AddToStack(hVBFJetEta2v[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotVBFJetEta2Ratio.SetYRange(0,2.0*(plotVBFJetEta2Ratio.GetStack()->GetMaximum()));
  plotVBFJetEta2Ratio.SetLegend(0.5,0.6,0.95,0.9);
  plotVBFJetEta2Ratio.DrawRatioStack(cr,hVBFJetEta2v[0],hVBFJetEta2_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotBDT("vbfmva","","VBF MVA value",ylabel);
  if(hBDTv[iewk] && hBDTv[izmm]) hBDTv[iewk]->Add(hBDTv[izmm]);
  if(hasData) { plotBDT.AddHist1D(hBDTv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hBDTv[isam]->Scale(10.);
      plotBDT.AddHist1D(hBDTv[isam],"(10#times) qq"+samplev[isam]->label,"hist",603,11,0)
;
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBDT.AddToStack(hBDTv[isam],samplev[isam]->label,samplev[isam]->color,1,1
,3);
  }
  //plotBDT.AddLine(0.82,0,0.82,140,kTeal,kDashed,"MVA cut");
  plotBDT.SetLegend(0.5,0.6,0.95,0.9);
  plotBDT.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBDTRatio("vbfmva_ratio","","VBF MVA value",ylabel);
  TH1F *hBDT_MC = new TH1F("hBDT_MC","",hBDTv[0]->GetNbinsX(),hBDTv[0]->GetXaxis()->GetXmin(),hBDTv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hBDTv[0]->GetNbinsX()+1; i++) {
    hBDT_MC->SetBinContent(i,(hBDTv[iztt]->GetBinContent(i)+hBDTv[ittbar]->GetBinContent(i)+hBDTv[iewk]->GetBinContent(i)+hBDTv[issfake]->GetBinContent(i)));
    hBDT_MC->SetBinError(i,(0.205*hBDTv[iztt]->GetBinContent(i)+0.215*hBDTv[ittbar]->GetBinContent(i)+0.252*hBDTv[iewk]->GetBinContent(i)+0.362*hBDTv[issfake]->GetBinContent(i)));
  }
  plotBDTRatio.AddHist1D(hBDTv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotBDTRatio.AddHist1D(hBDTv[isam],"(10#times) qq"+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBDTRatio.AddToStack(hBDTv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBDTRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotBDTRatio.DrawRatioStack(cr,hBDTv[0],hBDT_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotBDT_log("vbfmva_log","","VBF MVA value",ylabel);
  if(hasData) { plotBDT_log.AddHist1D(hBDTv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotBDT_log.AddHist1D(hBDTv[isam],"(10#times) qq"+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBDT_log.AddToStack(hBDTv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  //plotBDT_log.AddLine(0.82,0,0.82,140,kTeal,kDashed,"MVA cut");
  plotBDT_log.SetYRange(0.2,2000.0*(plotBDT_log.GetStack()->GetMaximum()));
  plotBDT_log.SetLogy();
  plotBDT_log.SetLegend(0.5,0.6,0.95,0.9);
  plotBDT_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBDT_log_Ratio("vbfmva_log_ratio","","VBF MVA value",ylabel);
  plotBDT_log_Ratio.AddHist1D(hBDTv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotBDT_log_Ratio.AddHist1D(hBDTv[isam],"(10#times) qq"+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBDT_log_Ratio.AddToStack(hBDTv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotBDT_log_Ratio.SetYRange(0.2,2000.0*(plotBDT_log_Ratio.GetStack()->GetMaximum()));
  plotBDT_log_Ratio.SetLogy();
  plotBDT_log_Ratio.SetLegend(0.5,0.6,0.95,0.9);
  plotBDT_log_Ratio.DrawRatioStack(cr,hBDTv[0],hBDT_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"Events");
  CPlot plotBDT_nobtag("vbfmva_nobtag","","VBF MVA value",ylabel);
  if(hBDT_nobtagv[iewk] && hBDT_nobtagv[izmm]) hBDT_nobtagv[iewk]->Add(hBDT_nobtagv[izmm]);
  if(hasData) { plotBDT_nobtag.AddHist1D(hBDT_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hBDT_nobtagv[isam]->Scale(10.);
      plotBDT_nobtag.AddHist1D(hBDT_nobtagv[isam],"(10#times) qq"+samplev[isam]->label,"hist",603,11,0)
;
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBDT_nobtag.AddToStack(hBDT_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1
,3);
  }
  //plotBDT_nobtag.AddLine(0.82,0,0.82,140,kTeal,kDashed,"MVA cut");
  plotBDT_nobtag.SetLegend(0.5,0.6,0.95,0.9);
  plotBDT_nobtag.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotBDT_nobtag_log("vbfmva_nobtag_log","","VBF MVA value",ylabel);
  if(hasData) { plotBDT_nobtag_log.AddHist1D(hBDT_nobtagv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotBDT_nobtag_log.AddHist1D(hBDT_nobtagv[isam],"(10#times) qq"+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotBDT_nobtag_log.AddToStack(hBDT_nobtagv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  //plotBDT_nobtag_log.AddLine(0.82,0,0.82,140,kTeal,kDashed,"MVA cut");
  plotBDT_nobtag_log.SetYRange(0.2,2000.0*(plotBDT_nobtag_log.GetStack()->GetMaximum()));
  plotBDT_nobtag_log.SetLogy();
  plotBDT_nobtag_log.SetLegend(0.5,0.6,0.95,0.9);
  plotBDT_nobtag_log.Draw(c,ecorr==kCenter,format);

  CPlot plotNPV("nvertices_reweighted","","N_{PV}","Events");
  if(hNPVv[iewk] && hNPVv[izmm]) hNPVv[iewk]->Add(hNPVv[izmm]);
  if(hasData) { plotNPV.AddHist1D(hNPVv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotNPV.AddToStack(hNPVv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNPV.SetYRange(0,1.3*(hNPVv[0]->GetMaximum()));
  plotNPV.SetLegend(0.5,0.65,0.95,0.9);
  plotNPV.Draw(c,ecorr==kCenter,format);

  CPlot plotNPVraw("nvertices_raw","","N_{PV}","Events");
  if(hNPVrawv[iewk] && hNPVrawv[izmm]) hNPVrawv[iewk]->Add(hNPVrawv[izmm]);
  if(hasData) { plotNPVraw.AddHist1D(hNPVrawv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    plotNPVraw.AddToStack(hNPVrawv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotNPVraw.SetYRange(0,1.3*(hNPVrawv[0]->GetMaximum()));
  plotNPVraw.SetLegend(0.5,0.65,0.95,0.9);
  plotNPVraw.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass("svfitmass_lept","","m_{#tau#tau} [GeV]",ylabel);
  if(hMassv[iewk] && hMassv[izmm]) hMassv[iewk]->Add(hMassv[izmm]);
  if(hMassv[ism_vbf] && hMassv[ism_gf] && hMassv[ism_vtth]) {
    hMassv[ism_vbf]->Add(hMassv[ism_gf]);
    hMassv[ism_vbf]->Add(hMassv[ism_vtth]);
  }
  if(hasData) {hMassv[0]->Scale(1.,"width"); plotMass.AddHist1D(hMassv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassv[isam]->Scale(1.,"width");
    if(snamev[isam].Contains("htt_vbf")) {
      hMassv[isam]->Scale(5.);
      plotMass.AddToStack(hMassv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass.AddToStack(hMassv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass.SetLegend(0.5,0.6,0.95,0.9);
  plotMass.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis("vismass_lept","","m_{vis} [GeV]",ylabel);
  if(hMassVisv[iewk] && hMassVisv[izmm]) hMassVisv[iewk]->Add(hMassVisv[izmm]);
  if(hMassVisv[ism_vbf] && hMassVisv[ism_gf] && hMassVisv[ism_vtth]) {
    hMassVisv[ism_vbf]->Add(hMassVisv[ism_gf]);
    hMassVisv[ism_vbf]->Add(hMassVisv[ism_vtth]);
  }
  if(hasData) {hMassVisv[0]->Scale(1.,"width"); plotMassVis.AddHist1D(hMassVisv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassVisv[isam]->Scale(1.,"width");
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVisv[isam]->Scale(5.);
      plotMassVis.AddToStack(hMassVisv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis.AddToStack(hMassVisv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotMassHigh("svfitmass_lept-high","","m_{#tau#tau} [GeV]",ylabel);
  if(hMassHighv[iewk] && hMassHighv[izmm]) hMassHighv[iewk]->Add(hMassHighv[izmm]);
  if(hasData) { plotMassHigh.AddHist1D(hMassHighv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    plotMassHigh.AddToStack(hMassHighv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassHigh.SetLogy();
  plotMassHigh.SetLegend(0.5,0.65,0.95,0.9);
  plotMassHigh.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"Events");
  CPlot plotMassVisHigh("vismass_lept-high","","m_{vis} [GeV]",ylabel);
  if(hMassVisHighv[iewk] && hMassVisHighv[izmm]) hMassVisHighv[iewk]->Add(hMassVisHighv[izmm]);
  if(hasData) { plotMassVisHigh.AddHist1D(hMassVisHighv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    plotMassVisHigh.AddToStack(hMassVisHighv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  } 
  plotMassVisHigh.SetLogy();
  plotMassVisHigh.SetLegend(0.5,0.65,0.95,0.9);
  plotMassVisHigh.Draw(c,ecorr==kCenter,format);

  /*CPlot plotMassVpmiss("massvpmiss","","m_{#tau#tau} [GeV]}","#slash{p}_{#zeta} [GeV]}");
  assert(hMassVpmissv[0]);
  plotMassVpmiss.AddHist2D((TH2D*)hMassVpmissv[0],"surf",kWhite,kBlue);*/

  /*CPlot plotProjVisVProjMet_ztt("pvisvpmiss_ztt","","p_{#zeta}^{vis} [GeV]","#slash{p}_{#zeta} [GeV]");
  assert(hProjVisVProjMetv[iemb]);
  plotProjVisVProjMet_ztt.AddHist2D((TH2D*)hProjVisVProjMetv[iemb],"COLZ",kWhite,kBlue);
  plotProjVisVProjMet_ztt.AddLine(0.0,-25,50,17.5,902);
  plotProjVisVProjMet_ztt.Draw(c,ecorr==kCenter,format);

  CPlot plotProjVisVProjMet_ttbar("pvisvpmiss_ttbar","","p_{#zeta}^{vis} [GeV]","#slash{p}_{#zeta} [GeV]");
  assert(hProjVisVProjMetv[ittbar]);
  plotProjVisVProjMet_ttbar.AddHist2D((TH2D*)hProjVisVProjMetv[ittbar],"COLZ",kWhite,kBlue);
  plotProjVisVProjMet_ttbar.AddLine(0.0,-25,50,17.5,902);
  plotProjVisVProjMet_ttbar.Draw(c,ecorr==kCenter,format);*/



  //----------------------------------------------------------------------------------------
  // inclusive
  //----------------------------------------------------------------------------------------  
  
  sprintf(ylabel,"Events");
  CPlot plotMass_i("svfitmass_incl","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_iv[iewk] && hMass_iv[izmm]) hMass_iv[iewk]->Add(hMass_iv[izmm]);
  if(hMass_iv[ism_vbf] && hMass_iv[ism_gf] && hMass_iv[ism_vtth]) {
    hMass_iv[ism_vbf]->Add(hMass_iv[ism_gf]);
    hMass_iv[ism_vbf]->Add(hMass_iv[ism_vtth]);
  }
  if(hasData) {//hMass_iv[0]->Scale(1.,"width"); 
    plotMass_i.AddHist1D(hMass_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    //hMass_iv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_iv[isam]->Scale(5.);
      plotMass_i.AddToStack(hMass_iv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_i.AddToStack(hMass_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_i.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_i.SetYRange(0.0,1.5*(plotMass_i.GetStack()->GetMaximum()));
  plotMass_i.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_iRatio("svfitmass_incl_ratio","","m_{#tau#tau} [GeV]",ylabel);
  TH1F *hMass_i_MC = new TH1F("hMass_i_MC","",hMass_iv[0]->GetNbinsX(),hMass_iv[0]->GetXaxis()->GetXmin(),hMass_iv[0]->GetXaxis()->GetXmax());
  for(Int_t i=1; i < hMass_iv[0]->GetNbinsX()+1; i++) {
    hMass_i_MC->SetBinContent(i,(hMass_iv[iztt]->GetBinContent(i)+hMass_iv[ittbar]->GetBinContent(i)+hMass_iv[iewk]->GetBinContent(i)+hMass_iv[issfake]->GetBinContent(i)));
    hMass_i_MC->SetBinError(i,(0.043*hMass_iv[iztt]->GetBinContent(i)+0.08*hMass_iv[ittbar]->GetBinContent(i)+0.153*hMass_iv[iewk]->GetBinContent(i)+0.301*hMass_iv[issfake]->GetBinContent(i)));
  }
  plotMass_iRatio.AddHist1D(hMass_iv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_iRatio.AddToStack(hMass_iv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_iRatio.AddToStack(hMass_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_iRatio.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_iRatio.DrawRatioStack(cr,hMass_iv[0],hMass_i_MC,ecorr==kCenter,format);
    
  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_i("vismass_incl","","m_{vis} [GeV]",ylabel);
  if(hMassVis_iv[iewk] && hMassVis_iv[izmm]) hMassVis_iv[iewk]->Add(hMassVis_iv[izmm]);
  if(hMassVis_iv[ism_vbf] && hMassVis_iv[ism_gf] && hMassVis_iv[ism_vtth]) {
    hMassVis_iv[ism_vbf]->Add(hMassVis_iv[ism_gf]);
    hMassVis_iv[ism_vbf]->Add(hMassVis_iv[ism_vtth]);
  }
  if(hasData) {hMassVis_iv[0]->Scale(1.,"width"); plotMassVis_i.AddHist1D(hMassVis_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassVis_iv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_iv[isam]->Scale(5.);
      plotMassVis_i.AddToStack(hMassVis_iv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_i.AddToStack(hMassVis_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_i.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_i.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_i_log("svfitmass_incl_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_i_log.AddHist1D(hMass_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_i_log.AddHist1D(hMass_iv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_i_log.AddToStack(hMass_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_i_log.SetLogy();
  plotMass_i_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_i_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_i_log("vismass_incl_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_i_log.AddHist1D(hMassVis_iv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_i_log.AddHist1D(hMassVis_iv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_i_log.AddToStack(hMassVis_iv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_i_log.SetLogy();
  plotMassVis_i_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_i_log.Draw(c,ecorr==kCenter,format);

  //----------------------------------------------------------------------------------------
  // 0-jet category
  //----------------------------------------------------------------------------------------
    
  /*sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_0jet("svfitmass_class_0jet","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_0jetv[iewk] && hMass_0jetv[izmm]) hMass_0jetv[iewk]->Add(hMass_0jetv[izmm]);
  if(hMass_0jetv[ism_vbf] && hMass_0jetv[ism_gf] && hMass_0jetv[ism_vtth]) {
    hMass_0jetv[ism_vbf]->Add(hMass_0jetv[ism_gf]);
    hMass_0jetv[ism_vbf]->Add(hMass_0jetv[ism_vtth]);
  }
  if(hasData) { plotMass_0jet.AddHist1D(hMass_0jetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_0jetv[isam]->Scale(5.);
      plotMass_0jet.AddToStack(hMass_0jetv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_0jet.AddToStack(hMass_0jetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_0jet.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_0jet.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_0jet("vismass_class_0jet","","m_{vis} [GeV]",ylabel);
  if(hMassVis_0jetv[iewk] && hMassVis_0jetv[izmm]) hMassVis_0jetv[iewk]->Add(hMassVis_0jetv[izmm]);
  if(hMassVis_0jetv[ism_vbf] && hMassVis_0jetv[ism_gf] && hMassVis_0jetv[ism_vtth]) {
    hMassVis_0jetv[ism_vbf]->Add(hMassVis_0jetv[ism_gf]);
    hMassVis_0jetv[ism_vbf]->Add(hMassVis_0jetv[ism_vtth]);
  }
  if(hasData) { plotMassVis_0jet.AddHist1D(hMassVis_0jetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_0jetv[isam]->Scale(5.);
      plotMassVis_0jet.AddToStack(hMassVis_0jetv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_0jet.AddToStack(hMassVis_0jetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_0jet.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_0jet.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_0jet_log("svfitmass_class_0jet_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_0jet_log.AddHist1D(hMass_0jetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_0jet_log.AddHist1D(hMass_0jetv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_0jet_log.AddToStack(hMass_0jetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_0jet_log.SetLogy();
  plotMass_0jet_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_0jet_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_0jet_log("vismass_class_0jet_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_0jet_log.AddHist1D(hMassVis_0jetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_0jet_log.AddHist1D(hMassVis_0jetv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_0jet_log.AddToStack(hMassVis_0jetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_0jet_log.SetLogy();
  plotMassVis_0jet_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_0jet_log.Draw(c,ecorr==kCenter,format);*/

  //----------------------------------------------------------------------------------------
  // 0-jet category, low-pt
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_0jet_lowpt("svfitmass_class_0jet_lowpt","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_0jet_lowptv[iewk_7TeV] && hMass_0jet_lowptv[izmm]) hMass_0jet_lowptv[iewk_7TeV]->Add(hMass_0jet_lowptv[izmm]);
  if(hMass_0jet_lowptv[ism_vbf] && hMass_0jet_lowptv[ism_gf] && hMass_0jet_lowptv[ism_vtth]) {
    hMass_0jet_lowptv[ism_vbf]->Add(hMass_0jet_lowptv[ism_gf]);
    hMass_0jet_lowptv[ism_vbf]->Add(hMass_0jet_lowptv[ism_vtth]);
  }
  if(hasData) {hMass_0jet_lowptv[0]->Scale(1.,"width"); plotMass_0jet_lowpt.AddHist1D(hMass_0jet_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMass_0jet_lowptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_0jet_lowptv[isam]->Scale(5.);
      plotMass_0jet_lowpt.AddToStack(hMass_0jet_lowptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_0jet_lowpt.AddToStack(hMass_0jet_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_0jet_lowpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_0jet_lowpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_0jet_lowpt("vismass_class_0jet_lowpt","","m_{vis} [GeV]",ylabel);
  if(hMassVis_0jet_lowptv[iewk_7TeV] && hMassVis_0jet_lowptv[izmm]) hMassVis_0jet_lowptv[iewk_7TeV]->Add(hMassVis_0jet_lowptv[izmm]);
  if(hMassVis_0jet_lowptv[ism_vbf] && hMassVis_0jet_lowptv[ism_gf] && hMassVis_0jet_lowptv[ism_vtth]) {
    hMassVis_0jet_lowptv[ism_vbf]->Add(hMassVis_0jet_lowptv[ism_gf]);
    hMassVis_0jet_lowptv[ism_vbf]->Add(hMassVis_0jet_lowptv[ism_vtth]);
  }
  if(hasData) {hMassVis_0jet_lowptv[0]->Scale(1.,"width");  plotMassVis_0jet_lowpt.AddHist1D(hMassVis_0jet_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassVis_0jet_lowptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_0jet_lowptv[isam]->Scale(5.);
      plotMassVis_0jet_lowpt.AddToStack(hMassVis_0jet_lowptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_0jet_lowpt.AddToStack(hMassVis_0jet_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_0jet_lowpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_0jet_lowpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_0jet_lowpt_log("svfitmass_class_0jet_lowpt_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_0jet_lowpt_log.AddHist1D(hMass_0jet_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_0jet_lowpt_log.AddHist1D(hMass_0jet_lowptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_0jet_lowpt_log.AddToStack(hMass_0jet_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_0jet_lowpt_log.SetLogy();
  plotMass_0jet_lowpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_0jet_lowpt_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_0jet_lowpt_log("vismass_class_0jet_lowpt_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_0jet_lowpt_log.AddHist1D(hMassVis_0jet_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_0jet_lowpt_log.AddHist1D(hMassVis_0jet_lowptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_0jet_lowpt_log.AddToStack(hMassVis_0jet_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_0jet_lowpt_log.SetLogy();
  plotMassVis_0jet_lowpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_0jet_lowpt_log.Draw(c,ecorr==kCenter,format);

  //----------------------------------------------------------------------------------------
  // 0-jet category, high-pt
  //----------------------------------------------------------------------------------------
    
  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_0jet_highpt("svfitmass_class_0jet_highpt","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_0jet_highptv[iewk_7TeV] && hMass_0jet_highptv[izmm]) hMass_0jet_highptv[iewk_7TeV]->Add(hMass_0jet_highptv[izmm]);
  if(hMass_0jet_highptv[ism_vbf] && hMass_0jet_highptv[ism_gf] && hMass_0jet_highptv[ism_vtth]) {
    hMass_0jet_highptv[ism_vbf]->Add(hMass_0jet_highptv[ism_gf]);
    hMass_0jet_highptv[ism_vbf]->Add(hMass_0jet_highptv[ism_vtth]);
  }
  if(hasData) {hMass_0jet_highptv[0]->Scale(1.,"width"); plotMass_0jet_highpt.AddHist1D(hMass_0jet_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMass_0jet_highptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_0jet_highptv[isam]->Scale(5.);
      plotMass_0jet_highpt.AddToStack(hMass_0jet_highptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_0jet_highpt.AddToStack(hMass_0jet_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_0jet_highpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_0jet_highpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{vis}");
  CPlot plotMassVis_0jet_highpt("vismass_class_0jet_highpt","","m_{vis} [GeV]",ylabel);
  if(hMassVis_0jet_highptv[iewk_7TeV] && hMassVis_0jet_highptv[izmm]) hMassVis_0jet_highptv[iewk_7TeV]->Add(hMassVis_0jet_highptv[izmm]);
  if(hMassVis_0jet_highptv[ism_vbf] && hMassVis_0jet_highptv[ism_gf] && hMassVis_0jet_highptv[ism_vtth]) {
    hMassVis_0jet_highptv[ism_vbf]->Add(hMassVis_0jet_highptv[ism_gf]);
    hMassVis_0jet_highptv[ism_vbf]->Add(hMassVis_0jet_highptv[ism_vtth]);
  }
  if(hasData) {hMassVis_0jet_highptv[0]->Scale(1.,"width"); plotMassVis_0jet_highpt.AddHist1D(hMassVis_0jet_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassVis_0jet_highptv[isam]->Scale(1.,"width");
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_0jet_highptv[isam]->Scale(5.);
      plotMassVis_0jet_highpt.AddToStack(hMassVis_0jet_highptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_0jet_highpt.AddToStack(hMassVis_0jet_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_0jet_highpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_0jet_highpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_0jet_highpt_log("svfitmass_class_0jet_highpt_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_0jet_highpt_log.AddHist1D(hMass_0jet_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_0jet_highpt_log.AddHist1D(hMass_0jet_highptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_0jet_highpt_log.AddToStack(hMass_0jet_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_0jet_highpt_log.SetLogy();
  plotMass_0jet_highpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_0jet_highpt_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_0jet_highpt_log("vismass_class_0jet_highpt_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_0jet_highpt_log.AddHist1D(hMassVis_0jet_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_0jet_highpt_log.AddHist1D(hMassVis_0jet_highptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_0jet_highpt_log.AddToStack(hMassVis_0jet_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_0jet_highpt_log.SetLogy();
  plotMassVis_0jet_highpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_0jet_highpt_log.Draw(c,ecorr==kCenter,format);

  //----------------------------------------------------------------------------------------
  // 1 boosted jet, no btag
  //----------------------------------------------------------------------------------------

  /*sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_boost("svfitmass_class_boost","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_boostv[iewk] && hMass_boostv[izmm]) hMass_boostv[iewk]->Add(hMass_boostv[izmm]);
  if(hMass_boostv[ism_vbf] && hMass_boostv[ism_gf] && hMass_boostv[ism_vtth]) {
    hMass_boostv[ism_vbf]->Add(hMass_boostv[ism_gf]);
    hMass_boostv[ism_vbf]->Add(hMass_boostv[ism_vtth]);
  }
  if(hasData) { plotMass_boost.AddHist1D(hMass_boostv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_boostv[isam]->Scale(5.);
      plotMass_boost.AddToStack(hMass_boostv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_boost.AddToStack(hMass_boostv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_boost.SetYRange(0.0,1.5*(plotMass_boost.GetStack()->GetMaximum()));
  plotMass_boost.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_boost.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_boost("vismass_class_boost","","m_{vis} [GeV]",ylabel);
  if(hMassVis_boostv[iewk] && hMassVis_boostv[izmm]) hMassVis_boostv[iewk]->Add(hMassVis_boostv[izmm]);
  if(hMassVis_boostv[ism_vbf] && hMassVis_boostv[ism_gf] && hMassVis_boostv[ism_vtth]) {
    hMassVis_boostv[ism_vbf]->Add(hMassVis_boostv[ism_gf]);
    hMassVis_boostv[ism_vbf]->Add(hMassVis_boostv[ism_vtth]);
  }
  if(hasData) { plotMassVis_boost.AddHist1D(hMassVis_boostv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_boostv[isam]->Scale(5.);
      plotMassVis_boost.AddToStack(hMassVis_boostv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_boost.AddToStack(hMassVis_boostv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_boost.SetYRange(0,1.3*(plotMassVis_boost.GetStack()->GetMaximum()));
  plotMassVis_boost.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_boost.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_boost_log("svfitmass_class_boost_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_boost_log.AddHist1D(hMass_boostv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_boost_log.AddHist1D(hMass_boostv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_boost_log.AddToStack(hMass_boostv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_boost_log.SetLogy();
  plotMass_boost_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_boost_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_boost_log("vismass_class_boost_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_boost_log.AddHist1D(hMassVis_boostv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_boost_log.AddHist1D(hMassVis_boostv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_boost_log.AddToStack(hMassVis_boostv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_boost_log.SetLogy();
  plotMassVis_boost_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_boost_log.Draw(c,ecorr==kCenter,format);*/

  //----------------------------------------------------------------------------------------
  // 1 jet, low pt
  //----------------------------------------------------------------------------------------

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_boost_lowpt("svfitmass_class_boost_lowpt","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_boost_lowptv[iewk] && hMass_boost_lowptv[izmm]) hMass_boost_lowptv[iewk]->Add(hMass_boost_lowptv[izmm]);
  if(hMass_boost_lowptv[ism_vbf] && hMass_boost_lowptv[ism_gf] && hMass_boost_lowptv[ism_vtth]) {
    hMass_boost_lowptv[ism_vbf]->Add(hMass_boost_lowptv[ism_gf]);
    hMass_boost_lowptv[ism_vbf]->Add(hMass_boost_lowptv[ism_vtth]);
  }
  if(hasData) {hMass_boost_lowptv[0]->Scale(1.,"width");  plotMass_boost_lowpt.AddHist1D(hMass_boost_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMass_boost_lowptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_boost_lowptv[isam]->Scale(5.);
      plotMass_boost_lowpt.AddToStack(hMass_boost_lowptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_boost_lowpt.AddToStack(hMass_boost_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_boost_lowpt.SetYRange(0.0,1.5*(plotMass_boost_lowpt.GetStack()->GetMaximum()));
  plotMass_boost_lowpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_boost_lowpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_boost_lowpt("vismass_class_boost_lowpt","","m_{vis} [GeV]",ylabel);
  if(hMassVis_boost_lowptv[iewk] && hMassVis_boost_lowptv[izmm]) hMassVis_boost_lowptv[iewk]->Add(hMassVis_boost_lowptv[izmm]);
  if(hMassVis_boost_lowptv[ism_vbf] && hMassVis_boost_lowptv[ism_gf] && hMassVis_boost_lowptv[ism_vtth]) {
    hMassVis_boost_lowptv[ism_vbf]->Add(hMassVis_boost_lowptv[ism_gf]);
    hMassVis_boost_lowptv[ism_vbf]->Add(hMassVis_boost_lowptv[ism_vtth]);
  }
  if(hasData) {hMassVis_boost_lowptv[0]->Scale(1.,"width");  plotMassVis_boost_lowpt.AddHist1D(hMassVis_boost_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassVis_boost_lowptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_boost_lowptv[isam]->Scale(5.);
      plotMassVis_boost_lowpt.AddToStack(hMassVis_boost_lowptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_boost_lowpt.AddToStack(hMassVis_boost_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_boost_lowpt.SetYRange(0,1.3*(plotMassVis_boost_lowpt.GetStack()->GetMaximum()));
  plotMassVis_boost_lowpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_boost_lowpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_boost_lowpt_log("svfitmass_class_boost_lowpt_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_boost_lowpt_log.AddHist1D(hMass_boost_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_boost_lowpt_log.AddHist1D(hMass_boost_lowptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_boost_lowpt_log.AddToStack(hMass_boost_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_boost_lowpt_log.SetLogy();
  plotMass_boost_lowpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_boost_lowpt_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_boost_lowpt_log("vismass_class_boost_lowpt_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_boost_lowpt_log.AddHist1D(hMassVis_boost_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_boost_lowpt_log.AddHist1D(hMassVis_boost_lowptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_boost_lowpt_log.AddToStack(hMassVis_boost_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_boost_lowpt_log.SetLogy();
  plotMassVis_boost_lowpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_boost_lowpt_log.Draw(c,ecorr==kCenter,format);

  //----------------------------------------------------------------------------------------
  // 1 jet, high pt
  //----------------------------------------------------------------------------------------

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_boost_highpt("svfitmass_class_boost_highpt","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_boost_highptv[iewk_7TeV] && hMass_boost_highptv[izmm]) hMass_boost_highptv[iewk_7TeV]->Add(hMass_boost_highptv[izmm]);
  if(hMass_boost_highptv[ism_vbf] && hMass_boost_highptv[ism_gf] && hMass_boost_highptv[ism_vtth]) {
    hMass_boost_highptv[ism_vbf]->Add(hMass_boost_highptv[ism_gf]);
    hMass_boost_highptv[ism_vbf]->Add(hMass_boost_highptv[ism_vtth]);
  }
  if(hasData) {hMass_boost_highptv[0]->Scale(1.,"width");  plotMass_boost_highpt.AddHist1D(hMass_boost_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMass_boost_highptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_boost_highptv[isam]->Scale(5.);
      plotMass_boost_highpt.AddToStack(hMass_boost_highptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_boost_highpt.AddToStack(hMass_boost_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_boost_highpt.SetYRange(0.0,1.5*(plotMass_boost_highpt.GetStack()->GetMaximum()));
  plotMass_boost_highpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_boost_highpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_boost_highpt("vismass_class_boost_highpt","","m_{vis} [GeV]",ylabel);
  if(hMassVis_boost_highptv[iewk_7TeV] && hMassVis_boost_highptv[izmm]) hMassVis_boost_highptv[iewk_7TeV]->Add(hMassVis_boost_highptv[izmm]);
  if(hMassVis_boost_highptv[ism_vbf] && hMassVis_boost_highptv[ism_gf] && hMassVis_boost_highptv[ism_vtth]) {
    hMassVis_boost_highptv[ism_vbf]->Add(hMassVis_boost_highptv[ism_gf]);
    hMassVis_boost_highptv[ism_vbf]->Add(hMassVis_boost_highptv[ism_vtth]);
  }
  if(hasData) {hMassVis_boost_highptv[0]->Scale(1.,"width");  plotMassVis_boost_highpt.AddHist1D(hMassVis_boost_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassVis_boost_highptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_boost_highptv[isam]->Scale(5.);
      plotMassVis_boost_highpt.AddToStack(hMassVis_boost_highptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_boost_highpt.AddToStack(hMassVis_boost_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_boost_highpt.SetYRange(0,1.3*(plotMassVis_boost_highpt.GetStack()->GetMaximum()));
  plotMassVis_boost_highpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_boost_highpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_boost_highpt_log("svfitmass_class_boost_highpt_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_boost_highpt_log.AddHist1D(hMass_boost_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_boost_highpt_log.AddHist1D(hMass_boost_highptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_boost_highpt_log.AddToStack(hMass_boost_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_boost_highpt_log.SetLogy();
  plotMass_boost_highpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_boost_highpt_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_boost_highpt_log("vismass_class_boost_highpt_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_boost_highpt_log.AddHist1D(hMassVis_boost_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_boost_highpt_log.AddHist1D(hMassVis_boost_highptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_boost_highpt_log.AddToStack(hMassVis_boost_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_boost_highpt_log.SetLogy();
  plotMassVis_boost_highpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_boost_highpt_log.Draw(c,ecorr==kCenter,format);

  //----------------------------------------------------------------------------------------
  // 1 btagged jet
  //----------------------------------------------------------------------------------------

  /*sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_b("svfitmass_class_b","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_bv[iewk] && hMass_bv[izmm]) hMass_bv[iewk]->Add(hMass_bv[izmm]);
  if(hMass_bv[ism_vbf] && hMass_bv[ism_gf] && hMass_bv[ism_vtth]) {
    hMass_bv[ism_vbf]->Add(hMass_bv[ism_gf]);
    hMass_bv[ism_vbf]->Add(hMass_bv[ism_vtth]);
  }
  if(hasData) { plotMass_b.AddHist1D(hMass_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_bv[isam]->Scale(5.);
      plotMass_b.AddToStack(hMass_bv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_b.AddToStack(hMass_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b.SetYRange(0.0,1.5*(plotMass_b.GetStack()->GetMaximum()));
  plotMass_b.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_b.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_b("vismass_class_b","","m_{vis} [GeV]",ylabel);
  if(hMassVis_bv[iewk] && hMassVis_bv[izmm]) hMassVis_bv[iewk]->Add(hMassVis_bv[izmm]);
  if(hMassVis_bv[ism_vbf] && hMassVis_bv[ism_gf] && hMassVis_bv[ism_vtth]) {
    hMassVis_bv[ism_vbf]->Add(hMassVis_bv[ism_gf]);
    hMassVis_bv[ism_vbf]->Add(hMassVis_bv[ism_vtth]);
  }
  if(hasData) { plotMassVis_b.AddHist1D(hMassVis_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_bv[isam]->Scale(5.);
      plotMassVis_b.AddToStack(hMassVis_bv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_b.AddToStack(hMassVis_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_b.SetYRange(0,1.3*(plotMassVis_b.GetStack()->GetMaximum()));
  plotMassVis_b.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_b.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_b_log("svfitmass_class_b_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_b_log.AddHist1D(hMass_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_b_log.AddHist1D(hMass_bv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_b_log.AddToStack(hMass_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b_log.SetLogy();
  plotMass_b_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_b_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_b_log("vismass_class_b_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_b_log.AddHist1D(hMassVis_bv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_b_log.AddHist1D(hMassVis_bv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_b_log.AddToStack(hMassVis_bv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_b_log.SetLogy();
  plotMassVis_b_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_b_log.Draw(c,ecorr==kCenter,format);*/

  //----------------------------------------------------------------------------------------
  // 1 btagged jet, low pT
  //----------------------------------------------------------------------------------------

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_b_lowpt("svfitmass_class_b_lowpt","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_b_lowptv[iewk] && hMass_b_lowptv[izmm]) hMass_b_lowptv[iewk]->Add(hMass_b_lowptv[izmm]);
  if(hMass_b_lowptv[ism_vbf] && hMass_b_lowptv[ism_gf] && hMass_b_lowptv[ism_vtth]) {
    hMass_b_lowptv[ism_vbf]->Add(hMass_b_lowptv[ism_gf]);
    hMass_b_lowptv[ism_vbf]->Add(hMass_b_lowptv[ism_vtth]);
  }
  if(hasData) {hMass_b_lowptv[0]->Scale(1.,"width"); plotMass_b_lowpt.AddHist1D(hMass_b_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMass_b_lowptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_b_lowptv[isam]->Scale(5.);
      plotMass_b_lowpt.AddToStack(hMass_b_lowptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_b_lowpt.AddToStack(hMass_b_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b_lowpt.SetYRange(0.0,1.5*(plotMass_b_lowpt.GetStack()->GetMaximum()));
  plotMass_b_lowpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_b_lowpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_b_lowpt("vismass_class_b_lowpt","","m_{vis} [GeV]",ylabel);
  if(hMassVis_b_lowptv[iewk] && hMassVis_b_lowptv[izmm]) hMassVis_b_lowptv[iewk]->Add(hMassVis_b_lowptv[izmm]);
  if(hMassVis_b_lowptv[ism_vbf] && hMassVis_b_lowptv[ism_gf] && hMassVis_b_lowptv[ism_vtth]) {
    hMassVis_b_lowptv[ism_vbf]->Add(hMassVis_b_lowptv[ism_gf]);
    hMassVis_b_lowptv[ism_vbf]->Add(hMassVis_b_lowptv[ism_vtth]);
  }
  if(hasData) {hMassVis_b_lowptv[0]->Scale(1.,"width");  plotMassVis_b_lowpt.AddHist1D(hMassVis_b_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassVis_b_lowptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_b_lowptv[isam]->Scale(5.);
      plotMassVis_b_lowpt.AddToStack(hMassVis_b_lowptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_b_lowpt.AddToStack(hMassVis_b_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_b_lowpt.SetYRange(0,1.3*(plotMassVis_b_lowpt.GetStack()->GetMaximum()));
  plotMassVis_b_lowpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_b_lowpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_b_lowpt_log("svfitmass_class_b_lowpt_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_b_lowpt_log.AddHist1D(hMass_b_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_b_lowpt_log.AddHist1D(hMass_b_lowptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_b_lowpt_log.AddToStack(hMass_b_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b_lowpt_log.SetLogy();
  plotMass_b_lowpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_b_lowpt_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_b_lowpt_log("vismass_class_b_lowpt_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_b_lowpt_log.AddHist1D(hMassVis_b_lowptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_b_lowpt_log.AddHist1D(hMassVis_b_lowptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_b_lowpt_log.AddToStack(hMassVis_b_lowptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_b_lowpt_log.SetLogy();
  plotMassVis_b_lowpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_b_lowpt_log.Draw(c,ecorr==kCenter,format);

  //----------------------------------------------------------------------------------------
  // 1 btagged jet, high pT
  //----------------------------------------------------------------------------------------

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_b_highpt("svfitmass_class_b_highpt","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_b_highptv[iewk] && hMass_b_highptv[izmm]) hMass_b_highptv[iewk]->Add(hMass_b_highptv[izmm]);
  if(hMass_b_highptv[ism_vbf] && hMass_b_highptv[ism_gf] && hMass_b_highptv[ism_vtth]) {
    hMass_b_highptv[ism_vbf]->Add(hMass_b_highptv[ism_gf]);
    hMass_b_highptv[ism_vbf]->Add(hMass_b_highptv[ism_vtth]);
  }
  if(hasData) {hMass_b_highptv[0]->Scale(1.,"width"); plotMass_b_highpt.AddHist1D(hMass_b_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMass_b_highptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_b_highptv[isam]->Scale(5.);
      plotMass_b_highpt.AddToStack(hMass_b_highptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_b_highpt.AddToStack(hMass_b_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b_highpt.SetYRange(0.0,1.5*(plotMass_b_highpt.GetStack()->GetMaximum()));
  plotMass_b_highpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_b_highpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_b_highpt("vismass_class_b_highpt","","m_{vis} [GeV]",ylabel);
  if(hMassVis_b_highptv[iewk] && hMassVis_b_highptv[izmm]) hMassVis_b_highptv[iewk]->Add(hMassVis_b_highptv[izmm]);
  if(hMassVis_b_highptv[ism_vbf] && hMassVis_b_highptv[ism_gf] && hMassVis_b_highptv[ism_vtth]) {
    hMassVis_b_highptv[ism_vbf]->Add(hMassVis_b_highptv[ism_gf]);
    hMassVis_b_highptv[ism_vbf]->Add(hMassVis_b_highptv[ism_vtth]);
  }
  if(hasData) {hMassVis_b_highptv[0]->Scale(1.,"width"); plotMassVis_b_highpt.AddHist1D(hMassVis_b_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    hMassVis_b_highptv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_b_highptv[isam]->Scale(5.);
      plotMassVis_b_highpt.AddToStack(hMassVis_b_highptv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_b_highpt.AddToStack(hMassVis_b_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_b_highpt.SetYRange(0,1.3*(plotMassVis_b_highpt.GetStack()->GetMaximum()));
  plotMassVis_b_highpt.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_b_highpt.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_b_highpt_log("svfitmass_class_b_highpt_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_b_highpt_log.AddHist1D(hMass_b_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMass_b_highpt_log.AddHist1D(hMass_b_highptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_b_highpt_log.AddToStack(hMass_b_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_b_highpt_log.SetLogy();
  plotMass_b_highpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_b_highpt_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_b_highpt_log("vismass_class_b_highpt_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_b_highpt_log.AddHist1D(hMassVis_b_highptv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==ifake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      plotMassVis_b_highpt_log.AddHist1D(hMassVis_b_highptv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_b_highpt_log.AddToStack(hMassVis_b_highptv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_b_highpt_log.SetLogy();
  plotMassVis_b_highpt_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_b_highpt_log.Draw(c,ecorr==kCenter,format);

  //----------------------------------------------------------------------------------------
  // VH
  //----------------------------------------------------------------------------------------

  /*sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_2jet("svfitmass_class_2jet","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_2jetv[iewk] && hMass_2jetv[izmm]) hMass_2jetv[iewk]->Add(hMass_2jetv[izmm]);
  if(hMass_2jetv[ism_vbf] && hMass_2jetv[ism_gf] && hMass_2jetv[ism_vtth]) {
    hMass_2jetv[ism_vbf]->Add(hMass_2jetv[ism_gf]);
    hMass_2jetv[ism_vbf]->Add(hMass_2jetv[ism_vtth]);
  }
  if(hasData) { plotMass_2jet.AddHist1D(hMass_2jetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_2jetv[isam]->Scale(5.);
      plotMass_2jet.AddToStack(hMass_2jetv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_2jet.AddToStack(hMass_2jetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_2jet.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_2jet.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_2jet("vismass_class_2jet","","m_{vis} [GeV]",ylabel);
  if(hMassVis_2jetv[iewk] && hMassVis_2jetv[izmm]) hMassVis_2jetv[iewk]->Add(hMassVis_2jetv[izmm]);
  if(hMassVis_2jetv[ism_vbf] && hMassVis_2jetv[ism_gf] && hMassVis_2jetv[ism_vtth]) {
    hMassVis_2jetv[ism_vbf]->Add(hMassVis_2jetv[ism_gf]);
    hMassVis_2jetv[ism_vbf]->Add(hMassVis_2jetv[ism_vtth]);
  }
  if(hasData) { plotMassVis_2jet.AddHist1D(hMassVis_2jetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_2jetv[isam]->Scale(5.);
      plotMassVis_2jet.AddToStack(hMassVis_2jetv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_2jet.AddToStack(hMassVis_2jetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_2jet.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_2jet.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_2jet_log("svfitmass_class_2jet_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_2jet_log.AddHist1D(hMass_2jetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) plotMass_2jet_log.AddHist1D(hMass_2jetv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_2jet_log.AddToStack(hMass_2jetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_2jet_log.SetLogy();
  plotMass_2jet_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_2jet_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_2jet_log("vismass_class_2jet_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_2jet_log.AddHist1D(hMassVis_2jetv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) plotMassVis_2jet_log.AddHist1D(hMassVis_2jetv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_2jet_log.AddToStack(hMassVis_2jetv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_2jet_log.SetLogy();
  plotMassVis_2jet_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_2jet_log.Draw(c,ecorr==kCenter,format);*/

  //----------------------------------------------------------------------------------------
  // vbf
  //----------------------------------------------------------------------------------------

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_vbf("svfitmass_class_vbf","","m_{#tau#tau} [GeV]",ylabel);
  if(hMass_vbfv[iewk] && hMass_vbfv[izmm]) hMass_vbfv[iewk]->Add(hMass_vbfv[izmm]);
  if(hMass_vbfv[ism_vbf] && hMass_vbfv[ism_gf] && hMass_vbfv[ism_vtth]) {
    hMass_vbfv[ism_vbf]->Add(hMass_vbfv[ism_gf]);
    hMass_vbfv[ism_vbf]->Add(hMass_vbfv[ism_vtth]);
  }
  if(hasData) {hMass_vbfv[0]->Scale(1.,"width"); plotMass_vbf.AddHist1D(hMass_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar) || (isam==issfake)) continue;
    hMass_vbfv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMass_vbfv[isam]->Scale(5.);
      plotMass_vbf.AddToStack(hMass_vbfv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_vbf.AddToStack(hMass_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_vbf.SetYRange(0.0,1.5*(plotMass_vbf.GetStack()->GetMaximum()));
  plotMass_vbf.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_vbf.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_vbf("vismass_class_vbf","","m_{vis} [GeV]",ylabel);
  if(hMassVis_vbfv[iewk] && hMassVis_vbfv[izmm]) hMassVis_vbfv[iewk]->Add(hMassVis_vbfv[izmm]);
  if(hMassVis_vbfv[ism_vbf] && hMassVis_vbfv[ism_gf] && hMassVis_vbfv[ism_vtth]) {
    hMassVis_vbfv[ism_vbf]->Add(hMassVis_vbfv[ism_gf]);
    hMassVis_vbfv[ism_vbf]->Add(hMassVis_vbfv[ism_vtth]);
  }
  if(hasData) {hMassVis_vbfv[0]->Scale(1.,"width");  plotMassVis_vbf.AddHist1D(hMassVis_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar) || (isam==issfake)) continue;
    hMassVis_vbfv[isam]->Scale(1.,"width"); 
    if(snamev[isam].Contains("htt_vbf")) {
      hMassVis_vbfv[isam]->Scale(5.);
      plotMassVis_vbf.AddToStack(hMassVis_vbfv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_vbf.AddToStack(hMassVis_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_vbf.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_vbf.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMass_vbf_log("svfitmass_class_vbf_log","","m_{#tau#tau} [GeV]",ylabel);
  if(hasData) { plotMass_vbf_log.AddHist1D(hMass_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) plotMass_vbf_log.AddHist1D(hMass_vbfv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMass_vbf_log.AddToStack(hMass_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMass_vbf_log.SetLogy();
  plotMass_vbf_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMass_vbf_log.Draw(c,ecorr==kCenter,format);

  sprintf(ylabel,"dN/dm_{#tau#tau}");
  CPlot plotMassVis_vbf_log("vismass_class_vbf_log","","m_{vis} [GeV]",ylabel);
  if(hasData) { plotMassVis_vbf_log.AddHist1D(hMassVis_vbfv[0],samplev[0]->label,"E"); }
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (isam==issfake)) continue;
    if(snamev[isam].Contains("htt_vbf")) plotMassVis_vbf_log.AddHist1D(hMassVis_vbfv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plotMassVis_vbf_log.AddToStack(hMassVis_vbfv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  plotMassVis_vbf_log.SetLogy();
  plotMassVis_vbf_log.SetLegend(0.5,0.6,0.95,0.9);
  plotMassVis_vbf_log.Draw(c,ecorr==kCenter,format);


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
  Double_t nTotal_0jet=0, nTotalVar_0jet=0;
  Double_t nTotal_0jet_lowpt=0, nTotalVar_0jet_lowpt=0;
  Double_t nTotal_0jet_highpt=0, nTotalVar_0jet_highpt=0;
  Double_t nTotal_boost=0, nTotalVar_boost=0;
  Double_t nTotal_boost_lowpt=0, nTotalVar_boost_lowpt=0;
  Double_t nTotal_boost_highpt=0, nTotalVar_boost_highpt=0;
  Double_t nTotal_b=0, nTotalVar_b=0;
  Double_t nTotal_b_lowpt=0, nTotalVar_b_lowpt=0;
  Double_t nTotal_b_highpt=0, nTotalVar_b_highpt=0;
  Double_t nTotal_2jet=0, nTotalVar_2jet=0;
  Double_t nTotal_vbf=0,   nTotalVar_vbf=0;

  if(ecorr==kCenter) {
    ofstream yieldfile;
    yieldfile.open(CPlot::sOutDir+"/yields.txt");
    
    yieldfile << setw(33) << "lepton sele." << setw(20) << "inclusive" << setw(20) << "0-jet" << setw(24) << "0-jet, low pT" << setw(22) << "0-jet, high pT" << setw(16) << "1-jet" << setw(24) << "1-jet, low pT" << setw(22) << "1-jet, high pT" << setw(16) << "btag" << setw(25) << "btag, low pT" << setw(23) << "btag, high pT" << setw(18) << "VBF" << endl;
  
    if(hMass_iv[iewk] && hMass_iv[izmm]) { // add zmm into the fakes
      nSelv[iewk]               += nSelv[izmm];                nSelVarv[iewk]             += nSelVarv[izmm];
      nSel_iv[iewk]             += nSel_iv[izmm];              nSelVar_iv[iewk]           += nSelVar_iv[izmm];
      nSel_0jetv[iewk]          += nSel_0jetv[izmm];           nSelVar_0jetv[iewk]        += nSelVar_0jetv[izmm];
      nSel_0jet_lowptv[iewk]    += nSel_0jet_lowptv[izmm];     nSelVar_0jet_lowptv[iewk]  += nSelVar_0jet_lowptv[izmm];
      nSel_0jet_highptv[iewk]   += nSel_0jet_highptv[izmm];    nSelVar_0jet_highptv[iewk] += nSelVar_0jet_highptv[izmm];
      nSel_boostv[iewk]         += nSel_boostv[izmm];          nSelVar_boostv[iewk]       += nSelVar_boostv[izmm];
      nSel_boost_lowptv[iewk]   += nSel_boost_lowptv[izmm];    nSelVar_boost_lowptv[iewk] += nSelVar_boost_lowptv[izmm];
      nSel_boost_highptv[iewk]  += nSel_boost_highptv[izmm];   nSelVar_boost_highptv[iewk]+= nSelVar_boost_highptv[izmm];
      nSel_bv[iewk]             += nSel_bv[izmm];              nSelVar_bv[iewk]           += nSelVar_bv[izmm];
      nSel_b_lowptv[iewk]       += nSel_b_lowptv[izmm];        nSelVar_b_lowptv[iewk]     += nSelVar_b_lowptv[izmm];
      nSel_b_highptv[iewk]      += nSel_b_highptv[izmm];       nSelVar_b_highptv[iewk]    += nSelVar_b_highptv[izmm];
      nSel_vbfv[iewk]             += nSel_vbfv[izmm];            nSelVar_vbfv[iewk]           += nSelVar_vbfv[izmm];
    }
  
    yieldfile << setw(15) << "fakes";
    yieldfile << setw(10) << setprecision(3) << fixed << nSelv[issfake]             << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[issfake]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jetv[issfake]        << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jetv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_lowptv[issfake]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_highptv[issfake] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boostv[issfake]       << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boostv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_lowptv[issfake] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_highptv[issfake]<< " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[issfake]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_bv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_lowptv[issfake]     << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_highptv[issfake]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[ifake]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[ifake]);
    yieldfile << endl;
    /*yieldfile << setw(15) << "ewk";
    yieldfile << setw(10) << setprecision(3) << fixed << nSelv[iewk]             << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[iewk]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[iewk]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[iewk]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jetv[iewk]        << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jetv[iewk]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_lowptv[iewk_7TeV]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_lowptv[iewk_7TeV]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_highptv[iewk_7TeV] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_highptv[iewk_7TeV]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boostv[iewk]       << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boostv[iewk]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_lowptv[iewk] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_lowptv[iewk]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_highptv[iewk_7TeV]<< " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_highptv[iewk_7TeV]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[iewk]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_bv[iewk]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_lowptv[iewk]     << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_lowptv[iewk]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_highptv[iewk]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_highptv[iewk]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[iewk]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[iewk]);
    yieldfile << endl;
    yieldfile << setw(15) << "ttbar";
    yieldfile << setw(10) << setprecision(3) << fixed << nSelv[ittbar]             << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[ittbar]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[ittbar]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[ittbar]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jetv[ittbar]        << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jetv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_lowptv[ittbar]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_highptv[ittbar] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boostv[ittbar]       << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boostv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_lowptv[ittbar] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_highptv[ittbar]<< " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[ittbar]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_bv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_lowptv[ittbar]     << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_highptv[ittbar]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[ittbar_7TeV]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[ittbar_7TeV]);
    yieldfile << endl;*/
  
  
  
    nTotal                 += nSelv[issfake] + nSelv[iewk] + nSelv[ittbar];
    nTotalVar              += nSelVarv[issfake] + nSelVarv[iewk] + nSelVarv[ittbar];
    nTotal_i               += nSel_iv[issfake] + nSel_iv[iewk] + nSel_iv[ittbar];
    nTotalVar_i            += nSelVar_iv[issfake] + nSelVar_iv[iewk] + nSelVar_iv[ittbar];
    nTotal_0jet            += nSel_0jetv[issfake] + nSel_0jetv[iewk] + nSel_0jetv[ittbar];
    nTotalVar_0jet         += nSelVar_0jetv[issfake] + nSelVar_0jetv[iewk] + nSelVar_0jetv[ittbar];
    nTotal_0jet_lowpt      += nSel_0jet_lowptv[issfake] + nSel_0jet_lowptv[iewk_7TeV] + nSel_0jet_lowptv[ittbar];
    nTotalVar_0jet_lowpt   += nSelVar_0jet_lowptv[issfake] + nSelVar_0jet_lowptv[iewk_7TeV] + nSelVar_0jet_lowptv[ittbar];
    nTotal_0jet_highpt     += nSel_0jet_highptv[issfake] + nSel_0jet_highptv[iewk_7TeV] + nSel_0jet_highptv[ittbar];
    nTotalVar_0jet_highpt  += nSelVar_0jet_highptv[issfake] + nSelVar_0jet_highptv[iewk_7TeV] + nSelVar_0jet_highptv[ittbar];
    nTotal_boost           += nSel_boostv[issfake] + nSel_boostv[iewk] + nSel_boostv[ittbar];
    nTotalVar_boost        += nSelVar_boostv[issfake] + nSelVar_boostv[iewk] + nSelVar_boostv[ittbar];
    nTotal_boost_lowpt     += nSel_boost_lowptv[issfake] + nSel_boost_lowptv[iewk] + nSel_boost_lowptv[ittbar];
    nTotalVar_boost_lowpt  += nSelVar_boost_lowptv[issfake] + nSelVar_boost_lowptv[iewk] + nSelVar_boost_lowptv[ittbar];
    nTotal_boost_highpt    += nSel_boost_highptv[issfake] + nSel_boost_highptv[iewk_7TeV] + nSel_boost_highptv[ittbar];
    nTotalVar_boost_highpt += nSelVar_boost_highptv[issfake] + nSelVar_boost_highptv[iewk_7TeV] + nSelVar_boost_highptv[ittbar];
    nTotal_b               += nSel_bv[issfake] + nSel_bv[iewk] + nSel_bv[ittbar];
    nTotalVar_b            += nSelVar_bv[issfake] + nSelVar_bv[iewk] + nSelVar_bv[ittbar];
    nTotal_b_lowpt         += nSel_b_lowptv[issfake] + nSel_b_lowptv[iewk] + nSel_b_lowptv[ittbar];
    nTotalVar_b_lowpt      += nSelVar_b_lowptv[issfake] + nSelVar_b_lowptv[iewk] + nSelVar_b_lowptv[ittbar];
    nTotal_b_highpt        += nSel_b_highptv[issfake] + nSel_b_highptv[iewk] + nSel_b_highptv[ittbar];
    nTotalVar_b_highpt     += nSelVar_b_highptv[issfake] + nSelVar_b_highptv[iewk] + nSelVar_b_highptv[ittbar];
    nTotal_2jet            += nSel_2jetv[ifake] + nSel_2jetv[iewk] + nSel_2jetv[ittbar];
    nTotalVar_2jet         += nSelVar_2jetv[ifake] + nSelVar_2jetv[iewk] + nSelVar_2jetv[ittbar];
    nTotal_vbf             += nSel_vbfv[ifake] + nSel_vbfv[iewk] + nSel_vbfv[ittbar_7TeV];
    nTotalVar_vbf          += nSelVar_vbfv[ifake] + nSelVar_vbfv[iewk] + nSelVar_vbfv[ittbar_7TeV];
 
  
    if(samplev.size()>1) { 
      for(UInt_t isam=1; isam<samplev.size(); isam++) {
        if((isam==izmm)) continue;
        //else if(snamev[isam].Contains("htt_") && snamev[isam].Contains("sm_"))  continue;
        else if((isam==ifake) || (isam==issfake)) continue;
        else {
          yieldfile << setw(15) << snamev[isam];
          yieldfile << setw(10) << setprecision(3) << fixed << nSelv[isam]             << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[isam]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jetv[isam]        << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jetv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_lowptv[isam]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_lowptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_highptv[isam] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_highptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_boostv[isam]       << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boostv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_lowptv[isam] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_lowptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_highptv[isam]<< " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_highptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[isam]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_bv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_lowptv[isam]     << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_lowptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_highptv[isam]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_highptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[isam]         << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[isam]);
          yieldfile << endl;
        }
  
        if(snamev[isam].Contains("mssm_"))  continue; // don't add higgs samples to total
        if(snamev[isam].Contains("sm_"))    continue;
        if(snamev[isam].Contains("ggH"))    continue;
        if(snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH"))    continue;
        if(snamev[isam].Contains("VH"))    continue;
        if(snamev[isam].Contains("bbH"))    continue;

        if((isam==iewk) || (isam==ittbar) || (isam==iewk_7TeV) || (isam==ittbar_7TeV)) continue;
  
        nTotal    	        += nSelv[isam];
        nTotalVar 	        += nSelVarv[isam];
        nTotal_i      	        += nSel_iv[isam];
        nTotalVar_i   	        += nSelVar_iv[isam];
        nTotal_0jet             += nSel_0jetv[isam];
        nTotalVar_0jet          += nSelVar_0jetv[isam];
        nTotal_0jet_lowpt       += nSel_0jet_lowptv[isam];
        nTotalVar_0jet_lowpt    += nSelVar_0jet_lowptv[isam];
        nTotal_0jet_highpt      += nSel_0jet_highptv[isam];
        nTotalVar_0jet_highpt   += nSelVar_0jet_highptv[isam];
        nTotal_boost            += nSel_boostv[isam];
        nTotalVar_boost         += nSelVar_boostv[isam];
        nTotal_boost_lowpt      += nSel_boost_lowptv[isam];
        nTotalVar_boost_lowpt   += nSelVar_boost_lowptv[isam];
        nTotal_boost_highpt     += nSel_boost_highptv[isam];
        nTotalVar_boost_highpt  += nSelVar_boost_highptv[isam];
        nTotal_b                += nSel_bv[isam];
        nTotalVar_b             += nSelVar_bv[isam];
        nTotal_b_lowpt          += nSel_b_lowptv[isam];
        nTotalVar_b_lowpt       += nSelVar_b_lowptv[isam];
        nTotal_b_highpt         += nSel_b_highptv[isam];
        nTotalVar_b_highpt      += nSelVar_b_highptv[isam];
        nTotal_2jet             += nSel_2jetv[isam];
        nTotalVar_2jet          += nSelVar_2jetv[isam];
        nTotal_vbf      	+= nSel_vbfv[isam];
        nTotalVar_vbf   	+= nSelVar_vbfv[isam];
      }
      yieldfile << endl;
      yieldfile << setw(15) << "bkg MC";
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal               << " +/- " << sqrt(nTotalVar);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_i             << " +/- " << sqrt(nTotalVar_i);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_0jet          << " +/- " << sqrt(nTotalVar_0jet);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_0jet_lowpt    << " +/- " << sqrt(nTotalVar_0jet_lowpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_0jet_highpt   << " +/- " << sqrt(nTotalVar_0jet_highpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_boost         << " +/- " << sqrt(nTotalVar_boost);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_boost_lowpt   << " +/- " << sqrt(nTotalVar_boost_lowpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_boost_highpt  << " +/- " << sqrt(nTotalVar_boost_highpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_b             << " +/- " << sqrt(nTotalVar_b);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_b_lowpt       << " +/- " << sqrt(nTotalVar_b_lowpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_b_highpt      << " +/- " << sqrt(nTotalVar_b_highpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_vbf           << " +/- " << sqrt(nTotalVar_vbf);
      yieldfile << endl;
    }
  
    if(hasData) {
      yieldfile << setw(15) << "Data";
      yieldfile << setw(10) << setprecision(3) << fixed << nSelv[0]             << " +/- " << sqrt(nSelVarv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[0]           << " +/- " << sqrt(nSelVar_iv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jetv[0]        << " +/- " << sqrt(nSelVar_0jetv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_lowptv[0]  << " +/- " << sqrt(nSelVar_0jet_lowptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_highptv[0] << " +/- " << sqrt(nSelVar_0jet_highptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_boostv[0]       << " +/- " << sqrt(nSelVar_boostv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_lowptv[0] << " +/- " << sqrt(nSelVar_boost_lowptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_highptv[0]<< " +/- " << sqrt(nSelVar_boost_highptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_bv[0]           << " +/- " << sqrt(nSelVar_bv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_lowptv[0]     << " +/- " << sqrt(nSelVar_b_lowptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_highptv[0]    << " +/- " << sqrt(nSelVar_b_highptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[0]         << " +/- " << sqrt(nSelVar_vbfv[0]);
      yieldfile << endl;
    }
  
    // cat the yields to stdout
    TString cmd("cat "+CPlot::sOutDir+"/yields.txt");
    system(cmd.Data());
  
    // write out the acceptances
    yieldfile << endl << endl;
    yieldfile << setw(25) << " " << setw(15) << "initial events";
    yieldfile << setw(20) << "lepton sele." <<setw(17) << "inclusive" << setw(19) << "0-jet" << setw(19) << "0-jet, low pT" << setw(19) << "0-jet, high pT" << setw(19) << "1-jet" << setw(19) << "1-jet, low pT" << setw(19) << "1-jet, high pT" << setw(19) << "btag" << setw(19) << "btag, low pT" << setw(19) << "btag, high pT" << setw(22) << "vbf" << endl;
  
    for(UInt_t isam=0; isam<samplev.size(); isam++) {
      if(!(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH"))) continue;
      if(snamev[isam].Contains("htt_")) continue;
      yieldfile << setw(25) << snamev[isam] << "  " <<
        setw(15) << nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSelv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_iv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_0jetv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_0jet_lowptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_0jet_highptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_boostv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_boost_lowptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_boost_highptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_bv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_b_lowptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_b_highptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_vbfv[isam]/nUnskEventsv[isam] << endl;
    }
    yieldfile.close();
  
    makeHTML(outputDir);
  }
  
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
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt.png\"><img src=\"plots/pt.png\" alt=\"plots/pt.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "</tr>" << endl;
  //htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtele.png\"><img src=\"plots/mtele.png\" alt=\"plots/mtele.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mtmu.png\"><img src=\"plots/mtmu.png\" alt=\"plots/mtmu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met.png\"><img src=\"plots/met.png\" alt=\"plots/met.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lepdeta.png\"><img src=\"plots/lepdeta.png\" alt=\"plots/lepdeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetamiss.png\"><img src=\"plots/pzetamiss.png\" alt=\"plots/pzetamiss.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetavis.png\"><img src=\"plots/pzetavis.png\" alt=\"plots/pzetavis.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetavar.png\"><img src=\"plots/pzetavar.png\" alt=\"plots/pzetavar.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetpt1.png\"><img src=\"plots/jetpt1.png\" alt=\"plots/jetpt1.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetpt2.png\"><img src=\"plots/jetpt2.png\" alt=\"plots/jetpt2.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jeteta1.png\"><img src=\"plots/jeteta1.png\" alt=\"plots/jeteta1.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jeteta2.png\"><img src=\"plots/jeteta2.png\" alt=\"plots/jeteta2.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "</tr>" << endl;
  //htmlfile << "<tr>" << endl;
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
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ptjj.png\"><img src=\"plots/ptjj.png\" alt=\"plots/ptjj.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ditaupt.png\"><img src=\"plots/ditaupt.png\" alt=\"plots/ditaupt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ditaujjdphi.png\"><img src=\"plots/ditaujjdphi.png\" alt=\"plots/ditaujjdphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetdphi_nobtag.png\"><img src=\"plots/jetdphi_nobtag.png\" alt=\"plots/jetdphi_nobtag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mjj_nobtag.png\"><img src=\"plots/mjj_nobtag.png\" alt=\"plots/mjj_nobtag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetdeta_nobtag.png\"><img src=\"plots/jetdeta_nobtag.png\" alt=\"plots/jetdeta_nobtag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetEtaProd_nobtag.png\"><img src=\"plots/jetEtaProd_nobtag.png\" alt=\"plots/jetEtaProd_nobtag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ptjj_nobtag.png\"><img src=\"plots/ptjj_nobtag.png\" alt=\"plots/ptjj_nobtag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ditaupt_nobtag.png\"><img src=\"plots/ditaupt_nobtag.png\" alt=\"plots/ditaupt_nobtag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/ditaujjdphi_nobtag.png\"><img src=\"plots/ditaujjdphi_nobtag.png\" alt=\"plots/ditaujjdphi_nobtag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
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
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vbfmva_log.png\"><img src=\"plots/vbfmva_log.png\" alt=\"plots/vbfmva_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vbfmva_nobtag_log.png\"><img src=\"plots/vbfmva_nobtag_log.png\" alt=\"plots/vbfmva_nobtag_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/btag.png\"><img src=\"plots/btag.png\" alt=\"plots/btag.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/btag_vbf.png\"><img src=\"plots/btag_vbf.png\" alt=\"plots/btag_vbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nvertices_reweighted.png\"><img src=\"plots/nvertices_reweighted.png\" alt=\"plots/nvertices_reweighted.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nvertices_raw.png\"><img src=\"plots/nvertices_raw.png\" alt=\"plots/nvertices_raw.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  //htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt1.png\"><img src=\"plots/pt1.png\" alt=\"plots/pt1.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt2.png\"><img src=\"plots/pt2.png\" alt=\"plots/pt2.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta1.png\"><img src=\"plots/eta1.png\" alt=\"plots/eta1.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta2.png\"><img src=\"plots/eta2.png\" alt=\"plots/eta2.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "</tr>" << endl;
  //htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi1.png\"><img src=\"plots/phi1.png\" alt=\"plots/phi1.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi2.png\"><img src=\"plots/phi2.png\" alt=\"plots/phi2.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "</tr>" << endl;
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
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_0jet.png\"><img src=\"plots/svfitmass_class_0jet.png\" alt=\"plots/svfitmass_class_0jet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_0jet_lowpt.png\"><img src=\"plots/svfitmass_class_0jet_lowpt.png\" alt=\"plots/svfitmass_class_0jet_lowpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_0jet_highpt.png\"><img src=\"plots/svfitmass_class_0jet_highpt.png\" alt=\"plots/svfitmass_class_0jet_highpt.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  //htmlfile << "</tr>" << endl;
  //htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_boost.png\"><img src=\"plots/svfitmass_class_boost.png\" alt=\"plots/svfitmass_class_boost.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_boost_lowpt.png\"><img src=\"plots/svfitmass_class_boost_lowpt.png\" alt=\"plots/svfitmass_class_boost_lowpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_boost_highpt.png\"><img src=\"plots/svfitmass_class_boost_highpt.png\" alt=\"plots/svfitmass_class_boost_highpt.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_2jet.png\"><img src=\"plots/svfitmass_class_2jet.png\" alt=\"plots/svfitmass_class_2jet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_vbf.png\"><img src=\"plots/svfitmass_class_vbf.png\" alt=\"plots/svfitmass_class_vbf.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  //htmlfile << "</tr>" << endl;
  //htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_b.png\"><img src=\"plots/svfitmass_class_b.png\" alt=\"plots/svfitmass_class_b.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_b_lowptlank\" href=\"plots/svfitmass_class_b_lowpt.png\"><img src=\"plots/svfitmass_class_b_lowpt.png\" alt=\"plots/svfitmass_class_b_lowpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_b_highptlank\" href=\"plots/svfitmass_class_b_highpt.png\"><img src=\"plots/svfitmass_class_b_highpt.png\" alt=\"plots/svfitmass_class_b_highpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_0jet_log.png\"><img src=\"plots/svfitmass_class_0jet_log.png\" alt=\"plots/svfitmass_class_0jet_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_0jet_lowpt_log.png\"><img src=\"plots/svfitmass_class_0jet_lowpt_log.png\" alt=\"plots/svfitmass_class_0jet_lowpt_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_0jet_highpt_log.png\"><img src=\"plots/svfitmass_class_0jet_highpt_log.png\" alt=\"plots/svfitmass_class_0jet_highpt_log.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  //htmlfile << "</tr>" << endl;
  //htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_boost_log.png\"><img src=\"plots/svfitmass_class_boost_log.png\" alt=\"plots/svfitmass_class_boost_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_boost_lowpt_log.png\"><img src=\"plots/svfitmass_class_boost_lowpt_log.png\" alt=\"plots/svfitmass_class_boost_lowpt_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_boost_highpt_log.png\"><img src=\"plots/svfitmass_class_boost_highpt_log.png\" alt=\"plots/svfitmass_class_boost_highpt_log.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_2jet_log.png\"><img src=\"plots/svfitmass_class_2jet_log.png\" alt=\"plots/svfitmass_class_2jet_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_vbf_log.png\"><img src=\"plots/svfitmass_class_vbf_log.png\" alt=\"plots/svfitmass_class_vbf_log.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_b_log.png\"><img src=\"plots/svfitmass_class_b_log.png\" alt=\"plots/svfitmass_class_b_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_b_lowpt_log.png\"><img src=\"plots/svfitmass_class_b_lowpt_log.png\" alt=\"plots/svfitmass_class_b_lowpt_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/svfitmass_class_b_highpt_log.png\"><img src=\"plots/svfitmass_class_b_highpt_log.png\" alt=\"plots/svfitmass_class_b_highpt_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
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
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_0jet.png\"><img src=\"plots/vismass_class_0jet.png\" alt=\"plots/vismass_class_0jet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_0jet_lowpt.png\"><img src=\"plots/vismass_class_0jet_lowpt.png\" alt=\"plots/vismass_class_0jet_lowpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_0jet_highpt.png\"><img src=\"plots/vismass_class_0jet_highpt.png\" alt=\"plots/vismass_class_0jet_highpt.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_boost.png\"><img src=\"plots/vismass_class_boost.png\" alt=\"plots/vismass_class_boost.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_boost_lowpt.png\"><img src=\"plots/vismass_class_boost_lowpt.png\" alt=\"plots/vismass_class_boost_lowpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_boost_highpt.png\"><img src=\"plots/vismass_class_boost_highpt.png\" alt=\"plots/vismass_class_boost_highpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_2jet.png\"><img src=\"plots/vismass_class_2jet.png\" alt=\"plots/vismass_class_2jet.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_vbf.png\"><img src=\"plots/vismass_class_vbf.png\" alt=\"plots/vismass_class_vbf.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_b.png\"><img src=\"plots/vismass_class_b.png\" alt=\"plots/vismass_class_b.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_b_lowptlank\" href=\"plots/vismass_class_b_lowpt.png\"><img src=\"plots/vismass_class_b_lowpt.png\" alt=\"plots/vismass_class_b_lowpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_b_highptlank\" href=\"plots/vismass_class_b_highpt.png\"><img src=\"plots/vismass_class_b_highpt.png\" alt=\"plots/vismass_class_b_highpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_0jet_log.png\"><img src=\"plots/vismass_class_0jet_log.png\" alt=\"plots/vismass_class_0jet_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_0jet_lowpt_log.png\"><img src=\"plots/vismass_class_0jet_lowpt_log.png\" alt=\"plots/vismass_class_0jet_lowpt_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_0jet_highpt_log.png\"><img src=\"plots/vismass_class_0jet_highpt_log.png\" alt=\"plots/vismass_class_0jet_highpt_log.png\" width=\"100%\"></a></td>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_boost_log.png\"><img src=\"plots/vismass_class_boost_log.png\" alt=\"plots/vismass_class_boost_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_boost_lowpt_log.png\"><img src=\"plots/vismass_class_boost_lowpt_log.png\" alt=\"plots/vismass_class_boost_lowpt_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_boost_highpt_log.png\"><img src=\"plots/vismass_class_boost_highpt_log.png\" alt=\"plots/vismass_class_boost_highpt_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  //htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_2jet_log.png\"><img src=\"plots/vismass_class_2jet_log.png\" alt=\"plots/vismass_class_2jet_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_vbf_log.png\"><img src=\"plots/vismass_class_vbf_log.png\" alt=\"plots/vismass_class_vbf_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_b_log.png\"><img src=\"plots/vismass_class_b_log.png\" alt=\"plots/vismass_class_b_log.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
            
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}
