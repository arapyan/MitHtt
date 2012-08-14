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
#include "MitHtt/Ntupler/interface/TSVfitter.h"
#include "EmuLimitInputs.hh"               // helper class for organizing inputs for limit calculation
#include "EmuPlots.hh"                     // helper class for producing plots
#include "EmuHTML.hh"                      // helper class for producing html page with plots

#endif

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

  const TString plotDir = outputDir + TString("/plots");
  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };  // final state type

  const Double_t kssFakeWgt = 1.3; //1.11;  // ratio of fake-rate fakes to same-sign fakes

  //
  // scale factors for each category (used only for madgraph ztt)
  //
  const Double_t kCat_0jet = 1.00;
  const Double_t kCat_boost = 1.00;
  const Double_t kCat_b = 1.00;
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
  vector<TH1F*> hProjMetv, hProjVisv, hProjVarv,hRawProjVarv;
  vector<TH1F*> hNjetsv, hNbjetsv;
  vector<TH1F*> hBdiscrv;
  vector<TH1F*> hDPhiv;				// delta phi between leptons
  vector<TH1F*> hPt1v, hEta1v, hPhi1v;   	// leading lepton
  vector<TH1F*> hPt2v, hEta2v, hPhi2v;   	// trailing lepton
  vector<TH1F*> hPtMuv, hEtaMuv, hPhiMuv;  	// muon
  vector<TH1F*> hPtElev, hEtaElev, hPhiElev; 	// electron
  vector<TH1F*> hMjjv, hDEtav, hEtaProdv, hJetDPhiv;       // kinematics of first two jets
  vector<TH1F*> hPtjjv, hPtHv, hHDiJetDPhiv;    // inputs for VBF mva
  vector<TH1F*> hBJetPtv, hBJetEtav, hBJetPhiv; // leading b-jet (can be same as leading two jets)
  vector<TH1F*> hNPVv, hNPVrawv;                // primary vertexes
  vector<TH1F*> hVBFJetPt1v, hVBFJetEta1v;      // leading jet
  vector<TH1F*> hVBFJetPt2v, hVBFJetEta2v;      // second pt jet
  vector<TH1F*> hBoostJetPtv, hBoostJetEtav;    // leading jet
  vector<TH1F*> hBDTv;                          // VBF bdtg values
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

  // Class vbf: two or more pt 30 jets and vbf cuts
  vector<TH1F*> hMass_vbfv;
  vector<TH1F*> hMassVis_vbfv;
  vector<TH1F*> hMassL_vbfv;

  vector<Double_t> nSelv,       nSelVarv;
  vector<Double_t> nSel_iv,     nSelVar_iv;
  vector<Double_t> nSel_0jet_lowptv,   nSelVar_0jet_lowptv;
  vector<Double_t> nSel_0jet_highptv,  nSelVar_0jet_highptv;
  vector<Double_t> nSel_boost_lowptv,  nSelVar_boost_lowptv;
  vector<Double_t> nSel_boost_highptv, nSelVar_boost_highptv;
  vector<Double_t> nSel_b_lowptv,      nSelVar_b_lowptv;
  vector<Double_t> nSel_b_highptv,     nSelVar_b_highptv;
  vector<Double_t> nSel_vbfv,   nSelVar_vbfv;
  
  vector<Double_t> nUnskEventsv;


  char hname[100];
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    // after lepton selection
    sprintf(hname,"hMet_%i",isam);     hMetv.push_back(new TH1F(hname,"",30,0,150));        hMetv[isam]->Sumw2();
    sprintf(hname,"hMetRaw_%i",isam);  hMetRawv.push_back(new TH1F(hname,"",30,0,150));     hMetRawv[isam]->Sumw2();
    sprintf(hname,"hProjMet_%i",isam);       hProjMetv.push_back(new TH1F(hname,"",30,-100,120));     hProjMetv[isam]->Sumw2();
    sprintf(hname,"hProjVis_%i",isam);       hProjVisv.push_back(new TH1F(hname,"",20,0,50));         hProjVisv[isam]->Sumw2();
    sprintf(hname,"hProjVar_%i",isam);       hProjVarv.push_back(new TH1F(hname,"",25,-200,100));     hProjVarv[isam]->Sumw2();
    sprintf(hname,"hRawProjVar_%i",isam);    hRawProjVarv.push_back(new TH1F(hname,"",25,-200,100));  hRawProjVarv[isam]->Sumw2();
    sprintf(hname,"hNjets_%i",isam);         hNjetsv.push_back(new TH1F(hname,"",11,-0.5,10.5));      hNjetsv[isam]->Sumw2();
    sprintf(hname,"hNbjets_%i",isam);        hNbjetsv.push_back(new TH1F(hname,"",6,-0.5,5.5));       hNbjetsv[isam]->Sumw2();
    sprintf(hname,"hBdiscr_%i",isam);        hBdiscrv.push_back(new TH1F(hname,"",20,0,1.0));          hBdiscrv[isam]->Sumw2();
    sprintf(hname,"hDPhi_%i",isam);    hDPhiv.push_back(new TH1F(hname,"",36,0,180));       hDPhiv[isam]->Sumw2();
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
    sprintf(hname,"hJetDPhi_%i",isam); hJetDPhiv.push_back(new TH1F(hname,"",15,0,180));    hJetDPhiv[isam]->Sumw2();
    sprintf(hname,"hMjj_%i",isam);     hMjjv.push_back(new TH1F(hname,"",6,0,1200));       hMjjv[isam]->Sumw2();
    sprintf(hname,"hPtjj_%i",isam);    hPtjjv.push_back(new TH1F(hname,"",20,0,500));       hPtjjv[isam]->Sumw2();
    sprintf(hname,"hPtH_%i",isam);     hPtHv.push_back(new TH1F(hname,"",20,0,500));        hPtHv[isam]->Sumw2();
    sprintf(hname,"hHDiJetDPhi_%i",isam); hHDiJetDPhiv.push_back(new TH1F(hname,"",15,0,180));    hHDiJetDPhiv[isam]->Sumw2();
    sprintf(hname,"hDEta_%i",isam);    hDEtav.push_back(new TH1F(hname,"",10,0,10));         hDEtav[isam]->Sumw2();
    sprintf(hname,"hEtaProd_%i",isam); hEtaProdv.push_back(new TH1F(hname,"",30,-7.5,7.5)); hEtaProdv[isam]->Sumw2();
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

    // lepton selection
    sprintf(hname,"hMass_%i",isam);          hMassv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));                                  hMassv[isam]->Sumw2();
    sprintf(hname,"hMassVis_%i",isam);       hMassVisv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));                               hMassVisv[isam]->Sumw2();
    sprintf(hname,"hMassHigh_%i",isam);      hMassHighv.push_back(new TH1F(hname,"",100,0,2000));                            hMassHighv[isam]->Sumw2();
    sprintf(hname,"hMassVisHigh_%i",isam);   hMassVisHighv.push_back(new TH1F(hname,"",100,0,1000));                         hMassVisHighv[isam]->Sumw2();
    sprintf(hname,"hMassL_%i",isam);         hMassLv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));           hMassLv[isam]->Sumw2();
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
    // VBF
    sprintf(hname,"hMass_vbf_%i",isam);    hMass_vbfv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));      hMass_vbfv[isam]->Sumw2();
    sprintf(hname,"hMassVis_vbf_%i",isam); hMassVis_vbfv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges));   hMassVis_vbfv[isam]->Sumw2();
    if(domssm) {sprintf(hname,"hMassL_vbf_%i",isam);   hMassL_vbfv.push_back(new TH1F(hname,"",massLNbins_MSSM-1,massLEdges_MSSM)); hMassL_vbfv[isam]->Sumw2();}
    else {sprintf(hname,"hMassL_vbf_%i",isam);   hMassL_vbfv.push_back(new TH1F(hname,"",massLNbins_SM-1,massLEdges)); hMassL_vbfv[isam]->Sumw2();}

    nSelv.push_back(0);       nSelVarv.push_back(0);
    nSel_iv.push_back(0);     nSelVar_iv.push_back(0);
    nSel_0jet_lowptv.push_back(0);    nSelVar_0jet_lowptv.push_back(0);
    nSel_0jet_highptv.push_back(0);   nSelVar_0jet_highptv.push_back(0);
    nSel_boost_lowptv.push_back(0);   nSelVar_boost_lowptv.push_back(0);
    nSel_boost_highptv.push_back(0);  nSelVar_boost_highptv.push_back(0);
    nSel_b_lowptv.push_back(0);       nSelVar_b_lowptv.push_back(0);
    nSel_b_highptv.push_back(0);      nSelVar_b_highptv.push_back(0);
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
    float lMSV         = 0;    eventTree->SetBranchAddress("m_svi"       ,&lMSV           );//SV Fit using integration method
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
    float lVBFMVA      = 0;    eventTree->SetBranchAddress("vbfmva"     ,&lVBFMVA        );//VBFMVA value
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
        if(lMJJ>0 && (0.85*lPZetaVis - lPZetaMiss <= 25) && lNBTag==0) {
	  if(isam>0 || (isam==0 && vbfMvaValue < 0.5)) {
	    hBDTv[isam]          ->Fill(vbfMvaValue,  wgt);
	  }
	}
      }


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

      //TLorentzVector svfit = fitter->fit(covM, dau1, dau2, lMet, lMetPhi);
      //lMSV = svfit.M();

      // fill plots after lepton selection
      if(isam!=ifake) {
	hMetv[isam]        ->Fill(lMVAMet,    	 wgt);
	hMetRawv[isam]     ->Fill(lMet,      	 wgt);
	hProjMetv[isam]    ->Fill(lPZetaMiss,     	 wgt);
	hProjVisv[isam]    ->Fill(lPZetaVis,     	 wgt);
	hProjVarv[isam]    ->Fill(-0.85*lPZetaVis + lPZetaMiss,  wgt);
	hNjetsv[isam]      ->Fill(lNJets,      wgt);
	hNbjetsv[isam]     ->Fill(lNBTag,     wgt);

	assert(npt20jets<kMaxPt20Jets);
	for(UInt_t ib=0;ib<npt20jets;ib++)
	  hBdiscrv[isam]   ->Fill((*btagArray)[ib],  wgt);

        hDPhiv[isam]   ->Fill(toolbox::deltaPhi(lPhi1,lPhi2)*180./pi, wgt);
	hPt1v[isam]    ->Fill((lPt1 > lPt2 ? lPt1 : lPt2),    wgt);
	hEta1v[isam]   ->Fill((lPt1 > lPt2 ? lEta1 : lEta2),   wgt);
	hPhi1v[isam]   ->Fill((lPt1 > lPt2 ? lPhi1 : lPhi2),   wgt);
	hPt2v[isam]    ->Fill((lPt1 > lPt2 ? lPt2 : lPt1),    wgt);
	hEta2v[isam]   ->Fill((lPt1 > lPt2 ? lEta2 : lEta1),   wgt);
	hPhi2v[isam]   ->Fill((lPt1 > lPt2 ? lPhi2 : lPhi1),   wgt);
	hPtMuv[isam]   ->Fill(lPt2,         wgt);
	hEtaMuv[isam]  ->Fill(lEta2,        wgt);
	hPhiMuv[isam]  ->Fill(lPhi2,        wgt);
	hPtElev[isam]  ->Fill(lPt1,         wgt);
	hEtaElev[isam] ->Fill(lEta1,        wgt);
	hPhiElev[isam] ->Fill(lPhi1,        wgt);
	if(lNJets>0 && (0.85*lPZetaVis - lPZetaMiss <= 25)) {
	  hBoostJetPtv[isam]	->Fill(lJPt1,     wgt);
	  hBoostJetEtav[isam]	->Fill(lJEta1,    wgt);
	}
	if(lNJets>1 && (0.85*lPZetaVis - lPZetaMiss <= 25)) {
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

	nSelv[isam]    += wgt;
	nSelVarv[isam] += wgt*wgt;
      }

      Bool_t bjets = lNBTag>0;
      Bool_t vbfcuts = lNJets>=2 && (lMJJ>mjjMin) && (fabs(lJEta1-lJEta2)>dEtaMin) && lNJetInGap==0;
      Bool_t vbfmvacuts = lNJets>=2 && lNJetInGap==0 && vbfMvaValue>0.5;
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

  // write hists to input file for limit inputs and pre-fit plots
  TString flimitstr(plotDir + "/limit-inputs.root");
  makeLimitInputs(domssm, flimitstr, ecorr, snamev, hMassL_0jet_lowptv, hMassL_0jet_highptv, hMassL_boost_lowptv, hMassL_boost_highptv, hMassL_b_lowptv, hMassL_b_highptv, hMassL_vbfv);

  
  // fix problem with fake template in VBF category
  double fakebin = hMass_vbfv[ifake]->GetBinContent(5);
  double fakebinerr = hMass_vbfv[ifake]->GetBinError(5);
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


  // make plots
  if(ecorr==kCenter) {
    makePlot(plotDir,format,"nvertices_raw","","N_{PV}","Events",snamev,samplev,hNPVrawv);
    makePlot(plotDir,format,"nvertices_reweighted","","N_{PV}","Events",snamev,samplev,hNPVv);
    makePlot(plotDir,format,"met","","#slash{E}_{T} [GeV]","Events",snamev,samplev,hMetv);
    makePlot(plotDir,format,"pzetamiss","","#slash{p}_{#zeta} [GeV]","Events",snamev,samplev,hProjMetv);
    makePlot(plotDir,format,"pzetavis","","p_{#zeta}^{vis} [GeV]","Events",snamev,samplev,hProjVisv);
    makePlot(plotDir,format,"pzetavar","","#slash{p}_{#zeta} - 0.85 #times p_{#zeta}^{vis} [GeV]","Events",snamev,samplev,hProjVarv);
    makePlot(plotDir,format,"njets","","Number of Jets","Events",snamev,samplev,hNjetsv);
    makePlot(plotDir,format,"njets_log","","Number of Jets","Events",snamev,samplev,hNjetsv,0,1);
    makePlot(plotDir,format,"nbjets","","Number of b-Tagged Jets","Events",snamev,samplev,hNbjetsv);
    makePlot(plotDir,format,"nbjets_log","","Number of b-Tagged Jets","Events",snamev,samplev,hNbjetsv,0,1,0.2,100.0);
    makePlot(plotDir,format,"btag","","b-Tag Discriminator","Events",snamev,samplev,hBdiscrv);
    makePlot(plotDir,format,"dphi_emu","","#Delta^{}#phi_{e#mu} [deg]","Events",snamev,samplev,hDPhiv);
    makePlot(plotDir,format,"pt_mu","","#mu p_{T} [GeV]","Events",snamev,samplev,hPtMuv);
    makePlot(plotDir,format,"eta_mu","","#mu #eta","Events",snamev,samplev,hEtaMuv);
    makePlot(plotDir,format,"phi_mu","","#mu #phi","Events",snamev,samplev,hPhiMuv);
    makePlot(plotDir,format,"pt_e","","e p_{T} [GeV]","Events",snamev,samplev,hPtElev);
    makePlot(plotDir,format,"eta_e","","e #eta","Events",snamev,samplev,hEtaElev);
    makePlot(plotDir,format,"phi_e","","e #phi","Events",snamev,samplev,hPhiElev);
    makePlot(plotDir,format,"pt1","","leading lepton p_{T} [GeV]","Events",snamev,samplev,hPt1v);
    makePlot(plotDir,format,"eta1","","leading lepton #eta","Events",snamev,samplev,hEta1v);
    makePlot(plotDir,format,"phi1","","leading lepton #phi","Events",snamev,samplev,hPhi1v);
    makePlot(plotDir,format,"pt2","","trailing lepton p_{T} [GeV]","Events",snamev,samplev,hPt2v);
    makePlot(plotDir,format,"eta2","","trailing lepton #eta","Events",snamev,samplev,hEta2v);
    makePlot(plotDir,format,"phi2","","trailing lepton #phi","Events",snamev,samplev,hPhi2v);
    makePlot(plotDir,format,"bjetpt","","lead b-jet p_{T} [GeV]","Events",snamev,samplev,hBJetPtv);
    makePlot(plotDir,format,"bjeteta","","lead b-jet #eta","Events",snamev,samplev,hBJetEtav);
    makePlot(plotDir,format,"bjetphi","","lead b-jet #phi","Events",snamev,samplev,hBJetPhiv);
    makePlot(plotDir,format,"boostjetpt","","Leading jet p_{T} [GeV]","Events",snamev,samplev,hBoostJetPtv);
    makePlot(plotDir,format,"boostjeteta","","Leading jet #eta","Events",snamev,samplev,hBoostJetEtav);
    makePlot(plotDir,format,"vbfjet1pt","","Leading jet p_{T} [GeV]","Events",snamev,samplev,hVBFJetPt1v);
    makePlot(plotDir,format,"vbfjet1eta","","Leading jet #eta","Events",snamev,samplev,hVBFJetEta1v);
    makePlot(plotDir,format,"vbfjet2pt","","Second jet p_{T} [GeV]","Events",snamev,samplev,hVBFJetPt2v);
    makePlot(plotDir,format,"vbfjet2eta","","Second jet #eta","Events",snamev,samplev,hVBFJetEta2v);
    makePlot(plotDir,format,"mjj","","M(jj) [GeV]","Events",snamev,samplev,hMjjv,1);
    makePlot(plotDir,format,"ptjj","","p_{T}(jj) [GeV]","Events",snamev,samplev,hPtjjv,1);
    makePlot(plotDir,format,"ditaupt","","p_{T}(#tau#tau) [GeV]","Events",snamev,samplev,hPtHv,1);
    makePlot(plotDir,format,"ditaujjdphi","","#Delta#phi (#tau#tau,jj)","Events",snamev,samplev,hHDiJetDPhiv,1);
    makePlot(plotDir,format,"jetdeta","","#Delta#eta (jj)","Events",snamev,samplev,hDEtav,1);
    makePlot(plotDir,format,"jetdphi","","#Delta#phi (jj)","Events",snamev,samplev,hJetDPhiv,1);
    makePlot(plotDir,format,"jetetaprod","","(#eta1*#eta2) (jj)","Events",snamev,samplev,hEtaProdv,1);
    makePlot(plotDir,format,"vbfmva","","VBF MVA value","Events",snamev,samplev,hBDTv,1);
    makePlot(plotDir,format,"vbfmva_log","","VBF MVA value","Events",snamev,samplev,hBDTv,1,1);
    makePlot(plotDir,format,"svfitmass_lept","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMassv);
    makePlot(plotDir,format,"svfitmass_lept_high","","m_{#tau#tau} [GeV]","Events",snamev,samplev,hMassHighv,0,1);
    makePlot(plotDir,format,"vismass_lept","","m_{vis} [gev]","dn/dm_{vis}",snamev,samplev,hMassVisv);
    makePlot(plotDir,format,"vismass_lept_high","","m_{vis} [GeV]","Events",snamev,samplev,hMassVisHighv,0,1);
    makePlot(plotDir,format,"svfitmass_incl","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_iv);
    makePlot(plotDir,format,"svfitmass_incl_log","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_iv,0,1);
    makePlot(plotDir,format,"vismass_incl","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_iv);
    makePlot(plotDir,format,"vismass_incl_log","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_iv,0,1);
    makePlot(plotDir,format,"svfitmass_0jet_lowpt","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_0jet_lowptv);
    makePlot(plotDir,format,"svfitmass_0jet_lowpt_log","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_0jet_lowptv,0,1);
    makePlot(plotDir,format,"vismass_0jet_lowpt","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_0jet_lowptv);
    makePlot(plotDir,format,"vismass_0jet_lowpt_log","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_0jet_lowptv,0,1);
    makePlot(plotDir,format,"svfitmass_0jet_highpt","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_0jet_highptv);
    makePlot(plotDir,format,"svfitmass_0jet_highpt_log","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_0jet_highptv,0,1);
    makePlot(plotDir,format,"vismass_0jet_highpt","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_0jet_highptv);
    makePlot(plotDir,format,"vismass_0jet_highpt_log","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_0jet_highptv,0,1);
    makePlot(plotDir,format,"svfitmass_boost_lowpt","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_boost_lowptv);
    makePlot(plotDir,format,"svfitmass_boost_lowpt_log","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_boost_lowptv,0,1);
    makePlot(plotDir,format,"vismass_boost_lowpt","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_boost_lowptv);
    makePlot(plotDir,format,"vismass_boost_lowpt_log","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_boost_lowptv,0,1);
    makePlot(plotDir,format,"svfitmass_boost_highpt","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_boost_highptv);
    makePlot(plotDir,format,"svfitmass_boost_highpt_log","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_boost_highptv,0,1);
    makePlot(plotDir,format,"vismass_boost_highpt","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_boost_highptv);
    makePlot(plotDir,format,"vismass_boost_highpt_log","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_boost_highptv,0,1);
    makePlot(plotDir,format,"svfitmass_b_lowpt","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_b_lowptv);
    makePlot(plotDir,format,"svfitmass_b_lowpt_log","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_b_lowptv,0,1);
    makePlot(plotDir,format,"vismass_b_lowpt","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_b_lowptv);
    makePlot(plotDir,format,"vismass_b_lowpt_log","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_b_lowptv,0,1);
    makePlot(plotDir,format,"svfitmass_b_highpt","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_b_highptv);
    makePlot(plotDir,format,"svfitmass_b_highpt_log","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_b_highptv,0,1);
    makePlot(plotDir,format,"vismass_b_highpt","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_b_highptv);
    makePlot(plotDir,format,"vismass_b_highpt_log","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_b_highptv,0,1);
    makePlot(plotDir,format,"svfitmass_vbf","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_vbfv,1);
    makePlot(plotDir,format,"svfitmass_vbf_log","","m_{#tau#tau} [GeV]","dN/dm_{#tau#tau}",snamev,samplev,hMass_vbfv,1,1);
    makePlot(plotDir,format,"vismass_vbf","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_vbfv,1);
    makePlot(plotDir,format,"vismass_vbf_log","","m_{vis} [GeV]","dN/dm_{vis}",snamev,samplev,hMassVis_vbfv,1,1);

    // generate html page with all the plots
    makeHTML(outputDir);
  }


  // Summary print out
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  Double_t nTotal=0,       nTotalVar=0;    // background MC totals
  Double_t nTotal_i=0,     nTotalVar_i=0;
  Double_t nTotal_0jet_lowpt=0, nTotalVar_0jet_lowpt=0;
  Double_t nTotal_0jet_highpt=0, nTotalVar_0jet_highpt=0;
  Double_t nTotal_boost_lowpt=0, nTotalVar_boost_lowpt=0;
  Double_t nTotal_boost_highpt=0, nTotalVar_boost_highpt=0;
  Double_t nTotal_b_lowpt=0, nTotalVar_b_lowpt=0;
  Double_t nTotal_b_highpt=0, nTotalVar_b_highpt=0;
  Double_t nTotal_vbf=0,   nTotalVar_vbf=0;

  if(ecorr==kCenter) {
    ofstream yieldfile;
    yieldfile.open(plotDir+"/yields.txt");
    
    yieldfile << setw(33) << "lepton sele." << setw(20) << "inclusive" << setw(24) << "0-jet, low pT" << setw(22) << "0-jet, high pT" << setw(24) << "1-jet, low pT" << setw(22) << "1-jet, high pT" << setw(25) << "btag, low pT" << setw(23) << "btag, high pT" << setw(18) << "VBF" << endl;
  
    if(hMass_iv[iewk] && hMass_iv[izmm]) { // add zmm to ewk contribution
      nSelv[iewk]               += nSelv[izmm];                nSelVarv[iewk]             += nSelVarv[izmm];
      nSel_iv[iewk]             += nSel_iv[izmm];              nSelVar_iv[iewk]           += nSelVar_iv[izmm];
      nSel_0jet_lowptv[iewk]    += nSel_0jet_lowptv[izmm];     nSelVar_0jet_lowptv[iewk]  += nSelVar_0jet_lowptv[izmm];
      nSel_0jet_highptv[iewk]   += nSel_0jet_highptv[izmm];    nSelVar_0jet_highptv[iewk] += nSelVar_0jet_highptv[izmm];
      nSel_boost_lowptv[iewk]   += nSel_boost_lowptv[izmm];    nSelVar_boost_lowptv[iewk] += nSelVar_boost_lowptv[izmm];
      nSel_boost_highptv[iewk]  += nSel_boost_highptv[izmm];   nSelVar_boost_highptv[iewk]+= nSelVar_boost_highptv[izmm];
      nSel_b_lowptv[iewk]       += nSel_b_lowptv[izmm];        nSelVar_b_lowptv[iewk]     += nSelVar_b_lowptv[izmm];
      nSel_b_highptv[iewk]      += nSel_b_highptv[izmm];       nSelVar_b_highptv[iewk]    += nSelVar_b_highptv[izmm];
      nSel_vbfv[iewk]             += nSel_vbfv[izmm];            nSelVar_vbfv[iewk]           += nSelVar_vbfv[izmm];
    }
  
    yieldfile << setw(15) << "fakes";
    yieldfile << setw(10) << setprecision(3) << fixed << nSelv[issfake]             << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVarv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[issfake]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_iv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_lowptv[issfake]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_highptv[issfake] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_lowptv[issfake] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_highptv[issfake]<< " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_lowptv[issfake]     << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_lowptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_highptv[issfake]    << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_b_highptv[issfake]);
    yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[ifake]           << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_vbfv[ifake]);
    yieldfile << endl;
  
  
    nTotal                 += nSelv[issfake] + nSelv[iewk] + nSelv[ittbar];
    nTotalVar              += nSelVarv[issfake] + nSelVarv[iewk] + nSelVarv[ittbar];
    nTotal_i               += nSel_iv[issfake] + nSel_iv[iewk] + nSel_iv[ittbar];
    nTotalVar_i            += nSelVar_iv[issfake] + nSelVar_iv[iewk] + nSelVar_iv[ittbar];
    nTotal_0jet_lowpt      += nSel_0jet_lowptv[issfake] + nSel_0jet_lowptv[iewk_7TeV] + nSel_0jet_lowptv[ittbar];
    nTotalVar_0jet_lowpt   += nSelVar_0jet_lowptv[issfake] + nSelVar_0jet_lowptv[iewk_7TeV] + nSelVar_0jet_lowptv[ittbar];
    nTotal_0jet_highpt     += nSel_0jet_highptv[issfake] + nSel_0jet_highptv[iewk_7TeV] + nSel_0jet_highptv[ittbar];
    nTotalVar_0jet_highpt  += nSelVar_0jet_highptv[issfake] + nSelVar_0jet_highptv[iewk_7TeV] + nSelVar_0jet_highptv[ittbar];
    nTotal_boost_lowpt     += nSel_boost_lowptv[issfake] + nSel_boost_lowptv[iewk] + nSel_boost_lowptv[ittbar];
    nTotalVar_boost_lowpt  += nSelVar_boost_lowptv[issfake] + nSelVar_boost_lowptv[iewk] + nSelVar_boost_lowptv[ittbar];
    nTotal_boost_highpt    += nSel_boost_highptv[issfake] + nSel_boost_highptv[iewk_7TeV] + nSel_boost_highptv[ittbar];
    nTotalVar_boost_highpt += nSelVar_boost_highptv[issfake] + nSelVar_boost_highptv[iewk_7TeV] + nSelVar_boost_highptv[ittbar];
    nTotal_b_lowpt         += nSel_b_lowptv[issfake] + nSel_b_lowptv[iewk] + nSel_b_lowptv[ittbar];
    nTotalVar_b_lowpt      += nSelVar_b_lowptv[issfake] + nSelVar_b_lowptv[iewk] + nSelVar_b_lowptv[ittbar];
    nTotal_b_highpt        += nSel_b_highptv[issfake] + nSel_b_highptv[iewk] + nSel_b_highptv[ittbar];
    nTotalVar_b_highpt     += nSelVar_b_highptv[issfake] + nSelVar_b_highptv[iewk] + nSelVar_b_highptv[ittbar];
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
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_lowptv[isam]  << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_lowptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_highptv[isam] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_0jet_highptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_lowptv[isam] << " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_lowptv[isam]);
          yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_highptv[isam]<< " +/- " << setw(6) << setprecision(3) << fixed << sqrt(nSelVar_boost_highptv[isam]);
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
        nTotal_0jet_lowpt       += nSel_0jet_lowptv[isam];
        nTotalVar_0jet_lowpt    += nSelVar_0jet_lowptv[isam];
        nTotal_0jet_highpt      += nSel_0jet_highptv[isam];
        nTotalVar_0jet_highpt   += nSelVar_0jet_highptv[isam];
        nTotal_boost_lowpt      += nSel_boost_lowptv[isam];
        nTotalVar_boost_lowpt   += nSelVar_boost_lowptv[isam];
        nTotal_boost_highpt     += nSel_boost_highptv[isam];
        nTotalVar_boost_highpt  += nSelVar_boost_highptv[isam];
        nTotal_b_lowpt          += nSel_b_lowptv[isam];
        nTotalVar_b_lowpt       += nSelVar_b_lowptv[isam];
        nTotal_b_highpt         += nSel_b_highptv[isam];
        nTotalVar_b_highpt      += nSelVar_b_highptv[isam];
        nTotal_vbf      	+= nSel_vbfv[isam];
        nTotalVar_vbf   	+= nSelVar_vbfv[isam];
      }
      yieldfile << endl;
      yieldfile << setw(15) << "bkg MC";
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal               << " +/- " << sqrt(nTotalVar);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_i             << " +/- " << sqrt(nTotalVar_i);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_0jet_lowpt    << " +/- " << sqrt(nTotalVar_0jet_lowpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_0jet_highpt   << " +/- " << sqrt(nTotalVar_0jet_highpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_boost_lowpt   << " +/- " << sqrt(nTotalVar_boost_lowpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_boost_highpt  << " +/- " << sqrt(nTotalVar_boost_highpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_b_lowpt       << " +/- " << sqrt(nTotalVar_b_lowpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_b_highpt      << " +/- " << sqrt(nTotalVar_b_highpt);
      yieldfile << setw(10) << setprecision(3) << fixed << nTotal_vbf           << " +/- " << sqrt(nTotalVar_vbf);
      yieldfile << endl;
    }
  
    if(hasData) {
      yieldfile << setw(15) << "Data";
      yieldfile << setw(10) << setprecision(3) << fixed << nSelv[0]             << " +/- " << sqrt(nSelVarv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_iv[0]           << " +/- " << sqrt(nSelVar_iv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_lowptv[0]  << " +/- " << sqrt(nSelVar_0jet_lowptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_0jet_highptv[0] << " +/- " << sqrt(nSelVar_0jet_highptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_lowptv[0] << " +/- " << sqrt(nSelVar_boost_lowptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_boost_highptv[0]<< " +/- " << sqrt(nSelVar_boost_highptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_lowptv[0]     << " +/- " << sqrt(nSelVar_b_lowptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_b_highptv[0]    << " +/- " << sqrt(nSelVar_b_highptv[0]);
      yieldfile << setw(10) << setprecision(3) << fixed << nSel_vbfv[0]         << " +/- " << sqrt(nSelVar_vbfv[0]);
      yieldfile << endl;
    }
  
    // cat the yields to stdout
    TString cmd("cat "+plotDir+"/yields.txt");
    system(cmd.Data());
  
    // write out the acceptances
    yieldfile << endl << endl;
    yieldfile << setw(25) << " " << setw(15) << "initial events";
    yieldfile << setw(20) << "lepton sele." <<setw(17) << "inclusive" << setw(19) << "0-jet, low pT" << setw(19) << "0-jet, high pT" << setw(19) << "1-jet, low pT" << setw(19) << "1-jet, high pT" << setw(19) << "btag, low pT" << setw(19) << "btag, high pT" << setw(22) << "vbf" << endl;
  
    for(UInt_t isam=0; isam<samplev.size(); isam++) {
      if(!(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH"))) continue;
      if(snamev[isam].Contains("htt_")) continue;
      yieldfile << setw(25) << snamev[isam] << "  " <<
        setw(15) << nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSelv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_iv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_0jet_lowptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_0jet_highptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_boost_lowptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_boost_highptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_b_lowptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_b_highptv[isam]/nUnskEventsv[isam] <<
        setw(20) << setprecision(6) <<nSel_vbfv[isam]/nUnskEventsv[isam] << endl;
    }
    yieldfile.close();
  
  }
  
  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  cout << endl;
  
  gBenchmark->Show("plotEmu");      
}

