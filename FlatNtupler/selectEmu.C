#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TSystem.h>                // interface to OS
#include <TFile.h>                  // file handle class
#include <TH1.h>                    // histogram base class
#include <TH2.h>                    // histogram base class
#include <TNtuple.h>                // class to access ntuples
#include <TTree.h>                  // class to access trees
#include <TRegexp.h>                // ROOT regexp class
#include <TClonesArray.h>           // ROOT array class
#include <TBenchmark.h>             // class to track macro running statistics
#include <TLorentzVector.h>         // 4-vector class
#include <TVector3.h>               // 3-vector class
#include <TMath.h>                  // ROOT math functions
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "MitHtt/Common/MitStyleRemix.hh"  // style settings for drawing
#include "MitHtt/Common/CSample.hh"        // helper class for organizing input ntuple files
#include "MitHtt/Common/MyTools.hh"        // miscellaneous helper functions

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh" 
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"   
#include "MitHtt/Ntupler/interface/TVertex.hh"   
#include "MitHtt/Ntupler/interface/TSVfit.h"
#include "MitHtt/Ntupler/interface/TSVfitter.h"

// recoil corrections
#include "MitHtt/Utils/RecoilCorrector.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh"

// event-based MVA
//#include "MitHtt/Utils/HttMVA.hh" 

#endif

//=== FUNCTION DECLARATIONS ======================================================================================
const Double_t pi = 3.14159265358979;
TRandom1 randm(0xDEADBEEF);
enum { kNo, kDown, kUp };                     // systematic variations 

// Initialize k-factors
TH1D* kfFHPInit(Int_t mH);

// Get k-factor
Double_t kfFHPValue(Double_t pt, TH1D* hKF);

// Is jet b-tagged?
Bool_t isbtagged(mithep::TJet *jet, Int_t isdata, UInt_t btageff, UInt_t mistag);

// Get higgs mass point from sample name
Int_t higgsmass(TString basename)
{
  stringstream ss(basename(TRegexp("[0-9][0-9]*"),3).Data());
  Int_t mass;
  ss >> mass;
  if(basename.Contains("-gf-")) assert(mass>85 && mass<1200);
  return mass;
}

// Get unfolding weights for embedded
Double_t embUnfoldWgt(Double_t pt1, Double_t eta1, Double_t pt2, Double_t eta2);

// Lepton id scale factors
Double_t eleIDscale(Double_t elept, Double_t eleeta);
Double_t muIDscale(Double_t mupt, Double_t mueta);

// Trigger scale factors/efficiencies
Double_t eleTrigScale(Double_t elept, Double_t eleeta);
Double_t muTrigScale(Double_t mupt, Double_t mueta);
Double_t eleTrigEff(Double_t elept, Double_t eleeta);
Double_t muTrigEff(Double_t mupt, Double_t mueta);

// Get number of entries in unskimmed tree
Double_t unskimmedEntries(TString skimname);

//=== MAIN MACRO =================================================================================================

void selectEmu(const TString conf,         // input config file
               const TString outputDir,    // output directory
	       const Double_t lumi,        // luminosity pb^-1
               const UInt_t btageff=0,     // b-tag efficiency scale factor uncertainty
               const UInt_t jetunc=0,      // jet energy uncertainties
               const UInt_t mistag=0,      // b mistag rate scale factor uncertainty
	       const UInt_t elescale=0     // electron energy scale/resolution uncertainty
) {
  gBenchmark->Start("selectEmu");
  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //============================================================================================================== 

  vector<TString>  snamev;      // sample name (for output file)  
  vector<CSample*> samplev;     // data/MC samples
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
      // fakes come from a separate macro
      if((TString(line).Contains("fake")) && (state>0)) continue;

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
      if(TString(fname).Contains("dummy",TString::kIgnoreCase)) continue;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(0);
      samplev.back()->xsecv.push_back(xsec);
    }
  }
  ifs.close();

  const TString ntupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(ntupDir,kTRUE);

  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };  // dataset type  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };      // final state type
  
  const Double_t kMuonPt1Min = 20;
  const Double_t kMuonPt2Min = 10;
  
  const Double_t kElePt1Min  = 20;
  const Double_t kElePt2Min  = 10;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 20;
  
  Bool_t doKFactors = kFALSE;      // not needed in Summer12

  Bool_t doNpuRwgt = kTRUE;

  // Access samples and fill histograms
  TFile *infile=0;
  TTree *eventTree=0;  

  //vbf MVA
  //HttMVA *vbfMVA = new HttMVA();
  //vbfMVA->Initialize("BDTG method", "/home/vdutta/cms/cmssw/new/CMSSW_4_4_1/src/MitHtt/Emu/Selection/MVAVBF_v3/weights/TMVA_BDTG.weights.xml", HttMVA::kVBF2);   // vbf mva
       
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVfit");

  Bool_t hasData = (samplev[0]->fnamev.size()>0);

  // loop over samples
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if(isam==0 && !hasData) continue;
  
    CSample* samp = samplev[isam];

    Double_t nSelEvents=0;

    // Set up output ntuple file for the sample
    TString outfname = ntupDir + TString("/") + snamev[isam] + TString("_select.root");
    TFile outfile(outfname,"RECREATE");
    TTree outtree("Events","Events");

    //Bookeeping
    int   lRun         = 0; outtree.Branch("run"        ,&lRun           ,"lRun/I"     );//Run
    int   lLumi        = 0; outtree.Branch("lumi"       ,&lLumi          ,"lLumi/I"    );//Lumi
    int   lEvt         = 0; outtree.Branch("evt"        ,&lEvt           ,"lEvt/I"     );//Evt

    //Event Variables
    int   lNPV         = 0; outtree.Branch("npv"        ,&lNPV           ,"lNPV/I"     );//NPV
    int   lNPU         = 0; outtree.Branch("npu"        ,&lNPU           ,"lNPU/I"     );//NPU
    float lRho         = 0; outtree.Branch("rho"        ,&lRho           ,"lRho/F"     );//Rho
    
    //Event Weights
    float lMCWeight    = 0; outtree.Branch("mcweight"   ,&lMCWeight      ,"lMCWeight/F");//MC Weight (xs/nevents * additional wieght (ie pt weight for gghiggs))
    float lPUWeight    = 0; outtree.Branch("puweight"   ,&lPUWeight      ,"lPUWeight/F");//Pielup Weight
    float lEffWeight   = 0; outtree.Branch("effweight"  ,&lEffWeight     ,"lEffWeight/F");//Effieiency Scale factor (all components multiplied in)
    float lWeight      = 0; outtree.Branch("weight"     ,&lWeight        ,"lWeight/F"  );//mcweight*puweight*effweight
   
    //Mass variables
    float lMass        = 0; outtree.Branch("mass"       ,&lMass          ,"lMass/F"     );//SV Fit using integration method
    float lMassUp      = 0; outtree.Branch("mass_Up"    ,&lMassUp        ,"lMassUp/F"   );//High Energy scale shape
    float lMassDown    = 0; outtree.Branch("mass_Down"  ,&lMassDown      ,"lMassDown/F" );//Low Energy Scale Shape
    float lMVis        = 0; outtree.Branch("m_vis"      ,&lMVis          ,"lMVis/F"    );//visible mass
 
    ///First lepton :  muon for mu Tau, electron for e Tau, electron for e mu, Leading (in pT) Tau for Tau Tau
    float lPt1         = 0; outtree.Branch("pt_1"       ,&lPt1           ,"lPt1/F"     ); //pT 
    float lPhi1        = 0; outtree.Branch("phi_1"      ,&lPhi1          ,"lPhi1/F"    ); //Phi 
    float lEta1        = 0; outtree.Branch("eta_1"      ,&lEta1          ,"lEta1/F"    ); //Eta 
    float lM1          = 0; outtree.Branch("m_1"        ,&lM1            ,"lM1/F"      ); //Mass 
    float lIso1        = 0; outtree.Branch("iso_1"      ,&lIso1          ,"lIso1/F"    ); //Delta Beta iso value 
    float lD01         = 0; outtree.Branch("d0_1"       ,&lD01           ,"lD01/F"     );//d0 with respect to primary vertex
    float lDZ1         = 0; outtree.Branch("dZ_1"       ,&lDZ1           ,"lDZ1/F"     );//dZ with respect to primary vertex
    bool  lPassIso1    = 0; outtree.Branch("passiso_1"  ,&lPassIso1      ,"lPassIso1/O");//Whether it passes iso 
    float lMt1         = 0; outtree.Branch("mt_1"       ,&lMt1           ,"lMt1/F"     );//mT of  first lepton wrt to met
    float lMVAMt1         = 0; outtree.Branch("mtMVA_1"       ,&lMVAMt1           ,"lMVAMt1/F"     );//mT of  first lepton wrt to MVA met

    ///Second lepton :  hadronic Tau for mu Tau had for e Tau, Muon for e mu, Trailing (in pT)  Tau for Tau Tau
    float lPt2         = 0; outtree.Branch("pt_2"       ,&lPt2           ,"lPt2/F"     );//pT
    float lPhi2        = 0; outtree.Branch("phi_2"      ,&lPhi2          ,"lPhi2/F"    );//Phi
    float lEta2        = 0; outtree.Branch("eta_2"      ,&lEta2          ,"lEta2/F"    );//Eta
    float lM2          = 0; outtree.Branch("m_2"        ,&lM2            ,"lM2/F"      );//Mass (visible mass for hadronic Tau)
    float lIso2        = 0; outtree.Branch("iso_2"      ,&lIso2          ,"lIso2/F"    );//MVA iso for hadronic Tau, Delta Beta for muon
    float lD02         = 0; outtree.Branch("d0_2"       ,&lD02           ,"lD02/F"     );//d0 with respect to primary vertex
    float lDZ2         = 0; outtree.Branch("dZ_2"       ,&lDZ2           ,"lDZ2/F"     );//dZ with respect to primary vertex
    bool  lPassIso2    = 0; outtree.Branch("passiso_2"  ,&lPassIso2      ,"lPassIso2/O");//Whether it passes iso (not necessarily id)
    float lMt2         = 0; outtree.Branch("mt_2"       ,&lMt2           ,"lMt2/F"     );//mT of 2nd lepton wrt to MVA met
    float lMVAMt2         = 0; outtree.Branch("mtMVA_2"       ,&lMVAMt2           ,"lMVAMt2/F"     );//mT of  first lepton wrt to MVA met
    
    //Met related variables
    float lMet         = 0; outtree.Branch("met"        ,&lMet           ,"lMet/F"      ); //pfmet
    float lMetPhi      = 0; outtree.Branch("metphi"     ,&lMetPhi        ,"lMetPhi/F"   ); //pfmet Phi
    float lMVAMet      = 0; outtree.Branch("mvamet"     ,&lMVAMet        ,"lMVAMet/F"   ); //mvamet
    float lMVAMetPhi   = 0; outtree.Branch("mvametphi"  ,&lMVAMetPhi     ,"lMVAMetPhi/F"); //mvamet Phi
    float lPZetaVis    = 0; outtree.Branch("pzetavis"   ,&lPZetaVis      ,"lPZetaVis/F" ); //pZeta Visible
    float lPZetaMiss   = 0; outtree.Branch("pzetamiss"  ,&lPZetaMiss     ,"lPZetaMiss/F"); //pZeta Missing
    float lPZetaMVAMiss = 0; outtree.Branch("pzetamvamiss"  ,&lPZetaMVAMiss     ,"lPZetaMVAMiss/F"); //pZeta MVA Missing
    //MET covariance matrices
    float lMetCov00    = 0; outtree.Branch("metcov00"   ,&lMetCov00      ,"lMetCov00/F"); //pf met covariance matrix 00 
    float lMetCov01    = 0; outtree.Branch("metcov01"   ,&lMetCov01      ,"lMetCov01/F"); //pf met covariance matrix 01 
    float lMetCov10    = 0; outtree.Branch("metcov10"   ,&lMetCov10      ,"lMetCov10/F"); //pf met covariance matrix 10 
    float lMetCov11    = 0; outtree.Branch("metcov11"   ,&lMetCov11      ,"lMetCov11/F"); //pf met covariance matrix 11 
    //MVAMet covariance matrices
    float lMVACov00    = 0; outtree.Branch("mvacov00"   ,&lMVACov00      ,"lMVACov00/F"); //mva met covariance matrix 00 
    float lMVACov01    = 0; outtree.Branch("mvacov01"   ,&lMVACov01      ,"lMVACov01/F"); //mva met covariance matrix 01 
    float lMVACov10    = 0; outtree.Branch("mvacov10"   ,&lMVACov10      ,"lMVACov10/F"); //mva met covariance matrix 10 
    float lMVACov11    = 0; outtree.Branch("mvacov11"   ,&lMVACov11      ,"lMVACov11/F"); //mva met covariance matrix 11 
 
    //First Jet   : leading jet after applying Jet energy corrections (excluding hadronic Tau)
    float lJPt1       = 0; outtree.Branch("jpt_1"      ,&lJPt1          ,"lJPt1/F"     );//Jet Pt after corrections
    float lJEta1      = 0; outtree.Branch("jeta_1"     ,&lJEta1         ,"lJEta1/F"    );//Jet Eta
    float lJPhi1      = 0; outtree.Branch("jphi_1"     ,&lJPhi1         ,"lJPhi1/F"    );//Jet Phi     
    float lJPtUnc1    = 0; outtree.Branch("jptunc_1"   ,&lJPtUnc1       ,"lJPtUnc1/F"  );//Jet Unc (relative to Jet corrected pT)
    float lJMVA1      = 0; outtree.Branch("jmva_1"     ,&lJMVA1         ,"lJMVA1/F"    );//Jet MVA id value
    bool  lJPass1     = 0; outtree.Branch("jpass_1"    ,&lJPass1        ,"lJPass1/O"   );//Whether Jet pass PU Id Loose WP

    //Second Jet  : 2nd leading jet (in pt) afer applying Jet energy corrections (excluding Tau)
    float lJPt2       = 0; outtree.Branch("jpt_2"      ,&lJPt2          ,"lJPt2/F"     );//Jet Pt after corrections
    float lJEta2      = 0; outtree.Branch("jeta_2"     ,&lJEta2         ,"lJEta2/F"    );//Jet Eta
    float lJPhi2      = 0; outtree.Branch("jphi_2"     ,&lJPhi2         ,"lJPhi2/F"    );//Jet Phi
    float lJPtUnc2    = 0; outtree.Branch("jptunc_2"   ,&lJPtUnc2       ,"lJPtUnc2/F"  );//Jet Unc (relative to Jet corrected pT)
    float lJMVA2      = 0; outtree.Branch("jmva_2"     ,&lJMVA2         ,"lJMVA2/F"    );//Jet MVA id value
    bool  lJPass2     = 0; outtree.Branch("jpass_2"    ,&lJPass2        ,"lJPass2/O"   );//Whether jet passes PU Id Loose WP 
    
    //B Tagged Jet : leading btagged jet (in pt) passing btag wp (pt > 20 + cvs medium)
    float lBTagPt     = 0; outtree.Branch("bpt"        ,&lBTagPt        ,"lBTagPt/F"   );//Corrected BTag Pt
    float lBTagEta    = 0; outtree.Branch("beta"       ,&lBTagEta       ,"lBTagEta/F"  );//Btag Eta
    float lBTagPhi    = 0; outtree.Branch("bphi"       ,&lBTagPhi       ,"lBTagPhi/F"  );//Btag Phi
  
    //Di Jet kinematic variables for VBF selection ==> Two leading pT Jets 
    float lMJJ        = 0; outtree.Branch("mjj"        ,&lMJJ           ,"lMJJ/F"      );//Mass Di Jet system  
    float lJDEta      = 0; outtree.Branch("jdeta"      ,&lJDEta         ,"lJDEta/F"    );//|jeta_1-jeta_2| 
    int   lNJetInGap  = 0; outtree.Branch("njetingap"  ,&lNJetInGap     ,"lNJetInGap/I");//# of Jets between two jets
    float lMVA        = 0; outtree.Branch("mva"        ,&lMVA           ,"lMVA/F"      );//VBF MVA value
    
    //Variables that go into the VBF MVA
    float lJDPhi      = 0; outtree.Branch("jdphi"      ,&lJDPhi         ,"lJDPhi/F"    );//Delta Phi between two leading jets
    float lDiJetPt    = 0; outtree.Branch("dijetpt"    ,&lDiJetPt       ,"lDiJetPt/F"  );//Pt of the di jet system
    float lDiJetPhi   = 0; outtree.Branch("dijetphi"   ,&lDiJetPhi      ,"lDiJetPhi/F" );//Phi of the di jet system
    float lHDJetPhi   = 0; outtree.Branch("hdijetphi"  ,&lHDJetPhi      ,"lHDJetPhi/F" );//Phi of the di jet system - Higgs system phi
    float lVisJetEta  = 0; outtree.Branch("visjeteta"  ,&lVisJetEta     ,"lVisJetEta/F");//TMath::Min(eta_vis - jeta,eta_vis,jeta2);
    float lPtVis      = 0; outtree.Branch("ptvis"      ,&lPtVis         ,"lPtVis/F"    );//Pt Vis
    float lPtH        = 0; outtree.Branch("pth"        ,&lPtH           ,"lPtH/F"      );//Pt of the higgs system

    //number of btags passing btag id ( pt > 20 )
    int   lNBTag      = 0; outtree.Branch("nbtag"      ,&lNBTag         ,"lNBTag/I");

    //number of jets passing jet id ( pt > 30 )
    int   lNJets      = 0; outtree.Branch("njets"      ,&lNJets         ,"lNJets/I");

    //generator lepton pt, eta, phi for embedded
    float lGenPt1     = 0; outtree.Branch("genlpt_1"   ,&lGenPt1        ,"lGenPt1/F"    );//pT leading
    float lGenPhi1    = 0; outtree.Branch("genlphi_1"  ,&lGenPhi1       ,"lGenPhi1/F"   );//Phi leading
    float lGenEta1    = 0; outtree.Branch("genleta_1"  ,&lGenEta1       ,"lGenEta1/F"   );//Eta leading
    float lGenPt2     = 0; outtree.Branch("genlpt_2"   ,&lGenPt2        ,"lGenPt2/F"    );//pT sub-leading
    float lGenPhi2    = 0; outtree.Branch("genlphi_2"  ,&lGenPhi2       ,"lGenPhi2/F"   );//Phi sub-leading
    float lGenEta2    = 0; outtree.Branch("genleta_2"  ,&lGenEta2       ,"lGenEta2/F"   );//Eta sub-leading

    UInt_t npt20jets;
    const UInt_t kMaxPt20Jets=50;
    TArrayF btagArray; btagArray.Set(kMaxPt20Jets); // array to hold b-tag values for pt-20 jets
    TArrayF jptArray; jptArray.Set(kMaxPt20Jets);   // array to hold jet pt values
    TArrayF jetaArray; jetaArray.Set(kMaxPt20Jets); // array to hold jet eta values

    // extra branches
    outtree.Branch("npt20jets",&npt20jets);
    outtree.Branch("btagArray",&btagArray);
    outtree.Branch("jptArray",&jptArray);
    outtree.Branch("jetaArray",&jetaArray);

    // loop through files
    cout <<  "processing " << snamev[isam] << ":" << endl;
    const UInt_t nfiles = samp->fnamev.size();
    for(UInt_t ifile=0; ifile<nfiles; ifile++) {
      printf("        %-55s",(samp->fnamev[ifile]+"...").Data()); fflush(stdout);
      infile = new TFile(samp->fnamev[ifile]); 
      assert(infile);

      TString sfname    = samp->fnamev[ifile];
      TString basename = sfname(sfname.Last('/')+1,sfname.Last('.') - sfname.Last('/') - 1);

      // which corrections to apply where
      Bool_t isdata     = !(samp->typev[ifile]==eMC);
      //Bool_t is52mc     = sfname.Contains("s12-");
      Bool_t isemb      = snamev[isam].Contains("emb");
      //Bool_t isfall11   = sfname.Contains("f11");
      Bool_t issamesign = snamev[isam].Contains("ss-fakes");
      Bool_t doRecoil   = (sfname.Contains("ztt") || sfname.Contains("-zll") || sfname.Contains("zjets") || snamev[isam].Contains("_sm_") || snamev[isam].Contains("_mssm_")) && !isemb;
      Bool_t reallyDoKf = doKFactors && sfname.Contains("-gf-");
      Bool_t ismadz     = sfname.Contains("-zll") || sfname.Contains("-zjets"); // madgraph z samples
      Bool_t ismadzmm   = snamev[isam].Contains("zmm") && (sfname.Contains("-zll") || sfname.Contains("-zjets")); // madgraph z samples
      //Bool_t istrainingsample = sfname.Contains("-zjets"); 
      Bool_t ismssm     = sfname.Contains("-ggh-") || sfname.Contains("-bbh-");
      Bool_t doIdScale  = !isdata;
      Bool_t doTrigScale= !isdata;
      Bool_t getGen     = doRecoil || reallyDoKf || ismadz ||isemb || ismssm;
      Bool_t doJetUnc   = (jetunc!=kNo);

      // PU reweighting
      TString pileupReweightFile;
      if(sfname.Contains("f11")) {
	cout << "Fall11 sample!" << endl;
	pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_Fall11toFull2011_PixelLumi_50bins.root"; //"/data/smurf/sixie/Pileup/weights/PileupReweighting.Fall11_To_Full2011.root";
      } else pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_S12To2012_5089ipb.root";
      TH1F *puWeights = 0;
      TFile *pufile = new TFile(pileupReweightFile.Data());
      puWeights = (TH1F*)pufile->Get("puWeights");

      // setup selecting with JSON file, if necessary
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if(isdata && (samp->jsonv[ifile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	ifstream jsonchk; jsonchk.open(samp->jsonv[ifile].Data()); assert(jsonchk.is_open()); jsonchk.close();
        rlrm.AddJSONFile(samp->jsonv[ifile].Data()); 
      }

      // recoil corrections
      RecoilCorrector *corrector=0;
      if(doRecoil) {
	cout << "doing recoil corrections" << endl;
        corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_datamm52X_v2_njet.root");
      }

      // k-factors
      TH1D *hKFactors = (reallyDoKf) ? kfFHPInit(higgsmass(basename)) : 0;

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",     &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("Muon",     &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
      eventTree->SetBranchAddress("Electron", &electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
      eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");      
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
      eventTree->SetBranchAddress("SVfitEMu", &svfitArr);    TBranch *svfitBr    = eventTree->GetBranch("SVfitEMu");
      TBranch *genBr=0;
      if(getGen) {
        eventTree->SetBranchAddress("Gen", &gen);
        genBr = eventTree->GetBranch("Gen");
      }

      // get weights for MC
      Double_t weight=1,treeEntries=-1; // (weight is only initialized for each *file*)
      if(!isdata) {
	if(sfname.Contains("_skim.root")) treeEntries = unskimmedEntries(sfname); // get entries from unskimmed file
	else                              treeEntries = (Double_t)eventTree->GetEntries();
	assert(treeEntries>0);
        weight = lumi*(samp->xsecv[ifile])/treeEntries;                           // (assumes you've merged filesets)
	if(isemb)  weight=1.0;
      }
      samp->weightv.push_back(weight);

      // counters
      Double_t nsel=0, nselvar=0; 
      Double_t nlowmass=0; // low mass z events (below 50)

      cout << eventTree->GetEntries() << " events" << endl;

      // loop over events
      for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
	if(ientry%1000000 == 0) cout << "processing " << ientry << endl;
        infoBr->GetEntry(ientry);
	//cout << info->runNum << " " << info->lumiSec << " " << info->evtNum << endl;

	if(getGen)  genBr->GetEntry(ientry);

	// skip non-tau events in madgraph sample
	if(ismadz && !ismadzmm && (fabs(gen->id_1_a)<15 || fabs(gen->id_1_a)>19)) continue;

        // skip non-mumu events in madgraph sample for zmm
        if(ismadzmm && (fabs(gen->id_1_a)>14 && fabs(gen->id_1_a)<20)) continue;

	// certified run selection
        mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);
        if(hasJSON && !rlrm.HasRunLumi(rl)) continue;

	// trigger
	if(!isemb && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] || info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL])) continue;

        // good primary vertex
        if(!info->hasGoodPV) continue;
	pvArr->Clear();
	pvBr->GetEntry(ientry);	

        // loop through muons
        vector<const mithep::TMuon*> goodMuonsv;
        vector<const mithep::TMuon*> looseMuonsv;
        muonArr->Clear();
        muonBr->GetEntry(ientry);

        for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
          const mithep::TMuon *muon = (mithep::TMuon*)((*muonArr)[i]);
	  looseMuonsv.push_back(muon);

	  // trigger matching
	  Bool_t trigmatch = ((info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && muon->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj]) || (info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && muon->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj]));
	  if(!isemb && !trigmatch)			continue;

	  if(!(muon->typeBits & kGlobal))	continue;
          if(muon->pt < kMuonPt2Min)		continue;
	  if(fabs(muon->eta) > 2.1)		continue;
          if(passTightPFMuonID(muon,0))  goodMuonsv.push_back(muon);
          //if(passTightPFMuonID(muon) && passMuonIsoPU(muon))  goodMuonsv.push_back(muon);
        }
	
        // loop through electrons 
        vector<const mithep::TElectron*> goodElectronsv;   
        electronArr->Clear();
        electronBr->GetEntry(ientry);

        for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
	  const mithep::TElectron *electron = (mithep::TElectron*)((*electronArr)[i]);

	  // trigger matching
	  Bool_t trigmatch = ((info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && electron->hltMatchBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj]) || (info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL] && electron->hltMatchBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj]));
	  if(!isemb && !trigmatch)                        continue;

          if(fabs(electron->eta) > 2.3)		continue;

	  // track matching against muons
          Bool_t hasMuonTrack=kFALSE;
          for(UInt_t imu=0; imu<goodMuonsv.size(); imu++) {
            if(electron->trkID == goodMuonsv[imu]->trkID) hasMuonTrack=kTRUE;
          }
          if(hasMuonTrack)			continue;

	  // clean against loose muons
          Bool_t matchLooseMuon=kFALSE;
	  for(UInt_t imu=0;imu<looseMuonsv.size();imu++) {
	    const mithep::TMuon *mu = looseMuonsv[imu];
	    if(toolbox::deltaR(electron->eta,electron->phi,mu->eta,mu->phi) < 0.3) matchLooseMuon=kTRUE;
	  }
	  if(matchLooseMuon)			continue;
	  
	  if(pass2012EleMVAID(electron, kMedium,0))  goodElectronsv.push_back(electron);
        }

	//----------------------------------------------------------------------------------------

	if(goodMuonsv.size()<1 || goodElectronsv.size()<1) continue;

	const mithep::TMuon *mu	   = goodMuonsv[0];
	const mithep::TElectron *ele = goodElectronsv[0];

        // SVFit
        svfitArr->Clear();
        svfitBr->GetEntry(ientry);

        Double_t met=info->pfMET, metphi=info->pfMETphi;
	Double_t mvamet =0, mvametphi = 0;
        Double_t cov_00=0, cov_01=0, cov_10=0, cov_11=0;
        Double_t mvacov_00=0, mvacov_01=0, mvacov_10=0, mvacov_11=0;
        mithep::FourVectorM dau1, dau2;

        for(Int_t i = 0; i < svfitArr->GetEntriesFast(); i++) {
          mithep::TSVfit *svfit = (mithep::TSVfit*) svfitArr->At(i);
          Int_t id = 0;
          if(toolbox::deltaR(ele->eta,ele->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.01           ) id = 1;
          if(toolbox::deltaR(mu->eta,mu->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi())   < 0.01 && id == 0) id = 2;
          if(id == 0) continue;
          if(toolbox::deltaR(ele->eta,ele->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.01 && id == 2) id = 3;
          if(toolbox::deltaR(mu->eta,mu->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi())   < 0.01 && id == 1) id = 4;
          if(id < 3) continue;
          mvamet    = svfit->mvaMET;
	  mvametphi = svfit->mvaMETphi;
          cov_00 = svfit->cov_00;
          cov_01 = svfit->cov_01;
          cov_10 = svfit->cov_10;
          cov_11 = svfit->cov_11;
          mvacov_00 = svfit->mvacov_00;
          mvacov_01 = svfit->mvacov_01;
          mvacov_10 = svfit->mvacov_10;
          mvacov_11 = svfit->mvacov_11;
          dau1 = svfit->daughter1;
          dau2 = svfit->daughter2;
        }
        if(cov_00==0 && cov_01==0 && cov_10==0 && cov_11==0) continue;

	Double_t elept = ele->pt;
	if(isemb && elescale==kNo) {
	  elept = 1.03*ele->pt;
	  met -= 0.03*ele->pt;
	}
	if(elescale==kUp) {
	  if(isemb) {
	    elept = 1.045*ele->pt;
	    met -= 0.045*ele->pt;
	  } else {
	    elept = 1.01*ele->pt;
	    met -= 0.01*ele->pt;
	  }
	}
	if(elescale==kDown) {
	  if(isemb) {
            elept = 1.015*ele->pt;
            met -= 0.015*ele->pt;
          } else {
            elept = 0.99*ele->pt;
            met += 0.01*ele->pt;
          }
        }

	if(mu->pt < kMuonPt2Min  || elept < kElePt2Min) continue;
	if(mu->pt < kMuonPt1Min  && elept < kElePt1Min) continue;

	// trigger requirements
	if(mu->pt  < kMuonPt1Min) {
          if(!isemb && !(info->triggerBits[kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL])) continue; // if failed trig1
	}
	else if(elept < kElePt1Min) {
	  if(!isemb && !(info->triggerBits[kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL])) continue; // if failed trig2
	}

	// same-sign requirements
	if(issamesign) {
	  if(mu->q != ele->q) continue;
	}
	else {
	  if(mu->q == ele->q) continue;
	} 

	// lepton 4-vectors
        TLorentzVector lep1, lep2, dilep;
	lep1.SetPtEtaPhiM(elept, ele->eta, ele->phi, 0.000511);
	lep2.SetPtEtaPhiM(mu->pt, mu->eta, mu->phi, 0.105658369);
	dilep = lep1+lep2;

        // loop through jets      
        jetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets = 0, nbjets = 0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	btagArray.Reset();	jptArray.Reset();	jetaArray.Reset();	npt20jets=0;
        for(Int_t i=0; i<jetArr->GetEntriesFast(); i++) {
	  mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);

	  if(doJetUnc) jet->pt *= (jetunc==kDown) ? (1-jet->unc) : (1+jet->unc);

          if(toolbox::deltaR(jet->eta,jet->phi,lep1.Eta(),lep1.Phi()) < 0.5) continue;
          if(toolbox::deltaR(jet->eta,jet->phi,lep2.Eta(),lep2.Phi()) < 0.5) continue;

          if(fabs(jet->eta) > 5) continue;
	  if(!jet->id) continue;

	  // look for b-jets
	  Int_t btagopt = 0;
	  if(isdata||isemb) btagopt = 1;
	  else btagopt = 2;
	  Bool_t btagged = isbtagged(jet,btagopt,btageff,mistag);
	  if((jet->pt > kBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
	    assert(npt20jets<kMaxPt20Jets);
	    btagArray.AddAt(jet->csv,npt20jets);
	    npt20jets++;
	    if(btagged) {
	      nbjets++;
	      if(!bjet || jet->pt > bjet->pt)
		bjet = jet; // leading b-jet
	    }
	  }

	  // look for jets
          if(jet->pt > kJetPtMin) {
            assert(njets<kMaxPt20Jets);
            jptArray.AddAt(jet->pt,njets);
            jetaArray.AddAt(jet->eta,njets);
  	    njets++;
	    if(!jet1 || jet->pt > jet1->pt) { // jet1 is highest pt, jet2 next highest
	      jet2 = jet1;
	      jet1 = jet;
	    } else if(!jet2 || jet->pt > jet2->pt) {
	      jet2 = jet;
	    }
          }		    
        }

	// dijet system
        TLorentzVector jv1, jv2, dijet;
	Int_t nCentralJets=0;
	if(njets>1) {
	  jv1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
	  jv2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
          dijet = jv1+jv2;
          for(Int_t i=2; i<jetArr->GetEntriesFast(); i++) {
            mithep::TJet *jet = (mithep::TJet*)((*jetArr)[i]);
	    if(!(jet->pt > kJetPtMin && fabs(jet->eta)<5 && jet->id==1)) continue;
	    if(jet1->eta > jet2->eta && jet->eta > jet2->eta && jet->eta < jet1->eta) nCentralJets++;
	    else if(jet2->eta > jet1->eta && jet->eta > jet1->eta && jet->eta < jet2->eta) nCentralJets++;
	  }
        }

	// recoil corrections
	Double_t pU1      = 0;  //--
	Double_t pU2      = 0;  //--
        if(doRecoil) corrector->CorrectAll(met, metphi, gen->vpt_a, gen->vphi_a, dilep.Pt(), dilep.Phi(), pU1, pU2, 0, 0, njets);

	// calculate projection variables
	TVector3 l1,l2,metv,mvametv;
	l1.SetPtEtaPhi(mu->pt,0,mu->phi);
	l2.SetPtEtaPhi(elept,0,ele->phi);
	metv.SetPtEtaPhi(met,0,metphi);
	mvametv.SetPtEtaPhi(mvamet,0,mvametphi);
	TVector3 bisector(l1.Unit() + l2.Unit());
	bisector = bisector.Unit();
	Double_t projVis  = (l1+l2).Dot(bisector);
	Double_t projMet  =  metv.Dot(bisector);
	Double_t projMVAMet = mvametv.Dot(bisector);

	// ditau system
	TLorentzVector met4v, higgs;
	met4v.SetPtEtaPhiM(met,0,metphi,0);
	higgs = dilep+met4v;

        // get k-factor if necessary
        Double_t kf=1;
        if(reallyDoKf) kf = kfFHPValue(gen->vpt_a, hKFactors);

	// do vertex reweighting
	Double_t npuWgt = 1;
	if(!isdata && !isemb && doNpuRwgt) {
	  assert(puWeights);
	  Int_t npuxbin = puWeights->GetXaxis()->FindFixBin(TMath::Min(double(info->nPU), 59.499));
	  npuWgt = puWeights->GetBinContent(npuxbin);
	}

	// lepton ID corrections
	Double_t idscale = 1;
	if(doIdScale) idscale = muIDscale(mu->pt,mu->eta)*eleIDscale(elept,ele->eta);

	// trigger scale factor for MC
	Double_t trigscale = 1;
	if(doTrigScale && !isemb) trigscale=muTrigScale(mu->pt,mu->eta)*eleTrigScale(elept,ele->eta);
	if(doTrigScale && isemb) trigscale=muTrigEff(mu->pt,mu->eta)*eleTrigEff(elept,ele->eta);

	// embedding weight for embedded sample
	Double_t embWgt = 1;
        Double_t pt1=0, eta1=0, phi1=0, pt2=0, eta2=0, phi2=0;
	if(doRecoil) {
	  if(gen->id_1_a == EGenType::kTauElectron && gen->id_2_a == EGenType::kTauMuon)
	    {
	      pt1 = gen->pt_1_a;
	      pt2 = gen->pt_2_a;
	      eta1 = gen->eta_1_a;
	      eta2 = gen->eta_2_a;
	      phi1 = gen->phi_1_a;
	      phi2 = gen->phi_2_a;
	    }
	  else
	    {
	      pt1 = gen->pt_2_a;
	      pt2 = gen->pt_1_a;
	      eta1 = gen->eta_2_a;
	      eta2 = gen->eta_1_a;
	      phi1 = gen->phi_2_a;
	      phi2 = gen->phi_1_a;
	    }
	}
	if(isemb)    {
	  if(gen->pt_1_a > gen->pt_2_a) {
	    pt1 = gen->pt_1_a;
	    eta1 = gen->eta_1_a;
	    phi1 = gen->phi_1_a;
            pt2 = gen->pt_2_a;
            eta2 = gen->eta_2_a;
	    phi2 = gen->phi_2_a;
	  } else {
            pt2 = gen->pt_1_a;
            eta2 = gen->eta_1_a;
            phi2 = gen->phi_1_a;
            pt1 = gen->pt_2_a;
            eta1 = gen->eta_2_a;
            phi1 = gen->phi_2_a;
	  }
	  embWgt=info->embWeight;
	  //embWgt=info->embWeight*embUnfoldWgt(pt1,eta1,pt2,eta2);
	}

	// events passing selection in this file
	nsel    += weight*kf*npuWgt*trigscale*idscale*embWgt;
	nselvar += weight*weight*kf*kf*npuWgt*npuWgt*trigscale*trigscale*idscale*idscale*embWgt*embWgt;
	if(corrector && (gen->vmass_a < 50)) nlowmass += weight*kf*npuWgt*trigscale*idscale*embWgt;

	// passing events in whole sample 
        nSelEvents += weight*kf*npuWgt*trigscale*idscale*embWgt;

	lRun		 = info->runNum;
	lLumi		 = info->lumiSec;
	lEvt		 = info->evtNum;
	lNPV		 = pvArr->GetEntriesFast();
	lNPU		 = info->nPU;
	lRho		 = info->rho;
	lMCWeight	 = weight*kf*embWgt/lumi;
	lPUWeight	 = npuWgt;
	lEffWeight	 = trigscale*idscale;
	lWeight		 = weight*kf*npuWgt*trigscale*idscale*embWgt/lumi;
	lMass		 = (ismssm)? gen->vmass_a: 0;
	lMassUp		 = 0;
	lMassDown	 = 0;
	lMVis		 = dilep.M();
	lPt1		 = elept;
	lPhi1		 = ele->phi;
	lEta1		 = ele->eta;
	lM1		 = 0.000511;
	lIso1		 = eleIsoPU(ele);
	lD01		 = ele->d0;
	lDZ1		 = ele->dz;
	lPassIso1	 = passEleIsoPU(ele);
	lMt1		 = sqrt(2.0*(elept*met*(1.0-cos(toolbox::deltaPhi(ele->phi,metphi)))));
	lMVAMt1		 = sqrt(2.0*(elept*mvamet*(1.0-cos(toolbox::deltaPhi(ele->phi,mvametphi)))));
        lPt2		 = mu->pt;
        lPhi2		 = mu->phi;
        lEta2		 = mu->eta;
        lM2		 = 0.105658369;
        lIso2		 = muonIsoPU(mu);
        lD02		 = mu->d0;
        lDZ2		 = mu->dz;
        lPassIso2	 = passMuonIsoPU(mu);
        lMt2		 = sqrt(2.0*(mu->pt*met*(1.0-cos(toolbox::deltaPhi(mu->phi,metphi)))));
	lMt2		 = sqrt(2.0*(mu->pt*mvamet*(1.0-cos(toolbox::deltaPhi(mu->phi,mvametphi)))));
	lMet		 = met;
	lMetPhi		 = metphi;
	lMVAMet		 = mvamet;
	lMVAMetPhi	 = mvametphi;
	lPZetaVis	 = projVis;
	lPZetaMiss	 = projMet;
	lPZetaMVAMiss    = projMVAMet;
	lMetCov00	 = cov_00;
        lMetCov01	 = cov_01;
        lMetCov10	 = cov_10;
        lMetCov11	 = cov_11;
        lMVACov00	 = mvacov_00;
        lMVACov01	 = mvacov_01;
        lMVACov10	 = mvacov_10;
        lMVACov11	 = mvacov_11;
	lJPt1		 = (jet1) ? jet1->pt  : 0;
        lJEta1		 = (jet1) ? jet1->eta : 0;
        lJPhi1		 = (jet1) ? jet1->phi : 0;
        lJPtUnc1	 = (jet1) ? jet1->unc : 0;
	lJMVA1		 = (jet1) ? jet1->mva : 0;
	lJPass1		 = (jet1) ? jet1->id  : 0;
        lJPt2		 = (jet2) ? jet2->pt  : 0;
        lJEta2		 = (jet2) ? jet2->eta : 0;
        lJPhi2		 = (jet2) ? jet2->phi : 0;
        lJPtUnc2	 = (jet2) ? jet2->unc : 0;
        lJMVA2		 = (jet2) ? jet2->mva : 0;
        lJPass2		 = (jet2) ? jet2->id  : 0;
	lBTagPt		 = (bjet) ? bjet->pt  : 0;
        lBTagEta	 = (bjet) ? bjet->eta : 0;
        lBTagPhi	 = (bjet) ? bjet->phi : 0;
	lMJJ		 = (njets>1) ? dijet.M() : 0;
	lJDEta		 = (njets>1) ? fabs(jet1->eta - jet2->eta) : 0;
	lNJetInGap	 = (njets>1) ? nCentralJets : 0;
	lMVA		 = 0.0;//vbfMVA->MVAValue(lMJJ, lJDEta, lJDPhi, lDiJetPt, lPtH, lHDJetPhi, lVisJetEta, lPtVis);
	lJDPhi		 = (njets>1) ? toolbox::deltaPhi(jet1->phi,jet2->phi) : 0;
	lDiJetPt	 = (njets>1) ? dijet.Pt()  : 0;
        lDiJetPhi	 = (njets>1) ? dijet.Phi() : 0;
        lHDJetPhi	 = (njets>1) ? toolbox::deltaPhi(dijet.Phi(),higgs.Phi()) : 0;
	lVisJetEta	 = (njets>1) ? TMath::Min(fabs(dilep.Eta()-jet1->eta),fabs(dilep.Eta()-jet2->eta)) : 0;
	lPtVis		 = dilep.Pt();
	lPtH		 = higgs.Pt();
	lNBTag		 = nbjets;
	lNJets		 = njets;
	lGenPt1		 = pt1;
        lGenEta1	 = eta1;
	lGenPhi1	 = phi1;
        lGenPt2          = pt2;
        lGenEta2         = eta2;
        lGenPhi2         = phi2;

	outtree.Fill();

      }

      printf("%8.2f +/- %-8.2f\n",nsel,sqrt(nselvar));
      if(nlowmass > 0) printf("           ---> selected events with z mass < 50:  %10.3f (out of %15.3f)\n",nlowmass,nsel);

      delete infile;
      if(corrector && doRecoil) {cout << "recoil corrections used" << endl; delete corrector;}
      infile=0, eventTree=0;    
    }
    outfile.Write();
    outfile.Close();

    if(samp->typev.size()>0 && samp->typev[0]==eMC)
      printf("    Yields for %1.2f/fb:",lumi/1000.);
    else
      printf("    Yields for data:    ");

    printf("%10.2f\n",nSelEvents);
    cout << endl;
  }

  delete info;
  delete gen;
  delete muonArr;
  delete electronArr;
  delete jetArr;
  delete pvArr;
  delete svfitArr;
  //delete fitter;


  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  if(!doNpuRwgt) cout << endl << endl << "Not doing npv reweight!" << endl;
  cout << endl; 
  
  gBenchmark->Show("selectEmu");
}


//=== FUNCTION DEFINITIONS ======================================================================================
//----------------------------------------------------------------------------------------
TH1D* kfFHPInit(Int_t mH)
{
  TH1D *kfhist=0; 
  char kfilename[100];
  sprintf(kfilename, "$CMSSW_BASE/src/MitHtt/Utils/HiggsKFactors/weight_ptH_%d.root", mH);
  cout << "Getting k-factors from " << kfilename << endl;
  TFile *kfile = TFile::Open(kfilename); assert(kfile->IsOpen());
  TDirectory *kfdir = (TDirectory*)kfile->FindObjectAny("powheg_weight");
  char kfhistname[100];
  sprintf(kfhistname, "weight_hqt_fehipro_fit_%d", mH);
  cout << "kfactor histogram: " << kfhistname << endl;
  kfhist = (TH1D*)(kfdir->Get(kfhistname)); assert(kfhist);
  return kfhist;
} 
//--------------------------------------------------------------------------------------------------
Double_t kfFHPValue(Double_t pt, TH1D* hKF)
{ 
  return hKF->Interpolate(pt);
}
//----------------------------------------------------------------------------------------
Bool_t isbtagged(mithep::TJet *jet, Int_t isdata, UInt_t btageff, UInt_t mistag)
{

  //          mistag                         scale factor
  // TCHEM  0.0175 \pm .0003 \pm .0038      1.21 \pm .02 \pm .17
  //          btag eff.                      scale factor
  // TCHEM  0.63 \pm 0.01                   0.93 \pm 0.02 \pm 0.07

  // new scale factors
  // TCHEM	btag eff: 0.96 \pm 0.04		mistag rate: 0.0286 \pm 0.0003		mistag scale factor: 1.20 \pm 0.14
  // CSVM	btag eff: 0.97 \pm 0.04		mistag rate: 0.0152 \pm 0.0002		mistag scale factor: 1.10 \pm 0.11

  Bool_t btagged;
  Double_t demoteProb=0; // ~probability to demote from tagged
  if(btageff==kNo)        demoteProb = fabs(1-0.97); //1-0.93;  // SF = 0.93 -> 0.07 = (prob to demote from tagged status)
  else if(btageff==kDown) demoteProb = fabs(1-0.97+0.04); //1-0.93+0.07;
  else if(btageff==kUp)   demoteProb = fabs(1-0.97-0.04); //1-0.93-0.07;
  Double_t promoteProb=0; // ~probability to promote to tagged
  if(mistag==kNo)         promoteProb = fabs(1.10-1)*0.0152/(1-0.0152); //(1.21-1)*0.0145/(1-0.0145);  // (1-SF)*mistag = (prob. to promote to tagged status)*(1-mistag)
  else if(mistag==kDown)  promoteProb = fabs(1.10-1+0.11)*0.0152/(1-0.0152);
  else if(mistag==kUp)    promoteProb = fabs(1.10-1-0.11)*0.0152/(1-0.0152);

  UInt_t jetflavor = 0;
                   
  if(isdata == 1) {
    if(jet->csv>0.679) btagged = kTRUE;
    else               btagged = kFALSE;
  } else { // MC
    //if(isdata == 0)jetflavor = abs(jet->mcFlavor);
    jetflavor = abs(jet->matchedId);
    if(jetflavor==5) {
      if(jet->csv>0.679) {
      if(randm.Uniform()>demoteProb) btagged = kTRUE;  // leave it tagged
      else                           btagged = kFALSE; // demote it
      } else                           btagged = kFALSE; // leave it untagged
    } else { // not bjet
      if(jet->csv>0.679)                   btagged = kTRUE;  // leave it tagged
      else if(randm.Uniform()<promoteProb) btagged = kTRUE;  // promote to tagged
      else                                 btagged = kFALSE; // leave it untagged
    }
  }

  return btagged;
}  
//----------------------------------------------------------------------------------------
Double_t embUnfoldWgt(Double_t pt1, Double_t eta1, Double_t pt2, Double_t eta2)
{
  TFile *unfFile1   = TFile::Open("data/unfold/v8/Unfold2D_1.root"); assert(unfFile1->IsOpen());
  TH2F  *unfWeight1 = (TH2F*) unfFile1->FindObjectAny("UnfoldDen1");
  TFile *unfFile2   = TFile::Open("data/unfold/v8/Unfold2D_2.root"); assert(unfFile2->IsOpen());
  TH2F  *unfWeight2 = (TH2F*) unfFile2->FindObjectAny("UnfoldDen2");
  double weight1 = unfWeight1->GetBinContent(unfWeight1->GetXaxis()->FindBin(eta1),unfWeight1->GetYaxis()->FindBin(pt1));
  double weight2 = unfWeight2->GetBinContent(unfWeight2->GetXaxis()->FindBin(eta2),unfWeight2->GetYaxis()->FindBin(pt2));
  unfFile1->Close();
  unfFile2->Close();
  return weight1*weight2;
}
//----------------------------------------------------------------------------------------
Double_t muIDscale(Double_t mupt, Double_t mueta)
{
  if((fabs(mueta) > 2.1) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  /*else if(mupt > 20) {
    if(fabs(mueta) < 1.479)     return 0.9930;
    else                          return 0.9981;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 1.479)     return 0.9455;
    else                          return 0.9604;
  }
  else {
    if(fabs(mueta) < 1.479)     return 0.9226;
    else                          return 0.9856;
  }*/
  /*if(mupt > 20) {
    if(fabs(mueta) < 1.479)   return 0.9918;
    else                      return 0.9942;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 1.479)   return 0.9951;
    else                      return 1.0024;
  }
  else {
    if(fabs(mueta) < 1.479)   return 0.9912;
    else                      return 1.0364;
  }*/
  if(mupt > 20) {
    if(fabs(mueta) < 1.5)   return 0.9900;
    else                      return 0.9924;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 1.5)   return 0.9875;
    else                      return 0.9891;
  }
  else {
    if(fabs(mueta) < 1.5)   return 0.9851;
    else                      return 0.9956;
  }
}
//----------------------------------------------------------------------------------------    
Double_t eleIDscale(Double_t elept, Double_t eleeta)
{
  if((fabs(eleeta) > 2.3) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  /*else if(elept > 20) {
    if(fabs(eleeta) < 1.479) return 0.9896;
    else                       return 1.0532;
  }
  else if(elept > 15) {
    if(fabs(eleeta) < 1.479) return 0.9783;
    else                       return 1.0623;
  }
  else {
    if(fabs(eleeta) < 1.479) return 1.1134;
    else                       return 1.1946;
  }*/
  /*if(elept > 20) {
    if(fabs(eleeta) < 0.8) return 0.9590;
    else if(fabs(eleeta) < 1.479) return 0.9544;
    else                     return 0.9684;
  }
  else if(elept > 15) {
    if(fabs(eleeta) < 0.8) return 0.9256;
    else if(fabs(eleeta) < 1.479) return 0.8530;
    else                     return 0.8376;
  }
  else {
    if(fabs(eleeta) < 0.8) return 0.8401;
    else if(fabs(eleeta) < 1.479) return 0.8374;
    else                     return 0.7217;
  }*/
  if(elept > 20) {
    if(fabs(eleeta) < 0.8) return 0.9562;
    else if(fabs(eleeta) < 1.479) return 0.9507;
    else                     return 0.9584;
  }
  else if(elept > 15) {
    if(fabs(eleeta) < 0.8) return 0.9030;
    else if(fabs(eleeta) < 1.479) return 0.8623;
    else                     return 0.7935;
  }
  else {
    if(fabs(eleeta) < 0.8) return 0.8500;
    else if(fabs(eleeta) < 1.479) return 0.8995;
    else                     return 0.6683;
  }
}
//----------------------------------------------------------------------------------------
Double_t eleTrigScale(Double_t elept, Double_t eleeta)
{
  if((fabs(eleeta) > 2.3) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  /*else if(elept > 30) {
    if(fabs(eleeta) < 1.479) return 1.0031;
    else                       return 1.0078;
  }
  else if(elept > 20) {
    if(fabs(eleeta) < 1.479) return 1.0012;
    else                       return 1.0040;
  }
  else if(elept > 15) {
    if(fabs(eleeta) < 1.479) return 1.0026;
    else                       return 1.0504;
  }
  else {
    if(fabs(eleeta) < 1.479) return 0.9769;
    else                       return 0.9696;
  }*/
  // MVA ID + DB Iso
  /*else if(elept > 30) {
    if(fabs(eleeta) < 0.8)        return 0.9992;
    else if(fabs(eleeta) < 1.479) return 0.9892;
    else                          return 0.9688;
  } 
  else if(elept > 25) {
    if(fabs(eleeta) < 0.8)        return 0.9697;
    else if(fabs(eleeta) < 1.479) return 0.9840;
    else                          return 1.0100;
  } 
  else if(elept > 20) {
    if(fabs(eleeta) < 0.8)        return 0.9773;
    else if(fabs(eleeta) < 1.479) return 0.9642;
    else                          return 1.0148;
  } 
  else if(elept > 15) {
    if(fabs(eleeta) < 0.8)        return 0.9915;
    else if(fabs(eleeta) < 1.479) return 0.9955;
    else                          return 1.0693;
  } 
  else {
    if(fabs(eleeta) < 0.8)        return 0.9910;
    else if(fabs(eleeta) < 1.479) return 0.8252;
    else                          return 0.9612;
  }*/
  else if(elept > 30) {
    if(fabs(eleeta) < 0.8)        return 1.0047;
    else if(fabs(eleeta) < 1.479) return 0.9981;
    else                          return 0.9829;
  } 
  else if(elept > 25) {
    if(fabs(eleeta) < 0.8)        return 1.0009;
    else if(fabs(eleeta) < 1.479) return 1.0202;
    else                          return 0.9920;
  } 
  else if(elept > 20) {
    if(fabs(eleeta) < 0.8)        return 0.9777;
    else if(fabs(eleeta) < 1.479) return 0.9646;
    else                          return 0.9663;
  } 
  else if(elept > 15) {
    if(fabs(eleeta) < 0.8)        return 1.0116;
    else if(fabs(eleeta) < 1.479) return 0.9866;
    else                          return 0.9742;
  } 
  else {
    if(fabs(eleeta) < 0.8)        return 0.9824;
    else if(fabs(eleeta) < 1.479) return 0.8795;
    else                          return 0.8638;
  }
}
//----------------------------------------------------------------------------------------
Double_t muTrigScale(Double_t mupt, Double_t mueta)
{
  if((fabs(mueta) > 2.1) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  /*else if(mupt > 30) {
    if(fabs(mueta) < 1.479) return 0.9922;
    else                      return 1.0550;
  }
  else if(mupt > 20) {
    if(fabs(mueta) < 1.479) return 0.9936;
    else                      return 1.0358;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 1.479) return 0.9918;
    else                      return 1.0712;
  }
  else {
    if(fabs(mueta) < 1.479)  return 1.0052;
    else                       return 1.0277;
  }*/
  //TPFID + DB iso
  /*else if(mupt > 30) {
    if(fabs(mueta) < 0.8)        return 1.0666;
    else if(fabs(mueta) < 1.2)   return 1.1177;
    else                         return 1.1247;
  }
  else if(mupt > 25) {
    if(fabs(mueta) < 0.8)        return 0.9958;
    else if(fabs(mueta) < 1.2)   return 1.0635;
    else                         return 1.0185;
  }
  else if(mupt > 20) {
    if(fabs(mueta) < 0.8)        return 1.0087;
    else if(fabs(mueta) < 1.2)   return 0.9794;
    else                         return 0.9610;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 0.8)        return 0.9977;
    else if(fabs(mueta) < 1.2)   return 1.0389;
    else                         return 1.0218;
  } 
  else {
    if(fabs(mueta) < 0.8)        return 1.0032;
    else if(fabs(mueta) < 1.2)   return 0.9939;
    else                         return 0.9828;
  }*/
  else if(mupt > 30) {
    if(fabs(mueta) < 0.8)        return 1.0670;
    else if(fabs(mueta) < 1.2)   return 1.0926;
    else                         return 1.1180;
  }
  else if(mupt > 25) {
    if(fabs(mueta) < 0.8)        return 0.9933;
    else if(fabs(mueta) < 1.2)   return 1.0224;
    else                         return 0.9755;
  }
  else if(mupt > 20) {
    if(fabs(mueta) < 0.8)        return 1.0014;
    else if(fabs(mueta) < 1.2)   return 0.9588;
    else                         return 0.9883;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 0.8)        return 0.9969;
    else if(fabs(mueta) < 1.2)   return 1.0325;
    else                         return 1.0063;
  } 
  else {
    if(fabs(mueta) < 0.8)        return 0.9928;
    else if(fabs(mueta) < 1.2)   return 1.0088;
    else                         return 0.9743;
  }
}
//----------------------------------------------------------------------------------------
Double_t muTrigEff(Double_t mupt, Double_t mueta)
{
  if((fabs(mueta) > 2.1) || (mupt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  /*else if(mupt > 30) {
    if(fabs(mueta) < 0.8)        return 0.9789;
    else if(fabs(mueta) < 1.2)   return 0.9278;
    else                         return 0.9215;
  }
  else if(mupt > 25) {
    if(fabs(mueta) < 0.8)        return 0.9607;
    else if(fabs(mueta) < 1.2)   return 0.9623;
    else                         return 0.9085;
  }
  else if(mupt > 20) {
    if(fabs(mueta) < 0.8)        return 0.9848;
    else if(fabs(mueta) < 1.2)   return 0.9231;
    else                         return 0.9023;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 0.8)        return 0.9738;
    else if(fabs(mueta) < 1.2)   return 0.9286;
    else                         return 0.9476;
  } 
  else {
    if(fabs(mueta) < 0.8)        return 0.9794;
    else if(fabs(mueta) < 1.2)   return 0.9609;
    else                         return 0.9363;
  }*/
  else if(mupt > 30) {
    if(fabs(mueta) < 0.8)        return 0.9726;
    else if(fabs(mueta) < 1.2)   return 0.9111;
    else                         return 0.9075;
  }
  else if(mupt > 25) {
    if(fabs(mueta) < 0.8)        return 0.9621;
    else if(fabs(mueta) < 1.2)   return 0.9489;
    else                         return 0.9035;
  }
  else if(mupt > 20) {
    if(fabs(mueta) < 0.8)        return 0.9794;
    else if(fabs(mueta) < 1.2)   return 0.9161;
    else                         return 0.9291;
  }
  else if(mupt > 15) {
    if(fabs(mueta) < 0.8)        return 0.9777;
    else if(fabs(mueta) < 1.2)   return 0.9319;
    else                         return 0.9459;
  } 
  else {
    if(fabs(mueta) < 0.8)        return 0.9716;
    else if(fabs(mueta) < 1.2)   return 0.9458;
    else                         return 0.9100;
  }
}
//----------------------------------------------------------------------------------------
Double_t eleTrigEff(Double_t elept, Double_t eleeta)
{
  if((fabs(eleeta) > 2.3) || (elept < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  /*else if(elept > 30) {
    if(fabs(eleeta) < 0.8)        return 0.9992;
    else if(fabs(eleeta) < 1.479) return 0.9892;
    else                          return 0.9688;
  } 
  else if(elept > 25) {
    if(fabs(eleeta) < 0.8)        return 0.9697;
    else if(fabs(eleeta) < 1.479) return 0.9840;
    else                          return 1.0100;
  } 
  else if(elept > 20) {
    if(fabs(eleeta) < 0.8)        return 0.9209;
    else if(fabs(eleeta) < 1.479) return 0.9580;
    else                          return 0.9464;
  } 
  else if(elept > 15) {
    if(fabs(eleeta) < 0.8)        return 0.9279;
    else if(fabs(eleeta) < 1.479) return 0.9568;
    else                          return 0.9565;
  } 
  else {
    if(fabs(eleeta) < 0.8)        return 0.9618;
    else if(fabs(eleeta) < 1.479) return 0.9813;
    else                          return 0.9677;
  }*/
  else if(elept > 30) {
    if(fabs(eleeta) < 0.8)        return 0.9661;
    else if(fabs(eleeta) < 1.479) return 0.9795;
    else                          return 0.9751;
  } 
  else if(elept > 25) {
    if(fabs(eleeta) < 0.8)        return 0.9468;
    else if(fabs(eleeta) < 1.479) return 0.9796;
    else                          return 0.9550;
  } 
  else if(elept > 20) {
    if(fabs(eleeta) < 0.8)        return 0.9213;
    else if(fabs(eleeta) < 1.479) return 0.9496;
    else                          return 0.9228;
  } 
  else if(elept > 15) {
    if(fabs(eleeta) < 0.8)        return 0.9023;
    else if(fabs(eleeta) < 1.479) return 0.9299;
    else                          return 0.8869;
  } 
  else {
    if(fabs(eleeta) < 0.8)        return 0.7925;
    else if(fabs(eleeta) < 1.479) return 0.7472;
    else                          return 0.7378;
  }
}
//----------------------------------------------------------------------------------------
Double_t unskimmedEntries(TString skimname)
{
  Double_t entries;
  
  skimname.ReplaceAll("_emu_skim.root","_ntuple.root");
  skimname.ReplaceAll("_emunod0_skim.root","_ntuple.root");
  TFile unskimmed(skimname);
  assert(unskimmed.IsOpen());
  TTree *tree = 0;
  unskimmed.GetObject("Events",tree);
  assert(tree);
  entries = (Double_t)tree->GetEntries();
  unskimmed.Close();

  return entries;
}
