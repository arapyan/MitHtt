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
#include "MitHtt/Ntupler/interface/TPFTau.hh"
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

// scale factros
#include "MitHtt/Utils/DataMC.hh"

#endif

const Double_t pi = 3.14159265358979;

//=== MAIN MACRO =================================================================================================

void selectTauTau(const TString conf,         // input config file
		  const TString outputDir,    // output directory
		  const Double_t lumi,        // luminosity pb^-1
		  const Int_t is2012,          //2012 or 2011 data
		  const UInt_t btageff=0,     // b-tag efficiency scale factor uncertainty
		  const UInt_t jetunc=0,      // jet energy uncertainties
		  const UInt_t mistag=0,      // b mistag rate scale factor uncertainty
		  const UInt_t elescale=0     // electron energy scale/resolution uncertainty
		  ) {
  gBenchmark->Start("selectTauTau");
  
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

  const Double_t kTauPtMin = 35;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 20;
  
  Bool_t doKFactors = kTRUE;
  if(is2012)
    doKFactors = kFALSE;      // not needed in Summer12
  
  Bool_t doNpuRwgt = kTRUE;

  // Access samples and fill histograms
  TFile *infile=0;
  TTree *eventTree=0;  
  TTree *lTree=0;
 
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  mithep::TGenInfo *gen     = new mithep::TGenInfo();
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVfit");
  TClonesArray *tauArr      = new TClonesArray("mithep::TPFTau"); 
  
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
    int   lq1          = 0; outtree.Branch("q_1"        ,&lq1            ,"lq1/I"      );  //charge
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
    int   lq2          = 0; outtree.Branch("q_2"        ,&lq2            ,"lq2/I"      );  //charge
    float lIso2        = 0; outtree.Branch("iso_2"      ,&lIso2          ,"lIso2/F"    );//MVA iso for hadronic Tau, Delta Beta for muon
    float lD02         = 0; outtree.Branch("d0_2"       ,&lD02           ,"lD02/F"     );//d0 with respect to primary vertex
    float lDZ2         = 0; outtree.Branch("dZ_2"       ,&lDZ2           ,"lDZ2/F"     );//dZ with respect to primary vertex
    bool  lPassIso2    = 0; outtree.Branch("passiso_2"  ,&lPassIso2      ,"lPassIso2/O");//Whether it passes iso (not necessarily id)
    float lMt2         = 0; outtree.Branch("mt_2"       ,&lMt2           ,"lMt2/F"     );//mT of 2nd lepton wrt to MVA met
    float lMVAMt2         = 0; outtree.Branch("mtMVA_2"       ,&lMVAMt2           ,"lMVAMt2/F"     );//mT of  second lepton wrt to MVA met

    
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
  
      // which corrections to apply where
      Bool_t isdata     = !(samp->typev[ifile]==eMC);
      Bool_t isemb      = snamev[isam].Contains("emb");
      Bool_t doRecoil   = (sfname.Contains("ztt") || sfname.Contains("-zll") || sfname.Contains("zjets") || snamev[isam].Contains("_sm_") || snamev[isam].Contains("_mssm_")) && !isemb;
      Bool_t reallyDoKf = doKFactors && sfname.Contains("-gf-");
      Bool_t ismadz     = sfname.Contains("-zll") || sfname.Contains("-zjets"); // madgraph z samples
      Bool_t ismadzmm   = snamev[isam].Contains("zmm") && (sfname.Contains("-zll") || sfname.Contains("-zjets")); // madgraph z samples 
      Bool_t ismssm     = sfname.Contains("-ggh-") || sfname.Contains("-bbh-");
      Bool_t doIdScale  = !isdata;
      Bool_t doTrigScale= !isdata;
      Bool_t getGen     = !isdata; //doRecoil || reallyDoKf || ismadz ||isemb || ismssm;
      Bool_t doJetUnc   = (jetunc!=kNo);

      // PU reweighting
      TString pileupReweightFile;
      if(!is2012) {
	cout << "Fall11 sample!" << endl;
	pileupReweightFile = "$CMSSW_BASE/src/MitHtt/data/pileup/PUWeights_Fall11toFull2011_PixelLumi_50bins.root";
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
      TH1D *hKFactors = (reallyDoKf) ? kfFHPInit(higgsmass(sfname)) : 0;

      // Get the TTree
      eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
      lTree =  (TTree*)infile->Get("hEvents"); assert(lTree);

      // Set branch address to structures that will store the info  
      eventTree->SetBranchAddress("Info",  &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
      eventTree->SetBranchAddress("HPSTau", &tauArr);   TBranch *tauBr = eventTree->GetBranch("HPSTau");
      eventTree->SetBranchAddress("PFJet",    &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");      
      eventTree->SetBranchAddress("PV",       &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
      eventTree->SetBranchAddress("SVfitTauTau", &svfitArr);    TBranch *svfitBr    = eventTree->GetBranch("SVfitTauTau");
      TBranch *genBr=0;
      if(getGen) {
        eventTree->SetBranchAddress("Gen", &gen);
        genBr = eventTree->GetBranch("Gen");
      }

      // get weights for MC
      Double_t weight=1,treeEntries=-1; // (weight is only initialized for each *file*)
      if(!isdata) {
	treeEntries = (Double_t)lTree->GetEntries();
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
	if(isdata && !(info->triggerBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30])) continue;

        // good primary vertex
        if(!info->hasGoodPV) continue;
	pvArr->Clear();
	pvBr->GetEntry(ientry);	

        // loop through HPSTaus
        
	tauArr->Clear();
	tauBr->GetEntry(ientry);

	vector<const mithep::TPFTau*> goodHPSTaus;
	const mithep::TPFTau *leadTau = NULL;
	const mithep::TPFTau *subTau  = NULL;
	for(Int_t i = 0; i < tauArr->GetEntries(); i++)
	  {
	    const mithep::TPFTau *tau = dynamic_cast<mithep::TPFTau *>(tauArr->At(i));
	    assert(tau);
          
	    // Tau ID
	    if(!(passtautauId(tau,0))) continue;
	      
		// Tau Kinematics
	    if(!(tau->pt > kTauPtMin && fabs(tau->eta) < 2.1)) continue;
	       
	    //Tau HLT
	    Bool_t trigmatch = tau->hltMatchBits[kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj];
	    if(isdata && !trigmatch)     continue;
	    
	    // Tau Isolation
	    //if(!(tau->ringIso > 0.795)) continue;
		     
	    if(!(tau->hcalOverP + tau->ecalOverP > 0.2 ||
	       tau->nSignalPFChargedHadrCands > 1 ||
		 tau->nSignalPFGammaCands > 0)) continue;
	    
	    goodHPSTaus.push_back(tau);
	    if(!leadTau || tau->pt > leadTau->pt)
	      {
		subTau = leadTau;
		leadTau = tau;
	      } else if(!subTau || tau->pt > subTau->pt) {
	      subTau = tau;
	    }
          }	
          
	if(goodHPSTaus.size()<2) continue;
	if(!(passtautauId(leadTau,1))) continue;
	if(!(tauIdElectronMVA(subTau,subTau->antiEleID))) continue;


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
          if(toolbox::deltaR(leadTau->eta,leadTau->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi()) < 0.01           ) id = 1;
          if(toolbox::deltaR(subTau->eta,subTau->phi,svfit->daughter1.Eta(),svfit->daughter1.Phi())   < 0.01 && id == 0) id = 2;
          if(id == 0) continue;
          if(toolbox::deltaR(leadTau->eta,leadTau->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi()) < 0.01 && id == 2) id = 3;
          if(toolbox::deltaR(subTau->eta,subTau->phi,svfit->daughter2.Eta(),svfit->daughter2.Phi())   < 0.01 && id == 1) id = 4;
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

	// lepton 4-vectors
        TLorentzVector lep1, lep2, dilep;
	lep1.SetPtEtaPhiM(leadTau->pt, leadTau->eta, leadTau->phi, leadTau->m);
	lep2.SetPtEtaPhiM(subTau->pt, subTau->eta, subTau->phi,subTau->m);
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
	l1.SetPtEtaPhi(leadTau->pt,0,leadTau->phi);
	l2.SetPtEtaPhi(subTau->pt,0,subTau->phi);
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
	
	// trigger scale factor for MC
	Double_t trigscale = 1;

	// embedding weight for embedded sample
	Double_t embWgt = 1;
        Double_t pt1=0, eta1=0, phi1=0, pt2=0, eta2=0, phi2=0;
	if(!isdata) {
	      pt1 = gen->pt_1_a;
	      pt2 = gen->pt_2_a;
	      eta1 = gen->eta_1_a;
	      eta2 = gen->eta_2_a;
	      phi1 = gen->phi_1_a;
	      phi2 = gen->phi_2_a;
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
	lPt1		 = lep1.Pt();
	lPhi1		 = lep1.Phi();
	lEta1		 = lep1.Eta();
	lM1		 = lep1.M();
	lq1              = leadTau->q;
	lIso1		 = leadTau->ringIso;
	lD01		 = leadTau->leadChargedHadronPFCand.d0;
	lDZ1		 = leadTau->leadChargedHadronPFCand.dz;
	lPassIso1	 = (leadTau->ringIso > 0.795);
	lMt1		 = sqrt(2.0*(lep1.Pt()*met*(1.0-cos(toolbox::deltaPhi(lep1.Phi(),metphi)))));
	lMVAMt1		 = sqrt(2.0*(lep1.Pt()*mvamet*(1.0-cos(toolbox::deltaPhi(lep1.Phi(),mvametphi)))));
        lPt2		 = lep2.Pt();
        lPhi2		 = lep2.Phi();
        lEta2		 = lep2.Eta();
        lM2		 = lep2.M();
	lq2              = subTau->q;
        lIso2		 = subTau->ringIso;
        lD02		 = subTau->leadChargedHadronPFCand.d0;
        lDZ2		 = subTau->leadChargedHadronPFCand.dz;
        lPassIso2	 = (subTau->ringIso > 0.795);
        lMt2		 = sqrt(2.0*(lep2.Pt()*met*(1.0-cos(toolbox::deltaPhi(lep2.Phi(),metphi)))));
	lMVAMt2		 = sqrt(2.0*(lep2.Pt()*mvamet*(1.0-cos(toolbox::deltaPhi(lep2.Phi(),mvametphi)))));
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
	lMVA		 = 0;
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
      infile=0, eventTree=0, lTree = 0;    
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
  delete jetArr;
  delete pvArr;
  delete tauArr;
  delete svfitArr;

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  if(!doNpuRwgt) cout << endl << endl << "Not doing npv reweight!" << endl;
  cout << endl; 
  
  gBenchmark->Show("selectTauTau");
}

