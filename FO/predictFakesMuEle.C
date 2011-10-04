#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                   // access to gROOT, entry point to ROOT system
#include <TSystem.h>                 // interface to OS
#include <TFile.h>                   // file handle class
#include <TTree.h>                   // class to access ntuples
#include <TClonesArray.h>            // ROOT array class
#include <TBenchmark.h>              // class to track macro running statistics
#include <TLorentzVector.h>          // class for 4-vector calculations
#include <TH2D.h>                    // 2D histogram class
#include <vector>                    // vector class
#include <iostream>                  // standard I/O
#include <iomanip>                   // functions to format standard I/O
#include <fstream>                   // functions for file I/O

#include "Common/MyTools.hh"         // custom helper functions
#include "Common/CEffUser2D.hh"

// define structures to read in ntuple
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Ntupler/interface/TSVFit.hh"

// lumi section selection with JSON files
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"

// lepton ID helper functions
#include "MitHtt/Utils/LeptonIDCuts.hh" 
#include "MitHtt/Emu/Selection/EmuData.hh"
#endif


//=== MAIN MACRO =================================================================================================

void predictFakesMuEle(const TString frname,     // fake rate data file
                       const TString outputDir   // output directory
) {
  gBenchmark->Start("predictFakesMuEle");


  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================   
  
  enum { eMC, eMuEl, eDiMu, eMu, eDiEl, eEl };
  
  vector<TString> fnamev;
  vector<TString> datatypev;
  vector<TString> jsonv;

  fnamev.push_back("/scratch/vdutta/htt/r11a-mueg-m10-v1_emu_skim.root"); datatypev.push_back(eMuEl);
  jsonv.push_back("/home/vdutta/cms/root/json/Cert_160404-163869_7TeV_May10ReReco_Collisions11_JSON_v2.txt");

  fnamev.push_back("/scratch/vdutta/htt/r11a-mueg-pr-v4_emu_skim.root");  datatypev.push_back(eMuEl);
  jsonv.push_back("/home/vdutta/cms/root/json/Cert_160404-167913_7TeV_PromptReco_Collisions11_JSON.txt");

  //fnamev.push_back("/scratch/vdutta/htt/lp/r11a-mueg-a05-v1_ntuple.root");  datatypev.push_back(eMuEl);
  //jsonv.push_back("/home/vdutta/cms/root/json/Cert_160404-172802_7TeV_PromptReco_Collisions11_JSON_v2.txt");

  //fnamev.push_back("/scratch/vdutta/htt/lp/r11a-mueg-pr-v6_ntuple.root");  datatypev.push_back(eMuEl);
  //jsonv.push_back("/home/vdutta/cms/root/json/Cert_160404-172802_7TeV_PromptReco_Collisions11_JSON_v2.txt");

  const Double_t kMuonPt1Min = 20;
  const Double_t kMuonPt2Min = 10;
  
  const Double_t kElePt1Min  = 20;
  const Double_t kElePt2Min  = 10;

  const Double_t kJetPtMin   = 30;
  const Double_t kBJetPtMin  = 20;
  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };      // final state type
    
  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //==============================================================================================================

  //
  // parse fake rate data
  //
  TFile frfile(frname);
  TH2D *hFR     = (TH2D*)frfile.Get("frEtaPt");
  TH2D *hFRErrl = (TH2D*)frfile.Get("errlEtaPt");
  TH2D *hFRErrh = (TH2D*)frfile.Get("errhEtaPt");
  CEffUser2D fr;
  fr.loadEff(hFR,hFRErrl,hFRErrh);
  
  // cout << fr.getEff(0.3,60) << endl; return;
    
  vector<Double_t> nEventsv;     // number of events per sample
  vector<Double_t> nEventsVarv;  // error^2 on number of events per sample
  Double_t nDenomTL  = 0;        // number of tight+loose candidates
  Double_t nTight    = 0;        // number of tight objects
  
  // electron + fake yield prediction
  Double_t nFOTL=0, nFOTLVarl=0, nFOTLVarh=0;

  // FILE *prs; prs = fopen("fake-print.txt","w");

  //
  // Set up output ntuple file for the sample
  //
  gSystem->mkdir(outputDir,kTRUE);
  TString outName = outputDir + TString("/") + TString("fakes-eps_select.root");
  TFile *outFile = new TFile(outName,"RECREATE");
  TTree outtree("Events","Events");

  EmuData data;
  Double_t rawMet,rawprojvar, trigeff=1;
  UInt_t npt20jets;
  const UInt_t kMaxPt20Jets=15;
  TArrayF btagArray; btagArray.Set(kMaxPt20Jets); // array to hold b-tag values for pt-20 jets
  Float_t varl,varh;
  outtree.Branch("Events",&data.runNum,
		 "runNum/i:evtNum:lumiSec:nPV:njets:nbjets:met/F:metphi:mass:scaledmass:dphi:mt:pt:phi:pmet:pvis:lpt1:leta1:lphi1:lpt2:leta2:lphi2:jpt1:jeta1:jphi1:jpt2:jeta2:jphi2:bjpt:bjeta:bjphi:mjj:svfmass:svflpt1:svfleta1:svflphi1:svflpt2:svfleta2:svflphi2:weight:state/I");

  // extra branches
  outtree.Branch("npt20jets",&npt20jets);
  outtree.Branch("btagArray",&btagArray);
  outtree.Branch("trigeff",&trigeff);
  outtree.Branch("rawMet",&rawMet);
  outtree.Branch("rawprojvar",&rawprojvar);
  outtree.Branch("varl",&varl);
  outtree.Branch("varh",&varh);

  TFile *infile=0;
  TTree *eventTree=0;
    
  // Data structures to store info from TTrees
  mithep::TEventInfo *info  = new mithep::TEventInfo();
  TClonesArray *muonArr     = new TClonesArray("mithep::TMuon");
  TClonesArray *electronArr = new TClonesArray("mithep::TElectron");
  TClonesArray *jetArr      = new TClonesArray("mithep::TJet");
  TClonesArray *pvArr       = new TClonesArray("mithep::TVertex");
  TClonesArray *svfitArr    = new TClonesArray("mithep::TSVFit");
  
  //
  // loop over data samples
  //
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {  
  
    // Read input file and get the TTrees
    cout << "Processing " << fnamev[ifile] << "..." << endl;
    infile = new TFile(fnamev[ifile]); 
    assert(infile);

    const TString json = jsonv[ifile];
    Bool_t hasJSON = kFALSE;
    mithep::RunLumiRangeMap rlrm;
    if(json.CompareTo("NONE")!=0) { 
      hasJSON = kTRUE;
      rlrm.AddJSONFile(json.Data()); 
    }
        
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);

    // Set branch address to structures that will store the info  
    eventTree->SetBranchAddress("Info",    &info);        TBranch *infoBr     = eventTree->GetBranch("Info");
    eventTree->SetBranchAddress("Muon",    &muonArr);     TBranch *muonBr     = eventTree->GetBranch("Muon");
    eventTree->SetBranchAddress("Electron",&electronArr); TBranch *electronBr = eventTree->GetBranch("Electron");
    eventTree->SetBranchAddress("PFJet",   &jetArr);      TBranch *jetBr      = eventTree->GetBranch("PFJet");
    eventTree->SetBranchAddress("PV",      &pvArr);       TBranch *pvBr       = eventTree->GetBranch("PV");
    eventTree->SetBranchAddress("SVFit",    &svfitArr);    TBranch *svfitBr    = eventTree->GetBranch("SVFit");
   
    Double_t weight=1;

    Double_t counter[30];         for(Int_t i=0; i<30; i++) { counter[i] = 0; }
   
    // perform fakeable objects extrapolation
    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      
      infoBr->GetEntry(ientry);

      // check for certified runs
      mithep::RunLumiRangeMap::RunLumiPairType rl(info->runNum, info->lumiSec);      
      if(hasJSON && !rlrm.HasRunLumi(rl)) continue;  
      
      // trigger requirement
      ULong_t trigger = 0;
      if(datatypev[ifile] == eMuEl) trigger = kHLT_Mu17_Ele8_CaloIdL | kHLT_Mu8_Ele17_CaloIdL;
//      if(datatypev[ifile] == eMuEl) trigger = kHLT_Mu17_Ele8_CaloIdL | kHLT_Mu8_Ele17_CaloIdL | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
      else {cout << "data type not defined" << endl; assert(0);}
      if(datatypev[ifile]!=eMC && !(info->triggerBits & trigger)) continue;
                
      if(!(info->hasGoodPV)) continue;
      
      pvArr->Clear();
      pvBr->GetEntry(ientry);
      
      // check for denominator objects and those passing tight selection
      vector<const mithep::TMuon*>     tightv;
      vector<const mithep::TElectron*> loosev;
      vector<Double_t> ratev;
      vector<Double_t> rateErrlv;
      vector<Double_t> rateErrhv;

      muonArr->Clear();
      muonBr->GetEntry(ientry);
      for(Int_t i=0; i<muonArr->GetEntriesFast(); i++) {
        const mithep::TMuon* muon = (mithep::TMuon*)((*muonArr)[i]);

        Bool_t trigmatch = muon->hltMatchBits & (kHLT_Mu17_Ele8_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdL_MuObj);        	
//	Bool_t trigmatch = muon->hltMatchBits & (kHLT_Mu17_Ele8_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj);
	// if(!trigmatch)            continue; // seems to be broken
        if(info->runNum<167000) // trigger matching broken after this run
          if(!trigmatch)                     continue;

	if(muon->pt < kMuonPt2Min)    continue;
	if(fabs(muon->eta) > 2.1)     continue;
	
	//if((info->triggerBits & kHLT_Mu17_Ele8_CaloIdL) && muon->pt < 20.0 ) continue;

	if(passMuonID(muon)) {
	  tightv.push_back(muon);
	  nTight+=weight;
	}
      }
      
      electronArr->Clear();
      electronBr->GetEntry(ientry);
      UInt_t ntightele=0;
      for(Int_t i=0; i<electronArr->GetEntriesFast(); i++) {
        const mithep::TElectron* electron = (mithep::TElectron*)((*electronArr)[i]);
	
	if(electron->pt        < kElePt2Min)  continue;
	if(fabs(electron->eta) > 2.5)         continue;
	//if(!(electron->isEcalDriven)) continue;

        Bool_t trigmatch = electron->hltMatchBits & (kHLT_Mu17_Ele8_CaloIdL_EGObj | kHLT_Mu8_Ele17_CaloIdL_EGObj);	
//	Bool_t trigmatch = electron->hltMatchBits & (kHLT_Mu17_Ele8_CaloIdL_EGObj | kHLT_Mu8_Ele17_CaloIdL_EGObj | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj);
        if(info->runNum<167000) // trigger matching broken after this run
	  if(!trigmatch)                          continue;
	// if(!(electron->hltMatchBits & trigger)) continue;
	//if((info->triggerBits & (kHLT_Mu8_Ele17_CaloIdL | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL)) && electron->pt < 20.0 ) continue;

	if(passEleID(electron)) {ntightele++; continue;}
	
	if(isEleFO(electron)) {
	  loosev.push_back(electron);
	  Double_t fopt = (electron->pt < 35) ? electron->pt : 34.99;
	  ratev.push_back(fr.getEff(fabs(electron->eta),fopt));
	  // assert(ratev.back() > 0);
	  rateErrlv.push_back(fr.getErrLow(fabs(electron->eta),fopt));
	  rateErrhv.push_back(fr.getErrHigh(fabs(electron->eta),fopt));
	}
      }
          
      // Tight+Loose
      if(tightv.size()==1) {
      // if((tightv.size() + ntightele)==1) {
	if((tightv.size() + ntightele)!=1) continue;

	const mithep::TMuon* mu = tightv[0];

	for(UInt_t i=0; i<loosev.size(); i++) {
	  const mithep::TElectron* ele = loosev[i];

          if(mu->q == ele->q) continue;                    
	  	  
	  if(mu->pt<kMuonPt1Min  && ele->pt<kElePt1Min) continue;

	  // if((info->triggerBits & kHLT_Mu8_Ele17_CaloIdL) && ele->pt < 20.0 )                     continue;

	  if(mu->pt < 20) {
            if(!(info->triggerBits & kHLT_Mu8_Ele17_CaloIdL)) continue; // if failed trig1
//	    if(!(info->triggerBits & (kHLT_Mu8_Ele17_CaloIdL | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL))) continue; // if failed trig1
	  }
	  else if(ele->pt < 20) {
	    if(!(info->triggerBits & kHLT_Mu17_Ele8_CaloIdL)) continue; // if failed trig2
	  }
	  
	  TLorentzVector lep1, lep2, dilep;  // lepton 4-vectors
	  Int_t finalState=-1;	           // final state type
          TLorentzVector svflep1, svflep2;  // lepton 4-vectors
          Double_t svfmass = -999;

	  if(mu->pt > ele->pt) {
	    lep1.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
	    lep2.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
	  } else {
	    lep1.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
	    lep2.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
	  }
	  dilep = lep1+lep2;

	  if(ele->pt > mu->pt) finalState=kEleMu; 
	  else                 finalState=kMuEle;

          // SVFit info
          svfitArr->Clear();
          svfitBr->GetEntry(ientry);
          mithep::TSVFit *svfit = (mithep::TSVFit*)((*svfitArr)[0]);
          svflep1.SetPtEtaPhiM(svfit->daughter1.Pt(), svfit->daughter1.Eta(), svfit->daughter1.Phi(), svfit->daughter1.M());
          svflep2.SetPtEtaPhiM(svfit->daughter2.Pt(), svfit->daughter2.Eta(), svfit->daughter2.Phi(), svfit->daughter2.M());
          svfmass = svfit->mass;

	  // loop through jets      
	  jetArr->Clear();
	  jetBr->GetEntry(ientry);
	  UInt_t njets = 0, nbjets = 0;
	  const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	  btagArray.Reset(); npt20jets=0;
	  for(Int_t ijet=0; ijet<jetArr->GetEntriesFast(); ijet++) {
	    const mithep::TJet *jet = (mithep::TJet*)((*jetArr)[ijet]);

	    if(toolbox::deltaR(jet->eta,jet->phi,lep1.Eta(),lep1.Phi()) < 0.3) continue;
	    if(toolbox::deltaR(jet->eta,jet->phi,lep2.Eta(),lep2.Phi()) < 0.3) continue;

	    if(fabs(jet->eta) > 5) continue;

	    // look for b-jets
	    if((jet->pt > kBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
	      assert(npt20jets<kMaxPt20Jets);
	      btagArray.AddAt(jet->tche,npt20jets); npt20jets++;
	      if(jet->tche > 3.3) {
		nbjets++;
		if(!bjet || jet->pt > bjet->pt)
		  bjet = jet; // leading b-jet
	      }
	    }

	    // look for vbf jets
	    if(jet->pt > kJetPtMin) {
	      njets++;
	      if(!jet1 || jet->pt > jet1->pt) { // jet1 is highest pt, jet2 next highest
		jet2 = jet1;
		jet1 = jet;
	      } else if(!jet2 || jet->pt > jet2->pt) {
		jet2 = jet;
	      }
	    }		    
	  }

	  TLorentzVector jv1, jv2, dijet;
	  if(njets>1) {
	    jv1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
	    jv2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
	    dijet = jv1+jv2;
	  } 

	  // calculate projection variables
	  TVector3 m,e,metv;
	  m.SetPtEtaPhi(mu->pt,0,mu->phi);
	  e.SetPtEtaPhi(ele->pt,0,ele->phi);
	  metv.SetPtEtaPhi(info->pfMET,0,info->pfMETphi); // uncorrected met
	  TVector3 bisector(m.Unit() + e.Unit());
	  bisector = bisector.Unit();
	  Double_t projVis  = (m+e).Dot(bisector);
	  Double_t projMet  =   metv.Dot(bisector);
	  rawprojvar  = 0.85*projVis - projMet;

	  rawMet = info->pfMET;
	  Double_t met=info->pfMET,metphi=info->pfMETphi;
	  metv.SetPtEtaPhi(met,0,metphi);
	  projMet  =   metv.Dot(bisector);

	  const Double_t r    = ratev[i]/(1. - ratev[i]);
	  if(ratev[i]<0) cout << "error: " << ele->pt << " " << ratev[i] << endl;
          const Double_t errl = rateErrlv[i]/(1.-ratev[i])/(1.-ratev[i]);
	  const Double_t errh = rateErrhv[i]/(1.-ratev[i])/(1.-ratev[i]);
	  
	  nFOTL     += weight*r;
	  nFOTLVarl += weight*weight*(r*r+errl*errl);
	  nFOTLVarh += weight*weight*(r*r+errh*errh);
	  
	  nDenomTL+=weight;
	  
	  varl    = weight*weight*(r*r+errl*errl);
	  varh    = weight*weight*(r*r+errh*errh);

	  // fprintf(prs,"%14d%14d%8.2f\n",info->runNum,info->evtNum,r);
 
	  data.runNum  = info->runNum;
	  data.evtNum  = info->evtNum;
	  data.lumiSec = info->lumiSec;
	  data.nPV     = pvArr->GetEntriesFast();
	  data.njets   = njets;
	  data.nbjets  = nbjets;
	  data.met     = met;
	  data.metphi  = metphi;
	  data.mass    = dilep.M();
          data.scaledmass = dilep.M();
	  data.dphi    = toolbox::deltaPhi(lep1.Phi(),lep2.Phi());
	  data.mt      = sqrt( 2.0 * (dilep.Pt()) * met * (1.0-cos(toolbox::deltaPhi(dilep.Phi(),metphi))) );
	  data.pt      = dilep.Pt();
	  data.phi     = dilep.Phi();
	  data.pmet    = projMet;
	  data.pvis    = projVis;
	  data.lpt1    = lep1.Pt();
	  data.leta1   = lep1.Eta();
	  data.lphi1   = lep1.Phi();
	  data.lpt2    = lep2.Pt();
	  data.leta2   = lep2.Eta();
	  data.lphi2   = lep2.Phi();
	  data.jpt1    = (jet1) ? jet1->pt  : 0;
	  data.jeta1   = (jet1) ? jet1->eta : 0;
	  data.jphi1   = (jet1) ? jet1->phi : 0;
	  data.jpt2    = (jet2) ? jet2->pt  : 0;
	  data.jeta2   = (jet2) ? jet2->eta : 0;
	  data.jphi2   = (jet2) ? jet2->phi : 0;
	  data.bjpt    = (bjet) ? bjet->pt  : 0;
	  data.bjeta   = (bjet) ? bjet->eta : 0;
	  data.bjphi   = (bjet) ? bjet->phi : 0;
	  data.mjj     = (njets>1) ? dijet.M() : 0;
          data.svfmass = svfmass;
          data.svflpt1 = svflep1.Pt();
          data.svfleta1= svflep1.Eta();
          data.svflphi1= svflep1.Phi();
          data.svflpt2 = svflep2.Pt();
          data.svfleta2= svflep2.Eta();
          data.svflphi2= svflep2.Phi();
	  data.weight  = weight*r;
	  data.state   = finalState;  	   

	  outtree.Fill();

	/*	  
          cout << endl;
          cout << ">> Run:" << info->runNum << ", Lumi:" << info->lumiSec << ", Event:" << info->evtNum << endl;  
          cout << ">> Tight: pt:" << tightv[0]->pt << ", eta:" << tightv[0]->eta << ", phi:" << tightv[0]->phi << endl;
          cout << ">> Loose: pt:" << loosev[i]->pt << ", eta:" << loosev[i]->eta << ", phi:" << loosev[i]->phi << endl;
          cout << ">> Weight: " << wgt << " (+" << sqrt(varh) << ", -" << sqrt(varl) << ")" << endl;
          cout << endl;
	*/	  
        }
      }      
    }    
    delete infile;
    infile=0, eventTree=0;
  }
  outFile->Write();
  outFile->Close();
  delete outFile;
  
  delete info;
  delete muonArr;
  delete electronArr;
  delete jetArr;
  delete pvArr;
  delete svfitArr;

  // fclose(prs);
  
  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;

  cout << "   data samples:" << endl;
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    cout << "     " << fnamev[ifile] << endl;
  }
  cout << endl;
  cout << "   fake rate data: " << frname << endl;
  cout << endl;
  cout << "  >>>   Events from TL: " << nFOTL << " (+" << sqrt(nFOTLVarh) << "/-" << sqrt(nFOTLVarl) << ")" << endl;
  cout << "  >>>      FOs from TL: " << nDenomTL << endl;  
  cout << "  >>>      Tight muons: " << nTight << endl;
  cout << endl;
  
  ofstream txtfile;
  char txtfname[100];
  sprintf(txtfname,"%s/summary.txt",outputDir.Data());
  txtfile.open(txtfname);
  txtfile << "*" << endl;
  txtfile << "* SUMMARY" << endl;
  txtfile << "*--------------------------------------------------" << endl;
  txtfile << endl;

  txtfile << "   data samples:" << endl;
  for(UInt_t ifile=0; ifile<fnamev.size(); ifile++) {
    txtfile << "     " << fnamev[ifile] << endl;
  }
  txtfile << endl;
  txtfile << "   fake rate data: " << frname << endl;
  txtfile << endl;
  txtfile << "  >>>   Events from TL: " << nFOTL << " (+" << sqrt(nFOTLVarh) << "/-" << sqrt(nFOTLVarl) << ")" << endl;
  txtfile << "  >>>      FOs from TL: " << nDenomTL << endl;  
  txtfile << "  >>>      Tight muons: " << nTight << endl;
  txtfile << endl;
  txtfile.close();
  
  cout << " <> Output saved in " << outputDir << "/" << endl;    
  cout << endl; 


  
  gBenchmark->Show("predictFakesMuEle"); 
} 
