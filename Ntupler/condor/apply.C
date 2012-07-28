#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TLorentzVector.h"
#include "TROOT.h"

#include "MitHtt/Ntupler/interface/TSVfit.h"
#include "MitHtt/Ntupler/interface/TSVfitter.h"

void apply(int id,int channel, const char* dataset,const char* infile,int nevents) 
 { 

   mithep::TSVfitter *fitter = new mithep::TSVfitter();
   std::cout << "Classifying" << std::endl;
   TFile *lFile = new TFile(infile);
   TTree *lTree = (TTree*) lFile->FindObjectAny("Events");

   float lMet = 0;
   float lMetPhi = 0;
   float lcov00 = 0;
   float lcov10 = 0;
   float lcov01 = 0;
   float lcov11 = 0;
   
   lTree->SetBranchAddress("mvamet"      ,&lMet);
   lTree->SetBranchAddress("mvametphi"   ,&lMetPhi);
   lTree->SetBranchAddress("mvacov00"    ,&lcov00);
   lTree->SetBranchAddress("mvacov10"    ,&lcov10);
   lTree->SetBranchAddress("mvacov01"    ,&lcov01);
   lTree->SetBranchAddress("mvacov11"    ,&lcov11);
 
   float lPt1 = 0; float lPhi1 = 0; float lEta1 = 0; float lM1 = 0; 
   float lPt2 = 0; float lPhi2 = 0; float lEta2 = 0; float lM2 = 0; 
   lTree->SetBranchAddress("pt_1"       ,&lPt1);
   lTree->SetBranchAddress("eta_1"      ,&lEta1);
   lTree->SetBranchAddress("phi_1"      ,&lPhi1);
   lTree->SetBranchAddress("m_1"        ,&lM1);
   lTree->SetBranchAddress("pt_2"       ,&lPt2);
   lTree->SetBranchAddress("eta_2"      ,&lEta2);
   lTree->SetBranchAddress("phi_2"      ,&lPhi2);
   lTree->SetBranchAddress("m_2"        ,&lM2);

   int lNEvents = lTree->GetEntries();
   char output[100];
   //char* idc = (char*)id;
   sprintf(output,"%s_%d_ntuple.root",dataset,id); 
   TFile *lOFile = new TFile(output,"RECREATE");
   TTree *lOTree = lTree->CloneTree(0);

   float lMSV        = 0;
   lOTree->Branch("m_sv"     ,&lMSV         ,"lMSV/F");
   
   for (Long64_t i0=0; i0<nevents; i0++) {
     if(i0 + id*nevents > lNEvents-1) continue;
     lTree->GetEntry(i0 + id*nevents);
     //mithep::TSVfit * svfit = new mithep::TSVfit();
     mithep::TSVfit svfit;
     if(lcov00*lcov11-lcov01*lcov01 == 0) continue;
     svfit.cov_00 = lcov00;
     svfit.cov_11 = lcov11;
     svfit.cov_01 = lcov01;
     svfit.cov_10 = lcov01;
     TLorentzVector lvec1; lvec1.SetPtEtaPhiM(lPt1,lEta1,lPhi1,lM1);
     TLorentzVector lvec2; lvec2.SetPtEtaPhiM(lPt2,lEta2,lPhi2,lM2);
     mithep::FourVectorM svlep1; svlep1.SetPxPyPzE(lvec1.Px(),lvec1.Py(),lvec1.Pz(),lvec1.E());
     mithep::FourVectorM svlep2; svlep2.SetPxPyPzE(lvec2.Px(),lvec2.Py(),lvec2.Pz(),lvec2.E());
     svfit.daughter1 = svlep1;
     svfit.daughter2 = svlep2;
     svfit.daughterId1 = 1;
     svfit.daughterId2 = 2;
     lMSV =  fitter->integrate(&svfit,lMet,lMetPhi,channel);
     lOTree->Fill();
     //delete svfit;
   }  
   
   lOTree->Write();
   lOFile->Close();
 }
