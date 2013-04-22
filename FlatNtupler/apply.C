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

double deltaPhi(double iPhi1,double iPhi2) { 
  double pDPhi = fabs(iPhi1-iPhi2);
  if( 2.*TMath::Pi()-pDPhi < pDPhi) pDPhi = 2.*TMath::Pi()-pDPhi;
  return pDPhi;
}

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
   int ngamma1 = 0; int ngamma2 = 0; int nprong1 = 0; int nprong2 = 0;
   
   lTree->SetBranchAddress("ngamma_1"       ,&ngamma1);
   lTree->SetBranchAddress("ngamma_2"       ,&ngamma2);
   lTree->SetBranchAddress("nprong_1"       ,&nprong1);
   lTree->SetBranchAddress("nprong_2"       ,&nprong2);

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

   float lMSV          = 0;
   float lPt           = 0;
   float lPhi          = 0;
   float lEta          = 0;
   float lMSVunc       = 0;
   float lPtunc        = 0;
   float lPhiunc       = 0;
   float lEtaunc       = 0;
   float lPtHigh       = 0;
   float lPhiHigh      = 0;
   float lEtaHigh      = 0;
   float lMSVuncHigh   = 0;
   float lPtuncHigh    = 0;
   float lPhiuncHigh   = 0;
   float lEtauncHigh   = 0;
   float lPtLow        = 0;
   float lPhiLow       = 0;
   float lEtaLow       = 0;
   float lMSVuncLow    = 0;
   float lPtuncLow     = 0;
   float lPhiuncLow    = 0;
   float lEtauncLow    = 0;
   float lMSVHigh      = 0;
   float lMVisHigh     = 0;
   float lMt1High      = 0;
   float lMt2High      = 0;
   float lMetHigh      = 0;
   float lMetPhiHigh   = 0;
   float lPt1High      = 0;
   float lPt2High      = 0;
   float lM1High       = 0;
   float lM2High       = 0;
   float lMSVLow       = 0;
   float lMVisLow      = 0;
   float lMt1Low       = 0;
   float lMt2Low       = 0;
   float lMetLow       = 0;
   float lMetPhiLow    = 0;
   float lPt1Low       = 0;
   float lPt2Low       = 0;
   float lM1Low        = 0;
   float lM2Low        = 0;
 

   lOTree->Branch("m_sv"          ,&lMSV         ,"lMSV/F");
   lOTree->Branch("pt_sv"          ,&lPt         ,"lPt/F");
   lOTree->Branch("phi_sv"          ,&lPhi         ,"lPhi/F");
   lOTree->Branch("eta_sv"          ,&lEta         ,"lEta/F");
   lOTree->Branch("munc_sv"          ,&lMSVunc         ,"lMSVunc/F");
   lOTree->Branch("ptunc_sv"          ,&lPtunc         ,"lPtunc/F");
   lOTree->Branch("phiunc_sv"          ,&lPhiunc         ,"lPhiunc/F");
   lOTree->Branch("munc_sv"          ,&lEtaunc         ,"lEtaunc/F");

   lOTree->Branch("pt_svhigh"          ,&lPtHigh         ,"lPtHigh/F");
   lOTree->Branch("phi_svhigh"          ,&lPhiHigh         ,"lPhiHigh/F");
   lOTree->Branch("eta_svhigh"          ,&lEtaHigh         ,"lEtaHigh/F");
   lOTree->Branch("munc_svhigh"          ,&lMSVuncHigh         ,"lMSVuncHigh/F");
   lOTree->Branch("ptunc_svhigh"          ,&lPtuncHigh         ,"lPtuncHigh/F");
   lOTree->Branch("phiunc_svhigh"          ,&lPhiuncHigh         ,"lPhiuncHigh/F");
   lOTree->Branch("munc_svhigh"          ,&lEtauncHigh         ,"lEtauncHigh/F");
   
   lOTree->Branch("pt_svlow"          ,&lPtLow         ,"lPtLow/F");
   lOTree->Branch("phi_svlow"          ,&lPhiLow         ,"lPhiLow/F");
   lOTree->Branch("eta_svlow"          ,&lEtaLow        ,"lEtaLow/F");
   lOTree->Branch("munc_svlow"          ,&lMSVuncLow         ,"lMSVuncLow/F");
   lOTree->Branch("ptunc_svlow"          ,&lPtuncLow         ,"lPtuncLow/F");
   lOTree->Branch("phiunc_svlow"          ,&lPhiuncLow         ,"lPhiuncLow/F");
   lOTree->Branch("munc_svlow"          ,&lEtauncLow         ,"lEtauncLow/F");

   lOTree->Branch("m_svhigh"      ,&lMSVHigh     ,"lMSVHigh/F");
   lOTree->Branch("mtMVA_1high"   ,&lMt1High     ,"lMt1High/F");
   lOTree->Branch("mtMVA_2high"   ,&lMt2High     ,"lMt2High/F");
   lOTree->Branch("mvamethigh"    ,&lMetHigh     ,"lMetHigh/F");
   lOTree->Branch("mvametphihigh" ,&lMetPhiHigh  ,"lMetPhiHigh/F");
   lOTree->Branch("pt_1high"      ,&lPt1High     ,"lPt1High/F");
   lOTree->Branch("pt_2high"      ,&lPt2High     ,"lPt2High/F");
   lOTree->Branch("m_1high"       ,&lM1High      ,"lM1High/F");
   lOTree->Branch("m_2high"       ,&lM2High      ,"lM2High/F");
   lOTree->Branch("m_vishigh"     ,&lMVisHigh    ,"lMVisHigh/F");  

   lOTree->Branch("m_svlow"       ,&lMSVLow      ,"lMSVLow/F");
   lOTree->Branch("mtMVA_1low"    ,&lMt1Low      ,"lMt1Low/F");
   lOTree->Branch("mtMVA_2low"    ,&lMt2Low      ,"lMt2Low/F");
   lOTree->Branch("mvametlow"     ,&lMetLow      ,"lMetLow/F");
   lOTree->Branch("mvametphilow"  ,&lMetPhiLow   ,"lMetPhiLow/F");
   lOTree->Branch("pt_1low"       ,&lPt1Low      ,"lPt1Low/F");
   lOTree->Branch("pt_2low"       ,&lPt2Low      ,"lPt2Low/F");
   lOTree->Branch("m_1low"       ,&lM1Low      ,"lM1Low/F");
   lOTree->Branch("m_2low"       ,&lM2Low      ,"lM2Low/F");
   lOTree->Branch("m_vislow"      ,&lMVisLow     ,"lMVisLow/F");      

   double lCorr1 = 0;
   double lCorr2 = 0;
   if(channel == 0) lCorr2  = 0.01;
   if(channel == 1) lCorr2  = 0.03;
   if(channel == 2)
     { 
       lCorr1 = 0.03;
       lCorr2 = 0.03;
     }
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
     lMSV =  fitter->integrateMarkov(&svfit,lMet,lMetPhi,channel);
     lPt  =  fitter->GetPt();
     lPhi = fitter->GetPhi();
     lEta = fitter->GetEta();
     lMSVunc = fitter->massUnc();
     lPtunc = fitter->GetPtUnc();
     lPhiunc = fitter->GetPhiUnc();
     lEtaunc = fitter->GetEtaUnc();

     TLorentzVector lVMetHigh; lVMetHigh.SetPtEtaPhiM(lMet,0,lMetPhi,0);
     TLorentzVector lVMetLow ; lVMetLow.SetPtEtaPhiM(lMet,0,lMetPhi,0);
     lPt1High = lPt1*(1.+lCorr1);
     lPt1Low  = lPt1*(1.-lCorr1);
     lPt2High = lPt2*(1.+lCorr2);
     lPt2Low  = lPt2*(1.-lCorr2);
     TLorentzVector ldvec1; 
     if(nprong1==3 || ngamma1>0)
       {
	 ldvec1.SetPtEtaPhiM(lPt1*lCorr1,lEta1,lPhi1,lM1*lCorr1);
	 lM1High = (1.+lCorr1)*lM1;
	 lM1Low  = (1.-lCorr1)*lM1;
       }
     else
       {
	 ldvec1.SetPtEtaPhiM(lPt1*lCorr1,lEta1,lPhi1,lM1*lCorr1);
	 lM1High = lM1;
	 lM1Low  = lM1;
       }
     TLorentzVector ldvec2; 
     if(nprong2==3 || ngamma2>0) 
       {
	 ldvec2.SetPtEtaPhiM(lPt2*lCorr2,lEta2,lPhi2,lM2*lCorr2);
	 lM2High = (1.+lCorr2)*lM2;
	 lM2Low  = (1.-lCorr2)*lM2;
       }
     else
       {
	 ldvec2.SetPtEtaPhiM(lPt2*lCorr2,lEta2,lPhi2,lM2); 
	 lM2High = lM2;
	 lM2Low  = lM2;
       }
	 
     lVMetHigh -= ldvec1;
     lVMetHigh -= ldvec2;
     lVMetLow  += ldvec1;
     lVMetLow  += ldvec2;
    
     lMetHigh     = lVMetHigh.Pt();
     lMetLow      = lVMetLow .Pt();

     lMetPhiHigh  = lVMetHigh.Phi();
     lMetPhiLow   = lVMetLow .Phi();

     TLorentzVector lvec1High; lvec1High.SetPtEtaPhiM(lPt1High,lEta1,lPhi1,lM1High);
     TLorentzVector lvec1Low ; lvec1Low.SetPtEtaPhiM(lPt1Low ,lEta1,lPhi1,lM1Low);
     TLorentzVector lvec2High; lvec2High.SetPtEtaPhiM(lPt2High,lEta2,lPhi2,lM2High);
     TLorentzVector lvec2Low ; lvec2Low.SetPtEtaPhiM(lPt2Low ,lEta2,lPhi2,lM2Low);

     lMVisHigh=(lvec1High+lvec2High).M();
     lMVisLow =(lvec1Low+lvec2Low).M();

     lMt1High = sqrt(2.0*(lvec1High.Pt()*lMetHigh*(1.0-cos(deltaPhi(lvec1High.Phi(),lMetPhiHigh)))));
     lMt1Low  = sqrt(2.0*(lvec1Low .Pt()*lMetLow *(1.0-cos(deltaPhi(lvec1Low .Phi(),lMetPhiLow )))));

     lMt2High = sqrt(2.0*(lvec2High.Pt()*lMetHigh*(1.0-cos(deltaPhi(lvec2High.Phi(),lMetPhiHigh)))));
     lMt2Low  = sqrt(2.0*(lvec2Low .Pt()*lMetLow *(1.0-cos(deltaPhi(lvec2Low .Phi(),lMetPhiLow )))));
     
     svlep1.SetPxPyPzE(lvec1High.Px(),lvec1High.Py(),lvec1High.Pz(),lvec1High.E());
     svlep2.SetPxPyPzE(lvec2High.Px(),lvec2High.Py(),lvec2High.Pz(),lvec2High.E());
     svfit.daughter1 = svlep1;
     svfit.daughter2 = svlep2;
     svfit.daughterId1 = 1;
     svfit.daughterId2 = 2;
     lMSVHigh =  fitter->integrateMarkov(&svfit,lMetHigh,lMetPhiHigh,channel);
     lPtHigh  =  fitter->GetPt();
     lPhiHigh = fitter->GetPhi();
     lEtaHigh = fitter->GetEta();
     lMSVuncHigh = fitter->massUnc();
     lPtuncHigh = fitter->GetPtUnc();
     lPhiuncHigh = fitter->GetPhiUnc();
     lEtauncHigh = fitter->GetEtaUnc();

     svlep1.SetPxPyPzE(lvec1Low.Px(),lvec1Low.Py(),lvec1Low.Pz(),lvec1Low.E());
     svlep2.SetPxPyPzE(lvec2Low.Px(),lvec2Low.Py(),lvec2Low.Pz(),lvec2Low.E());
     svfit.daughter1 = svlep1;
     svfit.daughter2 = svlep2;
     svfit.daughterId1 = 1;
     svfit.daughterId2 = 2;
     lMSVLow =  fitter->integrateMarkov(&svfit,lMetLow,lMetPhiLow,channel);
     lPtLow  =  fitter->GetPt();
     lPhiLow = fitter->GetPhi();
     lEtaLow = fitter->GetEta();
     lMSVuncLow = fitter->massUnc();
     lPtuncLow = fitter->GetPtUnc();
     lPhiuncLow = fitter->GetPhiUnc();
     lEtauncLow = fitter->GetEtaUnc();

     cout << "===> " << lMSV << " -- " << lMSVLow << " -- " << lMSVHigh << endl;
     lOTree->Fill();
     //delete svfit;
   }  
   
   lOTree->Write();
   lOFile->Close();
 }
