#include "RecoilCorrector.hh"
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"

void runRecoil()
{
  // ztt MC that we're going to correct
  TFile *fZMCFileNVtx = new TFile("/scratch/dkralph/Recoil/ztt_select_recoil_novtx.root");
  TTree *fZMCTreeNVtx = (TTree*) fZMCFileNVtx->FindObjectAny("Events");

  double lGenPt  = 0; double lGenPhi  = 0; 
  double lLepPt1 = 0; double lLepPhi1 = 0;  double lLepEta1 = 0;
  double lLepPt2 = 0; double lLepPhi2 = 0;  double lLepEta2 = 0;
  double lMet    = 0; double lMPhi    = 0;  double lWeight  = 0;
  fZMCTreeNVtx->SetBranchAddress("vpt"   ,&lGenPt); // gen z pt
  fZMCTreeNVtx->SetBranchAddress("vphi"  ,&lGenPhi); // gen z phi
  fZMCTreeNVtx->SetBranchAddress("lpt1"  ,&lLepPt1); // lepton 1 kinematics
  fZMCTreeNVtx->SetBranchAddress("lphi1" ,&lLepPhi1);
  fZMCTreeNVtx->SetBranchAddress("leta1" ,&lLepEta1);
  fZMCTreeNVtx->SetBranchAddress("lpt2"  ,&lLepPt2); // lepton 2 kinematics
  fZMCTreeNVtx->SetBranchAddress("lphi2" ,&lLepPhi2);
  fZMCTreeNVtx->SetBranchAddress("leta2" ,&lLepEta2);
  fZMCTreeNVtx->SetBranchAddress("met"   ,&lMet);    // reco met
  fZMCTreeNVtx->SetBranchAddress("metphi",&lMPhi);   // reco met phi

  RecoilCorrector corrector;

  for(int i0 = 0; i0 < fZMCTreeNVtx->GetEntries(); i0++) {
    fZMCTreeNVtx->GetEntry(i0);

    TLorentzVector m,e;
    m.SetPtEtaPhiM(lLepPt1,lLepEta1,lLepPhi1,0.105658369);
    e.SetPtEtaPhiM(lLepPt2,lLepEta2,lLepPhi2,0.000511);

    double pMet  = lMet;
    double pMPhi = lMPhi; 
    printf("====> Met Before ==> %10.2f%10.2f\n",pMet,pMPhi);
    corrector.Correct(pMet,pMPhi,lGenPt,lGenPhi,(m+e).Pt(),(m+e).Phi());
    printf("   => Met After  ==> %10.2f%10.2f\n",pMet,pMPhi);
    if(i0>100) break;
  }


}
