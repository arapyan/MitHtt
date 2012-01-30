#include "TLorentzVector.h"
#include "../interface/NSVfitStandaloneLikelihood.h"
#include "../interface/NSVfitStandaloneAlgorithm.h"
#include "../interface/TSVFitter.hh"
ClassImp(mithep::TSVFitter)

TLorentzVector mithep::TSVFitter::fit(mithep::TNSVFit *iFit,double iMet,double iMetPhi) {
  NSVfitStandalone::LorentzVector lL;
  TMatrixD lMM(2,2); 
  lMM(0,0) = iFit->cov_00; 
  lMM(0,1) = iFit->cov_01; 
  lMM(1,0) = iFit->cov_10; 
  lMM(1,1) = iFit->cov_11; 
  NSVfitStandalone::LorentzVector lLep1;
  NSVfitStandalone::LorentzVector lLep2;
  lLep1.SetPxPyPzE(iFit->daughter1.Px(),iFit->daughter1.Py(),iFit->daughter1.Pz(),iFit->daughter1.E());
  lLep2.SetPxPyPzE(iFit->daughter2.Px(),iFit->daughter2.Py(),iFit->daughter2.Pz(),iFit->daughter2.E());
  std::cout << "---> Mass Check " << iFit->daughter1.M() << " -- " << iFit->daughter2.M() << std::endl;
  NSVfitStandalone::kDecayType lId2 = NSVfitStandalone::kLepDecay; 
  if(iFit->daughter2.M() != 0.105658 || iFit->daughter2.M() != 0.000510999)  lId2 = NSVfitStandalone::kHadDecay;
  NSVfitStandalone::Vector         lMet;   lMet.SetPtThetaPhi(iMet,0,iMetPhi);
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, lLep1));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(lId2                       , lLep2));
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons,lMet, lMM, 0);
  algo.maxObjFunctionCalls(5000);
  algo.fit(); 
  lL  = algo.fittedDiTauSystem();
  //std::cout << "===> Fit check " << iFit->mass << " -- " << lL.M() << std::endl;
  return lL;
} 
