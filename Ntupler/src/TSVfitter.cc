#include "TVectorD.h"
#include "TLorentzVector.h"
#include "MitHtt/Ntupler/interface/TSVfitter.h"

#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneLikelihood.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"

ClassImp(mithep::TSVfitter)


TLorentzVector mithep::TSVfitter::fit(mithep::TSVfit* fit, double iMet, double iMetPhi) 
{
  // setup significance matrix
  TMatrixD lMM(2,2); 
  lMM(0,0) = fit->cov_00; 
  lMM(0,1) = fit->cov_01; 
  lMM(1,0) = fit->cov_10; 
  lMM(1,1) = fit->cov_11; 
  // setup lorentz vectors for input leptons
  NSVfitStandalone::LorentzVector lLep1; lLep1.SetPxPyPzE( fit->daughter1.Px(), fit->daughter1.Py(), fit->daughter1.Pz(), fit->daughter1.E() );
  NSVfitStandalone::LorentzVector lLep2; lLep2.SetPxPyPzE( fit->daughter2.Px(), fit->daughter2.Py(), fit->daughter2.Pz(), fit->daughter2.E() );
  //std::cout << "---> Mass Check " << fit->daughter1.M() << " -- " << fit->daughter2.M() << std::endl;

  // setup lepton types (depending on mass of the input particles)
  NSVfitStandalone::kDecayType lId2 = NSVfitStandalone::kLepDecay; 
  if(fit->daughter2.M() != 0.105658 || fit->daughter2.M() != 0.000510999){  
    lId2 = NSVfitStandalone::kHadDecay;
  }
  // setup MET vector
  TVector3 buffer; buffer.SetPtThetaPhi(iMet,0,iMetPhi);
  NSVfitStandalone::Vector lMet(buffer.x(), buffer.y(), buffer.z());
  // configure the fit
  std::vector<NSVfitStandalone::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(NSVfitStandalone::kLepDecay, lLep1));
  measuredTauLeptons.push_back(NSVfitStandalone::MeasuredTauLepton(lId2                       , lLep2));
  NSVfitStandaloneAlgorithm algo(measuredTauLeptons,lMet, lMM, 0);
  algo.maxObjFunctionCalls(5000);
  algo.fit(); 
  NSVfitStandalone::LorentzVector lL  = algo.fittedDiTauSystem(); TLorentzVector result;
  result.SetXYZM(lL.x(), lL.y(), lL.z(), lL.mass());
  //std::cout << "===> Fit check " << fit->mass << " -- " << lL.M() << std::endl;
  return result;
} 
