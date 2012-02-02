#ifndef MITHTT_NTUPLER_SIGNALGORESOLUTIONS_H
#define MITHTT_NTUPLER_SIGNALGORESOLUTIONS_H

#include <map>
#include <vector>
#include <iostream>

#include "DataFormats/Math/interface/LorentzVector.h"
#include "RecoMET/METAlgorithms/interface/SigInputObj.h"
#include "CondFormats/JetMETObjects/interface/JetResolution.h"

#include "MitAna/DataTree/interface/PFJet.h"
#include "MitAna/DataTree/interface/PFCandidate.h"

/**
   \class SignAlgoResolutions SignAlgoResolutions.h MitHtt/Ntupler/include/SignAlgoResolutions.h

   \brief Description: class to determine object resolutions for MET significance

   This is a class to determine object resolutions for the calculation of the MET significance 
   used for the svfit calculation. The input parameters for the resolution parametrizations are 
   taken from: 

   $CMSSW_RELEASE_BASE/src/RecoMET/METProducers/python/METSigParams_cfi.py

   This class is a mirrir of the same class in RecoMET/METAlgorithms replacing the official data 
   formats with Bambu data types.
*/

namespace mithep {
  /// type of resolution (for internal use)
  enum resolutionType { caloEE, caloEB, caloHE, caloHO, caloHF, caloHB, jet, electron, tau, muon,PFtype1,PFtype2, PFtype3, PFtype4, PFtype5, PFtype6, PFtype7 };
  /// type of resolution function (for internal use)
  enum resolutionFunc { ET, PHI, TRACKP,CONSTPHI };
  
  class SignAlgoResolutions{
    
  public:
    /// default constructor
    SignAlgoResolutions();

    /// add resolutions
    void addResolutions();
    /// evaluate resolution 
    double eval(const resolutionType& type, const resolutionFunc& func, const double& et, const double& phi, const double& eta, const double& p) const;
    /// evaluate resolution 
    double eval(const resolutionType& type, const resolutionFunc& func, const double& et, const double& phi, const double& eta) const;
    /// calculate significance (bypassed to use Bambu objects)
    metsig::SigInputObj evalPF(const mithep::PFCandidate* candidate) const;
    /// calculate significance (bypassed to use Bambu objects)
    metsig::SigInputObj evalPFJet(const mithep::PFJet *jet) const;
    /// check whether functionmap is filled
    bool isFilled() const {return functionmap_.size()>0;}

  private:
    /// get function for given functiuon and resolution type 
    double getfunc(const resolutionType& type,const resolutionFunc& func,  std::vector<double>& x) const;
    /// add a function of given function and resolution type with given parameters
    void addfunction(const resolutionType type, const resolutionFunc func, std::vector<double> parameters);
    /// initialize jet resolutions
    void initializeJetResolutions();
    
    /// combined resolution type and function type and 
    typedef std::pair<resolutionType, resolutionFunc> functionCombo;
    /// function parameters
    typedef std::vector<double> functionPars;
    /// mappgin of function parameters to a function of given resolution and function type 
    std::map<functionCombo,functionPars> functionmap_;

    /// resolution in Et
    double EtFunction( const functionPars &x,  const functionPars &  par) const;
    /// resolutino in phi
    double PhiFunction( const functionPars &x,  const functionPars & par) const;
    /// resolutino for particle flow candidates
    double PFunction( const functionPars &x, const functionPars &par) const;
    /// const resolution term in phi
    double PhiConstFunction(const functionPars &x, const functionPars &par) const;

    /// pt threshold for the pobjects in consideration    
    double ptResolThreshold_;
    //temporary fix for low pT jet resolutions
    //First index, eta bins, from 0 to 5;
    //Second index, pt bins, from 3 to 23 GeV;
    std::vector<double> jdpt[10];
    std::vector<double> jdphi[10];
    
    /// helper objects to determine the jet resolutions from JetMET in pt
    JetResolution *ptResol_;
    /// helper objects to determine the jet resolutions from JetMET in phi
    JetResolution *phiResol_;
    /// name of the input file for the function parameters
    std::string fName;
  };
}
#endif
