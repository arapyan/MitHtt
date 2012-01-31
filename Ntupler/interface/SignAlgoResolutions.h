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

/**\class SignAlgoResolutions SignAlgoResolutions.h MitHtt/Ntupler/include/SignAlgoResolutions.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/

namespace mithep {
  enum resolutionType { caloEE, caloEB, caloHE, caloHO, caloHF, caloHB, jet, electron, tau, muon,PFtype1,PFtype2, PFtype3, PFtype4, PFtype5, PFtype6, PFtype7 };
  enum resolutionFunc { ET, PHI,TRACKP,CONSTPHI };
  
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
    double getfunc(const resolutionType& type,const resolutionFunc& func,  std::vector<double>& x) const;
    void addfunction(const resolutionType type, const resolutionFunc func, std::vector<double> parameters);
    void initializeJetResolutions();
    
    typedef std::pair<resolutionType, resolutionFunc> functionCombo;
    typedef std::vector<double> functionPars;
    std::map<functionCombo,functionPars> functionmap_;
 
    double EtFunction( const functionPars &x,  const functionPars &  par) const;
    double PhiFunction( const functionPars &x,  const functionPars & par) const;
    double PFunction( const functionPars &x, const functionPars &par) const;
    double PhiConstFunction(const functionPars &x, const functionPars &par) const;

    double ptResolThreshold_;
    //temporary fix for low pT jet resolutions
    //First index, eta bins, from 0 to 5;
    //Second index, pt bins, from 3 to 23 GeV;
    std::vector<double> jdpt[10];
    std::vector<double> jdphi[10];
     
    JetResolution *ptResol_;
    JetResolution *phiResol_;
    std::string fName;
  };
}
#endif
