#ifndef MITHTT_NTUPLER_TNSVFIT_H
#define MITHTT_NTUPLER_TNSVFIT_H

#include <TObject.h>
#include <Math/Vector4D.h>

/**
   \class TSVfit TSVfit.h MitHtt/Ntupler/include/TSVfit.h

   \brief Description: <one line class summary>

   <Notes on implementation>
*/

namespace mithep 
{
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double_t> > FourVectorM;

  class TSVfit : public TObject
  {
  public:
    /// default contstructor
    TSVfit(){}
    /// default destructor
    ~TSVfit(){} 
  
    /// indicated whether has the fit been successful or not
    bool isValid;
    /// return fitted mass
    double mass;
    /// return fit uncertainty in up direction
    double massErrUp;
    /// return fit uncertainty in down direction
    double massErrDown;
    /// obsolete
    double massMean;
    /// obsolete
    double massMedian;
    /// obsolete
    double massMaximum;
    /// obsolete
    double massMaxInterpol;
    /// met significance matrix element [x|x]
    double cov_00;
    /// met significance matrix element [y|x]
    double cov_10;
    /// met significance matrix element [x|y]
    double cov_01;
    /// met significance matrix element [y|y]
    double cov_11;
    /// di-lepton resonance daughters. NOTE: daughter1 in doubt is always taken to be an electron/muon
    FourVectorM daughter1, daughter2;

    /// root class definition
    ClassDef(TSVfit,2)
  };  
}
#endif
