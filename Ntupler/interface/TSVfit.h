#ifndef MITHTT_NTUPLER_TNSVFIT_H
#define MITHTT_NTUPLER_TNSVFIT_H

#include <TObject.h>
#include <Math/Vector4D.h>

/**
   \class TSVfit TSVfit.h MitHtt/Ntupler/include/TSVfit.h

   \brief Description: Bacon svfit input

   All information that is available on Bacon to calculate svfit for a given lepton pair.
   In principle all this informatino can be recomputed on any level if the particle flow 
   candidates and the particle flow jets are kept.
*/

namespace mithep 
{
  /// four vector with (pt, eta, phi, m) as coordinates
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double_t> > FourVectorM;

  class TSVfit : public TObject
  {
  public:
    /// default contstructor
    TSVfit(){}
    /// default destructor
    ~TSVfit(){} 
  
    /// met significance matrix element [x|x]
    double cov_00;
    /// met significance matrix element [y|x]
    double cov_10;
    /// met significance matrix element [x|y]
    double cov_01;
    /// met significance matrix element [y|y]
    double cov_11;
    //MVA met 
    double mvaMET, mvaMETphi;
    //MVA met 53X
    double mvaMET53, mvaMETphi53;
    /// MVA met significance matrix element [x|x] (53)
    double mvacov_00,mvacov_0053;
    /// MVA met significance matrix element [y|x]
    double mvacov_10,mvacov_1053;
    /// MVA met significance matrix element [x|y]
    double mvacov_01,mvacov_0153;
    /// MVA met significance matrix element [y|y]
    double mvacov_11,mvacov_1153;
    /// EGenType Id of daughter1 and daughter2
    unsigned int daughterId1, daughterId2;
    /// di-lepton resonance daughters.
    FourVectorM daughter1, daughter2;

    ClassDef(TSVfit, 2);
  };  
}
#endif
