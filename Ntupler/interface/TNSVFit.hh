#ifndef EWKANA_NTUPLER_TNSVFIT_HH
#define EWKANA_NTUPLER_TNSVFIT_HH

#include <Math/Vector4D.h>
#include <TObject.h>
//#include "MitCommon/DataFormats/interface/Vect4M.h"

namespace mithep 
{
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double_t> > FourVectorM;

  class TNSVFit : public TObject
  {
  public:
    TNSVFit() {}
	
    ~TNSVFit() {} 
  
    Bool_t isValid;
    Double_t mass;
    Double_t massErrUp;
    Double_t massErrDown;
    Double_t massMean;
    Double_t massMedian;
    Double_t massMaximum;
    Double_t massMaxInterpol;
    Double_t cov_00;
    Double_t cov_10;
    Double_t cov_01;
    Double_t cov_11;
    FourVectorM daughter1, daughter2; //Note about daughter1 is always a lepton

    ClassDef(TNSVFit,2)
  };  
}
#endif
