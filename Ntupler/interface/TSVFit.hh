#ifndef MITHTT_NTUPLER_TSVFIT_HH
#define MITHTT_NTUPLER_TSVFIT_HH

#include <TObject.h>
#include "MitCommon/DataFormats/interface/Vect4M.h"

namespace mithep 
{
  class TSVFit : public TObject
  {
    public:
      TSVFit(){}
      ~TSVFit(){}

      Bool_t  isValid;
      Double_t mass, massErrUp, massErrDown;
      Double_t massMean, massMedian, massMaximum, massMaxInterpol;
      FourVectorM daughter1, daughter2;

    ClassDef(TSVFit,1)
  };
}
#endif
