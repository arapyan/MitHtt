#ifndef TSVFITTER_HH
#define TSVFITTER_HH

#include <TObject.h>
#include "TLorentzVector.h"
#include "TNSVFit.hh"

namespace mithep 
{
  class TSVFitter : public TObject
  {
  public:
    TSVFitter() {}
    ~TSVFitter() {} 
    TLorentzVector fit(TNSVFit *iFit,double iMet,double iMetPhi);
    
    ClassDef(TSVFitter,1)
  };  
}
#endif
