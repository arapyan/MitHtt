#ifndef MITHTT_NTUPLER_TSVFITTER_HH
#define MIT_HTT_NTUPLER_TSVFITTER_HH

#include "TSVfit.h"
#include "TObject.h"
#include "TLorentzVector.h"

/**
   \class TSVfit TSVfit.h MitHtt/Ntupler/include/TSVfit.h

   \brief Description: <one line class summary>

   <Notes on implementation>
*/

namespace mithep 
{
  class TSVfitter : public TObject
  {
  public:
    /// default contructor
    TSVfitter(){}
    ///default destructor
    ~TSVfitter(){}
    /// do the fit 
    TLorentzVector fit(TSVfit* fit, double met, double metPhi);
    
    /// root class definition
    ClassDef(TSVfitter,1)
  };  
}
#endif
