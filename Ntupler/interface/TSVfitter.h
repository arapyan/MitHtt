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
    TSVfitter(){
      fMassUnc = 0; 
      fFittedMET    = 0; fMeasMET    = 0; 
      fFittedMETPhi = 0; fMeasMETPhi = 0;
    }
    ///default destructor
    ~TSVfitter(){}
    /// do the fit 
    TLorentzVector fit(TSVfit* fit, double met, double metPhi);
    //return additional variables
    float          fittedMET()    {return fFittedMET;}
    float          measMET()      {return fMeasMET;}
    float          fittedMETPhi() {return fFittedMETPhi;}
    float          measMETPhi()   {return fMeasMETPhi;}
    float          massUnc()      {return fMassUnc;}

  protected:
    float fMassUnc;         //Uncertainty
    float fFittedMET;       //Best fit MET
    float fMeasMET;         //Original MET
    float fFittedMETPhi;    //Best fit MET phi
    float fMeasMETPhi;      //Original MET phi
    /// root class definition
    ClassDef(TSVfitter,1)
  };  
}
#endif