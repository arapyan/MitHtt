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
    double integrateMarkov(TSVfit* fit, double met, double metPhi, int id);
    double integrate(TSVfit* fit, double met, double metPhi, int id);
    //return additional variables
    float          fittedMET()    {return fFittedMET;}
    float          measMET()      {return fMeasMET;}
    float          fittedMETPhi() {return fFittedMETPhi;}
    float          measMETPhi()   {return fMeasMETPhi;}
    int            status()       {return fstatus;}

    //markov chain variables
    float      massUnc()      {return fMassUnc;}
    float      GetPt()   {return fPtMarkov;}
    float      GetPhi()  {return fPhiMarkov;}
    float      GetPtUnc() {return fPtMarkovUnc;}
    float      GetPhiUnc() {return fPhiMarkovUnc;}
    float      GetEta() {return fEtaMarkov;}
    float      GetEtaUnc() {return fEtaMarkovUnc;}

  protected:
    //Markov Chain DiTau system
    float fMassUnc;         //Uncertainty
    float fPtMarkov;   
    float fEtaMarkov;
    float fPhiMarkov;
    float fPhiMarkovUnc;
    float fEtaMarkovUnc;
    float fPtMarkovUnc;
    //**************************************************
    int fstatus;             //valid solution?
    float fFittedMET;       //Best fit MET
    float fMeasMET;         //Original MET
    float fFittedMETPhi;    //Best fit MET phi
    float fMeasMETPhi;      //Original MET phi
    /// root class definition
    ClassDef(TSVfitter,3)
  };  
}
#endif
