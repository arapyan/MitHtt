#ifndef MITHTT_NTUPLER_TMET_HH
#define MITHTT_NTUPLER_TMET_HH

#include "TObject.h"
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"

namespace mithep  {  
  class TMet : public TObject {
    public:
    TMet(){}
    ~TMet(){}
    Int_t   type;                       //0-All deposits ===> 1-Vtx related deposits
    Float_t etaMax,etaMin;
    Float_t phiMax,phiMin;
    UInt_t  nCharged,nNeut;		// Multiplicity
   
    Float_t met,      metPhi,      metm,      metpz,      sumEt;  
    Float_t chmet,    chmetPhi,    chmetm,    chmetpz,    chsumEt,   chErr;  
    Float_t chcl1met ,chcl1metPhi ,chcl1metm, chcl1metpz, chcl1sumEt;  
    Float_t chcl2met ,chcl2metPhi ,chcl2metm, chcl2metpz, chcl2sumEt;  
    Float_t neumet,   neumetPhi,   neumetm,   neumetpz,   neusumEt,  neErr;  
    Float_t neucl1met,neucl1metPhi,neucl1metm,neucl1metpz,neucl1sumEt;  
    Float_t neucl2met,neucl2metPhi,neucl2metm,neucl2metpz,neucl2sumEt;  

    Float_t emmet,      emmetPhi,      emmetm,      emmetpz,      emsumEt;  
    Float_t emchmet,    emchmetPhi,    emchmetm,    emchmetpz,    emchsumEt;  
    Float_t emchcl1met ,emchcl1metPhi ,emchcl1metm, emchcl1metpz, emchcl1sumEt;  
    Float_t emchcl2met ,emchcl2metPhi ,emchcl2metm, emchcl2metpz, emchcl2sumEt;  
    Float_t emneumet,   emneumetPhi,   emneumetm,   emneumetpz,   emneusumEt;  
    Float_t emneucl1met,emneucl1metPhi,emneucl1metm,emneucl1metpz,emneucl1sumEt;  
    Float_t emneucl2met,emneucl2metPhi,emneucl2metm,emneucl2metpz,emneucl2sumEt;  

    Float_t hadmet,      hadmetPhi,      hadmetm,      hadmetpz,      hadsumEt;  
    Float_t hadchmet,    hadchmetPhi,    hadchmetm,    hadchmetpz,    hadchsumEt;  
    Float_t hadchcl1met ,hadchcl1metPhi ,hadchcl1metm, hadchcl1metpz, hadchcl1sumEt;  
    Float_t hadchcl2met ,hadchcl2metPhi ,hadchcl2metm, hadchcl2metpz, hadchcl2sumEt;  
    Float_t hadneumet,   hadneumetPhi,   hadneumetm,   hadneumetpz,   hadneusumEt;  
    Float_t hadneucl1met,hadneucl1metPhi,hadneucl1metm,hadneucl1metpz,hadneucl1sumEt;  
    Float_t hadneucl2met,hadneucl2metPhi,hadneucl2metm,hadneucl2metpz,hadneucl2sumEt;  

    Float_t elmet,      elmetPhi,      elmetm,      elmetpz,      elsumEt   , elErr;      
    Float_t mumet,      mumetPhi,      mumetm,      mumetpz,      musumEt   , muErr;      
    Float_t gamet,      gametPhi,      gametm,      gametpz,      gasumEt   , gaErr;      
    Float_t hamet,      hametPhi,      hametm,      hametpz,      hasumEt   , haErr;      
    Float_t haneumet,   haneumetPhi,   haneumetm,   haneumetpz,   haneusumEt, haneuErr;      

    Int_t mcFlavor;			// PDG ID of matched parton flavor
    TriggerObjects hltTOMatchBits;		// bits from matching with HLT primitives
    TriggerBits hltMatchBits;		// bits from matching with HLT primitives
    
    ClassDef(TMet,5)
  };
}
#endif
