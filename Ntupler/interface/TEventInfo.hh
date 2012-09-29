#ifndef MITHTT_NTUPLER_TEVENTINFO_HH
#define MITHTT_NTUPLER_TEVENTINFO_HH

#include "TObject.h"
#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"

/**
   \class TEventInfo TEventInfo.h MitHtt/Ntupler/include/TEventInfo.h

   \brief Description: Common event information

   Common event information, filled in HttNtuple::fillCommon. Pileup information exists for MC events 
   only. It keeps the information of the number of simulated pileup events during the actual bunch 
   crossing (in time PU) and of the simulated number of event before and after the actual bunch 
   crossing (out of time PU). The beamspot is taken from Bambu w/o further processing. The reconstructed
   primary vertex is chosen to be the one with the highest sum pt**2, which fullfills a minimal set 
   of selection criteria. If no primary vertex is found, which fullfills the selection criteria the 
   selection falls back to the primary vertex with the highest sum pt**2, which does not fullfill the 
   selection criteria, but hasGoodPV will be set to false. For trkMet ony tracks from particle flow 
   candidates are considered, which have a DzCorrected to the selected primary vertex smaller than 0.1
   (hard coded in fillCommon()).
*/

namespace mithep 
{
  class TEventInfo : public TObject
  {
  public:
    /// default constructor
    TEventInfo(){}
    /// default destructor
    ~TEventInfo(){}
    
    /// run number
    unsigned int  runNum;
    /// event number
    unsigned int  evtNum;
    /// lumi section
    unsigned int  lumiSec;
    /// number of in-time pileup event per event (for MC only)
    unsigned int  nPU;  
    /// number of out of time pileup events in following event (for MC only)
    unsigned int  nPUPlus;
    /// number of out of time pileup events previous event (for MC only)
    unsigned int  nPUMinus;
    /// poisson mean from which npu was thrown
    float nPUTrue, nPUPlusTrue, nPUMinusTrue;
    /// HLT trigger bits that fired for this event
    TriggerBits triggerBits;
    /// coordiantes of the reconstructed primary vertex with highest sum pt squared
    float pvx, pvy, pvz;
    /// coordinates of the beam spot
    float bsx, bsy, bsz;
    /// particle flow MET (from Bambu w/o any further processing)
    float pfMET, pfMETphi, pfSumET;
    /// track MET 
    float trkMET, trkMETphi, trkSumET;
    /// MVA MET
    float mvaMET, mvaMETphi, mvaCov00, mvaCov11, mvaCov10, mvaCov01;
    /// mean energy density in the event (due to pileup, from L1Fastjet algorithm)
    float rho, rhoHighEta;
    /// does a reconstructed primary vertex exist that fullfills minimal selection criteria
    bool hasGoodPV;
    /// embedding weight (for embedding sample only) 
    float embWeight;

    ClassDef(TEventInfo, 2)
  };
}
#endif
