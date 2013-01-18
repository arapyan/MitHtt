#ifndef TRIGGERATIO_H
#define TRIGGERATIO_H

#include "TriggerEfficiency.h"

namespace mithep
{
  //
  // 2012 Trigger Efficiency Scale Factors
  //
  const double lumi_R12A = 0.14;
  const double lumi_R12B = 0.86;

  TrigEffRatio *getElectronTrigEffR12()
  {

    CBTurnOn *R12Barrel  = new CBTurnOn(22.9041, 1.04728 , 1.38544,  1.22576, 1.13019);
    CBTurnOn *S12Barrel  = new CBTurnOn(21.7243, 0.619015, 0.739301, 1.34903, 1.02594);
    
    CBTurnOn *R12Endcap  = new CBTurnOn(21.9941, 1.43419, 1.01152, 2.28622, 0.939872);
    CBTurnOn *S12Endcap  = new CBTurnOn(22.1217, 1.34054, 1.8885,  1.01855, 4.7241);

    CBTurnOnEta *R12  = new CBTurnOnEta(R12Barrel, R12Endcap, 1.5);
    CBTurnOnEta *S12  = new CBTurnOnEta(S12Barrel, S12Endcap, 1.5);

    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }

  TrigEffRatio *getMuonTrigEffR12()
  {

    CBTurnOn *R12EM  = new CBTurnOn(15.9825, 7.90724e-05, 5.49275e-08, 1.6403, 0.858285);
    CBTurnOn *S12EM  = new CBTurnOn(16.0051, 2.45144e-05, 4.3335e-09,  1.66134, 0.87045);
    
    CBTurnOn *R12TM  = new CBTurnOn(17.3283, 0.707103, 1.2047, 1.3732, 0.900519);
    CBTurnOn *S12TM  = new CBTurnOn(17.3135, 0.747636, 1.21803, 1.40611, 0.934983);

    CBTurnOn *R12BM  = new CBTurnOn(15.9828, 0.0412999, 0.0177441, 1.66934, 0.970097);
    CBTurnOn *S12BM  = new CBTurnOn(15.9556, 0.0236127, 0.00589832, 1.75409, 0.981338);

    CBTurnOn *R12BP  = new CBTurnOn(15.9802, 0.0548775, 0.020313, 1.79791, 0.968398);
    CBTurnOn *S12BP  = new CBTurnOn(15.9289, 0.0271317, 0.00448573, 1.92101, 0.978625);

    CBTurnOn *R12TP  = new CBTurnOn(16.8396, 0.458636, 0.633185, 1.5706, 0.8848);
    CBTurnOn *S12TP  = new CBTurnOn(16.5678, 0.328333, 0.354533, 1.67085, 0.916992);

    CBTurnOn *R12EP  = new CBTurnOn(15.9987, 8.94398e-05, 5.18549e-08, 1.8342, 0.854625);
    CBTurnOn *S12EP  = new CBTurnOn(15.997,  7.90069e-05, 4.40036e-08, 1.66272, 0.884502);

    CBTurnOn6Eta *R12  = new CBTurnOn6Eta(R12EM,R12TM,R12BM,R12BP,R12TP,R12EP,-1.2,-0.8,0.,0.8,1.2);
    CBTurnOn6Eta *S12  = new CBTurnOn6Eta(S12EM,S12TM,S12BM,S12BP,S12TP,S12EP,-1.2,-0.8,0.,0.8,1.2);


    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }

  TrigEffRatio *getTauMTrigEffR12()
  {
    //CBTurnOn *R12ABarrel = new CBTurnOn(18.52262128, 1.85879597, 3.48843815, 1.15491294, 1.02489024);
    //CBTurnOn *R12BBarrel = new CBTurnOn(17.92648563, 1.96846742, 4.46406075, 1.02023992, 1.52260575);
    //CBTurnOn *R12Barrel  = new CBTurnOn(18.50940288,1.62285299,2.73232995,1.79135412,0.91481432);
    //CBTurnOn *S12Barrel  = new CBTurnOn(18.80484409,0.19082817,0.19983010,1.81979820,0.93270649);
    //CBTurnOn *R12AEndcap = new CBTurnOn(18.90119559, 0.14025596, 0.14482632, 1.56126508, 0.81188198);
    //CBTurnOn *R12BEndcap = new CBTurnOn(18.59856420, 2.49132550, 10.99643595, 1.50651123, 0.87952970);

    //CBTurnOnEta *R12A = new CBTurnOnEta(R12ABarrel, R12AEndcap, 1.5);
    //CBTurnOnEta *R12B = new CBTurnOnEta(R12BBarrel, R12BEndcap, 1.5);

    CBTurnOn *R12Barrel  = new CBTurnOn(18.52036251, 1.47760312, 2.53574445, 1.71202550, 0.93019930 );
    CBTurnOn *S12Barrel  = new CBTurnOn(18.88740627, 0.10718873, 0.12277723, 1.60581265, 0.95041892 );
    
    CBTurnOn *R12Endcap  = new CBTurnOn(18.41225333, 0.76598912, 0.60544260, 5.38350881, 0.85870108 );
    CBTurnOn *S12Endcap  = new CBTurnOn(18.30439676, 1.44360240, 3.79358997, 1.07560564, 0.93103925 );

    CBTurnOnEta *R12  = new CBTurnOnEta(R12Barrel, R12Endcap, 1.5);
    CBTurnOnEta *S12  = new CBTurnOnEta(S12Barrel, S12Endcap, 1.5);

    //TrigEffMix *R12 = new TrigEffMix();
    //R12->add(lumi_R12A, R12A);
    //R12->add(lumi_R12B, R12B);
    //R12->add(1.0, R12);
    
    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }

  TrigEffRatio *getTauETrigEffR12()
  {


    CBTurnOn *R12Barrel  = new CBTurnOn(18.686211, 1.993524, 3.202713, 3.612693, 0.871640 );
    CBTurnOn *S12Barrel  = new CBTurnOn(18.431118, 1.572877, 3.301699, 4.760769, 0.899620 );
    
    CBTurnOn *R12Endcap  = new CBTurnOn(18.472954, 1.606388, 3.468975, 55.629620, 0.828977 );
    CBTurnOn *S12Endcap  = new CBTurnOn(18.257217, 1.632443, 9.283116, 40.219585, 0.858643 );

    CBTurnOnEta *R12  = new CBTurnOnEta(R12Barrel, R12Endcap, 1.5);
    CBTurnOnEta *S12  = new CBTurnOnEta(S12Barrel, S12Endcap, 1.5);

    //TrigEffMix *R12 = new TrigEffMix();
    //R12->add(lumi_R12A, R12A);
    //R12->add(lumi_R12B, R12B);
    //R12->add(1.0, R12);

    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator  (R12);
    ratio->setDenominator(S12);

    return ratio;
  }
}
#endif
