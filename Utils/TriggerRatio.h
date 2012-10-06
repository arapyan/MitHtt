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
    //CBTurnOn *R12ABarrel = new CBTurnOn(20.97643939, 1.15196354, 2.27544602, 1.01743868, 2.04391816);
    //CBTurnOn *R12BBarrel = new CBTurnOn(22.90752344, 1.32376429, 2.17813319, 1.03674051, 2.15454768);
    CBTurnOn *R12Barrel = new CBTurnOn(23.05556088,0.96047151,1.24782044,1.26042277,1.09675041);
    CBTurnOn *S12Barrel = new CBTurnOn(22.00666445,0.00036058,0.00000251,1.38456083,1.02640579);

    //CBTurnOn *R12AEndcap = new CBTurnOn(20.59874300, 1.25425435, 1.61098921, 1.00146962, 60.35067579);
    //CBTurnOn *R12BEndcap = new CBTurnOn(22.14553261, 1.19913124, 1.75642067, 1.00826962, 9.04331617);
    CBTurnOn *R12Endcap = new CBTurnOn(21.99911375,1.15806380,0.80675262,1.98765770,0.97138507);
    CBTurnOn *S12Endcap = new CBTurnOn(22.18226941,1.07762306,1.23712775,1.27324238,1.15312185);

    //CBTurnOnEta *R12B = new CBTurnOnEta(R12BBarrel, R12BEndcap, 1.479);
    CBTurnOnEta *R12  = new CBTurnOnEta(R12Barrel, R12Endcap, 1.479);
    CBTurnOnEta *S12  = new CBTurnOnEta(S12Barrel, S12Endcap,  1.479);

    //TrigEffMix *R12 = new TrigEffMix();
    //R12->add(lumi_R12A, R12A);
    //R12->add(lumi_R12B, R12B);

    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }

  TrigEffRatio *getMuonTrigEffR12()
  {
    //CBTurnOn *R12ABarrel = new CBTurnOn(15.99983195, -0.39072829, 0.28256338, 1.72861719, 0.95769408);
    //CBTurnOn *R12BBarrel = new CBTurnOn(17.21270264, 0.54997112, 1.02874912, 1.29646487, 0.96724273);
    CBTurnOn *R12Barrel =  new CBTurnOn(16.00061526,0.00737246,0.00029014,2.12854792,0.93371791);
    CBTurnOn *S12Barrel  = new CBTurnOn(16.00073094,0.00779095,0.00029834,2.13782323,0.95571348);
    
    //CBTurnOn *R12AEndcap = new CBTurnOn(18.49754887, -0.16941614, 0.26076717, 1.05494469, 1.53819978);
    //CBTurnOn *R12BEndcap = new CBTurnOn(15.98037640, 0.12062946, 0.02183977, 2.84751010, 0.83985656);
    CBTurnOn *R12Endcap  = new CBTurnOn(16.65093710,0.48774518,0.56076820,1.73768135,0.86107187);
    CBTurnOn *S12Endcap  = new CBTurnOn(17.03319591,0.73033173,1.02903291,1.46732719,0.89420534);

    //CBTurnOnEta *R12A = new CBTurnOnEta(R12ABarrel, R12AEndcap, 1.2);
    //CBTurnOnEta *R12B = new CBTurnOnEta(R12BBarrel, R12BEndcap, 1.2);
    CBTurnOnEta *R12  = new CBTurnOnEta(R12Barrel, R12Endcap, 1.2);
    CBTurnOnEta *S12  = new CBTurnOnEta(S12Barrel, S12Endcap, 1.2);

    //TrigEffMix *R12 = new TrigEffMix();
    //R12->add(lumi_R12A, R12A);
    //R12->add(lumi_R12B, R12B);

    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }

  TrigEffRatio *getTauMTrigEffR12()
  {
    //CBTurnOn *R12ABarrel = new CBTurnOn(18.52262128, 1.85879597, 3.48843815, 1.15491294, 1.02489024);
    //CBTurnOn *R12BBarrel = new CBTurnOn(17.92648563, 1.96846742, 4.46406075, 1.02023992, 1.52260575);
    CBTurnOn *R12Barrel  = new CBTurnOn(18.50940288,1.62285299,2.73232995,1.79135412,0.91481432);
    CBTurnOn *S12Barrel  = new CBTurnOn(18.80484409,0.19082817,0.19983010,1.81979820,0.93270649);
    
    //CBTurnOn *R12AEndcap = new CBTurnOn(18.90119559, 0.14025596, 0.14482632, 1.56126508, 0.81188198);
    //CBTurnOn *R12BEndcap = new CBTurnOn(18.59856420, 2.49132550, 10.99643595, 1.50651123, 0.87952970);
    CBTurnOn *R12Endcap  = new CBTurnOn(18.45678784,0.68697618,0.57008697,  3.73470825,0.84747211);
    CBTurnOn *S12Endcap  = new CBTurnOn(18.25975478,1.32745225,1.70380810,149.18410074,0.87377770);

    //CBTurnOnEta *R12A = new CBTurnOnEta(R12ABarrel, R12AEndcap, 1.5);
    //CBTurnOnEta *R12B = new CBTurnOnEta(R12BBarrel, R12BEndcap, 1.5);
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
    //CBTurnOn *R12A = new CBTurnOn(18.84658959, 0.25958704, 0.17300958, 2.43491208, 0.85872017);
    //CBTurnOn *R12B = new CBTurnOn(18.48663118, 1.63417147, 20.25695815, 138.55422224, 0.89456038);
    CBTurnOn *R12  = new CBTurnOn(18.48663118,1.63417147,20.25695815,138.55422224,0.89456038);
    CBTurnOn *S12  = new CBTurnOn(18.62733399,0.51301539, 0.38517573,  5.68099833,0.91536401);

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
