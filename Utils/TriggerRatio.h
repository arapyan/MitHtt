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
    CBTurnOn *R12ABarrel = new CBTurnOn(20.97643939, 1.15196354, 2.27544602, 1.01743868, 2.04391816);
    CBTurnOn *R12BBarrel = new CBTurnOn(22.90752344, 1.32376429, 2.17813319, 1.03674051, 2.15454768);
    CBTurnOn *S12Barrel  = new CBTurnOn(20.58604584, -1.89456806, 3.69311772, 1.05480046, 1.28655181);

    CBTurnOn *R12AEndcap = new CBTurnOn(20.59874300, 1.25425435, 1.61098921, 1.00146962, 60.35067579);
    CBTurnOn *R12BEndcap = new CBTurnOn(22.14553261, 1.19913124, 1.75642067, 1.00826962, 9.04331617);
    CBTurnOn *S12Endcap  = new CBTurnOn(20.15425918, 0.75449122, 1.06027513, 1.01106686, 7.01956561);

    CBTurnOnEta *R12A = new CBTurnOnEta(R12ABarrel, R12AEndcap, 1.479);
    CBTurnOnEta *R12B = new CBTurnOnEta(R12BBarrel, R12BEndcap, 1.479);
    CBTurnOnEta *S12  = new CBTurnOnEta(S12Barrel,  S12Endcap,  1.479);

    TrigEffMix *R12 = new TrigEffMix();
    R12->add(lumi_R12A, R12A);
    R12->add(lumi_R12B, R12B);

    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }

  TrigEffRatio *getMuonTrigEffR12()
  {
    CBTurnOn *R12ABarrel = new CBTurnOn(15.99983195, -0.39072829, 0.28256338, 1.72861719, 0.95769408);
    CBTurnOn *R12BBarrel = new CBTurnOn(17.21270264, 0.54997112, 1.02874912, 1.29646487, 0.96724273);
    CBTurnOn *S12Barrel  = new CBTurnOn(16.99389526, -0.04080190, 0.00794730, 1.60377906, 0.99626161);

    CBTurnOn *R12AEndcap = new CBTurnOn(18.49754887, -0.16941614, 0.26076717, 1.05494469, 1.53819978);
    CBTurnOn *R12BEndcap = new CBTurnOn(15.98037640, 0.12062946, 0.02183977, 2.84751010, 0.83985656);
    CBTurnOn *S12Endcap  = new CBTurnOn(16.99065795, -0.11993730, 0.01384991, 2.38867304, 0.86552275);

    CBTurnOnEta *R12A = new CBTurnOnEta(R12ABarrel, R12AEndcap, 1.2);
    CBTurnOnEta *R12B = new CBTurnOnEta(R12BBarrel, R12BEndcap, 1.2);
    CBTurnOnEta *S12  = new CBTurnOnEta(S12Barrel,  S12Endcap,  1.2);

    TrigEffMix *R12 = new TrigEffMix();
    R12->add(lumi_R12A, R12A);
    R12->add(lumi_R12B, R12B);

    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }

  TrigEffRatio *getTauMTrigEffR12()
  {
    CBTurnOn *R12ABarrel = new CBTurnOn(18.52262128, 1.85879597, 3.48843815, 1.15491294, 1.02489024);
    CBTurnOn *R12BBarrel = new CBTurnOn(17.92648563, 1.96846742, 4.46406075, 1.02023992, 1.52260575);
    CBTurnOn *S12Barrel  = new CBTurnOn(18.86257072, 0.25680380, 0.16916101, 2.42931257, 0.89590264);

    CBTurnOn *R12AEndcap = new CBTurnOn(18.90119559, 0.14025596, 0.14482632, 1.56126508, 0.81188198);
    CBTurnOn *R12BEndcap = new CBTurnOn(18.59856420, 2.49132550, 10.99643595, 1.50651123, 0.87952970);
    CBTurnOn *S12Endcap  = new CBTurnOn(18.74764561, 1.82036845, 701.46994969, 101.57913480, 0.82547043);

    CBTurnOnEta *R12A = new CBTurnOnEta(R12ABarrel, R12AEndcap, 1.5);
    CBTurnOnEta *R12B = new CBTurnOnEta(R12BBarrel, R12BEndcap, 1.5);
    CBTurnOnEta *S12  = new CBTurnOnEta(S12Barrel,  S12Endcap,  1.5);

    TrigEffMix *R12 = new TrigEffMix();
    R12->add(lumi_R12A, R12A);
    R12->add(lumi_R12B, R12B);

    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }

  TrigEffRatio *getTauETrigEffR12()
  {
    CBTurnOn *R12A = new CBTurnOn(18.84658959, 0.25958704, 0.17300958, 2.43491208, 0.85872017);
    CBTurnOn *R12B = new CBTurnOn(18.48663118, 1.63417147, 20.25695815, 138.55422224, 0.89456038);
    CBTurnOn *S12  = new CBTurnOn(18.77448606, 0.45765507, 0.26077509, 13.43372485, 0.88037836);

    TrigEffMix *R12 = new TrigEffMix();
    R12->add(lumi_R12A, R12A);
    R12->add(lumi_R12B, R12B);

    TrigEffRatio *ratio = new TrigEffRatio();
    ratio->setNumerator(R12);
    ratio->setDenominator(S12);

    return ratio;
  }
}
#endif
