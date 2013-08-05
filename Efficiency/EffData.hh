#ifndef EFF_DATA_HH
#define EFF_DATA_HH

struct EffData
{
  Float_t mass, pt, eta, phi, weight;
  Int_t q;
  UInt_t npv;
  Float_t npu;
  UInt_t pass;
  UInt_t runNum, lumiSec, evtNum;
  Float_t rho;
};

// "mass/F:pt:eta:phi:weight:q/I:npv/i:npu/F:pass/i:runNum:lumiSec:evtNum:rho/F"

#endif
