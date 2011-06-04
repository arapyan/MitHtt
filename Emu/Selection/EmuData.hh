#ifndef EMUDATA_HH
#define EMUDATA_HH

struct EmuData
{
  UInt_t  runNum;                   // run number in data
  UInt_t  evtNum;                   // event number in data
  UInt_t  lumiSec;                  // lumi section      
  UInt_t  nPV;                      // number of reconstructed primary vertices (with some requirements)                                          
  UInt_t  njets, nbjets;            // number of jets, b-jets
  Float_t met, metphi;              // MET
  Float_t mass, dphi, mt, pt, phi;  // dilepton mass, dPhi
  Float_t pmet, pvis;               // projected met, projected dilepton pt
  Float_t lpt1, leta1, lphi1;       // lepton 1 kinematics
  Float_t lpt2, leta2, lphi2;       // lepton 2 kinematics
  Float_t jpt1, jeta1, jphi1;       // leading jet kinematics
  Float_t jpt2, jeta2, jphi2;       // sub-leading jet kinematics
  Float_t bjpt, bjeta, bjphi;       // leading b-jet kinematics
  Float_t mjj;                      // dijet mass
  Float_t weight;                   // event weight per 1/pb
  Int_t   state;                    // dilepton final state: mm, ee, em, me
};

//"runNum/i:evtNum:lumiSec:nPV:njets:nbjets:met/F:metphi:mass:dphi:mt:pt:phi:pmet:pvis:
//lpt1:leta1:lphi1:lpt2:leta2:lphi2:
//jpt1:jeta1:jphi1:jpt2:jeta2:jphi2:
//mjj:weight:state/I"

#endif
