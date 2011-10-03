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
  Float_t mass, scaledmass, dphi, mt, pt, phi;  // dilepton mass, dPhi
  Float_t pmet, pvis;               // projected met, projected dilepton pt
  Float_t lpt1, leta1, lphi1;       // lepton 1 kinematics
  Float_t lpt2, leta2, lphi2;       // lepton 2 kinematics
  Float_t jpt1, jeta1, jphi1;       // leading jet kinematics
  Float_t jpt2, jeta2, jphi2;       // sub-leading jet kinematics
  Float_t bjpt, bjeta, bjphi;    // leading b-jet kinematics
//  Float_t bjpt2, bjeta2, bjphi2;    // sub-leading b-jet kinematics
  Float_t mjj;                 // dijet mass, dibjet mass
  Float_t svfmass;                  // svfit mass
  Float_t svflpt1, svfleta1, svflphi1; // lepton 1 kinematics
  Float_t svflpt2, svfleta2, svflphi2; // lepton 2 kinematics
  Float_t weight;                   // event weight per 1/pb
  Int_t   state;                    // dilepton final state: mm, ee, em, me
};

#endif
