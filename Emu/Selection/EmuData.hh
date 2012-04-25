#ifndef EMUDATA_HH
#define EMUDATA_HH

struct EmuData
{
  UInt_t  runNum;                     // run number in data
  UInt_t  evtNum;                     // event number in data
  UInt_t  lumiSec;                    // lumi section      
  UInt_t  nPV;                        // number of reconstructed primary vertices (with some requirements)
  UInt_t  njets, nbjets;              // number of jets, b-jets
  Float_t vpt, vphi;                  // generator boson pT/phi
  Float_t rawmet, rawmetphi;          // raw MET
  Float_t met, metphi;                // MET
  Float_t mass, dphi, mt, pt, phi;    // dilepton mass, dPhi
  Float_t pmet, pvis;                 // projected met, projected dilepton pt
  Float_t eleiso, muiso;              // lepton isolation
  Float_t eled0, eled0sig;            // electron d0
  Float_t eleip3d, eleip3dsig;        // electron 3D impact parameter
  Float_t mud0, mud0sig;              // muon d0
  Float_t muip3d, muip3dsig;          // muon 3D impact parameter
  Float_t lpt1, leta1, lphi1;         // lepton 1 kinematics
  Float_t lpt2, leta2, lphi2;         // lepton 2 kinematics
  Float_t jpt1, jeta1, jphi1;         // leading jet kinematics
  Float_t jpt2, jeta2, jphi2;         // sub-leading jet kinematics
  Float_t bjpt, bjeta, bjphi;         // leading b-jet kinematics
  Float_t mjj;                        // dijet mass
  Float_t svfmass, svfmassunc;        // svfit mass and mass uncertainty
  Float_t genlpt1, genleta1, genlphi1;// generator lepton 1 kinematics
  Float_t genlpt2, genleta2, genlphi2;// generator lepton 2 kinematics
  Float_t weight;                     // event weight per 1/pb
  Int_t   state;                      // dilepton final state: mm, ee, em, me
};

#endif
