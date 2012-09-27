#ifndef MITHTT_NTUPLER_HIGGSANADEFS_HH 
#define MITHTT_NTUPLER_HIGGSANADEFS_HH

#include "bitset"
#include "Rtypes.h"

/**
   \enum    ESampleType
   \brief   enumeration of special sample types which cause different filling of generator information in the Ntupler
*/
namespace ESampleType {
  enum {
    kH     = 1,
    kZ     = 2,
    kW     = 3,
    kVV    = 4,
    kHWW   = 5,
    kHZZ   = 6,
    kVttH  = 7,
    kEmbed = 8
  };
}

/**
   \enum    EGenType
   \brief   enumeration of usual pdgIds for generator particles
*/
namespace EGenType {
  enum {
    kElectron    = 11,
    kENeutrino   = 12,
    kMuon        = 13,
    kMNeutrino   = 14,
    kTau         = 15,
    kTNeutrino   = 16,
    kTauMuon     = 17,
    kTauElectron = 18,
    kTauHadr     = 19,
    kPhoton      = 22,
    kZ           = 23,
    kW           = 24,
    kHiggs       = 25,
    kttH         = 26,
    kUp          = 1,
    kDown        = 2,
    kStrange     = 3,
    kCharm       = 4,
    kBottom      = 5,
    kTop         = 6
  };
}

/**
   \enum    EMuType
   \brief   enumeration of muon reconstruction types
*/
enum EMuType 
{ 
  kGlobal     = 1, 
  kTracker    = 2, 
  kStandalone = 4
};

/**
   \enum    EQualityBit
   \brief   enumeration of muon quality bits (taken from DataFormats/MuonReco/interface/MuonSelectors.h)
*/
enum EQualityBit
{ 
  kAll  			          = 0x000001,  // dummy option - always true
  kAllGlobalMuons		          = 0x000002,  // checks isGlobalMuon flag
  kAllStandAloneMuons		          = 0x000004,  // checks isStandAloneMuon flag
  kAllTrackerMuons		          = 0x000008,  // checks isTrackerMuon flag
  kTrackerMuonArbitrated	          = 0x000010,  // resolve ambiguity of sharing segments
  kAllArbitrated		          = 0x000020,  // all muons with the tracker muon arbitrated
  kGlobalMuonPromptTight	          = 0x000040,  // global muons with tighter fit requirements
  kTMLastStationLoose		          = 0x000080,  // penetration depth loose selector
  kTMLastStationTight		          = 0x000100,  // penetration depth tight selector
  kTM2DCompatibilityLoose	          = 0x000200,  // likelihood based loose selector
  kTM2DCompatibilityTight	          = 0x000400,  // likelihood based tight selector
  kTMOneStationLoose		          = 0x000800,  // require one well matched segment
  kTMOneStationTight		          = 0x001000,  // require one well matched segment
  kTMLastStationOptimizedLowPtLoose       = 0x002000,  // combination of TMLastStation and TMOneStation
  kTMLastStationOptimizedLowPtTight       = 0x004000,  // combination of TMLastStation and TMOneStation
  kGMTkChiCompatibility 	          = 0x008000,  // require tk stub have good chi2 relative to glb track
  kGMStaChiCompatibility	          = 0x010000,  // require sta stub have good chi2 compatibility relative to glb track
  kGMTkKinkTight		          = 0x020000,  // require a small kink value in the tracker stub
  kTMLastStationAngLoose	          = 0x040000,  // TMLastStationLoose with additional angular cuts
  kTMLastStationAngTight	          = 0x080000,  // TMLastStationTight with additional angular cuts
  kTMOneStationAngLoose 	          = 0x100000,  // TMOneStationLoose with additional angular cuts
  kTMOneStationAngTight 	          = 0x200000,  // TMOneStationTight with additional angular cuts
  kTMLastStationOptimizedBarrelLowPtLoose = 0x400000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
  kTMLastStationOptimizedBarrelLowPtTight = 0x800000   // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
}; 

const unsigned int kNTrigBit = 128;
typedef std::bitset<kNTrigBit> TriggerBits;

/**
   \enum    ETriggerBit
   \brief   enumeration of trigger bits in the TriggerBits bitset
*/
enum ETriggerBit {    
  /* MuEG                       --------------------------------------------------------------------------------------------- */
  kHLT_Mu11_Ele8                                                                        = 0,  // MC
  kHLT_Mu17_Ele8_CaloIdL                                                                = 0,  // data
  kHLT_Mu8_Ele8                                                                         = 1,  // MC
  kHLT_Mu8_Ele17_CaloIdL                                                                = 1,  // data
  kHLT_Mu8_Photon20_CaloIdVT_IsoT                                                       = 2,  // data
  kHLT_Mu15_Photon20_CaloIdL                                                            = 3,  // data
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL                                                      = 4,
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL                                                      = 5,
  /* DoubleMu                   --------------------------------------------------------------------------------------------- */
  kHLT_DoubleMu5                                                                        = 6,  // MC
  kHLT_DoubleMu7                                                                        = 6,  // data
  kHLT_Mu13_Mu8                                                                         = 7,  // data
  kHLT_Mu17_Mu8                                                                         = 8,  // data
  kHLT_Mu5_Jet50U                                                                       = 9,  // MC
  kHLT_Mu8_Jet40                                                                        = 9,  // data
  kHLT_Mu17_TkMu8                                                                       = 10,
  /* SingleMu                   --------------------------------------------------------------------------------------------- */
  kHLT_Mu8                                                                              = 11, // data
  kHLT_Mu9                                                                              = 11, // MC
  kHLT_Mu12                                                                             = 12, // data
  kHLT_Mu11                                                                             = 13, // MC
  kHLT_Mu15                                                                             = 13, // MC, data
  kHLT_Mu24                                                                             = 14, // data
  kHLT_Mu21                                                                             = 15, // MC
  kHLT_Mu30                                                                             = 16, // data
  kHLT_IsoMu17                                                                          = 16, // MC, data
  kHLT_IsoMu24                                                                          = 17, // data
  kHLT_IsoMu15                                                                          = 18, // data
  /* DoubleElectron             ----------------------------------------------------------------------------------------------- */
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL                                   = 19, // data
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 20, // data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30                              = 21, // data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30                             = 22, // data
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                                                     = 23, // data
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17                                       = 24, // data
  /* SingleElectron             ----------------------------------------------------------------------------------------------- */
  kHLT_Ele8                                                                             = 25, // data
  kHLT_Ele8_CaloIdL_TrkIdVL                                                             = 26, // data
  kHLT_Ele8_CaloIdL_CaloIsoVL                                                           = 27, // data
  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL                                          = 28, // data
  kHLT_Ele17_CaloIdL_CaloIsoVL                                                          = 29, // data
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40                                                     = 30, // data
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL                                    = 31, // data
  kHLT_Ele17_SW_L1R                                                                     = 32, // MC
  kHLT_Ele22_SW_L1R                                                                     = 33, // MC
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 34, // data
  kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 35, // data
  kHLT_Ele45_CaloIdVT_TrkIdT                                                            = 36, // data
  kHLT_Ele52_CaloIdVT_TrkIdT                                                            = 37, // data
  kHLT_Ele65_CaloIdVT_TrkIdT                                                            = 38, // data
  kHLT_Ele80_CaloIdVT_TrkIdT                                                            = 39, // data
  /* Old Triggers               ---------------------------------------------------------------------------------------------- */
  kHLT_Photon10_L1R			                                                = 40,
  kHLT_Photon15_Cleaned_L1R		                                                = 41,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 42,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 43,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 44,
  kHLT_Ele22_SW_TighterCaloIdIsol_L1R	                                                = 45,
 
  /* JetMET                   ------------------------------------------------------------------------------------------------ */
  kHLT_CentralJet80_MET100		                                                = 46,
  kHLT_CentralJet80_MET160		                                                = 47,
  kHLT_CentralJet80_MET80		                                                = 48,
  kHLT_DiCentralJet20_BTagIP_MET65	                                                = 49,
  kHLT_DiCentralJet20_MET80     	                                                = 50,
  kHLT_MET100_HBHENoiseFiltered 	                                                = 51,
  kHLT_MET100   			                                                = 52,
  kHLT_MET120_HBHENoiseFiltered  	                                                = 53,
  kHLT_MET120   			                                                = 54,
  /* 2012 Triggers               ---------------------------------------------------------------------------------------------- */
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL                                     = 55,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL                                     = 56,
  kHLT_Ele27_WP80                                                                       = 57,
  kHLT_Mu17                                                                             = 58,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50                             = 59,
  kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50                              = 60,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50                                = 61,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30                                   = 62,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL                                         = 63,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30                                    = 64,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL                                          = 65,

  //Tau+Tau 2012
  kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30 = 66,
  //Tau+Tau 2011
  kHLT_DoubleIsoPFTau35_Trk5_eta2p1             = 67,
  kHLT_DoubleIsoPFTau25_Trk5_eta2p1             = 68,
  kHLT_DoubleIsoPFTau20_Trk5                    = 69,

  // Mu+Tau
  kHLT_Mu15_LooseIsoPFTau15           = 70,
  kHLT_Mu15_LooseIsoPFTau20           = 71,
  kHLT_Mu18_eta2p1_LooseIsoPFTau20    = 72,
  kHLT_IsoMu12_LooseIsoPFTau10        = 73,
  kHLT_IsoMu15_LooseIsoPFTau15        = 74,
  kHLT_IsoMu15_LooseIsoPFTau20        = 75,
  kHLT_IsoMu15_TightIsoPFTau20        = 76,
  kHLT_IsoMu15_eta2p1_LooseIsoPFTau20 = 77,
  kHLT_IsoMu15_eta2p1_MediumIsoPFTau20 = 78,
  kHLT_IsoMu15_eta2p1_TightIsoPFTau20 = 79,
  kHLT_IsoMu18_eta2p1_LooseIsoPFTau20 = 80,
  // E+Tau
  kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15 = 81,
  kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20 = 82,
  kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20 = 83,
  kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20 = 84,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15 = 85,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20 = 86,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20 = 87,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20 = 88,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20 = 89,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20 = 90,
  kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20 = 91,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20 = 92,
  kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20 = 93, 
  kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20 = 94,
  kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20 = 95,
  
  //additional MuTau 2012
  kHLT_IsoMu17_eta2p1_LooseIsoPFTau20 = 96,
  kHLT_Mu17_eta2p1_LooseIsoPFTau20   = 97,
   
  //misc
  kHLT_IsoMu15_L1ETM20 = 98,
  kHLT_Jet150          = 99,
  kHLT_Jet190          = 100,
  kHLT_Jet240          = 101,

  //additional 2012 Tau+Tau
  kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30 = 102,
  kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30 = 103
};

const unsigned int kNTrigObj = 256;
typedef std::bitset<kNTrigObj> TriggerObjects;

/**
   \enum    ETriggerObject
   \brief   enumeration of trigger bits in the TriggerObjects bitset
*/
enum ETriggerObject {
  /* MuEG                      ------------------------------------------------------------------------------------------------- */
  kHLT_Mu17_Ele8_CaloIdL_MuObj                                                                  =   0,
  kHLT_Mu17_Ele8_CaloIdL_EGObj                                                                  =   1,  
  kHLT_Mu8_Ele17_CaloIdL_MuObj                                                                  =   2,
  kHLT_Mu8_Ele17_CaloIdL_EGObj                                                                  =   3,
  kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj                                                         =   4,
  kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj                                                         =   5,
  kHLT_Mu15_Photon20_CaloIdL_MuObj                                                              =   6,
  kHLT_Mu15_Photon20_CaloIdL_EGObj                                                              =   7,
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj                                                        =   8,
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj                                                        =   9,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj                                                        =  10,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj                                                        =  11,
  /* DoubleMuon                ------------------------------------------------------------------------------------------------- */
  kHLT_DoubleMu7_MuObj                                                                          =  12,
  kHLT_Mu13_Mu8_Mu1Obj                                                                          =  13,
  kHLT_Mu13_Mu8_Mu2Obj                                                                          =  14,
  kHLT_Mu17_Mu8_Mu1Obj                                                                          =  15,
  kHLT_Mu17_Mu8_Mu2Obj                                                                          =  16,
  kHLT_Mu8_Jet40_MuObj                                                                          =  17,
  kHLT_Mu8_Jet40_JetObj                                                                         =  18,
  kHLT_Mu17_TkMu8_Mu1Obj                                                                        =  19,
  kHLT_Mu17_TkMu8_Mu2Obj                                                                        =  20,
  /* SingleMuon                ------------------------------------------------------------------------------------------------- */
  kHLT_Mu8_MuObj                                                                                =  21,
  kHLT_Mu9_MuObj                                                                                =  22,
  kHLT_Mu12_MuObj                                                                               =  23,
  kHLT_Mu15_MuObj                                                                               =  24,
  kHLT_Mu24_MuObj                                                                               =  25,
  kHLT_Mu30_MuObj                                                                               =  26,
  kHLT_IsoMu15_MuObj                                                                            =  27,
  kHLT_IsoMu17_MuObj                                                                            =  28,
  kHLT_IsoMu24_MuObj                                                                            =  29,
  /* DoubleElectron            ------------------------------------------------------------------------------------------------- */
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj                                   =  30,
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj                                   =  31,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj =  32,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj =  33,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj                               =  34,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj                                =  35,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj                             =  36,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj                             =  37,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj                                                      =  38,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj                                                       =  39,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj                                        =  40,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj                                         =  41,
  kHLT_Ele8_EleObj                                                                              =  42,
  kHLT_Ele8_CaloIdL_TrkIdVL_EleObj                                                              =  43,
  kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj                                                            =  44,
  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj                                           =  45,
  kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj                                                           =  46,
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj                                                      =  47,
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj                                                      =  48,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj                                     =  49,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj                                     =  50,
  /* SingleElectron             ------------------------------------------------------------------------------------------------- */  
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  51,
  kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  52,
  kHLT_Ele45_CaloIdVT_TrkIdT_EleObj                                                             =  53, // data
  kHLT_Ele52_CaloIdVT_TrkIdT_EleObj                                                             =  54, // data
  kHLT_Ele65_CaloIdVT_TrkIdT_EleObj                                                             =  55, // data
  kHLT_Ele80_CaloIdVT_TrkIdT_EleObj                                                             =  56, // data
  kHLT_Photon10_L1R_EleObj			                                                =  57,
  kHLT_Photon15_Cleaned_L1R_EleObj		                                                =  58,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  59,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  60,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  61,
  kHLT_Ele22_SW_TighterCaloIdIsol_L1R_EleObj	                                                =  62,
  /* JetMET                    ------------------------------------------------------------------------------------------------- */
  kHLT_Jet30_JetObj       = 63,
  kHLT_Jet150_JetObj                                                                            =  64,
  kHLT_Jet190_JetObj                                                                            =  65,
  kHLT_Jet240_JetObj                                                                            =  66,
  kHLT_Jet300_JetObj                                                                            =  67,
  kHLT_DiJetAve30_JetObj                                                                        =  68,
  kHLT_DiJetAve240_JetObj                                                                       =  69,
  kHLT_DiJetAve300_JetObj                                                                       =  70,
  /* Photon                    ------------------------------------------------------------------------------------------------- */
  kHLT_Photon20_CaloIdVL_IsoL_PhoObj                                                            =  71,
  kHLT_Photon30_CaloIdVL_IsoL_PhoObj                                                            = 72,
  kHLT_Photon50_CaloIdVL_IsoL_PhoObj                                                            = 73,
  kHLT_Photon75_CaloIdVL_PhoObj                                                                 = 74,
  /* 2012                      ------------------------------------------------------------------------------------------------- */
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj                                       = 75,
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj                                       = 76,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj                                       = 77,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj                                       = 78,
  kHLT_Ele27_WP80_EleObj                                                                        = 79,
  kHLT_Mu17_MuObj                                                                               = 80,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj                             = 81,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele2Obj                             = 82,
  kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj                               = 83,
  kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj                                = 84,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_EleObj                                 = 85,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_SCObj                                  = 86,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj                                    = 87,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj                                          = 88,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj                                     = 89,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj                                           = 90,
  // Mu + MET
  kHLT_IsoMu15_L1ETM20_MuObj = 91,
  
   // TauPlusX
  kHLT_IsoMu12_LooseIsoPFTau10_MuObj  = 92,
  kHLT_IsoMu12_LooseIsoPFTau10_TauObj = 93,
  kHLT_Mu15_LooseIsoPFTau15_MuObj     = 94,
  kHLT_Mu15_LooseIsoPFTau15_TauObj    = 95,
  kHLT_Mu15_LooseIsoPFTau20_MuObj     = 96,
  kHLT_Mu15_LooseIsoPFTau20_TauObj    = 97,
  kHLT_Mu18_eta2p1_LooseIsoPFTau20_MuObj = 98,
  kHLT_Mu18_eta2p1_LooseIsoPFTau20_TauObj = 99,
  kHLT_IsoMu15_LooseIsoPFTau15_MuObj  = 100,
  kHLT_IsoMu15_LooseIsoPFTau15_TauObj = 101,
  kHLT_IsoMu15_LooseIsoPFTau20_MuObj  = 102,
  kHLT_IsoMu15_LooseIsoPFTau20_TauObj = 103,
  kHLT_IsoMu15_TightIsoPFTau20_MuObj  = 104,
  kHLT_IsoMu15_TightIsoPFTau20_TauObj = 105,
  kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_MuObj = 106,
  kHLT_IsoMu15_eta2p1_LooseIsoPFTau20_TauObj = 107,
  kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_MuObj = 108,
  kHLT_IsoMu15_eta2p1_MediumIsoPFTau20_TauObj = 109,
  kHLT_IsoMu15_eta2p1_TightIsoPFTau20_MuObj = 110,
  kHLT_IsoMu15_eta2p1_TightIsoPFTau20_TauObj = 111,
  kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_MuObj = 112,
  kHLT_IsoMu18_eta2p1_LooseIsoPFTau20_TauObj = 113,
  kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_EleObj = 114,
  kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau15_TauObj = 115,
  kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj = 116,
  kHLT_Ele15_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj = 117,
  kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_EleObj = 118,
  kHLT_Ele15_CaloIdVT_TrkIdT_TightIsoPFTau20_TauObj = 119,
  kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_EleObj = 120,
  kHLT_Ele18_CaloIdVT_TrkIdT_MediumIsoPFTau20_TauObj = 121,
  kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_EleObj = 122,
  kHLT_Ele20_CaloIdVT_TrkIdT_LooseIsoPFTau20_TauObj = 123,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj = 124,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj = 125,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj = 126,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj = 127,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj = 128,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj = 129,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_EleObj = 130,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_TauObj = 131,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_EleObj = 132,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_TauObj = 133,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj = 134,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj = 135,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj = 136,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj = 137,
  kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj = 138,
  kHLT_Ele20_CaloIdVT_CaloIsoRhoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj = 139,

  //Tau+Tau
  kHLT_DoubleIsoPFTau35_Trk5_eta2p1Obj = 140,
  kHLT_DoubleIsoPFTau25_Trk5_eta2p1Obj = 141,
  kHLT_DoubleIsoPFTau20_Trk5Obj        = 142,
  kHLT_DoubleMediumIsoPFTau30_Trk1_eta2p1_Jet30Obj = 143,

  //MuTau, additional 2012
  kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_MuObj = 144,
  kHLT_IsoMu17_eta2p1_LooseIsoPFTau20_TauObj = 145,
  kHLT_Mu17_eta2p1_LooseIsoPFTau20_MuObj   = 146,
  kHLT_Mu17_eta2p1_LooseIsoPFTau20_TauObj  = 147,

  //ETau, additional 2012
  kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_EleObj = 148,
  kHLT_Ele22_eta2p1_WP90Rho_LooseIsoPFTau20_TauObj = 149,
  kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_EleObj = 150, 
  kHLT_Ele22_eta2p1_WP90NoIso_LooseIsoPFTau20_TauObj = 151,

  //Tau Tau, additional 2012
  kHLT_DoubleMediumIsoPFTau30_Trk5_eta2p1_Jet30Obj = 152,
  kHLT_DoubleMediumIsoPFTau25_Trk5_eta2p1_Jet30Obj = 153
};

#endif
