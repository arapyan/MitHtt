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
    kW           = 23,
    kZ           = 24,
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
  /* SingleMu CrossObject       ---------------------------------------------------------------------------------------------- */
  kHLT_IsoMu15_L1ETM20                                                                  = 19, // data
  kHLT_Mu15_LooseIsoPFTau15                                                             = 20,
  kHLT_Mu15_LooseIsoPFTau20                                                             = 21,
  kHLT_IsoMu12_LooseIsoPFTau10                                                          = 22,
  kHLT_IsoMu15_LooseIsoPFTau15                                                          = 23,
  kHLT_IsoMu15_LooseIsoPFTau20                                                          = 24,
  kHLT_IsoMu15_MediumIsoPFTau20                                                         = 25,
  kHLT_IsoMu15_TightIsoPFTau20                                                          = 26,
  /* DoubleElectron             ----------------------------------------------------------------------------------------------- */
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL                                   = 27, // data
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 28, // data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30                              = 29, // data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30                             = 30, // data
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                                                     = 31, // data
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17                                       = 32, // data
  /* SingleElectron             ----------------------------------------------------------------------------------------------- */
  kHLT_Ele8                                                                             = 33, // data
  kHLT_Ele8_CaloIdL_TrkIdVL                                                             = 34, // data
  kHLT_Ele8_CaloIdL_CaloIsoVL                                                           = 35, // data
  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL                                          = 36, // data
  kHLT_Ele17_CaloIdL_CaloIsoVL                                                          = 37, // data
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40                                                     = 38, // data
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL                                    = 39, // data
  kHLT_Ele17_SW_L1R                                                                     = 40, // MC
  kHLT_Ele22_SW_L1R                                                                     = 41, // MC
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 42, // data
  kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 43, // data
  kHLT_Ele45_CaloIdVT_TrkIdT                                                            = 44, // data
  kHLT_Ele52_CaloIdVT_TrkIdT                                                            = 45, // data
  kHLT_Ele65_CaloIdVT_TrkIdT                                                            = 46, // data
  kHLT_Ele80_CaloIdVT_TrkIdT                                                            = 47, // data
  /* Old Triggers               ---------------------------------------------------------------------------------------------- */
  kHLT_Photon10_L1R			                                                = 48,
  kHLT_Photon15_Cleaned_L1R		                                                = 49,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 50,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 51,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT                                           = 52,
  kHLT_Ele22_SW_TighterCaloIdIsol_L1R	                                                = 53,
  /* SingleElectron CrossObject ----------------------------------------------------------------------------------------------- */
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15                           = 54,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20                           = 55,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20                           = 56,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20                           = 57,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20                          = 58,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20                          = 59,
  /* JetMET                   ------------------------------------------------------------------------------------------------ */
  kHLT_CentralJet80_MET100		                                                = 60,
  kHLT_CentralJet80_MET160		                                                = 61,
  kHLT_CentralJet80_MET80		                                                = 62,
  kHLT_DiCentralJet20_BTagIP_MET65	                                                = 63,
  kHLT_DiCentralJet20_MET80     	                                                = 64,
  kHLT_MET100_HBHENoiseFiltered 	                                                = 65,
  kHLT_MET100   			                                                = 66,
  kHLT_MET120_HBHENoiseFiltered  	                                                = 67,
  kHLT_MET120   			                                                = 68,
  /* SingleTau                 ------------------------------------------------------------------------------------------------- */
  kHLT_LooseIsoPFTau10                                                                  = 69,
  kHLT_LooseIsoPFTau15                                                                  = 70,
  kHLT_LooseIsoPFTau20                                                                  = 71,
  kHLT_MediumIsoPFTau15                                                                 = 72,
  kHLT_MediumIsoPFTau20                                                                 = 74,
  kHLT_TightIsoPFTau15                                                                  = 75,
  kHLT_TightIsoPFTau20                                                                  = 76,
  /* 2012 Triggers               ---------------------------------------------------------------------------------------------- */
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL                                     = 77,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL                                     = 78,
  kHLT_Ele27_WP80                                                                       = 79,
  kHLT_Mu17                                                                             = 80,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50                             = 81,
  kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50                              = 82,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50                                = 83,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30                                   = 84,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL                                         = 85,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30                                    = 86,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL                                          = 87
};

const unsigned int kNTrigObj = 128;
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
  /* SingleMuon CrossObject    ------------------------------------------------------------------------------------------------- */
  kHLT_IsoMu15_L1ETM20_MuObj                                                                    =  30,
  kHLT_IsoMu15_L1ETM20_METObj                                                                   =  31,
  kHLT_Mu15_LooseIsoPFTau15_MuObj                                                               =  32,
  kHLT_Mu15_LooseIsoPFTau15_TauObj                                                              =  33,
  kHLT_Mu15_LooseIsoPFTau20_MuObj                                                               =  34,
  kHLT_Mu15_LooseIsoPFTau20_TauObj                                                              =  35,
  kHLT_IsoMu12_LooseIsoPFTau10_MuObj                                                            =  36,
  kHLT_IsoMu12_LooseIsoPFTau10_TauObj                                                           =  37,
  kHLT_IsoMu15_LooseIsoPFTau15_MuObj                                                            =  38,
  kHLT_IsoMu15_LooseIsoPFTau15_TauObj                                                           =  39,
  kHLT_IsoMu15_LooseIsoPFTau20_MuObj                                                            =  40,
  kHLT_IsoMu15_LooseIsoPFTau20_TauObj                                                           =  41,
  kHLT_IsoMu15_MediumIsoPFTau20_MuObj                                                           =  42,
  kHLT_IsoMu15_MediumIsoPFTau20_TauObj                                                          =  43,
  kHLT_IsoMu15_TightIsoPFTau20_MuObj                                                            =  44,
  kHLT_IsoMu15_TightIsoPFTau20_TauObj                                                           =  45,
  /* DoubleElectron            ------------------------------------------------------------------------------------------------- */
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj                                   =  46,
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj                                   =  47,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj =  48,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj =  49,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj                               =  50,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj                                =  51,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj                             =  52,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj                             =  53,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj                                                      =  54,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj                                                       =  55,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj                                        =  56,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj                                         =  57,
  kHLT_Ele8_EleObj                                                                              =  58,
  kHLT_Ele8_CaloIdL_TrkIdVL_EleObj                                                              =  59,
  kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj                                                            =  60,
  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj                                           =  61,
  kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj                                                           =  62,
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj                                                      =  63,
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj                                                      =  64,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj                                     =  65,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj                                     =  66,
  /* SingleElectron             ------------------------------------------------------------------------------------------------- */  
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  67,
  kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  68,
  kHLT_Ele45_CaloIdVT_TrkIdT_EleObj                                                             =  69, // data
  kHLT_Ele52_CaloIdVT_TrkIdT_EleObj                                                             =  70, // data
  kHLT_Ele65_CaloIdVT_TrkIdT_EleObj                                                             =  71, // data
  kHLT_Ele80_CaloIdVT_TrkIdT_EleObj                                                             =  72, // data
  kHLT_Photon10_L1R_EleObj			                                                =  73,
  kHLT_Photon15_Cleaned_L1R_EleObj		                                                =  74,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  75,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  76,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj                                            =  77,
  kHLT_Ele22_SW_TighterCaloIdIsol_L1R_EleObj	                                                =  78,
  /* SingleTau CorssObject     ------------------------------------------------------------------------------------------------- */
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_EleObj                            =  80,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau15_TauObj                            =  81,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj                            =  82,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj                            =  83,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_EleObj                            =  84,
  kHLT_Ele15_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_TightIsoPFTau20_TauObj                            =  85,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_EleObj                            =  86,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_LooseIsoPFTau20_TauObj                            =  87,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj                           =  88,
  kHLT_Ele18_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj                           =  89,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_EleObj                           =  90,
  kHLT_Ele20_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_MediumIsoPFTau20_TauObj                           =  91,
  /* JetMET                    ------------------------------------------------------------------------------------------------- */  
  kHLT_Jet150_JetObj                                                                            =  92,
  kHLT_Jet190_JetObj                                                                            =  93,
  kHLT_Jet240_JetObj                                                                            =  94,
  kHLT_Jet300_JetObj                                                                            =  95,
  kHLT_DiJetAve30_JetObj                                                                        =  96,
  kHLT_DiJetAve240_JetObj                                                                       =  97,
  kHLT_DiJetAve300_JetObj                                                                       =  98,
  /* Photon                    ------------------------------------------------------------------------------------------------- */
  kHLT_Photon20_CaloIdVL_IsoL_PhoObj                                                            =  99,
  kHLT_Photon30_CaloIdVL_IsoL_PhoObj                                                            = 100,
  kHLT_Photon50_CaloIdVL_IsoL_PhoObj                                                            = 101,
  kHLT_Photon75_CaloIdVL_PhoObj                                                                 = 102,
  /* Tau                       ------------------------------------------------------------------------------------------------- */
  kHLT_LooseIsoPFTau10_TauObj                                                                   = 103,
  kHLT_LooseIsoPFTau15_TauObj                                                                   = 104,
  kHLT_LooseIsoPFTau20_TauObj                                                                   = 105,
  kHLT_MediumIsoPFTau15_TauObj                                                                  = 106,
  kHLT_MediumIsoPFTau20_TauObj                                                                  = 107,
  kHLT_TightIsoPFTau15_TauObj                                                                   = 108,
  kHLT_TightIsoPFTau20_TauObj                                                                   = 109,
  /* 2012                      ------------------------------------------------------------------------------------------------- */
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj                                       = 110,
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj                                       = 111,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_MuObj                                       = 112,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EGObj                                       = 113,
  kHLT_Ele27_WP80_EleObj                                                                        = 114,
  kHLT_Mu17_MuObj                                                                               = 115,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele1Obj                             = 116,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_Ele2Obj                             = 117,
  kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_EleObj                               = 118,
  kHLT_Ele20_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC4_Mass50_SCObj                                = 119,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_EleObj                                 = 120,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_Mass50_SCObj                                  = 121,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj                                    = 122,
  kHLT_Ele17_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj                                          = 123,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_Jet30_EleObj                                     = 124,
  kHLT_Ele8_CaloIdT_CaloIsoVL_TrkIdVL_TrkIsoVL_EleObj                                           = 125
};

#endif


/* 
left overs from Valentina's implementation:

enum ETriggerBit
{  
  kHLT_Mu3                                                 = (ULong64_t)1<<48, // data
  kHLT_Mu5                                                 = (ULong64_t)1<<49, // data
  kHLT_Mu30                                                = (ULong64_t)1<<12, // data
  kHLT_Mu40                                                = (ULong64_t)1<<40, // data
  kHLT_IsoMu30                                             = (ULong64_t)1<<14, // data
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17         = (ULong64_t)1<<48,
  kHLT_Photon10_L1R                                        = (ULong64_t)1<<32,
  kHLT_Photon15_Cleaned_L1R                                = (ULong64_t)1<<33,
  kHLT_Ele15_SW_CaloEleId_L1R                              = (ULong64_t)1<<34,
  kHLT_Ele17_SW_CaloEleId_L1R                              = (ULong64_t)1<<35,
  kHLT_Ele17_SW_TightEleId_L1R                             = (ULong64_t)1<<36,
  kHLT_Ele22_SW_TighterCaloIdIsol_L1R                      = (ULong64_t)1<<37,
  kHLT_Jet30                                               = (ULong64_t)1<<44,
  kHLT_Jet60                                               = (ULong64_t)1<<45,
  kHLT_Jet80                                               = (ULong64_t)1<<46,
  kHLT_Jet110                                              = (ULong64_t)1<<47
};

enum ETriggerObjBit
{
  // SingleMu
  kHLT_Mu3_MuObj                                           = (ULong64_t)1<<58,
  kHLT_Mu5_MuObj                                           = (ULong64_t)1<<59,
  kHLT_Mu40_MuObj                                          = (ULong64_t)1<<49,
  kHLT_IsoMu30_MuObj                                       = (ULong64_t)1<<50,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele1Obj = (ULong64_t)1<<54,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele2Obj = (ULong64_t)1<<55,
  kHLT_Jet30_JetObj                                        = (ULong64_t)1<<54,
  kHLT_Jet60_JetObj                                        = (ULong64_t)1<<55,
  kHLT_Jet80_JetObj                                        = (ULong64_t)1<<56,
  kHLT_Jet110_JetObj                                       = (ULong64_t)1<<57
};
*/
