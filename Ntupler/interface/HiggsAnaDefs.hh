#ifndef MITHTT_NTUPLER_HIGGSANADEFS_HH 
#define MITHTT_NTUPLER_HIGGSANADEFS_HH

namespace EGenType {
enum {
  kMuon        = 1,
  kElectron    = 2,
  kTau         = 3,
  kTauMuon     = 4,
  kTauElectron = 5,
  kTauHadr     = 6,
  kW           = 7,
  kZ           = 8,
  kWW          = 9,
  kHiggs       = 10
};
}

enum EMuType 
{ 
  kGlobal     = 1, 
  kTracker    = 2, 
  kStandalone = 4
};

enum EQualityBit
{ 
  // descriptions from DataFormats/MuonReco/interface/MuonSelectors.h
  kAll  			    = 0x000001,  // dummy options - always true
  kAllGlobalMuons		    = 0x000002,  // checks isGlobalMuon flag
  kAllStandAloneMuons		    = 0x000004,  // checks isStandAloneMuon flag
  kAllTrackerMuons		    = 0x000008,  // checks isTrackerMuon flag
  kTrackerMuonArbitrated	    = 0x000010,  // resolve ambiguity of sharing segments
  kAllArbitrated		    = 0x000020,  // all muons with the tracker muon arbitrated
  kGlobalMuonPromptTight	    = 0x000040,  // global muons with tighter fit requirements
  kTMLastStationLoose		    = 0x000080,  // penetration depth loose selector
  kTMLastStationTight		    = 0x000100,  // penetration depth tight selector
  kTM2DCompatibilityLoose	    = 0x000200,  // likelihood based loose selector
  kTM2DCompatibilityTight	    = 0x000400,  // likelihood based tight selector
  kTMOneStationLoose		    = 0x000800,  // require one well matched segment
  kTMOneStationTight		    = 0x001000,  // require one well matched segment
  kTMLastStationOptimizedLowPtLoose = 0x002000,  // combination of TMLastStation and TMOneStation
  kTMLastStationOptimizedLowPtTight = 0x004000,  // combination of TMLastStation and TMOneStation
  kGMTkChiCompatibility 	    = 0x008000,  // require tk stub have good chi2 relative to glb track
  kGMStaChiCompatibility	    = 0x010000,  // require sta stub have good chi2 compatibility relative to glb track
  kGMTkKinkTight		    = 0x020000,  // require a small kink value in the tracker stub
  kTMLastStationAngLoose	    = 0x040000,  // TMLastStationLoose with additional angular cuts
  kTMLastStationAngTight	    = 0x080000,  // TMLastStationTight with additional angular cuts
  kTMOneStationAngLoose 	    = 0x100000,  // TMOneStationLoose with additional angular cuts
  kTMOneStationAngTight 	    = 0x200000,  // TMOneStationTight with additional angular cuts
  //The two algorithms that follow are identical to what were known as
  //TMLastStationOptimizedLowPt* (sans the Barrel) as late as revision
  //1.7 of this file. The names were changed because indeed the low pt
  //optimization applies only to the barrel region, whereas the sel-
  //ectors above are more efficient at low pt in the endcaps, which is
  //what we feel is more suggestive of the algorithm name. This will be
  //less confusing for future generations of CMS members, I hope...
  kTMLastStationOptimizedBarrelLowPtLoose = 0x400000,  // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
  kTMLastStationOptimizedBarrelLowPtTight = 0x800000   // combination of TMLastStation and TMOneStation but with low pT optimization in barrel only
}; 

enum ETriggerBit
{  

  // MuEG
  kHLT_Mu11_Ele8                  = (ULong64_t)1<<0,  // MC
  kHLT_Mu17_Ele8_CaloIdL          = (ULong64_t)1<<0,  // data
  kHLT_Mu8_Ele8                   = (ULong64_t)1<<1,  // MC
  kHLT_Mu8_Ele17_CaloIdL          = (ULong64_t)1<<1,  // data
  kHLT_Mu8_Photon20_CaloIdVT_IsoT = (ULong64_t)1<<2,  // data
  kHLT_Mu15_Photon20_CaloIdL      = (ULong64_t)1<<3,  // data
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL = (ULong64_t)1<<38,  // data
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL = (ULong64_t)1<<39,  // data

  
  // DoubleMu
  kHLT_DoubleMu5  = (ULong64_t)1<<4,  // MC
  kHLT_DoubleMu7  = (ULong64_t)1<<4,  // data
  kHLT_Mu13_Mu8   = (ULong64_t)1<<5,  // data
  kHLT_Mu17_Mu8   = (ULong64_t)1<<6,  // data
  kHLT_Mu5_Jet50U = (ULong64_t)1<<7,  // MC
  kHLT_Mu8_Jet40  = (ULong64_t)1<<7,  // data
  
  // SingleMu
  kHLT_Mu3     = (ULong64_t)1<<48,  // data
  kHLT_Mu5     = (ULong64_t)1<<49,  // data
  kHLT_Mu8     = (ULong64_t)1<<8,  // data
  kHLT_Mu9     = (ULong64_t)1<<8,  // MC
  kHLT_Mu11    = (ULong64_t)1<<9,  // MC
  kHLT_Mu12    = (ULong64_t)1<<9,  // data
  kHLT_Mu15    = (ULong64_t)1<<10, // MC, data
  kHLT_Mu21    = (ULong64_t)1<<11, // MC
  kHLT_Mu24    = (ULong64_t)1<<11, // data
  kHLT_Mu30    = (ULong64_t)1<<12, // data
  kHLT_Mu40    = (ULong64_t)1<<40, // data
  kHLT_IsoMu17 = (ULong64_t)1<<13, // MC, data
  kHLT_IsoMu24 = (ULong64_t)1<<14, // data
  kHLT_IsoMu30 = (ULong64_t)1<<14, // data
  
  // DoubleElectron
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL                                   = (ULong64_t)1<<15, // data
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = (ULong64_t)1<<16, // data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30                              = (ULong64_t)1<<17, // data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30                             = (ULong64_t)1<<18, // data
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                                                     = (ULong64_t)1<<19, // data
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17                                       = (ULong64_t)1<<20, // data
  kHLT_Ele8                                                                             = (ULong64_t)1<<21, // data
  kHLT_Ele8_CaloIdL_TrkIdVL                                                             = (ULong64_t)1<<22, // data
  kHLT_Ele8_CaloIdL_CaloIsoVL                                                           = (ULong64_t)1<<23, // data
  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL                                          = (ULong64_t)1<<24, // data
  kHLT_Ele17_CaloIdL_CaloIsoVL                                                          = (ULong64_t)1<<25, // data
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40                                                     = (ULong64_t)1<<26, // data
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL                                    = (ULong64_t)1<<27, // data
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17                                      = (ULong64_t)1<<48,

  // SingleElectron
  kHLT_Ele17_SW_L1R                           = (ULong64_t)1<<28, // MC
  kHLT_Ele22_SW_L1R                           = (ULong64_t)1<<29, // MC
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = (ULong64_t)1<<30, // data
  kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT = (ULong64_t)1<<31,  // data
  kHLT_Ele52_CaloIdVT_TrkIdT                  = (ULong64_t)1<<41,
  kHLT_Ele65_CaloIdVT_TrkIdT                  = (ULong64_t)1<<42,
  kHLT_Ele80_CaloIdVT_TrkIdT                  = (ULong64_t)1<<43,
  
  kHLT_Photon10_L1R                   = (ULong64_t)1<<32,
  kHLT_Photon15_Cleaned_L1R           = (ULong64_t)1<<33,
  kHLT_Ele15_SW_CaloEleId_L1R         = (ULong64_t)1<<34,
  kHLT_Ele17_SW_CaloEleId_L1R         = (ULong64_t)1<<35,
  kHLT_Ele17_SW_TightEleId_L1R        = (ULong64_t)1<<36,
  kHLT_Ele22_SW_TighterCaloIdIsol_L1R = (ULong64_t)1<<37,

  kHLT_Jet30                          = (ULong64_t)1<<44,
  kHLT_Jet60                          = (ULong64_t)1<<45,
  kHLT_Jet80                          = (ULong64_t)1<<46,
  kHLT_Jet110                         = (ULong64_t)1<<47
};

enum ETriggerObjBit
{
  // MuEG 
  kHLT_Mu17_Ele8_CaloIdL_MuObj          = (ULong64_t)1<<0,
  kHLT_Mu17_Ele8_CaloIdL_EGObj          = (ULong64_t)1<<1,  
  kHLT_Mu8_Ele17_CaloIdL_MuObj          = (ULong64_t)1<<2,
  kHLT_Mu8_Ele17_CaloIdL_EGObj          = (ULong64_t)1<<3,
  kHLT_Mu8_Photon20_CaloIdVT_IsoT_MuObj = (ULong64_t)1<<4,
  kHLT_Mu8_Photon20_CaloIdVT_IsoT_EGObj = (ULong64_t)1<<5,
  kHLT_Mu15_Photon20_CaloIdL_MuObj      = (ULong64_t)1<<6,
  kHLT_Mu15_Photon20_CaloIdL_EGObj      = (ULong64_t)1<<7,
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_MuObj= (ULong64_t)1<<45,
  kHLT_Mu17_Ele8_CaloIdT_CaloIsoVL_EGObj= (ULong64_t)1<<46,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj= (ULong64_t)1<<47,
  kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj= (ULong64_t)1<<48,


  // DoubleMu
  kHLT_DoubleMu7_MuObj  = (ULong64_t)1<<8,
  kHLT_Mu13_Mu8_Mu1Obj  = (ULong64_t)1<<9,
  kHLT_Mu13_Mu8_Mu2Obj  = (ULong64_t)1<<10,
  kHLT_Mu17_Mu8_Mu1Obj  = (ULong64_t)1<<11,
  kHLT_Mu17_Mu8_Mu2Obj  = (ULong64_t)1<<12,
  kHLT_Mu8_Jet40_MuObj  = (ULong64_t)1<<13,
  kHLT_Mu8_Jet40_JetObj = (ULong64_t)1<<14,

  // SingleMu
  kHLT_Mu3_MuObj     = (ULong64_t)1<<58,
  kHLT_Mu5_MuObj     = (ULong64_t)1<<59,
  kHLT_Mu8_MuObj     = (ULong64_t)1<<15,
  kHLT_Mu9_MuObj     = (ULong64_t)1<<45,
  kHLT_Mu12_MuObj    = (ULong64_t)1<<16,
  kHLT_Mu15_MuObj    = (ULong64_t)1<<17,
  kHLT_Mu24_MuObj    = (ULong64_t)1<<18,
  kHLT_Mu30_MuObj    = (ULong64_t)1<<19,
  kHLT_Mu40_MuObj    = (ULong64_t)1<<49,
  kHLT_IsoMu17_MuObj = (ULong64_t)1<<20,
  kHLT_IsoMu24_MuObj = (ULong64_t)1<<21,
  kHLT_IsoMu30_MuObj = (ULong64_t)1<<50,

  // DoubleElectron
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele1Obj                                   = (ULong64_t)1<<22,
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL_Ele2Obj                                   = (ULong64_t)1<<23,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele1Obj = (ULong64_t)1<<24,
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele2Obj = (ULong64_t)1<<25,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_EleObj                               = (ULong64_t)1<<26,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30_SCObj                                = (ULong64_t)1<<27,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele1Obj                             = (ULong64_t)1<<28,
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass30_Ele2Obj                             = (ULong64_t)1<<29,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_EleObj                                                      = (ULong64_t)1<<30,
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17_SCObj                                                       = (ULong64_t)1<<31,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_EleObj                                        = (ULong64_t)1<<32,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_SC17_SCObj                                         = (ULong64_t)1<<33,
  kHLT_Ele8_EleObj                                                                              = (ULong64_t)1<<34,
  kHLT_Ele8_CaloIdL_TrkIdVL_EleObj                                                              = (ULong64_t)1<<35,
  kHLT_Ele8_CaloIdL_CaloIsoVL_EleObj                                                            = (ULong64_t)1<<36,
  kHLT_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_EleObj                                           = (ULong64_t)1<<37,
  kHLT_Ele17_CaloIdL_CaloIsoVL_EleObj                                                           = (ULong64_t)1<<38,
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_EleObj                                                      = (ULong64_t)1<<39,
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40_JetObj                                                      = (ULong64_t)1<<40,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_PhoObj                                     = (ULong64_t)1<<41,
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL_EleObj                                     = (ULong64_t)1<<42,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele1Obj                                      = (ULong64_t)1<<54,
  kHLT_Ele32_CaloIdT_CaloIsoT_TrkIdT_TrkIsoT_Ele17_Ele2Obj                                      = (ULong64_t)1<<55,
  
  // SingleElectron
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj = (ULong64_t)1<<43,
  kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT_EleObj = (ULong64_t)1<<44,
  kHLT_Ele52_CaloIdVT_TrkIdT_EleObj                  = (ULong64_t)1<<51,
  kHLT_Ele65_CaloIdVT_TrkIdT_EleObj                  = (ULong64_t)1<<52,
  kHLT_Ele80_CaloIdVT_TrkIdT_EleObj                  = (ULong64_t)1<<53,

  // Jet
  kHLT_Jet30_JetObj      = (ULong64_t)1<<54,
  kHLT_Jet60_JetObj      = (ULong64_t)1<<55,
  kHLT_Jet80_JetObj      = (ULong64_t)1<<56,
  kHLT_Jet110_JetObj     = (ULong64_t)1<<57
};
#endif
