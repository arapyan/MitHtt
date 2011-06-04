#ifndef HIGGSANA_NTUPLER_HIGGSANADEFS_HH 
#define HIGGSANA_NTUPLER_HIGGSANADEFS_HH

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

enum ECiC
{
  kMedium      = 0x001,  
  kTight       = 0x002,
  kSuperTight  = 0x004,
  kHyperTight1 = 0x008,
  kHyperTight2 = 0x010,
  kHyperTight3 = 0x020,
  kHyperTight4 = 0x040
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
  kHLT_Mu8                                                 = 0x00000001, // data
  kHLT_Mu9                                                 = 0x00000001, // MC
  kHLT_Mu15                                                = 0x00000002, // data, MC
  kHLT_Mu21                                                = 0x00000004, // MC
  kHLT_Mu24                                                = 0x00000004, // data
  kHLT_IsoMu17                                             = 0x00000008, // data, MC

  kHLT_Ele8                                                = 0x00000010, // data
  kHLT_Ele8_CaloIdL_CaloIsoVL                              = 0x00000020, // data
  kHLT_Ele17_SW_L1R                                        = 0x00000040, // MC
  kHLT_Ele17_CaloIdL_CaloIsoVL                             = 0x00000040, // data
  kHLT_Ele22_SW_L1R                                        = 0x00000080, // MC
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT              = 0x00000100, // data
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL      = 0x00000200, // data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 = 0x00000400, // data
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 0x00000800, // data
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL       = 0x00001000, // data

  kHLT_DoubleMu5                                           = 0x00002000, // MC
  kHLT_DoubleMu7                                           = 0x00002000, // data
  kHLT_Mu5_Jet50U                                          = 0x00004000, // MC
  kHLT_Mu8_Jet40                                           = 0x00004000, // data      
  kHLT_Ele8_CaloIdL_CaloIsoVL_Jet40                        = 0x00008000, // data
  
//   kHLT_Mu11_Ele8                                           = 0x00010000, // MC
  kHLT_Mu17_Ele8_CaloIdL                                   = 0x00010000, // data
//   kHLT_Mu8_Ele8                                            = 0x00020000, // MC
  kHLT_Mu8_Ele17_CaloIdL                                   = 0x00020000, // data
  kHLT_Mu8_Photon20_CaloIdVT_IsoT                          = 0x00040000, // data

  kHLT_Jet15U                                              = 0x00080000, // MC
  kHLT_Jet30U                                              = 0x00100000, // MC
  kHLT_Jet30                                               = 0x00100000, // data, MC
  kHLT_DiJetAve15U                                         = 0x00200000, // MC
  kHLT_DiJetAve30U                                         = 0x00400000, // data, MC

  kHLT_Photon30_Cleaned_L1R                                = 0x00800000, // MC
  kHLT_Photon30_CaloIdVL_IsoL                              = 0x00800000, // data
  kHLT_Photon50_Cleaned_L1R                                = 0x01000000, // MC
  kHLT_Photon50_CaloIdVL_IsoL                              = 0x01000000, // data
  
  kHLT_DoubleMu0_Quarkonium                                = 0x02000000, // MC
  kHLT_DoubleMu3_Quarkonium                                = 0x04000000, // data
  kHLT_DoubleMu3_Jpsi                                      = 0x08000000, // data
  kHLT_Dimuon0_Barrel_Upsilon                              = 0x10000000, // data
  kHLT_Dimuon6p5_Barrel_Jpsi                               = 0x20000000, // data
  kHLT_Dimuon6p5_Jpsi                                      = 0x40000000, // data
  
  kHLT_Mu13_Mu8                                            = 0x80000000  // data
//kHLT_Ele32_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT4
//kHLT_Mu12
};
  
/*
enum ETriggerBit
{  
  kHLT_Mu8                                                 = 0x00000001, // data
  kHLT_Mu9                                                 = 0x00000001, // MC
  kHLT_Mu15                                                = 0x00000002, // data, MC
  kHLT_Mu19                                                = 0x00000004, // MC
  kHLT_Mu20                                                = 0x00000004, // data
  kHLT_Mu24                                                = 0x00000008, // data
  kHLT_IsoMu11                                             = 0x00000010, // MC
  kHLT_IsoMu12                                             = 0x00000010, // data
  kHLT_IsoMu15                                             = 0x00000020, // data, MC

  kHLT_Ele10_SW_L1R                                        = 0x00000040, // MC
  kHLT_Ele8                                                = 0x00000040, // data
  kHLT_Ele8_CaloIdL_CaloIsoVL                              = 0x00000080, // data
  kHLT_Ele17_SW_L1R                                        = 0x00000100, // MC
  kHLT_Ele17_CaloIdL_CaloIsoVL                             = 0x00000100, // data
  kHLT_Ele22_SW_L1R                                        = 0x00000200, // MC
  kHLT_Ele27_CaloIdVT_CaloIsoT_TrkIdT_TrkIsoT              = 0x00000200, // data
  kHLT_Ele17_CaloIdL_CaloIsoVL_Ele8_CaloIdL_CaloIsoVL      = 0x00000400, // data
  kHLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_SC8_Mass30 = 0x00000800, // data
  kHLT_Ele17_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL_Ele8_CaloIdT_TrkIdVL_CaloIsoVL_TrkIsoVL = 0x00001000, // data
  kHLT_Photon20_CaloIdVT_IsoT_Ele8_CaloIdL_CaloIsoVL       = 0x00002000, // data
  kHLT_Ele32_CaloIdL_CaloIsoVL_SC17                        = 0x00004000, // data

  kHLT_DoubleMu7                                           = 0x00008000, // data
  kHLT_Mu5_Jet50U                                          = 0x00010000, // MC
  kHLT_Mu8_Jet40                                           = 0x00010000, // data

  kHLT_Mu11_Ele8                                           = 0x00020000, // MC
  kHLT_Mu17_Ele8_CaloIdL                                   = 0x00020000, // data
  kHLT_Mu8_Ele8                                            = 0x00040000, // MC
  kHLT_Mu8_Ele17_CaloIdL                                   = 0x00040000, // data
  kHLT_Mu8_Photon20_CaloIdVT_IsoT                          = 0x00080000, // data

  kHLT_Jet30                                               = 0x00100000, // data, MC
  kHLT_Jet50U                                              = 0x00200000, // MC
  kHLT_Jet60                                               = 0x00200000, // data
  kHLT_Jet100U                                             = 0x00400000, // MC
  kHLT_Jet110                                              = 0x00400000, // data

  kHLT_DoubleMu0                                           = 0x00800000, // MC  
  kHLT_DoubleMu0_Quarkonium                                = 0x01000000, // MC
  kHLT_DoubleMu3_Quarkonium                                = 0x02000000, // data
  kHLT_DoubleMu3_Upsilon                                   = 0x04000000, // data

  kHLT_Photon30_Cleaned_L1R                                = 0x08000000, // MC
  kHLT_Photon30_CaloIdVL_IsoL                              = 0x08000000, // data
  kHLT_Photon50_Cleaned_L1R                                = 0x10000000, // MC
  kHLT_Photon50_CaloIdVL_IsoL                              = 0x10000000, // data
  kHLT_Photon70_Cleaned_L1R                                = 0x20000000, // MC
  kHLT_Photon75_CaloIdVL_IsoL                              = 0x20000000  // data
};
*/
#endif
