#ifndef MITHTT_NTUPLER_HTTNTUPLER_H
#define MITHTT_NTUPLER_HTTNTUPLER_H

#include "TFile.h"
#include "TTree.h"

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/Particle.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupInfoFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitAna/DataTree/interface/EmbedWeight.h"
#include "MitAna/DataTree/interface/EmbedWeightFwd.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataTree/interface/Names.h"

#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TPFTau.hh"
#include "MitHtt/Ntupler/interface/TPhoton.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TMet.hh"
#include "MitHtt/Ntupler/interface/TSVfit.h"

#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/MetSignificance.h"
#include "MitPhysics/Utils/interface/ElectronTools.h"
#include "MitPhysics/Utils/interface/MuonTools.h"
#include "MitPhysics/Utils/interface/ElectronIDMVA.h"
#include "MitPhysics/Utils/interface/TauIsoMVA.h"
#include "MitPhysics/Utils/interface/MuonIDMVA.h"
#include "MitPhysics/Utils/interface/JetIDMVA.h"
#include "MitPhysics/Utils/interface/MVAMet.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "MitHtt/Ntupler/interface/AntiElectronIDMVA.h"

/// forward declarations
class TString;

/**
   \class HttNtupler HttNtupler.h MitHtt/Ntupler/include/HttNtupler.h

   \brief Description: Main module to do the filling from Bambu to Bacon

   Implementation:
   <Notes on implementation>

   useGen :   1 : Fill generator information for Higgs
              2 : Fill generator information for Z Boson
              3 : Fill generator information for W Boson
              4 : Fill generator information for DiBoson
*/

namespace mithep
{  
  class HttNtupler : public BaseMod
  {    
    public:

    /// default contructor
    HttNtupler(const char *name="HttNtupler", const char *title="Bambu2Bacon");
    /// default destructor
    ~HttNtupler();	    
    
    /// name of the output file
    void SetOutputName(const char *f){ fOutputName = f; }
    /// set the data flag: 0=MC, 1=data, 2=embedded sample
    void SetIsData(const int flag){ fIsData = flag; }
    /// add generator information to the bacon tree
    void SetUseGen(const int useGen){ fUseGen = useGen; }
    /// consider trigger information  
    void Set2012(const int flag){ f2012 = flag; }
     /// set minimum pt for photons
    void SetSkipIfHLTFail(const bool flag){ fSkipIfHLTFail = flag; } 
    /// set minimum pt for muons
    void SetMuPtMin(const double pt){ fMuPtMin = pt; }
    /// set maximum pt for muons
    void SetMuPtMax(const double pt){ fMuPtMax = pt; }
    /// set minimum eta for muons
    void SetMuEtaMin(const double eta){ fMuEtaMin = eta; }
    /// set maximum eta for muons
    void SetMuEtaMax(const double eta){ fMuEtaMax = eta; }
    /// set minimum Et for electrons
    void SetEleEtMin(const double et){ fEleEtMin = et; }
    /// set maximum Et for electrons
    void SetEleEtMax(const double et){ fEleEtMax = et; }
    /// set minimum eta for electrons
    void SetEleEtaMin(const double eta){ fEleEtaMin = eta; }
    /// set maximum eta for electrons
    void SetEleEtaMax(const double eta){ fEleEtaMax = eta; }
    /// set minimum pt for taus
    void SetTauPtMin(const double pt){ fPFTauPtMin  = pt; }
    /// set maximum eta for taus
    void SetTauEtaMax(const double eta){ fPFTauEtaMax = eta; }
    /// set minimum pt for jets
    void SetJetPtMin(const double pt){ fJetPtMin = pt; }
    /// set minimum pt for photons
    void SetPhotonEtMin(const double et){ fPhotonEtMin = et; }
    /// set minimum number of tracks for primary vertex selection
    void SetMinNTracksFit(const unsigned int ntrks){ fMinNTracksFit = ntrks; }
    /// set minimum number of degrees of freedom for vertex fit for primary vertex selection
    void SetMinNdof(const double ndof){ fMinNdof = ndof; }
    /// set maximal distance of primary vertices along z-axis for primary vertex selection
    void SetMaxAbsZ(const double z){ fMaxAbsZ = z; }
    /// set maximal distance in rho-phi plane for primary vertex selection
    void SetMaxRho(const double rho){ fMaxRho = rho; }
    /// decide whether to print the HLT table at the beginning of each run or not
    void SetPrintHLT(const bool flag){ fPrintTable = flag; }
    /// add a given trigger path to the analysis and perform a given matching
    void AddTrigger(const char* name, const unsigned int id, const char* objName1="", const unsigned int objId1=0, const double minPt1=0, const char* objName2="", const unsigned int objId2=0, const double minPt2=0);
    /// add a given jet corrector to the list of jet correctors     
    void AddJetCorr(const char* name);
    /// add a given JSON file to the list of JSON files
    void AddJSON(const char* name);
    
    protected:

    /// everythin that needs to be done before starting the event loop
    void Begin();
    /// everything that needs to be done at the beginning of each run
    void BeginRun();
    /// everything that needs to be done at the beginning of each run
    void SlaveBegin();
    /// everything that needs to be done at the end of each run
    void EndRun();
    /// everything that needs to be done at the end of each run
    void SlaveTerminate();
    /// everything that needs to be done at the end of the event loop
    void Terminate();
    ///everything that needs to be done for each event
    void Process();
    /// load HLT trigger table from Bambu at the beginning of each event   
    void loadTriggerTable(TriggerBits& trigBits);
    /// load all relevant Bambu branches at the beginning of each event
    void loadBambuBranches();
     /// clear all arrays before filling
    void resetOutputArrays();
    /// fill generator infomration depending on the indicated sample type (distinguished by enumerator ESampleType)
    void fillGenerator();
    /// fill general event information from Bambu
    void fillCommon(TriggerBits& trigBits);
    /// fill muons from Bambu
    void fillMuons();
    /// fill electrons from Bambu
    void fillElecs();
    /// fill taus from Bambu
    void fillPFTaus();
    /// fill input information needed to run svfit
    void fillSVfit();
    /// fill particle flow jets from Bambu
    void fillJets();
    /// fill photons from Bambu
    void fillPhotons();
    /// request branches for input collections
    void setupInput();
    /// book output collections 
    void setupOutput();
    /// cleanup output collections before termination
    void cleanup();
    /// fill boson and daughters for the adding of gen information
    void fillMCParticles(const MCParticle*& boson, const MCParticle*& daughterA, const MCParticle*& daughterB);
    /// visible momentum of the tau lepton (on generator level)
    FourVectorM visibleMCMomentum(const MCParticle* tauLepton);
    /// step down the daughters of the particle until the given status is reached or until the particle does not have daughters any more (for the latter set status==0)
    void stableMCParticle(const MCParticle*& part, int status=0);
    /// determine the pdgId of the leptons, which is the PDG standard in most cases but can be special for leptonic/hadronic tau decays
    int pdgId(const MCParticle* part);
    /// match offline reconstructed object to object on HLT level (potentially including minimum pt cut)
    TriggerObjects matchHLT(const double eta, const double phi, const double pt=0.);
    /// check whether this electron is classified to come from a photon conversion according to MVA
    bool isConversion(const Electron* elec);
    /// get the vertex, which associated to the tau lepton for official tau isolation (not implemented)
    const Vertex* officialTauIsoAssociatedVertex(const PFTau* iTau);  
    /// quality cuts on the particle flow candidate for official tau isolation 
    bool officialTauIsoPFCandidateQuality(const PFCandidate* cand, const Vertex *vtx, bool invertDz=false, bool invertGeneral=false);
    /// separate PFCandidate collection into PileUp and NoPileUp (re-implementing PF2PAT)
    void separatePileup(const PFCandidateCol* pfCandidates, const Vertex* pv, const VertexCol* primVerts, PFCandidateOArr* pfPileup, PFCandidateOArr* pfNoPileup, bool checkClosestZVertex = kTRUE);
    /// compute particle flow isolation for muon with charged hadrons restricted to the primary vertex in use (no correction for neutrals yet though)
    float computePFMuonIso(const Muon* muon, const double dRMax);
    /// compute particle flow isolation for electron with charged hadrons restricted to the primary vertex in use (no correction for neutrals yet though)
    float computePFElecIso(const Electron* elec, const double dRMax);
    /// check whether PF candidate is a good lepton
    bool findLeptonFootprint(const PFCandidate* pfcand);
    /// compute particle flow isolation for tau with charged hadrons restricted to the primary vertex in use (no correction for neutrals yet though)
    float computePFTauIso(const PFTau* tau, const Track* track);
    /// compute official POG tau isolation
    float computeOfficialPFTauIso(const PFTau* tau);
    /// old particle flow isolation. The isoType can be: 1: (Charged)Hadron+Electron+Muon; 2: (Neutral)Hadron; 3: Photon
    double computeCommonIso(const ChargedParticle* p, PFCandidateCol* pfNoPileup, double ptMin, double extRadius, double intRadius, int isoType, double dzLep=0, double dzMax=-999.);
    /// old naive particle flow isolation (including charged and neutral pileup)
    double computeNaiveIso(const Particle* p, PFCandidateCol* pfPileup, double ptMin, double extRadius, double intRadius);
    /// does this tau fullfil the loose tau Id?
    bool looseTauId(const PFTau* iTau) ;
    /// does this muon fullfil the loose Muon Id?
    bool looseMuId(const Muon * imu);
    /// does this electron fullfil the loose electron Id? conv:check if it is a conversion electron
    bool looseEleId(const Electron *iElectron, bool conv);
    /// fill input information for svfit for a given svfit array
    void fillSVfit(TClonesArray*& iArr, Particle* lep1, unsigned int lepId1, Particle* lep2, unsigned int lepId2, TMatrixD iMatrix, double dcaSig3D, double dcaSig2D, double dca3DErr, double dca2DErr);
    /// fill input information for svfit for a given svfit array
    void fillSVfit(TClonesArray*& iArr, Particle* lep1, unsigned int lepId1, Particle* lep2, unsigned int lepId2, TMatrixD iMatrix);

  protected:

    /// output file
    TFile* fOutputFile;           
    /// output tree
    TTree* fEventTree;

    /// name of the output file
    TString fOutputName;     
    /// name of the beam spot branch in Bambu
    TString fBeamSpotName;
    /// name of the primary vertex collection in Bambu 
    TString fPrimVtxName;
    /// name of the trigger mask table in Bambu
    TString fTrigMaskName;
    /// name of pile-up info branch in Bambu
    TString fPileupName;
    /// name of the energy density rho for pileup corrections with fastjet in Bambu
    TString fPUEnergyDensityName;
    /// name of the MC particle collection in Bambu
    TString fPartName;
    /// name of the MC event info in Bambu
    TString fMCEvtInfoName;
    /// name of the muon collection in Bambu
    TString fMuonName;
    /// name of the electron collection in Bambu
    TString fElectronName;
    /// name of the tau collection in Bambu
    TString fHPSTauName;
    /// name of the particle flow jet collection in Bambu
    TString fPFJetName;
    /// name of the particle flow MET collection in Bambu
    TString fPFMetName;
    /// name of the photon collection in Bambu
    TString fPhotonName;
    /// name of the conversion collection in Bambu
    TString fConversionName;
    /// name of the particle flow candidate collection in Bambu
    TString fPFCandidateName;
    /// name of the embedded weight collection for the embedding sample in Bambu
    TString fEmbedWeightName;
    /// name of dca significance collection for lepton pairs in Bambu
    //TString fDCASigName;

    /// generator particles
    const MCParticleCol* fParticles;
    /// MC information       
    const MCEventInfo* fMCEvtInfo;
    /// generator jets 
    const GenJetCol* fGenJets;
    /// beam spot
    const BeamSpotCol* fBeamSpot;
    /// primary vertices
    const VertexCol* fPrimVerts;
    /// trigger mask
    const TriggerMask* fTrigMask;
    /// pileup information
    const PileupInfoCol* fPileup;
    /// energy density rho for fastjet pileup corrections
    const PileupEnergyDensityCol* fPUEnergyDensity;
    /// muons
    const MuonCol* fMuons;
    /// electrons
    const ElectronCol* fElectrons;
    /// taus
    const PFTauCol* fPFTaus;
    /// particle flow jets
    const PFJetCol* fPFJets;
    /// particle flow MET
    const PFMetCol* fPFMet;
    /// photons
    const PhotonCol* fPhotons;
    /// photon conversions
    const DecayParticleCol* fConversions;
    /// particle flow candidates
    const PFCandidateCol* fPFCandidates;
    /// weights for embedded sample(s)
    const EmbedWeightCol* fEmbedWeight;
    /// dca significance for lepton pairs
    //const DCASigCol* fDCASigs;
    
    /// flag to indicate if processing collision data
    bool fIsData;
    /// flag whether to look at generator info
    int fUseGen;
    /// flag whether to print out HLT table
    bool fPrintTable;
    /// flag whether to skip event processing if HLT does not accept
    bool fSkipIfHLTFail;
    /// flag whether 2012 or 2011 samples  
    bool f2012;
    /// minimum pt for muons
    double fMuPtMin;
    /// maximum pt for muons
    double fMuPtMax;
    /// minimum eta for muons
    double fMuEtaMin;
    /// maximum eta for muons
    double fMuEtaMax;
    /// minimum supercluster Et for electrons
    double fEleEtMin;
    /// maximum supercluster Et for electrons
    double fEleEtMax;
    /// minimum supercluster eta for electrons
    double fEleEtaMin;
    /// maximum supercluster eta for electrons
    double fEleEtaMax;
    /// minimum pt for taus
    double  fPFTauPtMin;
    /// maximum eta for taus
    double fPFTauEtaMax;
    /// minimum pt for jets
    double fJetPtMin;
    /// minimum supercluster Et for photons
    double fPhotonEtMin;
    /// minimum number of tracks used for a good primary vertex
    unsigned int fMinNTracksFit;
    /// minimum degrees of freedom for a good primary vertex
    double fMinNdof;
    /// maximum z displacement for a good primary vertex
    double fMaxAbsZ;
    /// maximum transverse displacement for a good primary vertex 
    double fMaxRho;
    /// maximum distance in dz for a good primary vertex
    double fDzMax;
    
    /// pfPileUp collection including the primary vertex association along the z-axis   
    PFCandidateOArr* fPFPileUp;
    /// pfNoPileUp collection including the primary vertex association along the z-axis   
    PFCandidateOArr* fPFNoPileUp;
    /// pfPileUp collection w/o the primary vertex association along the z-axis (needed?)  
    PFCandidateOArr* fPFPileUpNoZ;
    /// pfPileNOUp collection w/o the primary vertex association along the z-axis (needed?)   
    PFCandidateOArr* fPFNoPileUpNoZ;
    
    /// selected primary vertex
    const Vertex* fVertex;
    /// general event information
    TEventInfo fEventInfo;
    /// trigger bits
    TriggerBits fTriggerBits;
    /// generator information
    TGenInfo fGenInfo;
    /// muons
    TClonesArray* fMuonArr;
    /// electrons
    TClonesArray* fElectronArr;
    /// taus
    TClonesArray* fHPSTauArr;
    /// particle flow jets
    TClonesArray* fPFJetArr;
    /// met
    //TClonesArray* fMetArr;
    /// photons
    TClonesArray* fPhotonArr;
    /// valid primary vertices
    TClonesArray* fPVArr;
    /// input information for svfit for emu
    TClonesArray* fSVfitEMuArr;
    /// input information for svfit for mutau
    TClonesArray* fSVfitMuTauArr;
    /// input information for svfit for etau
    TClonesArray* fSVfitETauArr;
     /// input information for svfit for tautau
    TClonesArray* fSVfitTauTauArr;
    
    /// names of triggers relevant for the analysis
    std::vector<TString> fTriggerNamesv;
    /// corresponding trigger bits
    std::vector<unsigned int> fTriggerIdsv;
    /// trigger object names for first lepton 
    std::vector<TString> fTriggerObjNames1v;
    /// trigger object Ids for first lepton
    std::vector<unsigned int> fTriggerObjIds1v;
    /// minimal pt cut for trigger matching for first lepton 
    std::vector<double> fTriggerObjMinPt1v;
    /// trigger object names for second lepton 
    std::vector<TString> fTriggerObjNames2v;
    /// trigger object Ids for second lepton 
    std::vector<unsigned int> fTriggerObjIds2v;
    /// minimal pt cut for trigger matching for second lepton 
    std::vector<double> fTriggerObjMinPt2v;
    
    /// list of jet correction parameters corresponding to a given jet correction level
    std::vector<TString> fJetCorrParsv;
    /// factorized jet corrector 
    FactorizedJetCorrector* fJetCorrector;
    /// jet correction uncertainties
    JetCorrectionUncertainty* fJetUncertainties;
    /// electron tools
    ElectronTools* fEleTools;
    /// muon tools
    MuonTools* fMuonTools;
    //met significance
    MetSignificance* metSign;
    /// electron ID MVA
    ElectronIDMVA* fElectronMVAID;
    /// electron ID MVA for triggered electrons
    ElectronIDMVA* fElectronMVAIDTrig;
    /// muon ID MVA
    JetIDMVA* fJetIDMVA;
    /// MET MVA
    MVAMet* fMVAMet;
    /// Tau ring ISO
    TauIsoMVA * fTauMVAIso;
    //AntiElectron ID MVA
    AntiElectronIDMVA * fAntiElectronIDMVA;
    /// list JSON files to be applied
    std::vector<TString> fJSONv;
    /// map of certified runs and lumi sections (for internal use)
    RunLumiRangeMap frlrm;
    /// run and lumi information (for internal use)
    RunLumiSet fRunLumiSet;
    /// run and lumi information (for intenral use)
    RunLumiRangeMap::RunLumiPairType fLastRunLumi;
    
    /// root class definition for schema evolution
    ClassDef(HttNtupler, 1)
  };

  inline void
  HttNtupler::setupInput()
  {
    ReqBranch( fMuonName            , fMuons           );
    ReqBranch( fElectronName        , fElectrons       );
    ReqBranch( fHPSTauName          , fPFTaus          );
    ReqBranch( fPrimVtxName         , fPrimVerts       );
    ReqBranch( fBeamSpotName        , fBeamSpot        );
    ReqBranch( fPFJetName           , fPFJets          );
    ReqBranch( fTrigMaskName        , fTrigMask        );
    ReqBranch( fPFMetName           , fPFMet           );
    ReqBranch( fPhotonName          , fPhotons         );
    ReqBranch( fConversionName      , fConversions     );
    ReqBranch( fPUEnergyDensityName , fPUEnergyDensity );
    ReqBranch( fPFCandidateName     , fPFCandidates    );
    ReqBranch( fPartName            , fParticles       );       
    ReqBranch( fPileupName          , fPileup          );
    ReqBranch( Names::gkGenJetBrn   , fGenJets         );
    ReqBranch( fMCEvtInfoName       , fMCEvtInfo       );
    ReqBranch( fEmbedWeightName     , fEmbedWeight     ); 
    //ReqBranch( fDCASigName          , fDCASigs         );
  } 

  inline void
  HttNtupler::setupOutput()
  {
    // book output and helper arrays
    fPFPileUp      = new PFCandidateOArr; fPFNoPileUp    = new PFCandidateOArr;
    fPFPileUpNoZ   = new PFCandidateOArr; fPFNoPileUpNoZ = new PFCandidateOArr;
    fMuonArr       = new TClonesArray( "mithep::TMuon"     ); assert( fMuonArr       );
    fElectronArr   = new TClonesArray( "mithep::TElectron" ); assert( fElectronArr   );
    fHPSTauArr     = new TClonesArray( "mithep::TPFTau"    ); assert( fHPSTauArr     ); 
    fPFJetArr      = new TClonesArray( "mithep::TJet"      ); assert( fPFJetArr      );
    //fMetArr      = new TClonesArray( "mithep::TMet"      ); assert( fMetArr        );
    fPhotonArr     = new TClonesArray( "mithep::TPhoton"   ); assert( fPhotonArr     );
    fPVArr         = new TClonesArray( "mithep::TVertex"   ); assert( fPVArr         );
    fSVfitEMuArr   = new TClonesArray( "mithep::TSVfit"    ); assert( fSVfitEMuArr   );
    fSVfitETauArr  = new TClonesArray( "mithep::TSVfit"    ); assert( fSVfitETauArr  );
    fSVfitMuTauArr = new TClonesArray( "mithep::TSVfit"    ); assert( fSVfitMuTauArr );    
    fSVfitTauTauArr = new TClonesArray( "mithep::TSVfit"    ); assert( fSVfitTauTauArr );    
    // open file and configure branches
    fOutputFile    = new TFile( fOutputName, "RECREATE" );
    fEventTree     = new TTree( "Events"   , "Events"   );
    if( (fIsData!=1) && fUseGen ){ fEventTree->Branch( "Gen", &fGenInfo ); }
    fEventTree->Branch( "Info"       , &fEventInfo      );
    fEventTree->Branch( "Muon"       , &fMuonArr        );
    fEventTree->Branch( "Electron"   , &fElectronArr    );
    fEventTree->Branch( "HPSTau"     , &fHPSTauArr      );
    fEventTree->Branch( "PFJet"      , &fPFJetArr       );
    //fEventTree->Branch( "Met"      , &fMetArr         );
    fEventTree->Branch( "Photon"     , &fPhotonArr      );
    fEventTree->Branch( "PV"         , &fPVArr          );
    fEventTree->Branch( "SVfitEMu"   , &fSVfitEMuArr    );
    fEventTree->Branch( "SVfitETau"  , &fSVfitETauArr   );
    fEventTree->Branch( "SVfitMuTau" , &fSVfitMuTauArr  );
    fEventTree->Branch( "SVfitTauTau" , &fSVfitTauTauArr  );
  }

  inline void
  HttNtupler::cleanup()
  {
    delete fPFPileUp; 
    delete fPFNoPileUp; 
    delete fPFPileUpNoZ; 
    delete fPFNoPileUpNoZ; 
    delete fMuonArr; 
    delete fElectronArr; 
    delete fHPSTauArr; 
    delete fPFJetArr; 
    //delete fMetArr; 
    delete fPhotonArr; 
    delete fPVArr; 
    delete fSVfitEMuArr; 
    delete fSVfitETauArr; 
    delete fSVfitMuTauArr;     
    delete fSVfitTauTauArr;     
  }

  inline void 
  HttNtupler::AddTrigger(const char* name, const unsigned int id, const char* objName1, const unsigned int objId1, const double minPt1, const char* objName2, const unsigned int objId2, const double minPt2)
  {
    fTriggerNamesv    .push_back( name     );
    fTriggerIdsv      .push_back( id       );
    fTriggerObjNames1v.push_back( objName1 );
    fTriggerObjIds1v  .push_back( objId1   );
    fTriggerObjMinPt1v.push_back( minPt1   );
    fTriggerObjNames2v.push_back( objName2 );
    fTriggerObjIds2v  .push_back( objId2   );
    fTriggerObjMinPt2v.push_back( minPt2   );
  }

  inline void 
  HttNtupler::AddJetCorr(const char* name) 
  {
    ifstream jetcorrchk;
    jetcorrchk.open(name); assert(jetcorrchk.is_open()); jetcorrchk.close();
    fJetCorrParsv.push_back(name);
  }
  
  inline void 
  HttNtupler::AddJSON(const char* name) 
  {
    ifstream jsonchk;
    jsonchk.open(name); assert(jsonchk.is_open()); jsonchk.close();
    fJSONv.push_back(name);
    cout << "adding " << name << endl;
  }

}

#endif
