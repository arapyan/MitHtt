#ifndef MITHTT_NTUPLER_HTTNTUPLERMOD_H
#define MITHTT_NTUPLER_HTTNTUPLERMOD_H

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

#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TPhoton.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Ntupler/interface/TSVfit.h"
#include "MitHtt/Ntupler/interface/MetSignificance.h"

#include "MitPhysics/Utils/interface/ElectronTools.h"

#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

/// forward declarations
class TTree;
class TFile;
class TString;

/**
   \class HttNtuplerMod HttNtuplerMod.h MitHtt/Ntupler/include/HttNtuplerMod.h

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
  class HttNtuplerMod : public BaseMod
  {    
    public:

    /// default contructor
    HttNtuplerMod(const char *name="HttNtuplerMod", const char *title="Bambu2Bacon");
    /// default destructor
    ~HttNtuplerMod();	    
    
    /// name of the output file
    void SetOutputName(const char *f){ fOutputName = f; }
    /// set the data flag: 0=MC, 1=data, 2=embedded sample
    void SetIsData(const int flag){ fIsData = flag; }
    /// add generator information to the bacon tree
    void SetUseGen(const int useGen){ fUseGen = useGen; }
    /// consider trigger information  
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
    void AddTrigger(const char* name, const ULong64_t id, const char* objName1="", const ULong64_t objId1=0, const double minPt1=0, const char* objName2="", const ULong64_t objId2=0, const double minPt2=0);
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
      
    /// fill generator information for Z bosons
    void FillGenZ();
    /// fill generator information for Higgs bosons
    void FillGenH();
    /// fill generator information for W bosons
    void FillGenW();
    /// fill generator information for DiBoson (really needed?)
    void FillGenWW();
    /// fill input information for svfit
    void FillSVfit(TClonesArray* iarr, Particle* ipart0, Particle* ipart1, TMatrixD iMatrix);
    /// fill TMuon from Bambu
    void FillMuon(const Muon* muon);
    /// fill electron from Bambu
    void FillElectron(const Electron* elec);
    /// fill particle flow jet from Bambu
    void FillJet(const PFJet* jet);
    /// fill photon from Bambu
    void FillPhoton(const Photon* pho);
    /// fill primary vertex from Bambu
    void FillPV(const Vertex* pv);
    /// match offline reconstructed muon to muon on HLT level
    ULong64_t MatchHLT(const double eta, const double phi);
    /// match offline reconstructed muon to muon on HLT level
    ULong64_t MatchHLT(const double pt, const double eta, const double phi);
    /// check whether this electron is classified to come from a photon conversion according to MVA
    bool IsConversion(const Electron* elec);
    /// compute delta beta corrected particle flow isolation for muon 
    float computePFMuonIso(const Muon* muon, const double dRMax);
    /// compute delta beta corrected particle flow isolation for electron 
    float computePFElecIso(const Electron* elec, const double dRMax);
    /// old particle flow isolation. The isoType can be:
    /// 1: (Charged)Hadron+Electron+Muon; 2: (Neutral)Hadron; 3: Photon
    /// Return value the summed pt. (Deprecated?)
    double PFIsoNoPileup(const ChargedParticle* p, PFCandidateCol* pfNoPileup, double ptMin, double extRadius, double intRadius, int isoType);
    /// old naive particle flow isolation (Deprecated?)
    double PileupIsolation(const Particle* p, PFCandidateCol* pfPileup, double ptMin, double extRadius, double intRadius);
    /// separate PFCandidate collection into PileUp and NoPileUp (re-implementing PF2PAT)
    void separatePileUp(const PFCandidateCol* pfCandidates, const Vertex* pv, const VertexCol* primVerts, PFCandidateOArr* pfPileup, PFCandidateOArr* pfNoPileup, bool checkClosestZVertex = kTRUE);
    /// does this muon fullfil the loose Muon Id?
    bool looseMuId(const Muon * imu);
    /// does this electron fullfil the loose electron Id?
    bool looseEleId(const Electron *iElectron);

  protected:

    /*
      Names for input collection in Bambu
    */

    /// output file
    TFile* fOutputFile;           
    /// name of the output
    TString fOutputName;     
    /// name of the MC particle collection in Bambu
    TString fPartName;
    /// name of the MC event info in Bambu
    TString fMCEvtInfoName;
    /// name of the muon collection in Bambu
    TString fMuonName;
    /// name of the electron collection in Bambu
    TString fElectronName;
    /// name of hte primary vertex collection in Bambu 
    TString fPrimVtxName;
    /// name of the beam spot branch in Bambu
    TString fBeamSpotName;
    /// name of the particle flow jet collection in Bambu
    TString fPFJetName;
    /// name of the photon collection in Bambu
    TString fPhotonName;
    /// name of the trigger mask table in Bambu
    TString fTrigMaskName;
    /// name of the particle flow MET collection in Bambu
    TString fPFMetName;
    /// name of the conversion collection in Bambu
    TString fConversionName;
    /// name of pile-up info branch in Bambu
    TString fPileupName;
    /// name of the energy density rho for pileup corrections with fastjet in Bambu
    TString fPUEnergyDensityName;
    /// name of the particle flow candidate collection in Bambu
    TString fPFCandidateName;
    /// name of the embedded weight collection for the embedding sample in Bambu
    TString fEmbedWeightName;
    
    /*
      Pointer to input collection in Bambu
    */

    /// generator particles
    const MCParticleCol* fParticles;
    /// MC information       
    const MCEventInfo* fMCEvtInfo;
    /// generator jets 
    const GenJetCol* fGenJets;
    /// muons
    const MuonCol* fMuons;
    /// electrons
    const ElectronCol* fElectrons;
    /// primary vertices
    const VertexCol* fPrimVerts;
    /// beam spot
    const BeamSpotCol* fBeamSpot;
    /// particle flow jets
    const PFJetCol* fPFJets;
    /// photons
    const PhotonCol* fPhotons;
    /// trigger mask
    const TriggerMask* fTrigMask;
    /// particle flow MET
    const PFMetCol* fPFMet;
    /// photon conversions
    const DecayParticleCol* fConversions;
    /// pileup information
    const PileupInfoCol* fPileup;
    /// energy density rho for fastjet pileup corrections
    const PileupEnergyDensityCol* fPUEnergyDensity;
    /// particle flow candidates
    const PFCandidateCol* fPFCandidates;
    /// weights for embedded sample(s)
    const EmbedWeightCol* fEmbedWeight;
    
    /*
      Members for the output to the bacon
    */

    /// pfPileUp collection including the primary vertex association along the z-axis   
    PFCandidateOArr* fPFPileUp;
    /// pfNoPileUp collection including the primary vertex association along the z-axis   
    PFCandidateOArr* fPFNoPileUp;
    /// pfPileUp collection w/o the primary vertex association along the z-axis (needed?)  
    PFCandidateOArr* fPFPileUpNoZ;
    /// pfPileNOUp collection w/o the primary vertex association along the z-axis (needed?)   
    PFCandidateOArr* fPFNoPileUpNoZ;
    
    /// flag to indicate if processing collision data
    int fIsData;
    /// flag whether to look at generator info
    int fUseGen;
    /// flag whether to print out HLT table
    bool fPrintTable;
    /// flag whether to skip event processing if HLT does not accept
    bool fSkipIfHLTFail;
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
    
    /// output tree
    TTree* fEventTree;
    
    /// selected primary vertex
    const Vertex* fVertex;
    /// general event information
    TEventInfo fEventInfo;
    /// generator information
    TGenInfo fGenInfo;
    /// muons
    TClonesArray *fMuonArr;
    /// electrons
    TClonesArray *fElectronArr;
    /// particle flow jets
    TClonesArray *fPFJetArr;
    /// photons
    TClonesArray *fPhotonArr;
    /// valid primary vertices
    TClonesArray *fPVArr;
    /// input information for svfit
    TClonesArray *fSVfitEMuArr;
    
    /// names of triggers relevant for the analysis
    std::vector<TString> fTriggerNamesv;
    /// corresponding trigger bits
    std::vector<ULong64_t> fTriggerIdsv;
    /// trigger object names 
    std::vector<TString> fTriggerObjNames1v;
    /// trigger object Ids
    std::vector<ULong64_t> fTriggerObjIds1v;
    /// minimal pt cut for trigger matching
    std::vector<double> fTriggerObjMinPt1v;
    /// trigger object names
    std::vector<TString> fTriggerObjNames2v;
    /// trigger object Ids
    std::vector<ULong64_t> fTriggerObjIds2v;
    /// minimal pt cut for trigger matching
    std::vector<double> fTriggerObjMinPt2v;
    
    /// list of jet correction parameters corresponding to a given jet correction level
    std::vector<TString> fJetCorrParsv;
    /// factorized jet corrector 
    FactorizedJetCorrector* fJetCorrector;
    /// jet correction uncertainties
    JetCorrectionUncertainty* fJetUnc;
    /// Met significance as input for svfit    
    MetSignificance* fMetSignificance;
    /// electron tools
    ElectronTools* fEleTools;
    /// do selection based on JSON files (this could be a member function, which checks whether the list of JSON file is empty or not, no?)
    bool fhasJSON;
    /// list JSON files to be applied
    std::vector<TString> fJSONv;
    /// map of certified runs and lumi sections (for internal use)
    RunLumiRangeMap frlrm;
    /// run and lumi information (for internal use)
    RunLumiSet fRunLumiSet;
    /// run and lumi information (for intenral use)
    RunLumiRangeMap::RunLumiPairType fLastRunLumi;
    
    /// root class definition for schema evolution
    ClassDef(HttNtuplerMod,3)
  };
  
  inline void 
  HttNtuplerMod::AddTrigger(const char* name, const ULong64_t id, const char* objName1, const ULong64_t objId1, const double minPt1, const char* objName2, const ULong64_t objId2, const double minPt2)
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
  HttNtuplerMod::AddJetCorr(const char* name) 
  {
    ifstream jetcorrchk;
    jetcorrchk.open(name); assert(jetcorrchk.is_open()); jetcorrchk.close();
    fJetCorrParsv.push_back(name);
  }
  
  inline void 
  HttNtuplerMod::AddJSON(const char* name) 
  {
    ifstream jsonchk;
    jsonchk.open(name); assert(jsonchk.is_open()); jsonchk.close();
    fJSONv.push_back(name);
    cout << "adding " << name << endl;
  }

}

#endif
