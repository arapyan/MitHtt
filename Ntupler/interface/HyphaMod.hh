#ifndef MITHTT_NTUPLER_HYPHAMOD_H
#define MITHTT_NTUPLER_HYPHAMOD_H

#include "MitAna/TreeMod/interface/BaseMod.h"
#include "MitAna/DataTree/interface/MCEventInfo.h"
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/PileupInfoFwd.h"
#include "MitAna/DataTree/interface/PileupEnergyDensityFwd.h"
#include "MitAna/DataTree/interface/BaseVertex.h"
#include "MitAna/DataTree/interface/TriggerMask.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"

#include "HiggsAnaDefs.hh"
#include "TEventInfo.hh"
#include "TGenInfo.hh"
#include <TClonesArray.h>
#include "TMuon.hh"
#include "TElectron.hh"
#include "TJet.hh"
#include "TPhoton.hh"
#include "TVertex.hh"

#include <vector>

class TTree;
class TFile;
class TString;

namespace mithep
{  
  class HyphaMod : public BaseMod
  {    
    public:
      HyphaMod(const char *name="HyphaMod", const char *title="BAMBU to humus");
      ~HyphaMod();	    
      
      void SetOutputName(const char *f)         { fOutputName = f; }
       
      void SetIsData(const Bool_t flag)         { fIsData = flag; }
      void SetUseGen(const Int_t useGen)        { fUseGen = useGen; }
      void SetSkipIfHLTFail(const Bool_t flag)  { fSkipIfHLTFail = flag; } 
          
      void SetMuPtMin(const Double_t pt)        { fMuPtMin = pt; }
      void SetMuPtMax(const Double_t pt)        { fMuPtMax = pt; }
      void SetMuEtaMin(const Double_t eta)      { fMuEtaMin = eta; }
      void SetMuEtaMax(const Double_t eta)      { fMuEtaMax = eta; }
      void SetEleEtMin(const Double_t pt)       { fEleEtMin = pt; }
      void SetEleEtMax(const Double_t pt)       { fEleEtMax = pt; }
      void SetEleEtaMin(const Double_t eta)     { fEleEtaMin = eta; }
      void SetEleEtaMax(const Double_t eta)     { fEleEtaMax = eta; }
      void SetJetPtMin(const Double_t et)       { fJetPtMin = et; }
      void SetPhotonEtMin(const Double_t et)    { fPhotonEtMin = et; }
      void SetMinNTracksFit(const UInt_t ntrks) { fMinNTracksFit = ntrks; }
      void SetMinNdof(const Double_t ndof)      { fMinNdof = ndof; }
      void SetMaxAbsZ(const Double_t z)         { fMaxAbsZ = z; }
      void SetMaxRho(const Double_t rho)        { fMaxRho = rho; }
      void SetPrintHLT(const Bool_t flag)       { fPrintTable = flag; }
      
      void AddTrigger(const char* name,  UInt_t id, const Double_t minPt=0,
                      const char* firstObjectModuleName  = "", UInt_t firstObjectId=0,
                      const char* secondObjectModuleName = "", UInt_t secondObjectId=0) {
	fTriggerNamesv.push_back(name);
	fTriggerIdsv.push_back(id);
	fTriggerMinPtv.push_back(minPt);
	fFirstTriggerObjectModuleNamesv.push_back(firstObjectModuleName);
        fFirstTriggerObjectIdsv.push_back(firstObjectId);
	fSecondTriggerObjectModuleNamesv.push_back(secondObjectModuleName);
        fSecondTriggerObjectIdsv.push_back(secondObjectId);
      }

    protected:
      void Begin();
      void BeginRun();
      void EndRun();
      void SlaveBegin();
      void SlaveTerminate();
      void Terminate();
      void Process();
      
      // Fill MC info
      void FillGenZ();
      void FillGenH();
      void FillGenW();
      void FillGenWW();
      
      // Fill muon data object
      void FillMuon(const Muon *mu);
      
      // Fill electron data object
      void FillElectron(const Electron *ele);
      
      // Fill jet data object
      void FillJet(const PFJet *jet);
      
      // Fill photon data object
      void FillPhoton(const Photon *pho);
      
      // Fill vertex data object
      void FillPV(const Vertex *pv);
      
      // Match muon to HLT primitive
      UInt_t MatchHLT(const Double_t eta, const Double_t phi);
      UInt_t MatchHLT(const Double_t pt, const Double_t eta, const Double_t phi);
      
      // Check for conversion with MVF
      Bool_t IsConversion(const Electron *ele);
      
      // Compute PF isolation
      Float_t computePFMuonIso(const Muon *muon, const Double_t dRMax);
      Float_t computePFEleIso(const Electron *electron, const Double_t dRMax);
      
      TFile                        *fOutputFile;           // output file handle
      TString                       fOutputName;           // output file name
      
      TString                       fPartName;             // MC particle collection name
      TString                       fMCEvtInfoName;        // MC event info name
      TString                       fMuonName;             // muon collection name
      TString                       fElectronName;         // electron collection name
      TString                       fPrimVtxName;          // primary vertex collection name
      TString                       fBeamSpotName;         // pointer to beam spot branch
      TString                       fPFJetName;            // particle flow jet collection name
      TString                       fPhotonName;           // photon collection name
      TString                       fTrigMaskName;         // trigger mask name
      TString                       fTCMetName;            // track-corrected MET collection name
      TString                       fPFMetName;            // particle flow MET collection name
      TString                       fConversionName;       // conversion collection name
      TString                       fPileupName;           // pile-up info name
      TString                       fPUEnergyDensityName;  // Fastjet correction info name
      TString                       fPFCandidateName;      // particle flow candidates collection name           
      
      const MCParticleCol          *fParticles;       // MC particle collection handle
      const MCEventInfo            *fMCEvtInfo;       // MC event info handle
      const MuonCol                *fMuons;           // muon collection handle
      const ElectronCol            *fElectrons;       // electron collection handle
      const VertexCol              *fPrimVerts;       // primary vertex collection handle
      const BeamSpotCol            *fBeamSpot;        // pointer to beam spot branch
      const PFJetCol               *fPFJets;          // particle flow jet collection handle
      const PhotonCol              *fPhotons;         // photon collection handle
      const TriggerMask            *fTrigMask;        // trigger mask handle
      const MetCol                 *fTCMet;           // track-corrected MET handle
      const PFMetCol               *fPFMet;           // particle flow MET handle
      const DecayParticleCol       *fConversions;     // conversion collection handle
      const PileupInfoCol          *fPileup;          // pile-up info handle
      const PileupEnergyDensityCol *fPUEnergyDensity; // Fastjet correction info handle
      const PFCandidateCol         *fPFCandidates;    // particle flow candidates collection handle 
      
      Bool_t                  fIsData;          // flag to indicate if processing collision data
      Int_t                   fUseGen;          // flag whether to look at generator info
      Bool_t                  fPrintTable;      // flag whether to print out HLT table
      Bool_t                  fSkipIfHLTFail;   // flag whether to skip event processing if HLT does not accept
       
      Double_t                fMuPtMin;         // minimum reco muon pT
      Double_t                fMuPtMax;         // maximum reco muon pT
      Double_t                fMuEtaMin;        // minimum reco muon eta
      Double_t                fMuEtaMax;        // maximum reco muon eta
      Double_t                fEleEtMin;        // minimum electron supercluster ET
      Double_t                fEleEtMax;        // maximum electron supercluster ET
      Double_t                fEleEtaMin;       // minimum electron supercluster eta
      Double_t                fEleEtaMax;       // maximum electron supercluster eta
      Double_t                fJetPtMin;        // minimum jet ET
      Double_t                fPhotonEtMin;     // minimum photon supercluster ET
      UInt_t                  fMinNTracksFit;   // minimum number of tracks used for a good primary vertex
      Double_t                fMinNdof;         // minimum degrees of freedom for a good primary vertex
      Double_t                fMaxAbsZ;         // maximum z displacement for a good primary vertex
      Double_t                fMaxRho;          // maximum transverse displacement for a good primary vertex 
      
      TTree*                  fEventTree;       // event tree

      BaseVertex              fVertex;          // best primary vertex in the event
            
      TEventInfo              fEventInfo;       // general event information
      TGenInfo                fGenInfo;         // generator information
      TClonesArray           *fMuonArr;         // muon array
      TClonesArray           *fElectronArr;     // electron array
      TClonesArray           *fPFJetArr;        // particle flow jet array
      TClonesArray           *fPhotonArr;       // photon array
      TClonesArray           *fPVArr;           // valid primary vertex array
      
      vector<TString>         fTriggerNamesv;       // names of triggers we're interested in 
      vector<UInt_t>          fTriggerIdsv;         // corresponding ETriggerBit value
      vector<Double_t>        fTriggerMinPtv;       // minimum pT threshold for trigger object
      vector<TString>         fFirstTriggerObjectModuleNamesv; // 
      vector<TString>         fSecondTriggerObjectModuleNamesv; // 
      vector<UInt_t>          fFirstTriggerObjectIdsv;   // ETriggerObjectBit
      vector<UInt_t>          fSecondTriggerObjectIdsv;  // ETriggerObjectBit


      FactorizedJetCorrector *fJetCorrector;    // CMSSW class to handle jet corrections
      
    ClassDef(HyphaMod,1)
  };
}

#endif
