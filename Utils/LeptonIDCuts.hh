#ifndef LEPTONIDCUTS_HH
#define LEPTONIDCUTS_HH

#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Common/MyTools.hh"
#include <TFile.h>                  // file handle class
#include <TGraph.h>
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <cassert>
#include <iostream>


Bool_t passMuonID(const mithep::TMuon *muon);
Bool_t passLooseMuonID(const mithep::TMuon *muon);
Bool_t passMuonIso(const mithep::TMuon *muon);
Bool_t passMuonIsoPU(const mithep::TMuon *muon);
Bool_t passMuonIsoPUTauHad(const mithep::TMuon *muon);
Double_t muonIsoPU(const mithep::TMuon *muon);
Bool_t passEleID(const mithep::TElectron *electron);
Bool_t passLooseEleID(const mithep::TElectron *electron);
Bool_t passEleMVAID(const mithep::TElectron *electron, Double_t mvaValue);
Bool_t passEleIso(const mithep::TElectron *electron);
Bool_t passEleIsoPU(const mithep::TElectron *electron);
Bool_t passEleIsoPUTauHad(const mithep::TElectron *electron);
Double_t eleIsoPU(const mithep::TElectron *electron);
Bool_t isSoftMuon(const mithep::TMuon *muon);
Bool_t isMuonFO(const mithep::TMuon *muon, const Int_t ver=1);
Bool_t isEleFO(const mithep::TElectron *electron);
Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi);

//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
Bool_t passMuonID(const mithep::TMuon *muon)
{
  if(fabs(muon->eta) > 2.1)        return kFALSE;

  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->muNchi2	  > 10)    return kFALSE;
  if(muon->nMatch 	  < 2)     return kFALSE;
  if(muon->nValidHits	  < 1)     return kFALSE;
  if(muon->ptErr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;  
  if(!(muon->typeBits & kGlobal))  return kFALSE;
  if(!(muon->typeBits & kTracker)) return kFALSE;

  if(fabs(muon->d0)>0.02)         return kFALSE;

  return kTRUE;

}
//----------------------------------------------------------------------------------------
Bool_t passLooseMuonID(const mithep::TMuon *muon)
{
  return kTRUE;

  // if(muon->nSeg        <  2)                return kFALSE;
  // if(!(muon->typeBits & kStandalone))  return kFALSE;
  // return (muon->trkIso03 < 0.20*muon->pt);
}
//--------------------------------------------------------------------------------------------------
Bool_t passMuonIso(const mithep::TMuon *muon)
{
  if(muon->pt>20) {
    if(fabs(muon->eta)<1.479) return (muon->pfIso03<0.13*(muon->pt));
    else                      return (muon->pfIso03<0.09*(muon->pt));
  } else {
    if(fabs(muon->eta)<1.479) return (muon->pfIso03<0.06*(muon->pt));
    else                      return (muon->pfIso03<0.05*(muon->pt));
  }
}
//--------------------------------------------------------------------------------------------------
Bool_t passMuonIsoPU(const mithep::TMuon *muon)
{
  Double_t chargedIso = muon->pfIsoCharged;
  Double_t neutralIso = max(muon->pfIsoNeutral + muon->pfIsoGamma - 0.5 * muon->puIsoNoZ, 0.0);

  Double_t totalIso = chargedIso+neutralIso;

  if(fabs(muon->eta)<1.479) return (totalIso<0.15*(muon->pt));
  else                      return (totalIso<0.10*(muon->pt));
}
//--------------------------------------------------------------------------------------------------
Bool_t passMuonIsoPUTauHad(const mithep::TMuon *muon)
{
  Double_t chargedIso = muon->pfIsoCharged;
  Double_t neutralIso = max(muon->pfIsoNeutral + muon->pfIsoGamma - 0.5 * muon->puIsoNoZ, 0.0);

  Double_t totalIso = chargedIso+neutralIso;

  return (totalIso<0.10*(muon->pt));
}
//--------------------------------------------------------------------------------------------------
Double_t muonIsoPU(const mithep::TMuon *muon)
{
  Double_t chargedIso = muon->pfIsoCharged;
  Double_t neutralIso = max(muon->pfIsoNeutral + muon->pfIsoGamma - 0.5 * muon->puIsoNoZ, 0.0);
  
  return chargedIso+neutralIso;
}
//--------------------------------------------------------------------------------------------------
Bool_t passEleID(const mithep::TElectron *electron)
{

  if(fabs(electron->scEta) > 2.5) return kFALSE;

  if(fabs(electron->d0) > 0.02)   return kFALSE;
  if(fabs(electron->dz) > 0.1)    return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;

     
  // barrel/endcap dependent requirements      
  if(fabs(electron->scEta)<1.479) {
    //if(electron->pfIso04 > 0.13*(electron->pt)) return kFALSE;

    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.06)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(electron->HoverE	    > 0.04)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.01)  return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.004) return kFALSE;
      if(electron->HoverE	    > 0.025) return kFALSE;
    
    }

  } else {
    //if(electron->pfIso04 > 0.09*(electron->pt)) return kFALSE;
     
    if(electron->pt>20) {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.03)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;
    
    } else {
      if(electron->sigiEtaiEta	    > 0.03)  return kFALSE;
      if(fabs(electron->deltaPhiIn) > 0.02)  return kFALSE;
      if(fabs(electron->deltaEtaIn) > 0.005) return kFALSE;
      if(electron->HoverE	    > 0.10)  return kFALSE;
      
    }
  }
  
  if(electron->pt < 20)
    return ((electron->fBrem>0.15) || (fabs(electron->eta)<1 && electron->EoverP>0.95));
  
  return kTRUE;
}
//--------------------------------------------------------------------------------------------------
Bool_t passLooseEleID(const mithep::TElectron *electron)
{
  if(fabs(electron->d0) > 0.02)   return kFALSE;
  if(fabs(electron->dz) > 0.1)    return kFALSE;

  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;

  //Barrel 
  if (fabs(electron->scEta) < 1.479) {
    if (! ( (0==0)
            && electron->sigiEtaiEta < 0.01
            && fabs(electron->deltaEtaIn) < 0.007
            && fabs(electron->deltaPhiIn) < 0.15
            && electron->HoverE < 0.12
            && (electron->trkIso03) / electron->pt < 0.2
            && (TMath::Max(electron->emIso03 - 1.0, 0.0)) / electron->pt < 0.20
            && (electron->hadIso03) / electron->pt < 0.20

          )
      ) {
      return kFALSE;
    }
  }

  //Endcap
  else {
    if (! (  (0==0)
             && electron->sigiEtaiEta < 0.03
             && fabs(electron->deltaEtaIn) < 0.009
             && fabs(electron->deltaPhiIn) < 0.10
             && electron->HoverE < 0.10
             && (electron->trkIso03 ) / electron->pt < 0.2
             && (TMath::Max(electron->emIso03 - 1.0, 0.0)) / electron->pt < 0.20
             && (electron->hadIso03) / electron->pt < 0.20
  
          )
      ) {
      return kFALSE;
    }
  }
  return kTRUE; 
}
//--------------------------------------------------------------------------------------------------
Bool_t passEleMVAID(const mithep::TElectron *electron, Double_t mvaValue)
{
  if(fabs(electron->d0) > 0.02)   return kFALSE;
  if(fabs(electron->dz) > 0.1)    return kFALSE;

  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;

  // preselection
  // Barrel 
  if (fabs(electron->scEta) < 1.479) {
    if (! ( (0==0)
            && electron->sigiEtaiEta < 0.01
            && fabs(electron->deltaEtaIn) < 0.007
            && fabs(electron->deltaPhiIn) < 0.15
            && electron->HoverE < 0.12
            && (electron->trkIso03) / electron->pt < 0.2
            && (TMath::Max(electron->emIso03 - 1.0, 0.0)) / electron->pt < 0.20
            && (electron->hadIso03) / electron->pt < 0.20

          )
      ) {
      return kFALSE;
    }
  }

  // Endcap
  else {
    if (! (  (0==0)
             && electron->sigiEtaiEta < 0.03
             && fabs(electron->deltaEtaIn) < 0.009
             && fabs(electron->deltaPhiIn) < 0.10
             && electron->HoverE < 0.10
             && (electron->trkIso03 ) / electron->pt < 0.2
             && (TMath::Max(electron->emIso03 - 1.0, 0.0)) / electron->pt < 0.20
             && (electron->hadIso03) / electron->pt < 0.20

          )
      ) {
      return kFALSE;
    }
  }

  Int_t subdet = 0;
  if (fabs(electron->scEta) < 1.0) subdet = 0;
  else if (fabs(electron->scEta) < 1.479) subdet = 1;
  else subdet = 2;
  Int_t ptBin = 0;
  if (electron->pt > 20.0) ptBin = 1;
  Int_t MVABin = -1;
  if (subdet == 0 && ptBin == 0) MVABin = 0;
  if (subdet == 1 && ptBin == 0) MVABin = 1;
  if (subdet == 2 && ptBin == 0) MVABin = 2;
  if (subdet == 0 && ptBin == 1) MVABin = 3;
  if (subdet == 1 && ptBin == 1) MVABin = 4;
  if (subdet == 2 && ptBin == 1) MVABin = 5;

  Double_t MVACut = -9999;
  if (MVABin == 0) MVACut = 0.133;
  if (MVABin == 1) MVACut = 0.465;
  if (MVABin == 2) MVACut = 0.518; 
  if (MVABin == 3) MVACut = 0.942;
  if (MVABin == 4) MVACut = 0.947;
  if (MVABin == 5) MVACut = 0.878 ;

  if (mvaValue > MVACut) return kTRUE;
  return kFALSE;
}
//-------------------------------------------------------------------------------------------------
Bool_t passEleIso(const mithep::TElectron *electron)
{
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    if(electron->pfIso04 > 0.13*(electron->pt)) return kFALSE;
  } else {
    if(electron->pfIso04 > 0.09*(electron->pt)) return kFALSE;
  }
  return kTRUE;
}
//-------------------------------------------------------------------------------------------------
Bool_t passEleIsoPU(const mithep::TElectron *electron)
{
  Double_t chargedIso = electron->pfIsoCharged;
  Double_t neutralIso = max(electron->pfIsoNeutral + electron->pfIsoGamma - 0.5 * electron->puIsoNoZ, 0.0);

  Double_t totalIso = chargedIso+neutralIso;

  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    if(totalIso > 0.15*(electron->pt)) return kFALSE;
  } else {
    if(totalIso > 0.10*(electron->pt)) return kFALSE;
  }
  return kTRUE;
}
//-------------------------------------------------------------------------------------------------
Bool_t passEleIsoPUTauHad(const mithep::TElectron *electron)
{
  Double_t chargedIso = electron->pfIsoCharged;
  Double_t neutralIso = max(electron->pfIsoNeutral + electron->pfIsoGamma - 0.5 * electron->puIsoNoZ, 0.0);

  Double_t totalIso = chargedIso+neutralIso;

  if(totalIso > 0.10*(electron->pt)) return kFALSE;
  return kTRUE;
}
//-------------------------------------------------------------------------------------------------
Double_t eleIsoPU(const mithep::TElectron *electron)
{
  Double_t chargedIso = electron->pfIsoCharged;
  Double_t neutralIso = max(electron->pfIsoNeutral + electron->pfIsoGamma - 0.5 * electron->puIsoNoZ, 0.0);

  return chargedIso+neutralIso;
}
//--------------------------------------------------------------------------------------------------
Bool_t isSoftMuon(const mithep::TMuon *muon)
{
  if(muon->nTkHits  < 11)  return kFALSE;
  if(fabs(muon->d0) > 0.2) return kFALSE;
  if(fabs(muon->dz) > 0.1) return kFALSE;

  if(!(muon->typeBits & kTracker)) return kFALSE;  

  if(!(muon->qualityBits & kTMLastStationAngTight)) return kFALSE;
	  
  Double_t iso = (muon->trkIso03 + muon->emIso03 + muon->hadIso03)/muon->pt;
  if(muon->pt>20 && iso<0.1) return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t isMuonFO(const mithep::TMuon *muon, const Int_t ver)
{
  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->muNchi2	  > 10)    return kFALSE;
  if(muon->nMatch 	  < 2)     return kFALSE;
  if(muon->nValidHits	  < 1)     return kFALSE;
  if(muon->ptErr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;
  if(fabs(muon->d0)       > 0.2)   return kFALSE;  
  if(!(muon->typeBits & kGlobal))  return kFALSE;
  if(!(muon->typeBits & kTracker)) return kFALSE;

  Double_t iso = (muon->trkIso03 + muon->emIso03 + muon->hadIso03)/muon->pt;
  if(ver==1) return (iso<1.0);
  if(ver==2) return (iso<0.4);
  if(ver==3) return (muon->trkIso03/muon->pt<0.2 && muon->emIso03/muon->pt<0.2 && muon->hadIso03/muon->pt<0.2);
  
  return kFALSE;
}
//--------------------------------------------------------------------------------------------------
Bool_t isEleFO(const mithep::TElectron *electron)
{
  if(fabs(electron->dz) > 0.1)    return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;
  
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {  
    // barrel
    if(electron->sigiEtaiEta      > 0.01)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.15)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
    if(electron->HoverE           > 0.12)  return kFALSE;

    if(electron->trkIso03                         > 0.2*(electron->pt)) return kFALSE;
    if(TMath::Max(electron->emIso03-1,Float_t(0)) > 0.2*(electron->pt)) return kFALSE;
    if(electron->hadIso03                         > 0.2*(electron->pt)) return kFALSE;
        
  } else {
    // endcap
    if(electron->sigiEtaiEta      > 0.03)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.10)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.009) return kFALSE;
    if(electron->HoverE           > 0.10)  return kFALSE;

    if(electron->trkIso03 > 0.2*(electron->pt)) return kFALSE;
    if(electron->emIso03  > 0.2*(electron->pt)) return kFALSE;
    if(electron->hadIso03 > 0.2*(electron->pt)) return kFALSE;
  }
    
  return kTRUE;
}
//--------------------------------------------------------------------------------------------------
Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi) 
{
  const Double_t pi = 3.14159265358979;
  Double_t dphi = toolbox::deltaPhi(lepPhi,metPhi);
  if(dphi > 0.5*pi)
    return met;
    
  return met*sin(dphi);
}
#endif
