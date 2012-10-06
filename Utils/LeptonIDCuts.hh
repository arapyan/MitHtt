#ifndef LEPTONIDCUTS_HH
#define LEPTONIDCUTS_HH

#include "MitHtt/Ntupler/interface/HiggsAnaDefs.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TPFTau.hh"
#include "MitHtt/Common/MyTools.hh"
#include <TFile.h>                  // file handle class
#include <TGraph.h>
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <cassert>
#include <iostream>

enum EWorkingPoint {
        kLoose,
        kMedium,
        kTight,
	kTighter
};


Bool_t passMuonID(const mithep::TMuon *muon);
Bool_t passPFMuonID(const mithep::TMuon *muon);
Bool_t passTightPFMuonID(const mithep::TMuon *muon, Bool_t mutau);
Bool_t passMuonIsoPU(const mithep::TMuon *muon, Bool_t mutau);
Bool_t passMuonIsoPUTauHad(const mithep::TMuon *muon);
Bool_t passEleID(const mithep::TElectron *electron);
Bool_t passLooseEleID(const mithep::TElectron *electron);
Bool_t passEleMVAID(const mithep::TElectron *electron, Double_t mvaValue);
Bool_t pass2012EleMVAID(const mithep::TElectron *electron, EWorkingPoint WP, Bool_t etau);
Bool_t passEleIdVeto(const mithep::TElectron *ele);
Bool_t passEleNonTrigMVA(const mithep::TElectron *electron, EWorkingPoint WP);
Bool_t passEleIso(const mithep::TElectron *electron);
Bool_t passEleIsoPU(const mithep::TElectron *electron, Bool_t xtau);
Bool_t passEleIsoPUTauHad(const mithep::TElectron *electron);
Bool_t isSoftMuon(const mithep::TMuon *muon);
Bool_t isMuonFO(const mithep::TMuon *muon, const Int_t ver=1);
Bool_t isEleFO(const mithep::TElectron *electron, const Int_t ver=1);
Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi);
Bool_t passtauIdMu(const mithep::TPFTau *tau);
Bool_t tauIdElectron(const mithep::TPFTau *tau);
Bool_t tauIdElectronMVA(const mithep::TPFTau *tau);
Bool_t passtautauId(const mithep::TPFTau *tau,Bool_t ele);

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
//--------------------------------------------------------------------------------------------------
Bool_t passPFMuonID(const mithep::TMuon *muon)
{
  if(!(muon->typeBits & kGlobal || muon->typeBits & kTracker))  return kFALSE;
  if(!(muon->matchesPFCand && muon->matchedPFType==3)) return kFALSE;
  if(fabs(muon->dz)> 0.2)   return kFALSE;

  return kTRUE;

}
//--------------------------------------------------------------------------------------------------
Bool_t passTightPFMuonID(const mithep::TMuon *muon,Bool_t mutau)
{
  //std::cout << bool(muon->typeBits & kGlobal) << " -- " << fabs(muon->dz) << " -- " << fabs(muon->d0) << " -- " << muon->matchesPFCand << " -- " << muon->matchedPFType << std::endl; 
  if(!(muon->typeBits & kGlobal))  return kFALSE;

  if(mutau)
    {
      if(fabs(muon->dz)> 0.2)   return kFALSE;
      if(fabs(muon->d0)>0.045)          return kFALSE;
    }
  else
    {
      if(fabs(muon->dz)> 0.1)   return kFALSE;
      if(fabs(muon->d0)>0.02)          return kFALSE;
    }
  if(muon->muNchi2        > 10)    return kFALSE;
  if(muon->nValidHits     < 1)     return kFALSE;
  if(muon->nMatch         < 2)     return kFALSE;
  if(muon->nPixHits       < 1)     return kFALSE;
  if(muon->nTkHits        < 6)    return kFALSE;
  
  if(!(muon->matchesPFCand && muon->matchedPFType==3)) return kFALSE;
  //if(!(muon->matchesPFCand )) return kFALSE;

  return kTRUE;
  
}
//----------------------------------------------------------------------------------------
Bool_t passMuonIsoPU(const mithep::TMuon *muon,Int_t xtau)
{
  Double_t chargedIso = muon->pfIsoCharged;
  Double_t neutralIso = max(muon->pfIsoNeutral + muon->pfIsoGamma - 0.5 * muon->puIso, 0.0);
  
  Double_t totalIso = chargedIso+neutralIso;
  if(xtau == 2) return (totalIso<0.3*(muon->pt));
  
  if(fabs(muon->eta)<1.479 && !xtau) return (totalIso<0.15*(muon->pt));
  return (totalIso<0.10*(muon->pt));
}
//--------------------------------------------------------------------------------------------------
Bool_t passMuonIsoPUTauHad(const mithep::TMuon *muon)
{
  Double_t chargedIso = muon->pfIsoCharged;
  Double_t neutralIso = max(muon->pfIsoNeutral + muon->pfIsoGamma - 0.5 * muon->puIso, 0.0);

  Double_t totalIso = chargedIso+neutralIso;

  return (totalIso<0.10*(muon->pt));
}
//--------------------------------------------------------------------------------------------------
Double_t muonIsoPU(const mithep::TMuon *muon)
{
  Double_t chargedIso = muon->pfIsoCharged;
  Double_t neutralIso = max(muon->pfIsoNeutral + muon->pfIsoGamma - 0.5 * muon->puIso, 0.0);
  
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
  //if(fabs(electron->d0) > 0.2)   return kFALSE;
  if(fabs(electron->dz) > 0.1)    return kFALSE;

  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;

  //Barrel 
  if (fabs(electron->scEta) < 1.479) {
    if (! ( (0==0)
            //&& electron->sigiEtaiEta < 0.01
            //&& fabs(electron->deltaEtaIn) < 0.007
            //&& fabs(electron->deltaPhiIn) < 0.15
            //&& electron->HoverE < 0.12
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
             //&& electron->sigiEtaiEta < 0.03
             //&& fabs(electron->deltaEtaIn) < 0.009
             //&& fabs(electron->deltaPhiIn) < 0.10
             //&& electron->HoverE < 0.10
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
  /*if (fabs(electron->scEta) < 1.479) {
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
  }*/

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
//--------------------------------------------------------------------------------------------------
Bool_t pass2012EleMVAID(const mithep::TElectron *electron, EWorkingPoint WP, Bool_t etau)
{
   // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;
 
  if(etau)
    {
      if(fabs(electron->d0) > 0.045) return kFALSE;
      if(fabs(electron->dz) > 0.2)   return kFALSE;
    }
  else
    {
      if(fabs(electron->d0) > 0.02)   return kFALSE;
      if(fabs(electron->dz) > 0.1)    return kFALSE;
    }
 
  Int_t subdet = 0;
  if (fabs(electron->scEta) < 0.8) subdet = 0;
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
  if(WP==kMedium) {
    if (MVABin == 0) MVACut = 0.925;
    if (MVABin == 1) MVACut = 0.915; 
    if (MVABin == 2) MVACut = 0.965;  
    if(etau)
      {
	if (MVABin == 3) MVACut = 0.925;
	if (MVABin == 4) MVACut = 0.975;
	if (MVABin == 5) MVACut = 0.985;
      }
    else
      {
	if (MVABin == 3) MVACut = 0.905;
	if (MVABin == 4) MVACut = 0.955;
	if (MVABin == 5) MVACut = 0.975;
      }
  }

  if (electron->mvaValID > MVACut) return kTRUE;
  return kFALSE;
}  

////////////////////////////////////////////////////
// VBTF WP95 Electron ID
Bool_t passEleIdVeto(const mithep::TElectron *ele)
{
  if(fabs(ele->scEta) < 1.479)
  {
    // Barrel
    return (fabs(ele->sigiEtaiEta) < 0.01  &&
	    fabs(ele->deltaEtaIn)  < 0.007 &&
	    fabs(ele->deltaPhiIn)  < 0.8   &&
	    fabs(ele->HoverE)      < 0.15);
  }
  else
  {
    // Endcap
    return (fabs(ele->sigiEtaiEta) < 0.03 &&
	    fabs(ele->deltaEtaIn)  < 0.01 &&
	    fabs(ele->deltaPhiIn)  < 0.7  &&
	    fabs(ele->HoverE)      < 999);
  }
}
//--------------------------------------------------------------------------------------------------
Bool_t passEleNonTrigMVA(const mithep::TElectron *electron, EWorkingPoint WP)
{
  Int_t subdet = 0;
  if (fabs(electron->scEta) < 0.8) subdet = 0;
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
  if(WP==kMedium) {
    if (MVABin == 0) MVACut = 0.925;
    if (MVABin == 1) MVACut = 0.915;
    if (MVABin == 2) MVACut = 0.965;
    if (MVABin == 3) MVACut = 0.905;
    if (MVABin == 4) MVACut = 0.955;
    if (MVABin == 5) MVACut = 0.975;
  }

  if (electron->mvaValID > MVACut) return kTRUE;
  return kFALSE;
}
//-------------------------------------------------------------------------------------------------
Bool_t passEleIsoPU(const mithep::TElectron *electron, int xtau)
{
  Double_t chargedIso = electron->pfIsoCharged;
  Double_t neutralIso = max(electron->pfIsoNeutral + electron->pfIsoGamma - 0.5 * electron->puIso, 0.0);

  Double_t totalIso = chargedIso+neutralIso;
  if(xtau == 2) return (totalIso<0.3*(electron->pt));

  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479 && !xtau) {
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
  Double_t neutralIso = max(electron->pfIsoNeutral + electron->pfIsoGamma - 0.5 * electron->puIso, 0.0);

  Double_t totalIso = chargedIso+neutralIso;

  if(totalIso > 0.10*(electron->pt)) return kFALSE;
  return kTRUE;
}
//-------------------------------------------------------------------------------------------------
Double_t eleIsoPU(const mithep::TElectron *electron)
{
  Double_t chargedIso = electron->pfIsoCharged;
  Double_t neutralIso = max(electron->pfIsoNeutral + electron->pfIsoGammaNoZ - 0.5 * electron->puIso, 0.0);

  return chargedIso+neutralIso;
}
//--------------------------------------------------------------------------------------------------
Bool_t isSoftMuon(const mithep::TMuon *muon)
{
  if(fabs(muon->eta) > 2.1)        return kFALSE;
  //if(muon->nTkHits  < 11)  return kFALSE;
  if(fabs(muon->d0) > 0.2) return kFALSE;
  if(fabs(muon->dz) > 0.1) return kFALSE;

  if(!(muon->typeBits & kGlobal)) return kFALSE;  

  //if(!(muon->qualityBits & kTMLastStationAngTight)) return kFALSE;
	  
  //Double_t iso = (muon->trkIso03 + muon->emIso03 + muon->hadIso03)/muon->pt;
  //if(muon->pt>20 && iso<0.1) return kFALSE;

  return kTRUE;
}

//--------------------------------------------------------------------------------------------------
Bool_t isMuonFO(const mithep::TMuon *muon, const Int_t ver)
{
  if(ver<3) if(muon->nTkHits	  < 11)    return kFALSE;
  if(ver<3) if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->muNchi2	  > 10)    return kFALSE;
  if(muon->nMatch 	  < 2)     return kFALSE;
  if(muon->nValidHits	  < 1)     return kFALSE;
  if(muon->ptErr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;
  if(fabs(muon->d0)       > 0.2)   return kFALSE;  
  if(!(muon->typeBits & kGlobal))  return kFALSE;
  if(ver<3) if(!(muon->typeBits & kTracker)) return kFALSE;

  Double_t iso = (muon->trkIso03 + muon->emIso03 + muon->hadIso03)/muon->pt;
  if(ver==1) return (iso<1.0);
  if(ver==2) return (iso<0.4);
  if(ver==3) return (muon->trkIso03/muon->pt<0.2 && muon->emIso03/muon->pt<0.2 && muon->hadIso03/muon->pt<0.2);
  if(ver==4) return kTRUE;
  
  return kFALSE;
}
//--------------------------------------------------------------------------------------------------
Bool_t isEleFO(const mithep::TElectron *electron, const Int_t ver)
{
  if(fabs(electron->dz) > 0.1)    return kFALSE;
  if(fabs(electron->d0) > 0.02)    return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;
  
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {  
    // barrel
    if(ver==2) if(electron->sigiEtaiEta      > 0.01)  return kFALSE;
    if(ver==2) if(fabs(electron->deltaPhiIn) > 0.15)  return kFALSE;
    if(ver==2) if(fabs(electron->deltaEtaIn) > 0.007) return kFALSE;
    if(ver==2) if(electron->HoverE           > 0.12)  return kFALSE;

    if(electron->trkIso03                         > 0.2*(electron->pt)) return kFALSE;
    if(TMath::Max(electron->emIso03-1,Float_t(0)) > 0.2*(electron->pt)) return kFALSE;
    if(electron->hadIso03                         > 0.2*(electron->pt)) return kFALSE;
        
  } else {
    // endcap
    if(ver==2) if(electron->sigiEtaiEta      > 0.03)  return kFALSE;
    if(ver==2) if(fabs(electron->deltaPhiIn) > 0.10)  return kFALSE;
    if(ver==2) if(fabs(electron->deltaEtaIn) > 0.009) return kFALSE;
    if(ver==2) if(electron->HoverE           > 0.10)  return kFALSE;

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
//-----------------------------------------------------------------------------------------------------tau id
Bool_t passtauIdMu(const mithep::TPFTau *tau)
{
  return(tau->hpsDiscriminators & mithep::TPFTau::kLooseEle &&
	 tau->hpsDiscriminators & mithep::TPFTau::kTightMu &&
	 tau->hpsDiscriminators & mithep::TPFTau::kDecayMode);
}
Bool_t passtautauId(const mithep::TPFTau *tau,Bool_t ele)
{
  if(ele)
    return(tau->hpsDiscriminators & mithep::TPFTau::kLooseEle &&
	   tau->hpsDiscriminators & mithep::TPFTau::kLooseMu &&
	   tau->hpsDiscriminators & mithep::TPFTau::kDecayMode);
  else
    return(tau->hpsDiscriminators & mithep::TPFTau::kLooseMu &&
	   tau->hpsDiscriminators & mithep::TPFTau::kDecayMode);
}
Bool_t tauIdElectron(const mithep::TPFTau *tau)
{
  return(tau->hpsDiscriminators & mithep::TPFTau::kMediumEle &&
         tau->hpsDiscriminators & mithep::TPFTau::kLooseMu &&
         tau->hpsDiscriminators & mithep::TPFTau::kDecayMode);
}
Bool_t tauIdElectronMVA(const mithep::TPFTau *tau, Double_t MVAValue)
{
  Double_t TauEta = tau->eta;
  UInt_t   TauSignalPFGammaCands = tau->nSignalPFGammaCands;
  Bool_t   TauHasGsf = tau->hasGsf;

  if(tau->gammaDEta < -50 || tau->gammaDPhi < -50 ) return kFALSE;

  Bool_t pass = ( (tau->nSignalPFChargedHadrCands == 3) ||
		  (fabs(TauEta)<1.5 && TauSignalPFGammaCands == 0 &&               MVAValue > 0.054) ||
		  (fabs(TauEta)<1.5 && TauSignalPFGammaCands >  0 && TauHasGsf  && MVAValue > 0.060) ||
		  (fabs(TauEta)<1.5 && TauSignalPFGammaCands >  0 && !TauHasGsf && MVAValue > 0.054) ||
		  (fabs(TauEta)>1.5 && TauSignalPFGammaCands == 0 &&               MVAValue > 0.060) ||
		  (fabs(TauEta)>1.5 && TauSignalPFGammaCands >  0 && TauHasGsf  && MVAValue > 0.053) ||
		  (fabs(TauEta)>1.5 && TauSignalPFGammaCands >  0 && !TauHasGsf && MVAValue > 0.049) );

  return pass;
}
#endif
