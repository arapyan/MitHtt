#ifndef LEPTONIDCUTS_HH
#define LEPTONIDCUTS_HH

#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "Common/MyTools.hh"
#include <TFile.h>                  // file handle class
#include <TGraph.h>
#include <vector>                   // STL vector class
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings
#include <cassert>
#include <iostream>

Bool_t passMuonID(const mithep::TMuon *muon);
Bool_t passEleID(const mithep::TElectron *electron);

Bool_t isSoftMuon(const mithep::TMuon *muon);

Bool_t isMuonFO(const mithep::TMuon *muon, const Int_t ver=1);
Bool_t isEleFO(const mithep::TElectron *electron);

Double_t projectedMET(const Double_t met, const Double_t metPhi, const Double_t lepPhi);

Double_t trigEff(ETriggerBit trig, mithep::TMuon *muon, mithep::TElectron *ele);


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
Bool_t passMuonID(const mithep::TMuon *muon)
{
  if(muon->nTkHits	  < 11)    return kFALSE;
  if(muon->nPixHits	  < 1)     return kFALSE;
  if(muon->muNchi2	  > 10)    return kFALSE;
  if(muon->nMatch 	  < 2)     return kFALSE;
  if(muon->nValidHits	  < 1)     return kFALSE;
  if(muon->ptErr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.1)   return kFALSE;  
  if(!(muon->typeBits & kGlobal))  return kFALSE;
  if(!(muon->typeBits & kTracker)) return kFALSE;

/*
  Double_t iso = (muon->trkIso03 + muon->emIso03 + muon->hadIso03)/muon->pt;

  if(muon->pt>20) 
    return (iso<0.15 && fabs(muon->d0)<0.02);
  
  return (iso<0.10 && fabs(muon->d0)<0.01);
//*/
//*
  if(fabs(muon->d0)>0.02)         return kFALSE;

  if(muon->pt>20) {
    if(fabs(muon->eta)<1.479) return (muon->pfIso03<0.13*(muon->pt));
    else                      return (muon->pfIso03<0.09*(muon->pt));
  } else {
    if(fabs(muon->eta)<1.479) return (muon->pfIso03<0.06*(muon->pt));
    else                      return (muon->pfIso03<0.05*(muon->pt));
  }
//*/
}
//--------------------------------------------------------------------------------------------------
Bool_t passEleID(const mithep::TElectron *electron)
{
  if(fabs(electron->d0) > 0.02)   return kFALSE;
  if(fabs(electron->dz) > 0.1)    return kFALSE;
  
  // conversion rejection
  if(electron->nExpHitsInner > 0) return kFALSE;
  if(electron->isConv)            return kFALSE;
     
  // barrel/endcap dependent requirments      
  if(fabs(electron->scEta)<1.479) {
    // barrel
//    Double_t iso = (electron->trkIso03 + TMath::Max(electron->emIso03-1,Float_t(0)) + electron->hadIso03)/electron->pt;
//    if(iso > 0.1) return kFALSE;
    if(electron->pfIso04 > 0.13*(electron->pt)) return kFALSE;
     
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
    // endcap
//    Double_t iso = (electron->trkIso03 + electron->emIso03 + electron->hadIso03)/electron->pt;
//    if(iso > 0.1) return kFALSE;
    if(electron->pfIso04 > 0.09*(electron->pt)) return kFALSE;
     
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
  if(muon->ptErr/muon->pt > 0.1)   return kFALSE;
  if(fabs(muon->dz)       > 0.2)   return kFALSE;
/*
if(muon->muNchi2    > 10)        return kFALSE;
if(muon->nMatch     < 2)         return kFALSE;
if(muon->nValidHits < 1)         return kFALSE;
if(!(muon->typeBits & kTracker)) return kFALSE;
if(!(muon->typeBits & kGlobal))  return kFALSE;
*/
  Bool_t isGlobal  = (muon->typeBits & kGlobal) && (muon->muNchi2 < 10) && (muon->nMatch > 1) && (muon->nValidHits > 0);
  Bool_t isTracker = (muon->typeBits & kTracker) && (muon->qualityBits & kTMLastStationTight);
  if(!isGlobal && !isTracker) return kFALSE;

  if(fabs(muon->d0) > 0.1) return kFALSE;
  
  if(ver==1) return (muon->pfIso03/muon->pt<1.0);
  if(ver==2) return (muon->pfIso03/muon->pt<0.4);
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
    if(electron->HoverE	          > 0.12)  return kFALSE;

    if(electron->trkIso03                         > 0.2*(electron->pt)) return kFALSE;
    if(TMath::Max(electron->emIso03-1,Float_t(0)) > 0.2*(electron->pt)) return kFALSE;
    if(electron->hadIso03                         > 0.2*(electron->pt)) return kFALSE;
        
  } else {
    // endcap
    if(electron->sigiEtaiEta	  > 0.03)  return kFALSE;
    if(fabs(electron->deltaPhiIn) > 0.10)  return kFALSE;
    if(fabs(electron->deltaEtaIn) > 0.009) return kFALSE;
    if(electron->HoverE	          > 0.10)  return kFALSE;

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

using namespace std;
class TriggerEfficiency
{
public:
  TriggerEfficiency() { initEff(); }
  double trigEff(float pt1, float eta1, float pt2, float eta2, TString name1, TString name2);

protected:
  void initEff();
  void getEffGraphs(vector<TGraph*> &graphv, vector<float> &etaminv, vector<float> &etamaxv,  const char *subname);

  vector<const char*>      subnamev; // names of sub trigger elements
  vector<vector<TGraph*> > graphvv;  // vector of vector of efficiency graphs
  vector<vector<float> >   etaminvv; // eta range for each graph
  vector<vector<float> >   etamaxvv;

};

void TriggerEfficiency::getEffGraphs(vector<TGraph*> &graphv, vector<float> &etaminv, vector<float> &etamaxv,  const char *subname)
{
  if(!getenv("CMSSW_BASE")) {
    printf("error! TriggerEfficiency called without input files. Define CMSSW_BASE or add by hand.\n");
    assert(0);
  }
  TFile *infile = 0;
  if(TString(subname).Contains("Mu"))
    infile = TFile::Open("$CMSSW_BASE/src/MitHtt/Emu/Selection/data/Trig_Efficiencies_muon.root"); 
  else
    infile = TFile::Open("$CMSSW_BASE/src/MitHtt/Emu/Selection/data/Trig_Efficiencies_emutrig_electron.root");
  assert(infile);
  TIter next(infile->GetListOfKeys());
  for(UInt_t ikey=0;ikey<fabs(infile->GetNkeys());ikey++) {
    string name(next()->GetName());
    if(name.find(subname)  == string::npos) continue; // wrong subname
    if(name.find("Pt")     != string::npos) continue;
    if(name.find("alleta") != string::npos) continue;
    int ietamin = 15;
    int ietamax = name.find("to") + 2;
    string etaminstr(name.substr(ietamin,ietamax-2-ietamin));
    string etamaxstr(name.substr(ietamax,name.find_last_of("_")-ietamax));
    stringstream ssetamin(etaminstr);
    stringstream ssetamax(etamaxstr);
    float etamin=0,etamax=0;
    ssetamin >> etamin;
    ssetamax >> etamax;
    TGraph *gr=0;
    infile->GetObject(name.c_str(),gr); assert(gr);

    graphv.push_back(gr);
    etaminv.push_back(etamin);
    etamaxv.push_back(etamax);
  }

  infile->Close();
  return;
}

void TriggerEfficiency::initEff()
{
  subnamev.push_back("Mu8");
  subnamev.push_back("Mu15");
  subnamev.push_back("Ele8");
  subnamev.push_back("Ele17");

  for(UInt_t isub=0;isub<subnamev.size();isub++) {
    graphvv.push_back(  *(new vector<TGraph*>));
    etaminvv.push_back( *(new vector<float>));
    etamaxvv.push_back( *(new vector<float>));
    getEffGraphs(graphvv.back(),etaminvv.back(),etamaxvv.back(),subnamev[isub]);

    if(graphvv.back().size() < 5) cout << "error! fewer than five graphs found." << endl;
  }

  return;
}  

double TriggerEfficiency::trigEff(float pt1, float eta1, float pt2, float eta2, TString name1, TString name2)
{
  vector<float> ptv;   ptv.push_back(pt1);   ptv.push_back(pt2);
  vector<float> etav; etav.push_back(eta1); etav.push_back(eta2);
  vector<TString> twosubnamev; twosubnamev.push_back(name1); twosubnamev.push_back(name2);

  vector<Int_t> isubv; isubv.push_back(-1);isubv.push_back(-1); // indices of the subnames we want
  for(UInt_t isub=0;isub<subnamev.size();isub++) {
    if(twosubnamev[0]==subnamev[isub])   isubv[0] = isub;
    if(twosubnamev[1]==subnamev[isub])   isubv[1] = isub;
  }
  if(isubv[0]==-1 || isubv[1]==-1) {printf("error! subname not found\n"); return -1;}

  vector<float> effs;
  
  for(UInt_t ihalf=0;ihalf<isubv.size();ihalf++) {
    vector<TGraph*> graphv  =  graphvv[isubv[ihalf]];
    vector<float>   etaminv = etaminvv[isubv[ihalf]];
    vector<float>   etamaxv = etamaxvv[isubv[ihalf]];
    Bool_t found = kFALSE;
    for(UInt_t i=0;i<graphv.size();i++) {
      if(found) break;
      if((etav[ihalf]>etaminv[i] && etav[ihalf]<etamaxv[i]) || etav[ihalf]==etamaxv[i] || etav[ihalf]==etaminv[i]) {
	// cout << graphv[i]->GetName() << "---> ";
	// graphv[i]->Print();
        Double_t *x = graphv[i]->GetX();
        Double_t *y = graphv[i]->GetY();
        for(Int_t ipt=0;ipt<graphv[i]->GetN();ipt++) {
	  if(ptv[ihalf]<x[ipt]) {
	    // cout << x[ipt] << " " << y[ipt] << endl;
	    effs.push_back(y[ipt]);
	    found = kTRUE;
	    break;
	  }
	  if(ipt==graphv[i]->GetN()-1) { // fell throught to the end
	    // cout << "   pt too large: " << ptv[ihalf] << "  setting eff. to last bin: " << y[ipt] << endl;
	    effs.push_back(y[ipt]);
	    found = kTRUE;
	  }
        }
      }
    }

    if(!found) {
      printf("eff. not found! return 0: %s%20.15f%20.15f\n",twosubnamev[ihalf].Data(),ptv[ihalf],etav[ihalf]);
      return 0;
    }
  }

  return effs[0]*effs[1];
  
}

// void play() {

//   TriggerEfficiency effic;
//   cout << endl << effic.trigEff(27,-2.5,188,1.9,"Mu8","Ele17") << endl;
  
// }
#endif
