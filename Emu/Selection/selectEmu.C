#if !defined(__CINT__) || defined(__MAKECINT__)
#include "Selector.hh"
#endif
void selectEmu(const TString conf,         // input file
               const TString outputDir,    // output directory
	       const Double_t lumi        // luminosity pb^-1
) {
  
  Selector select(conf,outputDir,lumi);
  // Selector select("tmp.conf","/scratch/dkralph/foo",1000,0);
  select.fJetUnc  = kNo;
  select.fBtagEff = kNo;
  select.fMistag  = kNo;
  select.fMuonPt1Min = 20;
  select.fEMuMuonPt2Min  = 10; // min pt for second muon for emu channel
  select.fMuMuMuonPt2Min = 20; // min pt for second muon for mumu channel
  select.fElePt1Min  = 20;
  select.fElePt2Min  = 10;
  select.fJetPtMin   = 30;
  select.fBJetPtMin  = 20;

  
  select.SampleLoop();

}

