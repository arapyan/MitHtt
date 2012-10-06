#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TClonesArray.h>           // ROOT array class
#include <TCanvas.h>                // class for drawing
#include <TH1D.h>                   // 1D histograms
#include <TH2D.h>
#include <TNtuple.h>
#include <TEfficiency.h>
#include <TVector3.h>
#include <TF1.h>                    // 1D functions
#include <TRandom.h>
#include <TBenchmark.h>             // class to track macro running statistics
#include <vector>                   // STL vector class
#include <map>
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneLikelihood.h"
#include "TauAnalysis/CandidateTools/interface/NSVfitStandaloneAlgorithm.h"
#include "MitHtt/Ntupler/interface/TSVfit.h"
#include "MitHtt/Ntupler/interface/TSVfitter.h"
#endif

using namespace std;
using namespace mithep;

namespace mithep
{
  const Double_t kMuonMass = 105.658369e-3;
  const Double_t kElectronMass = 510.998928e-6;
  TSVfitter svfitter;

  vector<string> tokenize(string str, string delimiter)
  {
    size_t pos = 0, next_pos = 0;
    vector<string> tokens;

    while(next_pos != string::npos)
    {
      next_pos = str.find(delimiter, pos);
      tokens.push_back(str.substr(pos, next_pos - pos));
      pos = next_pos + 1;
    }

    return tokens;
  }

  void setBranches(TTree *tree, string vars, map<string, Double_t> &mvars)
  {
    vector<string> varnames = tokenize(vars, ":");

    for(vector<string>::iterator v = varnames.begin(); v != varnames.end(); v++)
    {
      mvars[*v] = 0;
      tree->SetBranchAddress(v->c_str(), &mvars[*v]);
    }
  }

  double mass(double lep1_pt, double lep1_eta, double lep1_phi, double lep1_m, double lep2_pt, double lep2_eta, double lep2_phi, double lep2_m, double met, double metphi, double cov_00, double cov_01, double cov_10, double cov_11, NSVfitStandalone::kDecayType lep1_decay, NSVfitStandalone::kDecayType lep2_decay)
  {
    if(cov_00 == 0 || cov_10 == 0 || cov_11 == 0) return 0;

    TSVfit svfit;
    svfit.cov_00 = cov_00;
    svfit.cov_01 = cov_01;
    svfit.cov_10 = cov_10;
    svfit.cov_11 = cov_11;

    mithep::FourVectorM dau1, dau2;
    dau1.SetPt(lep1_pt);
    dau1.SetEta(lep1_eta);
    dau1.SetPhi(lep1_phi);
    dau1.SetM(lep1_m);
    dau2.SetPt(lep2_pt);
    dau2.SetEta(lep2_eta);
    dau2.SetPhi(lep2_phi);
    dau2.SetM(lep2_m);

    svfit.daughter1 = dau1;
    svfit.daughter2 = dau2;

    int id = 0;

    if(lep1_decay == NSVfitStandalone::kLepDecay &&
       lep2_decay == NSVfitStandalone::kLepDecay)
    {
      id = 0;
    }
    else if(lep1_decay == NSVfitStandalone::kLepDecay &&
            lep2_decay == NSVfitStandalone::kHadDecay)
    {
      id = 1;
    }
    else if(lep1_decay == NSVfitStandalone::kHadDecay &&
            lep2_decay == NSVfitStandalone::kHadDecay)
    {
      id = 2;
    }
    else
    {
      assert(0);
    }
      
    double svfit_mass = svfitter.integrate(&svfit, met, metphi, id);

    // DEBUG
    // cout << "cov_00 " << covmatrix(0, 0)
    //      << " cov_01 " << covmatrix(0, 1)
    //      << " cov_11 " << covmatrix(1, 1)
    //      << " lep1_pt "  << dau1.Pt()
    //      << " lep1_eta " << dau1.Eta()
    //      << " lep1_phi " << dau1.Phi()
    //      << " lep1_m "   << dau1.M()
    //      << " lep2_pt "  << dau2.Pt()
    //      << " lep2_eta " << dau2.Eta()
    //      << " lep2_phi " << dau2.Phi()
    //      << " lep2_m "   << dau2.M()
    //      << " met " << met
    //      << " metphi " << metphi
    //      << " m " << svfit_mass.M()
    //      << endl;
    // END DEBUG

    return svfit_mass;
  }

  void fillTree(TTree *newTree, map<string, Double_t> &mvars, 
                string lep1_type, string lep2_type, bool zeroonly)
  {
    double lep1_m, lep2_m;
    NSVfitStandalone::kDecayType lep1_decay, lep2_decay;

    if(lep1_type == "electron")
    {
      lep1_m = kElectronMass;
      lep1_decay = NSVfitStandalone::kLepDecay;
    }
    else if(lep1_type == "muon")
    {
      lep1_m = kMuonMass;
      lep1_decay = NSVfitStandalone::kLepDecay;
    }
    else
    {
      lep1_m = mvars["lep1_m"];
      lep1_decay = NSVfitStandalone::kHadDecay;
    }

    if(lep2_type == "electron")
    {
      lep2_m = kElectronMass;
      lep2_decay = NSVfitStandalone::kLepDecay;
    }
    else if(lep2_type == "muon")
    {
      lep2_m = kMuonMass;
      lep2_decay = NSVfitStandalone::kLepDecay;
    }
    else
    {
      lep2_m = mvars["lep2_m"];
      lep2_decay = NSVfitStandalone::kHadDecay;
    }

    if(mvars["m"] == 0 || !zeroonly)
    {
      double svfit_mass = mass(mvars["lep1_pt"],
                               mvars["lep1_eta"],
                               mvars["lep1_phi"],
                               lep1_m,
                               mvars["lep2_pt"],
                               mvars["lep2_eta"],
                               mvars["lep2_phi"],
                               lep2_m,
                               mvars["met"],
                               mvars["metphi"],
                               mvars["mvacov_00"],
                               mvars["mvacov_01"],
                               mvars["mvacov_01"],
                               mvars["mvacov_11"],
                               lep1_decay,
                               lep2_decay);

      double svfit_mass_hi = mass(mvars["lep1_pt"],
                                  mvars["lep1_eta"],
                                  mvars["lep1_phi"],
                                  lep1_m,
                                  mvars["lep2_pt"]*1.03,
                                  mvars["lep2_eta"],
                                  mvars["lep2_phi"],
                                  lep2_m,
                                  mvars["met"],
                                  mvars["metphi"],
                                  mvars["mvacov_00"],
                                  mvars["mvacov_01"],
                                  mvars["mvacov_01"],
                                  mvars["mvacov_11"],
                                  lep1_decay,
                                  lep2_decay);

      double svfit_mass_lo = mass(mvars["lep1_pt"],
                                  mvars["lep1_eta"],
                                  mvars["lep1_phi"],
                                  lep1_m,
                                  mvars["lep2_pt"]*0.97,
                                  mvars["lep2_eta"],
                                  mvars["lep2_phi"],
                                  lep2_m,
                                  mvars["met"],
                                  mvars["metphi"],
                                  mvars["mvacov_00"],
                                  mvars["mvacov_01"],
                                  mvars["mvacov_01"],
                                  mvars["mvacov_11"],
                                  lep1_decay,
                                  lep2_decay);

      mvars["m"] = svfit_mass;
      mvars["mHi"] = svfit_mass_hi;
      mvars["mLo"] = svfit_mass_lo;
    }

    newTree->Fill();
  }

  void processTree(TTree *oldTree, TTree *newTree, string vars,
                   Int_t skipevents, Int_t events,
                   string lep1_type, string lep2_type,
                   bool zeroonly)
  {
    map<string, Double_t> mvars;

    setBranches(oldTree, vars, mvars);
    setBranches(newTree, vars, mvars);

    // Determine number of events to process
    UInt_t totalEvents = oldTree->GetEntries();
    printf("Processing %s\n", oldTree->GetName());
    printf("Sample contains %u events.\n", totalEvents);


    if(events >= 0 && skipevents + events <= (Int_t)totalEvents)
    {
      printf("Limiting analysis to %i.\n", events);
      totalEvents = skipevents + events;
    }
    else if(skipevents >= (Int_t)totalEvents)
    {
      printf("Skipping all events!\n");
    }

    //
    // Loop over events
    //
    for(UInt_t ev = skipevents; ev < totalEvents; ev++)
    {
      if(ev % 100 == 0) printf("--- Processing event %u ---\n", ev);

      oldTree->GetEntry(ev);
      fillTree(newTree, mvars, lep1_type, lep2_type, zeroonly);
    }
  }
}

////////////////
// Main Macro //
////////////////
void applySVFit(const char *input = "mutau_f11-zjets-v14b-pu_test.root",
                const char *output = "mutau_f11-zjets-v14b-pu_svfit.root",
                int skipevents = 0,
                int events = 10,
                string lep1_type = "muon",
                string lep2_type = "tau",
                bool zeroonly = false)
{
  gBenchmark->Start("macro");

  // Input File
  TFile *inFile = TFile::Open(input);
  assert(inFile);

  TTree *ntEvt = (TTree *)inFile->Get("ntEvt");
  assert(ntEvt);

  TTree *ntLooseEvt = (TTree *)inFile->Get("ntLooseEvt");
  assert(ntLooseEvt);

  // Output Tree
  TFile *outFile = new TFile(output, "recreate");
  TTree *ntEvtNew = ntEvt->CloneTree(0);
  TTree *ntLooseEvtNew = ntLooseEvt->CloneTree(0);
  string vars = "m:mHi:mLo:lep1_pt:lep1_eta:lep1_phi:lep2_pt:lep2_eta:lep2_phi:lep2_m:met:metphi:cov_00:cov_01:cov_11:mvacov_00:mvacov_01:mvacov_11";

  // Loop through trees
  processTree(ntEvt, ntEvtNew, vars, skipevents, events, lep1_type, lep2_type, zeroonly);
  processTree(ntLooseEvt, ntLooseEvtNew, vars, skipevents, events, lep1_type, lep2_type, zeroonly);

  // Write trees to file
  cout << "Writing output to " << output << endl;
  ntEvtNew->Write();
  ntLooseEvtNew->Write();

  // Also copy all histograms from the input file to the output file
  /*
    TKey *key;
    TIter it(inFile->GetListOfKeys());

    while((key = dynamic_cast<TKey *>(it.Next())))
    {
    TObject *obj = inFile->Get(key->GetName());

    if(obj == NULL) continue;

    if(obj->InheritsFrom("TH1"))
    {
    // cout << "Copying " << obj->GetName() << endl;
    obj->Write();
    }
    }
  */

  outFile->Close();
  inFile->Close();

  delete outFile;
  delete inFile;

  gBenchmark->Show("macro");
}
