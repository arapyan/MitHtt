#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>
#include <TBenchmark.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include "MitAna/DataCont/interface/RunLumiRangeMap.h"
#include "MitAna/DataCont/interface/RunLumiSet.h"
#endif

using namespace mithep;

// Main macro function
//--------------------------------------------------------------------------------------------------
void MergeJSON(const TString input) 
{
  // gBenchmark->Start("MergeJSON");

  TString outfilename;          // output of merged json files
  vector<TString> infilenames;  // list of input json files
  
  // 
  // parse input file
  //  
  ifstream ifs;
  ifs.open(input.Data()); 
  assert(ifs.is_open());
  string line;
  getline(ifs,line); 
  outfilename = line;
  outfilename.ReplaceAll("root","json");
  while(getline(ifs,line)) {
    infilenames.push_back(line);
    // use the same input file as for merging ntuples...
    infilenames.back().ReplaceAll("root","json");
  }
  ifs.close();

  //
  // Combine JSON files from each ntuples file
  //
  RunLumiRangeMap rlrm;
  for(UInt_t ifile=0; ifile<infilenames.size(); ifile++) {
    cout << "Adding " << infilenames[ifile] << endl;
    rlrm.AddJSONFile(infilenames[ifile].Data());
  }

  rlrm.DumpJSONFile(outfilename.Data());
  
  std::cout << outfilename << " created!" << std::endl;
  
  // gBenchmark->Show("MergeJSON");
}
