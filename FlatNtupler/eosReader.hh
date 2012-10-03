#include <iostream>
#include <iomanip>
#include <sstream>
#include "TChain.h"
#include "TFile.h"

TChain* loadReader(std::string iDir,std::string iName,int event) {
  TChain *lEvents = 0;
  if(event)  
    lEvents = new TChain("lEvents");
  else
    lEvents = new TChain("Events");
  for(int i0 = 0; i0 < 1000; i0++) {
    std::stringstream pSS; pSS << "root://eoscms//eos/cms"<< iDir  << iName << "/" << iName << "_" << std::setw(4) << std::setfill('0') << i0 << "_ntuple.root";
    //TFile *pFile = TFile::Open(pSS.str().c_str());
    //if(pFile == 0) continue;
    lEvents->Add(pSS.str().c_str());
  }
  std::cout << "loaded : " << iDir << "/" << iName << " --- Total Events ==> " << lEvents->GetEntries() << std::endl;
  return lEvents;
}

//void test() {
//  TChain *lTree = loadReader("/store/cmst3/user/pharris/production/028/","s12-zllm50-2-v9");
//}
