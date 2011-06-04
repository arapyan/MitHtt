{
//   TFile *_file0 = TFile::Open("/mnt/hadoop/cmsprod/filefi/020/p11-zttm20-v1g1-pu/58806FDA-6F4F-E011-B1B2-002618FDA211.root");
//   TFile *_file0 = TFile::Open("/mnt/hadoop/cmsprod/filefi/020/r11a-mueg-pr-v1/52BD2EA0-0454-E011-8952-0030486733B4.root");
  TFile *_file0 = TFile::Open("/mnt/hadoop/cmsprod/filefi/020/r11a-mueg-pr-v2/7EDD64CE-0674-E011-9B13-001D09F2462D.root");
  gSystem->Load("/home/dkralph/cms/cmssw/020/CMSSW_4_1_3/lib/slc5_amd64_gcc434/libMitAnaDataTree.so");
  gSystem->Load("/home/dkralph/cms/cmssw/020/CMSSW_4_1_3/lib/slc5_amd64_gcc434/libMitAnaDataCont.so");

  HLT->Print();
  TTree *hlt = (TTree*)_file0->FindObjectAny("HLT");
  TBranch *trigtable = hlt->GetBranch("HLTTriggerTable");
  vector<string> *strv;
  hlt->SetBranchAddress("HLTTriggerTable",&strv);
  hlt->GetEntry(0);
  for(UInt_t i=0;i<strv->size();i++) {
    cout << (*strv)[i] << endl;
  }
}
  
