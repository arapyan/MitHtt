#ifndef CSAMPLE_HH
#define CSAMPLE_HH

//
// helper class to handle sample inputs
//
class CSample
{
public:
  CSample(){}
  ~CSample(){}

  void print();
  
  TString          label;     // plot item label
  Int_t            color;     // plot item color
  vector<TString>  fnamev;    // ntuple files
  vector<Double_t> xsecv;     // per file cross section
  vector<TString>  jsonv;     // per file JSON file
  vector<Double_t> weightv;   // per file event weight
  map<string,bool> flags;     // flags to control configuration
  
  // data type
  //  0 : MC
  //  1 : mu-el
  //  2 : di-mu
  //  3 : mu
  //  4 : di-el
  //  5 : el
  vector<Int_t> typev;
};

void CSample::print()
{
  cout << endl;
  cout << label.Data() << endl;
  cout << color << endl;
  map<string,bool>::iterator it;
  for(it=flags.begin(); it!=flags.end(); it++)
    cout << setw(12) << (*it).first;
  cout << endl;
  for(it=flags.begin(); it!=flags.end(); it++)
    cout << setw(12) << (*it).second;
  cout << endl;
}

#endif
