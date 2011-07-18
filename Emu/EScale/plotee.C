#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>                  // access to gROOT, entry point to ROOT system
#include <TFile.h>                  // file handle class
#include <TTree.h>                  // class to access ntuples
#include <TCanvas.h>                // class for drawing
#include <TH1F.h>                   // 1D histograms
#include <TH2F.h>                   // 2D histograms
#include <TBenchmark.h>             // class to track macro running statistics
#include <TVector3.h>               // 3D vector class
#include <TRegexp.h>               // 3D vector class
#include <vector>                   // STL vector class
#include <iostream>                 // standard I/O
#include <iomanip>                  // functions to format standard I/O
#include <fstream>                  // functions for file I/O
#include <string>                   // C++ string class
#include <sstream>                  // class for parsing strings

#include "Common/CPlot.hh"          // helper class for plots
#include "Common/MitStyleRemix.hh"  // style settings for drawing
#include "Common/MyTools.hh"        // miscellaneous helper functions
#include "Common/CSample.hh"        // helper class for organizing input ntuple files

// define structure for output ntuple
#include "MitHtt/Emu/Selection/EmuData.hh"

#endif


//=== FUNCTION DECLARATIONS ======================================================================================

// generate web page
void makeHTML(const TString outDir);


//=== MAIN MACRO =================================================================================================

void plotee(const TString  conf,         // input file
             const TString  ntupleDir,    // directory of input ntuples
	     const TString  outputDir,    // output directory
             const TString  format,       // plot file format
	     const Double_t lumi          // luminosity (pb^-1)
	     ) {  
  gBenchmark->Start("plotee");

  
  //--------------------------------------------------------------------------------------------------------------
  // Settings 
  //==============================================================================================================
  
  vector<TString>  snamev;    // sample name (for output file)  
  vector<CSample*> samplev;   // data/MC samples
    
  //
  // parse .conf file
  //
  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  while(getline(ifs,line)) {
    if(line[0]=='#') continue;
    if(line[0]=='%') { 
      state++; 
      continue; 
    }
    if(line[0]=='$') {
      samplev.push_back(new CSample());
      stringstream ss(line);
      string chr;
      string sname;
      Int_t color;
      ss >> chr >> sname >> color;
      string label = line.substr(line.find('@')+1);
      snamev.push_back(sname);
      samplev.back()->label = label;
      samplev.back()->color = color;
      continue;
    }
    
    if(state==0) {  // define data sample
      string fname;
      string json;
      Int_t type;
      stringstream ss(line);
      ss >> fname >> type >> json;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(type);
      samplev.back()->xsecv.push_back(0);
      samplev.back()->jsonv.push_back(json);
    
    } else if(state==1) {  // define MC samples
      string fname;
      Double_t xsec;
      stringstream ss(line);
      ss >> fname >> xsec;
      samplev.back()->fnamev.push_back(fname);
      samplev.back()->typev.push_back(0);
      samplev.back()->xsecv.push_back(xsec);
    }
  }
  ifs.close();

  // get indices of fakes and zmm, so we can add them together later
  UInt_t ifake=9999, izmm=9999;
  for(UInt_t isam=0; isam<snamev.size(); isam++) {
    if(snamev[isam].Contains("fakes")) ifake = isam;
    if(snamev[isam].Contains("zmm"))   izmm  = isam;
  }
  // if(ifake>snamev.size() || izmm>snamev.size()) { cout << "error -- ifake: " << ifake << " izmm: " << izmm << endl << endl; return; }

  CPlot::sOutDir = outputDir + TString("/plots");
  
  enum { kMuMu, kEleEle, kEleMu, kMuEle };  // final state type

  //--------------------------------------------------------------------------------------------------------------
  // Main analysis code 
  //============================================================================================================== 
  
  // const Double_t pi = 3.14159265358979;
  
  Bool_t hasData = (samplev[0]->fnamev.size()>0);
  
  //
  // Set up histograms
  //

  EmuData data;
  Double_t trigeff,rawMet,rawprojvar;
  UInt_t npt20jets;
  const UInt_t kMaxPt20Jets=35;
  TArrayF *btagArray = new TArrayF;
  btagArray->Set(kMaxPt20Jets);
  
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  TTree *eventTree=0; 
  
  for(UInt_t isam=0; isam<samplev.size(); isam++) {
    if((isam==0) && !hasData) continue;
    
    const TString fname = ntupleDir + TString("/") + snamev[isam] + TString("_select.root");
    cout << "Processing " << fname << "..." << endl;   
    infile = new TFile(fname);
    assert(infile); 

    // Get the TTree and set branch address
    eventTree = (TTree*)infile->Get("Events"); assert(eventTree);
    eventTree->SetBranchAddress("Events",&data);
    // extra branches
    eventTree->SetBranchAddress("npt20jets",&npt20jets);
    eventTree->SetBranchAddress("btagArray",&btagArray);
    eventTree->SetBranchAddress("trigeff",&trigeff);
    eventTree->SetBranchAddress("rawMet",&rawMet);
    eventTree->SetBranchAddress("rawprojvar",&rawprojvar);

    // output flat ntuple
    TString flatfname(outputDir+"/"+snamev[isam]+"-flat.root");
    TFile flatfile(flatfname,"recreate");
    TTree flatree("Events","Events");
    Float_t pt1;  flatree.Branch("pt1"    ,&pt1);
    Float_t eta1; flatree.Branch("eta1"   ,&eta1);
    Float_t phi1; flatree.Branch("phi1"   ,&phi1);
    Float_t pt2;  flatree.Branch("pt2"    ,&pt2);
    Float_t eta2; flatree.Branch("eta2"   ,&eta2);
    Float_t phi2; flatree.Branch("phi2"   ,&phi2);
    Float_t mass; flatree.Branch("mass"   ,&mass);

    for(UInt_t ientry=0; ientry<eventTree->GetEntries(); ientry++) {
      eventTree->GetEntry(ientry);

      Double_t wgt = 1;
      if(isam!=0)
	wgt = data.weight*lumi;

      pt1  = data.lpt1;
      eta1 = data.leta1;
      phi1 = data.lphi1;
      pt2  = data.lpt2;
      eta2 = data.leta2;
      phi2 = data.lphi2;
      mass = data.mass;

      flatree.Fill();
    }
    delete infile;
    infile=0, eventTree=0;

    flatree.Print();
    flatfile.Write();
    flatfile.Close();
  }

  //--------------------------------------------------------------------------------------------------------------
  // Make plots
  //==============================================================================================================

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================

  cout << endl;
  cout << "*" << endl;
  cout << "* SUMMARY" << endl;
  cout << "*--------------------------------------------------" << endl;
  cout << endl;
  
  cout << endl;
  cout << " <-> Output saved in " << outputDir << "/" << endl;
  cout << endl;
  
  gBenchmark->Show("plotee");      
}


//=== FUNCTION DEFINITIONS ======================================================================================

//--------------------------------------------------------------------------------------------------
void makeHTML(const TString outDir)
{
  ofstream htmlfile;
  char htmlfname[100];
  sprintf(htmlfname,"%s/plots.html",outDir.Data());
  htmlfile.open(htmlfname);
  htmlfile << "<!DOCTYPE html" << endl;
  htmlfile << "    PUBLIC \"-//W3C//DTD HTML 3.2//EN\">" << endl;
  htmlfile << "<html>" << endl;
  htmlfile << "<head><title>No Need For a Title</title></head>" << endl;
  htmlfile << "<body bgcolor=\"000000\">" << endl;
  htmlfile << "<h3 style=\"text-align:left; color:DD6600;\">No Need For a Title</h3>" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt_mu.png\"><img src=\"plots/pt_mu.png\" alt=\"plots/pt_mu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt_e.png\"><img src=\"plots/pt_e.png\" alt=\"plots/pt_e.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta_mu.png\"><img src=\"plots/eta_mu.png\" alt=\"plots/eta_mu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta_e.png\"><img src=\"plots/eta_e.png\" alt=\"plots/eta_e.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phimu.png\"><img src=\"plots/phimu.png\" alt=\"plots/phimu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phiele.png\"><img src=\"plots/phiele.png\" alt=\"plots/phiele.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/dphi_emu.png\"><img src=\"plots/dphi_emu.png\" alt=\"plots/dphi_emu.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt.png\"><img src=\"plots/pt.png\" alt=\"plots/pt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/lepdeta.png\"><img src=\"plots/lepdeta.png\" alt=\"plots/lepdeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/met.png\"><img src=\"plots/met.png\" alt=\"plots/met.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metraw.png\"><img src=\"plots/metraw.png\" alt=\"plots/metraw.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mt.png\"><img src=\"plots/mt.png\" alt=\"plots/mt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/metdphi.png\"><img src=\"plots/metdphi.png\" alt=\"plots/metdphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetamiss.png\"><img src=\"plots/pzetamiss.png\" alt=\"plots/pzetamiss.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetavis.png\"><img src=\"plots/pzetavis.png\" alt=\"plots/pzetavis.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pzetavar.png\"><img src=\"plots/pzetavar.png\" alt=\"plots/pzetavar.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/rawProjVar.png\"><img src=\"plots/rawProjVar.png\" alt=\"plots/rawProjVar.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetpt1.png\"><img src=\"plots/jetpt1.png\" alt=\"plots/jetpt1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetpt2.png\"><img src=\"plots/jetpt2.png\" alt=\"plots/jetpt2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jeteta1.png\"><img src=\"plots/jeteta1.png\" alt=\"plots/jeteta1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jeteta2.png\"><img src=\"plots/jeteta2.png\" alt=\"plots/jeteta2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/bjetpt.png\"><img src=\"plots/bjetpt.png\" alt=\"plots/bjetpt.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/bjeteta.png\"><img src=\"plots/bjeteta.png\" alt=\"plots/bjeteta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/bjetphi.png\"><img src=\"plots/bjetphi.png\" alt=\"plots/bjetphi.png\" width=\"100%\"></a></td>" << endl;

  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetdphi.png\"><img src=\"plots/jetdphi.png\" alt=\"plots/jetdphi.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/mjj.png\"><img src=\"plots/mjj.png\" alt=\"plots/mjj.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetdeta.png\"><img src=\"plots/jetdeta.png\" alt=\"plots/jetdeta.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/jetEtaProd.png\"><img src=\"plots/jetEtaProd.png\" alt=\"plots/jetEtaProd.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/njets.png\"><img src=\"plots/njets.png\" alt=\"plots/njets.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nbjets.png\"><img src=\"plots/nbjets.png\" alt=\"plots/nbjets.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/btag.png\"><img src=\"plots/btag.png\" alt=\"plots/btag.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/nvertices_reweighted.png\"><img src=\"plots/nvertices_reweighted.png\" alt=\"plots/nvertices_reweighted.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;   
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_lept.png\"><img src=\"plots/vismass_lept.png\" alt=\"plots/vismass_lept.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_incl.png\"><img src=\"plots/vismass_incl.png\" alt=\"plots/vismass_incl.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_novbf.png\"><img src=\"plots/vismass_class_novbf.png\" alt=\"plots/vismass_class_novbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_vbf.png\"><img src=\"plots/vismass_class_vbf.png\" alt=\"plots/vismass_class_vbf.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;  
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_nob.png\"><img src=\"plots/vismass_class_nob.png\" alt=\"plots/vismass_class_nob.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_class_b.png\"><img src=\"plots/vismass_class_b.png\" alt=\"plots/vismass_class_b.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/vismass_lept-high.png\"><img src=\"plots/vismass_lept-high.png\" alt=\"plots/vismass_lept-high.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "<hr />" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt1.png\"><img src=\"plots/pt1.png\" alt=\"plots/pt1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/pt2.png\"><img src=\"plots/pt2.png\" alt=\"plots/pt2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta1.png\"><img src=\"plots/eta1.png\" alt=\"plots/eta1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/eta2.png\"><img src=\"plots/eta2.png\" alt=\"plots/eta2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi1.png\"><img src=\"plots/phi1.png\" alt=\"plots/phi1.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a target=\"_blank\" href=\"plots/phi2.png\"><img src=\"plots/phi2.png\" alt=\"plots/phi2.png\" width=\"100%\"></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "<hr />" << endl;

  htmlfile << "<table border=\"0\" cellspacing=\"5\" width=\"100%\">" << endl;   
  htmlfile << "<tr>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "<td width=\"25%\"><a></a></td>" << endl;
  htmlfile << "</tr>" << endl;
  htmlfile << "</table>" << endl;
  htmlfile << "<hr />" << endl;
            
  htmlfile << "</body>" << endl;
  htmlfile << "</html>" << endl;
  htmlfile.close();
}
