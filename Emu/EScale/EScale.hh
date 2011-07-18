#include <vector>
#include <sstream>
#include <string>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TProfile.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom1.h"

//
// Reads correections from the specified files for data and MC.
//

using namespace std;

enum { kCenter, kScDown, kScUp, kResDown, kResUp };

class EScale
{
public:
  EScale(TString datafname="$CMSSW_BASE/src/MitHtt/Emu/EScale/865-41x/data-EnergyScale.root",
	 TString mcfname="$CMSSW_BASE/src/MitHtt/Emu/EScale/865-41x/mc-EnergyScale.root");
  Double_t scale(Double_t eta, UInt_t updown);
  Double_t res(Double_t eta, UInt_t updown);
  Double_t pt(Double_t eta, Double_t pt, UInt_t updown);
  void     print();

protected:
  Int_t findbin(Double_t eta);
  
  vector<Double_t> lowedgev;
  TGraphErrors    *data_sc;
  TGraphErrors    *data_res;
  TGraphErrors    *mc_sc;
  TGraphErrors    *mc_res;
  TRandom         *rand;
};
  
EScale::EScale(TString datafname, TString mcfname)
{
  TFile datafile(datafname); assert(datafile.IsOpen());
  datafile.GetObject("ScalePlus",data_sc);
  datafile.GetObject("ResPlus",data_res);
  datafile.Close();

  TFile mcfile(mcfname); assert(mcfile.IsOpen());
  mcfile.GetObject("ScalePlus",mc_sc);
  mcfile.GetObject("ResPlus",mc_res);
  mcfile.Close();

  data_sc->Draw("AlP");
  data_sc->SetLineColor(kBlack);
  mc_sc->Draw("lP");

  assert(data_sc->GetN()==mc_sc->GetN());
  for(Int_t i=0; i<data_sc->GetN(); i++) {
    lowedgev.push_back(data_sc->GetX()[i] - data_sc->GetEX()[i]);
  }

  rand = new TRandom(1024);
}
//----------------------------------------------------------------------------------------
Int_t EScale::findbin(Double_t eta)
{
  Int_t ieta = -1;
  for(UInt_t i=0; i<lowedgev.size(); i++) {
    if((eta >= lowedgev[i]) && (eta <= lowedgev[i]+2*data_sc->GetEX()[i])) {
      ieta = i;
      break;
    }
  }
  if(ieta<0) { cout << "Error: bin not found: " << eta << endl; }

  return ieta;
}
//----------------------------------------------------------------------------------------  
Double_t EScale::scale(Double_t eta, UInt_t updown)
{
  Int_t ieta = findbin(eta);
  if(ieta < 0) return 1;                // eta out of range
  if(mc_sc->GetY()[ieta]==0) return 1;  // mc scale is zero

  Double_t data_delta, mc_delta;
  if(updown==kScDown) {
    data_delta = - data_sc->GetEY()[ieta];
    mc_delta   = - mc_sc->GetEY()[ieta];
  }
  else if(updown==kScUp) {
    data_delta = data_sc->GetEY()[ieta];
    mc_delta   = mc_sc->GetEY()[ieta];
  }
  else
    data_delta = mc_delta = 0;

  return (data_sc->GetY()[ieta] + data_delta)/(mc_sc->GetY()[ieta] + mc_delta);
}
//----------------------------------------------------------------------------------------  
Double_t EScale::res(Double_t eta, UInt_t updown)
{
  Int_t ieta = findbin(eta);
  if(ieta < 0) return 0;          // eta out of range

  Double_t data_delta, mc_delta;
  if(updown==kResDown) {
    data_delta = - data_res->GetEY()[ieta];
    mc_delta   = - mc_res->GetEY()[ieta];
  }
  else if(updown==kResUp) {
    data_delta = data_res->GetEY()[ieta];
    mc_delta   = mc_res->GetEY()[ieta];
  }
  else
    data_delta = mc_delta = 0;

  Double_t data = data_res->GetY()[ieta] + data_delta;
  Double_t mc   = mc_res->GetY()[ieta]   + mc_delta;
  Double_t diff = data*data - mc*mc;

  if(diff<0) return 0;           // currently can't implement if mc wider than data
  return sqrt(diff)/sqrt(2);
}
//----------------------------------------------------------------------------------------
Double_t EScale::pt(Double_t eta, Double_t ptraw, UInt_t updown)
{
  return scale(eta,updown)*ptraw + rand->Gaus(0,res(eta,updown));
}
//----------------------------------------------------------------------------------------
void EScale::print()
{
  printf("%8s%20s%20s\n","eta","scale","res.");
  for(Int_t i=0; i<data_sc->GetN(); i++) {
    Double_t eta = data_sc->GetX()[i];
    Double_t sc  = scale(eta,kCenter);
    Double_t re  = res(eta,kCenter);
    Double_t sc_err  = scale(eta,kScUp) - scale(eta,kCenter);
    Double_t re_err  = res(eta,kResUp)  - res(eta,kCenter);

    printf("%8.3f%13.4f +/-%7.4f%13.2f +/-%5.2f\n",eta,sc,fabs(sc_err),re,fabs(re_err));
  }
}
