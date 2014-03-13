#include <iostream>
#include <fstream>
#include <algorithm>
#include <string>
#include <sstream>
#include <TH1F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <Rtypes.h>

#include <TMath.h>
#include <TAxis.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TAttLine.h>
#include <TPaveText.h>
#include <TColor.h>

#include "TTree.h"
#include "TCanvas.h"


//#include "HttStyles.h"
#include "HttStyles.cc"

TTree * load(std::string iName) { 
  TFile *lFile = new TFile(iName.c_str());
  TTree *lTree  = (TTree*) lFile->FindObjectAny("Events");
  return lTree;
}

float maximum(TH1F* h, bool LOG=false){
  if(LOG){
    if(h->GetMaximum()>1000){ return 1000.*TMath::Nint(30*h->GetMaximum()/1000.); }
    if(h->GetMaximum()>  10){ return   10.*TMath::Nint(30*h->GetMaximum()/  10.); }
    return 50*h->GetMaximum(); 
  }
  else{
    if(h->GetMaximum()>  12){ return 10.*TMath::Nint((1.20*h->GetMaximum()/10.)); }
    if(h->GetMaximum()> 1.2){ return TMath::Nint((1.50*h->GetMaximum())); }
    //return 1.6*h->GetMaximum(); 
    return 1.7*h->GetMaximum(); 
  }
}

void blind(TH1F* iH,double low,double high) {
  for(int i = 0; i < iH->GetNbinsX(); i++) {
    double mass = iH->GetBinCenter(i);
    if(mass <= high && mass >= low)
      {
	iH->SetBinContent(i,0);
	iH->SetBinError(i,0);
      }
  }
}


void bbtt_tt(std::string var,int nbins, double xmin, double xmax,std::string xtitle, std::string ytitle)
{
  SetStyle(); gStyle->SetLineStyleString(11,"20 10");
  TH1::SetDefaultSumw2(1);

  std::string dir = "/data/blue/Bacon/029a/tautaunew/mhhvar/";
  double sigscale = 10;
  double sigscale1 = 10; 
  std::stringstream scale; scale << sigscale;
  std::stringstream scale1; scale1 << sigscale1;

  //Cut definitions
  double luminosity = 19712;
  std::stringstream lumi; lumi << luminosity;
  std::string objcutini = "(pt_1>45&&pt_2>45 && abs(eta_1) <2.1 && abs(eta_2)<2.1 && passAntiEleNewWPMVA3_2>0.5  && byCombinedIsolationDeltaBetaCorrRaw3Hits_1>1.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2>1.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<10.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<10.0 )";
  std::string objcut = "(pt_1>45&&pt_2>45 && abs(eta_1) <2.1 && abs(eta_2)<2.1 && passAntiEleNewWPMVA3_2>0.5 && byCombinedIsolationDeltaBetaCorrRaw3Hits_1<1.0 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.0)";
  //std::string jetcut = objcut+"*(mjj>100 && mjj<150 && jpt_1>30 && jpt_2>30 && ((jcsv_1>0.244 || jcsv_2>0.244) && (jcsv_1>0.244 && jcsv_2>0.244) && (jcsv_1>0.679 || jcsv_2>0.679) && (jcsv_1>0.679 && jcsv_2>0.679)))";
  //std::string jetcutini = objcutini+"*(mjj>100 && mjj<150 && jpt_1>30 && jpt_2>30 &&  ((jcsv_1>0.244 || jcsv_2>0.244) &&  (jcsv_1>0.244 && jcsv_2>0.244) && (jcsv_1>0.679 || jcsv_2>0.679) &&  (jcsv_1>0.679 && jcsv_2>0.679)))";
     std::string jetcut = objcut+"*(jpt_1>30 && jpt_2>30 && (jcsv_1>0.244 && jcsv_2>0.244))";
   std::string jetcutini = objcutini+"*(jpt_1>30 && jpt_2>30 && (jcsv_1>0.244 && jcsv_2>0.244))";
  //std::string jetcut = objcut+"*(jpt_1>30 && jpt_2>30 && mtMVA_1<30)";
  //signal region
  std::string datacut  = jetcut+"*(q_1*q_2<0)";
  std::string smhhcut = jetcut+"*weight*("+scale1.str()+"*decaymodeweight*0.714e-3/1.972e09)*(q_1*q_2<0)*"+lumi.str();
  std::string h300cut = jetcut+"*weight*("+scale.str()+"*decaymodeweight*0.0792/1.972e09)*(q_1*q_2<0)*"+lumi.str();
  std::string a300cut = jetcut+"*weight*("+scale.str()+"*decaymodeweight*0.1713/1.972e09)*(q_1*q_2<0)*"+lumi.str();
  std::string tthcut = jetcut+"*(weight*decaymodeweight*0.0789)*(q_1*q_2<0)*"+lumi.str();
  std::string vbfcut = jetcut+"*(weight*decaymodeweight*0.0997)*(q_1*q_2<0)*"+lumi.str();
  std::string gfcut = jetcut+"*(weight*decaymodeweight*1.2179)*(q_1*q_2<0)*"+lumi.str();
  std::string zttcut = jetcut+"*weight*(decaymodeweight)*emblid*(genmatch==3)*(q_1*q_2<0)";
  std::string ttbarcut = jetcut+"*1.11*0.96*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  std::string wjetcut = jetcut+"*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  std::string ewkcut = jetcut+"*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  std::string zttmccut = jetcut+"*weight*(decaymodeweight/1.972e09)*(genmatch==3)*(q_1*q_2<0)*"+lumi.str();
  std::string zllcut = jetcut+"*weight*(decaymodeweight)*(abs(genlid_1)<15 && abs(genlid_2)<15 && genmatch==1)*(q_1*q_2<0)*"+lumi.str();
  std::string zljcut = jetcut+"*weight*(decaymodeweight)*(genmatch!=3 && genmatch!=1)*(q_1*q_2<0)*"+lumi.str();
  std::string fakescut = jetcut+"*(q_1*q_2>0)";
  std::string fakescutosini = jetcutini+"*(q_1*q_2<0)";
  std::string fakescutini = jetcutini+"*(q_1*q_2>0)";
  //Ztt normalization
  std::string zttcutloose = objcut+"*weight*(decaymodeweight)*emblid*(genmatch==3)*(q_1*q_2<0)";
  std::string zttmccutloose = objcut+"*weight*(decaymodeweight)*(genmatch==3)*(q_1*q_2<0)*"+lumi.str();
  //QCD SS control region
  std::string zttmccutss = jetcut+"*weight*(decaymodeweight)*(genmatch==3)*(q_1*q_2>0)*"+lumi.str();
   std::string zllcutss = jetcut+"*weight*decaymodeweight*(abs(genlid_1)<15 && abs(genlid_2)<15 && genmatch==1)*(q_1*q_2>0)*"+lumi.str();
  std::string zljcutss = jetcut+"*weight*decaymodeweight*(genmatch!=3 && genmatch!=1)*(q_1*q_2>0)*"+lumi.str();
  std::string ttbarcutss = jetcut+"*1.11*0.96*(q_1*q_2>0)*(weight*decaymodeweight)*"+lumi.str();
  std::string ewkcutss = jetcut+"*(q_1*q_2>0)*(weight*decaymodeweight)*"+lumi.str();
  std::string wjetcutss = jetcut+"*(q_1*q_2>0)*(weight*decaymodeweight)*"+lumi.str();
  std::string zttmccutssini = jetcutini+"*weight*(decaymodeweight)*(genmatch==3)*(q_1*q_2>0)*"+lumi.str();
  std::string zllcutssini = jetcutini+"*weight*decaymodeweight*(abs(genlid_1)<15 && abs(genlid_2)<15 && genmatch==1)*(q_1*q_2>0)*"+lumi.str();
  std::string zljcutssini = jetcutini+"*weight*decaymodeweight*(genmatch!=3 && genmatch!=1)*(q_1*q_2>0)*"+lumi.str();
  std::string ttbarcutssini = jetcutini+"*1.11*0.96*(q_1*q_2>0)*(weight*decaymodeweight)*"+lumi.str();
  std::string ewkcutssini = jetcutini+"*(q_1*q_2>0)*(weight*decaymodeweight)*"+lumi.str();
  std::string wjetcutssini = jetcutini+"*(q_1*q_2>0)*(weight*decaymodeweight)*"+lumi.str();
  //QCD OS control region
   std::string zttmccutosini = jetcutini+"*weight*(decaymodeweight)*(genmatch==3)*(q_1*q_2<0)*"+lumi.str();
  std::string zllcutosini = jetcutini+"*weight*decaymodeweight*(abs(genlid_1)<15 && abs(genlid_2)<15 && genmatch==1)*(q_1*q_2<0)*"+lumi.str();
  std::string zljcutosini = jetcutini+"*weight*decaymodeweight*(genmatch!=3 && genmatch!=1)*(q_1*q_2<0)*"+lumi.str();
  std::string ttbarcutosini = jetcutini+"*1.11*0.96*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  std::string ewkcutosini = jetcutini+"*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  std::string wjetcutosini = jetcutini+"*(q_1*q_2<0)*(weight*decaymodeweight)*"+lumi.str();
  
  //--------------------------------------------------------------------------
  
  //Get the trees
  TTree *datatree = load(dir+"data_select.root"); 
  TTree *embtree = load(dir+"emb_select.root"); 
  TTree *ttbartree = load(dir+"ttbar-8TeV_select.root");
  TTree *wjettree = load(dir+"wjets_select.root");
  TTree *ewktree = load(dir+"ewk-8TeV_select.root");
  TTree *zttmctree = load(dir+"ztt-mad_select.root");
  TTree *tthtree = load(dir+"htt_vtth_sm_125_select.root");
  TTree *gftree = load(dir+"htt_gf_sm_125_select.root");
  TTree *vbftree = load(dir+"htt_vbf_sm_125_select.root");
  TTree *h300tree = load(dir+"s12-H3002h125-2t2b-8tev_select.root");
  TTree *a300tree = load(dir+"s12-A300zh125-2l2b-8tev_select.root");
  TTree *smhhtree = load(dir+"s12-hh125-bbtt-tt-8tev_select.root");
  TTree *zlltree = load(dir+"zll-mad_select.root");
  //-------------------------------------------------------------------------
  
  //Get histograms
  TCanvas *canv0 = MakeCanvas("canv", "histograms", 600, 600);
  canv0->cd();
  std::string vardraw;
  TH1F *data = new TH1F("Data","",nbins,xmin,xmax);
  vardraw = var+">>"+"Data";
  datatree->Draw(vardraw.c_str(),datacut.c_str());
  InitHist(data, xtitle.c_str(), ytitle.c_str()); InitData(data);
  TH1F *Ztt = new TH1F("Embedded","",nbins,xmin,xmax);
  vardraw = var+">>"+"Embedded";
  embtree->Draw(vardraw.c_str(),zttcut.c_str());
  InitHist(Ztt  , xtitle.c_str(), ytitle.c_str(), TColor::GetColor(248,206,104), 1001);
  TH1F *ttbar = new TH1F("TTbar","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbar";
  ttbartree->Draw(vardraw.c_str(),ttbarcut.c_str());
  InitHist(ttbar, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(155,152,204), 1001);
  TH1F *wjets = new TH1F("Wjets","",nbins,xmin,xmax);
  vardraw = var+">>"+"Wjets";
  wjettree->Draw(vardraw.c_str(),wjetcut.c_str());
  InitHist(wjets, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(222,90,106), 1001);
  TH1F *zll = new TH1F("Zll","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zll";
  zlltree->Draw(vardraw.c_str(),zllcut.c_str());
  InitHist(zll, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(222,90,106), 1001);
  TH1F *zlj = new TH1F("Zlj","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zlj";
  zlltree->Draw(vardraw.c_str(),zljcut.c_str());
  InitHist(zlj, xtitle.c_str(), ytitle.c_str(), TColor::GetColor(222,90,106), 1001);
  TH1F *ewk = new TH1F("Ewk","",nbins,xmin,xmax);
  vardraw = var+">>"+"Ewk";
  ewktree->Draw(vardraw.c_str(),ewkcut.c_str());
  InitHist(ewk, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(222,90,106), 1001);
  TH1F *zttmc = new TH1F("Zttmc","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zttmc";
  zttmctree->Draw(vardraw.c_str(),zttmccut.c_str());
  InitHist(zttmc, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(248,206,104), 1001);
  TH1F *fakes = new TH1F("Fakes","",nbins,xmin,xmax);
  vardraw = var+">>"+"Fakes";
  datatree->Draw(vardraw.c_str(),fakescut.c_str());
  InitHist(fakes, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(250,202,255), 1001);
  TH1F *tth = new TH1F("TTH","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTH";
  tthtree->Draw(vardraw.c_str(),tthcut.c_str());
  InitHist(tth, xtitle.c_str(), ytitle.c_str(), kGreen+2, 1001);
  TH1F *vbfh = new TH1F("VBFH","",nbins,xmin,xmax);
  vardraw = var+">>"+"VBFH";
  vbftree->Draw(vardraw.c_str(),vbfcut.c_str());
  TH1F *ggh = new TH1F("GGH","",nbins,xmin,xmax);
  vardraw = var+">>"+"GGH";
  tth->Add(vbfh); tth->Add(ggh);
  gftree->Draw(vardraw.c_str(),gfcut.c_str());
  TH1F *smhh = new TH1F("SMhh","",nbins,xmin,xmax);
  vardraw = var+">>"+"SMhh";
  smhhtree->Draw(vardraw.c_str(),smhhcut.c_str());
  InitSignal(smhh);
  smhh->SetLineColor(kBlack);
  TH1F *h300 = new TH1F("H300","",nbins,xmin,xmax);
  vardraw = var+">>"+"H300";
  h300tree->Draw(vardraw.c_str(),h300cut.c_str());
  InitSignal(h300);
  TH1F *a300 = new TH1F("A300","",nbins,xmin,xmax);
  vardraw = var+">>"+"A300";
  a300tree->Draw(vardraw.c_str(),a300cut.c_str());
  InitSignal(a300);
  a300->SetLineColor(6);
  
  //ztt normalization histograms
  TH1F *zttmcloose = new TH1F("Zttmcloose","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zttmcloose";
  zttmctree->Draw(vardraw.c_str(),zttmccutloose.c_str());
  TH1F *zttloose = new TH1F("Zttloose","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zttloose";
  embtree->Draw(vardraw.c_str(),zttcutloose.c_str());

  //qcd SS control region
  TH1F *ttbarss = new TH1F("TTbarss","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbarss";
  ttbartree->Draw(vardraw.c_str(),ttbarcutss.c_str());
  TH1F *wjetsss = new TH1F("Wjetsss","",nbins,xmin,xmax);
  vardraw = var+">>"+"Wjetsss";
  wjettree->Draw(vardraw.c_str(),wjetcutss.c_str());
  TH1F *ewkss = new TH1F("Ewkss","",nbins,xmin,xmax);
  vardraw = var+">>"+"Ewkss";
  ewktree->Draw(vardraw.c_str(),ewkcutss.c_str());
  TH1F *zttmcss = new TH1F("Zttmcss","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zttmcss";
  zttmctree->Draw(vardraw.c_str(),zttmccutss.c_str());
  TH1F *zllss = new TH1F("Zllss","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zllss";
  zlltree->Draw(vardraw.c_str(),zllcutss.c_str());
  TH1F *zljss = new TH1F("Zljss","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zljss";
  zlltree->Draw(vardraw.c_str(),zljcutss.c_str());
  TH1F *ttbarssini = new TH1F("TTbarssini","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbarssini";
  ttbartree->Draw(vardraw.c_str(),ttbarcutssini.c_str());
  TH1F *wjetsssini = new TH1F("Wjetsssini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Wjetsssini";
  wjettree->Draw(vardraw.c_str(),wjetcutssini.c_str());
  TH1F *ewkssini = new TH1F("Ewkssini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Ewkssini";
  ewktree->Draw(vardraw.c_str(),ewkcutssini.c_str());
  TH1F *zttmcssini = new TH1F("Zttmcssini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zttmcssini";
  zttmctree->Draw(vardraw.c_str(),zttmccutssini.c_str());
  TH1F *zllssini = new TH1F("Zllssini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zllssini";
  zlltree->Draw(vardraw.c_str(),zllcutssini.c_str());
  TH1F *zljssini = new TH1F("Zljssini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zljssini";
  zlltree->Draw(vardraw.c_str(),zljcutssini.c_str());
  TH1F *fakesini = new TH1F("Fakesini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Fakesini";
  datatree->Draw(vardraw.c_str(),fakescutini.c_str());
  //QCD OS control region
  TH1F *ttbarosini = new TH1F("TTbarosini","",nbins,xmin,xmax);
  vardraw = var+">>"+"TTbarosini";
  ttbartree->Draw(vardraw.c_str(),ttbarcutosini.c_str());
  TH1F *wjetsosini = new TH1F("Wjetsosini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Wjetsosini";
  wjettree->Draw(vardraw.c_str(),wjetcutosini.c_str());
  TH1F *ewkosini = new TH1F("Ewkosini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Ewkosini";
  ewktree->Draw(vardraw.c_str(),ewkcutosini.c_str());
  TH1F *zttmcosini = new TH1F("Zttmcosini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zttmcosini";
  zttmctree->Draw(vardraw.c_str(),zttmccutosini.c_str());
  TH1F *zllosini = new TH1F("Zllosini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zllosini";
  zlltree->Draw(vardraw.c_str(),zllcutosini.c_str());
  TH1F *zljosini = new TH1F("Zljosini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Zljosini";
  zlltree->Draw(vardraw.c_str(),zljcutosini.c_str());
  TH1F *fakesosini = new TH1F("Fakesosini","",nbins,xmin,xmax);
  vardraw = var+">>"+"Fakesosini";
  datatree->Draw(vardraw.c_str(),fakescutosini.c_str());
  InitHist(fakesosini, xtitle.c_str(), ytitle.c_str(),  TColor::GetColor(250,202,255), 1001);
  delete canv0;
  //----------------------------------------------------------------------------
  //Background estimations
  //embedded normalization from DY MC
  double sfactor = zttmcloose->Integral()/zttloose->Integral();
  Ztt->Scale(sfactor);
  //QCD SS subtraction
  TH1F* ssbkg = (TH1F*)ttbarss ->Clone("ssbkg");
  ssbkg->Add(wjetsss); ssbkg->Add(ewkss); ssbkg->Add(zttmcss); ssbkg->Add(zllss); ssbkg->Add(zljss);  
  fakes->Add(ssbkg,-1);
  for(int ibin = 0; ibin < fakes->GetNbinsX(); ibin++) if(fakes->GetBinContent(ibin) < 0) fakes->SetBinContent(ibin,0);
  TH1F* ssbkgini = (TH1F*)ttbarssini ->Clone("ssbkgini");
  ssbkgini->Add(wjetsssini); ssbkgini->Add(ewkssini); ssbkgini->Add(zttmcssini);  ssbkgini->Add(zllssini); ssbkgini->Add(zljssini); 
  fakesini->Add(ssbkgini,-1);
  for(int ibin = 0; ibin < fakesini->GetNbinsX(); ibin++) if(fakesini->GetBinContent(ibin) < 0) fakesini->SetBinContent(ibin,0);
  TH1F* ssbkgosini = (TH1F*)ttbarosini ->Clone("ssbkgosini");
  ssbkgosini->Add(wjetsosini); ssbkgosini->Add(ewkosini); ssbkgosini->Add(zttmcosini); ssbkgosini->Add(zllosini); ssbkgosini->Add(zljosini); 
  fakesosini->Add(ssbkgosini,-1);
  for(int ibin = 0; ibin < fakesosini->GetNbinsX(); ibin++) if(fakesosini->GetBinContent(ibin) < 0) fakesosini->SetBinContent(ibin,0);
  double sf = fakes->Integral()/fakesini->Integral();
  std::cout << "QCD SS scale factor  " << sf << std::endl;
  fakes = fakesosini;
  fakes->Scale(sf);
 
  //---------------------------------------------------------------------------
  //Print out the yields
  Double_t error=0.0;
  ofstream outfile;
  outfile.open("yields.txt");
  outfile << "Yields for the signal region." << std::endl;
  outfile << "Data    "  << data->IntegralAndError(0,data->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "SM hh   "  << smhh->IntegralAndError(0,smhh->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "H(300)->hh   "  << h300->IntegralAndError(0,h300->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "A(300)->hh   "  << a300->IntegralAndError(0,a300->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "TTH   "  << tth->IntegralAndError(0,tth->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "Fakes    "  << fakes->IntegralAndError(0,fakes->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "Ztt    "  << Ztt->IntegralAndError(0,Ztt->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "Ztt MC   "  << zttmc->IntegralAndError(0,zttmc->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "Zll    "  << zll->IntegralAndError(0,zll->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "Zlj    "  << zlj->IntegralAndError(0,zlj->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "ttbar    "  << ttbar->IntegralAndError(0,ttbar->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "ewk    "  << ewk->IntegralAndError(0,ewk->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "wjets    "  << wjets->IntegralAndError(0,wjets->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "S/sqrt(B)    "  << h300->Integral()/(wjets->Integral()+Ztt->Integral()+ttbar->Integral()+fakes->Integral()+ewk->Integral()) << endl;
  //--------------------------------------------------------------------------
  //stack the electroweak histtograms
  wjets->Add(zll); wjets->Add(zlj); ewk->Add(wjets);
  //continue outputing
  outfile << "Ewk total    "  << ewk->IntegralAndError(0,ewk->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << endl << endl << endl;
  outfile << "Same sign yields" << endl;
  outfile << "Ztt    "  << zttmcss->IntegralAndError(0,zttmcss->GetNbinsX(),error) << "+/-" << error << endl;
  outfile << "Zll    "  << zllss->IntegralAndError(0,zllss->GetNbinsX(),error) << "+/-" << error << endl;
   outfile << "Zlj    "  << zljss->IntegralAndError(0,zljss->GetNbinsX(),error) << "+/-" << error << endl;
   outfile << "ttbar    "  << ttbarss->IntegralAndError(0,ttbarss->GetNbinsX(),error) << "+/-" << error << endl;
   outfile << "ewk    "  << ewkss->IntegralAndError(0,ewkss->GetNbinsX(),error) << "+/-" << error << endl;
   outfile << "wjets    "  << wjetsss->IntegralAndError(0,wjetsss->GetNbinsX(),error) << "+/-" << error << endl;
   outfile << endl;
   outfile << "In the signal region (100,150GeV)  " <<endl;
   outfile << "H(300)->hh   "  << h300->IntegralAndError(5,6,error) << "+/-" << error << endl;
   outfile << "TTH   "  << tth->IntegralAndError(5,6,error) << "+/-" << error << endl;
   outfile << "Fakes    "  << fakes->IntegralAndError(5,6,error) << "+/-" << error << endl;
   outfile << "Ztt    "  << Ztt->IntegralAndError(5,6,error) << "+/-" << error << endl;
   outfile << "ttbar    "  << ttbar->IntegralAndError(5,6,error) << "+/-" << error << endl;
   outfile << "ewk    "  << ewk->IntegralAndError(5,6,error) << "+/-" << error << endl;
   outfile << "In the signal region (75,100GeV)  " <<endl;
   outfile << "A(300)->hh   "  << a300->IntegralAndError(4,4,error) << "+/-" << error << endl;
   outfile << "TTH   "  << tth->IntegralAndError(4,4,error) << "+/-" << error << endl;
   outfile << "Fakes    "  << fakes->IntegralAndError(4,4,error) << "+/-" << error << endl;
   outfile << "Ztt    "  << Ztt->IntegralAndError(4,4,error) << "+/-" << error << endl;
   outfile << "ttbar    "  << ttbar->IntegralAndError(4,4,error) << "+/-" << error << endl;
   outfile << "ewk    "  << ewk->IntegralAndError(4,4,error) << "+/-" << error << endl;
   outfile.close();
   //-----------------------------------------------------------------------
   //Draw the histograms
   TCanvas *canv = MakeCanvas("canv", "histograms", 600, 600);
   canv->cd();
   ewk->Add(fakes);  ttbar->Add(ewk); 
   //Ztt->Add(ttbar); tth->Add(Ztt);
   zttmc->Add(ttbar); tth->Add(zttmc);
   //Error band stat
   TH1F* errorBand = (TH1F*)tth ->Clone("errorBand");
   errorBand  ->SetMarkerSize(0);
   errorBand  ->SetFillColor(13);
   errorBand  ->SetFillStyle(3013);
   errorBand  ->SetLineWidth(1);
   //  for(int idx=0; idx<errorBand->GetNbinsX(); ++idx){
   //     if(errorBand->GetBinContent(idx)>0){
   //       std::cout << "Uncertainties on summed background samples: " << errorBand->GetBinError(idx)/errorBand->GetBinContent(idx) << std::endl;
   //       break;
   //     }
   //}
   data->SetMaximum(1.2*std::max(maximum(data, 0), maximum(tth, 0)));
   blind(data,75,150);
   data->Draw("e");
   tth->Draw("histsame");
   //Ztt->Draw("histsame");
   zttmc->Draw("histsame");
   ttbar->Draw("histsame");
   ewk->Draw("histsame");
   fakes->Draw("histsame");
   data->Draw("esame");
   errorBand->Draw("e2same");
   //smhh->Draw("histsame");
   a300->Draw("histsame");
   h300->Draw("histsame");
   canv->RedrawAxis();
   //---------------------------------------------------------------------------
   //Adding a legend
   TLegend* leg = new TLegend(0.53, 0.65, 0.95, 0.90);
   SetLegendStyle(leg);
   leg->AddEntry(h300  , TString::Format("%.0f#timesH(300 GeV)#rightarrowhh#rightarrow#tau#tau bb", sigscale) , "L" );
   leg->AddEntry(a300  , TString::Format("%.0f#timesA(300 GeV)#rightarrowZh#rightarrow#tau#tau bb", sigscale) , "L" );
   //leg->AddEntry(smhh , TString::Format("%.0f#timeshh#rightarrowhh#rightarrow#tau#tau bb", sigscale1) , "L" );
   leg->AddEntry(data , "Observed"                       , "LP");
   leg->AddEntry(tth  , "SM H#rightarrow#tau#tau"   , "F" );
   leg->AddEntry(Ztt  , "Z#rightarrow#tau#tau"           , "F" );
   leg->AddEntry(ttbar, "t#bar{t}"                       , "F" );
   leg->AddEntry(ewk  , "Electroweak"                    , "F" );
   leg->AddEntry(fakes, "QCD"                            , "F" );
   leg->AddEntry(errorBand,"bkg. uncertainty","F");
   leg->Draw();
   //---------------------------------------------------------------------------
   
   //CMS preliminary 
   const char* dataset = "CMS Preliminary,  H#rightarrow#tau#tau, 19.7 fb^{-1} at 8 TeV";
   const char* category = "";
   CMSPrelim(dataset, "#tau_{h}#tau_{h}", 0.17, 0.835);
   //CMSPrelim(dataset, "", 0.16, 0.835);
   TPaveText* chan     = new TPaveText(0.52, 0.35, 0.91, 0.55, "tlbrNDC");
   chan->SetBorderSize(   0 );
   chan->SetFillStyle(    0 );
   chan->SetTextAlign(   12 );
   chan->SetTextSize ( 0.05 );
   chan->SetTextColor(    1 );
   chan->SetTextFont (   62 );
   chan->AddText(category);
   chan->Draw();
   //-------------------------------------------------------------------------
   //Save histograms
   canv->Print((var+".png").c_str());
  
   /*
     Ratio Data over MC
   */
   TCanvas *canv1 = MakeCanvas("canv0", "histograms", 600, 400);
   canv1->SetGridx();
   canv1->SetGridy();
   canv1->cd();

   TH1F* model = (TH1F*)Ztt ->Clone("model");
   TH1F* test1 = (TH1F*)data->Clone("test1"); 
   for(int ibin=0; ibin<test1->GetNbinsX(); ++ibin){
     //the small value in case of 0 entries in the model is added to prevent the chis2 test from failing
     model->SetBinContent(ibin+1, model->GetBinContent(ibin+1)>0 ? model->GetBinContent(ibin+1)*model->GetBinWidth(ibin+1) : 0.01);
     //model->SetBinError  (ibin+1, CONVERVATIVE_CHI2 ? 0. : model->GetBinError  (ibin+1)*model->GetBinWidth(ibin+1));
     model->SetBinError  (ibin+1, 0);
     test1->SetBinContent(ibin+1, test1->GetBinContent(ibin+1)*test1->GetBinWidth(ibin+1));
     test1->SetBinError  (ibin+1, test1->GetBinError  (ibin+1)*test1->GetBinWidth(ibin+1));
   }
   double chi2prob = test1->Chi2Test      (model,"PUW");        std::cout << "chi2prob:" << chi2prob << std::endl;
   double chi2ndof = test1->Chi2Test      (model,"CHI2/NDFUW"); std::cout << "chi2ndf :" << chi2ndof << std::endl;
   double ksprob   = test1->KolmogorovTest(model);              std::cout << "ksprob  :" << ksprob   << std::endl;
   double ksprobpe = test1->KolmogorovTest(model,"DX");         std::cout << "ksprobpe:" << ksprobpe << std::endl;  

   std::vector<double> edges;
   TH1F* zero = (TH1F*)ttbar->Clone("zero"); zero->Clear();
   TH1F* rat1 = (TH1F*)data->Clone("rat1"); 
   for(int ibin=0; ibin<rat1->GetNbinsX(); ++ibin){
     rat1->SetBinContent(ibin+1, Ztt->GetBinContent(ibin+1)>0 ? data->GetBinContent(ibin+1)/Ztt->GetBinContent(ibin+1) : 0);
     rat1->SetBinError  (ibin+1, Ztt->GetBinContent(ibin+1)>0 ? data->GetBinError  (ibin+1)/Ztt->GetBinContent(ibin+1) : 0);
     zero->SetBinContent(ibin+1, 0.);
     zero->SetBinError  (ibin+1, Ztt->GetBinContent(ibin+1)>0 ? Ztt ->GetBinError  (ibin+1)/Ztt->GetBinContent(ibin+1) : 0);
   }
   for(int ibin=0; ibin<rat1->GetNbinsX(); ++ibin){
     if(rat1->GetBinContent(ibin+1)>0){
       edges.push_back(TMath::Abs(rat1->GetBinContent(ibin+1)-1.)+TMath::Abs(rat1->GetBinError(ibin+1)));
       // catch cases of 0 bins, which would lead to 0-alpha*0-1
       rat1->SetBinContent(ibin+1, rat1->GetBinContent(ibin+1)-1.);
     }
   }
   float range = 0.1;
   std::sort(edges.begin(), edges.end());
   if (edges[edges.size()-2]>0.1) { range = 0.2; }
   if (edges[edges.size()-2]>0.2) { range = 0.5; }
   if (edges[edges.size()-2]>0.5) { range = 1.0; }
   if (edges[edges.size()-2]>1.0) { range = 1.5; }
   if (edges[edges.size()-2]>1.5) { range = 2.0; }
   rat1->SetLineColor(kBlack);
   rat1->SetFillColor(kGray );
   rat1->SetMaximum(+range);
   rat1->SetMinimum(-range);
   rat1->GetYaxis()->CenterTitle();
   rat1->GetYaxis()->SetTitle("#bf{Data/MC-1}");
   rat1->GetXaxis()->SetTitle("#bf{m_{#tau#tau} [GeV]}");
   rat1->Draw();
   zero->SetFillStyle(  3013);
   zero->SetFillColor(kBlack);
   zero->SetLineColor(kBlack);
   zero->SetMarkerSize(0.1);
   zero->Draw("e2histsame");
   canv1->RedrawAxis();

   TPaveText* stat1 = new TPaveText(0.20, 0.76+0.061, 0.32, 0.76+0.161, "NDC");
   stat1->SetBorderSize(   0 );
   stat1->SetFillStyle(    0 );
   stat1->SetTextAlign(   12 );
   stat1->SetTextSize ( 0.05 );
   stat1->SetTextColor(    1 );
   stat1->SetTextFont (   62 );
   stat1->AddText(TString::Format("#chi^{2}/ndf=%.3f,  P(#chi^{2})=%.3f", chi2ndof, chi2prob));
   //stat1->AddText(TString::Format("#chi^{2}/ndf=%.3f,  P(#chi^{2})=%.3f, P(KS)=%.3f", chi2ndof, chi2prob, ksprob));
   //stat1->Draw();
   canv1->Print((var+"_ratio.png").c_str());
}
