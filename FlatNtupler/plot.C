#include <TROOT.h>
#include <TSystem.h>
//#include "plotLTau.C"
#include "plotTauTau.C"
//#include "plotTauTauMSSM.C"
//#include "plotLTau7TeV.C"
//#include "plotLTauMSSM.C"
//#include "plotETauMSSM.C"
//#include "plotETau.C"
void plotCategories(std::string iVar="m_vis",int makecard=0,bool secondjet=0) { 
  TFile *outfile=NULL;
  if(makecard)
    outfile = new TFile("htt_tt.inputs-sm-8TeV.root","RECREATE");
  //std::string cut="(antiele_2 == 1 && iso_1 > 0.5 && iso_2 > 0.5 && abs(eta_1) <2.1 && abs(eta_2)<2.1 && jpt_1 > 50.0  && abs(jeta_1) < 3.0 ";
  //std::string cut="(antiele_2 == 1 && iso_1 > 0.5 && iso_2 > 0.5 && abs(eta_1) <2.1 && jpt_1>50 && abs(eta_2)<2.1  && abs(jeta_1) < 3.0 ";
  //std::string cut="(antiele_2==1 && iso_1 > 0.5 && iso_2 > 0.5 && abs(eta_1) <2.1  && abs(eta_2)<2.1  ";
  //std::string cut="(abs(eta_1) <2.1 && abs(eta_2)<2.1 && ((pt_1>pt_2)*passAntiEleMVA3_2>0.5+(pt_1<pt_2)*passAntiEleMVA3_1>0.5)";
  
  std::string cut="(pt_1>45&&pt_2>45 && abs(eta_1) <2.1 && abs(eta_2)<2.1 && passAntiEleNewWPMVA3_2>0.5";
  
  //std::string cut="(pt_1>24&&pt_2>30 && abs(eta_1) <2.1 && abs(eta_2)<2.3 && abs(dZ_2)<0.2 && iso_1<0.5 && passAntiEleNewWPMVA3_2>1.5 && ((pvz+130./tan(2*TMath::ATan(exp(-1.0*eta_2)))) > 0.5 || (pvz+130./tan(2*TMath::ATan(exp(-1.0*eta_2))))<-1.5)";


  //std::string cut="(abs(eta_1) <2.1 && abs(eta_2)<2.3 && abs(dZ_2)<0.2 && iso_1<0.5 && pt_1>17&&pt_2>30 ";


  //   if(secondjet)
  //cut="(njets>1 && abs(eta_1) <2.1 && abs(eta_2)<2.3 && abs(dZ_2)<0.2 && iso_1<0.5 && byCombinedIsolationDeltaBetaCorrRaw3Hits_2<1.5 ";
  // std::string cut="(abs(eta_1) <2.1 && abs(eta_2)<2.3  && iso_1<0.5 && iso_2>0.795 ";
  //if(secondjet)
  //  cut="(antiele_2 == 1 && iso_1 > 0.5 && iso_2 > 0.5 && jpt_1 > 50  && abs(jeta_1) < 3.0 && abs(eta_1) <2.1 && abs(eta_2)<2.1 && jpt_2>0 ";
  
  //plotTauTauMSSM(iVar,6,makecard,outfile,cut,18252);
  //plotTauTauMSSM(iVar,4,makecard,outfile,cut,18252);
  //plotTauTauMSSM(iVar,5,makecard,outfile,cut,18252);
  

  //plotTauTau(iVar,1,makecard,outfile,cut,7274);
  //plotTauTau(iVar,1,makecard,outfile,cut,19490.6);
  //plotTauTau(iVar,2,makecard,outfile,cut,19490.6);
   //plotTauTau(iVar,0,makecard,outfile,cut,11450);
  //plotTauTau(iVar,1,makecard,outfile,cut,11450);
  
  plotTauTau(iVar,0,makecard,outfile,cut,19712);
  plotTauTau(iVar,1,makecard,outfile,cut,19712);
  plotTauTau(iVar,2,makecard,outfile,cut,19712);
  plotTauTau(iVar,3,makecard,outfile,cut,19712);
 
  //plotLTau(iVar,10,makecard,outfile,cut,19712);
  //plotLTau(iVar,0,makecard,outfile,cut,19712);
  //plotLTau(iVar,1,makecard,outfile,cut,19712);
  //plotLTau(iVar,2,makecard,outfile,cut,19712);
  //plotLTau(iVar,3,makecard,outfile,cut,19712);
  //plotLTau(iVar,4,makecard,outfile,cut,19712);
  //plotLTau(iVar,5,makecard,outfile,cut,19712);
  //plotLTau(iVar,6,makecard,outfile,cut,19712);
  //plotLTau(iVar,7,makecard,outfile,cut,19712);
  // plotLTau(iVar,10,makecard,outfile,cut,19800);


  //plotLTau7TeV(iVar,10,makecard,outfile,cut,4935);
  //plotLTau7TeV(iVar,1,makecard,outfile,cut,4935);
  //plotLTau7TeV(iVar,1,makecard,outfile,cut,19712);
  //plotLTau7TeV(iVar,2,makecard,outfile,cut,4935);
  //plotLTau7TeV(iVar,3,makecard,outfile,cut,4935);
  //plotLTau7TeV(iVar,4,makecard,outfile,cut,4935);
  //plotLTau7TeV(iVar,5,makecard,outfile,cut,4935);
  //plotLTau7TeV(iVar,6,makecard,outfile,cut,4935);
  //plotLTau7TeV(iVar,7,makecard,outfile,cut,4935);


  
  //plotLTauMSSM(iVar,8,makecard,outfile,cut,19800);
  //plotLTauMSSM(iVar,9,makecard,outfile,cut,19800);
  //plotLTauMSSM(iVar,10,makecard,outfile,cut,19800);
  
  //plotETauMSSM(iVar,8,makecard,outfile,cut,19800);
  // plotETauMSSM(iVar,9,makecard,outfile,cut,19800);
  //plotETauMSSM(iVar,10,makecard,outfile,cut,19800);

  
  //plotETau(iVar,10,makecard,outfile,cut,19712);
  //plotETau(iVar,0,makecard,outfile,cut,19712);
  //plotETau(iVar,1,makecard,outfile,cut,19712);
  //plotETau(iVar,2,makecard,outfile,cut,19712);
  //plotETau(iVar,3,makecard,outfile,cut,19712);
  //plotETau(iVar,4,makecard,outfile,cut,19712);
  //plotETau(iVar,5,makecard,outfile,cut,19712);
  //plotETau(iVar,6,makecard,outfile,cut,19712);
  //plotETau(iVar,7,makecard,outfile,cut,19712);
  //plotETau(iVar,10,makecard,outfile,cut,19800);


  //plotTauTau(iVar,7,makecard,outfile,cut,18361);
  //plotTauTau(iVar,3,makecard,outfile,cut,18361);
  //plotTauTau(iVar,4,makecard,outfile,cut,18361);
  //plotTauTau(iVar,1,makecard,outfile,cut,18361);
  //plotETau(iVar,0,makecard,outfile,cut,19800);
  //plotETau(iVar,0,makecard,outfile,cut,19800);
  //plotETau(iVar,1,makecard,outfile,cut,19800);
  //plotETau(iVar,2,makecard,outfile,cut,19800);
  //plotLTau(iVar,3,makecard,outfile,cut,19800);
  //plotETau(iVar,4,makecard,outfile,cut,19800);
  //plotTauTau(iVar,6,makecard,outfile,cut,19800);
  //plotTauTau(iVar,3,makecard,outfile,cut,19800);
  //plotTauTau(iVar,4,makecard,outfile,cut,19800);
  //plotETau(iVar,5,makecard,outfile,cut,19400);
  //plotLTau(iVar,5,makecard,outfile,cut,18400);
  //plotLTau(iVar,6,makecard,outfile,cut,18400);
  //plotLTau(iVar,7,makecard,outfile,cut,19400);
  // plotLTau(iVar,7,makecard,outfile,cut,19400);
  //plotETau(iVar,0,makecard,outfile,cut,19400);
  //plotETau(iVar,1,makecard,outfile,cut,19400);
  //7541.8
  //11985
  //19445
}
void plot(TString iDir,int makecard) {  
  gSystem->mkdir(iDir,true);
  gSystem->cd(iDir);
  //plotCategories("m_sv",makecard);
  //plotCategories("pzetavis",makecard);
  //plotCategories("pzetamvamiss",makecard);
  plotCategories("mvis",makecard);
  //plotCategories("pvz+130./tan(2*TMath::ATan(exp(-1.0*eta_2)))",makecard);
  //plotCategories("pt_1*(1/pt_2)",makecard);
  //plotCategories("pthmva",makecard);
  //plotCategories("drll",makecard);
  //plotCategories("mvamet",makecard);
  //plotCategories("npv");
  //plotCategories("pt_1",makecard);
  //plotCategories("btagArray.fArray[0]",makecard);
  //plotCategories("pt_1",makecard);
  //plotCategories("pt_2",makecard);
  //plotCategories("m_1",makecard);
  //plotCategories("m_2",makecard);
  //plotCategories("eta_1",makecard);
  //plotCategories("eta_2",makecard);
  //plotCategories("nprong_1",makecard);
  //plotCategories("nprong_2",makecard);
  //plotCategories("phi_1",makecard);
  //plotCategories("phi_2",makecard);
  //plotCategories("iso_1",makecard);
  //plotCategories("iso_2",makecard);
  //plotCategories("njets",makecard);
  //plotCategories("nbtag",makecard);
  //plotCategories("nbtag",makecard);
  //plotCategories("jpt_1",makecard,1);
  //plotCategories("jpt_2",makecard,1);
  //plotCategories("jeta_1",makecard,1);
  //plotCategories("bpt",makecard,1);
  //plotCategories("beta",makecard,1);
  //plotCategories("jeta_2",makecard,1);
  //plotCategories("jdeta",makecard,1);
  //plotCategories("mjj",makecard,1);
  //plotCategories("mtMVA_1",makecard);
  //plotCategories("mtMVA_2",makecard);
  //plotCategories("clmva",makecard);
  //plotCategories("(abs(jphi_1-phi_1)>3.14)*abs(abs(jphi_1-phi_1)-2.0*3.14)+(abs(jphi_1-phi_1)<3.14)*abs(jphi_1-phi_1)");
  //plotCategories("(abs(jphi_1-phi_1)>3.14)*abs(abs(jphi_1-phi_1)-2.0*3.14)+(abs(jphi_1-phi_1)<3.14)*abs(jphi_1-phi_1)");
  //plotCategories("(abs(jphi_1-phi_2)>3.14)*abs(abs(jphi_1-phi_2)-2.0*3.14)+(abs(jphi_1-phi_2)<3.14)*abs(jphi_1-phi_2)"); 
  //plotCategories("abs(eta_1-jeta_1)");
}
//njet,nbjet,met,pt_1,pt_2,eta_1,eta_2,phi_1,phi_2,drll,pthmva,metphi,npv,jpt_1,eta1jet,jphi,mjj,jdeta,_
