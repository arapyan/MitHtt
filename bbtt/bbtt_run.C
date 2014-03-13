#include "bbtt_mt.C"
//#include "bbtt_et.C"
//#include "bbtt_tt.C"

void bbtt_run(TString dir)
{
  //defintion of cuts
  gSystem->mkdir(dir,true);
  gSystem->cd(dir);
  //mt plotting
  //bbtt_tt("mjj",12,0,300,"m_{jj} [GeV]","Events");
  //bbtt_tt("ptvis",12,0,300,"p_{T}^{#tau#tau} [GeV]","Events");
  //bbtt_tt("pthmva",12,0,300,"p_{T}^{#tau#tau E_{MET}} [GeV]","Events");
  //bbtt_tt("dijetpt",12,0,300,"p_{T}^{jj} [GeV]","Events");
  //bbtt_tt("hdijetphi",15,1,3.14,"#Delta #phi(bb,#tau #tau)","Events");
  //bbtt_tt("mtMVA_1",12,0,300,"m_{T} [GeV]","Events");
  bbtt_mt("m_sv",12,0,300,"m_{#tau#tau} [GeV]","Events");
  

  // bbtt_mt("mt2new",16,0,400,"m_{2T} [GeV]","Events");
  //bbtt_mt("m_sv",12,0,300,"m_{2T} [GeV]","Events");
  //bbtt_et("m_sv",12,0,300,"m_{2T} [GeV]","Events");
  //bbtt_tt("m_H",24,100,700,"m_{H} [GeV]","Events");
  //bbtt_mt("mvis_H",24,100,700,"mvis_{H} [GeV]","Events");
  
  //et plotting
  // bbtt_et("mjj",12,0,300,"m_{jj} [GeV]","Events");
  //bbtt_et("mtMVA_1",12,0,300,"m_{T} [GeV]","Events");
  //bbtt_et("m_sv",12,0,300,"m_{#tau#tau} [GeV]","Events");

  //tt plotting
  //bbtt_tt("mjj",12,0,300,"m_{jj} [GeV]","Events");
  //bbtt_tt("m_sv",12,0,300,"m_{#tau#tau} [GeV]","Events");
  
}
