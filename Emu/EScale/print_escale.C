void print_escale()
{

  TFile data("data-EnergyScale.root");
  // cout << "data: " << ScalePlus->GetY()[0] << " " << ResPlus->GetY()[0] << endl;
  TGraph *sc=0; data->GetObject("ScalePlus",sc);
  TGraph *res=0; data->GetObject("ResPlus",res);
  ScalePlus->Print();
  ResPlus->Print();
  // TFile dataEE("data-EE-EnergyScale.root");
  // cout << "data EE: " << ScalePlus->GetY()[0] << " " << ResPlus->GetY()[0] << endl;
  // TFile mcEB("mc-EB-EnergyScale.root");
  // cout << "mc   EB: " << ScalePlus->GetY()[0] << " " << ResPlus->GetY()[0] << endl;
  TFile mc("mc-EnergyScale.root");
  TGraph *scmc=0; mc->GetObject("ScalePlus",scmc);
  TGraph *resmc=0; mc->GetObject("ResPlus",resmc);
  // cout << "mc  : " << ScalePlus->GetY()[0] << " " << ResPlus->GetY()[0] << endl;
  ScalePlus->Print();
  ResPlus->Print();

  // sc->Draw("AlP");
  // scmc->SetLineColor(kBlack);
  // scmc->Draw("lP");
  // res->Draw("AlP");
  // resmc->SetLineColor(kBlack);
  // resmc->Draw("lP");

  for(UInt_t bin=0;bin<sc->GetN();bin++) {
    cout << ((sc->GetY())[bin]) << " " << ((scmc->GetY())[bin]) << endl;
  }

}
