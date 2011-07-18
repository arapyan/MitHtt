#include "../NYStyle/test/NYStyle.h"

//   1  mpar0       -2.79958e-01   1.42387e-02   2.25611e-02  -5.60209e-02
//   2  mpar1       -1.09243e-02   5.64497e-03   7.00850e-04  -2.18486e-03
//   3  mpar2       -1.36237e+00   4.53507e-01   6.57210e-02  -2.75963e-01
//   4  mpar3        9.36308e-03   5.39088e-03   7.14352e-04   1.87262e-03
//   5  mpar4       -1.27907e+00   5.93203e-01   7.51255e-02  -2.58690e-01

void setMCPhi(TF1 *iF1,TF1 *iF2) {
  iF2->SetLineColor(kRed-7); iF1->SetLineWidth(3);
  iF1->SetLineColor(kBlue-7); iF2->SetLineWidth(3);
  iF1->SetParameter(0,-0.0109);
  iF2->SetParameter(0,0.00936);
  iF1->SetParameter(1,-1.279);
  iF2->SetParameter(1,-1.362);
  iF1->SetParameter(2,0.2799/91.3+1.);
  iF2->SetParameter(2,0.2799/91.3+1.);
 }
//   1  mpar0       -3.65278e-01   1.95033e-02  -2.10932e-02   2.18458e-02
//   2  mpar1       -5.17170e-03   9.47048e-03  -1.51808e-03   1.18615e-02
//   3  mpar2       -5.47163e-01   1.21826e+00  -2.68024e-01   2.85952e-01
//   4  mpar3        9.46163e-03   1.10365e-03  -1.48123e-03   1.38066e-03
//   5  mpar4        7.79331e-01   1.00701e+00  -1.53949e-01   1.64046e-01

//Data=>(fabs(iM)+0.37)/sqrt((1-0.00548*sin(iPhi1-0.5399))*(1+0.00937*sin(iPhi2+0.8056))); ///-0.0168 -0.681  1.54 1.66
void setDataPhi(TF1 *iF1,TF1 *iF2) {
  iF2->SetLineColor(kRed-7); iF1->SetLineWidth(3);
  iF1->SetLineColor(kBlue-7); iF2->SetLineWidth(3);
  iF1->SetParameter(0,-0.00547);
  iF2->SetParameter(0, 0.00946);
  iF1->SetParameter(1,-0.5172);
  iF2->SetParameter(1,+0.7793);
  iF1->SetParameter(2,0.37/91.3+1.);
  iF2->SetParameter(2,0.37/91.3+1.);
}
//   1  mpar0        5.30310e-02   3.88876e-02   2.58437e-02   1.06064e-02
//   2  mpar1       -9.24456e-03   1.33329e-03   5.30351e-04  -1.84891e-03
//   3  mpar2       -6.64266e-04   4.59133e-04   3.48520e-04  -1.32853e-04
//   4  mpar3        9.98830e-03   9.24758e-04   6.42868e-07   1.99766e-03
//   5  mpar4        1.43595e-04   4.39903e-04   4.24064e-05   2.87189e-05
//   6  mpar5        1.55671e-03   4.97718e-04   2.06109e-04   3.11342e-04
//   7  mpar6       -1.83657e-03   3.73729e-04   2.06383e-04  -3.67315e-04

//   1  mpar0       -3.40231e-01   2.92058e-02   6.91254e-04  -6.80988e-02
//   2  mpar1       -9.94123e-03   6.20209e-04   1.38480e-05  -1.98825e-03
//   3  mpar2       -9.59127e-04   3.22232e-04   8.53044e-06  -1.91825e-04
//   4  mpar3        1.05421e-02   5.66423e-04   1.44072e-05   2.10841e-03
//   5  mpar4        1.75044e-04   3.62295e-04   9.48697e-06   3.50088e-05
//   6  mpar5        1.78525e-03   2.01813e-04   4.76186e-06   3.57050e-04
//   7  mpar6       -2.09196e-03   2.08223e-04   5.60849e-06  -4.18392e-04

//   1  mpar0        6.95806e-02   3.44816e-02   7.11987e-04   1.39166e-02
//   2  mpar1       -9.33691e-03   8.21511e-04   1.51606e-05  -1.86738e-03
//   3  mpar2       -6.44017e-04   4.14640e-04   1.00012e-05  -1.28803e-04
//   4  mpar3        9.84519e-03   8.93495e-04   1.51356e-05   1.96904e-03
//   5  mpar4        2.28202e-04   4.17152e-04   9.88083e-06   4.56405e-05
//   6  mpar5        1.55733e-03   3.13774e-04   5.82024e-06   3.11465e-04
//   7  mpar6       -1.78372e-03   3.27399e-04   5.81250e-06  -3.56744e-04

void setDataEta(TF1 *iF1,TF1 *iF2) {
  iF2->SetLineColor(kRed-7); iF1->SetLineWidth(3);
  iF1->SetLineColor(kBlue-7); iF2->SetLineWidth(3);
  iF1->SetParameter(0,0.05/91.3+1);
  iF2->SetParameter(0,0.05/91.3+1);
  iF1->SetParameter(1,-0.00933691);
  iF2->SetParameter(1, 0.00984519);
  iF1->SetParameter(2,-0.000644017);
  iF2->SetParameter(2, 0.000228202);
  iF1->SetParameter(3, 0.00155733);
  iF2->SetParameter(3,-0.00178372);
}
//   1  mpar0       -2.78531e-01   2.92677e-02   2.65133e-02  -5.57350e-02
//   2  mpar1        1.41613e-03   7.75707e-04   5.41385e-04   2.83226e-04
//   3  mpar2        1.24622e-04   3.43319e-04   3.60777e-04   2.49245e-05
//   4  mpar3       -1.91784e-03   8.15336e-04   3.57637e-06  -3.83567e-04
//   5  mpar4        4.15290e-04   3.36938e-04   3.63808e-04   8.30580e-05
//   6  mpar5       -1.89737e-03   2.94334e-04   2.18543e-04  -3.79475e-04
//   7  mpar6        2.03956e-03   3.03364e-04   2.83126e-07   4.07912e-04
void setMCEta(TF1 *iF1,TF1 *iF2) {
  iF2->SetLineColor(kRed-7); iF1->SetLineWidth(3);
  iF1->SetLineColor(kBlue-7); iF2->SetLineWidth(3);
  iF1->SetParameter(0,0.27/91.3+1);
  iF2->SetParameter(0,0.27/91.3+1);
  iF1->SetParameter(1, 0.00141613);
  iF2->SetParameter(1,-0.00191784);
  iF1->SetParameter(2, 0.000124622);
  iF2->SetParameter(2, 0.000415290);
  iF1->SetParameter(3,-0.00189737);
  iF2->SetParameter(3, 0.00203956);
}


//   8  s1par0       1.29895e+00   1.97207e-02   2.34508e-02  -8.33383e-01
//   9  s1par1       3.51482e-03   4.92370e-02   1.62733e-02   7.02964e-04
//  10  s1par2      -1.19386e-02   2.07310e-02   1.09138e-02  -2.38773e-03
//  11  s1par3      -1.43569e-02   1.94919e-02   6.37825e-03  -2.87138e-03
//  12  s1par4       3.67112e-02   6.05739e-03   3.69476e-03   7.34232e-03

//   8  s1par0       1.40042e+00   1.70463e-02   1.46760e-04  -8.03682e-01
//   9  s1par1       8.54439e-02   2.67302e-02   1.68671e-02   1.70896e-02
//  10  s1par2      -9.10097e-02   1.47824e-02   1.11748e-02  -1.82029e-02
//  11  s1par3      -4.72824e-02   1.06480e-02   6.87880e-03  -9.45662e-03
//  12  s1par4       5.66641e-02   4.22902e-03   3.90164e-03   1.13331e-02

void setDataRes(TF1 *iF1,TF1 *iF2) {
  iF1->SetLineColor(kRed-7);  iF1->SetLineWidth(3);
  iF2->SetLineColor(kBlue-7); iF2->SetLineWidth(3);
  iF1->SetParameter(0, 1.29895);
  iF1->SetParameter(1, 0.00351482);
  iF1->SetParameter(2,-0.01193);
  iF1->SetParameter(3,-0.0143569);
  iF1->SetParameter(4, 0.0367);
  iF2->SetParameter(0, 1.40042);
  iF2->SetParameter(1, 0.0854439);
  iF2->SetParameter(2,-0.0910087);
  iF2->SetParameter(3,-0.0472824);
  iF2->SetParameter(4, 0.056664);
}

void plotResults(std::string iName="Data8Phi.root") {
  Prep();
  TFile *lMFile = new TFile("MC8Eta.root");
  TGraphErrors *lGM1 = (TGraphErrors*) lMFile->FindObjectAny("ResPlus"); lGM1->SetLineColor(kBlue); lGM1->SetMarkerColor(kBlue); lGM1->SetLineStyle(1);

  TFile *lFile = new TFile(iName.c_str());
  TGraphErrors *lGP0 = (TGraphErrors*) lFile->FindObjectAny("ScalePlus");
  TGraphErrors *lGM0 = (TGraphErrors*) lFile->FindObjectAny("ScaleMinus");
  TGraphErrors *lGP1 = (TGraphErrors*) lFile->FindObjectAny("ResPlus");
  //TGraphErrors *lGM1 = (TGraphErrors*) lFile->FindObjectAny("ResMinus");  

  TF1 *lF1 = new TF1("f1","[0]*sin(x+[1])+[2]",-10,10); 
  TF1 *lF2 = new TF1("f2","[0]*sin(x+[1])+[2]",-10,10); 

  TF1 *lF3 = new TF1("f3","[0]+[1]*x+[2]*x*x+[3]*x*x*x",-10,10);
  TF1 *lF4 = new TF1("f4","[0]+[1]*x+[2]*x*x+[3]*x*x*x",-10,10);

  TF1 *lF5 = new TF1("f5","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-10,10);
  TF1 *lF6 = new TF1("f6","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x",-10,10);

  double *lX = new double[5]; lX[0] = lGP0->GetXaxis()->GetXmin()*0.95;  lX[1]=lGP0->GetXaxis()->GetXmax()*0.95; lX[2]=lGP0->GetXaxis()->GetXmax()*0.95; lX[3]=lGP0->GetXaxis()->GetXmin()*0.95;
  double *lY = new double[5]; lY[0] = 1.000;                        lY[1]=1.000;                       lY[2]=0.992;                       lY[3]=0.992;
  TGraph *lG0 = new TGraph(4,lX,lY); lG0->SetFillColor(kGreen-10);  lG0->SetFillStyle(3000); 

  setDataPhi(lF1,lF2);
  //setMCEta(lF3,lF4);
  TLine *lLine = new TLine(lGP0->GetXaxis()->GetXmin()*0.95,1,lGM0->GetXaxis()->GetXmax()*0.95,1); lLine->SetLineWidth(1);  lLine->SetLineStyle(kDashed); 
  TLegend *lL = new TLegend(0.3,0.2,0.8,0.4); lL->SetFillColor(0); lL->SetBorderSize(0);
  lL->AddEntry(lF1,"Fitted plus trend","l");
  lL->AddEntry(lF2,"Fitted minus trend","l");
  lL->AddEntry(lGP0,"positive muon","lp");
  lL->AddEntry(lGM0,"negative muon","lp");
  lL->AddEntry(lG0 ,"uncertainty band","f");

  //lGP0->GetXaxis()->SetTitle("#eta");
  TCanvas *lC0 = new TCanvas("C0","C0",800,600); lC0->cd(); 
  lGP0->Draw("ape");
  //lG0->Draw("lF");
  lLine->Draw("sames");
  lGP0->Draw("pe");
  lGM0->Draw("pe");
  lF1->Draw("same");
  lF2->Draw("sames");
  lL->Draw();
  //return;
 
  setDataRes(lF5,lF6);
  double *lX1 = new double[100]; 
  double *lY1 = new double[100]; 
  for(int i0 = 0;  i0 < 50; i0++) {lX1[i0]=(lGP0->GetXaxis()->GetXmax()-lGP0->GetXaxis()->GetXmin())/50*i0+lGP0->GetXaxis()->GetXmin(); lY1[i0]=sqrt(lF6->Eval(lX1[i0])**2+0.25);}
  for(int i0 = 50; i0 > 0;  i0--) {lX1[50+50-i0]=(lGP0->GetXaxis()->GetXmax()-lGP0->GetXaxis()->GetXmin())/50*i0+lGP0->GetXaxis()->GetXmin(); lY1[50+50-i0]=sqrt(lF6->Eval(lX1[i0])**2-0.25);}
  TGraph *lG1 = new TGraph(100,lX1,lY1); lG1->SetFillColor(kGreen-10);  lG1->SetFillStyle(1001); 

  TLegend *lL1 = new TLegend(0.3,0.2,0.8,0.4); lL1->SetFillColor(0); lL1->SetBorderSize(0);
  lL1->AddEntry(lF5 ,"Fitted data","l");
  lL1->AddEntry(lF6 ,"Fitted simulation","l");
  lL1->AddEntry(lGP1,"Data","lp");
  lL1->AddEntry(lGM1,"Simulation","lp");
  lL1->AddEntry(lG1 ,"uncertainty band","f");

  lGP1->GetXaxis()->SetTitle("#eta");
  lGP1->GetYaxis()->SetTitle("Resolution(GeV/c^{2})");
  TCanvas *lC1 = new TCanvas("C1","C1",800,600); lC1->cd();
  lGP1->Draw("ape");
  lG1->Draw("F");
  lGP1->Draw("pe");
  lGM1->Draw("pe");
  lF5->Draw("same");
  lF6->Draw("sames");
  lL1->Draw();
}
