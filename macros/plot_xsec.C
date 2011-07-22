#include "mssm_xs_tools.h"

/**
   \class   plot_xsec plot_xsec.C "MitHtt/macros/plot_xsec.c"

   \brief   macro to create the cross section for Higgs production as function of mA/mH

   macro to retrieve most up to date cross sections and BR's for MSSM Higgs production and 
   decay. To run the tool the following files must be present in your working directory:

    + mssm_xs_tools.h
    + mssm_xs_tools.C
    + out.mhmax_7.root (expected in sub-directory MitHtt/macros/root)

   You can download these from the following web page: 

    + https://twiki.cern.ch/twiki/bin/view/LHCPhysics/MSSMNeutral


   Run this macro in root using CINT in the following way (sorry mssm_xs_tools.C is NOT 
   C++ compliant :-(...): 

   root -l
   .L .L mssm_xs_tools.C
   .x plot_xsec.C

   a plot will be created with the cross section times BR for ggH / qqH(VBF) for the SM 
   and bbH / ggH using the Santander matching for the MSSM. in addition and a set of 
   numbers will be printed to the screen including various flavour schemes for bbH 
   (MSSM). Sorry for the format, this should be cleaned up some time...

   For the SM cross sections and BR's the values have been taken from: 

    + https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt7TeV
    + https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageBR

   Please check for the root input file to be up to date from tiem to time (last check 
   for EPS 2011). 
) 
*/


int
plot_xsec()
{
  mssm_xs_tools mssm;
  mssm.SetInput("root/out.mhmax_7.root"); //mssm.help();

  /// SM Higgs boson
  TGraph* smHiggs = new TGraph();
  smHiggs->SetPoint( 0, 100., 24.02 * 8.36e-02);
  smHiggs->SetPoint( 1, 110., 19.84 * 8.03e-02);
  smHiggs->SetPoint( 2, 120., 16.63 * 7.11e-02);
  smHiggs->SetPoint( 3, 130., 14.12 * 5.49e-02);
  smHiggs->SetPoint( 4, 140., 12.13 * 3.54e-02);
  smHiggs->SetPoint( 5, 150., 10.50 * 1.79e-02);
  smHiggs->SetPoint( 6, 160., 9.080 * 3.97e-03);
  smHiggs->SetPoint( 7, 170., 7.729 * 9.20e-04);
  smHiggs->SetPoint( 8, 180., 6.739 * 5.87e-04);
  smHiggs->SetPoint( 9, 190., 5.896 * 3.76e-04);
  smHiggs->SetPoint(10, 200., 5.249 * 2.87e-04);
  smHiggs->SetPoint(11, 210., 4.723 * 2.34e-04);
  smHiggs->SetPoint(12, 220., 4.288 * 1.96e-04);
  smHiggs->SetPoint(13, 230., 3.908 * 1.68e-04);
  smHiggs->SetPoint(14, 240., 3.581 * 1.45e-04);
  smHiggs->SetPoint(15, 250., 3.312 * 1.27e-04);
  smHiggs->SetPoint(16, 260., 3.072 * 1.12e-04);
  smHiggs->SetPoint(17, 270., 2.864 * 1.00e-04);
  smHiggs->SetPoint(18, 280., 2.696 * 8.98e-05);
  smHiggs->SetPoint(19, 290., 2.546 * 8.09e-05);
  smHiggs->SetPoint(20, 300., 2.418 * 7.33e-05);
  smHiggs->SetPoint(21, 320., 2.248 * 6.12e-05);
  smHiggs->SetPoint(22, 340., 2.199 * 5.20e-05);
  smHiggs->SetPoint(23, 360., 2.359 * 4.23e-05);
  smHiggs->SetPoint(24, 380., 2.263 * 3.40e-05);
  smHiggs->SetPoint(25, 400., 2.035 * 2.84e-05);
  smHiggs->SetPoint(26, 450., 1.356 * 1.99e-05);
  smHiggs->SetPoint(27, 500., 0.8497* 1.53e-05);

  /// SM Higgs boson
  TGraph* vbfHiggs = new TGraph();
  vbfHiggs->SetPoint( 0, 100., 1.546 * 8.36e-02);
  vbfHiggs->SetPoint( 1, 110., 1.398 * 8.03e-02);
  vbfHiggs->SetPoint( 2, 120., 1.269 * 7.11e-02);
  vbfHiggs->SetPoint( 3, 130., 1.154 * 5.49e-02);
  vbfHiggs->SetPoint( 4, 140., 1.052 * 3.54e-02);
  vbfHiggs->SetPoint( 5, 150., .9617 * 1.79e-02);
  vbfHiggs->SetPoint( 6, 160., .8787 * 3.97e-03);
  vbfHiggs->SetPoint( 7, 170., .8173 * 9.20e-04);
  vbfHiggs->SetPoint( 8, 180., .7480 * 5.87e-04);
  vbfHiggs->SetPoint( 9, 190., .6925 * 3.76e-04);
  vbfHiggs->SetPoint(10, 200., .6371 * 2.87e-04);
  vbfHiggs->SetPoint(11, 210., .5869 * 2.34e-04);
  vbfHiggs->SetPoint(12, 220., .5420 * 1.96e-04);
  vbfHiggs->SetPoint(13, 230., .5011 * 1.68e-04);
  vbfHiggs->SetPoint(14, 240., .4641 * 1.45e-04);
  vbfHiggs->SetPoint(15, 250., .4304 * 1.27e-04);
  vbfHiggs->SetPoint(16, 260., .3988 * 1.12e-04);
  vbfHiggs->SetPoint(17, 270., .3715 * 1.00e-04);
  vbfHiggs->SetPoint(18, 280., .3461 * 8.98e-05);
  vbfHiggs->SetPoint(19, 290., .3226 * 8.09e-05);
  vbfHiggs->SetPoint(20, 300., .3011 * 7.33e-05);
  vbfHiggs->SetPoint(21, 320., .2627 * 6.12e-05);
  vbfHiggs->SetPoint(22, 340., .2286 * 5.20e-05);
  vbfHiggs->SetPoint(23, 360., .2018 * 4.23e-05);
  vbfHiggs->SetPoint(24, 380., .1808 * 3.40e-05);
  vbfHiggs->SetPoint(25, 400., .1620 * 2.84e-05);
  vbfHiggs->SetPoint(26, 440., .1304 * 1.88E-05);
  vbfHiggs->SetPoint(27, 500., .0949 * 1.53e-05);

  /// gg->A, tanB=30
  TGraph* ggfA30 = new TGraph();
  for(unsigned int i=0; i<41; ++i){
    float tanB = 30.; 
    float mass = 100+10*i;
    ggfA30->SetPoint(i, mass, mssm.Give_Xsec_ggFA(mass, tanB)*mssm.Give_BR_A_tautau(mass, tanB)/1000.);
  } 
  /// gg->A, tanB=25
  TGraph* ggfA25 = new TGraph();
  for(unsigned int i=0; i<41; ++i){
    float tanB = 25.; 
    float mass = 100+10*i;
    ggfA25->SetPoint(i, mass, mssm.Give_Xsec_ggFA(mass, tanB)*mssm.Give_BR_A_tautau(mass, tanB)/1000.);
  } 
  /// gg->h, tanB=20 
  TGraph* ggfA20 = new TGraph();
  for(unsigned int i=0; i<41; ++i){
    float tanB = 20.; 
    float mass = 100+10*i;
    ggfA20->SetPoint(i, mass, mssm.Give_Xsec_ggFA(mass, tanB)*mssm.Give_BR_A_tautau(mass, tanB)/1000.);
    std::cout << "ggfA (\tan\beta=20) " << std::endl;
    std::cout << "mass[GeV]: " << mass  << std::endl;
    std::cout << "xsec[pb ]: " << mssm.Give_Xsec_ggFA(mass, tanB)/1000. << std::endl;
    std::cout << "BR       : " << mssm.Give_BR_A_tautau(mass, tanB) << std::endl;
  } 


  /// bb->A, tanB=30
  TGraph* bbA5f30 = new TGraph();
  for(unsigned int i=0; i<41; ++i){
    float tanB = 30.; 
    float mass = 100+10*i;
    bbA5f30->SetPoint(i, mass, mssm.Give_Xsec_bbA5f(mass, tanB)*mssm.Give_BR_A_tautau(mass, tanB)/1000.);
  }  
  /// bb->A, tanB=25
  TGraph* bbA5f25 = new TGraph();
  for(unsigned int i=0; i<41; ++i){
    float tanB = 25.; 
    float mass = 100+10*i;
    bbA5f25->SetPoint(i, mass, mssm.Give_Xsec_bbA5f(mass, tanB)*mssm.Give_BR_A_tautau(mass, tanB)/1000.);
  }  
  /// bb->A, tanB=20
  TGraph* bbA5f20 = new TGraph();
  for(unsigned int i=0; i<41; ++i){
    float tanB = 20.; 
    float mass = 100+10*i;
    bbA5f20->SetPoint(i, mass, mssm.Give_Xsec_bbA5f(mass, tanB)*mssm.Give_BR_A_tautau(mass, tanB)/1000.);
    std::cout << "bbA$f (\tan\beta=20) " << std::endl;
    std::cout << "mass[GeV]: " << mass   << std::endl;
    std::cout << "xsec4[pb]: " << mssm.Give_Xsec_bbA4f(mass, tanB)/1000. << std::endl;
    std::cout << "xsec5[pb]: " << mssm.Give_Xsec_bbA5f(mass, tanB)/1000. << std::endl;
    std::cout << "xsecS[pb]: " << mssm.GiveXsec_Santander_A(mass, tanB)/1000. << std::endl;
    std::cout << "BR       : " << mssm.Give_BR_A_tautau(mass, tanB) << std::endl;
  }  

  TCanvas* canv1 = new TCanvas("canv1", "Cross Section * BR", 600, 600);
  canv1  ->cd();
  canv1  ->SetLogy(1);

  ggfA30 ->SetMaximum(10e3);
  ggfA30 ->SetMinimum(10e-7);
  ggfA30 ->SetLineColor(kRed-2 );
  ggfA30 ->SetLineWidth(  3.   );
  ggfA30 ->Draw("AC");
  ggfA30 ->GetXaxis()->SetTitle("m_{A/Higgs} [GeV]");
  ggfA30 ->GetXaxis()->SetLabelFont(62);
  ggfA30 ->GetXaxis()->SetTitleFont(62);
  ggfA30 ->GetXaxis()->SetTitleColor(1);
  ggfA30 ->GetXaxis()->SetTitleOffset(1.05);

  ggfA30 ->GetYaxis()->SetTitle("#sigma BR(Higgs/#phi#rightarrow#tau#tau) [pb]");
  ggfA30 ->GetYaxis()->SetLabelFont(62);
  ggfA30 ->GetYaxis()->SetTitleOffset(1.05);
  ggfA30 ->GetYaxis()->SetLabelSize(0.03);

  /*
  ggfA25 ->SetLineColor(kRed+0 );
  ggfA25 ->SetLineWidth(  3.   );
  ggfA25 ->Draw("Csame");
  */
  ggfA20 ->SetLineColor(kRed+2 );
  ggfA20 ->SetLineWidth(  3.   );
  ggfA20 ->Draw("Csame");

  bbA5f30->SetLineColor(kBlue-2);
  bbA5f30->SetLineWidth(    3. );
  bbA5f30->Draw("Csame");
  /*
  bbA5f25->SetLineColor(kBlue+0);
  bbA5f25->SetLineWidth(    3. );
  bbA5f25->Draw("Csame");
  */
  bbA5f20->SetLineColor(kBlue+2);
  bbA5f20->SetLineWidth(    3. );
  bbA5f20->Draw("Csame");

  smHiggs->SetLineColor(kBlack );
  smHiggs->SetLineWidth(    3. );
  smHiggs->Draw("Csame");

  vbfHiggs->SetLineColor(kGray+2);
  vbfHiggs->SetLineWidth(    3. );
  vbfHiggs->SetLineStyle(    2. );
  vbfHiggs->Draw("Csame");

  TPaveText* smH  = new TPaveText(0.24, 0.14, 0.40, 0.25, "NDC");
  smH->AddText("VBF Higgs(SM)");
  smH->SetBorderSize(    0    );
  smH->SetFillStyle (    0    );
  smH->SetTextAlign (   12    );
  smH->SetTextSize  ( 0.04    );
  smH->SetTextColor ( kGray+2 );
  smH->SetTextFont  (   62    );
  smH->Draw("same");

  TPaveText* smH  = new TPaveText(0.52, 0.68, 0.25, 0.35, "NDC");
  smH->AddText("gg#rightarrowHiggs(SM)");
  smH->SetBorderSize(    0    );
  smH->SetFillStyle (    0    );
  smH->SetTextAlign (   12    );
  smH->SetTextSize  ( 0.04    );
  smH->SetTextColor ( kBlack  );
  smH->SetTextFont  (   62    );
  smH->Draw("same");

  TPaveText* bbH  = new TPaveText(0.40, 0.53, 0.57, 0.57, "NDC");
  bbH->AddText("bb#rightarrow#phi");
  bbH->SetBorderSize(    0    );
  bbH->SetFillStyle (    0    );
  bbH->SetTextAlign (   12    );
  bbH->SetTextSize  ( 0.04    );
  bbH->SetTextColor ( kBlue+0 );
  bbH->SetTextFont  (   62    );
  bbH->Draw("same");

  TPaveText* ggF  = new TPaveText(0.33, 0.43, 0.40, 0.47, "NDC");
  ggF->AddText("gg#rightarrow#phi");
  ggF->SetBorderSize(    0    );
  ggF->SetFillStyle (    0    );
  ggF->SetTextAlign (   12    );
  ggF->SetTextSize  ( 0.04    );
  ggF->SetTextColor ( kRed+0  );
  ggF->SetTextFont  (   62    );
  ggF->Draw("same");

  TLegend* leg = new TLegend(0.51,0.60,0.88,0.88);
  leg->SetBorderSize( 0 );
  leg->SetFillStyle ( 0 );
  leg->AddEntry( bbA5f20 , "bb#rightarrow#phi (tan #beta=30)",  "L" );
  //leg->AddEntry( bbA5f25, "bb#rightarrow#phi (tan #beta=25)",  "L" );
  leg->AddEntry( bbA5f30 , "bb#rightarrow#phi (tan #beta=20)",  "L" );
  leg->AddEntry( ggfA30  , "gg#rightarrow#phi (tan #beta=30)",  "L" );
  //leg->AddEntry( ggfA25 , "gg#rightarrow#phi (tan #beta=25)",  "L" );
  leg->AddEntry( ggfA20  , "gg#rightarrow#phi (tan #beta=20)",  "L" );
  leg->AddEntry( smHiggs , "gg#rightarrow Higgs(SM)",           "L" );
  leg->AddEntry( vbfHiggs, "VBF Higgs(SM)",                     "L" );
  leg->Draw("same");

  return 0; 
}
