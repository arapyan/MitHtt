/**
   \class   plot_acc plot_acc.C "MitHtt/macros/plot_acc.c"

   \brief   macro to create the acceptance plots for the Inclusive Selection as descriobed in AN-11-153

   macro to creat the acceptance plots as function of mA/mH after the Inclusive Selection as desribed 
   in AN-11-153 (v20). The acceptances have been gathered from the analyses. Run this macro in root in
   the following way:

   root -l
   .x plot_acc.C

   A plot will be created that can be safed manually. 
*/

int
plot_acc()
{
  /// ggH (emu)
  TGraph* ggHemu   = new TGraph();
  ggHemu  ->SetPoint( 0,  90., 0.00492);
  ggHemu  ->SetPoint( 1, 100., 0.00590);
  ggHemu  ->SetPoint( 2, 120., 0.00855);
  ggHemu  ->SetPoint( 3, 140., 0.01069);
  ggHemu  ->SetPoint( 4, 160., 0.01198);
  ggHemu  ->SetPoint( 5, 180., 0.01371);
  ggHemu  ->SetPoint( 6, 300., 0.02118);
  ggHemu  ->SetPoint( 7, 400., 0.02259);
  ggHemu  ->SetPoint( 8, 450., 0.02289);
  ggHemu  ->SetPoint( 9, 500., 0.02191);

  /// ggH (etau)
  TGraph* ggHetau = new TGraph();
  ggHetau ->SetPoint( 0,  90., 0.008369);
  ggHetau ->SetPoint( 1, 100., 0.011047);
  ggHetau ->SetPoint( 2, 120., 0.015310);
  ggHetau ->SetPoint( 3, 140., 0.019884);
  ggHetau ->SetPoint( 4, 160., 0.023484);
  ggHetau ->SetPoint( 5, 180., 0.027789);
  ggHetau ->SetPoint( 6, 300., 0.037348);
  ggHetau ->SetPoint( 7, 350., 0.042010);
  ggHetau ->SetPoint( 8, 400., 0.043779);
  ggHetau ->SetPoint( 9, 500., 0.040514);

  /// ggH (mutau)
  TGraph* ggHmutau = new TGraph();
  ggHmutau->SetPoint( 0,  90., 0.013905);
  ggHmutau->SetPoint( 1, 100., 0.017807);
  ggHmutau->SetPoint( 2, 120., 0.024592);
  ggHmutau->SetPoint( 3, 140., 0.030957);
  ggHmutau->SetPoint( 4, 160., 0.037478);
  ggHmutau->SetPoint( 5, 180., 0.040827);
  ggHmutau->SetPoint( 6, 300., 0.056171);
  ggHmutau->SetPoint( 7, 350., 0.058985);
  ggHmutau->SetPoint( 8, 400., 0.057913);
  ggHmutau->SetPoint( 9, 500., 0.053962);

  /// bbH (emu)
  TGraph* bbHemu   = new TGraph();
  bbHemu  ->SetPoint( 0,  90., 0.00487);
  bbHemu  ->SetPoint( 1, 100., 0.00588);
  bbHemu  ->SetPoint( 2, 120., 0.00846);
  bbHemu  ->SetPoint( 3, 140., 0.01087);
  bbHemu  ->SetPoint( 4, 160., 0.01284);
  bbHemu  ->SetPoint( 5, 180., 0.01584);
  bbHemu  ->SetPoint( 6, 300., 0.02341);
  bbHemu  ->SetPoint( 7, 400., 0.02612);
  bbHemu  ->SetPoint( 8, 450., 0.02720);
  bbHemu  ->SetPoint( 9, 500., 0.02890);

  /// bbH (etau)
  TGraph* bbHetau = new TGraph();
  bbHetau ->SetPoint( 0,  90., 0.009058);
  bbHetau ->SetPoint( 1, 100., 0.011240);
  bbHetau ->SetPoint( 2, 120., 0.016419);
  bbHetau ->SetPoint( 3, 140., 0.022724);
  bbHetau ->SetPoint( 4, 160., 0.024868);
  bbHetau ->SetPoint( 5, 180., 0.030113);
  bbHetau ->SetPoint( 6, 300., 0.044448);
  bbHetau ->SetPoint( 7, 350., 0.046025);
  bbHetau ->SetPoint( 8, 400., 0.048670);
  bbHetau ->SetPoint( 9, 500., 0.050621);

  /// bbH (mutau)
  TGraph* bbHmutau = new TGraph();
  bbHmutau->SetPoint( 0,  90., 0.014981);
  bbHmutau->SetPoint( 1, 100., 0.018533);
  bbHmutau->SetPoint( 2, 120., 0.026935);
  bbHmutau->SetPoint( 3, 140., 0.033277);
  bbHmutau->SetPoint( 4, 160., 0.038454);
  bbHmutau->SetPoint( 5, 180., 0.045204);
  bbHmutau->SetPoint( 6, 300., 0.060625);
  bbHmutau->SetPoint( 7, 350., 0.050619);
  bbHmutau->SetPoint( 8, 400., 0.067533);
  bbHmutau->SetPoint( 9, 500., 0.069952);


  TCanvas* canv1 = new TCanvas("canv1", "Full Acceptance", 600, 600);
  canv1  ->cd();
  canv1  ->SetLogy(1);
  canv1  ->SetGridx(1);
  canv1  ->SetGridy(1);

  ggHemu ->SetMaximum(1);
  //ggHemu ->SetMinimum(10e-3);
  ggHemu  ->SetLineColor(kRed);
  ggHemu  ->SetLineWidth(  3.   );
  ggHemu  ->Draw("AC");
  ggHemu  ->GetXaxis()->SetTitle("m_{A} [GeV]");
  ggHemu  ->GetXaxis()->SetLabelFont(62);
  ggHemu  ->GetXaxis()->SetTitleFont(62);
  ggHemu  ->GetXaxis()->SetTitleColor(1);
  ggHemu  ->GetXaxis()->SetTitleOffset(1.05);

  ggHemu  ->GetYaxis()->SetTitle("A #times #epsilon");
  ggHemu  ->GetYaxis()->SetLabelFont(62);
  ggHemu  ->GetYaxis()->SetTitleOffset(1.05);
  ggHemu  ->GetYaxis()->SetLabelSize(0.03);

  bbHemu  ->SetLineColor(kRed   );
  bbHemu  ->SetLineWidth(  3.   );
  bbHemu  ->SetLineStyle(  2.   );
  bbHemu  ->Draw("Csame");

  ggHetau ->SetLineColor(kBlue  );
  ggHetau ->SetLineWidth(  3.   );
  ggHetau ->Draw("Csame");

  bbHetau ->SetLineColor(kBlue  );
  bbHetau ->SetLineWidth(  3.   );
  bbHetau ->SetLineStyle(  2.   );
  bbHetau ->Draw("Csame");

  ggHmutau->SetLineColor(kBlack );
  ggHmutau->SetLineWidth(  3.   );
  ggHmutau->Draw("Csame");

  bbHmutau->SetLineColor(kBlack );
  bbHmutau->SetLineWidth(  3.   );
  bbHmutau->SetLineStyle(  2.   );
  bbHmutau->Draw("Csame");

  TPaveText* sel  = new TPaveText(0.50, 0.14, 0.75, 0.25, "NDC");
  sel->AddText("Inclusive Selection");
  sel->SetBorderSize(    0    );
  sel->SetFillStyle (    0    );
  sel->SetTextAlign (   12    );
  sel->SetTextSize  ( 0.04    );
  sel->SetTextColor ( kBlack  );
  sel->SetTextFont  (   62    );
  sel->Draw("same");

  TLegend* leg0 = new TLegend(0.50,0.70,1.00,0.90);
  leg0->SetBorderSize( 0 );
  leg0->SetFillStyle ( 0 );
  leg0->SetHeader( "gg#rightarrow#phi (MSSM)" );
  leg0->AddEntry( ggHemu  , "e#mu-channel"   ,  "L" );
  leg0->AddEntry( ggHetau , "e#tau-channel"  ,  "L" );
  leg0->AddEntry( ggHmutau, "#mu#tau-channel",  "L" );
  leg0->Draw("same");

  TLegend* leg1 = new TLegend(0.10,0.70,0.60,0.90);
  leg1->SetBorderSize( 0 );
  leg1->SetFillStyle ( 0 );
  leg1->SetHeader( "bb#rightarrow#phi (MSSM)" );
  leg1->AddEntry( bbHemu  , "e#mu-channel"   ,  "L" );
  leg1->AddEntry( bbHetau , "e#tau-channel"  ,  "L" );
  leg1->AddEntry( bbHmutau, "#mu#tau-channel",  "L" );
  leg1->Draw("same");

  return 0; 
}
