int
plot_acc_nob()
{
  /// ggH (emu)
  TGraph* ggHemu   = new TGraph();
  ggHemu  ->SetPoint( 0,  90., 0.00429);
  ggHemu  ->SetPoint( 1, 100., 0.00531);
  ggHemu  ->SetPoint( 2, 120., 0.00769);
  ggHemu  ->SetPoint( 3, 140., 0.00937);
  ggHemu  ->SetPoint( 4, 160., 0.01063);
  ggHemu  ->SetPoint( 5, 180., 0.01371);
  ggHemu  ->SetPoint( 6, 300., 0.01751);
  ggHemu  ->SetPoint( 7, 400., 0.01813);
  ggHemu  ->SetPoint( 8, 450., 0.01815);
  ggHemu  ->SetPoint( 9, 500., 0.01709);

  /// ggH (etau)
  TGraph* ggHetau = new TGraph();
  ggHetau ->SetPoint( 0,  90., 0.008369);
  ggHetau ->SetPoint( 1, 100., 0.010946);
  ggHetau ->SetPoint( 2, 120., 0.015115);
  ggHetau ->SetPoint( 3, 140., 0.019884);
  ggHetau ->SetPoint( 4, 160., 0.022844);
  ggHetau ->SetPoint( 5, 180., 0.026670);
  ggHetau ->SetPoint( 6, 300., 0.037348);
  ggHetau ->SetPoint( 7, 400., 0.039514);
  ggHetau ->SetPoint( 8, 450., 0.038138);
  ggHetau ->SetPoint( 9, 500., 0.035500);

  /// ggH (mutau)
  TGraph* ggHmutau = new TGraph();
  ggHmutau->SetPoint( 0,  90., 0.013681);
  ggHmutau->SetPoint( 1, 100., 0.017738);
  ggHmutau->SetPoint( 2, 120., 0.024206);
  ggHmutau->SetPoint( 3, 140., 0.030261);
  ggHmutau->SetPoint( 4, 160., 0.036494);
  ggHmutau->SetPoint( 5, 180., 0.039271);
  ggHmutau->SetPoint( 6, 300., 0.051243);
  ggHmutau->SetPoint( 7, 400., 0.051832);
  ggHmutau->SetPoint( 8, 450., 0.052523);
  ggHmutau->SetPoint( 9, 500., 0.047083);


  /// bbH (emu)
  TGraph* bbHemu   = new TGraph();
  bbHemu  ->SetPoint( 0,  90., 0.00504);
  bbHemu  ->SetPoint( 1, 100., 0.00599);
  bbHemu  ->SetPoint( 2, 120., 0.00859);
  bbHemu  ->SetPoint( 3, 140., 0.01088);
  bbHemu  ->SetPoint( 4, 160., 0.01271);
  bbHemu  ->SetPoint( 5, 180., 0.01584);
  bbHemu  ->SetPoint( 6, 300., 0.02163);
  bbHemu  ->SetPoint( 7, 400., 0.02307);
  bbHemu  ->SetPoint( 8, 450., 0.02433);
  bbHemu  ->SetPoint( 9, 500., 0.02532);

  /// bbH (etau)
  TGraph* bbHetau = new TGraph();
  bbHetau ->SetPoint( 0,  90., 0.009058);
  bbHetau ->SetPoint( 1, 100., 0.009902);
  bbHetau ->SetPoint( 2, 120., 0.014198);
  bbHetau ->SetPoint( 3, 140., 0.019718);
  bbHetau ->SetPoint( 4, 160., 0.021350);
  bbHetau ->SetPoint( 5, 180., 0.025392);
  bbHetau ->SetPoint( 6, 300., 0.035614);
  bbHetau ->SetPoint( 7, 400., 0.038652);
  bbHetau ->SetPoint( 8, 450., 0.038943);
  bbHetau ->SetPoint( 9, 500., 0.038945);

  /// bbH (mutau)
  TGraph* bbHmutau = new TGraph();
  bbHmutau->SetPoint( 0,  90., 0.013663);
  bbHmutau->SetPoint( 1, 100., 0.016489);
  bbHmutau->SetPoint( 2, 120., 0.023989);
  bbHmutau->SetPoint( 3, 140., 0.028861);
  bbHmutau->SetPoint( 4, 160., 0.033172);
  bbHmutau->SetPoint( 5, 180., 0.038752);
  bbHmutau->SetPoint( 6, 300., 0.048483);
  bbHmutau->SetPoint( 7, 400., 0.053021);
  bbHmutau->SetPoint( 8, 450., 0.052450);
  bbHmutau->SetPoint( 9, 500., 0.054074);


  TCanvas* canv1 = new TCanvas("canv1", "Full Acceptance", 600, 600);
  canv1  ->cd()
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
  sel->AddText("Non-BTag Category");
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
