int
plot_acc_b()
{
  /// ggH (emu)
  TGraph* ggHemu   = new TGraph();
  ggHemu  ->SetPoint( 0,  90., 0.000019);
  ggHemu  ->SetPoint( 1, 100., 0.000031);
  ggHemu  ->SetPoint( 2, 120., 0.000105);
  ggHemu  ->SetPoint( 3, 140., 0.000166);
  ggHemu  ->SetPoint( 4, 160., 0.000117);
  ggHemu  ->SetPoint( 5, 180., 0.000228);
  ggHemu  ->SetPoint( 6, 300., 0.000306);
  ggHemu  ->SetPoint( 7, 400., 0.000402);
  ggHemu  ->SetPoint( 8, 450., 0.000226);
  ggHemu  ->SetPoint( 9, 500., 0.000336);
  ggHemu->Fit("pol3", "R", "",  90, 500);

  /// ggH (etau)
  TGraph* ggHetau = new TGraph();
  ggHetau ->SetPoint( 0,  90., 0.000232);
  ggHetau ->SetPoint( 1, 100., 0.000193);
  ggHetau ->SetPoint( 2, 120., 0.000275);
  ggHetau ->SetPoint( 3, 140., 0.000415);
  ggHetau ->SetPoint( 4, 160., 0.000454);
  ggHetau ->SetPoint( 5, 180., 0.000489);
  ggHetau ->SetPoint( 6, 300., 0.001119);
  ggHetau ->SetPoint( 7, 400., 0.001106);
  ggHetau ->SetPoint( 8, 450., 0.001163);
  ggHetau ->SetPoint( 9, 500., 0.001467);
  ggHetau->Fit("pol3", "R", "", 90, 500);

  /// ggH (mutau)
  TGraph* ggHmutau = new TGraph();
  ggHmutau->SetPoint( 0,  90., 0.000247);
  ggHmutau->SetPoint( 1, 100., 0.000307);
  ggHmutau->SetPoint( 2, 120., 0.000451);
  ggHmutau->SetPoint( 3, 140., 0.000743);
  ggHmutau->SetPoint( 4, 160., 0.000614);
  ggHmutau->SetPoint( 5, 180., 0.000872);
  ggHmutau->SetPoint( 6, 300., 0.001673);
  ggHmutau->SetPoint( 7, 400., 0.001591);
  ggHmutau->SetPoint( 8, 450., 0.001588);
  ggHmutau->SetPoint( 9, 500., 0.001771);
  ggHmutau->Fit("pol3", "R", "", 90,500);

  /// bbH (emu)
  TGraph* bbHemu   = new TGraph();
  bbHemu  ->SetPoint( 0,  90., 0.000619);
  bbHemu  ->SetPoint( 1, 100., 0.000766);
  bbHemu  ->SetPoint( 2, 120., 0.001148);
  bbHemu  ->SetPoint( 3, 140., 0.001551);
  bbHemu  ->SetPoint( 4, 160., 0.001838);
  bbHemu  ->SetPoint( 5, 180., 0.002518);
  bbHemu  ->SetPoint( 6, 300., 0.003179);
  bbHemu  ->SetPoint( 7, 400., 0.004527);
  bbHemu  ->SetPoint( 8, 450., 0.004158);
  bbHemu  ->SetPoint( 9, 500., 0.004468);

  /// bbH (etau)
  TGraph* bbHetau = new TGraph();
  bbHetau ->SetPoint( 0,  90., 0.002189);
  bbHetau ->SetPoint( 1, 100., 0.002741);
  bbHetau ->SetPoint( 2, 120., 0.004656);
  bbHetau ->SetPoint( 3, 140., 0.006371);
  bbHetau ->SetPoint( 4, 160., 0.007031);
  bbHetau ->SetPoint( 5, 180., 0.009215);
  bbHetau ->SetPoint( 6, 300., 0.013464);
  bbHetau ->SetPoint( 7, 400., 0.014848);
  bbHetau ->SetPoint( 8, 450., 0.016270);
  bbHetau ->SetPoint( 9, 500., 0.015556);

  /// bbH (mutau)
  TGraph* bbHmutau = new TGraph();
  bbHmutau->SetPoint( 0,  90., 0.003360);
  bbHmutau->SetPoint( 1, 100., 0.004575);
  bbHmutau->SetPoint( 2, 120., 0.006725);
  bbHmutau->SetPoint( 3, 140., 0.009464);
  bbHmutau->SetPoint( 4, 160., 0.010915);
  bbHmutau->SetPoint( 5, 180., 0.012411);
  bbHmutau->SetPoint( 6, 300., 0.019682);
  bbHmutau->SetPoint( 7, 400., 0.019740);
  bbHmutau->SetPoint( 8, 450., 0.020162);
  bbHmutau->SetPoint( 9, 500., 0.021271);


  TCanvas* canv1 = new TCanvas("canv1", "Full Acceptance", 600, 600);
  canv1  ->cd()
  canv1  ->SetLogy(1);
  canv1  ->SetGridx(1);
  canv1  ->SetGridy(1);

  bbHemu ->SetMaximum(1);
  bbHemu ->SetMinimum(10e-6);
  bbHemu  ->GetXaxis()->SetTitle("m_{A} [GeV]");
  bbHemu  ->GetXaxis()->SetLabelFont(62);
  bbHemu  ->GetXaxis()->SetTitleFont(62);
  bbHemu  ->GetXaxis()->SetTitleColor(1);
  bbHemu  ->GetXaxis()->SetTitleOffset(1.05);

  bbHemu  ->GetYaxis()->SetTitle("A #times #epsilon");
  bbHemu  ->GetYaxis()->SetLabelFont(62);
  bbHemu  ->GetYaxis()->SetTitleOffset(1.05);
  bbHemu  ->GetYaxis()->SetLabelSize(0.03);

  bbHemu  ->SetLineColor(kRed);
  bbHemu  ->SetLineWidth(  3.   );
  bbHemu  ->SetLineStyle(  2.   );
  bbHemu  ->Draw("AC");

  ggHetau ->SetLineColor(kBlue  );
  ggHetau ->SetLineWidth(  3.   );
  ggHetau ->GetFunction("pol3")->SetLineColor(kBlue);
  ggHetau ->GetFunction("pol3")->SetLineWidth(  3. );
  ggHetau ->GetFunction("pol3")->Draw("same");
  //ggHetau ->Draw("Csame");

  ggHmutau->SetLineColor(kBlack );
  ggHmutau->SetLineWidth(  3.   );
  ggHmutau ->GetFunction("pol3")->SetLineColor(kBlack);
  ggHmutau ->GetFunction("pol3")->SetLineWidth(   3. );
  ggHmutau ->GetFunction("pol3")->Draw("same");
  //ggHmutau->Draw("Csame");

  ggHemu  ->SetLineColor(kRed   );
  ggHemu  ->SetLineWidth(  3.   );
  //ggHemu  ->SetLineStyle(  2.   );
  ggHemu  ->GetFunction("pol3")->SetLineColor(kRed);
  ggHemu  ->GetFunction("pol3")->SetLineWidth(  3. );
  ggHemu  ->GetFunction("pol3")->Draw("same");
  //ggHemu  ->Draw("Csame");


  bbHetau ->SetLineColor(kBlue  );
  bbHetau ->SetLineWidth(  3.   );
  bbHetau ->SetLineStyle(  2.   );
  bbHetau ->Draw("Csame");


  bbHmutau->SetLineColor(kBlack );
  bbHmutau->SetLineWidth(  3.   );
  bbHmutau->SetLineStyle(  2.   );
  bbHmutau->Draw("Csame");

  TPaveText* sel  = new TPaveText(0.50, 0.14, 0.75, 0.25, "NDC");
  sel->AddText("BTag Category");
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
