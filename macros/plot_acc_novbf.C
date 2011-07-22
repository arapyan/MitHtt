int
plot_acc_novbf()
{
  /// ggH (emu)
  TGraph* ggHemu   = new TGraph();
  ggHemu  ->SetPoint( 0, 110., 0.007216);
  ggHemu  ->SetPoint( 1, 115., 0.007547);
  ggHemu  ->SetPoint( 2, 120., 0.008011);
  ggHemu  ->SetPoint( 3, 125., 0.008591);
  ggHemu  ->SetPoint( 4, 130., 0.009488);
  ggHemu  ->SetPoint( 5, 135., 0.009890);
  ggHemu  ->SetPoint( 6, 145., 0.010621);

  /// ggH (etau)
  TGraph* ggHetau = new TGraph();
  ggHetau ->SetPoint( 0, 110., 0.012894);
  ggHetau ->SetPoint( 1, 115., 0.014369);
  ggHetau ->SetPoint( 2, 120., 0.015495);
  ggHetau ->SetPoint( 3, 125., 0.015173);
  ggHetau ->SetPoint( 4, 130., 0.016117);
  ggHetau ->SetPoint( 5, 135., 0.019066);
  ggHetau ->SetPoint( 6, 145., 0.019967);

  /// ggH (mutau)
  TGraph* ggHmutau = new TGraph();
  ggHmutau->SetPoint( 0, 110., 0.019957);
  ggHmutau->SetPoint( 1, 115., 0.023482);
  ggHmutau->SetPoint( 2, 120., 0.023996);
  ggHmutau->SetPoint( 3, 125., 0.025581);
  ggHmutau->SetPoint( 4, 130., 0.027448);
  ggHmutau->SetPoint( 5, 135., 0.028570);
  ggHmutau->SetPoint( 6, 145., 0.031137);

  /// qqH (emu)
  TGraph* qqHemu   = new TGraph();
  qqHemu  ->SetPoint( 0, 110., 0.007128);
  qqHemu  ->SetPoint( 1, 115., 0.007588);
  qqHemu  ->SetPoint( 2, 120., 0.008373);
  qqHemu  ->SetPoint( 3, 125., 0.008454);
  qqHemu  ->SetPoint( 4, 130., 0.009126);
  qqHemu  ->SetPoint( 5, 135., 0.009457);
  qqHemu  ->SetPoint( 6, 145., 0.010225);

  /// qqH (etau)
  TGraph* qqHetau = new TGraph();
  qqHetau ->SetPoint( 0, 110., 0.012506);
  qqHetau ->SetPoint( 1, 115., 0.013459);
  qqHetau ->SetPoint( 2, 120., 0.013678);
  qqHetau ->SetPoint( 3, 125., 0.015173);
  qqHetau ->SetPoint( 4, 130., 0.016117);
  qqHetau ->SetPoint( 5, 135., 0.017092);
  qqHetau ->SetPoint( 6, 145., 0.019217);

  /// qqH (mutau)
  TGraph* qqHmutau = new TGraph();
  qqHmutau->SetPoint( 0, 110., 0.019134);
  qqHmutau->SetPoint( 1, 115., 0.021443);
  qqHmutau->SetPoint( 2, 120., 0.022787);
  qqHmutau->SetPoint( 3, 125., 0.023985);
  qqHmutau->SetPoint( 4, 130., 0.025318);
  qqHmutau->SetPoint( 5, 135., 0.026565);
  qqHmutau->SetPoint( 6, 145., 0.027886);


  TCanvas* canv1 = new TCanvas("canv1", "Full Acceptance", 600, 600);
  canv1  ->cd()
  canv1  ->SetLogy(1);
  canv1  ->SetGridx(1);
  canv1  ->SetGridy(1);

  ggHemu ->SetMaximum(1);
  ggHemu ->SetMinimum(10e-4);
  ggHemu  ->SetLineColor(kRed);
  ggHemu  ->SetLineWidth(  3.   );
  ggHemu  ->Draw("AC");
  ggHemu  ->GetXaxis()->SetTitle("m_{H} [GeV]");
  ggHemu  ->GetXaxis()->SetLabelFont(62);
  ggHemu  ->GetXaxis()->SetTitleFont(62);
  ggHemu  ->GetXaxis()->SetTitleColor(1);
  ggHemu  ->GetXaxis()->SetTitleOffset(1.05);

  ggHemu  ->GetYaxis()->SetTitle("A #times #epsilon");
  ggHemu  ->GetYaxis()->SetLabelFont(62);
  ggHemu  ->GetYaxis()->SetTitleOffset(1.05);
  ggHemu  ->GetYaxis()->SetLabelSize(0.03);

  qqHemu  ->SetLineColor(kRed   );
  qqHemu  ->SetLineWidth(  3.   );
  qqHemu  ->SetLineStyle(  2.   );
  qqHemu  ->Draw("Csame");

  ggHetau ->SetLineColor(kBlue  );
  ggHetau ->SetLineWidth(  3.   );
  ggHetau ->Draw("Csame");

  qqHetau ->SetLineColor(kBlue  );
  qqHetau ->SetLineWidth(  3.   );
  qqHetau ->SetLineStyle(  2.   );
  qqHetau ->Draw("Csame");

  ggHmutau->SetLineColor(kBlack );
  ggHmutau->SetLineWidth(  3.   );
  ggHmutau->Draw("Csame");

  qqHmutau->SetLineColor(kBlack );
  qqHmutau->SetLineWidth(  3.   );
  qqHmutau->SetLineStyle(  2.   );
  qqHmutau->Draw("Csame");

  TPaveText* sel  = new TPaveText(0.50, 0.14, 0.75, 0.25, "NDC");
  sel->AddText("Non-VBF Category");
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
  leg0->SetHeader( "gg#rightarrow H (SM)" );
  leg0->AddEntry( ggHemu  , "e#mu-channel"   ,  "L" );
  leg0->AddEntry( ggHetau , "e#tau-channel"  ,  "L" );
  leg0->AddEntry( ggHmutau, "#mu#tau-channel",  "L" );
  leg0->Draw("same");

  TLegend* leg1 = new TLegend(0.10,0.70,0.60,0.90);
  leg1->SetBorderSize( 0 );
  leg1->SetFillStyle ( 0 );
  leg1->SetHeader( "qq#rightarrow H (SM)" );
  leg1->AddEntry( qqHemu  , "e#mu-channel"   ,  "L" );
  leg1->AddEntry( qqHetau , "e#tau-channel"  ,  "L" );
  leg1->AddEntry( qqHmutau, "#mu#tau-channel",  "L" );
  leg1->Draw("same");

  return 0; 
}
