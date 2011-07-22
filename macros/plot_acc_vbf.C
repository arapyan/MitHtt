int
plot_acc_vbf()
{
  /// ggH (emu)
  TGraph* ggHemu   = new TGraph();
  ggHemu  ->SetPoint( 0, 110., 0.000068);
  ggHemu  ->SetPoint( 1, 115., 0.000072);
  ggHemu  ->SetPoint( 2, 120., 0.000010);
  ggHemu  ->SetPoint( 3, 125., 0.000017);
  ggHemu  ->SetPoint( 4, 130., 0.000049);
  ggHemu  ->SetPoint( 5, 135., 0.000034);
  ggHemu  ->SetPoint( 6, 145., 0.000073);
  ggHemu->Fit("pol1", "R", "", 110, 145);

  /// ggH (etau)
  TGraph* ggHetau = new TGraph();
  ggHetau ->SetPoint( 0, 110., 0.000062);
  ggHetau ->SetPoint( 1, 115., 0.000045);
  ggHetau ->SetPoint( 2, 120., 0.000061);
  ggHetau ->SetPoint( 3, 125., 0.000050);
  ggHetau ->SetPoint( 4, 130., 0.000079);
  ggHetau ->SetPoint( 5, 135., 0.000058);
  ggHetau ->SetPoint( 6, 145., 0.000085);
  ggHetau->Fit("pol1", "R", "",110, 145);

  /// ggH (mutau)
  TGraph* ggHmutau = new TGraph();
  ggHmutau->SetPoint( 0, 110., 0.000074);
  ggHmutau->SetPoint( 1, 115., 0.000039);
  ggHmutau->SetPoint( 2, 120., 0.000109);
  ggHmutau->SetPoint( 3, 125., 0.000128);
  ggHmutau->SetPoint( 4, 130., 0.000115);
  ggHmutau->SetPoint( 5, 135., 0.000180);
  ggHmutau->SetPoint( 6, 145., 0.000134);
  ggHmutau->Fit("pol1", "R", "",110,145);

  /// qqH (emu)
  TGraph* qqHemu   = new TGraph();
  qqHemu  ->SetPoint( 0, 110., 0.001305);
  qqHemu  ->SetPoint( 1, 115., 0.001201);
  qqHemu  ->SetPoint( 2, 120., 0.001337);
  qqHemu  ->SetPoint( 3, 125., 0.001320);
  qqHemu  ->SetPoint( 4, 130., 0.001477);
  qqHemu  ->SetPoint( 5, 135., 0.001603);
  qqHemu  ->SetPoint( 6, 145., 0.001743);

  /// qqH (etau)
  TGraph* qqHetau = new TGraph();
  qqHetau ->SetPoint( 0, 110., 0.001543);
  qqHetau ->SetPoint( 1, 115., 0.001718);
  qqHetau ->SetPoint( 2, 120., 0.001704);
  qqHetau ->SetPoint( 3, 125., 0.001878);
  qqHetau ->SetPoint( 4, 130., 0.002129);
  qqHetau ->SetPoint( 5, 135., 0.001992);
  qqHetau ->SetPoint( 6, 145., 0.002279);

  /// qqH (mutau)
  TGraph* qqHmutau = new TGraph();
  qqHmutau->SetPoint( 0, 110., 0.002394);
  qqHmutau->SetPoint( 1, 115., 0.002462);
  qqHmutau->SetPoint( 2, 120., 0.002692);
  qqHmutau->SetPoint( 3, 125., 0.002520);
  qqHmutau->SetPoint( 4, 130., 0.002706);
  qqHmutau->SetPoint( 5, 135., 0.003144);
  qqHmutau->SetPoint( 6, 145., 0.003423);


  TCanvas* canv1 = new TCanvas("canv1", "Full Acceptance", 600, 600);
  canv1  ->cd();
  canv1  ->SetLogy(1);
  canv1  ->SetGridx(1);
  canv1  ->SetGridy(1);

  qqHemu ->SetMaximum(1);
  qqHemu ->SetMinimum(10e-7);
  qqHemu  ->SetLineColor(kRed);
  qqHemu  ->SetLineWidth(  3.   );
  qqHemu  ->SetLineStyle(  2.   );
  qqHemu  ->Draw("AC");
  qqHemu  ->GetXaxis()->SetTitle("m_{H} [GeV]");
  qqHemu  ->GetXaxis()->SetLabelFont(62);
  qqHemu  ->GetXaxis()->SetTitleFont(62);
  qqHemu  ->GetXaxis()->SetTitleColor(1);
  qqHemu  ->GetXaxis()->SetTitleOffset(1.05);

  qqHemu  ->GetYaxis()->SetTitle("A #times #epsilon");
  qqHemu  ->GetYaxis()->SetLabelFont(62);
  qqHemu  ->GetYaxis()->SetTitleOffset(1.05);
  qqHemu  ->GetYaxis()->SetLabelSize(0.03);

  ggHemu  ->SetLineColor(kRed   );
  ggHemu  ->SetLineWidth(  3.   );
  ggHemu  ->GetFunction("pol1")->SetLineColor(kRed);
  ggHemu  ->GetFunction("pol1")->SetLineWidth(  3. );
  ggHemu  ->GetFunction("pol1")->Draw("same");
  //ggHemu  ->Draw("Csame");

  ggHetau ->SetLineColor(kBlue  );
  ggHetau ->SetLineWidth(  3.   );
  ggHetau ->GetFunction("pol1")->SetLineColor(kBlue);
  ggHetau ->GetFunction("pol1")->SetLineWidth(  3. );
  ggHetau ->GetFunction("pol1")->Draw("same");
  //ggHetau ->Draw("Csame");

  qqHetau ->SetLineColor(kBlue  );
  qqHetau ->SetLineWidth(  3.   );
  qqHetau ->SetLineStyle(  2.   );
  qqHetau ->Draw("Csame");

  ggHmutau->SetLineColor(kBlack );
  ggHmutau->SetLineWidth(  3.   );
  ggHmutau->GetFunction("pol1")->SetLineColor(kBlack);
  ggHmutau->GetFunction("pol1")->SetLineWidth(  3. );
  ggHmutau->GetFunction("pol1")->Draw("same");
  //ggHmutau->Draw("Csame");

  qqHmutau->SetLineColor(kBlack );
  qqHmutau->SetLineWidth(  3.   );
  qqHmutau->SetLineStyle(  2.   );
  qqHmutau->Draw("Csame");

  TPaveText* sel  = new TPaveText(0.50, 0.14, 0.75, 0.25, "NDC");
  sel->AddText("VBF Category");
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
