#include "CPlot.hh"
#include <TLatex.h>

CPlot::CPlot()
{
  TString name = "plot"; 
  name += sCount;
  CPlot(name,"","",""); 
}

CPlot::CPlot(TString name, TString title, TString xtitle, TString ytitle):
fStack(0),
fName(name),
fTitle(title),
fXTitle(xtitle),
fYTitle(ytitle),
fXmin(0),
fXmax(0),
fYmin(0),
fYmax(0),
fLogx(0),
fLogy(0),
fGridx(0),
fGridy(0),
fRebin(1),
fLeg(0),
fShowStats(0),
fStatsX(0.68),
fStatsY(0.90)
{
  sCount++;
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddHist1D(TH1F *h, TString drawopt, int color, int linesty, int fillsty, int linecolor)
{
  if(!h)
    return;
   
  if(linecolor==0)
    h->SetLineColor(color);
  else
    h->SetLineColor(linecolor);
  h->SetLineStyle(linesty);
  h->SetFillColor(color);
  h->SetFillStyle(fillsty);
  
  if(drawopt.CompareTo("E",TString::kIgnoreCase)==0)
    h->SetMarkerSize(1.3);
  
  CPlotItem item;
  item.hist1D = h;
  item.drawopt = drawopt;
  fItems.push_back(item);
}

void CPlot::AddHist1D(TH1F *h, TString label, TString drawopt, int color, int linesty, int fillsty)
{
  if(!h)
    return;

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
 
  if(drawopt.CompareTo("E",TString::kIgnoreCase)==0) {
    //fLeg->AddEntry(h,label,"P");
    fLeg->AddEntry(h,label,"PL");
    h->SetMarkerSize(1.3);
    h->SetMarkerStyle(20);
    h->SetLineWidth(3);
  } else {
    if(fillsty>0) fLeg->AddEntry(h,label,"F");
    else          fLeg->AddEntry(h,label,"L");
  } 
  
  AddHist1D(h,drawopt,color,linesty,fillsty);
}

void CPlot::AddHist1D(TFile *f, TString histName, TString drawopt, int color, int linesty, int fillsty)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddHist1D(h,drawopt,color,linesty,fillsty);
}

void CPlot::AddHist1D(TFile *f, TString histName, TString label, TString drawopt, int color, int linesty, int fillsty)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddHist1D(h,label,drawopt,color,linesty,fillsty);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddToStack(TH1F *h, int color)
{
  if(!h)
    return;
    
  if(!fStack)
    fStack = new THStack(fName+TString("_stack"),"");
  
  fStack->Add(h);
  AddHist1D(h,"",color,1,1001);
}

void CPlot::AddToStack(TH1F *h, TString label, int color, int linecolor,
		       int linestyle, int linewidth, int fillstyle)
{
  if(!h)
    return;
  
  if(!fStack)
    fStack = new THStack(fName+TString("_stack"),"");
  
  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);  
    
  // make legend entries appear in reverse of the order the histograms are added
  if(label.Length()>0) {
    fStackEntries.push_back(fLeg->AddEntry(h,label,"F"));
    for(Int_t ientry=(fStackEntries.size()-2); ientry>=0; ientry--) {
      TObject* hh = fStackEntries[ientry]->GetObject();
      TString ll  = fStackEntries[ientry]->GetLabel();
      fStackEntries[ientry+1]->SetObject(hh);
      fStackEntries[ientry+1]->SetLabel(ll);
    }
    fStackEntries[0]->SetObject(h);
    fStackEntries[0]->SetLabel(label);
  }
     
  fStack->Add(h);  

  h->SetLineColor(linecolor);
  h->SetLineStyle(linestyle);
  h->SetLineWidth(linewidth);
  h->SetFillStyle(fillstyle);
  h->SetLabelSize(0.05);

  AddHist1D(h,"",color,linestyle,fillstyle,linecolor);
  // cout << h->GetLineColor() << " " << h->GetLineStyle() << " " << h->GetLineWidth() << endl;
}

void CPlot::AddToStack(TFile *f, TString histName, int color)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddToStack(h,color);
}

void CPlot::AddToStack(TFile *f, TString histName, TString label, int color)
{
  if(!f)
    return;
  
  TH1F *h = (TH1F*)f->FindObjectAny(histName);
  AddToStack(h,label,color);
}  
  
//--------------------------------------------------------------------------------------------------
void CPlot::AddHist2D(TH2D *h, TString drawopt, int fillcolor, int linecolor)
{
  if(!h)
    return;
  
  h->SetLineColor(linecolor);
  h->SetFillColor(fillcolor);
  h->SetMarkerStyle(kFullDotMedium); 
  
  CPlotItem item;
  item.hist2D = h;
  item.drawopt = drawopt;
  fItems.push_back(item);
}

void CPlot::AddHist2D(TFile *f, TString histName, TString drawopt, int fillcolor, int linecolor)
{
  if(!f)
    return;
  
  TH2D *h = (TH2D*)f->FindObjectAny(histName);
  AddHist2D(h,drawopt,linecolor,fillcolor);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddGraph(TGraph *gr, TString drawopt, int color, int marksty, int linesty)
{
  if(!gr)
    return;
  
  gr->SetMarkerColor(color);
  gr->SetLineColor(color);
  gr->SetLineStyle(linesty);
  gr->SetLineWidth(2);
  gr->SetMarkerStyle(marksty);
  gr->SetMarkerSize(1.2);

if(gr->GetN()==1) gr->SetMarkerSize(3.5);
  
  CPlotItem item;
  item.graph = gr;
  item.drawopt = drawopt;  
  fItems.push_back(item);
}

void CPlot::AddGraph(TGraph *gr, TString label, TString drawopt, int color, int marksty, int linesty)
{
  if(!gr)
    return;

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  if( (drawopt.Contains("L",TString::kIgnoreCase) ||
      drawopt.Contains("C",TString::kIgnoreCase))
      && !(drawopt.Contains("P",TString::kIgnoreCase)) ) {    
    fLeg->AddEntry(gr,label,"L");
  } else {  
    fLeg->AddEntry(gr,label,"LP"); 
  }
  
  AddGraph(gr,drawopt,color,marksty,linesty);
}

void CPlot::AddGraph(TFile *f, TString grName, TString drawopt, int color, int marksty, int linesty)
{
  if(!f)
    return;
  
  TGraph *gr = (TGraph*)f->FindObjectAny(grName);
  AddGraph(gr,drawopt,color,marksty,linesty);
}

void CPlot::AddGraph(TFile *f, TString grName, TString label, TString drawopt, int color, int marksty, int linesty)
{
  if(!f)
    return;
  
  TGraph *gr = (TGraph*)f->FindObjectAny(grName);
  AddGraph(gr,label,drawopt,color,marksty,linesty);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddProfile(TProfile *pr, TString drawopt, int color, int marksty, int linesty)
{
  if(!pr)
    return;
  
  pr->SetMarkerColor(color);
  pr->SetLineColor(color);
  pr->SetLineStyle(linesty);
  pr->SetLineWidth(2);
  pr->SetMarkerStyle(marksty);
  pr->SetMarkerSize(1.2);
  
  CPlotItem item;
  item.prof = pr;
  item.drawopt = drawopt;  
  fItems.push_back(item);
}

void CPlot::AddProfile(TProfile *pr, TString label, TString drawopt, int color, int marksty, int linesty)
{
  if(!pr)
    return;

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
    
  fLeg->AddEntry(pr,label,"LP");
  
  AddProfile(pr,drawopt,color,marksty,linesty);
}

void CPlot::AddProfile(TFile *f, TString prName, TString drawopt, int color, int marksty, int linesty)
{
  if(!f)
    return;
  
  TProfile *pr = (TProfile*)f->FindObjectAny(prName);
  AddProfile(pr,drawopt,color,marksty,linesty);
}

void CPlot::AddProfile(TFile *f, TString prName, TString label, TString drawopt, int color, int marksty, int linesty)
{
  if(!f)
    return;
  
  TProfile *pr = (TProfile*)f->FindObjectAny(prName);
  AddProfile(pr,label,drawopt,color,marksty,linesty);
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddTextBox(TString text, double x1, double y1, double x2, double y2,
                       int bordersize, int textcolor, int fillcolor)
{
  TPaveText *tb = new TPaveText(x1,y1,x2,y2,"NDC");
  tb->SetTextColor(textcolor);
  if(fillcolor==-1)
    tb->SetFillStyle(0);
  else
    tb->SetFillColor(fillcolor);
  tb->SetBorderSize(bordersize);
  tb->AddText(text);
  fTextBoxes.push_back(tb);
}

void CPlot::AddTextBox(double x1, double y1, double x2, double y2, 
                       int bordersize, int textcolor, int fillcolor, int nlines,...)
{
  TPaveText *tb = new TPaveText(x1,y1,x2,y2,"NDC");
  tb->SetTextColor(textcolor);
  if(fillcolor==-1)
    tb->SetFillStyle(0);
  else
    tb->SetFillColor(fillcolor);
  tb->SetBorderSize(bordersize);
  tb->SetTextAlign(12);
  
  va_list ap;
  va_start(ap,nlines);
  for(int i=0; i<nlines; i++) {
    TString textline(va_arg(ap,char*));
    tb->AddText(textline);
  }
  va_end(ap);
  
  fTextBoxes.push_back(tb);  
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddLine(double x1, double y1, double x2, double y2, int color, int style)
{
  TLine *line = new TLine(x1,y1,x2,y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(2);
  fLines.push_back(line);
}

void CPlot::AddLine(double x1, double y1, double x2, double y2, 
                    int color, int style, TString label)
{
  TLine *line = new TLine(x1,y1,x2,y2);
  line->SetLineColor(color);
  line->SetLineStyle(style);
  line->SetLineWidth(2);
  fLines.push_back(line);

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  fLeg->AddEntry(line,label,"L");
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddBox(double x1, double y1, double x2, double y2, 
                   int linecolor, int linesty, int fillcolor)
{
  TBox *box = new TBox(x1,y1,x2,y2);
  box->SetLineColor(linecolor);
  box->SetLineStyle(linesty);
  if(fillcolor==-1)
    box->SetFillStyle(0);
  else
    box->SetFillColor(fillcolor);
  box->SetLineWidth(2);
  fBoxes.push_back(box);
}

void CPlot::AddBox(double x1, double y1, double x2, double y2, 
                   int linecolor, int linesty, int fillcolor, TString label)
{
  TBox *box = new TBox(x1,y1,x2,y2);
  box->SetLineColor(linecolor);
  box->SetLineStyle(linesty);
  if(fillcolor==-1)
    box->SetFillStyle(0);
  else
    box->SetFillColor(fillcolor);
  box->SetLineWidth(2);
  fBoxes.push_back(box);

  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  
  if(fillcolor<0) fLeg->AddEntry(box,label,"L");
  else            fLeg->AddEntry(box,label,"F");
}

//--------------------------------------------------------------------------------------------------
void CPlot::AddFcn(TF1* fcn, int color, int linesty)
{
  if(!fcn)
    return;
  
  fcn->SetLineColor(color);
  fcn->SetLineStyle(linesty);
  fFcns.push_back(fcn);
}

void CPlot::AddFcn(TF1* fcn, TString label, int color, int linesty)
{
  if(!fcn)
    return;
    
  if(!fLeg)
    fLeg = new TLegend(0.6,0.84,0.93,0.9);
  else
    fLeg->SetY1(fLeg->GetY1()-0.06);
  fLeg->AddEntry(fcn,label,"L");
  
  AddFcn(fcn,color,linesty);
}

//-------------------------------------------------------------------------------------------------
void CPlot::DrawRatio(TCanvas *c, TH1F *h1, TH1F *h2, bool doSave, TString format)
{
  c->cd();
  TPad *p1 = new TPad("p1","p1",0,0.3,1,1);
  p1->SetLeftMargin  (0.18);
  p1->SetTopMargin   (0.10);
  p1->SetRightMargin (0.07);
  p1->SetBottomMargin(0.03);
  p1->Draw();
  p1->cd();
  c->SetLogx(fLogx);
  c->SetLogy(fLogy);
  double ymax = h1->GetMaximum();
  if(h2->GetMaximum()>ymax) ymax=h2->GetMaximum();
  h1->SetMarkerSize(1.3);
  h1->SetMarkerStyle(20);
  h1->SetLineWidth(3);
  h2->SetLineWidth(3);
  h2->GetXaxis()->SetTitle(fXTitle);
  h2->GetXaxis()->SetLabelOffset(0.20);
  h2->GetYaxis()->SetTitle(fYTitle);
  h2->GetYaxis()->SetRangeUser(0,1.1*ymax);
  h2->Draw("hist");
  h1->DrawCopy("same");
  if(fLeg) {
    fLeg->SetFillStyle(0);
    fLeg->SetBorderSize(0);
    fLeg->Draw();
  }
  c->cd();
  TPad *p2 = new TPad("p2","p2",0,0,1,0.3);
  p2->SetLeftMargin  (0.18);
  p2->SetTopMargin   (0.00);
  p2->SetRightMargin (0.07);
  p2->SetBottomMargin(0.30);
  p2->SetGridy(1);
  p2->Draw();
  p2->cd();
  h1->Divide(h2);
  h1->SetLineColor(kBlack);
  h1->GetXaxis()->SetTitle(fXTitle);
  h1->SetTitleSize  (0.12,"Y");
  h1->SetTitleOffset(0.60,"Y");
  h1->SetTitleSize  (0.12,"X");
  h1->SetLabelSize  (0.10,"X");
  h1->SetLabelSize  (0.08,"Y");
  h1->GetYaxis()->SetTitleFont(42);
  h1->GetYaxis()->CenterTitle(kTRUE);
  h1->GetXaxis()->SetTitleFont(42);
  h1->GetYaxis()->SetNdivisions(305);
  h1->GetXaxis()->SetTitle(fXTitle);
  h1->GetYaxis()->SetTitle("Ratio");
  h1->GetYaxis()->SetRangeUser(0.0,2.0);
  h1->Draw("EP");
  double xmin = h1->GetXaxis()->GetXmin();
  double xmax = h1->GetXaxis()->GetXmax();
  TLine *l = new TLine(xmin,1,xmax,1);
  l->SetLineWidth(3);
  l->Draw("same");
  c->cd();

  CMSPrelim("#tau_{e}#tau_{#mu}", 0.18, 0.93);


  if(doSave) {
    gSystem->mkdir(sOutDir,true);
    TString outname = sOutDir+TString("/")+fName+TString(".");
    if(format.CompareTo("all",TString::kIgnoreCase)==0) {
      c->SaveAs(outname+TString("png"));
      c->SaveAs(outname+TString("eps"));
      c->SaveAs(outname+TString("C"));
    } else {
      c->SaveAs(outname+format);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void CPlot::DrawRatioStack(TCanvas *c, TH1F *hData, TH1F *hMC, bool doSave, TString format)
{
  c->cd();

  TPad *p1 = new TPad("p1","p1",0,0.3,1,1);
  p1->SetLeftMargin  (0.18);
  p1->SetTopMargin   (0.10);
  p1->SetRightMargin (0.07);
  p1->SetBottomMargin(0.03);  
  p1->Draw();
  if(fLogy) p1->SetLogy();
  p1->cd();

  double ymax = hData->GetMaximum();
  if(hMC->GetMaximum()>ymax) ymax=hMC->GetMaximum();
  hData->SetMarkerSize(1.3);
  hData->SetMarkerStyle(20);
  hData->SetLineWidth(3);

  TH1F *h00 = (TH1F*)hData->Clone("data0");

  if(fStack) {
    if(!fLogy) {
      if(fYmin < fYmax) {
        fStack->SetMaximum(fYmax);
        fStack->SetMinimum(fYmin);
      } else {
        fStack->SetMaximum(1.1*ymax);
        fStack->SetMinimum(0);
      }
    }
    fStack->SetTitle(fTitle);
    fStack->Draw("hist");
    fStack->GetXaxis()->SetLabelOffset(0.20);
    fStack->GetYaxis()->SetTitle(fYTitle);
    if(fYmin < fYmax) fStack->GetYaxis()->SetRangeUser(fYmin,fYmax);
    else fStack->GetYaxis()->SetRangeUser(0,1.1*ymax);
    fStack->GetYaxis()->SetTitleSize(0.055);
    fStack->Draw("histsame");
  } 

  h00->Draw("same");


  if(fLeg) {
    fLeg->SetFillStyle(0);
    fLeg->SetBorderSize(0);
    fLeg->Draw();
  }

  c->cd();
  TPad *p2 = new TPad("p2","p2",0,0,1,0.3);
  p2->SetLeftMargin  (0.18);
  p2->SetTopMargin   (0.00);
  p2->SetRightMargin (0.07);
  p2->SetBottomMargin(0.30);
  p2->SetGridy(1);
  p2->Draw();
  p2->cd();

  TH1F *h1 = (TH1F*)hData->Clone("data1"); 
  TH1F *h2 = (TH1F*)hData->Clone("data2"); 
  TH1F *h3 = (TH1F*)hData->Clone("data3"); 
  TH1F *h4 = (TH1F*)hData->Clone("data4");

  h3->Add(hMC,-1);
  h3->SetLineColor(kBlack);

  h3->SetTitleSize  (0.12,"Y");
  h3->SetTitleOffset(0.60,"Y");
  h3->SetTitleSize  (0.12,"X");
  h3->SetLabelSize  (0.10,"X");
  h3->SetLabelSize  (0.08,"Y");
  h3->GetYaxis()->SetTitleFont(42);
  h3->GetYaxis()->CenterTitle(kTRUE);
  h3->GetXaxis()->SetTitleFont(42);
  h3->GetYaxis()->SetRangeUser(-5.0,5.0);
  h3->GetYaxis()->SetNdivisions(305);
  h3->GetXaxis()->SetTitle(fXTitle);
  h3->GetYaxis()->SetTitle("#frac{(N_{obs}-N_{exp})}{#sqrt{N_{obs}}}");


  double xmin = hData->GetXaxis()->GetXmin();
  double xmax = hData->GetXaxis()->GetXmax();
  TLine *l = new TLine(xmin,0,xmax,0);
  l->SetLineWidth(3);
  l->SetLineColor(kBlack);

  Int_t nbins = hData->GetNbinsX();
  
  for(int i = 0; i < nbins; i++) {
    double err=0;
    if(h1->GetBinContent(i+1)==0 ) err=0;
    else err=(0.5*h1->GetBinError(i+1)/sqrt(h1->GetBinContent(i+1)));
    h2->SetBinContent(i+1,sqrt(h1->GetBinContent(i+1)));
    h2->SetBinError(i+1,err);
  }

  h3->Divide(h2);

  for(int i = 0; i < nbins; i++) {
    h4->SetBinContent(i+1,0);
    h4->SetBinError(i+1,h3->GetBinError(i+1));
    h3->SetBinError(i+1,0);
    if(h3->GetBinContent(i+1)==0) h3->SetBinContent(i+1,-10);      // hack, because I can't figure out how to suppress the points for bins with zeroes
  }

  TGraphErrors *grErr = new TGraphErrors(h4);
  grErr->SetFillColor(kPink-5);
  grErr->SetFillStyle(3001);


  h3->DrawCopy("P");  
  grErr->Draw("same2");
  h3->Draw("Psame");
  l->Draw("same");

  c->cd();

  CMSPrelim("#tau_{e}#tau_{#mu}", 0.18, 0.93);

  if(doSave) {
    gSystem->mkdir(sOutDir,true);
    TString outname = sOutDir+TString("/")+fName+TString(".");
    if(format.CompareTo("all",TString::kIgnoreCase)==0) {
      c->SaveAs(outname+TString("png"));
      c->SaveAs(outname+TString("eps"));
      c->SaveAs(outname+TString("C"));
    } else {
      c->SaveAs(outname+format);
    }
  }
}


//--------------------------------------------------------------------------------------------------
void CPlot::Draw(TCanvas *c, bool doSave, TString format)
{ 
  c->cd();
  
  c->SetLogy(0);
  c->SetLogx(0);
  
  if(!fItems.size())
    return;   
  
  int nHist1D=0, nHist2D=0, nGraph=0, nProf=0;
  for(uint i=0; i<fItems.size(); i++) {
    if(fItems[i].hist1D != 0) nHist1D++;
    if(fItems[i].hist2D != 0) nHist2D++;
    if(fItems[i].graph  != 0) nGraph++;
    if(fItems[i].prof != 0) nProf++;
  }
  
  //
  // Draw 2D histogram, save if necessary, then exit
  //   Suggested options for:
  //     contour plot -> "CONT4Z"
  //     lego plot    -> "LEGO1 0"
  //     color plot   -> "COLZ"
  //   Default is scatter plot
  //  
  if(nHist2D>0) {
    for(uint i=0; i<fItems.size(); i++) {
      if(fItems[i].hist2D==0) continue;
      
      fItems[i].hist2D->Draw(fItems[i].drawopt);
      fItems[i].hist2D->SetTitle(fTitle);
      fItems[i].hist2D->GetXaxis()->SetTitle(fXTitle);
      fItems[i].hist2D->GetYaxis()->SetTitle(fYTitle);
    
      //
      // Set log scale if necessary
      // 
      c->SetLogx(fLogx);
      c->SetLogy(fLogy);
      
      for(uint k=0; k<fLines.size(); k++)
        fLines[k]->Draw();

      for(uint k=0; k<fBoxes.size(); k++)
        fBoxes[k]->Draw();
      
      for(uint j=0; j<fTextBoxes.size(); j++)
        fTextBoxes[j]->Draw();
            
      if(doSave) {
        gSystem->mkdir(sOutDir,true);
        TString outname = sOutDir+TString("/")+fName+TString(".");
	if(format.CompareTo("all",TString::kIgnoreCase)==0) {
	  c->SaveAs(outname+TString("pdf"));
	  c->SaveAs(outname+TString("png"));
	  c->SaveAs(outname+TString("eps"));	  
	  c->SaveAs(outname+TString("C"));
	} else {
	  c->SaveAs(outname+format);
	}
      }
      
      return;
    }
  }    
  
  // 
  // Draw 1D histograms
  //   Histograms are cloned so that content and properties 
  //   of the original histograms are not changed
  //
  std::vector<TH1F*> vHists;
  std::vector<TString> vHistOpts;
  if(nHist1D>0) {   
    
    double ymax=0;
    uint ifirst=0;
    
    for(uint i=0; i<fItems.size(); i++) {
      if(fItems[i].hist1D==0) continue;
      if(fStack && fStack->GetHists()->Contains(fItems[i].hist1D)) continue;
      
      TString hname = fName;
      hname += "_h_";
      hname += i;
    
      TH1F *h;     
      if(fRebin>1)      
        h = (TH1F*)fItems[i].hist1D->Rebin(fRebin,hname);
      else
        h = (TH1F*)fItems[i].hist1D->Clone(hname);
  
      if(fXmin < fXmax) {
        h->GetXaxis()->SetRangeUser(fXmin,fXmax);
      }
      
      if(fYmin < fYmax) { 
        h->GetYaxis()->SetRangeUser(fYmin,fYmax);
      } else {
        if(ymax < h->GetMaximum()) {
	  ymax = h->GetMaximum();
	  ifirst = vHists.size();
	}
      }
      
      vHists.push_back(h);
      vHistOpts.push_back(fItems[i].drawopt);
    }
    
    if(vHists.size()>0) {
      vHists[ifirst]->SetTitle(fTitle);
      vHists[ifirst]->GetXaxis()->SetTitle(fXTitle);
      vHists[ifirst]->GetYaxis()->SetTitle(fYTitle);
      vHists[ifirst]->SetLineWidth(2);
      vHists[ifirst]->Draw(vHistOpts[ifirst].Data());
  }
   
    //
    // Draw histogram stack
    //
    if(fStack) {
      if(vHists.size()>0) { 
        if(fYmin < fYmax) {
	  fStack->Draw("hist same");
	} else {
	  if(fStack->GetMaximum() > vHists[ifirst]->GetMaximum()) {
	    fStack->SetTitle(fTitle);
	    fStack->Draw();
	    fStack->GetXaxis()->SetTitle(fXTitle);
	    fStack->GetYaxis()->SetTitle(fYTitle);
	    fStack->Draw("hist");
	  } else {	    
	    fStack->Draw("hist same"); 
	  }
        }
	
      } else {
        // NOTE: Must draw first before accessing axes
	fStack->Draw("hist"); 

	if(fXmin < fXmax)
          fStack->GetXaxis()->SetRangeUser(fXmin,fXmax);
	 
	if(fYmin < fYmax) {
          fStack->SetMaximum(fYmax);  
	  fStack->SetMinimum(fYmin);
        }
	 
        fStack->SetTitle(fTitle);
	fStack->GetXaxis()->SetTitle(fXTitle);
	fStack->GetYaxis()->SetTitle(fYTitle);	
        fStack->Draw("hist");	 
      } 
    }
        
    for(uint i=0; i<vHists.size(); i++) {
      TH1F *h = vHists[i];              
      h->SetLineWidth(2);
      char opt[100];
      sprintf(opt,"same%s",vHistOpts[i].Data());
      h->Draw(opt);
    }
  }  
  c->RedrawAxis();
  
  //
  // Draw graphs
  //
  std::vector<TGraph*> vGraphs;
  std::vector<TString> vGraphOpts;
  if(nGraph>0) {    
    for(uint i=0; i<fItems.size(); i++) {
      if(fItems[i].graph==0) continue;
    
      TString grName = fName;
      grName += "_gr_";
      grName += i;
      
      TGraph *gr = (TGraph*)fItems[i].graph->Clone(grName);
      
      if(fXmin < fXmax)
        gr->GetXaxis()->SetLimits(fXmin,fXmax);
//        gr->GetXaxis()->SetRangeUser(fXmin,fXmax);
    
      if(fYmin < fYmax)
        gr->GetYaxis()->SetRangeUser(fYmin,fYmax);
	
      vGraphs.push_back(gr);
      vGraphOpts.push_back(fItems[i].drawopt);
    }
    
    if(vHists.size()==0) {
      vGraphs[0]->SetTitle(fTitle);
      vGraphs[0]->GetXaxis()->SetTitle(fXTitle);
      vGraphs[0]->GetYaxis()->SetTitle(fYTitle);
    }
    
    for(uint i=0; i<vGraphs.size(); i++) {
      TGraph *gr = vGraphs[i];
      char opt[100];
      (i==0 && nHist1D==0) ? sprintf(opt,"AP%s",vGraphOpts[i].Data()) : sprintf(opt,"P%s",vGraphOpts[i].Data());
      gr->Draw(opt);
    }
  }

  //
  // Draw profile histograms
  //
  std::vector<TProfile*> vProfiles;
  std::vector<TString> vProfileOpts;
  if(nProf>0) {    
    for(uint i=0; i<fItems.size(); i++) {
      if(fItems[i].prof==0) continue;
    
      TString prName = fName;
      prName += "_pr_";
      prName += i;
      
      TProfile *pr = (TProfile*)fItems[i].prof->Clone(prName);
      
      if(fXmin < fXmax)
        pr->GetXaxis()->SetLimits(fXmin,fXmax);
//        pr->GetXaxis()->SetRangeUser(fXmin,fXmax);
    
      if(fYmin < fYmax)
        pr->GetYaxis()->SetRangeUser(fYmin,fYmax);
	
      vProfiles.push_back(pr);
      vProfileOpts.push_back(fItems[i].drawopt);
    }
    
    if(vHists.size()==0) {
      vProfiles[0]->SetTitle(fTitle);
      vProfiles[0]->GetXaxis()->SetTitle(fXTitle);
      vProfiles[0]->GetYaxis()->SetTitle(fYTitle);
    }
    
    for(uint i=0; i<vProfiles.size(); i++) {
      TProfile *pr = vProfiles[i];
      char opt[100];
      if(i>0 || nHist1D>0 || nGraph>0) 
        sprintf(opt,"same%s",vProfileOpts[i].Data());
      else
        sprintf(opt,"%s",vProfileOpts[i].Data());
      pr->Draw(opt);
    }
  }
      
  //
  // Draw legend
  //
  if(fLeg) {
    fLeg->SetFillStyle(0);
    fLeg->SetFillColor(0);
    fLeg->SetBorderSize(0);
    fLeg->Draw();
  }
  
  //
  // Draw statistics box
  //
  TLatex *stat=0, *sval=0;
  if(fShowStats) {
    char buffer[20];
    stat = new TLatex[3*vHists.size()];
    sval = new TLatex[3*vHists.size()];
    for(uint i=0; i<vHists.size(); i++) {
      int x = fShowStats;
      
      // number of entries
      if(x / 100) {
        stat[3*i].SetNDC(); stat[3*i].SetTextAlign(13); stat[3*i].SetTextSize(0.03);
        stat[3*i].SetText(fStatsX,fStatsY-0.04*(3*i)-0.005*i,"Entries");
	stat[3*i].SetTextColor(vHists[i]->GetLineColor());
        stat[3*i].Draw();   
        sprintf(buffer,"%i",int(vHists[i]->Integral()));
        sval[3*i].SetNDC(); sval[3*i].SetTextAlign(33); sval[3*i].SetTextSize(0.03);
        sval[3*i].SetText(fStatsX+0.25,fStatsY-0.04*(3*i)-0.005*i,buffer);
	sval[3*i].SetTextColor(vHists[i]->GetLineColor());
        sval[3*i].Draw();
      }
      
      // mean
      x = x % 100;
      if(x / 10) {
        stat[3*i+1].SetNDC(); stat[3*i+1].SetTextAlign(13); stat[3*i+1].SetTextSize(0.03);
        stat[3*i+1].SetText(fStatsX,fStatsY-0.04*(3*i+1)-0.005*i,"Mean");
	stat[3*i+1].SetTextColor(vHists[i]->GetLineColor());
        stat[3*i+1].Draw();   
        sprintf(buffer,"%g",vHists[i]->GetMean());
        sval[3*i+1].SetNDC(); sval[3*i+1].SetTextAlign(33); sval[3*i+1].SetTextSize(0.03);
        sval[3*i+1].SetText(fStatsX+0.25,fStatsY-0.04*(3*i+1)-0.005*i,buffer);
	sval[3*i+1].SetTextColor(vHists[i]->GetLineColor());
        sval[3*i+1].Draw();
      }
      
      // RMS
      x = x % 10;
      if(x) {
        stat[3*i+2].SetNDC(); stat[3*i+2].SetTextAlign(13); stat[3*i+2].SetTextSize(0.03);
        stat[3*i+2].SetText(fStatsX,fStatsY-0.04*(3*i+2)-0.005*i,"RMS");
	stat[3*i+2].SetTextColor(vHists[i]->GetLineColor());
        stat[3*i+2].Draw();   
        sprintf(buffer,"%g",vHists[i]->GetRMS());
        sval[3*i+2].SetNDC(); sval[3*i+2].SetTextAlign(33); sval[3*i+2].SetTextSize(0.03);
        sval[3*i+2].SetText(fStatsX+0.25,fStatsY-0.04*(3*i+2)-0.005*i,buffer);
	sval[3*i+2].SetTextColor(vHists[i]->GetLineColor());
        sval[3*i+2].Draw();
      }
    }
  }
 
  //
  // Draw functions
  //
  for(uint i=0; i<fFcns.size(); i++)
    (i==0 && vHists.size()==0 && vGraphs.size()==0) ? fFcns[i]->Draw() : fFcns[i]->Draw("sameC");
  
  //
  // Draw lines
  //
  for(uint i=0; i<fLines.size(); i++)
    fLines[i]->Draw();
  
  //
  // Draw Boxes
  //
  for(uint i=0; i<fBoxes.size(); i++)
    fBoxes[i]->Draw();
  
  //
  // Draw textboxes
  //
  for(uint i=0; i<fTextBoxes.size(); i++)
    fTextBoxes[i]->Draw();    
  
  //
  // Set log scale if necessary
  // 
  c->SetLogx(fLogx);
  c->SetLogy(fLogy);
  
  //
  // Set grid lines if necessary
  //
  c->SetGridx(fGridx);
  c->SetGridy(fGridy);

  //
  // Add CMS Preliminary and lumi text
  //
  CMSPrelim("#tau_{e}#tau_{#mu}", 0.18, 0.92);

  //
  // Save plot if necessary
  //  
  if(doSave) {
    gSystem->mkdir(sOutDir,true);
    TString outname = sOutDir+TString("/")+fName+TString(".");
    if(format.CompareTo("all",TString::kIgnoreCase)==0) {
      c->SaveAs(outname+TString("png"));
      c->SaveAs(outname+TString("eps"));
      c->SaveAs(outname+TString("C"));
    } else {
      c->SaveAs(outname+format);
    }
    
    delete [] stat;
    delete [] sval;
//    for(uint i=0; i<vHists.size(); i++)
//      delete vHists[i];
  }
}

void CPlot::CMSPrelim(const char* channel, double lowX, double lowY)
{
  TPaveText* cmsprel  = new TPaveText(lowX, lowY, lowX+0.25, lowY+0.04, "NDC");
  cmsprel->SetBorderSize(   0 );
  cmsprel->SetFillStyle(    0 );
  cmsprel->SetTextAlign(   11 );
  cmsprel->SetTextSize ( 0.035);
  cmsprel->SetTextColor(    1 );
  cmsprel->SetTextFont (   42 );
  cmsprel->AddText("CMS Data 2011, 4.6 fb^{-1}, #sqrt{s} = 7 TeV");
  cmsprel->Draw();

  TPaveText* lumi     = new TPaveText(lowX+0.35, lowY, lowX+0.55, lowY+0.04, "NDC");
  lumi->SetBorderSize(   0 );
  lumi->SetFillStyle(    0 );
  lumi->SetTextAlign(   11 );
  lumi->SetTextSize ( 0.035);
  lumi->SetTextColor(    1 );
  lumi->SetTextFont (   42 );
  lumi->AddText("4.6 fb^{-1}  #sqrt{s}=7 TeV");
  //lumi->Draw();

  TPaveText* chan     = new TPaveText(lowX+0.68, lowY, lowX+0.73, lowY+0.04, "NDC");
  chan->SetBorderSize(   0 );
  chan->SetFillStyle(    0 );
  chan->SetTextAlign(   11 );
  chan->SetTextSize ( 0.035);
  chan->SetTextColor(    1 );
  chan->SetTextFont (   42 );
  chan->AddText(channel);
  chan->Draw();
}
