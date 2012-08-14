#ifndef PLOTSEMU_HH
#define PLOTSEMU_HH

void makePlot(TString outdirName, TString format, TString name, TString title, TString xtitle, TString ytitle, vector<TString> snamev, vector<CSample*> samplev, vector<TH1F*> histv, Bool_t useFR=kFALSE, Bool_t makeLog=kFALSE, Double_t yMin = 0.2, Double_t yMax = 30.0) {

  // get indices of samples
  UInt_t ifake=9999, izmm=9999, iztt=9999, ittbar=9999, imssm_gg=9999, imssm_bb=9999, ism_vbf=9999, ism_gf=9999, ism_vtth=9999, issfake=9999, iemb=9999, iewk=9999, iewk_7TeV=9999, ittbar_7TeV=9999;
  for(UInt_t isam=0; isam<snamev.size(); isam++) {
    if(snamev[isam].Contains("fakes") && !snamev[isam].Contains("ss-fakes")) ifake = isam;
    if(snamev[isam].Contains("ss-fakes")) issfake = isam;
    if(snamev[isam].Contains("zmm"))   izmm  = isam;
    if(snamev[isam].Contains("ztt"))   iztt = isam;
    if(snamev[isam].Contains("emb"))   iemb = isam;
    if(snamev[isam].Contains("ttbar-8TeV")) ittbar = isam;
    if(snamev[isam].Contains("ttbar-7TeV")) ittbar_7TeV = isam;
    if(snamev[isam].Contains("ewk-8TeV"))   iewk = isam;
    if(snamev[isam].Contains("ewk-7TeV"))   iewk_7TeV = isam;
    if(snamev[isam].Contains("htt_gg_mssm")) imssm_gg=isam;
    if(snamev[isam].Contains("htt_bb_mssm")) imssm_bb=isam;
    if(snamev[isam].Contains("htt_vbf_sm")) ism_vbf=isam;
    if(snamev[isam].Contains("htt_gf_sm")) ism_gf=isam;
    if(snamev[isam].Contains("htt_vtth_sm")) ism_vtth=isam;
  }

  SetStyle();
  TCanvas *c = MakeCanvas("c","c",700,700);
  
  CPlot *plot = new CPlot(name,title,xtitle,ytitle);
  plot->sOutDir = outdirName;
  gSystem->mkdir(outdirName,kTRUE);
  if(!makeLog && histv[iewk] && histv[izmm]) histv[iewk]->Add(histv[izmm]);
  if(!makeLog && histv[ism_vbf] && histv[ism_gf] && histv[ism_vtth]) {
    histv[ism_vbf]->Add(histv[ism_gf]);
    histv[ism_vbf]->Add(histv[ism_vtth]);
  }
  if(ytitle.Contains("dN")) histv[0]->Scale(1.,"width");
  plot->AddHist1D(histv[0],samplev[0]->label,"E");
  for(UInt_t isam=1; isam<samplev.size(); isam++) {
    if((isam==izmm) || (isam==iewk_7TeV) || (isam==ittbar_7TeV) || (useFR && isam==issfake) || (!useFR && isam==ifake)) continue;
    if(ytitle.Contains("dN")) histv[isam]->Scale(1.,"width");
    if(snamev[isam].Contains("htt_vbf")) {
      if(!makeLog) {
	histv[isam]->Scale(5.);
	plot->AddToStack(histv[isam],"(5#times) "+samplev[isam]->label,samplev[isam]->color,603,11,3,0);
      } else {
	plot->AddHist1D(histv[isam],"(5#times) "+samplev[isam]->label,"hist",603,11,0);
      }
    }
    else if(snamev[isam].Contains("mssm_") || snamev[isam].Contains("sm_") || snamev[isam].Contains("ggH") || snamev[isam].Contains("bbH") || snamev[isam].Contains("qqH") || snamev[isam].Contains("VH")) continue;
    else plot->AddToStack(histv[isam],samplev[isam]->label,samplev[isam]->color,1,1,3);
  }
  if(makeLog) {
    plot->SetYRange(yMin,yMax*plot->GetStack()->GetMaximum());
    plot->SetLogy();
  }
  plot->SetLegend(0.5,0.65,0.95,0.9);
  plot->Draw(c,kTRUE,format);

  delete c;
  delete plot;

}



#endif
