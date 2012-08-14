#ifndef LIMITINPUTSEMU_HH
#define LIMITINPUTSEMU_HH

Double_t rescale(TH1F* hist, TString catname, TString histname, TFile* file){

  TH1F* hcenter = 0;
  
  TDirectory *dir = (TDirectory*)file->FindObjectAny(catname);
  hcenter = (TH1F*)(dir->Get(histname)); assert(hcenter);

  Double_t center = hcenter->Integral(0,hist->GetNbinsX()+1);
  Double_t integral = hist->Integral(0,hist->GetNbinsX()+1);

  Double_t scale = 1.0*center/integral;

  return scale;
}

void makeLimitInputs(Bool_t domssm, TString output, UInt_t ecorr, vector<TString> snamev, vector<TH1F*> hMassL_0jet_lowptv, vector<TH1F*> hMassL_0jet_highptv, vector<TH1F*> hMassL_boost_lowptv, vector<TH1F*> hMassL_boost_highptv, vector<TH1F*> hMassL_b_lowptv, vector<TH1F*> hMassL_b_highptv, vector<TH1F*> hMassL_vbfv) {

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

  TFile *flimit;
  if (ecorr == 0) {
    flimit = new TFile(output,"recreate");
    flimit->mkdir("emu_0jet_low");
    flimit->mkdir("emu_0jet_high");
    flimit->mkdir("emu_boost_low");
    flimit->mkdir("emu_boost_high");
    flimit->mkdir("emu_btag_low");
    flimit->mkdir("emu_btag_high");
    flimit->mkdir("emu_vbf");
  }
  else flimit = new TFile(output,"update");

  // add zmm into the fakes histogram
    hMassL_0jet_lowptv[iewk_7TeV]   ->Add(hMassL_0jet_lowptv[izmm]);
    hMassL_0jet_highptv[iewk_7TeV]  ->Add(hMassL_0jet_highptv[izmm]);
    hMassL_boost_lowptv[iewk]  ->Add(hMassL_boost_lowptv[izmm]);
    hMassL_boost_highptv[iewk_7TeV] ->Add(hMassL_boost_highptv[izmm]);
    hMassL_b_lowptv[iewk]      ->Add(hMassL_b_lowptv[izmm]);
    hMassL_b_highptv[iewk]     ->Add(hMassL_b_highptv[izmm]);
    hMassL_vbfv[iewk]		  ->Add(hMassL_vbfv[izmm]);

  double fakebin = hMassL_vbfv[ifake]->GetBinContent(5);
  double fakebinerr = hMassL_vbfv[ifake]->GetBinError(5);
  hMassL_vbfv[ifake]->SetBinContent(5,0.3*fakebin);
  hMassL_vbfv[ifake]->SetBinError(5,0.3*fakebinerr);
  hMassL_vbfv[ifake]->SetBinContent(3,hMassL_vbfv[ifake]->GetBinContent(3)+0.15*fakebin);
  hMassL_vbfv[ifake]->SetBinError(3,hMassL_vbfv[ifake]->GetBinError(3)+0.15*fakebinerr);
  hMassL_vbfv[ifake]->SetBinContent(4,hMassL_vbfv[ifake]->GetBinContent(4)+0.2*fakebin);
  hMassL_vbfv[ifake]->SetBinError(4,hMassL_vbfv[ifake]->GetBinError(4)+0.2*fakebinerr);
  hMassL_vbfv[ifake]->SetBinContent(6,hMassL_vbfv[ifake]->GetBinContent(6)+0.25*fakebin);
  hMassL_vbfv[ifake]->SetBinError(6,hMassL_vbfv[ifake]->GetBinError(6)+0.25*fakebinerr);
  hMassL_vbfv[ifake]->SetBinContent(7,hMassL_vbfv[ifake]->GetBinContent(7)+0.1*fakebin);
  hMassL_vbfv[ifake]->SetBinError(7,hMassL_vbfv[ifake]->GetBinError(7)+0.1*fakebinerr);

  TString histname;
  for(UInt_t isam=0;isam<snamev.size();isam++) {

    if(isam==izmm) continue;
    if(!domssm && snamev[isam].Contains("mssm")) continue;
    else if(domssm && snamev[isam].Contains("_sm")) continue;

    if(snamev[isam].Contains("ewk",   TString::kIgnoreCase))            histname = "EWK";
    else if(snamev[isam].Contains("fakes", TString::kIgnoreCase))       histname = "Fakes";
    else if(snamev[isam].Contains("ttbar", TString::kIgnoreCase))       histname = "ttbar";
    else if(snamev[isam].Contains("Zmm", TString::kIgnoreCase))         histname = "Zmm";
    else if(snamev[isam].Contains("Ztt",   TString::kIgnoreCase))       histname = "Ztt";
    else if(snamev[isam].Contains("emb",   TString::kIgnoreCase))       histname = "Ztt";
    else if((snamev[isam].Contains("data",  TString::kIgnoreCase)) && !(snamev[isam].Contains("emb",   TString::kIgnoreCase)))          {histname = "data_obs";}
    else if(snamev[isam].Contains("htt_")) continue;
    else if(snamev[isam].Contains("gf_sm")) histname = snamev[isam].ReplaceAll("gf_sm_","ggH");
    else if(snamev[isam].Contains("vbf_sm")) histname = snamev[isam].ReplaceAll("vbf_sm_","qqH");
    else if(snamev[isam].Contains("vtth_sm")) histname = snamev[isam].ReplaceAll("vtth_sm_","VH");
    else if(snamev[isam].Contains("gg_mssm")) histname = snamev[isam].ReplaceAll("gg_mssm_","ggH");
    else if(snamev[isam].Contains("bb_mssm")) histname = snamev[isam].ReplaceAll("bb_mssm_","bbH");
    else { cout << "error! name not found" << endl; assert(0); }

    /*if(hMassL_vbfv[isam]->Integral(0,hMassL_vbfv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_vbfv[isam]->SetBinContent(ibin,0.00001);
      hMassL_vbfv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_bv[isam]->Integral(0,hMassL_bv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_bv[isam]->SetBinContent(ibin,0.00001);
      hMassL_bv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_b_lowptv[isam]->Integral(0,hMassL_b_lowptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_b_lowptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_b_lowptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_b_highptv[isam]->Integral(0,hMassL_b_highptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_b_highptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_b_highptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_boostv[isam]->Integral(0,hMassL_boostv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_boostv[isam]->SetBinContent(ibin,0.00001);
      hMassL_boostv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_boost_lowptv[isam]->Integral(0,hMassL_boost_lowptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_boost_lowptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_boost_lowptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_boost_highptv[isam]->Integral(0,hMassL_boost_highptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_boost_highptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_boost_highptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_0jetv[isam]->Integral(0,hMassL_0jetv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_0jetv[isam]->SetBinContent(ibin,0.00001);
      hMassL_0jetv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_0jet_lowptv[isam]->Integral(0,hMassL_0jet_lowptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_0jet_lowptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_0jet_lowptv[isam]->SetBinError(ibin,0.00001);
    }
    if(hMassL_0jet_highptv[isam]->Integral(0,hMassL_0jet_highptv[isam]->GetNbinsX()+1) == 0) {
      int ibin = int(randm.Uniform(1,massLNbins_SM));
      hMassL_0jet_highptv[isam]->SetBinContent(ibin,0.00001);
      hMassL_0jet_highptv[isam]->SetBinError(ibin,0.00001);
    }*/

    if (isam>0 && isam!=iewk && isam!=iewk_7TeV && isam !=ittbar && isam!=ittbar_7TeV && isam!=ifake && isam!=issfake) {

      if (ecorr==1) {
        hMassL_0jet_lowptv[isam]->Scale(rescale(hMassL_0jet_lowptv[isam],"emu_0jet_low",histname,flimit));
        hMassL_0jet_highptv[isam]->Scale(rescale(hMassL_0jet_highptv[isam],"emu_0jet_high",histname,flimit));
        hMassL_boost_lowptv[isam]->Scale(rescale(hMassL_boost_lowptv[isam],"emu_boost_low",histname,flimit));
        hMassL_boost_highptv[isam]->Scale(rescale(hMassL_boost_highptv[isam],"emu_boost_high",histname,flimit));
        hMassL_b_lowptv[isam]->Scale(rescale(hMassL_b_lowptv[isam],"emu_btag_low",histname,flimit));
        hMassL_b_highptv[isam]->Scale(rescale(hMassL_b_highptv[isam],"emu_btag_high",histname,flimit));
        hMassL_vbfv[isam]->Scale(rescale(hMassL_vbfv[isam],"emu_vbf",histname,flimit));

        histname += "_CMS_scale_e_8TeVDown";
      }
      if (ecorr==2) {
        hMassL_0jet_lowptv[isam]->Scale(rescale(hMassL_0jet_lowptv[isam],"emu_0jet_low",histname,flimit));
        hMassL_0jet_highptv[isam]->Scale(rescale(hMassL_0jet_highptv[isam],"emu_0jet_high",histname,flimit));
        hMassL_boost_lowptv[isam]->Scale(rescale(hMassL_boost_lowptv[isam],"emu_boost_low",histname,flimit));
        hMassL_boost_highptv[isam]->Scale(rescale(hMassL_boost_highptv[isam],"emu_boost_high",histname,flimit));
        hMassL_b_lowptv[isam]->Scale(rescale(hMassL_b_lowptv[isam],"emu_btag_low",histname,flimit));
        hMassL_b_highptv[isam]->Scale(rescale(hMassL_b_highptv[isam],"emu_btag_high",histname,flimit));
        hMassL_vbfv[isam]->Scale(rescale(hMassL_vbfv[isam],"emu_vbf",histname,flimit));

        histname += "_CMS_scale_e_8TeVUp";
      }

    } else {if (ecorr!=0) continue; }

    if (isam!=ifake && isam!=ittbar_7TeV) {
      if(isam!=iewk) {
	flimit->cd("emu_0jet_low");
	hMassL_0jet_lowptv[isam]->SetName(histname);
	hMassL_0jet_lowptv[isam]->Write();
	flimit->cd("emu_0jet_high");
	hMassL_0jet_highptv[isam]->SetName(histname);
	hMassL_0jet_highptv[isam]->Write();
      }
      if(isam!=iewk_7TeV) {
	flimit->cd("emu_boost_low");
	hMassL_boost_lowptv[isam]->SetName(histname);
	hMassL_boost_lowptv[isam]->Write();
      }
      if(isam!=iewk) {
	flimit->cd("emu_boost_high");
	hMassL_boost_highptv[isam]->SetName(histname);
	hMassL_boost_highptv[isam]->Write();
      }
      if(isam!=iewk_7TeV) {
	flimit->cd("emu_btag_low");
	hMassL_b_lowptv[isam]->SetName(histname);
	hMassL_b_lowptv[isam]->Write();
	flimit->cd("emu_btag_high");
	hMassL_b_highptv[isam]->SetName(histname);
	hMassL_b_highptv[isam]->Write();
      }
    }
    if (isam!=issfake && isam!=ittbar && isam!=iewk_7TeV) {
      flimit->cd("emu_vbf");
      hMassL_vbfv[isam]->SetName(histname);
      hMassL_vbfv[isam]->Write();
    }

  }

  flimit->Close();
}

#endif
