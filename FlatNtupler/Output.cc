#include "Output.hh"
#include "MitHtt/Common/MyTools.hh"        // miscellaneous helper functions


Output::Output(TString name):
  fRun(0),
  fLumi(0),
  fEvt(0),
  fpvx(0),
  fpvy(0),
  fpvz(0),
  fNPV(0),
  fNPU(0),
  fRho(0),
  fMCWeight(0),
  fPUWeight(0),
  fEffWeight(0),
  fWeight(0),
  fEmbWeight(0),
  fEmbKinWeight(0),
  fEmbIdScale(0),
  fEmbGenWeight(0), 
  fEmbSpinnerWeight(0),
  fEmbMuEff(0),
  fEmbMuRad(0),
  fEmbKinMassWeight(0),
  fMass(0),
  fMVis(0),
  fPt1(0),
  fPhi1(0),
  fEta1(0),
  fM1(0),
  fq1(0),
  fIso1(0),
  fD01(0),
  fDZ1(0),
  fPassIso1(0),
  fMt1(0),
  fMVAMt1(0),
  fngamma1(0),
  fnprong1(0),
  fantiele1(0),
  fPt2(0),
  fPhi2(0),
  fEta2(0),
  fM2(0),
  fq2(0),
  fIso2(0),
  fD02(0),
  fDZ2(0),
  fPassIso2(0),
  fMt2(0),
  fMVAMt2(0),
  fngamma2(0),
  fnprong2(0),
  fantiele2(0),
  fdrll(0),
  fMet(0),
  fMetPhi(0),
  fMVAMet(0),
  fMVAMetPhi(0),
  fPZetaVis(0),
  fPZetaMiss(0),
  fPZetaMVAMiss(0),
  fMetCov00(0),
  fMetCov01(0),
  fMetCov10(0),
  fMetCov11(0),
  fMVACov00(0),
  fMVACov01(0),
  fMVACov10(0),
  fMVACov11(0),
  fJPt1(0),
  fJEta1(0),
  fJPhi1(0),
  fJM1(0),
  fJPtUnc1(0),
  fJMVA1(0),
  fJcsv1(0),
  fJPass1(0),
  fJPt2(0),
  fJEta2(0),
  fJPhi2(0),
  fJM2(0),
  fJPtUnc2(0),
  fJMVA2(0),
  fJcsv2(0),
  fJPass2(0),
  fBTagPt1(0),
  fBTagEta1(0),
  fBTagPhi1(0),
  fBTagM1(0),
  fbcsv1(0),
  fBTagPt2(0),
  fBTagEta2(0),
  fBTagPhi2(0),
  fBTagM2(0),
  fbcsv2(0),
  fMJJ(0),
  fJDEta(0),
  fNJetInGap(0),
  fMVA(0),
  fJDPhi(0),
  fDiJetPt(0),
  fDiJetPhi(0),
  fHDJetPhi(0),
  fVisJetEta(0),
  fPtVis(0),
  fPtH(0),
  fPtHMVA(0),
  fNBTag(0),
  fNJets(0),
  fNJetsClean(0),
  fGenPt1(0),
  fGenPhi1(0),
  fGenEta1(0),
  fGenId1(0),
  fGenPt2(0),
  fGenPhi2(0),
  fGenEta2(0),
  fGenId2(0),
  fVisGenPt1(0),
  fVisGenPhi1(0),
  fVisGenEta1(0),
  fVisGenId1(0),
  fVisGenPt2(0),
  fVisGenPhi2(0),
  fVisGenEta2(0),
  fVisGenId2(0),
  byCombinedIsolationDeltaBetaCorrRaw3Hits1(0),
  byCombinedIsolationDeltaBetaCorrRaw3Hits2(0),
  againstElectronMVA3raw1(0),
  againstElectronMVA3raw2(0),
  byIsolationMVA2raw1(0),
  byIsolationMVA2raw2(0),
  againstMuonLoose21(0),
  againstMuonLoose22(0),
  againstMuonMedium21(0),
  againstMuonMedium22(0),
  againstMuonTight21(0),
  againstMuonTight22(0),
  passAntiEleMVA31(0),
  passAntiEleMVA32(0),
  passAntiEleNewWPMVA31(0),
  passAntiEleNewWPMVA32(0),
  antitightele1(0),
  antitightele2(0),
  fpartons(0),
  fGenMatch(0),
  fLepMatchAny(0),
  fGenMatchLep(0),
  doRecoil(0),
  doEmu(0),
  npt20jets(0),
  fOutputFile(0),
  fEventTree(0),
  lshift1(0),
  lshift2(0),
  corrector(0)
{
  btagArray.Set(50);
  jptArray.Set(50);
  jetaArray.Set(50);

  setupOutput(name);
}
void Output::cd() { 
  fOutputFile->cd();
}
void Output::setupRecoil(int doRec, bool is2012, bool isEmu)
{
  doEmu = isEmu;
  if(doRec) {
    doRecoil = doRec;
    cout << "doing recoil corrections" << endl;
    if(!isEmu) {
      if(doRecoil == 1) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_zmm53X_2012_njet.root");
      if(doRecoil == 2) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_wjets53X_20pv_njet.root");
      if(doRecoil == 3) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_higgs53X_20pv_njet.root");
      corrector->addMCFile      ("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_zmm53XRR_2012_njet.root");
      corrector->addDataFile    ("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_datamm53XRR_2012_njet.root");
    } else {
      if(is2012) {
	if(doRecoil == 1) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_higgsem53X_20pv_njet.root");
	if(doRecoil == 2) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_zmm53X_20pv_njet.root");
	corrector->addMCFile      ("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_zmm53XRR_2012_njet.root");
	corrector->addDataFile    ("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_datamm53XRR_2012_njet.root");
      } else {
	if(doRecoil == 1) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_higgsem42X_20pv_njet.root");
	if(doRecoil == 2) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_zmm42X_20pv_njet.root");
	corrector->addMCFile      ("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_zmm42X_20pv_njet.root");
	corrector->addDataFile    ("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_datamm42X_20pv_njet.root");
      }
    }
  }
}
void Output::setupOutput(TString name) {
  fOutputFile = new TFile( name, "RECREATE" );
  fEventTree = new TTree("Events","Events");
  fEventTree->Branch("run",&fRun,"fRun/I");
  fEventTree->Branch("lumi",&fLumi,"fLumi/I");
  fEventTree->Branch("evt",&fEvt,"fEvt/I");
  fEventTree->Branch("pvx",&fpvx,"fpvx/f");
  fEventTree->Branch("pvy",&fpvy,"fpvy/f");
  fEventTree->Branch("pvz",&fpvz,"fpvz/f");
  fEventTree->Branch("npv",&fNPV,"fNPV/I");
  fEventTree->Branch("npu",&fNPU,"fNPU/I");
  fEventTree->Branch("rho",&fRho,"fRho/F");
  fEventTree->Branch("mcweight",&fMCWeight,"fMCWeight/F");
  fEventTree->Branch("puweight",&fPUWeight,"fPUWeight/F");
  fEventTree->Branch("effweight",&fEffWeight,"fEffWeight/F");
  fEventTree->Branch("weight",&fWeight,"fWeight/F");
  fEventTree->Branch("embweight",&fEmbWeight,"fEmbWeight/F");
  fEventTree->Branch("embkinweight",&fEmbKinWeight,"fEmbKinWeight/F");
  fEventTree->Branch("embgenweight",&fEmbGenWeight,"fEmbGenWeight/F");
  fEventTree->Branch("embspinnerweight",&fEmbSpinnerWeight,"fEmbSpinnerWeight/F");
  fEventTree->Branch("embmueffeweight",&fEmbMuEff,"fEmbMuEff/F");
  fEventTree->Branch("embmuradweight",&fEmbMuRad,"fEmbMuRad/F");
  fEventTree->Branch("embkinmassweight",&fEmbKinMassWeight,"fEmbKinMassWeight/F");
  fEventTree->Branch("emblid",&fEmbIdScale,"fEmbIdScale/F");
  fEventTree->Branch("genmass",&fMass,"fMass/F");
  fEventTree->Branch("genpt",&fgenpt,"fgenpt/F");
  fEventTree->Branch("genphi",&fgenphi,"fgenphi/F");
  fEventTree->Branch("mvis",&fMVis,"fMVis/F");
  fEventTree->Branch("pt_1",&fPt1,"fPt1/F");
  fEventTree->Branch("phi_1",&fPhi1,"fPhi1/F");
  fEventTree->Branch("eta_1",&fEta1,"fEta1/F");
  fEventTree->Branch("m_1",&fM1,"fM1/F");
  fEventTree->Branch("q_1",&fq1,"fq1/I");
  fEventTree->Branch("iso_1",&fIso1,"fIso1/F");
  fEventTree->Branch("d0_1",&fD01,"fD01/F");
  fEventTree->Branch("dZ_1",&fDZ1,"fDZ1/F");
  fEventTree->Branch("passiso_1",&fPassIso1,"fPassIso1/O");
  fEventTree->Branch("mt_1",&fMt1,"fMt1/F");
  fEventTree->Branch("mtMVA_1",&fMVAMt1,"fMVAMt1/F");
  fEventTree->Branch("ngamma_1",&fngamma1,"fngamma1/I");
  fEventTree->Branch("nprong_1",&fnprong1,"fnprong1/I");
  fEventTree->Branch("antiele_1",&fantiele1,"fantiele1/O");
  fEventTree->Branch("pt_2",&fPt2,"fPt2/F");
  fEventTree->Branch("phi_2",&fPhi2,"fPhi2/F");
  fEventTree->Branch("eta_2",&fEta2,"fEta2/F");
  fEventTree->Branch("m_2",&fM2,"fM2/F");
  fEventTree->Branch("q_2",&fq2,"fq2/I");
  fEventTree->Branch("iso_2",&fIso2,"fIso2/F");
  fEventTree->Branch("d0_2",&fD02,"fD02/F");
  fEventTree->Branch("dZ_2",&fDZ2,"fDZ2/F");
  fEventTree->Branch("passiso_2",&fPassIso2,"fPassIso2/O");
  fEventTree->Branch("mt_2",&fMt2,"fMt2/F");
  fEventTree->Branch("mtMVA_2",&fMVAMt2,"fMVAMt2/F");
  fEventTree->Branch("ngamma_2",&fngamma2,"fngamma2/I");
  fEventTree->Branch("nprong_2",&fnprong2,"fnprong2/I");
  fEventTree->Branch("antiele_2",&fantiele2,"fantiele2/O");
  fEventTree->Branch("drll",&fdrll,"fdrll/F");
  fEventTree->Branch("met",&fMet,"fMet/F");
  fEventTree->Branch("metphi",&fMetPhi,"fMetPhi/F");
  fEventTree->Branch("mvamet",&fMVAMet,"fMVAMet/F");
  fEventTree->Branch("mvametphi",&fMVAMetPhi,"fMVAMetPhi/F");
  fEventTree->Branch("pzetavis",&fPZetaVis,"fPZetaVis/F");
  fEventTree->Branch("pzetamiss",&fPZetaMiss,"fPZetaMiss/F");
  fEventTree->Branch("pzetamvamiss",&fPZetaMVAMiss,"fPZetaMVAMiss/F");
  fEventTree->Branch("metcov00",&fMetCov00,"fMetCov00/F");
  fEventTree->Branch("metcov01",&fMetCov01,"fMetCov01/F");
  fEventTree->Branch("metcov10",&fMetCov10,"fMetCov10/F");
  fEventTree->Branch("metcov11",&fMetCov11,"fMetCov11/F");
  fEventTree->Branch("mvacov00",&fMVACov00,"fMVACov00/F");
  fEventTree->Branch("mvacov01",&fMVACov01,"fMVACov01/F");
  fEventTree->Branch("mvacov10",&fMVACov10,"fMVACov10/F");
  fEventTree->Branch("mvacov11",&fMVACov11,"fMVACov11/F");
  fEventTree->Branch("jpt_1",&fJPt1,"fJPt1/F");
  fEventTree->Branch("jeta_1",&fJEta1,"fJEta1/F");
  fEventTree->Branch("jphi_1",&fJPhi1,"fJPhi1/F");
  fEventTree->Branch("jm_1",&fJM1,"fJM1/F");
  fEventTree->Branch("jptunc_1",&fJPtUnc1,"fJPtUnc1/F");
  fEventTree->Branch("jmva_1",&fJMVA1,"fJMVA1/F");
  fEventTree->Branch("jcsv_1",&fJcsv1,"fJcsv1/F");
  fEventTree->Branch("jpass_1",&fJPass1,"fJPass1/O");
  fEventTree->Branch("jpt_2",&fJPt2,"fJPt2/F");
  fEventTree->Branch("jeta_2",&fJEta2,"fJEta2/F");
  fEventTree->Branch("jphi_2",&fJPhi2,"fJPhi2/F");
  fEventTree->Branch("jm_2",&fJM2,"fJM2/F");
  fEventTree->Branch("jptunc_2",&fJPtUnc2,"fJPtUnc2/F");
  fEventTree->Branch("jmva_2",&fJMVA2,"fJMVA2/F");
  fEventTree->Branch("jcsv_2",&fJcsv2,"fJcsv2/F");
  fEventTree->Branch("jpass_2",&fJPass2,"fJPass2/O");
  fEventTree->Branch("bpt",&fBTagPt1,"fBtagPt1/F");
  fEventTree->Branch("beta",&fBTagEta1,"fBTagEta1/F");
  fEventTree->Branch("bphi",&fBTagPhi1,"fBTagPhi1/F");
  fEventTree->Branch("bm",&fBTagM1,"fBTagM1/F");
  fEventTree->Branch("bcsv",&fbcsv1,"fbcsv1/F");
  fEventTree->Branch("bpt_2",&fBTagPt2,"fBtagPt2/F");
  fEventTree->Branch("beta_2",&fBTagEta2,"fBTagEta2/F");
  fEventTree->Branch("bphi_2",&fBTagPhi2,"fBTagPhi2/F");
  fEventTree->Branch("bm_2",&fBTagM2,"fBTagM2/F");
  fEventTree->Branch("bcsv_2",&fbcsv2,"fbcsv2/F");
  fEventTree->Branch("mjj",&fMJJ,"fMJJ/F");
  fEventTree->Branch("jdeta",&fJDEta,"fJDEta/F");
  fEventTree->Branch("jdphi",&fJDPhi,"fJDPhi/F");
  fEventTree->Branch("dijetpt",&fDiJetPt,"fDiJetPt/F");
  fEventTree->Branch("dijetphi",&fDiJetPhi,"fDiJetPhi/F");
  fEventTree->Branch("hdijetphi",&fHDJetPhi,"fHDJetPhi/F");
  fEventTree->Branch("visjeteta",&fVisJetEta,"fVisJetEta/F");
  fEventTree->Branch("ptvis",&fPtVis,"fPtVis/F");
  fEventTree->Branch("pth",&fPtH,"fPtH/F");
  fEventTree->Branch("pthmva",&fPtHMVA,"fPtHMVA/F");
  fEventTree->Branch("nbtag",&fNBTag,"fNBTag/I");
  fEventTree->Branch("njets",&fNJetsClean,"fNJetsClean/I");
  fEventTree->Branch("njetsRaw",&fNJets,"fNJets/I");
  fEventTree->Branch("genlpt_1",&fGenPt1,"fJGenPt1/F");
  fEventTree->Branch("genlphi_1",&fGenPhi1,"fGenPhi1/F");
  fEventTree->Branch("genleta_1",&fGenEta1,"fGenEta1/F");
  fEventTree->Branch("genlid_1",&fGenId1 ,"fGenId1/I");
  fEventTree->Branch("genlpt_2",&fGenPt2,"fJGenPt2/F");
  fEventTree->Branch("genlphi_2",&fGenPhi2,"fGenPhi2/F");
  fEventTree->Branch("genleta_2",&fGenEta2,"fGenEta2/F");
  fEventTree->Branch("genlid_2",&fGenId2 ,"fGenId2/I");
  fEventTree->Branch("genvislpt_1",&fVisGenPt1,"fVisGenPt1/F");
  fEventTree->Branch("genvislphi_1",&fVisGenPhi1,"fVisGenPhi1/F");
  fEventTree->Branch("genvisleta_1",&fVisGenEta1,"fVisGenEta1/F");
  fEventTree->Branch("genvislid_1",&fVisGenId1 ,"fVisGenId1/I");
  fEventTree->Branch("genvislpt_2",&fVisGenPt2,"fJVisGenPt2/F");
  fEventTree->Branch("genvislphi_2",&fVisGenPhi2,"fVisGenPhi2/F");
  fEventTree->Branch("genvisleta_2",&fVisGenEta2,"fVisGenEta2/F");
  fEventTree->Branch("genvislid_2",&fVisGenId2 ,"fVisGenId2/I");
  fEventTree->Branch("genmatch",&fGenMatch ,"fGenMatch/I");
  fEventTree->Branch("genlepmatch",&fLepMatchAny ,"fLepMatchAny/I");
  fEventTree->Branch("genmatchlep",&fGenMatchLep ,"fGenMatchLep/I");
  fEventTree->Branch("vbfmva",&fMVA,"fMVA/F");
  fEventTree->Branch("njetingap",&fNJetInGap,"fNJetInGap/I");
  fEventTree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_1",&byCombinedIsolationDeltaBetaCorrRaw3Hits1,"byCombinedIsolationDeltaBetaCorrRaw3Hits1/F");
  fEventTree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2",&byCombinedIsolationDeltaBetaCorrRaw3Hits2,"byCombinedIsolationDeltaBetaCorrRaw3Hits2/F");
  fEventTree->Branch("againstElectronMVA3raw_1",&againstElectronMVA3raw1,"againstElectronMVA3raw1/F");
  fEventTree->Branch("againstElectronMVA3raw_2",&againstElectronMVA3raw2,"againstElectronMVA3raw2/F");
  fEventTree->Branch("byIsolationMVA2raw_1",&byIsolationMVA2raw1,"byIsolationMVA2raw1/F");
  fEventTree->Branch("byIsolationMVA2raw_2",&byIsolationMVA2raw2,"byIsolationMVA2raw2/F");
  fEventTree->Branch("againstMuonLoose2_1",&againstMuonLoose21,"againstMuonLoose21/F");
  fEventTree->Branch("againstMuonLoose2_2",&againstMuonLoose22,"againstMuonLoose22/F");
  fEventTree->Branch("againstMuonMedium2_1",&againstMuonMedium21,"againstMuonMedium21/F");
  fEventTree->Branch("againstMuonMedium2_2",&againstMuonMedium22,"againstMuonMedium22/F");
  fEventTree->Branch("againstMuonTight2_1",&againstMuonTight21,"againstMuonTight21/F");
  fEventTree->Branch("againstMuonTight2_2",&againstMuonTight22,"againstMuonTight22/F");
  fEventTree->Branch("passAntiEleMVA3_1",&passAntiEleMVA31,"passAntiEleMVA31/F");
  fEventTree->Branch("passAntiEleMVA3_2",&passAntiEleMVA32,"passAntiEleMVA32/F");
  fEventTree->Branch("passAntiEleNewWPMVA3_1",&passAntiEleNewWPMVA31,"passAntiEleNewWPMVA31/F");
  fEventTree->Branch("passAntiEleNewWPMVA3_2",&passAntiEleNewWPMVA32,"passAntiEleNewWPMVA32/F");
  fEventTree->Branch("Npartons",&fpartons,"fpartons/F");
  fEventTree->Branch("antitightele_1",&antitightele1,"antitightele1/O");
  fEventTree->Branch("antitightele_2",&antitightele2,"antitightele2/O");
  fEventTree->Branch("npt20jets",&npt20jets);
  fEventTree->Branch("btagArray",&btagArray);
  fEventTree->Branch("jptArray",&jptArray);
  fEventTree->Branch("jetaArray",&jetaArray);
}

void Output::save()
{
  fOutputFile->Write(); 
  delete fEventTree;
  delete fOutputFile;
  if(doRecoil) delete corrector;
}

void Output::fillMuon(const mithep::TMuon *muon, bool first, double iso, bool passiso)
{
  if(first) {
    fPt1 = muon->pt;
    fPhi1 = muon->phi;
    fEta1 = muon->eta;
    fM1  = 105.658369e-3;
    fq1 = muon->q;
    fIso1 = iso;
    fD01 = muon->d0;
    fDZ1 = muon->dz;
    fPassIso1 = passiso;
  } else {
    fPt2 = muon->pt;
    fPhi2 = muon->phi;
    fEta2 = muon->eta;
    fM2  = 105.658369e-3;
    fq2 = muon->q;
    fIso2 = iso;
    fD02 = muon->d0;
    fDZ2 = muon->dz;
    fPassIso2 = passiso;
  }
}

void Output::fillElectron(const mithep::TElectron *ele, bool first, double iso, bool passiso, unsigned int scale)
{
  if(first)
    {
      fPt1 = ele->pt;
      if(scale==1) { fPt1 *= 0.99; }
      else if(scale==2) { fPt1 *= 1.01; }
      fPhi1 = ele->phi;
      fEta1 = ele->eta;
      fM1  = 0.000511;
      fq1 = ele->q;
      fIso1 = iso;
      fD01 = ele->d0;
      fDZ1 = ele->dz;
      fPassIso1 = passiso;
    }
  else
    {
      fPt2 = ele->pt;
      fPhi2 = ele->phi;
      fEta2 = ele->eta;
      fM2  = 0.000511;
      fq2 = ele->q;
      fIso2 = iso/ele->pt;
      fD02 = ele->d0;
      fDZ2 = ele->dz;
      fPassIso2 = passiso;
    }
}

void Output::fillTau(const mithep::TPFTau *tau, bool first, bool passiso, bool tauescorr, float antielemva3nwp)
{
  if(first)
    {
      fPhi1 = tau->phi;
      fEta1 = tau->eta;
      fq1 = tau->q;
      fIso1 = tau->ringIso;
      fD01 =  tau->leadChargedHadronPFCand.d0;
      fDZ1 = tau->leadChargedHadronPFCand.dz;
      fPassIso1 = passiso;
      fngamma1  = tau->nSignalPFGammaCands;
      fnprong1  = tau->nSignalPFChargedHadrCands;
      fantiele1   = (tau->hpsDiscriminators & mithep::TPFTau::kMVAEle);
      antitightele1 = (tau->hpsDiscriminators & mithep::TPFTau::kTightEle);
      passAntiEleMVA31 = tau->passAntiEleMVA3;
      passAntiEleNewWPMVA31 = antielemva3nwp;
      if(tauescorr)
	{
	  if(fnprong1==1 && fngamma1==0)
	    lshift1 = 0.0;
	  else if(fnprong1==1 && fngamma1>0)
	    lshift1 = prong1(tau->pt);
	  else
	    lshift1 = prong3(tau->pt);
	}     
      fPt1=(1.0+lshift1)*tau->pt;
      fM1 = (1.0+lshift1)*tau->m;
      byCombinedIsolationDeltaBetaCorrRaw3Hits1 = tau->rawIso3Hits;
      againstElectronMVA3raw1 = tau->antiEleMVA3;
      byIsolationMVA2raw1  = tau->ringIso2;
      if(tau->passAntiMu2>0) againstMuonLoose21=1;
      if(tau->passAntiMu2>1.5) againstMuonMedium21=1;
      if(tau->passAntiMu2>2.5) againstMuonTight21=1;
    }
  else
    {
      fPhi2 = tau->phi;
      fEta2 = tau->eta;
      fq2 = tau->q;
      fIso2 = tau->ringIso;
      fD02 =  tau->leadChargedHadronPFCand.d0;
      fDZ2 = tau->leadChargedHadronPFCand.dz;
      fPassIso2 = passiso;
      fngamma2  = tau->nSignalPFGammaCands;
      fnprong2  = tau->nSignalPFChargedHadrCands;
      fantiele2   = (tau->hpsDiscriminators & mithep::TPFTau::kMVAEle);	
      antitightele2 = (tau->hpsDiscriminators & mithep::TPFTau::kTightEle);
      passAntiEleMVA32 = tau->passAntiEleMVA3;
      passAntiEleNewWPMVA32 = antielemva3nwp;
      if(tauescorr)
	{
	  if(fnprong2==1 && fngamma2==0)
	    lshift2 =0.0;
	  else if(fnprong2==1 && fngamma2>0)
	    lshift2 = prong1(tau->pt);
	  else
	    lshift2 = prong3(tau->pt);
	}     
      fPt2=(1.0+lshift2)*tau->pt;
      fM2 = (1.0+lshift2)*tau->m;
      byCombinedIsolationDeltaBetaCorrRaw3Hits2 = tau->rawIso3Hits;
      againstElectronMVA3raw2 = tau->antiEleMVA3;
      byIsolationMVA2raw2  = tau->ringIso2;
      if(tau->passAntiMu2>0) againstMuonLoose22=1;
      //againstMuonLoose22 = tau->hpsDiscriminators & mithep::TPFTau::kTightMu;
      if(tau->passAntiMu2>1.5) againstMuonMedium22=1;
      if(tau->passAntiMu2>2.5) againstMuonTight22=1;
    }
}

void Output::fillCov(mithep::TSVfit *svfit)
{
  fMVAMet    = svfit->mvaMET;
  fMVAMetPhi = svfit->mvaMETphi;
  fMetCov00 = svfit->cov_00;
  fMetCov01 = svfit->cov_01;
  fMetCov10 = svfit->cov_10;
  fMetCov11 = svfit->cov_11;
  fMVACov00 = svfit->mvacov_00;
  fMVACov01 = svfit->mvacov_01;
  fMVACov10 = svfit->mvacov_10;
  fMVACov11 = svfit->mvacov_11;
}

void Output::fillJets(const mithep::TJet *jet1,const mithep::TJet *jet2, const mithep::TJet *bjet1, const mithep::TJet *bjet2, int njetsclean, int njets, int nbjets, int npt20,int nCentralJets)
{
  fJPt1		 = (jet1) ? jet1->pt  : 0;
  fJEta1	 = (jet1) ? jet1->eta : 0;
  fJPhi1	 = (jet1) ? jet1->phi : 0;
  fJM1	         = (jet1) ? jet1->mass : 0;
  fJPtUnc1	 = (jet1) ? jet1->unc : 0;
  fJMVA1	 = (jet1) ? jet1->mva : 0;
  fJcsv1         = (jet1) ? jet1->csv : 0;
  fJPass1	 = (jet1) ? jet1->id  : 0;
  fJPt2		 = (jet2) ? jet2->pt  : 0;
  fJEta2	 = (jet2) ? jet2->eta : 0;
  fJPhi2	 = (jet2) ? jet2->phi : 0;
  fJM2	         = (jet2) ? jet2->mass : 0;
  fJPtUnc2	 = (jet2) ? jet2->unc : 0;
  fJMVA2	 = (jet2) ? jet2->mva : 0;
  fJcsv2         = (jet2) ? jet2->csv : 0;
  fJPass2	 = (jet2) ? jet2->id  : 0;
  fBTagPt1	 = (bjet1) ? bjet1->pt  : 0;
  fBTagEta1	 = (bjet1) ? bjet1->eta : 0;
  fBTagPhi1	 = (bjet1) ? bjet1->phi : 0; 
  fBTagM1	 = (bjet1) ? bjet1->mass : 0; 
  fbcsv1         = (bjet1) ? bjet1->csv : 0;
  fBTagPt2	 = (bjet2) ? bjet2->pt  : 0;
  fBTagEta2	 = (bjet2) ? bjet2->eta : 0;
  fBTagPhi2	 = (bjet2) ? bjet2->phi : 0; 
  fBTagM2	 = (bjet2) ? bjet2->mass : 0; 
  fbcsv2         = (bjet2) ? bjet2->csv : 0;
  fNJetsClean    = njetsclean;
  fNJets         = njets;
  fNBTag         = nbjets;
  fNJetInGap	 = (njets>1) ? nCentralJets : 0;
  npt20jets      = npt20;
}

void Output::fillGen(mithep::TGenInfo *gen)
{

  if(gen->pt_1_a > gen->pt_2_a) {
    fGenPt1  = gen->pt_1_a;
    fGenEta1 = gen->eta_1_a;
    fGenPhi1 = gen->phi_1_a;
    fGenId1  = gen->id_1_a; 
    fGenPt2  = gen->pt_2_a;
    fGenEta2 = gen->eta_2_a;
    fGenPhi2 = gen->phi_2_a;
    fGenId2  = gen->id_2_a; 
  } else {
    fGenPt2  = gen->pt_1_a;
    fGenEta2 = gen->eta_1_a;
    fGenPhi2 = gen->phi_1_a;
    fGenId2  = gen->id_1_a; 
    fGenPt1  = gen->pt_2_a;
    fGenEta1 = gen->eta_2_a;
    fGenPhi1 = gen->phi_2_a;
    fGenId1  = gen->id_2_a; 
  } 
  if(gen->pt_1_b > gen->pt_2_b) {
    fVisGenPt1  = gen->pt_1_b;
    fVisGenEta1 = gen->eta_1_b;
    fVisGenPhi1 = gen->phi_1_b;
    fVisGenId1  = gen->id_1_b; 
    fVisGenPt2  = gen->pt_2_b;
    fVisGenEta2 = gen->eta_2_b;
    fVisGenPhi2 = gen->phi_2_b;
    fVisGenId2  = gen->id_2_b; 
  } else {
    fVisGenPt2  = gen->pt_1_b;
    fVisGenEta2 = gen->eta_1_b;
    fVisGenPhi2 = gen->phi_1_b;
    fVisGenId2  = gen->id_1_b; 
    fVisGenPt1  = gen->pt_2_b;
    fVisGenEta1 = gen->eta_2_b;
    fVisGenPhi1 = gen->phi_2_b;
    fVisGenId1  = gen->id_2_b; 
  } 
  int lNTauMatch = 0;
  int lNLepMatch = 0;
  fLepMatchAny = 0;
  fGenMatch = 0;
  fGenMatchLep = 0;
  
  if(toolbox::deltaR(fEta1,fPhi1,fGenEta1,fGenPhi1) < 0.3) 
    {
      if(fabs(fGenId1) > 14)  
	{
	  lNTauMatch++;
	  if(fabs(fGenId1) == 17 || fabs(fGenId1) == 18)
	    fGenMatchLep++;
	} 
      else
	{lNLepMatch++;}
    }
  if(toolbox::deltaR(fEta1,fPhi1,fGenEta2,fGenPhi2) < 0.3) 
    {
      if(fabs(fGenId2) > 14)  
	{
	  lNTauMatch++;
	  if(fabs(fGenId2) == 17 || fabs(fGenId2) == 18)
	    fGenMatchLep++;
	} 
      else 
	{lNLepMatch++;}
    }
  if(toolbox::deltaR(fEta2,fPhi2,fGenEta1,fGenPhi1) < 0.3) 
    {
      if(fabs(fGenId1) > 14)  
	{
	  lNTauMatch++;
	  if(fabs(fGenId1) == 17 || fabs(fGenId1) == 18)
	    fGenMatchLep++;
	} 
      else {
	lNLepMatch++;
      if(fGenPt1>8)
	fLepMatchAny++;
      }
    }
  if(toolbox::deltaR(fEta2,fPhi2,fGenEta2,fGenPhi2) < 0.3) 
    {
      if(fabs(fGenId2) > 14)  
	{
	  lNTauMatch++;
	  if(fabs(fGenId2) == 17 || fabs(fGenId2) == 18)
	    fGenMatchLep++;
	} 
      else {
	lNLepMatch++;
	if(fGenPt2>8)
	  fLepMatchAny++;
      }
    }
  if(lNLepMatch == 2                   ) fGenMatch = 1;
  if(lNTauMatch == 1 && lNLepMatch == 1) fGenMatch = 2;
  if(lNTauMatch == 2                   ) fGenMatch = 3;
  if(lNLepMatch == 1 && lNTauMatch == 0) fGenMatch = 4;
  if(lNLepMatch == 0 && lNTauMatch == 1) fGenMatch = 5;
 
  fMass  = gen->vmass_a;
  fgenpt = gen->vpt_a;
  fgenphi = gen->vphi_a;
  fpartons = gen->npartons;
}

void Output::fillEvent(mithep::TEventInfo *info, HttMVA *vbfmva, int npv, double scalecorr)
{
  fRun = info->runNum;
  fLumi = info->lumiSec;
  fEvt = info->evtNum;
  fNPV = npv;
  fpvx = info->pvx;
  fpvy = info->pvy;
  fpvz = info->pvz;
  fNPU = info->nPUTrue;
  fRho = info->rho;

  fMet = info->pfMET + scalecorr;
  fMVAMet += scalecorr;
  fMetPhi = info->pfMETphi;
  //fMVAMet     = info->mvaMET;
  //MVAMetPhi  = info->mvaMETphi;
  //fMVACov00 = info->mvaCov00;
  //fMVACov01 = info->mvaCov01;
  //fMVACov10 = info->mvaCov10;
  //fMVACov11 = info->mvaCov11;
  TLorentzVector ldvecr1; ldvecr1.SetPtEtaPhiM(fPt1*lshift1/(lshift1+1),fEta1,fPhi1,fM1*lshift1/(1.0+lshift1));
  TLorentzVector ldvecr2; ldvecr2.SetPtEtaPhiM(fPt2*lshift2/(lshift2+1),fEta2,fPhi2,fM2*lshift2/(1.0+lshift2));
  TLorentzVector lVMet;  lVMet.SetPtEtaPhiM(fMet,0,fMetPhi,0);
  TLorentzVector lVMVAMet;  lVMVAMet.SetPtEtaPhiM(fMVAMet,0,fMVAMetPhi,0);
  lVMet  -= ldvecr1;
  lVMet  -= ldvecr2;
  lVMVAMet  -= ldvecr1;
  lVMVAMet  -= ldvecr2;
  fMet=lVMet.Pt();
  fMetPhi=lVMet.Phi();
  fMVAMet=lVMVAMet.Pt();
  fMVAMetPhi=lVMVAMet.Phi();

  TLorentzVector lep1, lep2, dilep;
  lep1.SetPtEtaPhiM(fPt1, fEta1, fPhi1,fM1);
  lep2.SetPtEtaPhiM(fPt2, fEta2, fPhi2,fM2);
  dilep = lep1+lep2;
 
  TLorentzVector jv1, jv2,dijet;
  jv1.SetPtEtaPhiM(fJPt1,fJEta1,fJPhi1,fJM1);
  jv2.SetPtEtaPhiM(fJPt2,fJEta2,fJPhi2,fJM2);
  dijet = jv1+jv2;
  
  // recoil corrections
  double pU1         = 0;  //--
  double pU2         = 0;  //--
  double lMVAMet     = fMVAMet;
  double lMVAMetPhi  = fMVAMetPhi;
  
  if(corrector && doEmu) corrector->CorrectType1(lMVAMet, lMVAMetPhi, fgenpt, fgenphi, dilep.Pt(), dilep.Phi(), pU1, pU2, 0, 0, fNJets);
  else if(corrector && doRecoil != 2) corrector->CorrectType1(lMVAMet, lMVAMetPhi, fgenpt, fgenphi, dilep.Pt(), dilep.Phi(), pU1, pU2, 0, 0, fNJetsClean);
  else if(corrector && doRecoil == 2) corrector->CorrectType1(lMVAMet, lMVAMetPhi, fgenpt, fgenphi, lep1 .Pt(), lep1 .Phi(), pU1, pU2, 0, 0, fNJets);
  fMVAMet            = lMVAMet;
  fMVAMetPhi         = lMVAMetPhi;

  // calculate projection variables
  TVector3 l1,l2,metv,mvametv;
  l1.SetPtEtaPhi(fPt1,0,fPhi1);
  l2.SetPtEtaPhi(fPt2,0,fPhi2);
  metv.SetPtEtaPhi(fMet,0,fMetPhi);
  mvametv.SetPtEtaPhi(fMVAMet,0,fMVAMetPhi);
  TVector3 bisector(l1.Unit() + l2.Unit());
  bisector = bisector.Unit();
  
  // ditau system
  TLorentzVector met4v,mvamet4v, higgs,higgsmvamet;
  met4v.SetPtEtaPhiM(fMet,0,fMetPhi,0);
  mvamet4v.SetPtEtaPhiM(fMVAMet,0,fMVAMetPhi,0);
  higgs = dilep+met4v;
  higgsmvamet = dilep+mvamet4v;

  fMVis = dilep.M();
  fdrll = toolbox::deltaR(fEta1,fPhi1,fEta2,fPhi2);
  fMt1	= sqrt(2.0*(lep1.Pt()*fMet*(1.0-cos(toolbox::deltaPhi(lep1.Phi(),fMetPhi)))));
  fMVAMt1 = sqrt(2.0*(lep1.Pt()*fMVAMet*(1.0-cos(toolbox::deltaPhi(lep1.Phi(),fMVAMetPhi)))));
  fMt2	 = sqrt(2.0*(lep2.Pt()*fMet*(1.0-cos(toolbox::deltaPhi(lep2.Phi(),fMetPhi)))));
  fMVAMt2 = sqrt(2.0*(lep2.Pt()*fMVAMet*(1.0-cos(toolbox::deltaPhi(lep2.Phi(),fMVAMetPhi)))));
  fPZetaVis	 = (l1+l2).Dot(bisector);
  fPZetaMiss	 = metv.Dot(bisector);
  fPZetaMVAMiss  = mvametv.Dot(bisector);
  fMJJ		 = (fNJetsClean>1) ? dijet.M() : 0;
  fJDEta	 = (fNJetsClean>1) ? fabs(fJEta1 - fJEta2) : 0;
  fJDPhi	 = (fNJetsClean>1) ? toolbox::deltaPhi(fJPhi1,fJPhi2) : 0;
  fDiJetPt	 = (fNJetsClean>1) ? dijet.Pt()  : 0;
  fDiJetPhi	 = (fNJetsClean>1) ? dijet.Phi() : 0;
  fHDJetPhi	 = (fNJetsClean>1) ? toolbox::deltaPhi(dijet.Phi(),higgs.Phi()) : 0;
  fVisJetEta	 = (fNJetsClean>1) ? TMath::Min(fabs(dilep.Eta()-fJEta1),fabs(dilep.Eta()-fJEta2)) : 0;
  fPtVis	 = dilep.Pt();
  fPtH		 = higgs.Pt();
  fPtHMVA        = higgsmvamet.Pt();
  fMVA		 = (vbfmva)  ? vbfmva->MVAValue(fMJJ, fJDEta, fJDPhi, fDiJetPt, fPtH, fHDJetPhi, fVisJetEta, fPtVis) : 0;
  fEventTree->Fill();
} 
//double Output::prong1(double pt)
//{
// return 0.025+0.001*TMath::Min(TMath::Max(pt-45,0.0),10.0);
//}
//double Output::prong3(double pt)
//{
//  return 0.012+0.001*TMath::Min(TMath::Max(pt-32,0.0),18.0);
//}

double Output::prong1(double pt)
{
  return 0.012;
}
double Output::prong3(double pt)
{
  return 0.012;
}


        

