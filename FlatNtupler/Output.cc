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
  fJPtUnc1(0),
  fJMVA1(0),
  fJcsv1(0),
  fJPass1(0),
  fJPt2(0),
  fJEta2(0),
  fJPhi2(0),
  fJPtUnc2(0),
  fJMVA2(0),
  fJcsv2(0),
  fJPass2(0),
  fBTagPt(0),
  fBTagEta(0),
  fBTagPhi(0),
  fbcsv(0),
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
  fNBTag(0),
  fNJets(0),
  fGenPt1(0),
  fGenPhi1(0),
  fGenEta1(0),
  fGenId1(0),
  fGenPt2(0),
  fGenPhi2(0),
  fGenEta2(0),
  fGenId2(0),
  doRecoil(0),
  npt20jets(0),
  fOutputFile(0),
  fEventTree(0)
{
  btagArray.Set(50);
  jptArray.Set(50);
  jetaArray.Set(50);

  if(doRecoil > 0)
     {
       cout << "doing recoil corrections" << endl;
       if(doRecoil == 1) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_zmm53X_20pv_njet.root");
       if(doRecoil == 2) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_wjets53X_20pv_njet.root");
       if(doRecoil == 3) corrector = new RecoilCorrector("$CMSSW_BASE/src/MitHtt/Utils/recoilfits/recoilfit_higgs53X_20pv_njet.root");
       corrector->addMCFile      ("files/recoilfit_zmm53X_20pv_njet.root");
       corrector->addDataFile    ("files/recoilfit_datamm53X_20pv_njet.root");
     }
  setupOutput(name);
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
  fEventTree->Branch("genmass",&fMass,"fMass/F");
  fEventTree->Branch("genpt",&fgenpt,"fgenpt/F");
  fEventTree->Branch("genphi",&fgenphi,"fgenphi/F");
  fEventTree->Branch("m_vis",&fMVis,"fMVis/F");
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
  fEventTree->Branch("ngamma_1",&fngamma1,"fngamma1/F");
  fEventTree->Branch("nprong_1",&fnprong1,"fnprong1/F");
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
  fEventTree->Branch("ngamma_2",&fngamma2,"fngamma2/F");
  fEventTree->Branch("nprong_2",&fnprong2,"fnprong2/F");
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
  fEventTree->Branch("jptunc_1",&fJPtUnc1,"fJPtUnc1/F");
  fEventTree->Branch("jmva_1",&fJMVA1,"fJMVA1/F");
  fEventTree->Branch("jcsv_1",&fJcsv1,"fJcsv1/F");
  fEventTree->Branch("jpass_1",&fJPass1,"fJPass1/O");
  fEventTree->Branch("jpt_2",&fJPt2,"fJPt2/F");
  fEventTree->Branch("jeta_2",&fJEta2,"fJEta2/F");
  fEventTree->Branch("jphi_2",&fJPhi2,"fJPhi2/F");
  fEventTree->Branch("jptunc_2",&fJPtUnc2,"fJPtUnc2/F");
  fEventTree->Branch("jmva_2",&fJMVA2,"fJMVA2/F");
  fEventTree->Branch("jcsv_2",&fJcsv2,"fJcsv2/F");
  fEventTree->Branch("jpass_2",&fJPass2,"fJPass2/O");
  fEventTree->Branch("bpt",&fBTagPt,"fBtagPt/F");
  fEventTree->Branch("beta",&fBTagEta,"fBTagEta/F");
  fEventTree->Branch("bphi",&fBTagPhi,"fBTagPhi/F");
  fEventTree->Branch("bcsv",&fbcsv,"fbcsv/F");
  fEventTree->Branch("mjj",&fMJJ,"fMJJ/F");
  fEventTree->Branch("jdeta",&fJDEta,"fJDEta/F");
  fEventTree->Branch("jdphi",&fJDPhi,"fJDPhi/F");
  fEventTree->Branch("dijetpt",&fDiJetPt,"fDiJetPt/F");
  fEventTree->Branch("dijetphi",&fDiJetPhi,"fDiJetPhi/F");
  fEventTree->Branch("hdijetphi",&fHDJetPhi,"fHDJetPhi/F");
  fEventTree->Branch("visjeteta",&fVisJetEta,"fVisJetEta/F");
  fEventTree->Branch("ptvis",&fPtVis,"fPtVis/F");
  fEventTree->Branch("pth",&fPtH,"fPtH/F");
  fEventTree->Branch("nbtag",&fNBTag,"fNBTag/I");
  fEventTree->Branch("njets",&fNJets,"fNJets/I");
  fEventTree->Branch("genlpt_1",&fGenPt1,"fJGenPt1/F");
  fEventTree->Branch("genlphi_1",&fGenPhi1,"fGenPhi1/F");
  fEventTree->Branch("genleta_1",&fGenEta1,"fGenEta1/F");
  fEventTree->Branch("genlid_1",&fGenId1 ,"fGenId1/I");
  fEventTree->Branch("genlpt_2",&fGenPt2,"fJGenPt2/F");
  fEventTree->Branch("genlphi_2",&fGenPhi2,"fGenPhi2/F");
  fEventTree->Branch("genleta_2",&fGenEta2,"fGenEta2/F");
  fEventTree->Branch("genlid_2",&fGenId2 ,"fGenId2/I");
  fEventTree->Branch("genmatch",&fGenMatch ,"fGenMatch/I");
  fEventTree->Branch("vbfmva",&fMVA,"fMVA/F");
  fEventTree->Branch("njetingap",&fNJetInGap,"fNJetInGap/F");
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

void Output::fillMuon(const mithep::TMuon *muon, double iso, bool passiso)
{
  fPt1 = muon->pt;
  fPhi1 = muon->phi;
  fEta1 = muon->eta;
  fM1  = 105.658369e-3;
  fq1 = muon->q;
  fIso1 = iso/muon->pt;
  fD01 = muon->d0;
  fDZ1 = muon->dz;
  fPassIso1 = passiso;
}

void Output::fillElectron(const mithep::TElectron *ele, bool first, double iso, bool passiso)
{
  if(first)
    {
      fPt1 = ele->pt;
      fPhi1 = ele->phi;
      fEta1 = ele->eta;
      fM1  = 0.000511;
      fq1 = ele->q;
      fIso1 = iso/ele->pt;
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

void Output::fillTau(const mithep::TPFTau *tau, bool first, bool passiso)
{
  if(first)
    {
      fPt1 = tau->pt;
      fPhi1 = tau->phi;
      fEta1 = tau->eta;
      fM1  = tau->m;
      fq1 = tau->q;
      fIso1 = tau->ringIso;
      fD01 =  tau->leadChargedHadronPFCand.d0;
      fDZ1 = tau->leadChargedHadronPFCand.dz;
      fPassIso1 = passiso;
      fngamma1  = tau->nSignalPFGammaCands;
      fnprong1  = tau->nSignalPFChargedHadrCands;
    }
  else
    {
      fPt2 = tau->pt;
      fPhi2 = tau->phi;
      fEta2 = tau->eta;
      fM2  = tau->m;
      fq2 = tau->q;
      fIso2 = tau->ringIso;
      fD02 =  tau->leadChargedHadronPFCand.d0;
      fDZ2 = tau->leadChargedHadronPFCand.dz;
      fPassIso2 = passiso;
      fngamma2  = tau->nSignalPFGammaCands;
      fnprong2  = tau->nSignalPFChargedHadrCands;
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

void Output::fillJets(const mithep::TJet *jet1,const mithep::TJet *jet2, const mithep::TJet *bjet, int njets, int nbjets, int npt20,int nCentralJets)
{
  fJPt1		 = (jet1) ? jet1->pt  : 0;
  fJEta1	 = (jet1) ? jet1->eta : 0;
  fJPhi1	 = (jet1) ? jet1->phi : 0;
  fJPtUnc1	 = (jet1) ? jet1->unc : 0;
  fJMVA1	 = (jet1) ? jet1->mva : 0;
  fJcsv1         = (jet1) ? jet1->csv : 0;
  fJPass1	 = (jet1) ? jet1->id  : 0;
  fJPt2		 = (jet2) ? jet2->pt  : 0;
  fJEta2	 = (jet2) ? jet2->eta : 0;
  fJPhi2	 = (jet2) ? jet2->phi : 0;
  fJPtUnc2	 = (jet2) ? jet2->unc : 0;
  fJMVA2	 = (jet2) ? jet2->mva : 0;
  fJcsv2         = (jet2) ? jet2->csv : 0;
  fJPass2	 = (jet2) ? jet2->id  : 0;
  fBTagPt	 = (bjet) ? bjet->pt  : 0;
  fBTagEta	 = (bjet) ? bjet->eta : 0;
  fBTagPhi	 = (bjet) ? bjet->phi : 0; 
  fbcsv         = (bjet) ? bjet->csv : 0;
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
  int lNTauMatch = 0;
  int lNLepMatch = 0;
  fGenMatch = 0;
  
  if(toolbox::deltaR(fEta1,fPhi1,fGenEta1,fGenPhi1) < 0.3) {if(fabs(fGenId1) > 14)  {lNTauMatch++;} else{lNLepMatch++;}}
  if(toolbox::deltaR(fEta1,fPhi1,fGenEta2,fGenPhi2) < 0.3) {if(fabs(fGenId2) > 14)  {lNTauMatch++;} else {lNLepMatch++;}}
  if(toolbox::deltaR(fEta2,fPhi2,fGenEta1,fGenPhi1) < 0.3) {if(fabs(fGenId1) > 14)  {lNTauMatch++;} else {lNLepMatch++;}}
  if(toolbox::deltaR(fEta2,fPhi2,fGenEta2,fGenPhi2) < 0.3) {if(fabs(fGenId2) > 14)  {lNTauMatch++;} else {lNLepMatch++;}}
  if(lNLepMatch == 2                   ) fGenMatch = 1;
  if(lNTauMatch == 1 && lNLepMatch == 1) fGenMatch = 2;
  if(lNTauMatch == 2                   ) fGenMatch = 3;
  if(lNLepMatch == 1 && lNTauMatch == 0) fGenMatch = 4;
  if(lNLepMatch == 0 && lNTauMatch == 1) fGenMatch = 5;
 
  fMass  = gen->vmass_a;
  fgenpt = gen->vpt_a;
  fgenphi = gen->vphi_a;
}

void Output::fillEvent(mithep::TEventInfo *info, HttMVA *vbfmva, int npv)
{
  fRun = info->runNum;
  fLumi = info->lumiSec;
  fEvt = info->evtNum;
  fNPV = npv;
  fpvx = info->pvx;
  fpvy = info->pvy;
  fpvz = info->pvz;
  fNPU = info->nPU;
  fRho = info->rho;

  fMet = info->pfMET;
  fMetPhi = info->pfMETphi;
  
  TLorentzVector lep1, lep2, dilep;
  lep1.SetPtEtaPhiM(fPt1, fEta1, fPhi1,fM1);
  lep2.SetPtEtaPhiM(fPt2, fEta2, fPhi2,fM2);
  dilep = lep1+lep2;
 
  TLorentzVector jv1, jv2, dijet;
  jv1.SetPtEtaPhiM(fJPt1,fJEta1,fJPhi1,fJPhi1);
  jv2.SetPtEtaPhiM(fJPt2,fJEta2,fJPhi2,fJPhi2);
  dijet = jv1+jv2;
  
  // recoil corrections
  double pU1         = 0;  //--
  double pU2         = 0;  //--
  double lMVAMet     = fMVAMet;
  double lMVAMetPhi  = fMVAMetPhi;
  if(corrector && doRecoil != 2) corrector->CorrectType2(lMVAMet, lMVAMetPhi, fgenpt, fgenphi, dilep.Pt(), dilep.Phi(), pU1, pU2, 0, 0, fNJets);
  if(corrector && doRecoil == 2) corrector->CorrectType2(lMVAMet, lMVAMetPhi, fgenpt, fgenphi, lep1 .Pt(), lep1 .Phi(), pU1, pU2, 0, 0, fNJets);
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
  TLorentzVector met4v, higgs;
  met4v.SetPtEtaPhiM(fMVAMet,0,fMVAMetPhi,0);
  higgs = dilep+met4v;

  fMVis = dilep.M();
  fdrll = toolbox::deltaR(fEta1,fPhi1,fEta2,fPhi2);
  fMt1	= sqrt(2.0*(lep1.Pt()*fMet*(1.0-cos(toolbox::deltaPhi(lep1.Phi(),fMetPhi)))));
  fMVAMt1 = sqrt(2.0*(lep1.Pt()*fMVAMet*(1.0-cos(toolbox::deltaPhi(lep1.Phi(),fMVAMetPhi)))));
  fMt2	 = sqrt(2.0*(lep2.Pt()*fMet*(1.0-cos(toolbox::deltaPhi(lep2.Phi(),fMetPhi)))));
  fMVAMt2 = sqrt(2.0*(lep2.Pt()*fMVAMet*(1.0-cos(toolbox::deltaPhi(lep2.Phi(),fMVAMetPhi)))));
  fPZetaVis	 = (l1+l2).Dot(bisector);
  fPZetaMiss	 = metv.Dot(bisector);
  fPZetaMVAMiss  = mvametv.Dot(bisector);
  fMJJ		 = (fNJets>1) ? dijet.M() : 0;
  fJDEta	 = (fNJets>1) ? fabs(fJEta1 - fJEta2) : 0;
  fJDPhi	 = (fNJets>1) ? toolbox::deltaPhi(fJPhi1,fJPhi2) : 0;
  fDiJetPt	 = (fNJets>1) ? dijet.Pt()  : 0;
  fDiJetPhi	 = (fNJets>1) ? dijet.Phi() : 0;
  fHDJetPhi	 = (fNJets>1) ? toolbox::deltaPhi(dijet.Phi(),higgs.Phi()) : 0;
  fVisJetEta	 = (fNJets>1) ? TMath::Min(fabs(dilep.Eta()-fJEta1),fabs(dilep.Eta()-fJEta2)) : 0;
  fPtVis	 = dilep.Pt();
  fPtH		 = higgs.Pt();
  fMVA		 = (vbfmva)  ? vbfmva->MVAValue(fMJJ, fJDEta, fJDPhi, fDiJetPt, fPtH, fHDJetPhi, fVisJetEta, fPtVis) : 0;
  fEventTree->Fill();
} 




        

