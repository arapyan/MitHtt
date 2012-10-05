std::vector<std::string> fMap;
std::vector<std::string> fXAxis;
std::vector<std::string> fYAxis;
std::vector<int>         fNBins;
std::vector<double>      fXMin;
std::vector<double>      fXMax;

void loadfMap() {
  fMap.push_back("metcov00"); fNBins.push_back(25);  fXMin.push_back(0);    fXMax.push_back(250);   fXAxis.push_back("cov_{00}  (GeV/c^{2})"); fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("metcov10"); fNBins.push_back(25);  fXMin.push_back(0);    fXMax.push_back(250);   fXAxis.push_back("cov_{10}  (GeV/c^{2})"); fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("metcov01"); fNBins.push_back(25);  fXMin.push_back(0);    fXMax.push_back(250);   fXAxis.push_back("cov_{01}  (GeV/c^{2})"); fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("metcov11"); fNBins.push_back(25);  fXMin.push_back(0);    fXMax.push_back(23350); fXAxis.push_back("cov_{11}  (GeV/c^{2})"); fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("sumet");    fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(2000);  fXAxis.push_back("#Sigma E_{T}  (GeV/c)"); fYAxis.push_back("Events/20 GeV/c^{2}");
  fMap.push_back("sumet-pt_1-pt_2");       fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(2000);  fXAxis.push_back("#Sigma E_{T}  (GeV/c)"); fYAxis.push_back("Events/20 GeV/c^{2}");
  fMap.push_back("chsumet-pt_1-pt_2");     fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(2000);  fXAxis.push_back("#Sigma E_{T}  (GeV/c)"); fYAxis.push_back("Events/20 GeV/c^{2}");
  fMap.push_back("nopusumet");   fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(2000);  fXAxis.push_back("#Sigma E_{T}  (GeV/c)"); fYAxis.push_back("Events/20 GeV/c^{2}");
  fMap.push_back("pucsumet");   fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(2000);  fXAxis.push_back("#Sigma E_{T}  (GeV/c)"); fYAxis.push_back("Events/20 GeV/c^{2}");
  fMap.push_back("pusumet");   fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(2000);  fXAxis.push_back("#Sigma E_{T}  (GeV/c)"); fYAxis.push_back("Events/20 GeV/c^{2}");
  fMap.push_back("chsumet");   fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(2000);  fXAxis.push_back("#Sigma E_{T}  (GeV/c)"); fYAxis.push_back("Events/20 GeV/c^{2}");

  fMap.push_back("rho");      fNBins.push_back(250); fXMin.push_back(0);   fXMax.push_back(75);    fXAxis.push_back("#rho  (GeV)");           fYAxis.push_back("Events/GeV");
  fMap.push_back("TMath::Min(abs(tkmet-pt1_l),abs(met-pt1_l))");      fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(50); fXAxis.push_back("#Delta #slash{E_{T}} (true-reco)  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("TMath::Min(abs(nopumet-pt1_l),abs(met-pt1_l))");      fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(50); fXAxis.push_back("#Delta #slash{E_{T}} (true-reco)  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("TMath::Min(abs(tkmet-pt1_l),TMath::Min(abs(nopumet-pt1_l),abs(met-pt1_l)))");      fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(50); fXAxis.push_back("#Delta #slash{E_{T}} (true-reco)  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("abs(met-pt1_l)");      fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(50); fXAxis.push_back("|#Delta #slash{E_{T}}| (true-reco)  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("met-pt1_l");      fNBins.push_back(100); fXMin.push_back(-50);   fXMax.push_back(50); fXAxis.push_back("#Delta #slash{E_{T}} (true-reco)  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("tkmet-pt1_l");      fNBins.push_back(100); fXMin.push_back(-50);   fXMax.push_back(50); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("nopumet-pt1_l");      fNBins.push_back(100); fXMin.push_back(-50);   fXMax.push_back(50); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("TMath::Min(tkmet,met)");      fNBins.push_back(1000); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("TMath::Min(nopumet,met)");      fNBins.push_back(1000); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("TMath::Min(TMath::Min(nopumet,met),tkmet)");      fNBins.push_back(1000); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");

  fMap.push_back("(rpt_z*sqrt(cos(rphi_z)*cos(uphi_mva) + sin(rphi_z)*sin(uphi_mva))-ux_mva*u)/sqrt(u1_cov)");         fNBins.push_back(50); fXMin.push_back(-5.);   fXMax.push_back(5.); fXAxis.push_back("Pull #parallel"); fYAxis.push_back("Pull");
  fMap.push_back("(rpt_z*(cos(rphi_z)*sin(uphi_mva) - sin(rphi_z)*cos(uphi_mva)))/sqrt(u2_cov)");         fNBins.push_back(50); fXMin.push_back(-5.);   fXMax.push_back(5.); fXAxis.push_back("Pull #perpl"); fYAxis.push_back("Pull #perp");
  fMap.push_back("pt_z");        fNBins.push_back(30); fXMin.push_back(0);   fXMax.push_back(240); fXAxis.push_back("pT_{W}^{Gen}"); fYAxis.push_back("Events/8 GeV");
  fMap.push_back("rpt_z");        fNBins.push_back(250); fXMin.push_back(0);   fXMax.push_back(250); fXAxis.push_back("#slash{E_{T}}  (GeV)"); fYAxis.push_back("Events/6 GeV");
  fMap.push_back("met");         fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV)"); fYAxis.push_back("Events/6 GeV");
  fMap.push_back("mvamet");      fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV)"); fYAxis.push_back("Events/6 GeV");
  fMap.push_back("tkmet");       fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events/2 GeV");
  fMap.push_back("nopumet");     fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events/2 GeV");
  fMap.push_back("met_mva");     fNBins.push_back(25); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events/2 GeV");
  fMap.push_back("mety_mva");     fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("metx_mva");     fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("u1_cov");       fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#sigma_{u_{1}} (GeV)");   fYAxis.push_back("Events GeV");
  fMap.push_back("u2_cov");       fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#sigma_{u_{1}} (GeV)");   fYAxis.push_back("Events GeV");
  fMap.push_back("sqrt(unc_cov)");       fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(30); fXAxis.push_back("#sigma_{cov} (GeV)");   fYAxis.push_back("Events/GeV");
  fMap.push_back("metx_mva/sqrt(u1_cov)");   fNBins.push_back(50); fXMin.push_back(0);   fXMax.push_back(30); fXAxis.push_back("#slash{E_{T}}/#sigma_{u_{1}} (GeV)");   fYAxis.push_back("Events/GeV");
  fMap.push_back("pumet");       fNBins.push_back(25); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("pucmet");      fNBins.push_back(25); fXMin.push_back(0);   fXMax.push_back(150); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("nopumet-pt2_l");      fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(100); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("pumet-pt2_l");      fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(100); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("pucmet-pt2_l");      fNBins.push_back(100); fXMin.push_back(0);   fXMax.push_back(100); fXAxis.push_back("#slash{E_{T}}  (GeV/c)"); fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("m_sv");       fNBins.push_back(30); fXMin.push_back(0);    fXMax.push_back(600); fXAxis.push_back("m_{#tau#tau}  (GeV/c^{2})");  fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("m_svipf");    fNBins.push_back(30); fXMin.push_back(0);    fXMax.push_back(600); fXAxis.push_back("m_{#tau#tau}  (GeV/c^{2})");  fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("m_sv2");      fNBins.push_back(30); fXMin.push_back(0);    fXMax.push_back(300); fXAxis.push_back("m_{#tau#tau}  (GeV/c^{2})");  fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("nopumass");   fNBins.push_back(60); fXMin.push_back(0);    fXMax.push_back(300); fXAxis.push_back("m_{#tau#tau}  (GeV/c^{2})");  fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("pumass");     fNBins.push_back(60); fXMin.push_back(0);    fXMax.push_back(300); fXAxis.push_back("m_{#tau#tau}  (GeV/c^{2})");  fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("m_vis");    fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(200); fXAxis.push_back("m_{vis}  (GeV/c^{2})");  fYAxis.push_back("Events/7.5 GeV/c^{2}");
  fMap.push_back("mjj");      fNBins.push_back(25); fXMin.push_back(0);    fXMax.push_back(1500); fXAxis.push_back("m_{jj}  (GeV/c^{2})");  fYAxis.push_back("Events/60 GeV/c^{2}");
  fMap.push_back("pt_vis/u1");      fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(2);   fXAxis.push_back("m_{#tau} (GeV/c^{2})");  fYAxis.push_back("Events/0.05 GeV/c^{2}");
  fMap.push_back("jpfiga06_1");      fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(10.);   fXAxis.push_back("jet 04 charged isolation");                    fYAxis.push_back("Events/0.25");
  fMap.push_back("jnehaf_1");      fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(1.);   fXAxis.push_back("jet area");                         fYAxis.push_back("Events/0.02");
  fMap.push_back("lLfr_1");      fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(2.0);   fXAxis.push_back("p_{T} fraction of leading cand");    fYAxis.push_back("Events/0.02");
  fMap.push_back("jnehaf_1+jneemf_1");      fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(3.);   fXAxis.push_back("jet em/had");                    fYAxis.push_back("Events/0.02");
  fMap.push_back("jneemf_1");      fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(1.);   fXAxis.push_back("jet area");                    fYAxis.push_back("Events/0.02");
  fMap.push_back("genfrac_1");      fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(1.6);   fXAxis.push_back("gen fraction from hard scatter");                    fYAxis.push_back("Events/0.04");
  fMap.push_back("higgswidth"); fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(0.6);   fXAxis.push_back("Relative");                    fYAxis.push_back("Events/0.04");
  fMap.push_back("genfrac_2");      fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(3.);   fXAxis.push_back("jet area");                    fYAxis.push_back("Events/0.02");
  fMap.push_back("genfrac_3");      fNBins.push_back(40); fXMin.push_back(0.);   fXMax.push_back(3.);   fXAxis.push_back("jet area");                    fYAxis.push_back("Events/0.02");
  fMap.push_back("jarea_1");      fNBins.push_back(40); fXMin.push_back(0.6);   fXMax.push_back(1.2);   fXAxis.push_back("jet area");                    fYAxis.push_back("Events/0.02");
  fMap.push_back("jbeta_1");      fNBins.push_back(40); fXMin.push_back(-10);   fXMax.push_back(10);   fXAxis.push_back("jet #beta");                    fYAxis.push_back("Events/0.02");
  fMap.push_back("u_mva");      fNBins.push_back(20); fXMin.push_back(0);   fXMax.push_back(3);   fXAxis.push_back("MET Scale corrections");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("jmva_1");      fNBins.push_back(50); fXMin.push_back(-1);   fXMax.push_back(1);   fXAxis.push_back("leading jet mva");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("JetID");      fNBins.push_back(40); fXMin.push_back(-1);   fXMax.push_back(1);   fXAxis.push_back("leading jet mva");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("jmva_2");      fNBins.push_back(40); fXMin.push_back(-1);   fXMax.push_back(1);   fXAxis.push_back("mva");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("mva");      fNBins.push_back(40); fXMin.push_back(-1);   fXMax.push_back(1);   fXAxis.push_back("mva");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("ggmva");    fNBins.push_back(40); fXMin.push_back(-1);   fXMax.push_back(1);   fXAxis.push_back("mva");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("vttmva");   fNBins.push_back(40); fXMin.push_back(-1);   fXMax.push_back(1);   fXAxis.push_back("mva");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("vbfmva");   fNBins.push_back(40); fXMin.push_back(-1);   fXMax.push_back(1);   fXAxis.push_back("mva");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("cos(phi_1-phi_2)");   fNBins.push_back(40); fXMin.push_back(-1);   fXMax.push_back(1);   fXAxis.push_back("cos(#Delta #phi)");                    fYAxis.push_back("Events/0.05");
  fMap.push_back("m_2");      fNBins.push_back(40); fXMin.push_back(0);    fXMax.push_back(2);   fXAxis.push_back("m_{#tau} (GeV/c^{2})");   fYAxis.push_back("Events/0.1 GeV/c^{2}");
  fMap.push_back("pmet");     fNBins.push_back(40); fXMin.push_back(-50);  fXMax.push_back(50);  fXAxis.push_back("#slash{p_{#zeta}} (GeV)"); fYAxis.push_back("Events/2.5 GeV");
  fMap.push_back("pzeta");    fNBins.push_back(50); fXMin.push_back(-50);  fXMax.push_back(200); fXAxis.push_back("p_{#zeta} (GeV)");        fYAxis.push_back("Events/5 GeV");
  fMap.push_back("-pzeta");   fNBins.push_back(30); fXMin.push_back(-300);  fXMax.push_back(200); fXAxis.push_back("p_{#zeta} (GeV)");        fYAxis.push_back("Events/5 GeV");
  fMap.push_back("emfrac_2"); fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(1);   fXAxis.push_back("emfrac_{#tau} ");         fYAxis.push_back("Events/0.05 ");
  fMap.push_back("mt_1");     fNBins.push_back(25);  fXMin.push_back(0);    fXMax.push_back(150); fXAxis.push_back("m_{T}    (GeV/c^{2})");   fYAxis.push_back("Events/5 GeV/c^{2}");
  fMap.push_back("mvamt1");     fNBins.push_back(25);  fXMin.push_back(0);    fXMax.push_back(150); fXAxis.push_back("m_{T}    (GeV/c^{2})");   fYAxis.push_back("Events/5 GeV/c^{2}");
  fMap.push_back("npv");      fNBins.push_back(60); fXMin.push_back(-0.5); fXMax.push_back(59.5); fXAxis.push_back("# pv");                   fYAxis.push_back("Events");
  fMap.push_back("npu");      fNBins.push_back(48); fXMin.push_back(-0.5); fXMax.push_back(47.5); fXAxis.push_back("# pu");                   fYAxis.push_back("Events");
  fMap.push_back("npumean");      fNBins.push_back(6); fXMin.push_back(-0.5); fXMax.push_back(47.5); fXAxis.push_back("# pu mean");                   fYAxis.push_back("Events");
  fMap.push_back("npuplus");      fNBins.push_back(6); fXMin.push_back(-0.5); fXMax.push_back(47.5); fXAxis.push_back("# +pu");                   fYAxis.push_back("Events");
  fMap.push_back("pt_vis");   fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(100); fXAxis.push_back("p^{vis}_{T} (GeV/c)");    fYAxis.push_back("Events/5 GeV/c^{2}");
  fMap.push_back("pt_sys");   fNBins.push_back(25); fXMin.push_back(0);    fXMax.push_back(100); fXAxis.push_back("p^{sys}_{T} (GeV/c)");    fYAxis.push_back("Events/10 GeV/c");
  fMap.push_back("pt_1");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(60);  fXAxis.push_back("p_{T}    (GeV/c)");       fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("lnept_1");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(20);  fXAxis.push_back("leading neutral p_{T} (GeV/c)");       fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("lempt_1");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(20);  fXAxis.push_back("p_{T}    (GeV/c)");       fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("lchpt_1");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(20);  fXAxis.push_back("p_{T}    (GeV/c)");       fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("spt_1/lpt_1");    fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(3.);  fXAxis.push_back("p_{T}    (GeV/c)");       fYAxis.push_back("Events/10 GeV/c^{2}");
  fMap.push_back("jm_1");     fNBins.push_back(30); fXMin.push_back(0);    fXMax.push_back(30);  fXAxis.push_back("mass leading Jet (GeV/c^{2})");       fYAxis.push_back("Events GeV/c^{2}");
  fMap.push_back("jpt1");     fNBins.push_back(80); fXMin.push_back(0);    fXMax.push_back(200);  fXAxis.push_back("p_{T} leading Jet (GeV/c)");       fYAxis.push_back("Events/10 GeV/c");
  fMap.push_back("jsptraw_1");     fNBins.push_back(90); fXMin.push_back(10);    fXMax.push_back(100);  fXAxis.push_back("p_{T} leading Jet (GeV/c)");       fYAxis.push_back("Events/10 GeV/c");
  fMap.push_back("(jspt_1*(abs(jseta_1) > 2.5)+jsptraw_1*(abs(jseta_1) < 2.5))");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(100);  fXAxis.push_back("p_{T} leading Jet (GeV/c)");       fYAxis.push_back("Events/10 GeV/c");
  fMap.push_back("jspt_1");     fNBins.push_back(30); fXMin.push_back(10);    fXMax.push_back(100);  fXAxis.push_back("p_{T} leading Jet (GeV/c)");       fYAxis.push_back("Events/3 GeV/c");
  fMap.push_back("jspt_3");     fNBins.push_back(30); fXMin.push_back(10);    fXMax.push_back(100);  fXAxis.push_back("p_{T} leading Jet (GeV/c)");       fYAxis.push_back("Events/3 GeV/c");
  fMap.push_back("lchpt_1");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(200);  fXAxis.push_back("p_{T} leading charged cand (GeV/c)");       fYAxis.push_back("Events/10 GeV/c");
 fMap.push_back("jspt_2");     fNBins.push_back(80); fXMin.push_back(0);    fXMax.push_back(200);  fXAxis.push_back("p_{T} leading Jet (GeV/c)");       fYAxis.push_back("Events/10 GeV/c");
  fMap.push_back("mt_2");     fNBins.push_back(30); fXMin.push_back(0);    fXMax.push_back(150); fXAxis.push_back("m^{#tau}_{T} (GeV/c^{2})"); fYAxis.push_back("Events/5 GeV/c^{2}");
  fMap.push_back("pt_2");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(60); fXAxis.push_back("p^{#tau}_{T} (GeV/c)");    fYAxis.push_back("Events/5 GeV/c^{2}");
  fMap.push_back("pt_2*emfrac_2");  fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(50); fXAxis.push_back("p^{em-#tau}_{T} (GeV/c)");     fYAxis.push_back("Events/2.5 GeV/c");
  fMap.push_back("jd0_1");          fNBins.push_back(25); fXMin.push_back(-0.10); fXMax.push_back(0.10); fXAxis.push_back("d^{jet}_{0}(cm)");        fYAxis.push_back("Events/0.002");
  fMap.push_back("jdZ_1");          fNBins.push_back(25); fXMin.push_back(-0.50); fXMax.push_back(0.50); fXAxis.push_back("d^{jet}_{Z}(cm)");        fYAxis.push_back("Events/0.002");
  fMap.push_back("abs(d0_1)");     fNBins.push_back(25); fXMin.push_back(0.00); fXMax.push_back(0.05); fXAxis.push_back("d^{l}_{0}(cm)");        fYAxis.push_back("Events/0.002");
  fMap.push_back("abs(d0_2)");     fNBins.push_back(25); fXMin.push_back(0.00); fXMax.push_back(25.0); fXAxis.push_back("d^{#tau}_{0}(sig)");     fYAxis.push_back("Events/1");
  fMap.push_back("jceta");    fNBins.push_back(11); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta^{J}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("nopumeteta");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta^{J}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("jeta1");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta^{J}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("jseta_1");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta^{J}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("leta_1");      fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta leading cand");        fYAxis.push_back("Events/0.3");
  fMap.push_back("lemeta_1");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta leading em cand");        fYAxis.push_back("Events/0.3");
  fMap.push_back("lcheta_1");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta leading charged cand");   fYAxis.push_back("Events/0.3");
  fMap.push_back("lneeta_1");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta leading neut cand");   fYAxis.push_back("Events/0.3");
  fMap.push_back("seta_1");      fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta 2nd leading cand");   fYAxis.push_back("Events/0.3");
  fMap.push_back("jseta_2");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta^{J}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("jseta_3");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta^{J}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("min(abs(jeta1),abs(jeta2))");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta^{J}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("jeta2");    fNBins.push_back(31); fXMin.push_back(-5.); fXMax.push_back(5.); fXAxis.push_back("#eta^{J}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("eta_1");    fNBins.push_back(25); fXMin.push_back(-2.5); fXMax.push_back(2.5); fXAxis.push_back("#eta^{l}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("y_z");      fNBins.push_back(25); fXMin.push_back(-2.5); fXMax.push_back(2.5); fXAxis.push_back("#eta^{W}_{Gen}");     fYAxis.push_back("Events/0.3");
  fMap.push_back("eta_2");    fNBins.push_back(30); fXMin.push_back(-2.5); fXMax.push_back(2.5); fXAxis.push_back("#eta^{#tau}");     fYAxis.push_back("Events/0.5");
  fMap.push_back("iso_1");    fNBins.push_back(30); fXMin.push_back(-0.2); fXMax.push_back(0.5); fXAxis.push_back("#eta^{#tau}");     fYAxis.push_back("Events/0.5");
  fMap.push_back("njet");     fNBins.push_back(10); fXMin.push_back(-0.5); fXMax.push_back(9.5); fXAxis.push_back("njet");             fYAxis.push_back("Events/1.0");
  fMap.push_back("nbtag");    fNBins.push_back(6); fXMin.push_back(-0.5); fXMax.push_back(5.5); fXAxis.push_back("njet");             fYAxis.push_back("Events/1.0");
  fMap.push_back("TMath::Min(abs(phi_sys-jphi1),6.28-abs(phi_sys-jphi1))");    fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("#phi^{#ell}");     fYAxis.push_back("Events/0.32");
  fMap.push_back("TMath::Min(abs(jphi1-jphi2),6.28-abs(jphi1-jphi2))"); fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("#phi^{#ell}");     fYAxis.push_back("Events/0.32");
  fMap.push_back("sqrt(TMath::Min(abs(phi_1-phi_2),2.*TMath::Pi()-abs(phi_1-phi_2))**2+abs(eta_1-eta_2)**2)"); fNBins.push_back(20); fXMin.push_back(0.);fXMax.push_back(6.); fXAxis.push_back("#Deta_{R}");     fYAxis.push_back("Events/0.32");
  fMap.push_back("jphi1");    fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("#phi^{#ell}");     fYAxis.push_back("Events/0.32");
  fMap.push_back("phi_1");    fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("#phi^{#ell}");     fYAxis.push_back("Events/0.32");
  fMap.push_back("jsphi_1");  fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("leading jet #phi");     fYAxis.push_back("Events/0.32");
  fMap.push_back("lchphi_1");  fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("leading charged cand #phi");     fYAxis.push_back("Events/0.32");
  fMap.push_back("lemphi_1");  fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("leading em cand #phi");     fYAxis.push_back("Events/0.32");
  fMap.push_back("lnephi_1");  fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("leading neut cand #phi");    fYAxis.push_back("Events/0.32");
  fMap.push_back("lphi_1");    fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("leading #phi");              fYAxis.push_back("Events/0.32");
  fMap.push_back("sphi_1");    fNBins.push_back(20); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("2nd leading #phi");          fYAxis.push_back("Events/0.32");
  fMap.push_back("phi_2");     fNBins.push_back(200); fXMin.push_back(-3.15);fXMax.push_back(3.15); fXAxis.push_back("#phi^{#tau}");           fYAxis.push_back("Events/0.32");
 fMap.push_back("-u1_mva/pt_z");        fNBins.push_back(60); fXMin.push_back(-2); fXMax.push_back(2);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
 fMap.push_back("-u1x_mva/pt_vis");        fNBins.push_back(60); fXMin.push_back(-2); fXMax.push_back(2);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
 fMap.push_back("-u1y_mva/pt_z");        fNBins.push_back(60); fXMin.push_back(-2); fXMax.push_back(2);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
 fMap.push_back("-u1/pt_vis");        fNBins.push_back(60); fXMin.push_back(-2); fXMax.push_back(2);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
 fMap.push_back("-tku1/pt_vis");        fNBins.push_back(60); fXMin.push_back(-2); fXMax.push_back(2);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
 fMap.push_back("-nopuu1/pt_vis");        fNBins.push_back(60); fXMin.push_back(-2); fXMax.push_back(2);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");

  fMap.push_back("u1x_mva+pt_vis");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u1y_mva+pt_vis");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u1x_mva+rpt_z");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u1_mva+pt_vis");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u1+pt_vis");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("#Delta u_{1}-p_{T}^{Z}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("tku1+pt_vis");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("nopuu1+pt_vis");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u_mva*nopuu-rpt_z");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u-rpt_z");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("#Delta u-p_{T}^{Z}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("tku-rpt_z");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("nopuu-rpt_z");        fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u1+rpt_z+pucu1*0.5");     fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("tku1+rpt_z+pucu1*0.5");   fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("nopuu1+rpt_z+pucu1*1.5"); fNBins.push_back(60); fXMin.push_back(-75); fXMax.push_back(75);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u1x_mva");    fNBins.push_back(125); fXMin.push_back(-75); fXMax.push_back(50);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u1_mva");    fNBins.push_back(125); fXMin.push_back(-75); fXMax.push_back(50);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u1");        fNBins.push_back(125); fXMin.push_back(-75); fXMax.push_back(50);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("nopuu2*u_mva");    fNBins.push_back(100); fXMin.push_back( -50); fXMax.push_back(50);  fXAxis.push_back("u_{2}    (GeV)");          fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u2x_mva");    fNBins.push_back(100); fXMin.push_back( -50); fXMax.push_back(50);  fXAxis.push_back("u_{2}    (GeV)");          fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u2_mva");    fNBins.push_back(100); fXMin.push_back( -50); fXMax.push_back(50);  fXAxis.push_back("u_{2}    (GeV)");          fYAxis.push_back("Events/1 GeV");
  fMap.push_back("u2");        fNBins.push_back(100); fXMin.push_back( -50); fXMax.push_back(50);  fXAxis.push_back("u_{2}    (GeV)");          fYAxis.push_back("Events/1 GeV");
  fMap.push_back("tku1");      fNBins.push_back(125); fXMin.push_back(-75); fXMax.push_back(50);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("tku2");      fNBins.push_back(100); fXMin.push_back( -50); fXMax.push_back(50);  fXAxis.push_back("u_{2}    (GeV)");          fYAxis.push_back("Events/1 GeV");
  fMap.push_back("nopuu1");    fNBins.push_back(125); fXMin.push_back(-75); fXMax.push_back(50);  fXAxis.push_back("u_{1}    (GeV)");           fYAxis.push_back("Events/1 GeV");
  fMap.push_back("nopuu2");    fNBins.push_back(100); fXMin.push_back( -50); fXMax.push_back(50);  fXAxis.push_back("u_{2}    (GeV)");          fYAxis.push_back("Events/1 GeV");
  fMap.push_back("TMath::Min(abs(phi_1-metphi),2.*TMath::Pi()-abs(phi_1-metphi))");
  fNBins.push_back(10); fXMin.push_back(0); fXMax.push_back(TMath::Pi()); fXAxis.push_back("#Delta #phi_{#ell#slash{E_{T}}}");       fYAxis.push_back("Events/#pi/10");
  fMap.push_back("npart_1");   fNBins.push_back(50); fXMin.push_back(-0.5); fXMax.push_back(49.5);  fXAxis.push_back("jet n particles");       fYAxis.push_back("Events");
  fMap.push_back("lpt_1");   fNBins.push_back(50); fXMin.push_back(-0.5); fXMax.push_back(49.5);  fXAxis.push_back("leading cand p_{T} (GeV)");       fYAxis.push_back("Events");
  fMap.push_back("spt_1");   fNBins.push_back(50); fXMin.push_back(-0.5); fXMax.push_back(49.5);  fXAxis.push_back("2nd leading cand p_{T} (GeV)");       fYAxis.push_back("Events");
  fMap.push_back("drlc_1");   fNBins.push_back(50); fXMin.push_back(-0); fXMax.push_back(0.5);  fXAxis.push_back("#Delta R leading to jet");       fYAxis.push_back("Events");
  fMap.push_back("drls_1");   fNBins.push_back(50); fXMin.push_back(-0); fXMax.push_back(0.5);  fXAxis.push_back("#Delta R leading to 2nd");       fYAxis.push_back("Events");
  fMap.push_back("drm_1");   fNBins.push_back(30); fXMin.push_back(-0); fXMax.push_back(0.5);  fXAxis.push_back("<#Delta R of deposits>");       fYAxis.push_back("Events");
  fMap.push_back("drm_3");   fNBins.push_back(30); fXMin.push_back(-0); fXMax.push_back(0.5);  fXAxis.push_back("<#Delta R of deposits>");       fYAxis.push_back("Events");
  fMap.push_back("drch_1");     fNBins.push_back(50); fXMin.push_back(-0); fXMax.push_back(0.5);  fXAxis.push_back("<#Delta R of charges particles>");       fYAxis.push_back("Events");
  fMap.push_back("jdrmin_1");   fNBins.push_back(210); fXMin.push_back(-0); fXMax.push_back(10.5);  fXAxis.push_back(">8 GeV gen jet #Delta R Min");      fYAxis.push_back("Events");
  fMap.push_back("jmcf_1");   fNBins.push_back(50); fXMin.push_back(-24.5); fXMax.push_back(25.5);  fXAxis.push_back(">8 GeV gen jet #Delta R Min");      fYAxis.push_back("Events");
  fMap.push_back("drem_1");     fNBins.push_back(50); fXMin.push_back(-0); fXMax.push_back(0.35);  fXAxis.push_back("<#Delta R of em particles>");       fYAxis.push_back("Events");
  fMap.push_back("drmne_1");     fNBins.push_back(50); fXMin.push_back(-0); fXMax.push_back(0.35);  fXAxis.push_back("<#Delta R of neut particles>");       fYAxis.push_back("Events");
  ////
  fMap.push_back("lpt_1/jspt_1");       fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(1);  fXAxis.push_back("leading p_{T}/Jet p_{T} ");           fYAxis.push_back("Events/0.05");
  fMap.push_back("spt_1/jspt_1");       fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(1);  fXAxis.push_back("2nd leading p_{T}/Jet p_{T} ");       fYAxis.push_back("Events/0.05");
  fMap.push_back("lnept_1/jspt_1");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(1);  fXAxis.push_back("leading neutral p_{T}/Jet p_{T} ");   fYAxis.push_back("Events/0.05");
  fMap.push_back("lempt_1/jspt_1");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(1);  fXAxis.push_back("leading em p_{T}/Jet p_{T}");         fYAxis.push_back("Events/0.05");
  fMap.push_back("lchpt_1/jspt_1");     fNBins.push_back(20); fXMin.push_back(0);    fXMax.push_back(1);  fXAxis.push_back("leading ch p_{T}/Jet p_{T}");         fYAxis.push_back("Events/0.05");
}

int getId(std::string iStr) { 
  for(unsigned int i0 = 0; i0 < fMap.size(); i0++) { 
    if(iStr == fMap[i0]) return i0;
  }
  return 0;
}
const char*       getXAxis(std::string iStr) { return fXAxis[getId(iStr)].c_str();}
const char*       getYAxis(std::string iStr) { return fYAxis[getId(iStr)].c_str();}
int         getNBins(std::string iStr) { return fNBins[getId(iStr)];}
double      getXMin (std::string iStr) { return fXMin [getId(iStr)];}
double      getXMax (std::string iStr) { return fXMax [getId(iStr)];}
