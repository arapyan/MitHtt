#include "Selector.hh"

//----------------------------------------------------------------------------------------
Selector::Selector(TString conf, TString outputDir, Double_t lumival):
  fJetUnc		(kNo),
  fBtagEff		(kNo),
  fMistag		(kNo),
  fMuonPt1Min		(20),
  fEMuMuonPt2Min	(10),
  fMuMuMuonPt2Min       (20),
  fElePt1Min		(20),
  fElePt2Min		(10),
  fJetPtMin		(30),
  fBJetPtMin		(20),
  fFailHist             (0),
  fEvtFail              (0),
  fLumi			(lumival),
  fDoNpuRwgt		(kTRUE),
  fCheckNpuHists	(kTRUE),
  fMakeNpuHists		(kTRUE),
  fDebug                (kFALSE),
  // fData		(),
  fRawMet		(0),
  fRawprojvar		(0),
  fNpuWgt		(1),
  kMaxPt20Jets		(35)
  // fBtagArray		(),
  
{
  fRandm = TRandom(234);

  fCorrector = RecoilCorrector("data/recoilfitZDat_800pb.root",
			       "data/recoilfitZMC.root", 0xDEADBEEF);
  
  fNtupDir = outputDir + TString("/ntuples");
  gSystem->mkdir(fNtupDir,kTRUE);

  fBtagArray.Set(kMaxPt20Jets);

  ParseConfig(conf);

  fInfo  = new mithep::TEventInfo();
  fGen   = new mithep::TGenInfo();
  fMuonArr     = new TClonesArray("mithep::TMuon");
  fElectronArr = new TClonesArray("mithep::TElectron");
  fJetArr      = new TClonesArray("mithep::TJet");
  fPvArr       = new TClonesArray("mithep::TVertex");
  fSvfitArr    = new TClonesArray("mithep::TSVFit");

  fDataNPVfname = "data/Pileup_2011_EPS_8_jul.root"; // officially produced predicted npu distribution

}
//----------------------------------------------------------------------------------------
Selector::~Selector()
{
  delete fInfo;
  delete fGen;
  delete fMuonArr;
  delete fElectronArr;
  delete fJetArr;
  delete fPvArr;
  delete fSvfitArr;
}
//----------------------------------------------------------------------------------------
void Selector::ParseConfig(TString conf)
{
  //
  // parse .conf file
  //

  ifstream ifs;
  ifs.open(conf.Data());
  assert(ifs.is_open());
  string line;
  Int_t state=0;
  vector<string> flagnamev; // names of the different corrections
  map<string,string> varmap;
  while(getline(ifs,line)) {
    //
    // skip comment lines (and define correction flags)
    //
    if(line[0]=='#') {
      // if line begins with "#flags" then read in correction flags
      if(line.find("#flags") != string::npos) {
	stringstream ss(line);
	string flagname;
	while(!ss.eof()) {
	  ss >> flagname;
	  if(flagname=="#flags") continue;
	  flagnamev.push_back(flagname);
	}
	continue;
      }
      // regular comment line: skip it
      else
	continue;
    }

    //
    // read in variables
    //
    if(line[0]=='^') {
      int ieq = line.find('=');
      string var = line.substr(1,ieq-1);
      string val = line.substr(ieq+1,line.size()-ieq);
      varmap[var] = val;
      continue;
    }

    //
    // switch from data -> mc mode
    //
    if(line[0]=='%') {
      state++; 
      continue; 
    }

    //
    // push back a new CSample
    //
    if(line[0]=='$') {
      if((TString(line).Contains("fake")) && (state>0)) continue; // fakes come from a separate macro: skip

      fSamplev.push_back(new CSample());
      stringstream ss(line);
      string chr,sname;
      Int_t color;
      ss >> chr >> sname >> color;
      string label = line.substr(line.find('@')+1);
      fSnamev.push_back(sname);
      fSamplev.back()->label = label;
      fSamplev.back()->color = color;
      continue;
    }

    //
    // initialize the new CSample
    //
    if(state==0) {  // define data sample
      stringstream ss(line);
      string fname;
      Int_t type;
      string json;
      ss >> fname >> type >> json;
      map<string,bool> *flags = new map<string,bool>;
      UInt_t iflag=0;
      bool flagvalue;
      while(!ss.eof()) {
	ss >> flagvalue;
	(*flags)[flagnamev[iflag]] = flagvalue;
	iflag++;
      }
      // replace the "variables" in file paths
      map<string,string>::iterator it;
      for(it=varmap.begin(); it!=varmap.end(); it++) {
      	ULong64_t pos = fname.find((*it).first);
      	if(pos!=string::npos)
      	  fname.replace(pos,(*it).first.size(),(*it).second);
	pos = json.find((*it).first);
      	if(pos!=string::npos)
      	  json.replace(pos,(*it).first.size(),(*it).second);
      }
      fSamplev.back()->fnamev.push_back(fname);
      fSamplev.back()->typev.push_back(type);
      fSamplev.back()->xsecv.push_back(0);
      fSamplev.back()->jsonv.push_back(json);
      fSamplev.back()->flags = *flags;

    } else if(state==1) {  // define MC samples
      stringstream ss(line);
      string fname;
      Double_t xsec;
      ss >> fname >> xsec;
      if(TString(fname).Contains("dummy",TString::kIgnoreCase)) continue; // skip the fakes
      map<string,bool> *flags = new map<string,bool>;
      UInt_t iflag=0;
      bool flagvalue;
      while(!ss.eof()) {
	ss >> flagvalue;
	(*flags)[flagnamev[iflag]] = flagvalue;
	iflag++;
      }
      // replace the "variables" in file paths
      map<string,string>::iterator it;
      for(it=varmap.begin(); it!=varmap.end(); it++) {
      	ULong64_t pos = fname.find((*it).first);
      	if(pos!=string::npos)
      	  fname.replace(pos,(*it).first.size(),(*it).second);
      }
      fSamplev.back()->fnamev.push_back(fname);
      fSamplev.back()->typev.push_back(0);
      fSamplev.back()->xsecv.push_back(xsec);
      fSamplev.back()->flags = *flags;
    }
  }
  ifs.close();
}

//????????????????????????????????????????????????????????????????????????????????????????
// // set up energy scale/smearing
// UInt_t escale = kCenter; // enums defined in EScale.hh
// EScale scaler("data/data-EnergyScale.root","data/mc-EnergyScale.root");
//????????????????????????????????????????????????????????????????????????????????????????

void Selector::SampleLoop()
{
  //
  // Access samples and fill histograms
  //  
  TFile *infile=0;
  fEventTree=0;  
  
  Bool_t hasData = (fSamplev[0]->fnamev.size()>0);

  //
  // loop over samples
  //
  for(fIsam=0; fIsam<fSamplev.size(); fIsam++) {
    if(fIsam==0 && !hasData) continue;
  
    CSample* samp = fSamplev[fIsam];

    Double_t nSelEvents[3]; for(Int_t i=0; i<3; i++) nSelEvents[i]=0; // events in this sample

    InitOutput();
    
    //
    // loop through files
    //
    cout <<  "processing " << fSnamev[fIsam] << ":" << endl;
    const UInt_t nfiles = samp->fnamev.size();
    for(fIfile=0; fIfile<nfiles; fIfile++) {
      printf("        %-65s",(samp->fnamev[fIfile]+"...").Data()); fflush(stdout);
      infile = new TFile(samp->fnamev[fIfile]); 
      assert(infile);

      // configuration and names
      TString sfname   = samp->fnamev[fIfile];
      fIsdata          = !(samp->typev[fIfile]==eMC);
      TString basename = sfname(sfname.Last('/')+1,sfname.Last('.') - sfname.Last('/') - 1);	// filename for npu hists
      fmcNPVfname = "npu/"+basename+"-npu.root";						// files from which to extract MC npu distribs.
      // get rid of skim label in filename -- we need the npu distrib. from the unskimmed sample
      RemoveSkimName(fmcNPVfname);

      // set up npu reweighting
      TH1D *hpu=0, *hpuRwgt=0;									// npu before/after reweighting
      vector<Double_t> puwgtv;
      InitNPU(sfname,hpu,hpuRwgt,puwgtv);

      // set up selecting with JSON file
      Bool_t hasJSON = kFALSE;
      mithep::RunLumiRangeMap rlrm;
      if((samp->jsonv.size()>0) && (samp->jsonv[fIfile].CompareTo("NONE")!=0)) { 
        hasJSON = kTRUE;
	ifstream jsonchk; jsonchk.open(samp->jsonv[fIfile].Data()); assert(jsonchk.is_open()); jsonchk.close();
        rlrm.AddJSONFile(samp->jsonv[fIfile].Data()); 
      }

      // Set up NNLO-NNLL k-factor reweighting (if necessary) [ not implemented ]
      TH1D *hKFactors = (samp->flags["doKFactors"]) ? kfInit("data/ggHWW_KFactors_PowhegToHQT.root",higgsmass(basename)) : 0;

      // Get the TTree
      fEventTree = (TTree*)infile->Get("Events"); assert(fEventTree);

      // Set branch address to structures that will store the info  
      fEventTree->SetBranchAddress("Info",     &fInfo);        TBranch *infoBr     = fEventTree->GetBranch("Info");
      fEventTree->SetBranchAddress("Muon",     &fMuonArr);     TBranch *muonBr     = fEventTree->GetBranch("Muon");
      fEventTree->SetBranchAddress("Electron", &fElectronArr); TBranch *electronBr = fEventTree->GetBranch("Electron");
      fEventTree->SetBranchAddress("PFJet",    &fJetArr);      TBranch *jetBr      = fEventTree->GetBranch("PFJet");      
      fEventTree->SetBranchAddress("PV",       &fPvArr);       TBranch *pvBr       = fEventTree->GetBranch("PV");
      fEventTree->SetBranchAddress("SVFit",    &fSvfitArr);    TBranch *svfitBr    = fEventTree->GetBranch("SVFit");
      TBranch *genBr=0;
      if(samp->flags["getGen"]) {
        fEventTree->SetBranchAddress("Gen",    &fGen); genBr = fEventTree->GetBranch("Gen");
      }

      // ~(xs/entries)
      Double_t weight = 1;
      if(!fIsdata) weight = GetWeight(sfname,samp);

      // counters
      Double_t nsel[3], nselvar[3]; for(Int_t i=0; i<3; i++)  { nsel[i] = nselvar[i] = 0; }	// events in this file
      Double_t nlowmass=0;									// low mass z events (below 50)

      // loop over events
      for(UInt_t ientry=0; ientry</*100000*/ fEventTree->GetEntries(); ientry++) {

	// if(!(ientry%10000)) if(fDebug) cout << "(" << ientry << "/" << fEventTree->GetEntries() << ")" << endl;

	fEvtFail = 0;

        infoBr->GetEntry(ientry);

	// fill pu histograms
	if(!fIsdata && (fDoNpuRwgt || fMakeNpuHists)) FillNPU(hpu,hpuRwgt,puwgtv);

	if(samp->flags["getGen"])  genBr->GetEntry(ientry);

	// madgraph z sample needs some selection for the xs to be correct
	if(samp->flags["ismadz"]) {
	  Int_t id = fabs(fGen->id_1_a);

	  if(samp->flags["MuMuSel"]) {
	    // mumu selection: skip non-mu events
	    if(id != EGenType::kMuon) { EvtFail(1); continue; }
	  } else {
	    // emu selection: skip non-tau events
	    Bool_t istau = id==EGenType::kTau         || id==EGenType::kTauMuon ||
	                   id==EGenType::kTauElectron || id==EGenType::kTauHadr;
	    if(!istau) { EvtFail(1); continue; }
	  }
	}
	  
	// not certified run? Skip to next event...
	mithep::RunLumiRangeMap::RunLumiPairType rl(fInfo->runNum, fInfo->lumiSec);
        if(hasJSON && !rlrm.HasRunLumi(rl)) {
	  EvtFail(2);
	  continue;
	}

        // No good primary vertex? Skip to next event...
        if(!fInfo->hasGoodPV) {
	  EvtFail(3);
	  continue;
	}

	//----------------------------------------------------------------------------------------
        // loop through muons
	//----------------------------------------------------------------------------------------
	vector<const mithep::TMuon*> goodMuonsv,looseMuonsv;
	fMuonArr->Clear();
	muonBr->GetEntry(ientry);
	MuonLoop(samp->flags["MuMuSel"],goodMuonsv,looseMuonsv);
	
	//----------------------------------------------------------------------------------------
        // loop through electrons
	//----------------------------------------------------------------------------------------
        vector<mithep::TElectron*> goodElectronsv;   
        fElectronArr->Clear();
        electronBr->GetEntry(ientry);
	ElectronLoop(goodElectronsv,looseMuonsv,goodMuonsv);

	//----------------------------------------------------------------------------------------
	// Require two leptons of the desired flavor
	//----------------------------------------------------------------------------------------
	TLorentzVector lep1, lep2, dilep; // lepton 4-vectors
	Int_t finalState=-1;	          // final state type
	const mithep::TMuon *mu=0;
	mithep::TElectron *ele=0;
	const mithep::TMuon *mu_2=0;
	if(samp->flags["MuMuSel"]) {
	  // apply mu mu selection
	  if(!PassMuMu(goodMuonsv,goodElectronsv,lep1,lep2,dilep,finalState,mu,mu_2)) {
	    continue;
	  }
	  assert(mu); assert(mu_2);
	}
	else {
	  // apply ele mu selection
	  if(!PassEmu(goodMuonsv,goodElectronsv,lep1,lep2,dilep,finalState,mu,ele)) {
	    continue;
	  }
	  assert(mu); assert(ele);
	}

	if(fDebug) cout << "------------> candidate! <-------------" << endl;

	//----------------------------------------------------------------------------------------
        // loop through jets      
	//----------------------------------------------------------------------------------------
        fJetArr->Clear();
        jetBr->GetEntry(ientry);
        UInt_t njets=0, nbjets=0;
        const mithep::TJet *jet1=0, *jet2=0, *bjet=0;
	fBtagArray.Reset(); fNpt20jets=0;
	JetLoop(lep1,lep2,njets,nbjets,jet1,jet2,bjet);

	//----------------------------------------------------------------------------------------
	// corrections
	//----------------------------------------------------------------------------------------

	// calculate and correct the projection variables
	Double_t projVis, projMet, met, metphi;
	Projections(lep1,lep2,dilep,samp->flags["doRecoil"],projVis,projMet,met,metphi);

        // get k-factor
        Double_t kf=1;
        if(samp->flags["doKFactors"]) kf = kfValue(fGen->vpt_a, hKFactors);

	// do vertex reweighting
	fNpuWgt = 1;
	if(!fIsdata && fDoNpuRwgt) fNpuWgt = (fInfo->nPU >= puwgtv.size()) ? 0 : puwgtv[fInfo->nPU];
	
	// multiply by trigger effic. in MC
	Double_t trigeff=1;
	if(samp->flags["doTrigEff"] && !samp->flags["MuMuSel"]) trigeff = GetTrigEff(mu,ele);

	//----------------------------------------------------------------------------------------
	// Get full event weight and add it up
	//----------------------------------------------------------------------------------------

	Double_t fullEvtWgt = weight*kf*fNpuWgt*trigeff; // event weight including corrections

	nsel[0]       += fullEvtWgt;    // passing events in this file
	nselvar[0]    += fullEvtWgt*fullEvtWgt;
        nSelEvents[0] += fullEvtWgt;    // passing events in whole sample 
	if(samp->flags["doRecoil"] && (fGen->vmass_a < 50)) nlowmass += fullEvtWgt;

	//----------------------------------------------------------------------------------------
	// get a few things that we only need to fill the tree
	//----------------------------------------------------------------------------------------

	// npv
	fPvArr->Clear(); pvBr->GetEntry(ientry);	

	// dijet mass
        TLorentzVector jv1, jv2, dijet;
	if(njets>1) {
	  jv1.SetPtEtaPhiM(jet1->pt,jet1->eta,jet1->phi,jet1->mass);
	  jv2.SetPtEtaPhiM(jet2->pt,jet2->eta,jet2->phi,jet2->mass);
          dijet = jv1+jv2;
        } 

	// SVFit info
        fSvfitArr->Clear();
        svfitBr->GetEntry(ientry);
        mithep::TSVFit *svfit = (mithep::TSVFit*)((*fSvfitArr)[0]);

	//----------------------------------------------------------------------------------------
	// fill the tree
	//----------------------------------------------------------------------------------------

        fData.runNum   = fInfo->runNum;
        fData.evtNum   = fInfo->evtNum;
        fData.lumiSec  = fInfo->lumiSec;
        fData.nPV      = fPvArr->GetEntriesFast();
        fData.njets    = njets;
        fData.nbjets   = nbjets;
        fData.met      = met;
	fData.metphi   = metphi;
        fData.mass     = dilep.M();
	fData.dphi     = toolbox::deltaPhi(lep1.Phi(),lep2.Phi());
	fData.mt       = sqrt( 2.0 * (dilep.Pt()) * met * (1.0-cos(toolbox::deltaPhi(dilep.Phi(),metphi))) );
	fData.pt       = dilep.Pt();
	fData.phi      = dilep.Phi();
	fData.pmet     = projMet;
	fData.pvis     = projVis;
        fData.lpt1     = lep1.Pt();
	fData.leta1    = lep1.Eta();
	fData.lphi1    = lep1.Phi();
        fData.lpt2     = lep2.Pt();
	fData.leta2    = lep2.Eta();
	fData.lphi2    = lep2.Phi();
        fData.jpt1     = (jet1) ? jet1->pt  : 0;
	fData.jeta1    = (jet1) ? jet1->eta : 0;
	fData.jphi1    = (jet1) ? jet1->phi : 0;
        fData.jpt2     = (jet2) ? jet2->pt  : 0;
	fData.jeta2    = (jet2) ? jet2->eta : 0;
	fData.jphi2    = (jet2) ? jet2->phi : 0;
        fData.bjpt     = (bjet) ? bjet->pt  : 0;
	fData.bjeta    = (bjet) ? bjet->eta : 0;
	fData.bjphi    = (bjet) ? bjet->phi : 0;
        fData.mjj      = (njets>1) ? dijet.M() : 0;
        fData.svfmass  = svfit->mass;
        fData.svflpt1  = (svfit->mass==0) ? 0 : svfit->daughter1.Pt();
        fData.svfleta1 = (svfit->mass==0) ? 0 : svfit->daughter1.Eta();
        fData.svflphi1 = (svfit->mass==0) ? 0 : svfit->daughter1.Phi();
        fData.svflpt2  = (svfit->mass==0) ? 0 : svfit->daughter2.Pt();
        fData.svfleta2 = (svfit->mass==0) ? 0 : svfit->daughter2.Eta();
        fData.svflphi2 = (svfit->mass==0) ? 0 : svfit->daughter2.Phi();
        fData.weight   = (fIsam==0) ? 1 : fullEvtWgt/fLumi;
        fData.state    = finalState;

	// reinitialize mass so we can use mass=0 as a criterion to see if it was filled for the next event
	svfit->mass = 0;

	fOuttree->Fill();
      } // event loop

      printf("%8.2f +/- %-8.2f\n",nsel[0],sqrt(nselvar[0]));
      if(nlowmass > 0) printf("           ---> selected events with z mass < 50:  %10.3f (out of %15.3f)\n",nlowmass,nsel[0]);

      if(!fIsdata && (fDoNpuRwgt || fMakeNpuHists)) WriteNPU(hpu,hpuRwgt,basename);
      
      delete infile;
      infile=0, fEventTree=0;    
    } // file loop

    fOutfile->Write();
    fOutfile->Close();

    delete fOutfile;

    if(samp->typev.size()>0 && samp->typev[0]==eMC)
      printf("    Yields for %1.2f/fb:",fLumi/1000.);
    else
      printf("    Yields for data:    ");

    printf("%10.2f\n",nSelEvents[0]);
    cout << endl;
  } // CSample loop

  //--------------------------------------------------------------------------------------------------------------
  // Summary print out
  //==============================================================================================================
  cout << endl;
  cout << " <-> Output saved in " << fNtupDir << "/" << endl;
  if(!fDoNpuRwgt) cout << endl << endl << "Note: not doing npv reweight!" << endl;
  cout << endl; 
}
//----------------------------------------------------------------------------------------
void Selector::InitOutput()
{
  //
  // Set up output ntuple file for this CSample
  //
  TString outfname = fNtupDir + TString("/") + fSnamev[fIsam] + TString("_select.root");
  fOutfile = TFile::Open(outfname,"RECREATE");

  fFailHist = new TH1F("fail","fail",100,-0.5,99.5);

  fOuttree = new TTree("Events","Events");     // output ntuple
  fOuttree->Branch("Events",&fData.runNum,
		 "runNum/i:evtNum:lumiSec:nPV:njets:nbjets:met/F:metphi:mass:dphi:mt:pt:phi:pmet:pvis:lpt1:leta1:lphi1:lpt2:leta2:lphi2:jpt1:jeta1:jphi1:jpt2:jeta2:jphi2:bjpt:bjeta:bjphi:mjj:svfmass:svflpt1:svfleta1:svflphi1:svflpt2:svfleta2:svflphi2:weight:state/I");
  // extra branches
  fOuttree->Branch("npt20jets",&fNpt20jets);
  fOuttree->Branch("btagArray",&fBtagArray);
  fOuttree->Branch("rawMet",&fRawMet);
  fOuttree->Branch("rawprojvar",&fRawprojvar);
  fOuttree->Branch("npuWgt",&fNpuWgt);
}
//----------------------------------------------------------------------------------------
void Selector::InitNPU(TString sfname, TH1D *&hpu, TH1D *&hpuRwgt, vector<Double_t> &puwgtv)
{
  //
  // make the raw and reweighted npu hists
  //
  UInt_t npubins = 55;
  if(fDoNpuRwgt || fMakeNpuHists) {
    if(fMakeNpuHists && sfname.Contains("skim")) { cout << "make npu distribs *before* selection" << endl; assert(0); }
    hpu     = new TH1D("hpu","hpu",npubins,-0.5,npubins-0.5); hpu->Sumw2();		    // raw npu distribution in MC
    hpuRwgt = new TH1D("hpuRwgt","hpuRwgt",npubins,-0.5,npubins-0.5); hpuRwgt->Sumw2();     // reweighted
  }
  if(!fIsdata && (fDoNpuRwgt || fMakeNpuHists)) puwgtv = generate_flat10_weights(fDataNPVfname, fmcNPVfname);
}
//----------------------------------------------------------------------------------------
Double_t Selector::GetWeight(TString sfname, CSample *samp)
{
  //
  // get uncorrected weight (~ xs/entries)
  //
  Double_t treeEntries=-1; // (weight is only initialized for each *file*)
  if(sfname.Contains("_skim.root")) treeEntries = UnskimmedEntries(sfname);		// get entries from unskimmed file
  else                              treeEntries = (Double_t)fEventTree->GetEntries();
  assert(treeEntries>0);
  return fLumi*(samp->xsecv[fIfile])/treeEntries;					// (assumes you've merged filesets)
}
//----------------------------------------------------------------------------------------
void Selector::RemoveSkimName(TString &name)
{
  //
  // remove the skim-label from file name
  //
  TRegexp reg("_[a-zA-Z][a-zA-Z]*_skim");
  if(name.Contains(reg)) {
    TString substr = name(name.Index(reg),name.Length()-name.Index("_skim")); // substr is, eg, "_emu_skim"
    name.ReplaceAll(substr,"_ntuple");
  }
}
//----------------------------------------------------------------------------------------
Double_t Selector::UnskimmedEntries(TString skimname)
{
  //
  // get number of entries in unskimmed tree
  //
  Double_t entries;
  
  RemoveSkimName(skimname);
  TFile unskimmed(skimname);
  assert(unskimmed.IsOpen());
  TTree *tree = 0;
  unskimmed.GetObject("Events",tree);
  assert(tree);
  entries = (Double_t)tree->GetEntries();
  unskimmed.Close();

  return entries;
}
//----------------------------------------------------------------------------------------
void Selector::FillNPU(TH1D *hpu, TH1D *hpuRwgt,vector<Double_t> puwgtv)
{
  //
  // fill the raw and reweighted npu hists with the number of pu interactions in this event
  //
  assert(hpu);
  hpu->Fill(fInfo->nPU); // do this before *any* selections
  Double_t tmpwgt = (fInfo->nPU >= puwgtv.size()) ? 0 : puwgtv[fInfo->nPU];
  assert(hpuRwgt);
  hpuRwgt->Fill(fInfo->nPU,tmpwgt);
}
//----------------------------------------------------------------------------------------
void Selector::WriteNPU(TH1D *hpu, TH1D *hpuRwgt, TString basename)
{
  //
  // rescale and write out npu hists
  //
  assert(hpu); assert(hpuRwgt);

  hpu    ->Scale(1./    hpu->Integral(0,     hpu->GetNbinsX()+1));
  hpuRwgt->Scale(1./hpuRwgt->Integral(0, hpuRwgt->GetNbinsX()+1));

  // write out root files of npu distributions for later use
  if(fMakeNpuHists) {
    TFile puoutfile("npu/"+basename+"-npu.root","recreate");
    hpu->Write();
    puoutfile.Close();
  }

  //
  // make pngs of (estimated npu in data) vs. (npu in mc)
  //
  TCanvas c3("c3","c3");
  // get estimated-for-data npu histogram
  TH1D* data_npu=0; TFile *foofile = TFile::Open(fDataNPVfname); foofile->GetObject("pileup",data_npu); assert(data_npu);
  data_npu->SetDirectory(0); data_npu->Sumw2(); foofile->Close(); delete foofile;
  // scale data estimated npu
  data_npu->Scale(1./data_npu->Integral(0,data_npu->GetNbinsX()+1));
  data_npu->SetMarkerStyle(20);
  data_npu->SetMarkerSize(0.9);
  data_npu->Draw("EP");
  // scale reweighted mc npu
  hpuRwgt->SetLineColor(kBlue);
  hpuRwgt->Draw("histsame");
  hpu->SetLineColor(kRed);
  hpu->Draw("histsame");
  if(fMakeNpuHists || fCheckNpuHists) c3.SaveAs("npu/"+basename+"-npu.png");
	
  delete hpu; delete hpuRwgt;
}
//----------------------------------------------------------------------------------------
Bool_t Selector::EvtFail(ULong64_t failshift)
{
  fEvtFail = (1<<failshift);
  fFailHist->Fill(failshift);
  return kFALSE;
}
//----------------------------------------------------------------------------------------
void Selector::MuonLoop(Bool_t muMu, vector<const mithep::TMuon*> &goodMuonsv, vector<const mithep::TMuon*> &looseMuonsv)
{
  //
  // fill "good" and "loose" muon vectors
  //
  for(Int_t i=0; i<fMuonArr->GetEntriesFast(); i++) {
    const mithep::TMuon *muon = (mithep::TMuon*)((*fMuonArr)[i]);

    looseMuonsv.push_back(muon);

    Double_t minpt = muMu ? fMuMuMuonPt2Min : fEMuMuonPt2Min;
    if(muon->pt < minpt)                           continue;
    if(fabs(muon->eta) > 2.1)                      continue;

    if(passMuonID(muon))  goodMuonsv.push_back(muon);
  }
}
//----------------------------------------------------------------------------------------
void Selector::ElectronLoop(vector<mithep::TElectron*> &goodElectronsv, vector<const mithep::TMuon*> looseMuonsv,
			    vector<const mithep::TMuon*> goodMuonsv)
{
  //
  // fill "good" electron vector
  //
  for(Int_t i=0; i<fElectronArr->GetEntriesFast(); i++) {
    mithep::TElectron *electron = (mithep::TElectron*)((*fElectronArr)[i]);

    if(electron->pt < fElePt2Min)                   continue;
    if(fabs(electron->eta) > 2.5)   	            continue;

    // if(!fIsdata) electron->pt = scaler.pt(electron->eta,electron->pt,escale); // not really proper: applies 42x corrections to 41x

    Bool_t hasMuonTrack=kFALSE;
    for(UInt_t imu=0; imu<goodMuonsv.size(); imu++) {
      if(electron->trkID == goodMuonsv[imu]->trkID) hasMuonTrack=kTRUE;
    }
    if(hasMuonTrack) continue;

    // clean against loose muons
    Bool_t matchLooseMuon=kFALSE;
    for(UInt_t imu=0;imu<looseMuonsv.size();imu++) {
      const mithep::TMuon *mu = looseMuonsv[imu];
      if(toolbox::deltaR(electron->eta,electron->phi,mu->eta,mu->phi) < 0.3) matchLooseMuon=kTRUE;
    }
    if(matchLooseMuon) continue;

    if(passEleID(electron))    goodElectronsv.push_back(electron);
  }
}
//----------------------------------------------------------------------------------------
Bool_t Selector::PassMuMu(vector<const mithep::TMuon*> goodMuonsv, vector<mithep::TElectron*> goodElectronsv,
			  TLorentzVector &lep1, TLorentzVector &lep2, TLorentzVector &dilep, Int_t &finalState,
			  const mithep::TMuon *&mu, const mithep::TMuon *&mu_2)
{
  //
  // select mu mu events and make some useful vectors/flags
  //
  if(goodMuonsv.size()<2 || goodElectronsv.size()>0)		return EvtFail(5);
  if(fDebug) cout << "nmu: " << goodMuonsv.size() << " nele: " << goodElectronsv.size() << endl;

  mu   = goodMuonsv[0];	  
  mu_2 = goodMuonsv[1];	  
							  
  if(mu->pt < fMuonPt1Min || mu_2->pt < fMuMuMuonPt2Min)	return EvtFail(6);
  if(fDebug) cout << "mu pt: " << setw(15) << mu->pt << setw(15) << mu_2->pt << endl;

  //
  // trigger requirements
  //
  ULong64_t trigbits = kHLT_DoubleMu7 | kHLT_Mu13_Mu8;
  if(!(fInfo->triggerBits & trigbits))				return EvtFail(7);
  if(fDebug) cout << "trigbits: " << setw(25) << hex << trigbits << setw(25) << hex << fInfo->triggerBits << endl;

  ULong64_t trigmatch=0;
  if(!fIsdata) trigmatch = (fInfo->triggerBits & kHLT_DoubleMu7) && (mu->hltMatchBits   & kHLT_DoubleMu7_MuObj)
		                                                 && (mu_2->hltMatchBits & kHLT_DoubleMu7_MuObj);
  else {
    if(fInfo->runNum<165085) {
      trigmatch = (fInfo->triggerBits & kHLT_DoubleMu7) && (mu->hltMatchBits   & kHLT_DoubleMu7_MuObj)
                                                        && (mu_2->hltMatchBits & kHLT_DoubleMu7_MuObj);
    }
    else {
      trigmatch = (fInfo->triggerBits & kHLT_Mu13_Mu8) && (mu->hltMatchBits   & (kHLT_Mu13_Mu8_Mu1Obj | kHLT_Mu13_Mu8_Mu2Obj))
                                                       && (mu_2->hltMatchBits & (kHLT_Mu13_Mu8_Mu1Obj | kHLT_Mu13_Mu8_Mu2Obj));
    }
  }
  if(!trigmatch)						return EvtFail(8);
  if(fDebug) cout << "trigmatch: " << hex << setw(25) << mu->hltMatchBits << endl;

  //
  // OS
  //
  if(mu->q == mu_2->q)						return EvtFail(9);
  if(fDebug) cout << "OS: " << endl;

  lep1.SetPtEtaPhiM(mu->pt,   mu->eta,   mu->phi,   0.105658369);
  lep2.SetPtEtaPhiM(mu_2->pt, mu_2->eta, mu_2->phi, 0.105658369);
  dilep = lep1+lep2;

  //
  // inside the z peak
  //
  if(dilep.M()<60.0 || dilep.M()>120.0)				return EvtFail(10);
  if(fDebug) cout << "z mass: " << dilep.M() << endl;

  finalState=kMuMu; 

  return kTRUE;
}
//----------------------------------------------------------------------------------------
Bool_t Selector::PassEmu(vector<const mithep::TMuon*> goodMuonsv, vector<mithep::TElectron*> goodElectronsv,
			 TLorentzVector &lep1, TLorentzVector &lep2, TLorentzVector &dilep, Int_t &finalState,
			 const mithep::TMuon *&mu, mithep::TElectron *&ele)
{
  //
  // select e mu events, and make some useful vectors/flags
  //
  if(goodMuonsv.size()<1 || goodElectronsv.size()<1)		return EvtFail(5);
  if(fDebug) cout << "nmu: " << goodMuonsv.size() << "nele: " << goodElectronsv.size() << endl;

  mu  = goodMuonsv[0];	  
  ele = goodElectronsv[0];	  

  //
  // kinematics
  //
  if(mu->pt < fEMuMuonPt2Min  || ele->pt < fElePt2Min)		return EvtFail(6);
  if(mu->pt < fMuonPt1Min     && ele->pt < fElePt1Min)		return EvtFail(7);
  if(fDebug) cout << "emu kinematics" << endl;

  //
  // trigger requirements
  //
  ULong64_t trigbits = kHLT_Mu17_Ele8_CaloIdL | kHLT_Mu8_Ele17_CaloIdL | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL;
  if(!(fInfo->triggerBits & trigbits))				return EvtFail(8);
  if(fDebug) cout << "trigger: " << hex << setw(25) << fInfo->triggerBits << hex << setw(25) << trigbits << endl;

  ULong64_t muObjBits = kHLT_Mu17_Ele8_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdL_MuObj | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_MuObj;
  if(fInfo->runNum<167000) // trigger matching broken after this run
    if(fIsdata && !(mu->hltMatchBits & muObjBits))		return EvtFail(9);
  if(fDebug) cout << "mu trig match" << endl;
  
  ULong64_t eleObjBits = kHLT_Mu17_Ele8_CaloIdL_EGObj | kHLT_Mu8_Ele17_CaloIdL_EGObj | kHLT_Mu8_Ele17_CaloIdT_CaloIsoVL_EGObj;
  if(fInfo->runNum<167000)
    if(fIsdata && !(ele->hltMatchBits & eleObjBits))		return EvtFail(10);
  if(fDebug) cout << "ele trig match" << endl;
  
  if(fIsdata) {						  
    if(mu->pt  < fMuonPt1Min) {				  
      if(!(fInfo->triggerBits & kHLT_Mu8_Ele17_CaloIdL))	return EvtFail(11);
    }							  
    else if(ele->pt < fElePt1Min) {			  
      if(!(fInfo->triggerBits & kHLT_Mu17_Ele8_CaloIdL))	return EvtFail(11);
    }							  
  }
  if(fDebug) cout << "trig and pt" << endl;
							  
  if(mu->q == ele->q)						return EvtFail(12);
  if(fDebug) cout << "OS: " << endl;

  // assign lead and next-to-lead lepton vectors
  if(mu->pt > ele->pt) {
    lep1.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
    lep2.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
  } else {
    lep1.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
    lep2.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
  }
  dilep = lep1+lep2;

  // assign final state type
  if(ele->pt > mu->pt) finalState=kEleMu; 
  else                 finalState=kMuEle;

  return kTRUE;
}
//----------------------------------------------------------------------------------------
void Selector::JetLoop(TLorentzVector lep1, TLorentzVector lep2, UInt_t &njets, UInt_t &nbjets,
		       const mithep::TJet *&jet1, const mithep::TJet *&jet2, const mithep::TJet *&bjet)
{
  //
  // select jets for vbf and b-tagging use
  //
  for(Int_t i=0; i<fJetArr->GetEntriesFast(); i++) {
    mithep::TJet *jet = (mithep::TJet*)((*fJetArr)[i]);

    // move jet pt up/down by its uncertainty
    if(fJetUnc!=kNo) jet->pt *= (fJetUnc==kDown) ? (1-jet->unc) : (1+jet->unc);

    if(toolbox::deltaR(jet->eta,jet->phi,lep1.Eta(),lep1.Phi()) < 0.3) continue;
    if(toolbox::deltaR(jet->eta,jet->phi,lep2.Eta(),lep2.Phi()) < 0.3) continue;

    if(fabs(jet->eta) > 5) continue;

    Bool_t btagged = IsBtagged(jet);

    // look for b-jets
    if((jet->pt > fBJetPtMin) && (fabs(jet->eta) < 2.4)) { // note: bjet can be the same as jet1 or jet2
      assert(fNpt20jets<kMaxPt20Jets);
      fBtagArray.AddAt(jet->tche,fNpt20jets); fNpt20jets++;
      if(btagged) {
	nbjets++;
	if(!bjet || jet->pt > bjet->pt)
	  bjet = jet; // leading b-jet
      }
    }

    // look for vbf jets
    if(jet->pt > fJetPtMin) {
      njets++;
      if(!jet1 || jet->pt > jet1->pt) { // jet1 is highest pt, jet2 next highest
	jet2 = jet1;
	jet1 = jet;
      } else if(!jet2 || jet->pt > jet2->pt) {
	jet2 = jet;
      }
    }		    
  }
}
//----------------------------------------------------------------------------------------
void Selector::Projections(TLorentzVector lep1, TLorentzVector lep2, TLorentzVector dilep,
			   Bool_t doRecoil, Double_t &projVis, Double_t &projMet, Double_t &met,
			   Double_t &metphi)
{
  //
  // calculate projection variables
  //
  TVector3 lep1v,lep2v,metv;
  lep1v = lep1.Vect();
  lep2v = lep2.Vect();
  metv.SetPtEtaPhi(fInfo->pfMET,0,fInfo->pfMETphi); // uncorrected met
  TVector3 bisector(lep1v.Unit() + lep2v.Unit());
  bisector = bisector.Unit();
  // calculate uncorrected variables
  projVis  = (lep1v+lep2v).Dot(bisector);
  projMet  = metv.Dot(bisector);
  fRawprojvar  = 0.85*projVis - projMet;
  fRawMet = fInfo->pfMET;
	  
  met = fInfo->pfMET;
  metphi = fInfo->pfMETphi;
  // apply recoil corrections
  if(doRecoil) fCorrector.Correct(met,metphi,fGen->vpt_a,fGen->vphi_a,dilep.Pt(),dilep.Phi());
  metv.SetPtEtaPhi(met,0,metphi); // corrected met vector
  projMet  =  metv.Dot(bisector); // corrected projected met
}
//----------------------------------------------------------------------------------------
Double_t Selector::GetTrigEff(const mithep::TMuon *mu, const mithep::TElectron *ele)
{
  //
  // get trigger efficiency corrections
  //
  Double_t t1eff=1,t2eff=1,trigeff=1;
  if(!fIsdata) { // this is a scale factor, not an efficiency, for 42x
    t1eff = t2eff = 0.991*0.991;
  }
  else { cout << "Error: no trigger efficiency defined." << endl; assert(0); }

  if(mu->pt < fMuonPt1Min)        trigeff = t1eff;
  else if(ele->pt > fElePt1Min)   trigeff = t1eff + t2eff*(1-t1eff);
  else                            trigeff = t2eff;

  return trigeff;
}
//----------------------------------------------------------------------------------------    
vector<Double_t> Selector::generate_flat10_weights(TString datafname, TString mcfname)
{
  //
  // Initialize npu weights
  //

  TH1D* data_npu = 0;
  TH1D* mc_npu = 0;

  TFile *infile = TFile::Open(datafname); assert(infile->IsOpen());
  infile->GetObject("pileup",data_npu);   assert(data_npu);
  data_npu->SetDirectory(0);
  infile->Close();
  data_npu->Scale(1./data_npu->Integral(0,data_npu->GetNbinsX()+1));

  infile = TFile::Open(mcfname);
  if(!infile) {
    infile = TFile::Open("npu/p11-vvj-v1g1-pu_ntuple-npu.root");
    cout << endl << endl << "Warning: using NPU from -vvj- sample. Run again to use npu from this file." << endl << endl;
    assert(infile->IsOpen());
  }
  infile->GetObject("hpu",mc_npu); assert(mc_npu);
  mc_npu->SetDirectory(0);
  infile->Close();


  const UInt_t nbins = TMath::Min(data_npu->GetNbinsX(),mc_npu->GetNbinsX());
  
  vector<Double_t> result(nbins);
  Double_t sum = 0;
  for(UInt_t npu=0; npu<nbins; ++npu){
    Double_t data_wgt = data_npu->GetBinContent(data_npu->GetXaxis()->FindBin(npu));                              
    Double_t   mc_wgt =   mc_npu->GetBinContent(mc_npu->GetXaxis()->FindBin(npu));                              
    result[npu] = (mc_wgt==0) ? 0 : data_wgt / mc_wgt;
    sum += result[npu];
  }

  return result;
}
//----------------------------------------------------------------------------------------
// trigger efficiency correction functions (numbers from kevin)
//----------------------------------------------------------------------------------------
Double_t Selector::eleTrigEff(const mithep::TElectron *ele)
{
  if((fabs(ele->eta) > 2.5) || (ele->pt < 10)) { cout << "ele kinematics out of range" << endl; assert(0); }
  else if(ele->pt > 20) {
    if(fabs(ele->eta) < 1.479) return 0.9970;
    else                       return 0.999;
  }
  else if(ele->pt > 15) {
    if(fabs(ele->eta) < 1.479) return 0.9947;
    else                       return 1;
  }
  else {
    if(fabs(ele->eta) < 1.479) return 0.978;
    else                       return 1;
  }
  
}
//----------------------------------------------------------------------------------------
Double_t Selector::muTrigEff(const mithep::TMuon *mu)
{
  if((fabs(mu->eta) > 2.4) || (mu->pt < 10)) { cout << "mu kinematics out of range" << endl; assert(0); }
  else if(fabs(mu->eta) > 1.2) {
    if(mu->pt > 20)  return 0.9548;
    else             return 0.9478;
  }
  else if(fabs(mu->eta) > 0.8) {
    if(mu->pt > 20)  return 0.9488;
    else             return 0.9609;
  }
  else {
    if(mu->pt > 20)  return 0.9784;
    else             return 0.9674;
  }
}
//----------------------------------------------------------------------------------------
Bool_t Selector::IsBtagged(mithep::TJet *jet)
{
  //
  // Calculate if we will treat this jet as b-tagged (includes uncertainties)
  //

  //          mistag                         scale factor
  // TCHEM  0.0175 \pm .0003 \pm .0038      1.21 \pm .02 \pm .17
  //          btag eff.                      scale factor
  // TCHEM  0.63 \pm 0.01                   0.93 \pm 0.02 \pm 0.07

  Bool_t btagged;
  Double_t demoteProb=0; // ~probability to demote from tagged
  if(fBtagEff==kNo)        demoteProb = 1-0.93;  // SF = 0.93 -> 0.07 = (prob to demote from tagged status)
  else if(fBtagEff==kDown) demoteProb = 1-0.93+0.07;
  else if(fBtagEff==kUp)   demoteProb = 1-0.93-0.07;
  Double_t promoteProb=0; // ~probability to promote to tagged
  if(fMistag==kNo)         promoteProb = (1.21-1)*0.0145/(1-0.0145);  // (1-SF)*mistag = (prob. to promote to tagged status)*(1-mistag)
  else if(fMistag==kDown)  promoteProb = (1.21-1+0.17)*0.0145/(1-0.0145);
  else if(fMistag==kUp)    promoteProb = (1.21-1-0.17)*0.0145/(1-0.0145);
                   
  if(fIsdata) {
    if(jet->tche>3.3)				btagged = kTRUE;
    else					btagged = kFALSE;
  } else { // MC
    if(abs(jet->mcFlavor)==5) {
      if(jet->tche>3.3) {
      if(fRandm.Uniform()>demoteProb)		btagged = kTRUE;  // leave it tagged
      else					btagged = kFALSE; // demote it
      } else					btagged = kFALSE; // leave it untagged
    } else { // not bjet
      if(jet->tche>3.3)				btagged = kTRUE;  // leave it tagged
      else if(fRandm.Uniform()<promoteProb)	btagged = kTRUE;  // promote to tagged
      else					btagged = kFALSE; // leave it untagged
    }
  }

  return btagged;
}
//----------------------------------------------------------------------------------------
// K-factors [NOT fully implemented]
//----------------------------------------------------------------------------------------
TH1D* Selector::kfInit(TString kfilename, Int_t mH)
{
  TFile kfile(kfilename);
  TH1D *kf=0;
  char kfname[100];
  sprintf(kfname, "KFactor_PowhegToHQT_mH%d", mH); // file doesn't have all of our mass points
  kf = (TH1D*)(kfile.Get(kfname)); assert(kf);
  kf->SetDirectory(0);
  kfile.Close();

  return kf;
//weightFactor *= HiggsPtKFactor->GetBinContent( HiggsPtKFactor->GetXaxis()->FindFixBin(higgsPt));
}
//--------------------------------------------------------------------------------------------------
Double_t Selector::kfValue(const Double_t pt, const TH1D* hKF)
{
  cout << "needs to be checked -- only halfway implemented" << endl; assert(0);
  if(pt < hKF->GetBinLowEdge(1)) {
    return hKF->GetBinContent(0);
  
  } else if(pt > hKF->GetBinLowEdge(hKF->GetNbinsX())) { // I think this is the wrong bin, but it won't really matter
    return hKF->GetBinContent(hKF->GetNbinsX());
  
  } else {
    for(Int_t ibin=1; ibin<=hKF->GetNbinsX(); ibin++) {
      if(pt >= hKF->GetBinLowEdge(ibin) && pt < hKF->GetBinLowEdge(ibin+1)) {
        return hKF->GetBinContent(ibin);
      }
    }
  }
  return 1;
}
Int_t Selector::higgsmass(TString basename)
{
  stringstream ss(basename(TRegexp("[0-9][0-9]*"),3).Data());
  Int_t mass;
  ss >> mass;
  if(basename.Contains("-gf-")) assert(mass>85 && mass<1200);
  return mass;
}
	// //----------------------------------------------------------------------------------------

	// // boost into the rest frame
	// TLorentzVector nlep,plep,bos;
	// Double_t boseta = eta(fGen->vpt,fGen->vy,fGen->vphi,fGen->vmass);
	// bos.SetPtEtaPhiM(fGen->vpt,boseta,fGen->vphi,fGen->vmass);
	// if(mu->q > 0) {
	//   plep.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
	//   nlep.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
	// } else {
	//   plep.SetPtEtaPhiM(ele->pt, ele->eta, ele->phi, 0.000511);
	//   nlep.SetPtEtaPhiM(mu->pt,  mu->eta,  mu->phi,  0.105658369);
	// }
	// // plep.Boost(0,0,v(lep1.Theta(),2*pi-lep2.Theta()));
	// // nlep.Boost(0,0,v(lep1.Theta(),2*pi-lep2.Theta()));
	// plep.Boost(-bos.BoostVector());
	// nlep.Boost(-bos.BoostVector());

	// // bos.Boost(0,0,v(lep1.Theta(),2*pi-lep2.Theta()));
	// // Double_t v_me = v(lep1.Theta(),2*pi-lep2.Theta());
	// // Double_t v_re = (-bos.BoostVector()).Pz();
	// // cout << "       " << v(lep1.Theta(),2*pi-lep2.Theta()) << endl;
	// // cout << "       "; (-bos.BoostVector()).Print();
	// // cout << endl;
	// Double_t cths = bos.Vect()*plep.Vect()/(bos.P()*plep.P());
	// // Double_t e1e2 = plep.E()/nlep.E();
	// hb.Fill(cths);
	// if(ientry>200000) break;
	// //----------------------------------------------------------------------------------------


    // // write out cutflow
    // FILE *fcut = fopen(outputDir+"/cutflow.txt","a");
    // Double_t kssFakeWgt = 1;
    // if(fSnamev[fIsam].Contains("ss-fakes")) kssFakeWgt = 1.312;
    // fprintf(fcut,"%45s%20.2f%20.2f%20.2f%20.2f%20.2f%20.2f%20.2f\n",
    // 	    fSnamev[fIsam].Data(),kssFakeWgt*counter[0],kssFakeWgt*counter[1],
    // 	    kssFakeWgt*counter[2],kssFakeWgt*counter[3],kssFakeWgt*counter[4],kssFakeWgt*counter[5],kssFakeWgt*counter[6]);
    // fclose(fcut);
