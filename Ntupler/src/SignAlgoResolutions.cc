#include <math.h>
#include <string>
#include <cstdlib>
#include <sstream>
#include <iostream>
#include <sys/stat.h>

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/PythonParameterSet/interface/MakeParameterSets.h"

#include "MitHtt/Ntupler/interface/SignAlgoResolutions.h"


mithep::SignAlgoResolutions::SignAlgoResolutions() :
  functionmap_(), ptResol_(0), phiResol_(0),
  fName("/build/pharris/CMSSW_4_2_4_patch1/src/RecoMET/METProducers/python/METSigParams_cfi.py")
{
  addResolutions();
}

double 
mithep::SignAlgoResolutions::eval(const mithep::resolutionType& type, const mithep::resolutionFunc& func, const double& et, const double& phi, const double& eta) const 
{
  // derive p from et and eta;
  double theta = 2*atan(exp(-eta));
  // rough assumption: take e and p equivalent...
  double p     = et / sin(theta); 
  return eval(type,func,et,phi,eta,p);
}

double 
mithep::SignAlgoResolutions::eval(const mithep::resolutionType& type, const mithep::resolutionFunc& func, const double& et, const double& phi, const double& eta, const double& p) const 
{
  functionPars x(4);
  x[0]=et;
  x[1]=phi;
  x[2]=eta;
  x[3]=p;
  //  std::cout << "getting function of type " << type << " " << func << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << std::endl;
  return getfunc(type,func,x);
}

metsig::SigInputObj  
mithep::SignAlgoResolutions::evalPF(const mithep::PFCandidate* candidate) const 
{
  double eta = candidate->Eta();
  double phi = candidate->Phi();
  double et  = candidate->Et();//*sin(candidate->Theta());
  mithep::resolutionType thetype;
  std::string name;
  int type = candidate->ObjType();
  switch (type) 
    {
    case mithep::PFCandidate::eHadron: 
      thetype=PFtype1; name="PFChargedHadron"; 
      break;
    case mithep::PFCandidate::eElectron: 
      thetype=PFtype2; name="PFChargedEM"; 
      break;
    case mithep::PFCandidate::eMuon: 
      thetype=PFtype3; name="PFMuon"; 
      break;
    case mithep::PFCandidate::eGamma: 
      thetype=PFtype4; name="PFNeutralEM"; 
      break;
    case mithep::PFCandidate::eNeutralHadron: 
      thetype=PFtype5; name="PFNeutralHadron"; 
      break;
    case mithep::PFCandidate::eHadronHF: 
      thetype=PFtype6; name="PFtype6"; 
      break;
    case mithep::PFCandidate::eEGammaHF:
      thetype=PFtype7; name="PFtype7"; 
      break;
    default:
      thetype=PFtype7; name="PFunknown"; 
      break;
    }
  //d_phi here is the error on the phi component of et
  double d_et=0, d_phi=0; 
  const mithep::Track* trackRef = candidate->Trk();
  if(trackRef != 0 && type != 2){
    d_et = trackRef->PtErr();
    d_phi = et*trackRef->Phi0Err();
  }
  else{
    d_et = eval(thetype,ET,et,phi,eta);
    d_phi = eval(thetype,PHI,et,phi,eta);
  }
  
  metsig::SigInputObj resultingobj(name,et,phi,d_et,d_phi);
  return resultingobj;
}

metsig::SigInputObj
mithep::SignAlgoResolutions::evalPFJet(const mithep::PFJet* jet) const
{
  double jpt  = jet->Pt();
  double jphi = jet->Phi();
  double jeta = jet->Eta();
  double jdeltapt = 999.;
  double jdeltapphi = 999.;
  
  if(jpt<ptResolThreshold_ && jpt<20.){ //use temporary fix for low pT jets
    double feta = TMath::Abs(jeta);
    int ieta = feta<5.? int(feta/0.5) : 9; //bin size = 0.5 
    int ipt  = jpt>3. ? int(jpt-3./2) : 0; //bin size =2, starting from ptmin=3GeV
    jdeltapt   = jdpt[ieta][ipt];
    jdeltapphi = jpt*jdphi[ieta][ipt];
  }
  else{
    TF1* fPtEta  = ptResol_->parameterEta("sigma",jeta);
    TF1* fPhiEta = phiResol_->parameterEta("sigma",jeta);
    jdeltapt   = jpt>ptResolThreshold_ ? jpt*fPtEta->Eval(jpt)  : jpt*fPtEta->Eval(ptResolThreshold_);
    jdeltapphi = jpt>ptResolThreshold_ ? jpt*fPhiEta->Eval(jpt) : jpt*fPhiEta->Eval(ptResolThreshold_);
    delete fPtEta;
    delete fPhiEta;
  }
  
  std::string inputtype = "jet";
  metsig::SigInputObj obj_jet(inputtype,jpt,jphi,jdeltapt,jdeltapphi);
  //std::cout << "RESOLUTIONS JET: " << jpt << "   " << jphi<< "   " <<jdeltapt << "   " << jdeltapphi << std::endl;
  return obj_jet;
}

void 
mithep::SignAlgoResolutions::addResolutions()
{
  //using namespace std;
  edm::ParameterSet iConfig = edm::readPSetsFrom(fName)->getParameter<edm::ParameterSet>("METSignificance_params");
  // Jet Resolutions - for now load from the files. Migrate to EventSetup asap.
  mithep::SignAlgoResolutions::initializeJetResolutions( );
  ptResolThreshold_ = iConfig.getParameter<double>("ptresolthreshold");
  //get temporary low pT pfjet resolutions
  for (int ieta=0; ieta<10; ieta++){
    jdpt[ieta] = iConfig.getParameter<std::vector<double> >(Form("jdpt%d", ieta));
    jdphi[ieta] = iConfig.getParameter<std::vector<double> >(Form("jdphi%d", ieta));
  }
  // for now: do this by hand - this can obviously also be done via ESSource etc.
  functionPars etparameters(3,0);
  functionPars phiparameters(1,0);
  // set the parameters per function:
  // ECAL, BARREL:
  std::vector<double> ebet = iConfig.getParameter<std::vector<double> >("EB_EtResPar");
  std::vector<double> ebphi = iConfig.getParameter<std::vector<double> >("EB_PhiResPar");
  
  etparameters[0]=ebet[0];
  etparameters[1]=ebet[1];
  etparameters[2]=ebet[2];
  phiparameters[0]=ebphi[0];
  addfunction(caloEB,ET,etparameters);
  addfunction(caloEB,PHI,phiparameters);
  // ECAL, ENDCAP:
  std::vector<double> eeet = iConfig.getParameter<std::vector<double> >("EE_EtResPar");
  std::vector<double> eephi = iConfig.getParameter<std::vector<double> >("EE_PhiResPar");
  
  etparameters[0]=eeet[0];
  etparameters[1]=eeet[1];
  etparameters[2]=eeet[2];
  phiparameters[0]=eephi[0];
  addfunction(caloEE,ET,etparameters);
  addfunction(caloEE,PHI,phiparameters);
  // HCAL, BARREL:
  std::vector<double> hbet = iConfig.getParameter<std::vector<double> >("HB_EtResPar");
  std::vector<double> hbphi = iConfig.getParameter<std::vector<double> >("HB_PhiResPar");
  
  etparameters[0]=hbet[0];
  etparameters[1]=hbet[1];
  etparameters[2]=hbet[2];
  phiparameters[0]=hbphi[0];
  addfunction(caloHB,ET,etparameters);
  addfunction(caloHB,PHI,phiparameters);
  // HCAL, ENDCAP:
  std::vector<double> heet = iConfig.getParameter<std::vector<double> >("HE_EtResPar");
  std::vector<double> hephi = iConfig.getParameter<std::vector<double> >("HE_PhiResPar");
  
  etparameters[0]=heet[0];
  etparameters[1]=heet[1];
  etparameters[2]=heet[2];
  phiparameters[0]=hephi[0];
  addfunction(caloHE,ET,etparameters);
  addfunction(caloHE,PHI,phiparameters);
  // HCAL, Outer
  std::vector<double> hoet = iConfig.getParameter<std::vector<double> >("HO_EtResPar");
  std::vector<double> hophi = iConfig.getParameter<std::vector<double> >("HO_PhiResPar");
  
  etparameters[0]=hoet[0];
  etparameters[1]=hoet[1];
  etparameters[2]=hoet[2];
  phiparameters[0]=hophi[0];
  addfunction(caloHO,ET,etparameters);
  addfunction(caloHO,PHI,phiparameters);
  // HCAL, Forward
  std::vector<double> hfet = iConfig.getParameter<std::vector<double> >("HF_EtResPar");
  std::vector<double> hfphi = iConfig.getParameter<std::vector<double> >("HF_PhiResPar");
  
  etparameters[0]=hfet[0];
  etparameters[1]=hfet[1];
  etparameters[2]=hfet[2];
  phiparameters[0]=hfphi[0];
  addfunction(caloHF,ET,etparameters);
  addfunction(caloHF,PHI,phiparameters);
  
  // PF objects:
  // type 1:
  std::vector<double> pf1et = iConfig.getParameter<std::vector<double> >("PF_EtResType1");
  std::vector<double> pf1phi = iConfig.getParameter<std::vector<double> >("PF_PhiResType1");
  etparameters[0]=pf1et[0];
  etparameters[1]=pf1et[1];
  etparameters[2]=pf1et[2];
  phiparameters[0]=pf1phi[0];
  addfunction(PFtype1,ET,etparameters);
  addfunction(PFtype1,PHI,phiparameters);
  
  // PF objects:
  // type 2:
  std::vector<double> pf2et = iConfig.getParameter<std::vector<double> >("PF_EtResType2");
  std::vector<double> pf2phi = iConfig.getParameter<std::vector<double> >("PF_PhiResType2");
  etparameters[0]=pf2et[0];
  etparameters[1]=pf2et[1];
  etparameters[2]=pf2et[2];
  phiparameters[0]=pf2phi[0];
  addfunction(PFtype2,ET,etparameters);
  addfunction(PFtype2,PHI,phiparameters);
  
  // PF objects:
  // type 3:
  std::vector<double> pf3et = iConfig.getParameter<std::vector<double> >("PF_EtResType3");
  std::vector<double> pf3phi = iConfig.getParameter<std::vector<double> >("PF_PhiResType3");
  etparameters[0]=pf3et[0];
  etparameters[1]=pf3et[1];
  etparameters[2]=pf3et[2];
  phiparameters[0]=pf3phi[0];
  addfunction(PFtype3,ET,etparameters);
  addfunction(PFtype3,PHI,phiparameters);
  
  // PF objects:
  // type 4:
  std::vector<double> pf4et = iConfig.getParameter<std::vector<double> >("PF_EtResType4");
  std::vector<double> pf4phi = iConfig.getParameter<std::vector<double> >("PF_PhiResType4");
  etparameters[0]=pf4et[0];
  etparameters[1]=pf4et[1];
  etparameters[2]=pf4et[2];
  //phiparameters[0]=pf4phi[0];
  addfunction(PFtype4,ET,etparameters);
  addfunction(PFtype4,PHI,pf4phi); //use the same functional form for photon phi error as for pT, pass whole vector
  
  // PF objects:
  // type 5:
  std::vector<double> pf5et = iConfig.getParameter<std::vector<double> >("PF_EtResType5");
  std::vector<double> pf5phi = iConfig.getParameter<std::vector<double> >("PF_PhiResType5");
  etparameters[0]=pf5et[0];
  etparameters[1]=pf5et[1];
  etparameters[2]=pf5et[2];
  phiparameters[0]=pf5phi[0];
  addfunction(PFtype5,ET,etparameters);
  addfunction(PFtype5,PHI,pf5phi);
  
  // PF objects:
  // type 6:
  std::vector<double> pf6et = iConfig.getParameter<std::vector<double> >("PF_EtResType6");
  std::vector<double> pf6phi = iConfig.getParameter<std::vector<double> >("PF_PhiResType6");
  etparameters[0]=pf6et[0];
  etparameters[1]=pf6et[1];
  etparameters[2]=pf6et[2];
  phiparameters[0]=pf6phi[0];
  addfunction(PFtype6,ET,etparameters);
  addfunction(PFtype6,PHI,phiparameters);
  
  // PF objects:
  // type 7:
  std::vector<double> pf7et = iConfig.getParameter<std::vector<double> >("PF_EtResType7");
  std::vector<double> pf7phi = iConfig.getParameter<std::vector<double> >("PF_PhiResType7");
  etparameters[0]=pf7et[0];
  etparameters[1]=pf7et[1];
  etparameters[2]=pf7et[2];
  phiparameters[0]=pf7phi[0];
  addfunction(PFtype7,ET,etparameters);
  addfunction(PFtype7,PHI,phiparameters);
  
  return;
}

void 
mithep::SignAlgoResolutions::addfunction(mithep::resolutionType type, mithep::resolutionFunc func, functionPars parameters)
{
  functionCombo mypair(type,func);
  functionmap_[mypair]=parameters;
}

double 
mithep::SignAlgoResolutions::getfunc(const mithep::resolutionType& type, const mithep::resolutionFunc& func, functionPars& x) const
{  
  double result=0;
  functionCombo mypair(type,func);
  if(functionmap_.count(mypair)==0){
    return result;
  }
  
  functionPars values = (functionmap_.find(mypair))->second;
  switch ( func ){
  case mithep::ET :
    return EtFunction(x, values);
  case mithep::PHI :
    return PhiFunction(x, values);
  case mithep::TRACKP :
    return PFunction(x, values);
  case mithep::CONSTPHI :
    return PhiConstFunction(x, values);
  } 
  //  std::cout << "returning function " << type << " " << func << " " << result << " " << x[0] << std::endl; 
  return result;
}

double 
mithep::SignAlgoResolutions::EtFunction(const functionPars& x, const functionPars& par) const
{
  if(par.size()<3)
    return 0.;
  if(x.size()<1)
    return 0.;
  double et=x[0];
  if(et<=0.)
    return 0.;
  double result = et*sqrt((par[2]*par[2])+(par[1]*par[1]/et)+(par[0]*par[0]/(et*et)));
  return result;
}

double 
mithep::SignAlgoResolutions::PhiFunction(const functionPars& x, const  functionPars& par) const
{
  double et=x[0];
  if(et<=0.){
    return 0.;
  }
  //if 1 parameter is C provided, returns C*pT, if three parameters N, S, C are provided, it returns the usual resolution value, as for sigmaPt
  if(par.size()!=1 && par.size()!=3){//only 1 or 3 parameters supported for phi function
    return 0.;
  }
  else if(par.size()==1){
    return par[0]*et;
  }
  else{
    return et*sqrt((par[2]*par[2])+(par[1]*par[1]/et)+(par[0]*par[0]/(et*et)));
  }
}

double 
mithep::SignAlgoResolutions::PFunction(const functionPars& x, const functionPars& par) const
{
  // not currently implemented
  return 0;
}

double 
mithep::SignAlgoResolutions::PhiConstFunction(const functionPars& x, const functionPars &par) const
{
  return par[0];
}

void
mithep::SignAlgoResolutions::initializeJetResolutions()
{  
  //using namespace std;
  edm::ParameterSet iConfig = edm::readPSetsFrom(fName)->getParameter<edm::ParameterSet>("METSignificance_params");
  // only reinitialize the resolution if the pointers are zero
  if ( ptResol_ == 0 ) {
    std::string resolutionsAlgo  = iConfig.getParameter<std::string>("resolutionsAlgo");     
    std::string resolutionsEra   = iConfig.getParameter<std::string>("resolutionsEra");     
    
    std::string cmssw_base(getenv("CMSSW_BASE"));
    std::string cmssw_release_base(getenv("CMSSW_RELEASE_BASE"));
    std::string path = cmssw_base + "/src/CondFormats/JetMETObjects/data";
    struct stat st;
    if (stat(path.c_str(),&st)!=0) {
      path = cmssw_release_base + "/src/CondFormats/JetMETObjects/data";
    }
    if (stat(path.c_str(),&st)!=0) {
      std::cerr << "ERROR: tried to set path but failed, abort." << std::endl;
    }    
    std::string era(resolutionsEra);
    std::string alg(resolutionsAlgo);
    std::string ptFileName  = path + "/" + era + "_PtResolution_" +alg+".txt";
    std::string phiFileName = path + "/" + era + "_PhiResolution_"+alg+".txt";
    
    ptResol_ = new JetResolution(ptFileName,false);
    phiResol_ = new JetResolution(phiFileName,false);
  }
}
