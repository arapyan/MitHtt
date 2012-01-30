#ifndef MITHTT_NTUPLER_LINKDEF_H
#define MITHTT_NTUPLER_LINKDEF_H
#include "MitHtt/Ntupler/interface/BambuGenDumperMod.hh"
#include "MitHtt/Ntupler/interface/HttNtuplerMod.hh"
#include "MitHtt/Ntupler/interface/TEventInfo.hh"
#include "MitHtt/Ntupler/interface/TGenInfo.hh"
#include "MitHtt/Ntupler/interface/TMuon.hh"
#include "MitHtt/Ntupler/interface/TElectron.hh"
#include "MitHtt/Ntupler/interface/TJet.hh"
#include "MitHtt/Ntupler/interface/TPhoton.hh"
#include "MitHtt/Ntupler/interface/TVertex.hh"
#include "MitHtt/Ntupler/interface/TNSVFit.hh"
#include "MitHtt/Ntupler/interface/MetSignificance.hh"
#include "MitHtt/Ntupler/interface/TSVFitter.hh"
#include "MitHtt/Ntupler/interface/svFitAuxFunctions.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;
#pragma link C++ namespace SVfit_namespace;

#pragma link C++ class mithep::BambuGenDumperMod+;
#pragma link C++ class mithep::HttNtuplerMod+;
#pragma link C++ class mithep::TEventInfo+;
#pragma link C++ class mithep::TGenInfo+;
#pragma link C++ class mithep::TMuon+;
#pragma link C++ class mithep::TElectron+;
#pragma link C++ class mithep::TJet+;
#pragma link C++ class mithep::TPhoton+;
#pragma link C++ class mithep::TVertex+;
#pragma link C++ class mithep::TNSVFit+;
#pragma link C++ class mithep::MetSignificance+;
#pragma link C++ class mithep::TSVFitter+;
#endif
