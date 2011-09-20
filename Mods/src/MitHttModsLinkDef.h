
// $Id: MitHttModsLinkDef.h,v 1.6 2011/02/18 14:58:22 dkralph Exp $

#ifndef MITHTT_MODS_LINKDEF_H
#define MITHTT_MODS_LINKDEF_H

#include "MitHtt/Mods/interface/AcceptanceMod.h"
#include "MitHtt/Mods/interface/MatchingMod.h"
#include "MitHtt/Mods/interface/MCPartonColPubMod.h"
#include "MitHtt/Mods/interface/ZeeAnalysis.h"
#include "MitHtt/Mods/interface/ZmmAnalysis.h"
#include "MitHtt/Mods/interface/ZttAnalysis.h"
#include "MitHtt/Mods/interface/EMUAnalysis.h"
#include "MitHtt/Mods/interface/EMUAnalysis2.h"
#include "MitHtt/Mods/interface/TauPlayAnalysis.h"
#include "MitHtt/Mods/interface/PlotInfo.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::AcceptanceMod+;
#pragma link C++ class mithep::MatchingMod+;
#pragma link C++ class mithep::MCPartonColPubMod+;
#pragma link C++ class mithep::ZeeAnalysis+;
#pragma link C++ class mithep::ZmmAnalysis+;
#pragma link C++ class mithep::ZttAnalysis+;
#pragma link C++ class mithep::EMUAnalysis+;
#pragma link C++ class mithep::EMUAnalysis2+;
#pragma link C++ class mithep::TauPlayAnalysis+;
#pragma link C++ class mithep::PlotInfo+;
#endif