
// $Id: $

#ifndef MITHTT_MODS_LINKDEF_H
#define MITHTT_MODS_LINKDEF_H
#include "MitHtt/Mods/interface/HttAnalysis.h"
#include "MitHtt/Mods/interface/ZeeAnalysis.h"
#include "MitHtt/Mods/interface/ZmmAnalysis.h"
#include "MitHtt/Mods/interface/ZttAnalysis.h"
#include "MitHtt/Mods/interface/EMUAnalysis.h"
#include "MitHtt/Mods/interface/TauPlayAnalysis.h"
#include "MitHtt/Mods/interface/WWTauAnalysis.h"
#include "MitHtt/Mods/interface/WWTauAnalysisRepeat.h"
#include "MitHtt/Mods/interface/MergeLeptonsFromAnywhereMod.h"
#include "MitHtt/Mods/interface/MessMod.h"
#include "MitHtt/Mods/interface/MessMod2.h"
#include "MitHtt/Mods/interface/PlotInfo.h"
#endif

#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;
#pragma link C++ nestedtypedef;
#pragma link C++ namespace mithep;

#pragma link C++ class mithep::HttAnalysis+;
#pragma link C++ class mithep::ZeeAnalysis+;
#pragma link C++ class mithep::ZmmAnalysis+;
#pragma link C++ class mithep::ZttAnalysis+;
#pragma link C++ class mithep::EMUAnalysis+;
#pragma link C++ class mithep::TauPlayAnalysis+;
#pragma link C++ class mithep::WWTauAnalysis+;
#pragma link C++ class mithep::WWTauAnalysisRepeat+;
#pragma link C++ class mithep::MergeLeptonsFromAnywhereMod+;
#pragma link C++ class mithep::MessMod+;
#pragma link C++ class mithep::MessMod2+;
#pragma link C++ class mithep::PlotInfo+;
#endif
