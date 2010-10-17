// $Id: $

#ifndef MITHTT_MODS_LINKDEF_H
#define MITHTT_MODS_LINKDEF_H
#include "MitHtt/Mods/interface/HttAnalysis.h"
#include "MitHtt/Mods/interface/ZeeAnalysis.h"
#include "MitHtt/Mods/interface/ZmmAnalysis.h"
#include "MitHtt/Mods/interface/ZttAnalysis.h"
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
#endif
