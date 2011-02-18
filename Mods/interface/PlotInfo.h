//--------------------------------------------------------------------------------------------------
// $Id: $
// holds some trivial plot information to avoid typing the stupid info twice in your module.
//--------------------------------------------------------------------------------------------------


#ifndef MITHTT_MODS_PLOTINFO_H
#define MITHTT_MODS_PLOTINFO_H

#include <string>
#include <TObject.h>

#include "MitAna/TreeMod/interface/BaseMod.h" 
#include "MitAna/DataTree/interface/CollectionsFwd.h"
#include "MitAna/DataTree/interface/CaloJetCol.h"


namespace mithep
{
  class PlotInfo : public TObject
  {
  public:
    PlotInfo(const char* name,const char *title,int nbins,double xmin,double xmax) :
      fname(name),ftitle(title),fnbins(nbins),fxmin(xmin),fxmax(xmax)
      {};
    
    string fname;
    string ftitle;
    int fnbins;
    double fxmin;
    double fxmax;
    ClassDef(PlotInfo, 1)
  };
}
#endif
