/*
 *  http://babar-hn.slac.stanford.edu:5090/cgi-bin/internal/cvsweb.cgi/RooFitBabar/RooErf.cc?rev=1.5;content-type=text%2Fplain;only_with_tag=MAIN
 */
 
#include "RooFit.h"
#include "Riostream.h"
#include <math.h>

#include "RooErf.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooAbsFunc.h"

ClassImp(RooErf);

RooErf::RooErf(const char *name, const char *title,
		       RooAbsReal& m, RooAbsReal& m0,
		       RooAbsReal& width) :
  RooAbsPdf(name,title),
  _m("dm","Dstar-D0 Mass Diff",this,m),
  _m0("dm0","Threshold",this,m0),
  _width("c","Shape Parameter",this,width)
{
}

RooErf::RooErf(const RooErf& other, const char *name) :
  RooAbsPdf(other,name), _m("m",this,other._m), _m0("m0",this,other._m0),
  _width("width",this,other._width)
{
}

Double_t RooErf::evaluate() const
{
  Double_t func = 0.5 * ( 1 - erf((_m - _m0)/_width*sqrt(2.)) );
  return ( func > 1.0e-06 ? func : 1.0e-06 ) ;
}
