/*
 *  http://www.slac.stanford.edu/BFROOT/www/doc/workbook_backup_010108/examples/ex1/workdir/RELEASE/workdir/PARENT/shtmp/Linux24SL3_i386_gcc323/RooFitBabar/RooErf.hh
 */
#ifndef ROO_ERF
#define ROO_ERF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"

class RooRealVar;

class RooErf : public RooAbsPdf {
public:
  RooErf(){}
  RooErf(const char *name, const char *title,
	     RooAbsReal& m, RooAbsReal& m0, RooAbsReal& width);

  RooErf(const RooErf& other, const char *name=0) ;
  virtual TObject *clone(const char *newname) const { 
    return new RooErf(*this,newname); }
  inline virtual ~RooErf() { }

protected:

  RooRealProxy _m ;
  RooRealProxy _m0 ;
  RooRealProxy _width;

  Double_t evaluate() const;
  
private:
  
  ClassDef(RooErf,1) // D*-D0 mass difference bg PDF
};

#endif
