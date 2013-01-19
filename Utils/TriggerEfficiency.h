#ifndef TRIGGEREFFICIENCY_H
#define TRIGGEREFFICIENCY_H

#include <vector>
#include <TMath.h>

using namespace std;

namespace mithep
{
  class TrigEff
  {
  public:
    virtual double eff(double pt, double eta) = 0;
  };


  class CBTurnOn : public TrigEff
  {
  public:
    CBTurnOn(double m0_, double sigma_, double alpha_, double n_, double norm_) :
      m0(m0_), sigma(sigma_), alpha(alpha_), n(n_), norm(norm_)
    {}

    double eff(double pt, double eta)
    {
      double m = pt;

      const double sqrtPiOver2 = 1.2533141373;
      const double sqrt2 = 1.4142135624;
      double sig = fabs((double) sigma);
      double t = (m - m0)/sig;
      if(alpha < 0)
        t = -t;
      double absAlpha = fabs(alpha/sig);
      double a = TMath::Power(n/absAlpha,n)*exp(-0.5*absAlpha*absAlpha);
      double b = absAlpha - n/absAlpha;
      double ApproxErf;
      double arg = absAlpha / sqrt2;
      if (arg > 5.) ApproxErf = 1;
      else if (arg < -5.) ApproxErf = -1;
      else ApproxErf = TMath::Erf(arg);
      double leftArea = (1 + ApproxErf) * sqrtPiOver2;
      double rightArea = ( a * 1/TMath::Power(absAlpha - b,n-1)) / (n - 1);
      double area = leftArea + rightArea;
      if( t <= absAlpha ){
        arg = t / sqrt2;
        if(arg > 5.) ApproxErf = 1;
        else if (arg < -5.) ApproxErf = -1;
        else ApproxErf = TMath::Erf(arg);
        return norm * (1 + ApproxErf) * sqrtPiOver2 / area;
      }
      else{
        return norm * (leftArea + a * (1/TMath::Power(t-b,n-1) -
                                       1/TMath::Power(absAlpha - b,n-1)) / (1 - n)) / area;
      }
    }

    void setParam(double m0_, double sigma_, double alpha_, double n_, double norm_)
    {
      m0 = m0_;
      sigma = sigma_;
      alpha = alpha_;
      n = n_;
      norm = norm_;
    }

    double getM0()    { return m0; }
    double getSigma() { return sigma; }
    double getAlpha() { return alpha; }
    double getN()     { return n; }
    double getNorm()  { return norm; }

  private:
    double m0, sigma, alpha, n, norm;
  };


  class CBTurnOnEta : public TrigEff
  {
  public:
    CBTurnOnEta() : barrel(0), endcap(0), boundary(1.479) {}
    CBTurnOnEta(CBTurnOn *barrel_, CBTurnOn *endcap_, double boundary_) :
      barrel(barrel_), endcap(endcap_), boundary(boundary_)
    {}

    double eff(double pt, double eta)
    {
      if(fabs(eta) < boundary)
      {
        return barrel->eff(pt, eta);
      }
      else
      {
        return endcap->eff(pt, eta);
      }
    }

    void setBoundary(double b) { boundary = b; }
    double getBoundary() { return boundary; }

  private:
    CBTurnOn *barrel;
    CBTurnOn *endcap;
    
    double boundary;
  };

  
  
  class CBTurnOn6Eta : public TrigEff
  {
  public:
    CBTurnOn6Eta() : ecM(0),trM(0),bM(0),bP(0),trP(0),ecP(0),b1(-1.2),b2(-0.8),b3(0.),b4(0.8),b5(1.2) {}
    CBTurnOn6Eta(CBTurnOn *ecM_, CBTurnOn *trM_, CBTurnOn *bM_,CBTurnOn *bP_,CBTurnOn *trP_,CBTurnOn *ecP_,double b1_,double b2_,double b3_,double b4_,double b5_) :
      ecM(ecM_),trM(trM_),bM(bM_),bP(bP_),trP(trP_),ecP(ecP_),b1(b1_),b2(b2_),b3(b3_),b4(b4_),b5(b5_) 
      {}
	
   double eff(double pt, double eta)
    {
      if(eta < b1)             return ecM->eff(pt, eta);
      if(eta > b1 && eta < b2) return trM->eff(pt, eta);
      if(eta > b2 && eta < b3) return bM ->eff(pt, eta);
      if(eta > b3 && eta < b4) return ecP->eff(pt, eta);
      if(eta > b4 && eta < b5) return trP->eff(pt, eta);
      return                          ecP->eff(pt, eta);
    }
										   
  private:
   CBTurnOn *ecM;
   CBTurnOn *trM;
   CBTurnOn *bM;
   CBTurnOn *bP;
   CBTurnOn *trP;
   CBTurnOn *ecP;

   double b1,b2,b3,b4,b5;
  };


  class TrigEffMix : public TrigEff
  {
  public:
    double eff(double pt, double eta)
    {
      double sum = 0;
      double mix = 0;

      for(unsigned int i = 0; i < subTrigEffs.size(); i++)
      {
        sum += subTrigEffs[i].first;
        mix += subTrigEffs[i].first * subTrigEffs[i].second->eff(pt, eta);
      }

      mix /= sum;

      return mix;
    }

    void add(double weight, TrigEff *trigEff)
    {
      subTrigEffs.push_back(pair<double, TrigEff *>(weight, trigEff));
    }

    vector<pair<double, TrigEff *> > &getTrigEffs() { return subTrigEffs; }

  private:
    vector<pair<double, TrigEff *> > subTrigEffs;
  };


  class TrigEffRatio : public TrigEff
  {
  public:
    double eff(double pt, double eta)
    {
      return (num->eff(pt, eta) / den->eff(pt, eta));
    }

    void setNumerator(TrigEff *num_) { num = num_; }
    void setDenominator(TrigEff *den_) { den = den_; }

    TrigEff *getNumerator()   { return num; }
    TrigEff *getDenominator() { return den; }

  private:
    TrigEff *num, *den;
  };
}
#endif
