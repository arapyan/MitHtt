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
