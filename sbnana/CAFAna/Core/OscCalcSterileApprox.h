#pragma once

#include "OscLib/IOscCalc.h"

namespace ana
{
  enum class SterileOscAngles {
    kNone = 0,
    kSinSq2ThetaMuMu = 1,
    kSinSq2ThetaMuE  = (1 << 1),
    kSinSq2ThetaEE   = (1 << 2)
  };
  
  class OscCalcSterileApprox: public osc::IOscCalc
  {
  public:
    // if flavAfter == 0, give the active fraction
    virtual double P(int from, int to, double E) override;
    double P_range(int from, int to, double Elo, double Ehi);

    double P_LoverE(int from, int to, double LElo, double LEhi);

    using osc::IOscCalc::P;

    virtual OscCalcSterileApprox* Copy() const override;

    virtual void Print(const std::string& prefix = "") const override;

    void SetDmsq(double d) {fDmsq = d;}
    double GetDmsq() const {return fDmsq;}

    void SetSinSq2ThetaMuMu(double t);
    double GetSinSq2ThetaMuMu() const;

    void SetSinSq2ThetaMuE(double t);
    double GetSinSq2ThetaMuE() const;

    void SetSinSq2ThetaEE(double t);
    double GetSinSq2ThetaEE() const;

    // TODO - potentially remove L in the brave new L/E future
    void SetL(double L) {fL = L;}
    double GetL() const {return fL;}

    TMD5* GetParamsHash() const override;
  protected:
    double PFromDelta(int from, int to, double Delta) const;

    double fDmsq;
    double fSinSq2ThetaMuMu;
    double fSinSq2ThetaMuE;
    double fSinSq2ThetaEE;
    bool   fSinSq2ThetaMuMuSet = false;
    bool   fSinSq2ThetaMuESet  = false;
    bool   fSinSq2ThetaEESet   = false;
    double fL;
  };

  class OscCalcSterileApproxAdjustable: public osc::IOscCalcAdjustable
  {
  public:
    OscCalcSterileApprox calc;

    virtual double P(int from, int to, double E) override
    {
      return calc.P(from, to, E);
    }

    using osc::IOscCalc::P;

    virtual OscCalcSterileApproxAdjustable* Copy() const override
    {
      auto ret = new OscCalcSterileApproxAdjustable;
      auto c = calc.Copy();
      ret->calc = *c;
      delete c;
      return ret;
    }

    virtual void SetL     (double L     ) override {calc.SetL(L);}
    virtual void SetRho   (double rho   ) override {}
    virtual void SetDmsq21(const double& dmsq21) override {}
    virtual void SetDmsq32(const double& dmsq32) override {}
    virtual void SetTh12  (const double& th12  ) override {}
    virtual void SetTh13  (const double& th13  ) override {}
    virtual void SetTh23  (const double& th23  ) override {}
    virtual void SetdCP   (const double& dCP   ) override {}

    TMD5* GetParamsHash() const override {return calc.GetParamsHash();}
  };

  inline SterileOscAngles operator|(const SterileOscAngles a, const SterileOscAngles b)
  {
    int a_int = static_cast<int>(a);
    int b_int = static_cast<int>(b);
    return static_cast<SterileOscAngles>(a_int | b_int);
  }

  inline SterileOscAngles operator&(const SterileOscAngles a, const SterileOscAngles b)
  {
    int a_int = static_cast<int>(a);
    int b_int = static_cast<int>(b);
    return static_cast<SterileOscAngles>(a_int & b_int);
  }

  OscCalcSterileApproxAdjustable* DefaultSterileApproxCalc(SterileOscAngles angles = SterileOscAngles::kSinSq2ThetaMuMu | SterileOscAngles::kSinSq2ThetaMuE);

  const OscCalcSterileApprox* DowncastToSterileApprox(const osc::IOscCalc* calc, bool allowFail = false);
  OscCalcSterileApprox* DowncastToSterileApprox(osc::IOscCalc* calc, bool allowFail = false);
}
