#pragma once

#include "sbnana/CAFAna/Core/IFitVar.h"
#include "TMath.h"

namespace ana
{
  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{32} \f$
  class FitDmSq32Sterile: public IFitVar
  {
  public:
    FitDmSq32Sterile() : IFitVar("dmsq32", "#Deltam^{2}_{32} (eV^{2})") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };

  /// \f$ \Delta m^2_{32} \f$
  const FitDmSq32Sterile kFitDmSq32Sterile;
  
  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{41} \f$
  class FitDmSq41Sterile: public IFitVar
  {
  public:
    FitDmSq41Sterile() : IFitVar("dmsq41", "#Deltam^{2}_{41} (eV^{2})") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };

  /// \f$ \Delta m^2_{41} \f$
  const FitDmSq41Sterile kFitDmSq41Sterile;

  //----------------------------------------------------------------------

  /// \f$ \Delta m^2_{43} \f$
  class FitDmSq43Sterile: public IFitVar
  {
  public:
    FitDmSq43Sterile() : IFitVar("dmsq43", "#Deltam^{2}_{43} (eV^{2})") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };
  
  /// \f$ \Delta m^2_{43} \f$
  const FitDmSq43Sterile kFitDmSq43Sterile;
  
  //----------------------------------------------------------------------

  /// \f$ \delta_{13}/\pi \f$
  class FitDelta13InPiUnitsSterile: public IFitVar
  {
  public:
    FitDelta13InPiUnitsSterile() : IFitVar("delta13(pi)", "#delta_{13} / #pi") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };

  /// \f$ \delta_{CP}/\pi \f$
  const FitDelta13InPiUnitsSterile kFitDelta13InPiUnitsSterile;

  //----------------------------------------------------------------------

  /// \f$ \delta_{13}/\pi \f$
  class FitDelta14InPiUnitsSterile: public IFitVar
  {
  public:
    FitDelta14InPiUnitsSterile() : IFitVar("delta14(pi)", "#delta_{14} / #pi") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };

  /// \f$ \delta_{14}/\pi \f$
  const FitDelta14InPiUnitsSterile kFitDelta14InPiUnitsSterile;

  //----------------------------------------------------------------------

  /// \f$ \delta_{24}/\pi \f$
  class FitDelta24InPiUnitsSterile: public IFitVar
  {
  public:
    FitDelta24InPiUnitsSterile() : IFitVar("delta24(pi)", "#delta_{24} / #pi") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;
  };

  /// \f$ \delta_{24}/\pi \f$
  const FitDelta24InPiUnitsSterile kFitDelta24InPiUnitsSterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{13} \f$
  class FitTheta13Sterile: public IConstrainedFitVar
  {
  public:
    FitTheta13Sterile() : IConstrainedFitVar("th13", "#theta_{13}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/2;}
  };

  /// \f$ \theta_{13} \f$
  const FitTheta13Sterile kFitTheta13Sterile;

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{13} \f$
  class FitSinSqTheta13Sterile: public IConstrainedFitVar
  {
  public:
    FitSinSqTheta13Sterile() : IConstrainedFitVar("ssth13", "sin^{2}#theta_{13}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{13} \f$
  const FitSinSqTheta13Sterile kFitSinSqTheta13Sterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{23} \f$
  class FitTheta23Sterile: public IConstrainedFitVar
  {
  public:
    FitTheta23Sterile() : IConstrainedFitVar("th23", "#theta_{23}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/2;}
  };

  /// \f$ \theta_{23} \f$
  const FitTheta23Sterile kFitTheta23Sterile;

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{23} \f$
  class FitSinSqTheta23Sterile: public IConstrainedFitVar
  {
  public:
    FitSinSqTheta23Sterile() : IConstrainedFitVar("ssth23", "sin^{2}#theta_{23}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{23} \f$
  const FitSinSqTheta23Sterile kFitSinSqTheta23Sterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{14} \f$
  class FitTheta14Sterile: public IConstrainedFitVar
  {
  public:
    FitTheta14Sterile() : IConstrainedFitVar("th14", "#theta_{14}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/2;}
  };

  /// \f$ \theta_{14} \f$
  const FitTheta14Sterile kFitTheta14Sterile;

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{14} \f$
  class FitSinSqTheta14Sterile: public IConstrainedFitVar
  {
  public:
    FitSinSqTheta14Sterile() : IConstrainedFitVar("ssth14", "sin^{2}#theta_{14}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{14} \f$
  const FitSinSqTheta14Sterile kFitSinSqTheta14Sterile;


 //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{14} \f$
  class FitSinSq2Theta14Sterile: public IConstrainedFitVar
  {
  public:
    FitSinSq2Theta14Sterile() : IConstrainedFitVar("ss2th14", "sin^{2}2#theta_{14}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{14} \f$
  const FitSinSq2Theta14Sterile kFitSinSq2Theta14Sterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{24} \f$
  class FitTheta24Sterile: public IConstrainedFitVar
  {
  public:
    FitTheta24Sterile() : IConstrainedFitVar("th24", "#theta_{24}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/4;}
  };

  /// \f$ \theta_{24} \f$
  const FitTheta24Sterile kFitTheta24Sterile;

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{24} \f$
  class FitSinSqTheta24Sterile: public IConstrainedFitVar
  {
  public:
    FitSinSqTheta24Sterile() : IConstrainedFitVar("ssth24", "sin^{2}#theta_{24}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{24} \f$
  const FitSinSqTheta24Sterile kFitSinSqTheta24Sterile;

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{24} \f$
  class FitSinSq2Theta24Sterile: public IConstrainedFitVar
  {
  public:
    FitSinSq2Theta24Sterile() : IConstrainedFitVar("ss2th24", "sin^{2}2#theta_{24}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{24} \f$
  const FitSinSq2Theta24Sterile kFitSinSq2Theta24Sterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{34} \f$
  class FitTheta34Sterile: public IConstrainedFitVar
  {
  public:
    FitTheta34Sterile() : IConstrainedFitVar("th34", "#theta_{34}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return TMath::Pi()/2;}
  };

  /// \f$ \theta_{34} \f$
  const FitTheta34Sterile kFitTheta34Sterile;

  //----------------------------------------------------------------------

  /// \f$ \sin^2\theta_{34} \f$
  class FitSinSqTheta34Sterile: public IConstrainedFitVar
  {
  public:
    FitSinSqTheta34Sterile() : IConstrainedFitVar("ssth34", "sin^{2}#theta_{34}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^2\theta_{34} \f$
  const FitSinSqTheta34Sterile kFitSinSqTheta34Sterile;

  //----------------------------------------------------------------------

  /// \f$ \sin^22\theta_{34} \f$
  class FitSinSq2Theta34Sterile: public IConstrainedFitVar
  {
  public:
    FitSinSq2Theta34Sterile() : IConstrainedFitVar("ss2th34", "sin^{2}2#theta_{34}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 1;}
  };

  /// \f$ \sin^22\theta_{34} \f$
  const FitSinSq2Theta34Sterile kFitSinSq2Theta34Sterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{13} \f$
  class FitTheta13InDegreesSterile: public IConstrainedFitVar
  {
  public:
    FitTheta13InDegreesSterile() : IConstrainedFitVar("th13(degrees)", "#theta_{13}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 90;}
  };

  /// \f$ \theta_{13} \f$
  const FitTheta13InDegreesSterile kFitTheta13InDegreesSterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{23} \f$
  class FitTheta23InDegreesSterile: public IConstrainedFitVar
  {
  public:
    FitTheta23InDegreesSterile() : IConstrainedFitVar("th23(degrees)", "#theta_{23}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 90;}
  };

  /// \f$ \theta_{23} \f$
  const FitTheta23InDegreesSterile kFitTheta23InDegreesSterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{14} \f$
  class FitTheta14InDegreesSterile: public IConstrainedFitVar
  {
  public:
    FitTheta14InDegreesSterile() : IConstrainedFitVar("th14(degrees)", "#theta_{14}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 90;}
  };

  /// \f$ \theta_{14} \f$
  const FitTheta14InDegreesSterile kFitTheta14InDegreesSterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{24} \f$
  class FitTheta24InDegreesSterile: public IConstrainedFitVar
  {
  public:
    FitTheta24InDegreesSterile() : IConstrainedFitVar("th24(degrees)", "#theta_{24}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 45;}
  };

  /// \f$ \theta_{24} \f$
  const FitTheta24InDegreesSterile kFitTheta24InDegreesSterile;

  //----------------------------------------------------------------------

  /// \f$ \theta_{34} \f$
  class FitTheta34InDegreesSterile: public IConstrainedFitVar
  {
  public:
    FitTheta34InDegreesSterile() : IConstrainedFitVar("th34(degrees)", "#theta_{34}") {}

    virtual double GetValue(const osc::IOscCalcAdjustable* osc) const;
    virtual void SetValue(osc::IOscCalcAdjustable* osc, double val) const;

    virtual double LowLimit() const {return 0;}
    virtual double HighLimit() const {return 90;}
  };

  /// \f$ \theta_{34} \f$
  const FitTheta34InDegreesSterile kFitTheta34InDegreesSterile;

  //----------------------------------------------------------------------


} // namespace
