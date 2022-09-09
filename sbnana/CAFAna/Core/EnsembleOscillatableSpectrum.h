#pragma once

#include "sbnana/CAFAna/Core/EnsembleReweightableSpectrum.h"

#include "cafanacore/FwdDeclare.h"
#include "cafanacore/ThreadLocal.h"

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/HistAxis.h"
#include "sbnana/CAFAna/Core/SpectrumLoaderBase.h"
#include "sbnana/CAFAna/Core/StanTypedefs.h"
#include "sbnana/CAFAna/Core/IRecordSource.h"

#include "TMD5.h"

#include <string>

namespace osc
{
  template<class T> class _IOscCalc;
  typedef _IOscCalc<double> IOscCalc;
}

namespace ana
{
  class Binning;

  struct OscCache
  {
    std::unique_ptr<TMD5> hash;
    std::optional<EnsembleSpectrum> spect;
  };

  /// %Spectrum with true energy information, allowing it to be oscillated
  class EnsembleOscillatableSpectrum: public EnsembleReweightableSpectrum
  {
  public:
    EnsembleOscillatableSpectrum(ana::ISliceEnsembleSource& src,
                                 const HistAxis& axis);

    ~EnsembleOscillatableSpectrum();

    /// Copy constructor
    EnsembleOscillatableSpectrum(const EnsembleOscillatableSpectrum& rhs);
    EnsembleOscillatableSpectrum(EnsembleOscillatableSpectrum&& rhs);
    /// Assignment operator
    EnsembleOscillatableSpectrum& operator=(const EnsembleOscillatableSpectrum& rhs);
    EnsembleOscillatableSpectrum& operator=(EnsembleOscillatableSpectrum&& rhs);

    // These under a different name
    EnsembleSpectrum Unoscillated() const {return UnWeighted();}
    EnsembleSpectrum TrueEnergy() const {return WeightingVariable();}

    EnsembleSpectrum Oscillated(osc::IOscCalc* calc, int from, int to) const;

    EnsembleOscillatableSpectrum& operator+=(const EnsembleOscillatableSpectrum& rhs);
    EnsembleOscillatableSpectrum operator+(const EnsembleOscillatableSpectrum& rhs) const;

    EnsembleOscillatableSpectrum& operator-=(const EnsembleOscillatableSpectrum& rhs);
    EnsembleOscillatableSpectrum operator-(const EnsembleOscillatableSpectrum& rhs) const;

    void SaveTo(TDirectory* dir, const std::string& name) const;
    static std::unique_ptr<EnsembleOscillatableSpectrum> LoadFrom(TDirectory* dir, const std::string& name);

  protected:
    EnsembleOscillatableSpectrum(const FitMultiverse* multiverse,
                                 const Eigen::MatrixXd&& mat,
                                 const HistAxis& recoAxis,
                                 double pot, double livetime);

    template<class T> EnsembleSpectrum _Oscillated(osc::_IOscCalc<T>* calc, int from, int to) const;

    mutable ThreadLocal<OscCache> fCache;
  };
}
