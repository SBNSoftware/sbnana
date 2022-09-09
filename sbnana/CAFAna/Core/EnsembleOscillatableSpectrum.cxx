#include "sbnana/CAFAna/Core/EnsembleOscillatableSpectrum.h"

#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/FitMultiverse.h"
#include "sbnana/CAFAna/Core/OscCurve.h"
#include "cafanacore/Ratio.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnana/CAFAna/Analysis/ExpInfo.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "OscLib/IOscCalc.h"

#include "TDirectory.h"
#include "TH2.h"
#include "TMD5.h"
#include "TObjString.h"

#include <cassert>
#include <memory>

namespace ana
{
  // Duplicate here because we can't include Vars.h
  namespace{
    const Var kTrueE([](const caf::SRSliceProxy* slc) -> double {return slc->truth.E;});

    const Var kBaseline([](const caf::SRSliceProxy* slc) -> double {return slc->truth.baseline * 1e-3;}); // m -> km

    const Var kTrueLOverE = kBaseline / kTrueE;

    const Cut kHasNu([](const caf::SRSliceProxy* sr)
                     {
                       return sr->truth.index >= 0;
                     });
  } // end anonymous namespace

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum::
  EnsembleOscillatableSpectrum(ISliceEnsembleSource& src,
                               const HistAxis& axis)
    : EnsembleReweightableSpectrum(src[kHasNu], axis, HistAxis("True L / E (km / GeV)", kTrueLOverEBins, kTrueLOverE))
  {
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum::
  EnsembleOscillatableSpectrum(const FitMultiverse* multiverse,
                               const Eigen::MatrixXd&& mat,
                               const HistAxis& recoAxis,
                               double pot, double livetime)
    : EnsembleReweightableSpectrum(multiverse,
                                   std::move(mat),
                                   recoAxis,
                                   HistAxis("True L / E (km / GeV)",
                                            kTrueLOverEBins,
                                            kTrueLOverE),
                                   pot, livetime)
  {
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum::~EnsembleOscillatableSpectrum()
  {
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum::EnsembleOscillatableSpectrum(const EnsembleOscillatableSpectrum& rhs)
    : EnsembleReweightableSpectrum(rhs)
  {
    if(rhs.fCache->hash){
      fCache->spect = rhs.fCache->spect;
      fCache->hash = std::make_unique<TMD5>(*rhs.fCache->hash);
    }
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum::EnsembleOscillatableSpectrum(EnsembleOscillatableSpectrum&& rhs)
    : EnsembleReweightableSpectrum(rhs)
  {
    if(rhs.fCache->hash){
      fCache->spect = std::move(rhs.fCache->spect);
      fCache->hash = std::move(rhs.fCache->hash);
    }
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum& EnsembleOscillatableSpectrum::operator=(const EnsembleOscillatableSpectrum& rhs)
  {
    if(this == &rhs) return *this;

    EnsembleReweightableSpectrum::operator=(rhs);

    if(rhs.fCache->hash){
      fCache->spect = rhs.fCache->spect;
      fCache->hash = std::make_unique<TMD5>(*rhs.fCache->hash);
    }
    else{
      fCache->hash.reset();
    }

    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum& EnsembleOscillatableSpectrum::operator=(EnsembleOscillatableSpectrum&& rhs)
  {
    if(this == &rhs) return *this;

    EnsembleReweightableSpectrum::operator=(rhs);

    if(rhs.fCache->hash){
      fCache->spect = std::move(rhs.fCache->spect);
      fCache->hash = std::move(rhs.fCache->hash);
    }
    else{
      fCache->hash.reset();
    }

    return *this;
  }

  //----------------------------------------------------------------------
  template<class T> EnsembleSpectrum EnsembleOscillatableSpectrum::
  _Oscillated(osc::_IOscCalc<T>* calc, int from, int to) const
  {
    // POT = 0 implies empty spectrum and oscillated result will also always be
    // empty
    if(fCache->hash && fPOT == 0) return *fCache->spect;

    TMD5* hash = calc->GetParamsHash();
    if(hash && fCache->hash && *hash == *fCache->hash){
      delete hash;
      return *fCache->spect;
    }

    const OscCurve curve(calc, from, to);
    const EnsembleSpectrum ret = WeightedBy(curve);
    if(hash){
      fCache->spect = ret;
      fCache->hash.reset(hash);
    }

    return ret;
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum EnsembleOscillatableSpectrum::Oscillated(osc::IOscCalc* calc,
                                                            int from, int to) const
  {
    return _Oscillated(calc, from, to);
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum& EnsembleOscillatableSpectrum::operator+=(const EnsembleOscillatableSpectrum& rhs)
  {
    EnsembleReweightableSpectrum::operator+=(rhs);

    // invalidate
    fCache->hash.reset(nullptr);

    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum EnsembleOscillatableSpectrum::operator+(const EnsembleOscillatableSpectrum& rhs) const
  {
    EnsembleOscillatableSpectrum ret = *this;
    ret += rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum& EnsembleOscillatableSpectrum::operator-=(const EnsembleOscillatableSpectrum& rhs)
  {
    EnsembleReweightableSpectrum::operator-=(rhs);

    // invalidate
    fCache->hash.reset(nullptr);

    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleOscillatableSpectrum EnsembleOscillatableSpectrum::operator-(const EnsembleOscillatableSpectrum& rhs) const
  {
    EnsembleOscillatableSpectrum ret = *this;
    ret -= rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  void EnsembleOscillatableSpectrum::SaveTo(TDirectory* dir, const std::string& name) const
  {
    _SaveTo(dir, name, "EnsembleOscillatableSpectrum");
  }

  //----------------------------------------------------------------------
  std::unique_ptr<EnsembleOscillatableSpectrum> EnsembleOscillatableSpectrum::LoadFrom(TDirectory* dir, const std::string& name)
  {
    dir = dir->GetDirectory(name.c_str()); // switch to subdir
    assert(dir);

    DontAddDirectory guard;

    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "EnsembleOscillatableSpectrum");
    delete tag;

    TH2D* spect = (TH2D*)dir->Get("hist");
    assert(spect);
    TH1* hPot = (TH1*)dir->Get("pot");
    assert(hPot);
    TH1* hLivetime = (TH1*)dir->Get("livetime");
    assert(hLivetime);

    std::vector<std::string> labels;
    std::vector<Binning> bins;

    for(int i = 0; ; ++i){
      const std::string subname = TString::Format("bins%d", i).Data();
      TDirectory* subdir = dir->GetDirectory(subname.c_str());
      if(!subdir) break;
      delete subdir;
      bins.push_back(*Binning::LoadFrom(dir, subname));
      TObjString* label = (TObjString*)dir->Get(TString::Format("label%d", i));
      labels.push_back(label ? label->GetString().Data() : "");
      delete label;
    }

    delete dir;

    const HistAxis recoAxis(labels, bins);

    const FitMultiverse* multiverse = FitMultiverse::LoadFrom(dir, "multiverse");

    // ROOT histogram storage is row-major, but Eigen is column-major by
    // default
    typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen:: Dynamic, Eigen::RowMajor> MatRowMajor;

    std::unique_ptr<EnsembleOscillatableSpectrum> ret(new EnsembleOscillatableSpectrum(
      multiverse,    
      Eigen::Map<MatRowMajor>(spect->GetArray(),
                              kTrueLOverEBins.NBins()+2,
                              (recoAxis.GetBins1D().NBins() + 2) * multiverse->NUniv()),
      recoAxis,
      hPot->Integral(0, -1),
      hLivetime->Integral(0, -1)));

    delete spect;

    delete hPot;
    delete hLivetime;

    return ret;
  }
}
