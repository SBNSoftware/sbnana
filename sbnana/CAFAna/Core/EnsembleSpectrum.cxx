#include "sbnana/CAFAna/Core/EnsembleSpectrum.h"

#include "sbnana/CAFAna/Core/EnsembleRatio.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  //----------------------------------------------------------------------
  EnsembleSpectrum::EnsembleSpectrum(SpectrumLoaderBase& loader,
                                     const HistAxis& axis,
                                     const SpillCut& spillcut,
                                     const Cut& cut,
                                     const std::vector<SystShifts>& univ_shifts,
                                     const Var& cv_wei)
    : fNom(loader, axis, spillcut, cut, kNoShift, cv_wei)
  {
    fUnivs.reserve(univ_shifts.size());
    for(const SystShifts& ss: univ_shifts){
      fUnivs.emplace_back(loader, axis, spillcut, cut, ss, cv_wei);
    }
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum::EnsembleSpectrum(SpectrumLoaderBase& loader,
                                     const HistAxis& axis,
                                     const SpillCut& spillcut,
                                     const Cut& cut,
                                     const std::vector<Var>& univ_weis,
                                     const Var& cv_wei)
    : fNom(loader, axis, spillcut, cut, kNoShift, cv_wei)
  {
    fUnivs.reserve(univ_weis.size());
    for(const Var& w: univ_weis){
      fUnivs.emplace_back(loader, axis, spillcut, cut, kNoShift, cv_wei * w);
    }
  }

  //----------------------------------------------------------------------
  void EnsembleSpectrum::Scale(double c)
  {
    fNom.Scale(c);
    for(Spectrum& s: fUnivs) s.Scale(c);
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum& EnsembleSpectrum::operator+=(const EnsembleSpectrum& rhs)
  {
    fNom += rhs.fNom;
    assert(fUnivs.size() == rhs.fUnivs.size());
    for(unsigned int i = 0; i < fUnivs.size(); ++i) fUnivs[i] += rhs.fUnivs[i];
    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum EnsembleSpectrum::operator+(const EnsembleSpectrum& rhs) const
  {
    EnsembleSpectrum ret = *this;
    ret += rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum& EnsembleSpectrum::operator-=(const EnsembleSpectrum& rhs)
  {
    fNom -= rhs.fNom;
    assert(fUnivs.size() == rhs.fUnivs.size());
    for(unsigned int i = 0; i < fUnivs.size(); ++i) fUnivs[i] -= rhs.fUnivs[i];
    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum EnsembleSpectrum::operator-(const EnsembleSpectrum& rhs) const
  {
    EnsembleSpectrum ret = *this;
    ret -= rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum& EnsembleSpectrum::operator*=(const EnsembleRatio& rhs)
  {
    fNom *= rhs.Nominal();
    assert(rhs.NUniverses() == fUnivs.size());
    for(unsigned int i = 0; i < fUnivs.size(); ++i) fUnivs[i] *= rhs.Universe(i);
    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum EnsembleSpectrum::operator*(const EnsembleRatio& rhs) const
  {
    EnsembleSpectrum ret = *this;
    ret *= rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum& EnsembleSpectrum::operator/=(const EnsembleRatio& rhs)
  {
    fNom /= rhs.Nominal();
    assert(rhs.NUniverses() == fUnivs.size());
    for(unsigned int i = 0; i < fUnivs.size(); ++i) fUnivs[i] /= rhs.Universe(i);
    return *this;
  }

  //----------------------------------------------------------------------
  EnsembleSpectrum EnsembleSpectrum::operator/(const EnsembleRatio& rhs) const
  {
    EnsembleSpectrum ret = *this;
    ret /= rhs;
    return ret;
  }

  //----------------------------------------------------------------------
  void EnsembleSpectrum::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    TObjString("EnsembleSpectrum").Write("type");

    fNom.SaveTo(dir->mkdir("nom"));
    for(unsigned int i = 0; i < fUnivs.size(); ++i){
      fUnivs[i].SaveTo(dir->mkdir(("univ"+std::to_string(i)).c_str()));
    }

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<EnsembleSpectrum> EnsembleSpectrum::LoadFrom(TDirectory* dir)
  {
    DontAddDirectory guard;

    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "EnsembleSpectrum");
    delete tag;

    std::unique_ptr<EnsembleSpectrum> ret(new EnsembleSpectrum(*Spectrum::LoadFrom(dir->GetDirectory("nom"))));

    for(unsigned int i = 0; ; ++i){
      TDirectory* d = dir->GetDirectory(("univ"+std::to_string(i)).c_str());
      if(!d) break;
      ret->fUnivs.push_back(*Spectrum::LoadFrom(d));
    }

    return ret;
  }
}
