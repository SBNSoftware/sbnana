#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/StandardRecord/Proxy/SRProxy.h"

#include <vector>

namespace ana
{
  class MECSyst: public ISyst
  {
  public:
    MECSyst() : ISyst("mec", "MEC Syst") {}
    void Shift(double sigma, caf::SRSliceProxy* sr, double& weight) const override
    {
      if(sigma < -1) sigma = -1.0;
      if(sr->truth.genie_intcode == 10) weight *= 1.0 + sigma;
    }
  };

  const MECSyst& GetMECSyst()
  {
    static const MECSyst mSyst;
    return mSyst;
  }

  enum class NormSystTerm {
    Constant,
    Sqrt,
    InvSqrt
  };

  enum class NormSystDet {
    All,
    SBND,
    MicroBooNE,
    ICARUS
  };

  template<NormSystTerm term, NormSystDet detector = NormSystDet::All>
  class NormalizationSyst: public ISyst
  {
  public:
    NormalizationSyst(double _uncertainty, std::string name) : 
      ISyst(name, name), uncertainty(_uncertainty) {}
    void Shift(double sigma, caf::SRSliceProxy *sr, double& weight) const override
    {
      auto baseline_cut = [](double baseline){ 
        auto all = (detector == NormSystDet::All);
        auto nd = (detector == NormSystDet::SBND && baseline < 120);
        auto ub = (detector == NormSystDet::MicroBooNE && baseline > 800);
        auto fd = (detector == NormSystDet::ICARUS && baseline > 500);
        return all || nd || ub || fd; 
      };
      if(sr->truth.iscc && abs(sr->truth.pdg) == 14 && baseline_cut(sr->truth.baseline) && !std::isnan(sr->fake_reco.nuE)) {
        double energy_term = 1;
        switch(term) {
        case NormSystTerm::Constant:
          energy_term = 1;
          break;
        case NormSystTerm::Sqrt:
          energy_term = std::sqrt(sr->fake_reco.nuE);
          break;
        case NormSystTerm::InvSqrt:
          energy_term = 1.0/std::sqrt(sr->fake_reco.nuE + 0.1);
          break;
        }
        weight *= (1.0 + uncertainty*sigma*energy_term);
      }
    }   
  private:
    double uncertainty;
  };

  template<NormSystTerm term>
  using NormSystCorr = NormalizationSyst<term, NormSystDet::All>;
  template<NormSystTerm term>
  using NormSystUncorrND = NormalizationSyst<term, NormSystDet::SBND>;
  template<NormSystTerm term>
  using NormSystUncorrUB = NormalizationSyst<term, NormSystDet::MicroBooNE>;
  template<NormSystTerm term>
  using NormSystUncorrFD = NormalizationSyst<term, NormSystDet::ICARUS>;

  extern const NormSystCorr<NormSystTerm::Constant> kNormSyst;
  extern const NormSystCorr<NormSystTerm::Sqrt>     kNormSystSqrt;
  extern const NormSystCorr<NormSystTerm::InvSqrt>  kNormSystInvSqrt;
 
  extern const NormSystUncorrND<NormSystTerm::Constant> kNormSystND;
  extern const NormSystUncorrND<NormSystTerm::Sqrt>     kNormSystSqrtND;
  extern const NormSystUncorrND<NormSystTerm::InvSqrt>  kNormSystInvSqrtND;

  extern const NormSystUncorrUB<NormSystTerm::Constant> kNormSystUB;
  extern const NormSystUncorrUB<NormSystTerm::Sqrt>     kNormSystSqrtUB;
  extern const NormSystUncorrUB<NormSystTerm::InvSqrt>  kNormSystInvSqrtUB;

  extern const NormSystUncorrFD<NormSystTerm::Constant> kNormSystFD;
  extern const NormSystUncorrFD<NormSystTerm::Sqrt>     kNormSystSqrtFD;
  extern const NormSystUncorrFD<NormSystTerm::InvSqrt>  kNormSystInvSqrtFD;

  const std::vector<const ISyst*> GetNormSysts();
}
