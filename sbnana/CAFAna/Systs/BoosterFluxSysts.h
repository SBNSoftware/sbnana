#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"

#include "TString.h"

class TH1;

namespace ana {

class BoosterFluxHadronSyst : public ISyst
{
public:
  virtual ~BoosterFluxHadronSyst();

  void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override;

protected:
  friend const BoosterFluxHadronSyst* GetBoosterFluxHadronSyst(const std::string &particle, unsigned int);

  BoosterFluxHadronSyst(const std::string &particle, int i)
      : ISyst((particle + std::to_string(i)).c_str(),
              (particle + " flux systematic #" + std::to_string(i)).c_str()),
        fIdx(i), particle(particle), fScale() {}

  int fIdx;
  std::string particle;

  mutable TH1* fScale[3];
};

const BoosterFluxHadronSyst* GetBoosterFluxHadronSyst(const std::string &particle, unsigned int i);

// Because vector<T*> won't automatically convert to vector<U*> even when U
// inherits from V.
struct BoosterFluxHadronSystVector : public std::vector<const BoosterFluxHadronSyst*> {
  operator std::vector<const ISyst*>() {
    return std::vector<const ISyst*>(begin(), end());
  }
};

/// \param N total number of systematics
BoosterFluxHadronSystVector GetBoosterFluxHadronSysts(const std::string &particle, unsigned int N);

} // namespace ana
