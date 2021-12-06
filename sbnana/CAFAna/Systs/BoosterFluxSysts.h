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
  friend const BoosterFluxHadronSyst* GetBoosterFluxHadronSyst(unsigned int);

  BoosterFluxHadronSyst(int i)
      : ISyst(TString::Format("booster_%i", i).Data(),
              TString::Format("Booster flux hadron syst #%i", i).Data()),
        fIdx(i), fScale() {}

  int fIdx;

  mutable TH1* fScale[3];
};

const BoosterFluxHadronSyst* GetBoosterFluxHadronSyst(unsigned int i);

// Because vector<T*> won't automatically convert to vector<U*> even when U
// inherits from V.
struct BoosterFluxHadronSystVector : public std::vector<const BoosterFluxHadronSyst*> {
  operator std::vector<const ISyst*>() {
    return std::vector<const ISyst*>(begin(), end());
  }
};

/// \param N total number of systematics
BoosterFluxHadronSystVector GetBoosterFluxHadronSysts(unsigned int N);

} // namespace ana
