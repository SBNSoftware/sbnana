#include "sbnana/CAFAna/Systs/Systs.h"

namespace ana
{
  const NormSystCorr<NormSystTerm::Constant> kNormSyst(0.02, "normsyst");
  const NormSystCorr<NormSystTerm::Sqrt>     kNormSystSqrt(0.02, "normsystsqrt");
  const NormSystCorr<NormSystTerm::InvSqrt>  kNormSystInvSqrt(0.02, "normsystinvsqrt");

  const NormSystUncorrND<NormSystTerm::Constant> kNormSystND(0.01, "normsystnd");
  const NormSystUncorrND<NormSystTerm::Sqrt>     kNormSystSqrtND(0.01, "normsystsqrtnd");
  const NormSystUncorrND<NormSystTerm::InvSqrt>  kNormSystInvSqrtND(0.01, "normsystinvsqrtnd");

  const NormSystUncorrUB<NormSystTerm::Constant> kNormSystUB(0.01, "normsystub");
  const NormSystUncorrUB<NormSystTerm::Sqrt>     kNormSystSqrtUB(0.01, "normsystsqrtub");
  const NormSystUncorrUB<NormSystTerm::InvSqrt>  kNormSystInvSqrtUB(0.01, "normsystinvsqrtub");

  const NormSystUncorrFD<NormSystTerm::Constant> kNormSystFD(0.01, "normsystfd");
  const NormSystUncorrFD<NormSystTerm::Sqrt>     kNormSystSqrtFD(0.01, "normsystsqrtfd");
  const NormSystUncorrFD<NormSystTerm::InvSqrt>  kNormSystInvSqrtFD(0.01, "normsystinvsqrtfd");

  const std::vector<const ISyst*> GetNormSysts()
  {
    std::vector<const ISyst*> ret;
    ret.push_back(&kNormSyst);
    ret.push_back(&kNormSystSqrt);
    ret.push_back(&kNormSystInvSqrt);
    ret.push_back(&kNormSystND);
    ret.push_back(&kNormSystSqrtND);
    ret.push_back(&kNormSystInvSqrtND);
    ret.push_back(&kNormSystUB);
    ret.push_back(&kNormSystSqrtUB);
    ret.push_back(&kNormSystInvSqrtUB);
    ret.push_back(&kNormSystFD);
    ret.push_back(&kNormSystSqrtFD);
    ret.push_back(&kNormSystInvSqrtFD);
    return ret;
  }
}

