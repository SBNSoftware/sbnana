#include "sbnana/CAFAna/Prediction/PredictionIncDirt.h"

#include "sbnana/CAFAna/Core/LoadFromFile.h"

#include "TDirectory.h"
#include "TObjString.h"

namespace ana
{
  // --------------------------------------------------------------------------
  PredictionIncDirt::PredictionIncDirt(ISliceSource& srcNonswap,
                                       ISliceSource& srcNue,
                                       ISliceSource& srcNuTau,
                                       ISliceSource& srcIntrinsic,
                                       ISliceSource& srcDirt,
                                       const HistAxis& axis)
    : fDet(srcNonswap, srcNue, srcNuTau, srcIntrinsic, axis),
      fDirt(srcDirt, kNullSliceSource, kNullSliceSource, kNullSliceSource, axis)
  {
  }

  // --------------------------------------------------------------------------
  PredictionIncDirt::PredictionIncDirt(SliceSources& srcs,
                                       ISliceSource& srcDirt,
                                       const HistAxis& axis)
    : fDet(srcs, axis),
      fDirt(srcDirt, kNullSliceSource, kNullSliceSource, kNullSliceSource, axis)
  {
  }

  // --------------------------------------------------------------------------
  PredictionIncDirt::~PredictionIncDirt()
  {
  }

  // --------------------------------------------------------------------------
  std::unique_ptr<PredictionIncDirt>
  PredictionIncDirt::LoadFrom(TDirectory* dir)
  {
    assert(dir->GetDirectory("det") && dir->GetDirectory("dirt"));

    return std::unique_ptr<PredictionIncDirt>(new PredictionIncDirt(ana::LoadFrom<PredictionNoExtrap>(dir->GetDirectory("det")),
                                                                    ana::LoadFrom<PredictionNoExtrap>(dir->GetDirectory("dirt"))));
  }

  // --------------------------------------------------------------------------
  void PredictionIncDirt::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("PredictionIncDirt").Write("type");

    fDet.SaveTo(dir->mkdir("det"));
    fDirt.SaveTo(dir->mkdir("dirt"));

    tmp->cd();
  }
}
