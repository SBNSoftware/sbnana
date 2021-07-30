#include "sbnana/CAFAna/Prediction/PredictionNoExtrap.h"

#include "sbnana/CAFAna/Extrap/IExtrap.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"

#include "sbnana/CAFAna/Core/Loaders.h"
#include "sbnana/CAFAna/Extrap/TrivialExtrap.h"

#include "TDirectory.h"
#include "TObjString.h"


namespace ana
{
  //----------------------------------------------------------------------
  PredictionNoExtrap::PredictionNoExtrap(ISliceSource& srcNonswap,
                                         ISliceSource& srcNue,
                                         ISliceSource& srcNuTau,
					 ISliceSource& srcIntrinsic,
                                         const HistAxis& axis)
    : PredictionExtrap(new TrivialExtrap(srcNonswap, srcNue, srcNuTau, srcIntrinsic, axis))
  {
  }

  //----------------------------------------------------------------------
  PredictionNoExtrap::PredictionNoExtrap(SliceSources& srcs,
                                         const HistAxis& axis)
    : PredictionExtrap(new TrivialExtrap(srcs, axis))
  {
  }

  //----------------------------------------------------------------------
  void PredictionNoExtrap::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();

    TObjString("PredictionNoExtrap").Write("type");

    fExtrap->SaveTo(dir->mkdir("extrap"));

    tmp->cd();
  }


  //----------------------------------------------------------------------
  std::unique_ptr<PredictionNoExtrap> PredictionNoExtrap::LoadFrom(TDirectory* dir)
  {
    assert(dir->GetDirectory("extrap"));
    TrivialExtrap* extrap = ana::LoadFrom<TrivialExtrap>(dir->GetDirectory("extrap")).release();
    return std::unique_ptr<PredictionNoExtrap>(new PredictionNoExtrap(extrap));
  }


  //----------------------------------------------------------------------
  PredictionNoExtrap::~PredictionNoExtrap()
  {
    // We created this in the constructor so it's our responsibility
    delete fExtrap;
  }
}
