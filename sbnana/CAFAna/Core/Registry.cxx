#include "CAFAna/Core/Registry.txx"

#include "sbnana/CAFAna/Core/IFitVar.h"
#include "sbnana/CAFAna/Core/ISyst.h"

namespace ana
{
  // Instantiate the registries we need
  template class Registry<IFitVar>;
  template class Registry<_ISyst<caf::SRSliceProxy>>;
  template class Registry<ISyst>;
}
