#include "cafanacore/Registry.txx"

#include "sbnana/CAFAna/Core/IFitVar.h"
#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Multiverse.h"

namespace ana
{
  // Instantiate the registries we need
  template class Registry<IFitVar>;
  template class Registry<ISyst>;
  template class Registry<Multiverse>;
}
