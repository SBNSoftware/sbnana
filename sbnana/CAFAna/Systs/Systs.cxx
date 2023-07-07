#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Systs/Systs.h"

namespace ana
{

const MECSyst& GetMECSyst()
{
  static const MECSyst msyst;
  return msyst;
}

const POTSyst& GetPOTSyst()
{
  static const POTSyst psyst;
  return psyst;
}

}
