#include "sbnana/CAFAna/Systs/Systs.h"

namespace ana
{
  const std::vector<const ISyst*>& GetMiscSysts()
  {
    static std::vector<const ISyst*> ret(5);
    if(ret.empty()) {
      ret.push_back(&GetMECSyst());
      ret.push_back(&GetPOTSyst());
      ret.push_back(&GetNormSyst());
      ret.push_back(&GetNormSystND());
      ret.push_back(&GetNormSystFD());
    }
    return ret;
  }
}

