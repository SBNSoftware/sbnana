#include "sbnana/CAFAna/Experiment/IExperiment.h"

#include "TFile.h"
#include "TObjString.h"

#include <cassert>
#include <iostream>

// To implement LoadFrom()
#include "sbnana/CAFAna/Experiment/CountingExperiment.h"
#include "sbnana/CAFAna/Experiment/SingleSampleExperiment.h"
#include "sbnana/CAFAna/Experiment/SolarConstraints.h"
#include "sbnana/CAFAna/Experiment/MultiExperiment.h"
#include "sbnana/CAFAna/Experiment/ReactorExperiment.h"

#include "sbnana/CAFAna/Core/Utilities.h"

namespace ana
{
  //----------------------------------------------------------------------
  // Definition to satisfy declaration in Core/LoadFromFile.h
  template<> std::unique_ptr<IExperiment> LoadFrom<IExperiment>(TDirectory* dir, const std::string& name)
  {
    TObjString* ptag = (TObjString*)dir->GetDirectory(name.c_str())->Get("type");
    assert(ptag);

    const TString tag = ptag->GetString();

    if(tag == "CountingExperiment") return CountingExperiment::LoadFrom(dir, name);
    if(tag == "ReactorExperiment") return ReactorExperiment::LoadFrom(dir, name);
    if(tag == "SingleSampleExperiment") return SingleSampleExperiment::LoadFrom(dir, name);
    if(tag == "SolarConstraints") return SolarConstraints::LoadFrom(dir, name);
    if(tag == "MultiExperiment") return MultiExperiment::LoadFrom(dir, name);

    std::cerr << "Unknown Experiment type '" << tag << "'" << std::endl;
    abort();
  }

  //----------------------------------------------------------------------
  void IExperiment::SaveTo(TDirectory* dir, const std::string& name) const
  {
    assert(0 && "Not implemented");
  }
}
