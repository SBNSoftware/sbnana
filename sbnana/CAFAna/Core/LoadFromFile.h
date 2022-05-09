#pragma once

#include <cassert>
#include <iostream>
#include <memory>
#include <string>

#include "TFile.h"

class TDirectory;

namespace osc
{
  template<class T> class _IOscCalc;
  typedef _IOscCalc<double> IOscCalc;
}

namespace ana
{
  //----------------------------------------------------------------------
  // Most classes are happy to load themselves
  template<class T> std::unique_ptr<T> LoadFrom(TDirectory* dir, const std::string& name)
  {
    return T::LoadFrom(dir, name);
  }

  //----------------------------------------------------------------------
  // But if you're trying to load a base class we need to figure out which
  // derived class is actually in the file and hand off to that. The
  // implementations of these are in the cxx files for the base classes in
  // question.
  class IExtrap;
  template<> std::unique_ptr<IExtrap> LoadFrom<IExtrap>(TDirectory* dir, const std::string& name);
  class IPrediction;
  template<> std::unique_ptr<IPrediction> LoadFrom<IPrediction>(TDirectory* dir, const std::string& name);
  class IExperiment;
  template<> std::unique_ptr<IExperiment> LoadFrom<IExperiment>(TDirectory* dir, const std::string& name);

  // This one is actually implemented in LoadFromFile.cxx to avoid polluting
  // OscLib with CAFAna conventions.
  template<> std::unique_ptr<osc::IOscCalc> LoadFrom<osc::IOscCalc>(TDirectory* dir, const std::string& name);

  //----------------------------------------------------------------------
  // For symmetry
  template<class T> void SaveTo(const T& x, TDirectory* dir, const std::string& name)
  {
    x.SaveTo(dir, name);
  }

  // Also in the cxx, to avoid having to put this logic into OscLib
  template<> void SaveTo(const osc::IOscCalc& x, TDirectory* dir, const std::string& name);

  //----------------------------------------------------------------------
  template<class T> std::unique_ptr<T> LoadFromFile(const std::string& fname,
                                                    const std::string& label)
  {
    TFile fin(fname.c_str());
    assert(!fin.IsZombie());
    TDirectory* dir = fin.GetDirectory(label.c_str());
    if(!dir){
      std::cerr << "Didn't find '" << label << "' in " << fname << std::endl;
      abort();
    }
    return LoadFrom<T>(&fin, label);
  }

  //----------------------------------------------------------------------
  template<class T> void SaveToFile(const T& x,
                                    const std::string& fname,
                                    const std::string& label)
  {
    TFile fout(fname.c_str(), "RECREATE");
    x.SaveTo(&fout, label.c_str());
  }
}
