#pragma once

#include <functional>
#include <set>
#include <string>
#include <vector>

#include "sbnana/CAFAna/Core/Binning.h"

#include "sbnanaobj/StandardRecord/Proxy/FwdDeclare.h"

namespace ana
{
  /// A Var that returns multiple results for each slice. eg the properties of
  /// multiple prongs. All results will be filled into the Spectrum.
  template<class T> class _MultiVar
  {
  public:
    /// The type of the function part of a var
    typedef std::vector<double> (VarFunc_t)(const T* sr);

    /// std::function can wrap a real function, function object, or lambda
    _MultiVar(const std::function<VarFunc_t>& fun);

    /// Allows a variable to be called with double value = myVar(sr) syntax
    std::vector<double> operator()(const T* sr) const
    {
      return fFunc(sr);
    }

    /// Vars with the same definition will have the same ID
    int ID() const {return fID;}

    static int MaxID() {return fgNextID-1;}
  protected:
    _MultiVar(const std::function<VarFunc_t>& fun, int id)
      : fFunc(fun), fID(id)
    {
    }

    std::function<VarFunc_t> fFunc;

    int fID;
    /// The next ID that hasn't yet been assigned
    static int fgNextID;
  };

  typedef _MultiVar<caf::SRSliceProxy> MultiVar;
  typedef _MultiVar<caf::SRSpillProxy> SpillMultiVar;
  typedef _MultiVar<caf::SRTrueInteractionProxy> TruthMultiVar;

  template<class T> _MultiVar<T>
  MultiVar2D(const _MultiVar<T>& a, const Binning& binsa,
             const _MultiVar<T>& b, const Binning& binsb);

} // namespace
