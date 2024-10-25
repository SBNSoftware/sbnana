#include "sbnana/CAFAna/Core/MultiVar.h"

#include <algorithm>
#include <map>
#include <set>

#include <iostream>

namespace ana
{
  /// std::function can wrap a real function, function object, or lambda
  template<class T> _MultiVar<T>::_MultiVar(const std::function<VarFunc_t>& fun)
    : fFunc(fun), fID(fgNextID--)
  {
  }

  //----------------------------------------------------------------------
  /// Helper for \ref MultiVar2D
  template<class T> class MultiVar2DFunc
  {
  public:
    MultiVar2DFunc(const _MultiVar<T>& a, const Binning binsa,
                   const _MultiVar<T>& b, const Binning binsb)
      : fA(a), fBinsA(binsa),
        fB(b), fBinsB(binsb)
    {
    }

    std::vector<double> operator()(const T* sr) const
    {
      std::vector<double> returnVec;

      const std::vector<double> vaVec = fA(sr);
      const std::vector<double> vbVec = fB(sr);

      if(vaVec.size() != vbVec.size())
	{
	  std::cout << "MultiVars need to be same size, these two are size " 
		    << vaVec.size() << " " << vbVec.size() << "respectively." << std::endl;
	  std::abort();
	}
      
      for(unsigned n = 0; n < vaVec.size(); ++n){
        const double va = vaVec.at(n);
        const double vb = vbVec.at(n);
        // Since there are no overflow/underflow bins, check the range
        if(va < fBinsA.Min() || vb < fBinsB.Min()){ returnVec.push_back(-1); continue;}
        if(va > fBinsA.Max() || vb > fBinsB.Max()){ returnVec.push_back(fBinsA.NBins() * fBinsB.NBins()); continue;}

        // FindBin uses root convention, first bin is bin 1, bin 0 is underflow
        const int ia = fBinsA.FindBin(va) - 1;
        const int ib = fBinsB.FindBin(vb) - 1;

        const int i = ia*fBinsB.NBins()+ib;

        returnVec.push_back(i+.5);
      }
      return returnVec;
    }
  protected:
    const _MultiVar<T> fA;
    const Binning fBinsA;
    const _MultiVar<T> fB;
    const Binning fBinsB;
  };

  //----------------------------------------------------------------------
  template<class T> _MultiVar<T>
  MultiVar2D(const _MultiVar<T>& a, const Binning& binsa,
             const _MultiVar<T>& b, const Binning& binsb)
  {
    return _MultiVar<T>(MultiVar2DFunc<T>(a, binsa, b, binsb));
  }

  // explicitly instantiate the template for the types we know we have
  template MultiVar MultiVar2D(const MultiVar&, const Binning&, const MultiVar&, const Binning&);
  template SpillMultiVar MultiVar2D(const SpillMultiVar&, const Binning&, const SpillMultiVar&, const Binning&);
  template TruthMultiVar MultiVar2D(const TruthMultiVar&, const Binning&, const TruthMultiVar&, const Binning&);

  // explicitly instantiate the templates for the types we know we have
  template class _MultiVar<caf::SRSpillProxy>;
  template class _MultiVar<caf::SRSliceProxy>;
  template class _MultiVar<caf::SRTrueInteractionProxy>;

  // Stupid hack to avoid colliding with the IDs of actual Vars. Just count
  // down through negative numbers.
  template<class T> int _MultiVar<T>::fgNextID = -1;

}
