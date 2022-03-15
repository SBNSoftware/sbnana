#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "sbnana/CAFAna/Core/Progress.h"
#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"
#include "sbnana/CAFAna/Core/SAMProjectSource.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnana/CAFAna/Core/GenieWeightList.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <cassert>
#include <iostream>
#include <cmath>

#include "TFile.h"
#include "TH2.h"
#include "TTree.h"

namespace ana
{
  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::string& wildcard, DataSource src, int max)
    : SpectrumLoaderBase(wildcard, src), max_entries(max)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(const std::vector<std::string>& fnames,
                                 DataSource src, int max)
    : SpectrumLoaderBase(fnames, src), max_entries(max)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader::SpectrumLoader(DataSource src)
    : SpectrumLoaderBase(src)
  {
  }

  //----------------------------------------------------------------------
  SpectrumLoader SpectrumLoader::FromSAMProject(const std::string& proj,
                                                DataSource src,
                                                int fileLimit)
  {
    SpectrumLoader ret;
    ret.fSource = src;
    ret.fWildcard = "project "+proj;
    ret.fFileSource = std::unique_ptr<IFileSource>(new SAMProjectSource(proj, fileLimit));
    return ret;
  }

  //----------------------------------------------------------------------
  SpectrumLoader::~SpectrumLoader()
  {
  }

  struct CompareByID
  {
    bool operator()(const Cut& a, const Cut& b) const
    {
      return a.ID() < b.ID();
    }
  };

  //----------------------------------------------------------------------
  void SpectrumLoader::Go()
  {
    if(fGone){
      std::cerr << "Error: can only call Go() once on a SpectrumLoader" << std::endl;
      abort();
    }
    fGone = true;

    // Find all the unique cuts
    //    std::set<Cut, CompareByID> cuts;
    //    for(auto& shiftdef: fHistDefs)
    //      for(auto& cutdef: shiftdef.second)
    //        cuts.insert(cutdef.first);
    //    for(const Cut& cut: cuts) fAllCuts.push_back(cut);

    //    fLivetimeByCut.resize(fAllCuts.size());
    //    fPOTByCut.resize(fAllCuts.size());


    const int Nfiles = NFiles();

    Progress* prog = 0;

    caf::SRBranchRegistry::clear();

    int fileIdx = -1;
    while(TFile* f = GetNextFile()){
      ++fileIdx;

      if(Nfiles >= 0 && !prog) prog = new Progress(TString::Format("Filling %lu spectra from %d files matching '%s'", fHistDefs.TotalSize(), Nfiles, fWildcard.c_str()).Data());

      HandleFile(f, Nfiles == 1 ? prog : 0);

      if(Nfiles > 1 && prog) prog->SetProgress((fileIdx+1.)/Nfiles);
    } // end for fileIdx

    if(prog){
      prog->Done();
      delete prog;
    }

    StoreExposures(); // also triggers the POT printout

    fHistDefs.RemoveLoader(this);
    fHistDefs.Clear();
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleFile(TFile* f, Progress* prog)
  {
    assert(!f->IsZombie());

    // Test for multi (trees are inside a directory) or single (tree is at top
    // level) tree cases.
    TDirectory* dir = 0;
    TTree* tr = 0;

    TObject* obj = f->Get("recTree");
    assert(obj); // Must exist in one form or the other

    // It might seem like you could use GetObject() to do the type-checking
    // directly, but that method seems to have a memory leak in the case that
    // the types don't match.
    if(obj->ClassName() == std::string("TTree")){
      tr = (TTree*)obj;
    }
    else{
      dir = (TDirectory*)obj;
      tr = (TTree*)dir->Get("rec");
      assert(tr);
    }

    // We try to access this field for every record. It was only added to the
    // files in late 2021, and we don't want to render all earlier files
    // unusable at a stroke. This logic can safely be removed once all extant
    // files have such a field (estimate mid-2022?)
    const bool has_husk = tr->GetLeaf("rec.hdr.husk");

    const caf::CAFType type = caf::GetCAFType(dir, tr);

    long n;
    caf::SRSpillProxy sr(dir, tr, "rec", n, 0);

    //    FloatingExceptionOnNaN fpnan;

    long Nentries = tr->GetEntries();
    if (max_entries != 0 && max_entries < Nentries) Nentries = max_entries;

    for(n = 0; n < Nentries; ++n){
      if(type != caf::kFlatMultiTree) tr->LoadTree(n); // for all single-tree modes

      // If there is no husk field there is no concept of husk events
      if(!has_husk) sr.hdr.husk = false;

      HandleRecord(&sr);

      if(prog) prog->SetProgress(double(n)/Nentries);
    } // end for n
  }

  //----------------------------------------------------------------------
  /// Helper for \ref HandleRecord
  template<class T, class U, class V> class CutVarCache
  {
  public:
    CutVarCache() : fVals(U::MaxID()+1), fValsSet(U::MaxID()+1, false) {}

    inline T Get(const U& var, const V* sr)
    {
      const unsigned int id = var.ID();

      if(fValsSet[id]){
        return fVals[id];
      }
      else{
        const T val = var(sr);
        fVals[id] = val;
        fValsSet[id] = true;
        return val;
      }
    }

  protected:
    // Seems to be faster to do this than [unordered_]map
    std::vector<T> fVals;
    std::vector<bool> fValsSet;
  };

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleRecord(caf::SRSpillProxy* sr)
  {
    if(sr->hdr.first_in_subrun){
      fPOT += sr->hdr.pot;
      // TODO think about if this should be gated behind first_in_file. At the
      // moment I think these will be synonymous. And despite the comment on
      // hdr.pot, I think it may be file-based in practice too.
      fNGenEvt += sr->hdr.ngenevt;
    }

    // This record was only kept as a receptacle for exposure information. It
    // shouldn't be included in any selected spectra.
    if(sr->hdr.husk) return;

    // Do the spill-level spectra first. Keep this very simple because we
    // intend to change it.
    for(auto& spillcutdef: fSpillHistDefs){
      if(!spillcutdef.first(sr)) continue;
      for(auto& spillweidef: spillcutdef.second){
        const double wei = spillweidef.first(sr);
        if(wei == 0) continue;
        for(auto spillvardef: spillweidef.second){
          if(spillvardef.first.IsMulti()){
            for(double val: spillvardef.first.GetMultiVar()(sr)){
              for(Spectrum* s: spillvardef.second.spects) s->Fill(val, wei);
            }
          }
          else{
            const double val = spillvardef.first.GetVar()(sr);
            for(Spectrum* s: spillvardef.second.spects) s->Fill(val, wei);
          }
        }
      }
    }


    for(auto& spillcutdef: fHistDefs){
      const SpillCut& spillcut = spillcutdef.first;

      const bool spillpass = spillcut(sr); // nomSpillCutCache.Get(spillcut, sr);
      // Cut failed, skip all the histograms that depended on it
      if(!spillpass) continue;

      for(caf::SRSliceProxy& slc: sr->slc){

        // Some shifts only adjust the weight, so they're effectively nominal,
        // but aren't grouped with the other nominal histograms. Keep track of
        // the results for nominals in these caches to speed those systs up.
        CutVarCache<bool, Cut, caf::SRSliceProxy> nomCutCache;
        CutVarCache<double, Var, caf::SRSliceProxy> nomWeiCache;
        CutVarCache<double, Var, caf::SRSliceProxy> nomVarCache;

        for(auto& shiftdef: spillcutdef.second){
          const SystShifts& shift = shiftdef.first;

          // Need to provide a clean slate for each new set of systematic
          // shifts to work from. Copying the whole StandardRecord is pretty
          // expensive, so modify it in place and revert it afterwards.

          caf::SRProxySystController::BeginTransaction();

          bool shifted = false;

          double systWeight = 1;
          // Can special-case nominal to not pay cost of Shift()
          if(!shift.IsNominal()){
            shift.Shift(&slc, systWeight);
            // If there were only weighting systs applied then the cached
            // nominal values are still valid.
            shifted = caf::SRProxySystController::AnyShifted();
          }

          for(auto& cutdef: shiftdef.second){
            const Cut& cut = cutdef.first;

            const bool pass = shifted ? cut(&slc) : nomCutCache.Get(cut, &slc);
            // Cut failed, skip all the histograms that depended on it
            if(!pass) continue;

            for(auto& weidef: cutdef.second){
              const Var& weivar = weidef.first;

              double wei = shifted ? weivar(&slc) : nomWeiCache.Get(weivar, &slc);

              wei *= systWeight;
              if(wei == 0) continue;

              for(auto& vardef: weidef.second){
                if(vardef.first.IsMulti()){
                  for(double val: vardef.first.GetMultiVar()(&slc)){
                    for(Spectrum* s: vardef.second.spects)
                      s->Fill(val, wei);
                  }
                  continue;
                }

                const Var& var = vardef.first.GetVar();

                const double val = shifted ? var(&slc) : nomVarCache.Get(var, &slc);

                if(std::isnan(val) || std::isinf(val)){
                  std::cerr << "Warning: Bad value: " << val
                            << " returned from a Var. The input variable(s) could "
                            << "be NaN in the CAF, or perhaps your "
                            << "Var code computed 0/0?";
                  std::cout << " Not filling into this histogram for this slice." << std::endl;
                  continue;
                }

                for(Spectrum* s: vardef.second.spects) s->Fill(val, wei);

                for(ReweightableSpectrum* rw: vardef.second.rwSpects){
                  const double yval = rw->ReweightVar()(&slc);

                  if(std::isnan(yval) || std::isinf(yval)){
                    std::cerr << "Warning: Bad value: " << yval
                              << " for reweighting Var";
                    std::cout << ". Not filling into histogram." << std::endl;
                    continue;
                  }

                  rw->fHist->Fill(val, yval, wei);
                } // end for rw
              } // end for vardef
            } // end for weidef
          } // end for cutdef

          // Return StandardRecord to its unshifted form ready for the next
          // histogram.
          caf::SRProxySystController::Rollback();
        } // end for shiftdef
      } // end for slc
    } // end for spillcutdef
  }

  //----------------------------------------------------------------------
  void SpectrumLoader::StoreExposures()
  {
    if(fabs(fPOT - fPOTFromHist)/std::min(fPOT, fPOTFromHist) > 0.001){
      std::cout << fPOT << " POT from hdr differs from " << fPOTFromHist << " POT from the TotalPOT histogram!" << std::endl;
      abort();
    }

    std::cout << fPOT << " POT over " << fNGenEvt << " readouts" << std::endl;

    for(auto& shiftdef: fHistDefs){
      for(auto& spillcutdef: shiftdef.second){
        for(auto& cutdef: spillcutdef.second){
          for(auto& weidef: cutdef.second){
            for(auto& vardef: weidef.second){
              for(Spectrum* s: vardef.second.spects){
                s->fPOT += fPOT;
                s->fLivetime += fNGenEvt;
              }
              for(ReweightableSpectrum* rw: vardef.second.rwSpects){
                rw->fPOT += fPOT;
                rw->fLivetime += fNGenEvt;
              }
            }
          }
        }
      }
    }


    for(auto& spillcutdef: fSpillHistDefs){
      for(auto& spillweidef: spillcutdef.second){
        for(auto spillvardef: spillweidef.second){
          for(Spectrum* s: spillvardef.second.spects){
            s->fPOT += fPOT;
            s->fLivetime += fNGenEvt;
          }
        }
      }
    }

  }
} // namespace
