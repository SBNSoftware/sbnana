#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "sbnana/CAFAna/Core/Progress.h"
#include "sbnana/CAFAna/Core/ReweightableSpectrum.h"
#include "sbnana/CAFAna/Core/SAMProjectSource.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/Tree.h"

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

    fSpillHistDefs.RemoveLoader(this);
    fSpillHistDefs.Clear();

    fTruthHistDefs.RemoveLoader(this);
    fTruthHistDefs.Clear();

    fTruthHistWithCutDefs.RemoveLoader(this);
    fTruthHistWithCutDefs.Clear();

  }

  //----------------------------------------------------------------------
  void SpectrumLoader::HandleFile(TFile* f, Progress* prog)
  {
    assert(!f->IsZombie());

    TTree* tr = (TTree*)f->Get("recTree");
    assert(tr);

    // We try to access this field for every record. It was only added to the
    // files in late 2021, and we don't want to render all earlier files
    // unusable at a stroke. This logic can safely be removed once all extant
    // files have such a field (estimate mid-2022?)
    const bool has_husk = tr->GetLeaf("rec.hdr.husk");

    caf::SRSpillProxy sr(tr, "rec");

    //    FloatingExceptionOnNaN fpnan;

    long Nentries = tr->GetEntries();
    if (max_entries != 0 && max_entries < Nentries) Nentries = max_entries;

    for(long n = 0; n < Nentries; ++n){
      tr->LoadTree(n);

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
      if(sr->hdr.ismc) fNReadouts += sr->hdr.ngenevt;
    }

    if(!sr->hdr.ismc){
      const int nbnb  = sr->hdr.bnbinfo.size();
      const int nnumi = sr->hdr.numiinfo.size();
      if(nbnb > 0 && nnumi > 0){
        std::cout << "SpectrumLoader: nonzero number of both BNB (" << nbnb
                  << ") and NuMI (" << nnumi << ") triggers. I'm confused"
                  << std::endl;
        abort();
      }

      fNReadouts += nbnb + nnumi;
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

    // TruthVar without (Slice)Cut

    for(auto& spillcutdef: fTruthHistDefs){
      const SpillCut& spillcut = spillcutdef.first;

      const bool spillpass = spillcut(sr); // nomSpillCutCache.Get(spillcut, sr);
      // Cut failed, skip all the histograms that depended on it
      if(!spillpass) continue;

      // now start the nu loop
      for(caf::SRTrueInteractionProxy& nu: sr->mc.nu){

        // Some shifts only adjust the weight, so they're effectively nominal,
        // but aren't grouped with the other nominal histograms. Keep track of
        // the results for nominals in these caches to speed those systs up.
        CutVarCache<bool, TruthCut, caf::SRTrueInteractionProxy> nomTruthCutCache;
        CutVarCache<double, TruthVar, caf::SRTrueInteractionProxy> nomTruthWeiCache;
        CutVarCache<double, TruthVar, caf::SRTrueInteractionProxy> nomTruthVarCache;

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
            shift.Shift(&nu, systWeight);
            // If there were only weighting systs applied then the cached
            // nominal values are still valid.
            shifted = caf::SRProxySystController::AnyShifted();
          }

          for(auto& truthcutdef: shiftdef.second){

            const TruthCut& truthcut = truthcutdef.first;

            const bool truthpass = shifted ? truthcut(&nu) : nomTruthCutCache.Get(truthcut, &nu);

            // TruthCut failed, skip all the histograms that depended on it
            if(!truthpass) continue;

            for(auto& truthweidef: truthcutdef.second){

              const TruthVar& truthweivar = truthweidef.first;

              double truthwei = shifted ? truthweivar(&nu) : nomTruthWeiCache.Get(truthweivar, &nu);

              truthwei *= systWeight;
              if(truthwei == 0) continue;

              for(auto& truthvardef: truthweidef.second){

                // if TruthMultiVar
                if(truthvardef.first.IsMulti()){
                  for(double truthval: truthvardef.first.GetMultiVar()(&nu)){
                    for(Spectrum* s: truthvardef.second.spects)
                      s->Fill(truthval, truthwei);
                  }
                }
                // if TruthVar
                else{

                  const TruthVar& truthvar = truthvardef.first.GetVar();
                  const double truthval = shifted ? truthvar(&nu) : nomTruthVarCache.Get(truthvar, &nu);

                  if(std::isnan(truthval) || std::isinf(truthval)){
                    std::cerr << "Warning: Bad value: " << truthval
                              << " returned from a TruthVar. The input variable(s) could "
                              << "be NaN in the CAF, or perhaps your "
                              << "Var code computed 0/0?";
                    std::cout << " Not filling into this histogram for this slice." << std::endl;
                    continue;
                  }

                  for(Spectrum* s: truthvardef.second.spects) s->Fill(truthval, truthwei);

                }

              } // end for truthvardef



            } // end for truthweidef

          } // end for truthcutdef

          // Return StandardRecord to its unshifted form ready for the next
          // histogram.
          caf::SRProxySystController::Rollback();

        } // end for shiftdef

      } // end for nu loop

    } // end for spillcutdef


    // TruthVar with (Slice)Cut by truth-matching

    for(auto& spillcutdef: fTruthHistWithCutDefs){
      const SpillCut& spillcut = spillcutdef.first;

      const bool spillpass = spillcut(sr); // nomSpillCutCache.Get(spillcut, sr);
      // Cut failed, skip all the histograms that depended on it
      if(!spillpass) continue;

      for(auto& cutdef: spillcutdef.second){

        const Cut& cut = cutdef.first;

        // now start the nu loop
        for(caf::SRTrueInteractionProxy& nu: sr->mc.nu){

          // Some shifts only adjust the weight, so they're effectively nominal,
          // but aren't grouped with the other nominal histograms. Keep track of
          // the results for nominals in these caches to speed those systs up.
          CutVarCache<bool, TruthCut, caf::SRTrueInteractionProxy> nomTruthCutCache;
          CutVarCache<double, TruthVar, caf::SRTrueInteractionProxy> nomTruthWeiCache;
          CutVarCache<double, TruthVar, caf::SRTrueInteractionProxy> nomTruthVarCache;

          for(auto& shiftdef: cutdef.second){
            const SystShifts& shift = shiftdef.first;

            // Loop over reco slices, and check if the truth-matched slice pass the (Slice)Cut
            // We have to shift Slice, and rollback for the actual neutrino shifts
            // NB) Here, the "reweighting" shifts is not considered but only the lateral shifts on reco selection
            bool HasMatchedSlicePassCut = false;
            for(caf::SRSliceProxy& slc: sr->slc){
              caf::SRProxySystController::BeginTransaction();
              double dummy_systWeight = 1;
              if(!shift.IsNominal()){
                shift.Shift(&slc, dummy_systWeight);
              }
              if ( slc.truth.index < 0 ) continue;
              else if ( slc.truth.index != nu.index ) continue;
              if( cut(&slc) ){
                HasMatchedSlicePassCut = true;
                break;
              }
              caf::SRProxySystController::Rollback();
            }
            if(!HasMatchedSlicePassCut) continue;


            // Need to provide a clean slate for each new set of systematic
            // shifts to work from. Copying the whole StandardRecord is pretty
            // expensive, so modify it in place and revert it afterwards.

            caf::SRProxySystController::BeginTransaction();

            bool shifted = false;
            double systWeight = 1;
            // Can special-case nominal to not pay cost of Shift()
            if(!shift.IsNominal()){
              shift.Shift(&nu, systWeight);
              // If there were only weighting systs applied then the cached
              // nominal values are still valid.
              shifted = caf::SRProxySystController::AnyShifted();
            }

            for(auto& truthcutdef: shiftdef.second){

              const TruthCut& truthcut = truthcutdef.first;
              const bool truthpass = shifted ? truthcut(&nu) : nomTruthCutCache.Get(truthcut, &nu);

              // TruthCut failed, skip all the histograms that depended on it
              if(!truthpass) continue;

              for(auto& truthweidef: truthcutdef.second){

                const TruthVar& truthweivar = truthweidef.first;

                double truthwei = shifted ? truthweivar(&nu) : nomTruthWeiCache.Get(truthweivar, &nu);

                truthwei *= systWeight;
                if(truthwei == 0) continue;

                for(auto& truthvardef: truthweidef.second){

                  // if TruthMultiVar
                  if(truthvardef.first.IsMulti()){
                    for(double truthval: truthvardef.first.GetMultiVar()(&nu)){
                      for(Spectrum* s: truthvardef.second.spects)
                        s->Fill(truthval, truthwei);
                    }
                  }
                  // if TruthVar
                  else{

                    const TruthVar& truthvar = truthvardef.first.GetVar();
                    const double truthval = shifted ? truthvar(&nu) : nomTruthVarCache.Get(truthvar, &nu);

                    if(std::isnan(truthval) || std::isinf(truthval)){
                      std::cerr << "Warning: Bad value: " << truthval
                                << " returned from a TruthVar. The input variable(s) could "
                                << "be NaN in the CAF, or perhaps your "
                                << "Var code computed 0/0?";
                      std::cout << " Not filling into this histogram for this slice." << std::endl;
                      continue;
                    }

                    for(Spectrum* s: truthvardef.second.spects) s->Fill(truthval, truthwei);

                  }

                } // end for truthvardef

              } // end for truthweidef

            } // end for truthcutdef

            // Return StandardRecord to its unshifted form ready for the next
            // histogram.
            caf::SRProxySystController::Rollback();

          } // end for shiftdef 

        } // end for nu loop

      } // end for cutdef

    } // end for spillcutdef



    // Trees
    //unsigned int idxSpillCut = 0; // testing
    for ( auto& [spillcut, shiftmap] : fTreeDefs ) {
      const bool spillpass = spillcut(sr);
      // Cut failed, skip all the histograms that depend on it
      if(!spillpass) continue;

      unsigned int idxSlice = 0; // in case we want to save the slice number to the tree
      for( caf::SRSliceProxy& slc: sr->slc ) {
        // Some shifts only adjust the weight, so they're effectively nominal,
        // but aren't grouped with the other nominal histograms. Keep track of
        // the results for nominals in these caches to speed those systs up.
        CutVarCache<bool, Cut, caf::SRSliceProxy> nomCutCache;
        CutVarCache<double, Var, caf::SRSliceProxy> nomVarCache;

        //unsigned int idxShift = 0; // testing
        for ( auto& [shift, cutmap] : shiftmap ) {
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

          //unsigned int idxCut = 0; // testing
          for ( auto& [cut, treemap] : cutmap ) {
            const bool pass = shifted ? cut(&slc) : nomCutCache.Get(cut, &slc);
            // Cut failed, skip all the histograms that depended on it
            if(!pass) continue;

            //unsigned int idxTree = 0;
            for ( std::map<Tree*, std::map<VarOrMultiVar, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
              //unsigned int idxVar = 0;
              std::map<std::string, std::vector<double>> recordVals;
              unsigned int numEntries=0;
              for ( auto& [varormulti, varname] : treemapIt->second ) {
                //std::cout << "SpillCut " << idxSpillCut << " Slice " << idxSlice << " Shift " << idxShift << " Cut " << idxCut << " Tree " << idxTree << " Var " << idxVar << std::endl;
                if(varormulti.IsMulti()){
                  auto const& vals = varormulti.GetMultiVar()(&slc);
                  for(double val: vals) recordVals[varname].push_back(val);
                  if (numEntries==0)    numEntries = vals.size();
                  continue;
                }

                const Var& var = varormulti.GetVar();
                const double val = shifted ? var(&slc) : nomVarCache.Get(var, &slc);

                //std::cout << "    VAL = " << val << std::endl;

                if(std::isnan(val) || std::isinf(val)){
                  std::cerr << "Warning: Bad value: " << val
                            << " returned from a Var. The input variable(s) could "
                            << "be NaN in the CAF, or perhaps your "
                            << "Var code computed 0/0?";
                  std::cout << " Still filling into the ''branch'' for this slice." << std::endl;
                }

                recordVals[varname].push_back(val);
                if( numEntries==0 ) numEntries = 1;
                //idxVar+=1;
              } // end for var/varname
              // If fSaveRunSubrunEvt then fill these entries...
              if ( treemapIt->first->SaveRunSubEvent() || treemapIt->first->SaveSliceNum() ) {
                for ( unsigned int idxRun=0; idxRun<numEntries; ++idxRun ) {
                  if ( treemapIt->first->SaveRunSubEvent() ) {
                    recordVals["Run/i"].push_back( sr->hdr.run );
                    recordVals["Subrun/i"].push_back( sr->hdr.subrun );
                    recordVals["Evt/i"].push_back( sr->hdr.evt );
                  }
                  if ( treemapIt->first->SaveSliceNum() ) {
                    recordVals["Slice/i"].push_back( idxSlice );
                  }
                }
              }
              treemapIt->first->UpdateEntries(recordVals);
              //idxTree+=1;
            } // end for tree
            //idxCut+=1;
          } // end for cut

          // Return StandardRecord to its unshifted form ready for the next
          // histogram.
          caf::SRProxySystController::Rollback();

          //idxShift+=1;
        } // end for shift
        idxSlice+=1;
      } // end for slice
      //idxSpillCut+=1;
    } // end for spillcut

    // Spill Trees
    for ( auto& [spillcut, treemap] : fSpillTreeDefs ) {
      const bool spillpass = spillcut(sr);
      // Cut failed, skip all the histograms that depend on it
      if(!spillpass) continue;

      for ( std::map<Tree*, std::map<SpillVarOrMultiVar, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
        std::map<std::string, std::vector<double>> recordVals;
        unsigned int numEntries = 0;
        for ( auto& [varormulti, varname] : treemapIt->second ) {
          if(varormulti.IsMulti()){
            auto const& vals = varormulti.GetMultiVar()(sr);
            for(double val: vals) recordVals[varname].push_back(val);
            if ( numEntries==0 )  numEntries = vals.size();
            continue;
          }

          const double val = varormulti.GetVar()(sr);

          if(std::isnan(val) || std::isinf(val)){
            std::cerr << "Warning: Bad value: " << val
                      << " returned from a Var. The input variable(s) could "
                      << "be NaN in the CAF, or perhaps your "
                      << "Var code computed 0/0?";
            std::cout << " Still filling into the ''branch'' for this slice." << std::endl;
          }

          recordVals[varname].push_back(val);
          if ( numEntries==0 ) numEntries=1;
        } // end for var/varname
        if ( treemapIt->first->SaveRunSubEvent() ) {
          for ( unsigned int idxRun=0; idxRun<numEntries; ++idxRun ) {
            recordVals["Run/i"].push_back( sr->hdr.run );
            recordVals["Subrun/i"].push_back( sr->hdr.subrun );
            recordVals["Evt/i"].push_back( sr->hdr.evt );
          }
        }
        treemapIt->first->UpdateEntries(recordVals);
      } // end for tree
    } // end for spillcut

    // TruthTrees
    //unsigned int idxSpillCut = 0; // testing
    for ( auto& [spillcut, shiftmap] : fTruthTreeDefs ) {
      const bool spillpass = spillcut(sr);

      for(caf::SRTrueInteractionProxy& nu: sr->mc.nu){
        // Some shifts only adjust the weight, so they're effectively nominal,
        // but aren't grouped with the other nominal histograms. Keep track of
        // the results for nominals in these caches to speed those systs up.
        CutVarCache<bool, TruthCut, caf::SRTrueInteractionProxy> nomTruthCutCache;
        CutVarCache<double, TruthVar, caf::SRTrueInteractionProxy> nomTruthVarCache;

        //unsigned int idxShift = 0; // testing
        for ( auto& [shift, truthcutmap] : shiftmap ) {
          // Need to provide a clean slate for each new set of systematic
          // shifts to work from. Copying the whole StandardRecord is pretty
          // expensive, so modify it in place and revert it afterwards.
          caf::SRProxySystController::BeginTransaction();

          bool shifted = false;

          double systWeight = 1;
          // Can special-case nominal to not pay cost of Shift()
          if(!shift.IsNominal()){
            shift.Shift(&nu, systWeight);
            // If there were only weighting systs applied then the cached
            // nominal values are still valid.
            shifted = caf::SRProxySystController::AnyShifted();
          }

          //unsigned int idxCut = 0; // testing
          for ( auto& [truthcut, treemap] : truthcutmap ) {
            const bool pass = shifted ? truthcut(&nu) : nomTruthCutCache.Get(truthcut, &nu);
            // Cut failed, skip all the histograms that depended on it
            if(!pass) continue;

            //unsigned int idxTree = 0;
            for ( std::map<Tree*, std::map<TruthVarOrMultiVar, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
              //unsigned int idxVar = 0;
              std::map<std::string, std::vector<double>> recordVals;
              unsigned int numEntries=0;
              for ( auto& [truthvarormulti, truthvarname] : treemapIt->second ) {
                //std::cout << "SpillCut " << idxSpillCut << " Shift " << idxShift << " Cut " << idxCut << " Tree " << idxTree << " Var " << idxVar << std::endl;
                if(truthvarormulti.IsMulti()){
                  auto const& truthvals = truthvarormulti.GetMultiVar()(&nu);
                  for(double truthval: truthvals) recordVals[truthvarname].push_back(truthval);
                  if (numEntries==0)    numEntries = truthvals.size();
                  continue;
                }

                const TruthVar& truthvar = truthvarormulti.GetVar();
                const double truthval = shifted ? truthvar(&nu) : nomTruthVarCache.Get(truthvar, &nu);

                //std::cout << "    VAL = " << turthval << std::endl;

                if(std::isnan(truthval) || std::isinf(truthval)){
                  std::cerr << "Warning: Bad value: " << truthval
                            << " returned from a Var. The input variable(s) could "
                            << "be NaN in the CAF, or perhaps your "
                            << "Var code computed 0/0?";
                  std::cout << " Still filling into the ''branch'' for this slice." << std::endl;
                }

                recordVals[truthvarname].push_back(truthval);
                if( numEntries==0 ) numEntries = 1;
                //idxVar+=1;
              } // end for truthvar/truthvarname
              // If fSaveRunSubrunEvt then fill these entries...
              if ( treemapIt->first->SaveRunSubEvent() ) {
                for ( unsigned int idxRun=0; idxRun<numEntries; ++idxRun ) {
                  recordVals["Run/i"].push_back( sr->hdr.run );
                  recordVals["Subrun/i"].push_back( sr->hdr.subrun );
                  recordVals["Evt/i"].push_back( sr->hdr.evt );
                }
              }
              // Adding CutType

              if ( treemapIt->first->SaveTruthCutType() ){
                // Loop over reco slices, and check if the truth-matched slice pass the (Slice)Cut
                bool HasMatchedSlicePassCut = false;
                for ( auto const& slc : sr->slc ) {
                  if ( slc.truth.index < 0 ) continue;
                  else if ( slc.truth.index != nu.index ) continue;
                  if( treemapIt->first->GetSignalSelectionCut()(&slc) ){
                    HasMatchedSlicePassCut = true;
                    break;
                  }
                }
                int tmp_CutType = HasMatchedSlicePassCut ? 1 : 0;
                recordVals["CutType/i"].push_back( tmp_CutType );

                int tmp_SpillCutType = spillpass ? 1 : 0;
                recordVals["SpillCutType/i"].push_back( tmp_SpillCutType );

              }



              treemapIt->first->UpdateEntries(recordVals);
              //idxTree+=1;
            } // end for tree
            //idxCut+=1;
          } // end for cut

          // Return StandardRecord to its unshifted form ready for the next
          // histogram.
          caf::SRProxySystController::Rollback();

          //idxShift+=1;
        } // end for shift
      } // end for slice
      //idxSpillCut+=1;
    } // end for spillcut

    // Weights trees
    // Sigma knobs
    for ( auto& [spillcut, shiftmap] : fNSigmasTreeDefs ) {
      const bool spillpass = spillcut(sr);
      // Cut failed, skip all the histograms that depend on it
      if(!spillpass) continue;

      unsigned int idxSlice = 0; // in case we want to save the slice number to the tree
      // NB: We DON'T keep track of Nominal Cut/Var/etc. because we want to shift/reset shifts for the weight saving... We sacrifice potential speed here by choice.
      for( caf::SRSliceProxy& slc: sr->slc ) {
        for ( auto& [shift, cutmap] : shiftmap ) {
          for ( auto& [cut, treemap] : cutmap ) {
            const bool pass = cut(&slc);
            // Cut failed, skip all the histograms that depended on it
            if(!pass) continue;

            for ( std::map<NSigmasTree*, std::map<const ISyst*, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {

              std::map<std::string, std::vector<double>> headerVals;
              std::map<std::string, std::vector<double>> recordVals;
              for ( auto& [syst, systname] : treemapIt->second ) {
                for ( int sigma=treemapIt->first->NSigmaLo(systname); sigma<=treemapIt->first->NSigmaHi(systname); ++sigma ) {

                  // Need to provide a clean slate for each new set of systematic
                  // shifts to work from. Copying the whole StandardRecord is pretty
                  // expensive, so modify it in place and revert it afterwards.
                  caf::SRProxySystController::BeginTransaction();

                  double systWeight = 1;
                  // Can special-case nominal to not pay cost of Shift()
                  if(!shift.IsNominal()){
                    shift.Shift(&slc, systWeight);
                  }

                  // Now shift for the weight we want to save
                  const SystShifts& shiftSigma = SystShifts(syst,sigma);
                  double systWeightSigma = 1;
                  shiftSigma.Shift(&slc, systWeightSigma);

                  recordVals[systname].push_back(systWeightSigma);

                  // Reset shifts to get the next sigma
                  caf::SRProxySystController::Rollback();
                }
              }
              // If fSaveRunSubrunEvt then fill these entries...
              if ( treemapIt->first->SaveRunSubEvent() ) {
                headerVals["Run/i"].push_back( sr->hdr.run );
                headerVals["Subrun/i"].push_back( sr->hdr.subrun );
                headerVals["Evt/i"].push_back( sr->hdr.evt );
              }
              if ( treemapIt->first->SaveSliceNum() ) {
                headerVals["Slice/i"].push_back( idxSlice );
              }

              treemapIt->first->UpdateEntries(headerVals,recordVals);
            } // end for tree
          } // end for cut
        } // end for shift
        idxSlice+=1;
      } // end for slice
    } // end for spillcut

    // Truth Sigma knobs
    // NB: We DON'T keep track of Nominal Cut/Var/etc. because we want to shift/reset shifts for the weight saving... We sacrifice potential speed here by choice.
    for(caf::SRTrueInteractionProxy& nu: sr->mc.nu){
      for ( auto& [shift, truthcutmap] : fTruthNSigmasTreeDefs ) {
        for ( auto& [truthcut, treemap] : truthcutmap ) {
          const bool pass = truthcut(&nu);
          // Cut failed, skip all the histograms that depended on it
          if(!pass) continue;

          for ( std::map<NSigmasTree*, std::map<const ISyst*, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {

            std::map<std::string, std::vector<double>> headerVals;
            std::map<std::string, std::vector<double>> recordVals;
            for ( auto& [syst, systname] : treemapIt->second ) {
              for ( int sigma=treemapIt->first->NSigmaLo(systname); sigma<=treemapIt->first->NSigmaHi(systname); ++sigma ) {

                // Need to provide a clean slate for each new set of systematic
                // shifts to work from. Copying the whole StandardRecord is pretty
                // expensive, so modify it in place and revert it afterwards.
                caf::SRProxySystController::BeginTransaction();

                double systWeight = 1;
                // Can special-case nominal to not pay cost of Shift()
                if(!shift.IsNominal()){
                  shift.Shift(&nu, systWeight);
                }

                // Now shift for the weight we want to save
                const SystShifts& shiftSigma = SystShifts(syst,sigma);
                double systWeightSigma = 1;
                shiftSigma.Shift(&nu, systWeightSigma);

                recordVals[systname].push_back(systWeightSigma);

                // Reset shifts to get the next sigma
                caf::SRProxySystController::Rollback();
              }
            }
            // If fSaveRunSubrunEvt then fill these entries...
            if ( treemapIt->first->SaveRunSubEvent() ) {
              headerVals["Run/i"].push_back( sr->hdr.run );
              headerVals["Subrun/i"].push_back( sr->hdr.subrun );
              headerVals["Evt/i"].push_back( sr->hdr.evt );
            }

            treemapIt->first->UpdateEntries(headerVals,recordVals);
          } // end for tree
        } // end for truthcut
      } // end for shift
    } // end for nu

    // Universe knobs
    for ( auto& [spillcut, shiftmap] : fNUniversesTreeDefs ) {
      const bool spillpass = spillcut(sr);
      // Cut failed, skip all the histograms that depend on it
      if(!spillpass) continue;

      unsigned int idxSlice = 0; // in case we want to save the slice number to the tree
      for( caf::SRSliceProxy& slc: sr->slc ) {
        // Some shifts only adjust the weight, so they're effectively nominal,
        // but aren't grouped with the other nominal histograms. Keep track of
        // the results for nominals in these caches to speed those systs up.
        CutVarCache<bool, Cut, caf::SRSliceProxy> nomCutCache;

        for ( auto& [shift, cutmap] : shiftmap ) {
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

          for ( auto& [cut, treemap] : cutmap ) {
            const bool pass = shifted ? cut(&slc) : nomCutCache.Get(cut, &slc);
            // Cut failed, skip all the histograms that depended on it
            if(!pass) continue;

            for ( std::map<NUniversesTree*, std::map<std::vector<VarOrMultiVar>, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
              std::map<std::string, std::vector<double>> headerVals;
              std::map<std::string, std::vector<double>> recordVals;
              for ( auto& [universes, systname] : treemapIt->second ) {
                for ( auto const& var : universes ) {
                  double val = var.GetVar()(&slc);
                  recordVals[systname].push_back(val);
                }
              }
              // If fSaveRunSubrunEvt then fill these entries...
              if ( treemapIt->first->SaveRunSubEvent() ) {
                headerVals["Run/i"].push_back( sr->hdr.run );
                headerVals["Subrun/i"].push_back( sr->hdr.subrun );
                headerVals["Evt/i"].push_back( sr->hdr.evt );
              }
              if ( treemapIt->first->SaveSliceNum() ) {
                headerVals["Slice/i"].push_back( idxSlice );
              }

              treemapIt->first->UpdateEntries(headerVals,recordVals);
            } // end for tree
          } // end for cut

          // Return StandardRecord to its unshifted form ready for the next
          // histogram.
          caf::SRProxySystController::Rollback();
        } // end for shift
        idxSlice+=1;
      } // end for slice
    } // end for spillcut

    // Universe knobs, Truth
    for(caf::SRTrueInteractionProxy& nu: sr->mc.nu){
      // Some shifts only adjust the weight, so they're effectively nominal,
      // but aren't grouped with the other nominal histograms. Keep track of
      // the results for nominals in these caches to speed those systs up.
      CutVarCache<bool, TruthCut, caf::SRTrueInteractionProxy> nomTruthCutCache;

      for ( auto& [shift, truthcutmap] : fTruthNUniversesTreeDefs ) {
        // Need to provide a clean slate for each new set of systematic
        // shifts to work from. Copying the whole StandardRecord is pretty
        // expensive, so modify it in place and revert it afterwards.
        caf::SRProxySystController::BeginTransaction();

        bool shifted = false;

        double systWeight = 1;
        // Can special-case nominal to not pay cost of Shift()
        if(!shift.IsNominal()){
          shift.Shift(&nu, systWeight);
          // If there were only weighting systs applied then the cached
          // nominal values are still valid.
          shifted = caf::SRProxySystController::AnyShifted();
        }

        for ( auto& [truthcut, treemap] : truthcutmap ) {
          const bool pass = shifted ? truthcut(&nu) : nomTruthCutCache.Get(truthcut, &nu);
          // Cut failed, skip all the histograms that depended on it
          if(!pass) continue;

          for ( std::map<NUniversesTree*, std::map<std::vector<TruthVarOrMultiVar>, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
            std::map<std::string, std::vector<double>> headerVals;
            std::map<std::string, std::vector<double>> recordVals;
            for ( auto& [universes, systname] : treemapIt->second ) {
              for ( auto const& truthvar : universes ) {
                double truthval = truthvar.GetVar()(&nu);
                recordVals[systname].push_back(truthval);
              }
            }
            // If fSaveRunSubrunEvt then fill these entries...
            if ( treemapIt->first->SaveRunSubEvent() ) {
              headerVals["Run/i"].push_back( sr->hdr.run );
              headerVals["Subrun/i"].push_back( sr->hdr.subrun );
              headerVals["Evt/i"].push_back( sr->hdr.evt );
            }

            treemapIt->first->UpdateEntries(headerVals,recordVals);
          } // end for tree
        } // end for truthcut

        // Return StandardRecord to its unshifted form ready for the next
        // histogram.
        caf::SRProxySystController::Rollback();
      } // end for shift

    } // end for nu


  }

  //----------------------------------------------------------------------
  void SpectrumLoader::StoreExposures()
  {
    if(fabs(fPOT - fPOTFromHist)/std::min(fPOT, fPOTFromHist) > 0.001){
      std::cout << fPOT << " POT from hdr differs from " << fPOTFromHist << " POT from the TotalPOT histogram!" << std::endl;
      abort();
    }

    std::cout << fPOT << " POT over " << fNReadouts << " readouts" << std::endl;

    for(auto& shiftdef: fHistDefs){
      for(auto& spillcutdef: shiftdef.second){
        for(auto& cutdef: spillcutdef.second){
          for(auto& weidef: cutdef.second){
            for(auto& vardef: weidef.second){
              for(Spectrum* s: vardef.second.spects){
                s->fPOT += fPOT;
                s->fLivetime += fNReadouts;
              }
              for(ReweightableSpectrum* rw: vardef.second.rwSpects){
                rw->fPOT += fPOT;
                rw->fLivetime += fNReadouts;
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
            s->fLivetime += fNReadouts;
          }
        }
      }
    }

    for(auto& spillcutdef: fTruthHistDefs){
      for(auto& shiftdef: spillcutdef.second){
        for(auto& truthcutdef: shiftdef.second){
          for(auto& truthweidef: truthcutdef.second){
            for(auto& truthvardef: truthweidef.second){
              for(Spectrum* s: truthvardef.second.spects){
                s->fPOT += fPOT;
                s->fLivetime += fNReadouts;
              }
            }
          }
        }
      }
    }

    for(auto& spillcutdef: fTruthHistWithCutDefs){
      for(auto& cutdef: spillcutdef.second){
        for(auto& shiftdef: cutdef.second){
          for(auto& truthcutdef: shiftdef.second){
            for(auto& truthweidef: truthcutdef.second){
              for(auto& truthvardef: truthweidef.second){
                for(Spectrum* s: truthvardef.second.spects){
                  s->fPOT += fPOT;
                  s->fLivetime += fNReadouts;
                }
              }
            }
          }
        }
      }
    }

    // Trees
    for ( auto& [spillcut, shiftmap] : fTreeDefs ) {
      for ( auto& [shift, cutmap] : shiftmap ) {
        for ( auto& [cut, treemap] : cutmap ) {
          for ( std::map<Tree*, std::map<VarOrMultiVar, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
            treemapIt->first->UpdateExposure(fPOT,fNReadouts);
          }
        }
      }
    }

    // Spill Trees
    for ( auto& [spillcut, treemap] : fSpillTreeDefs ) {
      for ( std::map<Tree*, std::map<SpillVarOrMultiVar, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
        treemapIt->first->UpdateExposure(fPOT,fNReadouts);
      }
    }

    // Truth Trees
    for ( auto& [spillcut, shiftmap] : fTruthTreeDefs ) {
      for ( auto& [shift, truthcutmap] : shiftmap ) {
        for ( auto& [truthcut, treemap] : truthcutmap ) {
          for ( std::map<Tree*, std::map<TruthVarOrMultiVar, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
            treemapIt->first->UpdateExposure(fPOT,fNReadouts);
          }
        }
      }
    }

    // NSigmasTrees
    for ( auto& [spillcut, shiftmap] : fNSigmasTreeDefs ) {
      for ( auto& [shift, cutmap] : shiftmap ) {
        for ( auto& [cut, treemap] : cutmap ) {
          for ( std::map<NSigmasTree*, std::map<const ISyst*, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
            treemapIt->first->UpdateExposure(fPOT,fNReadouts);
          }
        }
      }
    }

    // NSigmasTrees, Truth
    for ( auto& [shift, truthcutmap] : fTruthNSigmasTreeDefs ) {
      for ( auto& [truthcut, treemap] : truthcutmap ) {
        for ( std::map<NSigmasTree*, std::map<const ISyst*, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
          treemapIt->first->UpdateExposure(fPOT,fNReadouts);
        }
      } 
    } 

    // NUniversesTrees
    for ( auto& [spillcut, shiftmap] : fNUniversesTreeDefs ) {
      for ( auto& [shift, cutmap] : shiftmap ) {
        for ( auto& [cut, treemap] : cutmap ) {
          for ( std::map<NUniversesTree*, std::map<std::vector<VarOrMultiVar>, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
            treemapIt->first->UpdateExposure(fPOT,fNReadouts);
          }
        }
      }
    }

    // NUniversesTrees, Truth
    for ( auto& [shift, truthcutmap] : fTruthNUniversesTreeDefs ) {
      for ( auto& [truthcut, treemap] : truthcutmap ) {
        for ( std::map<NUniversesTree*, std::map<std::vector<TruthVarOrMultiVar>, std::string>>::iterator treemapIt=treemap.begin(); treemapIt!=treemap.end(); ++treemapIt ) {
            treemapIt->first->UpdateExposure(fPOT,fNReadouts);
        }     
      }     
    }    


  }
} // namespace
