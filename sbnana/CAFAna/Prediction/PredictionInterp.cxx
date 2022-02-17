#include "sbnana/CAFAna/Prediction/PredictionInterp.h"
#include "sbnana/CAFAna/Core/LoadFromFile.h"
#include "CAFAna/Core/Ratio.h"
#include "CAFAna/Core/Registry.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "TDirectory.h"
#include "TH2.h"
#include "TMatrixD.h"
#include "TObjString.h"
#include "TVectorD.h"

// For debug plots
#include "TGraph.h"
#include "TCanvas.h"

#include "OscLib/IOscCalc.h"
#include "sbnana/CAFAna/Core/MathUtil.h"

#include "sbnana/CAFAna/Core/Loaders.h"

#include <algorithm>
#include <malloc.h>

namespace ana
{
  //----------------------------------------------------------------------
  PredictionInterp::PredictionInterp(std::vector<const ISyst*> systs,
                                     osc::IOscCalc* osc,
                                     const IPredictionGenerator& predGen,
                                     Loaders& loaders,
                                     const SystShifts& shiftMC,
                                     EMode_t mode)
    : fOscOrigin(osc ? osc->Copy() : 0),
      fBinning(Spectrum::Uninitialized()),
      fSplitBySign(mode == kSplitBySign)
  {
    for(const ISyst* syst: systs){
      ShiftedPreds sp;
      sp.systName = syst->ShortName();

      for(int x = -syst->PredInterpMaxNSigma(); x <= +syst->PredInterpMaxNSigma(); ++x){
        sp.shifts.push_back(x);
      }

      for(int sigma: sp.shifts){
        SystShifts shiftHere = shiftMC;
        shiftHere.SetShift(syst, sigma);
        sp.preds.push_back(predGen.Generate(loaders, shiftHere).release());
      }

      fPreds.emplace(syst, sp);
    } // end for syst

    fPredNom = predGen.Generate(loaders, shiftMC);
  }

  //----------------------------------------------------------------------
  PredictionInterp::~PredictionInterp()
  {
    //    for(auto it: fPreds) for(IPrediction* p: it.second.preds) delete p;
    //    delete fOscOrigin;

    // It isn't really a unique ptr when we use PredictionInterpTemplates
    fPredNom.release();
  }

  //----------------------------------------------------------------------
  std::vector<std::vector<PredictionInterp::Coeffs>> PredictionInterp::
  FitRatios(const std::vector<double>& shifts,
            const std::vector<Eigen::ArrayXd>& ratios) const
  {
    if(ratios.size() < 2){
      std::cout << "PredictionInterp::FitRatios(): ratios.size() = " << ratios.size() << " - how did that happen?" << std::endl;

      abort();
    }

    assert(shifts.size() == ratios.size());

    std::vector<std::vector<Coeffs>> ret;

    const int binMax = ratios[0].size();

    for(int binIdx = 0; binIdx < binMax; ++binIdx){
      ret.push_back({});

      // This is cubic interpolation. For each adjacent set of four points we
      // determine coefficients for a cubic which will be the curve between the
      // center two. We constrain the function to match the two center points
      // and to have the right mean gradient at them. This causes this patch to
      // match smoothly with the next one along. The resulting function is
      // continuous and first and second differentiable. At the ends of the
      // range we fit a quadratic instead with only one constraint on the
      // slope. The coordinate conventions are that point y1 sits at x=0 and y2
      // at x=1. The matrices are simply the inverses of writing out the
      // constraints expressed above.

      // Special-case for linear interpolation
      if(ratios.size() == 2){
        const double y0 = ratios[0][binIdx];
        const double y1 = ratios[1][binIdx];

        ret.back().emplace_back(0, 0, y1-y0, y0);
        continue;
      }

      {
        const double y1 = ratios[0][binIdx];
        const double y2 = ratios[1][binIdx];
        const double y3 = ratios[2][binIdx];
        const double v[3] = {y1, y2, (y3-y1)/2};
        const double m[9] = { 1, -1,  1,
                             -2,  2, -1,
                              1,  0,  0};
        const TVectorD res = TMatrixD(3, 3, m) * TVectorD(3, v);
        ret.back().emplace_back(0, res(0), res(1), res(2));
      }

      // We're assuming here that the shifts are separated by exactly 1 sigma.
      for(unsigned int shiftIdx = 1; shiftIdx < ratios.size()-2; ++shiftIdx){
        const double y0 = ratios[shiftIdx-1][binIdx];
        const double y1 = ratios[shiftIdx  ][binIdx];
        const double y2 = ratios[shiftIdx+1][binIdx];
        const double y3 = ratios[shiftIdx+2][binIdx];

        const double v[4] = {y1, y2, (y2-y0)/2, (y3-y1)/2};
        const double m[16] = { 2, -2,  1,  1,
                              -3,  3, -2, -1,
                               0,  0,  1,  0,
                               1,  0,  0,  0};
        const TVectorD res = TMatrixD(4, 4, m) * TVectorD(4, v);
        ret.back().emplace_back(res(0), res(1), res(2), res(3));
      } // end for shiftIdx

      {
        const int N = ratios.size()-3;
        const double y0 = ratios[N  ][binIdx];
        const double y1 = ratios[N+1][binIdx];
        const double y2 = ratios[N+2][binIdx];
        const double v[3] = {y1, y2, (y2-y0)/2};
        const double m[9] = {-1,  1, -1,
                              0,  0,  1,
                              1,  0,  0};
        const TVectorD res = TMatrixD(3, 3, m) * TVectorD(3, v);
        ret.back().emplace_back(0, res(0), res(1), res(2));
      }
    } // end for binIdx

    double stride = -1;
    for(unsigned int i = 0; i < shifts.size()-1; ++i){
      const double newStride = shifts[i+1]-shifts[i];
      assert((stride < 0 || fabs(stride-newStride) < 1e-3) &&
             "Variably-spaced syst templates are unsupported");
      stride = newStride;
    }

    // If the stride is actually not 1, need to rescale all the coefficients
    for(std::vector<Coeffs>& cs: ret)
      for(Coeffs& c: cs){
        c = Coeffs(c.a/util::cube(stride),
                   c.b/util::sqr(stride),
                   c.c/stride,
                   c.d);}
    return ret;
  }

  //----------------------------------------------------------------------
  std::vector<std::vector<PredictionInterp::Coeffs>> PredictionInterp::
  FitComponent(const std::vector<double>& shifts,
               const std::vector<IPrediction*>& preds,
               Flavors::Flavors_t flav,
               Current::Current_t curr,
               Sign::Sign_t sign,
               const std::string& shortName) const
  {
    IPrediction* pNom = 0;
    for(unsigned int i = 0; i < shifts.size(); ++i){
      if(shifts[i] == 0) pNom = preds[i];
    }
    assert(pNom);

    // Do it this way rather than via fPredNom so that systematics evaluated
    // relative to some alternate nominal (eg Birks C where the appropriate
    // nominal is no-rock) can work.
    const Spectrum nom = pNom->PredictComponent(fOscOrigin,
                                                flav, curr, sign);

    std::vector<Eigen::ArrayXd> ratios;
    ratios.reserve(preds.size());
    for(auto& p: preds){
      ratios.emplace_back(Ratio(p->PredictComponent(fOscOrigin,
                                                    flav, curr, sign),
                                nom).GetEigen());

      // Check none of the ratio values is crazy
      Eigen::ArrayXd& r = ratios.back();
      for(int i = 0; i < r.size(); ++i){
        if (r[i] > 2){
          // std::cout << "PredictionInterp: WARNING, ratio in bin "
          // 	    << i << " for " << shifts[&p-&preds.front()]
          //           << " sigma shift of " << systName << " is " << y
          //           << " which exceeds limit of 2. Capping." << std::endl;
          r[i] = 2;
        }
      }
    }

    return FitRatios(shifts, ratios);
  }

  //----------------------------------------------------------------------
  void PredictionInterp::InitFitsHelper(ShiftedPreds& sp,
                                        std::vector<std::vector<std::vector<Coeffs>>>& fits,
                                        Sign::Sign_t sign) const
  {
    fits.resize(kNCoeffTypes);

    fits[kNueApp]   = FitComponent(sp.shifts, sp.preds, Flavors::kNuMuToNuE,  Current::kCC, sign, sp.systName);
    fits[kNueSurv]  = FitComponent(sp.shifts, sp.preds, Flavors::kNuEToNuE,   Current::kCC, sign, sp.systName);
    fits[kNumuSurv] = FitComponent(sp.shifts, sp.preds, Flavors::kNuMuToNuMu, Current::kCC, sign, sp.systName);

    fits[kNC]       = FitComponent(sp.shifts, sp.preds, Flavors::kAll, Current::kNC, sign, sp.systName);

    fits[kOther] = FitComponent(sp.shifts, sp.preds, Flavors::kNuEToNuMu | Flavors::kAllNuTau, Current::kCC, sign, sp.systName);
  }

  //----------------------------------------------------------------------
  void PredictionInterp::InitFits() const
  {
    // No systs
    if(fPreds.empty()){
      if(fBinning.POT() > 0 || fBinning.Livetime() > 0) return;
    }
    // Already initialized
    else if(!fPreds.empty() && !fPreds.begin()->second.fits.empty()) return;

    for(auto& it: fPreds){
      ShiftedPreds& sp = it.second;

      if(fSplitBySign){
        InitFitsHelper(sp, sp.fits, Sign::kNu);
        InitFitsHelper(sp, sp.fitsNubar, Sign::kAntiNu);
      }
      else{
        InitFitsHelper(sp, sp.fits, Sign::kBoth);
      }
      sp.nCoeffs = sp.fits[0][0].size();
    }

    // Predict something, anything, so that we can know what binning to use
    fBinning = fPredNom->Predict(fOscOrigin);
    fBinning.Clear();
  }

  //----------------------------------------------------------------------
  void PredictionInterp::SetOscSeed(osc::IOscCalc* oscSeed){
    fOscOrigin = oscSeed->Copy();
    for(auto& it: fPreds) it.second.fits.clear();
  }

  //----------------------------------------------------------------------
  Spectrum PredictionInterp::Predict(osc::IOscCalc* calc) const
  {
    return fPredNom->Predict(calc);
  }

  //----------------------------------------------------------------------
  Spectrum PredictionInterp::PredictComponent(osc::IOscCalc* calc,
                                              Flavors::Flavors_t flav,
                                              Current::Current_t curr,
                                              Sign::Sign_t sign) const
  {
    return fPredNom->PredictComponent(calc, flav, curr, sign);
  }

  //----------------------------------------------------------------------
  Spectrum PredictionInterp::PredictSyst(osc::IOscCalc* calc,
                                         const SystShifts& shift) const
  {
    InitFits();

    return PredictComponentSyst(calc, shift,
                                Flavors::kAll,
                                Current::kBoth,
                                Sign::kBoth);
  }

  //----------------------------------------------------------------------
  Spectrum PredictionInterp::ShiftSpectrum(const Spectrum& s, CoeffsType type,
                                           bool nubar,
                                           const SystShifts &shift) const
  {
    // Save time by not shifting a spectrum that is all zeros anyway
    if(s.POT() == 0 && s.Livetime() == 0) return s;

    Eigen::ArrayXd vec = s.GetEigen(s.POT());
    ShiftBins(vec.size(), vec.data(), type, nubar, shift);
    return Spectrum(std::move(vec), HistAxis(s.GetLabels(), s.GetBinnings()), s.POT(), s.Livetime());
  }

  //----------------------------------------------------------------------
  void PredictionInterp::ShiftBins(unsigned int N,
                                   double* arr,
                                   CoeffsType type,
                                   bool nubar,
                                   const SystShifts& shift) const
  {
    if(nubar) assert(fSplitBySign);

    std::vector<double> corr(N, 1);

    for(auto& it: fPreds){
      const ISyst* syst = it.first;
      const ShiftedPreds& sp = it.second;

      double x = shift.GetShift(syst);

      if(x == 0) continue;

      int shiftBin = (x - sp.shifts[0])/sp.Stride();
      shiftBin = std::max(0, shiftBin);
      shiftBin = std::min(shiftBin, sp.nCoeffs - 1);

      auto& fits = nubar ? sp.fitsNubar : sp.fits;

      x -= sp.shifts[shiftBin];

      const double x_cube = util::cube(x);
      const double x_sqr = util::sqr(x);

      for(unsigned int n = 0; n < N; ++n){
        // Uncomment to debug crashes in this function
        // assert(type < fits.size());
        // assert(n < sp.fits[type].size());
        // assert(shiftBin < int(sp.fits[type][n].size()));
        const Coeffs& f = fits[type][n][shiftBin];

        corr[n] *= f.a*x_cube + f.b*x_sqr + f.c*x + f.d;
      } // end for n

    } // end for syst

    for(unsigned int n = 0; n < N; ++n) {
      arr[n] *= (corr[n] > 0.) ? corr[n] : 0.;
    }
  }

  //----------------------------------------------------------------------
  Spectrum PredictionInterp::
  ShiftedComponent(osc::IOscCalc* calc,
                   const TMD5* hash,
                   const SystShifts& shift,
                   Flavors::Flavors_t flav,
                   Current::Current_t curr,
                   Sign::Sign_t sign,
                   CoeffsType type) const
  {
    if(fSplitBySign && sign == Sign::kBoth){
      return (ShiftedComponent(calc, hash, shift, flav, curr, Sign::kAntiNu, type)+
              ShiftedComponent(calc, hash, shift, flav, curr, Sign::kNu,     type));
    }

    // Must be the base case of the recursion to use the cache. Otherwise we
    // can cache systematically shifted versions of our children, which is
    // wrong. Also, some calcs won't hash themselves.
    const bool canCache = (hash != 0);

    const Key_t key = {flav, curr, sign};
    auto it = fNomCache.find(key);

    // Should the interpolation use the nubar fits?
    const bool nubar = (fSplitBySign && sign == Sign::kAntiNu);

    // We have the nominal for this exact combination of flav, curr, sign, calc
    // stored.  Shift it and return.
    if(canCache && it != fNomCache.end() && it->second.hash == *hash){
      return ShiftSpectrum(it->second.nom, type, nubar, shift);
    }

    // We need to compute the nominal again for whatever reason
    const Spectrum nom = fPredNom->PredictComponent(calc, flav, curr, sign);

    if(canCache){
      // Insert into the cache if not already there, or update if there but
      // with old oscillation parameters.
      if(it == fNomCache.end())
        fNomCache.emplace(key, Val_t({*hash, nom}));
      else
        it->second = {*hash, nom};
    }

    return ShiftSpectrum(nom, type, nubar, shift);
  }

  //----------------------------------------------------------------------
  Spectrum PredictionInterp::PredictComponentSyst(osc::IOscCalc* calc,
                                                  const SystShifts& shift,
                                                  Flavors::Flavors_t flav,
                                                  Current::Current_t curr,
                                                  Sign::Sign_t sign) const
  {
    InitFits();

    Spectrum& ret = fBinning;
    ret.Clear();

    if(ret.POT()==0) ret.OverridePOT(1e24);

    // Check that we're able to handle all the systs we were passed
    for(const ISyst* syst: shift.ActiveSysts()){
      if(fPreds.find(syst) == fPreds.end()){
        std::cerr << "This PredictionInterp is not set up to handle the requested systematic: " << syst->ShortName() << std::endl;
        abort();
      }
    } // end for syst


    const TMD5* hash = calc ? calc->GetParamsHash() : 0;

    if(curr & Current::kCC){
      if(flav & Flavors::kNuEToNuE)    ret += ShiftedComponent(calc, hash, shift, Flavors::kNuEToNuE,    Current::kCC, sign, kNueSurv);
      if(flav & Flavors::kNuEToNuMu)   ret += ShiftedComponent(calc, hash, shift, Flavors::kNuEToNuMu,   Current::kCC, sign, kOther  );
      if(flav & Flavors::kNuEToNuTau)  ret += ShiftedComponent(calc, hash, shift, Flavors::kNuEToNuTau,  Current::kCC, sign, kOther  );

      if(flav & Flavors::kNuMuToNuE)   ret += ShiftedComponent(calc, hash, shift, Flavors::kNuMuToNuE,   Current::kCC, sign, kNueApp  );
      if(flav & Flavors::kNuMuToNuMu)  ret += ShiftedComponent(calc, hash, shift, Flavors::kNuMuToNuMu,  Current::kCC, sign, kNumuSurv);
      if(flav & Flavors::kNuMuToNuTau) ret += ShiftedComponent(calc, hash, shift, Flavors::kNuMuToNuTau, Current::kCC, sign, kOther   );
    }
    if(curr & Current::kNC){
      assert(flav == Flavors::kAll); // Don't know how to calculate anything else

      ret += ShiftedComponent(calc, hash, shift, Flavors::kAll, Current::kNC, sign, kNC);
    }

    delete hash;

    return ret;
  }

  //----------------------------------------------------------------------
  void PredictionInterp::SaveTo(TDirectory* dir) const
  {
    TDirectory* tmp = gDirectory;

    dir->cd();
    TObjString("PredictionInterp").Write("type");


    fPredNom->SaveTo(dir->mkdir("pred_nom"));


    for(auto it: fPreds){
      const ShiftedPreds& sp = it.second;

      for(unsigned int i = 0; i < sp.shifts.size(); ++i){
        if(!sp.preds[i]){
          std::cout << "Can't save a PredictionInterp after MinimizeMemory()" << std::endl;
          abort();
        }
        sp.preds[i]->SaveTo(dir->mkdir(TString::Format("pred_%s_%+d",
                                                       sp.systName.c_str(),
                                                       int(sp.shifts[i])).Data()));
      } // end for i
    } // end for it

    ana::SaveTo(*fOscOrigin, dir->mkdir("osc_origin"));

    if(!fPreds.empty()){
      TH1F hSystNames("syst_names", ";Syst names", fPreds.size(), 0, fPreds.size());
      int binIdx = 1;
      for(auto it: fPreds){
        hSystNames.GetXaxis()->SetBinLabel(binIdx++, it.second.systName.c_str());
      }
      hSystNames.Write("syst_names");
    }

    TObjString split_sign = fSplitBySign ? "yes" : "no";
    split_sign.Write("split_sign");

    tmp->cd();
  }

  //----------------------------------------------------------------------
  std::unique_ptr<PredictionInterp> PredictionInterp::LoadFrom(TDirectory* dir)
  {
    TObjString* tag = (TObjString*)dir->Get("type");
    assert(tag);
    assert(tag->GetString() == "PredictionInterp");

    std::unique_ptr<PredictionInterp> ret(new PredictionInterp);

    LoadFromBody(dir, ret.get());

    TObjString* split_sign = (TObjString*)dir->Get("split_sign");
    // Can be missing from old files
    ret->fSplitBySign = (split_sign && split_sign->String() == "yes");

    return ret;
  }

  //----------------------------------------------------------------------
  void PredictionInterp::LoadFromBody(TDirectory* dir, PredictionInterp* ret,
                                      std::vector<const ISyst*> veto)
  {
    ret->fPredNom = ana::LoadFrom<IPrediction>(dir->GetDirectory("pred_nom"));

    TH1* hSystNames = (TH1*)dir->Get("syst_names");
    if(hSystNames){
      for(int systIdx = 0; systIdx < hSystNames->GetNbinsX(); ++systIdx){
        ShiftedPreds sp;
        sp.systName = hSystNames->GetXaxis()->GetBinLabel(systIdx+1);

        const ISyst* syst = Registry<ISyst>::ShortNameToPtr(sp.systName);

        if(std::find(veto.begin(), veto.end(), syst) != veto.end()) continue;

        for(int shift = -3; shift <= +3; ++shift){
          TDirectory* preddir = dir->GetDirectory(TString::Format("pred_%s_%+d", sp.systName.c_str(), shift).Data());
          if(!preddir) continue; // Can happen for genie systs

          IPrediction* pred = ana::LoadFrom<IPrediction>(preddir).release();

          sp.shifts.push_back(shift);
          sp.preds.push_back(pred);
        } // end for shift

        ret->fPreds.emplace(syst, sp);
      } // end for systIdx
    } // end if hSystNames

    ret->fOscOrigin = ana::LoadFrom<osc::IOscCalc>(dir->GetDirectory("osc_origin")).release();
  }

  //----------------------------------------------------------------------
  void PredictionInterp::MinimizeMemory()
  {
    InitFits();

    std::set<IPrediction*> todel;
    for(auto& it: fPreds){
      std::vector<IPrediction*>& preds = it.second.preds;
      for(unsigned int i = 0; i < preds.size(); ++i){
        if(preds[i] != fPredNom.get()){
          todel.insert(preds[i]);
          preds[i] = 0;
        }
      }
    }

    for(IPrediction* p: todel) delete p;

    // We probably just freed up a lot of memory, but malloc by default hangs
    // on to all of it as cache.
    #ifndef DARWINBUILD
    malloc_trim(0);
    #endif
  }

  //----------------------------------------------------------------------
  void PredictionInterp::DebugPlot(const ISyst* syst,
                                   osc::IOscCalc* calc,
                                   Flavors::Flavors_t flav,
                                   Current::Current_t curr,
                                   Sign::Sign_t sign) const
  {
    DontAddDirectory guard;

    InitFits();

    auto it = fPreds.find(syst);
    if(it == fPreds.end()){
      std::cout << "PredictionInterp::DebugPlots(): "
                << syst->ShortName() << " not found" << std::endl;
      return;
    }

    std::unique_ptr<TH1> nom(fPredNom->PredictComponent(calc, flav, curr, sign).ToTH1(18e20));
    const int nbins = nom->GetNbinsX();

    std::vector<TGraph*> curves(nbins);
    std::vector<TGraph*> points(nbins);

    for(int i = 0; i <= 80; ++i){
      const double x = .1*i-4;
      const SystShifts ss(it->first, x);
      std::unique_ptr<TH1> h(PredictComponentSyst(calc, ss, flav, curr, sign).ToTH1(18e20));

      for(int bin = 0; bin < nbins; ++bin){
        if(i == 0){
          curves[bin] = new TGraph;
          points[bin] = new TGraph;
        }

        const double ratio = h->GetBinContent(bin+1)/nom->GetBinContent(bin+1);

        if(!std::isnan(ratio)) curves[bin]->SetPoint(curves[bin]->GetN(), x, ratio);
        else curves[bin]->SetPoint(curves[bin]->GetN(), x, 1);
      } // end for bin
    } // end for i (x)

    // As elswhere, to allow BirksC etc that need a different nominal to plot
    // right.
    IPrediction* pNom = 0;
    for(unsigned int shiftIdx = 0; shiftIdx < it->second.shifts.size(); ++shiftIdx){
      if(it->second.shifts[shiftIdx] == 0) pNom = it->second.preds[shiftIdx];
    }
    assert(pNom);
    std::unique_ptr<TH1> hnom(pNom->PredictComponent(calc, flav, curr, sign).ToTH1(18e20));

    for(unsigned int shiftIdx = 0; shiftIdx < it->second.shifts.size(); ++shiftIdx){
      if(!it->second.preds[shiftIdx]) continue; // Probably MinimizeMemory()
      std::unique_ptr<TH1> h(it->second.preds[shiftIdx]->PredictComponent(calc, flav, curr, sign).ToTH1(18e20));

      for(int bin = 0; bin < nbins; ++bin){
        const double ratio = h->GetBinContent(bin+1)/hnom->GetBinContent(bin+1);
        if(!std::isnan(ratio)) points[bin]->SetPoint(points[bin]->GetN(), it->second.shifts[shiftIdx], ratio);
        else points[bin]->SetPoint(points[bin]->GetN(), it->second.shifts[shiftIdx], 1);
      }
    } // end for shiftIdx


    int nx = int(sqrt(nbins));
    int ny = int(sqrt(nbins));
    if(nx*ny < nbins) ++nx;
    if(nx*ny < nbins) ++ny;

    TCanvas* c = new TCanvas;
    c->Divide(nx, ny);

    for(int bin = 0; bin < nbins; ++bin){
      c->cd(bin+1);
      (new TH2F(UniqueName().c_str(),
                TString::Format("%s %g < %s < %g;Shift;Ratio",
                                it->second.systName.c_str(),
                                nom->GetXaxis()->GetBinLowEdge(bin+1),
                                nom->GetXaxis()->GetTitle(),
                                nom->GetXaxis()->GetBinUpEdge(bin+1)),
                100, -4, +4, 100, .5, 1.5))->Draw();
      curves[bin]->Draw("l same");
      points[bin]->SetMarkerStyle(kFullDotMedium);
      points[bin]->Draw("p same");
    } // end for bin

    c->cd(0);
  }

  //----------------------------------------------------------------------
  void PredictionInterp::DebugPlots(osc::IOscCalc* calc,
				    const std::string& savePattern,
				    Flavors::Flavors_t flav,
				    Current::Current_t curr,
				    Sign::Sign_t sign) const
  {
    const bool save = !savePattern.empty();
    const bool multiFile = savePattern.find("%s") != std::string::npos;

    if(save && !multiFile){
      new TCanvas;
      gPad->Print((savePattern+"[").c_str());
    }

    for(auto& it: fPreds){
      DebugPlot(it.first, calc, flav, curr, sign);

      if(save){
        if(multiFile){
          gPad->Print(TString::Format(savePattern.c_str(), it.second.systName.c_str()).Data());
        }
        else{
          gPad->Print(savePattern.c_str());
        }
      }
    } // end for it

    if(save && !multiFile){
      gPad->Print((savePattern+"]").c_str());
    }
  }

  //----------------------------------------------------------------------
  void PredictionInterp::DebugPlotColz(const ISyst* syst,
                                       osc::IOscCalc* calc,
                                       Flavors::Flavors_t flav,
                                       Current::Current_t curr,
                                       Sign::Sign_t sign) const
  {
    InitFits();

    std::unique_ptr<TH1> nom(fPredNom->PredictComponent(calc, flav, curr, sign).ToTH1(18e20));
    const int nbins = nom->GetNbinsX();

    TH2* h2 = new TH2F("", (syst->LatexName()+";;#sigma").c_str(),
                       nbins, nom->GetXaxis()->GetXmin(), nom->GetXaxis()->GetXmax(),
                       80, -4, +4);
    h2->GetXaxis()->SetTitle(nom->GetXaxis()->GetTitle());

    for(int i = 1; i <= 80; ++i){
      const double y = h2->GetYaxis()->GetBinCenter(i);
      const SystShifts ss(syst, y);
      std::unique_ptr<TH1> h(PredictComponentSyst(calc, ss, flav, curr, sign).ToTH1(18e20));

      for(int bin = 0; bin < nbins; ++bin){
        const double ratio = h->GetBinContent(bin+1)/nom->GetBinContent(bin+1);

        if(!isnan(ratio) && !isinf(ratio))
          h2->Fill(h2->GetXaxis()->GetBinCenter(bin), y, ratio);
      } // end for bin
    } // end for i (x)

    h2->Draw("colz");
    h2->SetMinimum(0.5);
    h2->SetMaximum(1.5);
  }

  //----------------------------------------------------------------------
  void PredictionInterp::DebugPlotsColz(osc::IOscCalc* calc,
                                        const std::string& savePattern,
                                        Flavors::Flavors_t flav,
                                        Current::Current_t curr,
                                        Sign::Sign_t sign) const
  {
    const bool save = !savePattern.empty();
    const bool multiFile = savePattern.find("%s") != std::string::npos;

    if(save && !multiFile){
      new TCanvas;
      gPad->Print((savePattern+"[").c_str());
    }

    for(auto it: fPreds){
      new TCanvas;
      DebugPlotColz(it.first, calc, flav, curr, sign);

      if(save){
        if(multiFile){
          gPad->Print(TString::Format(savePattern.c_str(), it.second.systName.c_str()).Data());
        }
        else{
          gPad->Print(savePattern.c_str());
        }
      }
    } // end for it

    if(save && !multiFile){
      gPad->Print((savePattern+"]").c_str());
    }
  }

} // namespace
