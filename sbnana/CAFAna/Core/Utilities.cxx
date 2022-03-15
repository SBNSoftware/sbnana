#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/Ratio.h"
#include "sbnana/CAFAna/Core/MathUtil.h"

#include "TArrayD.h"
#include "TClass.h"
#include "TDirectory.h"
#include "TH2.h"
#include "TH3.h"
#include "TF1.h"
#include "TMatrixD.h"
#include "TObjString.h"
#include "TString.h"
#include "TTree.h"
#include "TVector3.h"
#include "TVectorD.h"

#include <vector>
#include <cassert>
#include <cmath>
#include <fstream>
#include "sys/stat.h"
#include "wordexp.h"

namespace ana
{
  double LLPerBinFracSystErr::fgErr = -1;

  //----------------------------------------------------------------------
  std::string UniqueName()
  {
    static int N = 0;
    return TString::Format("cafanauniq%d", N++).Data();
  }

  //----------------------------------------------------------------------
  DontAddDirectory::DontAddDirectory()
  {
    fBackup = TH1::AddDirectoryStatus();
    TH1::AddDirectory(false);
  }

  //----------------------------------------------------------------------
  DontAddDirectory::~DontAddDirectory()
  {
    TH1::AddDirectory(fBackup);
  }

  //----------------------------------------------------------------------
  IFDHSilent::IFDHSilent()
  {
    const char* s = getenv("IFDH_SILENT");
    fSet = (s && s == std::string("0"));
    if(!fSet) setenv("IFDH_SILENT", "1", 1);
    if(s) std::cout << "IFDH_SILENT=" << s << std::endl;
  }

  //----------------------------------------------------------------------
  IFDHSilent::~IFDHSilent()
  {
    if(!fSet) unsetenv("IFDH_SILENT");
  }

  //----------------------------------------------------------------------
  FloatingExceptionOnNaN::FloatingExceptionOnNaN(bool enable)
  {
    // Don't want any pending FPEs to trigger when we flip exceptions
    // on. Whoever had them off previously had a reason.
    feclearexcept(FE_INVALID);

    fegetexceptflag(&fBackup, FE_INVALID);

#ifndef DARWINBUILD
    if(enable)
      feenableexcept(FE_INVALID);
    else
      fedisableexcept(FE_INVALID);
#else
    std::cerr << "WARNING: CAFAna/Core/Utilities.cxx built on OS X, no feenableexcept available" << std::endl;
#endif
  }

  //----------------------------------------------------------------------
  FloatingExceptionOnNaN::~FloatingExceptionOnNaN()
  {
    fesetexceptflag(&fBackup, FE_INVALID);
  }

  //----------------------------------------------------------------------
  std::string Experiment()
  {
    const char* ret = getenv("EXPERIMENT");

    if(!ret){
      std::cout << "\nERROR: Environment variable $EXPERIMENT not set." << std::endl;
      exit(1);
    }

    return ret;
  }

  //----------------------------------------------------------------------
  std::string SAMExperiment()
  {
    const char* ret = getenv("SAM_EXPERIMENT");

    if(!ret){
      std::cout << "\nERROR: Environment variable $SAM_EXPERIMENT not set.\nThis is required for various ifdh/sam functionality.\nYou likely want it to be set to 'sbn', though it should have been set automatically by setup scripts." << std::endl;
      exit(1);
    }

    return ret;
  }

  //----------------------------------------------------------------------
  std::unique_ptr<TMatrixD> CalcCovMx(const std::vector<TArrayD*> & binSets, int firstBin, int lastBin)
  {
    if (binSets.size() < 2)
      return std::unique_ptr<TMatrixD>(nullptr);

    if (lastBin < 0)
      lastBin = binSets[0]->GetSize() - 1;  // indexed from 0

    int nBins = lastBin - firstBin + 1;  // firstBin and lastBin are inclusive

    std::vector<double> binMeans(nBins);
    for( const auto & binSet : binSets )
    {
      for ( decltype(lastBin) binIdx = firstBin; binIdx <= lastBin; binIdx++ )
        binMeans[binIdx] += (*binSet)[binIdx];
    }
    for (decltype(lastBin) binIdx = firstBin; binIdx <= lastBin; binIdx++)
      binMeans[binIdx] /= binSets.size();


    auto covmx = std::make_unique<TMatrixD>(nBins, nBins);

    for( unsigned int hist_idx = 0; hist_idx < binSets.size(); ++hist_idx )
    {
      // first calculate the weighted sum of squares of the deviations
      for( decltype(nBins) i = 0; i < nBins; i++ )
      {
        double xi = (*(binSets[hist_idx]))[i];
        for( decltype(nBins) k = i; k < nBins; k++ )
        {
          double xk = (*(binSets[hist_idx]))[k];
          (*covmx)[i][k] += (xi - binMeans[i]) * (xk - binMeans[k]);
          if (i != k)
            (*covmx)[k][i] = (*covmx)[i][k];  // covariance matrices are always symmetric
        }
      }
    } // for (hist_idx)

    // now divide by N-1 to get sample covariance
    (*covmx) *= 1./(binSets.size()-1);

    return covmx;
  }

  //----------------------------------------------------------------------
  double LogLikelihood(double e, double o)
  {
    // http://www.wolframalpha.com/input/?i=d%2Fds+m*(1%2Bs)+-d+%2B+d*ln(d%2F(m*(1%2Bs)))%2Bs%5E2%2FS%5E2%3D0
    // http://www.wolframalpha.com/input/?i=solve+-d%2F(s%2B1)%2Bm%2B2*s%2FS%5E2%3D0+for+s
    const double S = LLPerBinFracSystErr::GetError();
    if(S > 0){
      const double S2 = util::sqr(S);
      const double s = .25*(sqrt(8*o*S2+util::sqr(e*S2-2))-e*S2-2);
      e *= 1+s;
    }

    // With this value, negative expected events and one observed
    // event gives a chisq from this one bin of 182.
    const double minexp = 1e-40; // Don't let expectation go lower than this

    assert(o >= 0);
    if(e < minexp){
      if(!o) return 0;
      e = minexp;
    }

    if(o*1000 > e){
      // This strange form is for numerical stability when e~o
      return 2*o*((e-o)/o + log1p((o-e)/e));
    }
    else{
      // But log1p doesn't like arguments near -1 (observation much smaller
      // than expectation), and it's better to use the usual formula in that
      // case.
      if(o){
        return 2*(e-o + o*log(o/e));
      }
      else{
        return 2*e;
      }
    }
  }

  //----------------------------------------------------------------------
  double LogLikelihood(const TH1* eh, const TH1* oh, bool useOverflow)
  {
    assert(eh->GetNbinsX() == oh->GetNbinsX());

    double chi = 0;

    int bufferBins = useOverflow? 2 : 1;

    for(int i = 0; i < eh->GetNbinsX()+bufferBins; ++i){
      const double e = eh->GetBinContent(i);
      const double o = oh->GetBinContent(i);

      chi += LogLikelihood(e, o);
    }

    return chi;
  }

  //----------------------------------------------------------------------
  double Chi2CovMx(const TVectorD* e, const TVectorD* o, const TMatrixD* covmxinv)
  {
    assert (e->GetNrows() == o->GetNrows());

    TVectorD diff = *o - *e;
    return diff * ((*covmxinv) * diff);  // operator* for two TVectorDs is the "dot product" (i.e., v1 * v2 = v1^{trans}v1)
  }

  //----------------------------------------------------------------------
  double Chi2CovMx(const TH1* e, const TH1* o, const TMatrixD* covmxinv)
  {
    TVectorD eVec(e->GetNbinsX());
    TVectorD oVec(o->GetNbinsX());
    for (int bin = 1; bin <= e->GetNbinsX(); bin++)
      eVec[bin-1] = e->GetBinContent(bin);
    for (int bin = 1; bin <= o->GetNbinsX(); bin++)
      oVec[bin-1] = o->GetBinContent(bin);

    return Chi2CovMx(&eVec, &oVec, covmxinv);
  }

  //----------------------------------------------------------------------
  TH2F* ExpandedHistogram(const std::string& title,
                          int nbinsx, double xmin, double xmax, bool xlog,
                          int nbinsy, double ymin, double ymax, bool ylog)
  {
    DontAddDirectory guard;

    if(xlog){xmin = log(xmin); xmax = log(xmax);}
    if(ylog){ymin = log(ymin); ymax = log(ymax);}

    // How wide the bins will be once we're done
    const double xwidth = (xmax-xmin)/(nbinsx-1);
    const double ywidth = (ymax-ymin)/(nbinsy-1);

    // Move the bin edges so that the limits occur at the centres
    xmin -= xwidth/2; ymin -= ywidth/2;
    xmax += xwidth/2; ymax += ywidth/2;

    std::vector<double> xedges(nbinsx+1);
    std::vector<double> yedges(nbinsy+1);

    for(int i = 0; i <= nbinsx; ++i){
      xedges[i] = xmin + (xmax-xmin)*i/double(nbinsx);
      if(xlog) xedges[i] = exp(xedges[i]);
    }
    for(int i = 0; i <= nbinsy; ++i){
      yedges[i] = ymin + (ymax-ymin)*i/double(nbinsy);
      if(ylog) yedges[i] = exp(yedges[i]);
    }

    return new TH2F(UniqueName().c_str(), title.c_str(),
                    nbinsx, &xedges.front(),
                    nbinsy, &yedges.front());
  }


  //----------------------------------------------------------------------
  std::unique_ptr<TMatrixD> SymmMxInverse(const TMatrixD& mx)
  {
    // check if there are any null rows/columns.
    // if there are, they make the matrix singular.
    // we will remove them temporarily,
    // invert the matrix, then put them back afterwards.
    std::set<int> nullRows;
    for (auto row = mx.GetRowLwb(); row <= mx.GetRowUpb(); row++)
    {
      bool rowIsNull = true;
      for (auto col = mx.GetColLwb(); col <= mx.GetColUpb(); col++)
      {
        if (mx[row][col] != 0.)
        {
          rowIsNull = false;
          break;
        }
      }

      if (rowIsNull)
        nullRows.insert(row);
    }

    std::cerr << " Notice: covariance matrix '" << mx.GetName() << "' has " << nullRows.size() << " null rows.\n"
        << "They will be removed before inverting and added back afterwards." << std::endl;

    // create a new matrix for inverting, skipping the null rows
    auto invMx = std::make_unique<TMatrixD>(mx.GetRowLwb(), mx.GetRowUpb() - nullRows.size(),
                                           mx.GetColLwb(), mx.GetColUpb() - nullRows.size());
    unsigned int skippedRows = 0;
    for (auto row = mx.GetRowLwb(); row <= mx.GetRowUpb(); row++)
    {
      if (nullRows.find(row) != nullRows.end())
      {
        skippedRows++;
        continue;
      }
      unsigned int skippedCols = 0;
      for (auto col = mx.GetColLwb(); col <= mx.GetColUpb(); col++)
      {
        // since we assumed the matrix is symmetric,
        // we can just use the null rows list here
        if (nullRows.find(col) != nullRows.end())
        {
          skippedCols++;
          continue;
        }

        (*invMx)[col-skippedCols][row-skippedRows] = (*invMx)[row-skippedRows][col-skippedCols] = mx[row][col];
      }
    }

    invMx->Invert();

    // put back the empty rows if there were any
    if (nullRows.size())
    {
      skippedRows = 0;
      auto retMx = std::make_unique<TMatrixD>(mx.GetRowLwb(), mx.GetRowUpb(),
                                              mx.GetColLwb(), mx.GetColUpb());
      for (auto row = mx.GetRowLwb(); row <= mx.GetRowUpb(); row++)
      {
        if (nullRows.find(row) != nullRows.end())
        {
          skippedRows++;
          continue;
        }

        unsigned int skippedCols = skippedRows;
        for (auto col = row; col <= mx.GetColUpb(); col++)
        {
          if (nullRows.find(col) != nullRows.end())
          {
            skippedCols++;
            continue;
          }

          (*retMx)[col][row] = (*retMx)[row][col] = (*invMx)[row-skippedRows][col-skippedCols];
        }
      }

      return retMx;
    }

    return invMx;
  }

  // Helper functions for MakeTHND().
  namespace{
    // Eventually the bin parameters will all be unpacked and we just pass them
    // on to the regular constructor.
    template<class T, class... A> T* MakeHist(A... a)
    {
      DontAddDirectory guard;
      return new T(a...);
    }

    // This function consumes bins from the start of the argument list and
    // pushes their translations onto the list of arguments at the end.
    template<class T, class... A> T* MakeHist(const Binning& firstBin,
                                              A... args)
    {
      if(firstBin.IsSimple())
        return MakeHist<T>(args...,
                           firstBin.NBins(), firstBin.Min(), firstBin.Max());
      else
        return MakeHist<T>(args...,
                           firstBin.NBins(), &firstBin.Edges().front());
    }
  }

  // Concrete instantiations. MakeHist() requires us to put the bin arguments
  // first...
  //----------------------------------------------------------------------
  TH1D* MakeTH1D(const char* name, const char* title, const Binning& bins)
  {
    return MakeHist<TH1D>(bins, name, title);
  }

  //----------------------------------------------------------------------
  TH2D* MakeTH2D(const char* name, const char* title,
                 const Binning& binsx,
                 const Binning& binsy)
  {
    return MakeHist<TH2D>(binsx, binsy, name, title);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2(const Spectrum& s, double exposure, ana::EExposureType expotype,
             const Binning& binsx, const Binning& binsy, ana::EBinType bintype)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(s.ToTH1(exposure, expotype));
    return ToTH2Helper(std::move(h1), binsx, binsy, bintype);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2(const Ratio& r,
             const Binning& binsx, const Binning& binsy)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(r.ToTH1());
    return ToTH2Helper(std::move(h1), binsx, binsy);
  }

  //----------------------------------------------------------------------
  TH2* ToTH2Helper(std::unique_ptr<TH1> h1,
		   const Binning& binsx, const Binning& binsy,
		   ana::EBinType bintype)
  {
    // Make sure it's compatible with having been made with this binning
    assert(h1->GetNbinsX() == binsx.NBins()*binsy.NBins());

    TH2* ret = MakeTH2D("", UniqueName().c_str(), binsx, binsy);

    for(int i = 0; i < h1->GetNbinsX(); ++i){
      const double val = h1->GetBinContent(i+1);
      const double err = h1->GetBinError(i+1);

      const int ix = i / binsy.NBins();
      const int iy = i % binsy.NBins();

      ret->SetBinContent(ix+1, iy+1, val);
      ret->SetBinError  (ix+1, iy+1, err);
    }

    if(bintype == ana::EBinType::kBinDensity) ret->Scale(1, "width");

    return ret;
  }

  //----------------------------------------------------------------------

  TH3* ToTH3(const Spectrum& s, double exposure, ana::EExposureType expotype,
             const Binning& binsx, const Binning& binsy, const Binning& binsz,
	     ana::EBinType bintype)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(s.ToTH1(exposure, expotype));

    return ToTH3Helper(std::move(h1), binsx, binsy, binsz, bintype);
  }

  //----------------------------------------------------------------------

  TH3* ToTH3(const Ratio& r,
             const Binning& binsx, const Binning& binsy, const Binning& binsz)
  {
    DontAddDirectory guard;

    std::unique_ptr<TH1> h1(r.ToTH1());

    return ToTH3Helper(std::move(h1), binsx, binsy, binsz);
  }

  //----------------------------------------------------------------------
  TH3* ToTH3Helper(std::unique_ptr<TH1> h1,
		   const Binning& binsx,
		   const Binning& binsy,
		   const Binning& binsz,
		   ana::EBinType bintype)
  {

    const int nx = binsx.NBins();
    const int ny = binsy.NBins();
    const int nz = binsz.NBins();

    // Make sure it's compatible with having been made with this binning
    assert(h1->GetNbinsX() == nx*ny*nz);

    TH3* ret;

    // If all three axes are simple, we can call a simpler constructor
    if(binsx.IsSimple() && binsy.IsSimple() && binsz.IsSimple()){
      ret = new TH3F(UniqueName().c_str(), "",
                     nx, binsx.Min(), binsx.Max(),
                     ny, binsy.Min(), binsy.Max(),
                     nz, binsz.Min(), binsz.Max());

      if(!binsx.IsSimple() || !binsy.IsSimple() || !binsz.IsSimple()){
        // TH3 doesn't have the constructors for mixed simple and custom
        std::cerr << "ToTH3: one or more axes is custom, but not all three. Applying Simple binning to all three axes" << std::endl;
      }
    }
    else{
      ret = new TH3F(UniqueName().c_str(), "",
                     nx, &binsx.Edges().front(),
                     ny, &binsy.Edges().front(),
                     nz, &binsz.Edges().front());
    }

    for(int i = 0; i < h1->GetNbinsX(); ++i){
      const double val = h1->GetBinContent(i+1);
      const double err = h1->GetBinError(i+1);

      const int nynz = ny*nz;
      const int nmodnynz = i%nynz;
      const int ix = i/nynz;
      const int iy = nmodnynz/nz;
      const int iz = i%nz;

      ret->SetBinContent(ix+1, iy+1, iz+1, val);
      ret->SetBinError  (ix+1, iy+1, iz+1, err);
    }

    if(bintype == ana::EBinType::kBinDensity) ret->Scale(1, "width");

    return ret;

  }

  //----------------------------------------------------------------------
  std::vector<std::string> Wildcard(const std::string& wildcardString)
  {
    // Expand environment variables and wildcards like the shell would
    wordexp_t p;
    const int status = wordexp(wildcardString.c_str(), &p, WRDE_SHOWERR);

    if(status != 0){
      std::cerr << "Wildcard string '" << wildcardString
                << "' returned error " << status << " from wordexp()."
                << std::endl;
      return {};
    }

    std::vector<std::string> fileList;

    for(unsigned int i = 0; i < p.we_wordc; ++i){
      // Check the file exists before adding it
      struct stat sb;
      if(stat(p.we_wordv[i], &sb) == 0)
        fileList.push_back(p.we_wordv[i]);
    }

    wordfree(&p);

    return fileList;
  }

  //----------------------------------------------------------------------
  std::string FindCAFAnaDir()
  {
    return std::string(getenv("MRB_SOURCE"))+"/sbnana/sbnana/CAFAna";
  }

  //----------------------------------------------------------------------
  std::vector<std::string> LoadFileList(const std::string& listfile)
  {
    std::vector<std::string> ret;

    std::ifstream is(listfile);
    if(!is.good()){
      std::cerr << "Can't open file list '" << listfile << "'. Aborting." << std::endl;
      abort();
    }

    while(!is.eof()){
      std::string fname;
      is >> fname;
      if(!fname.empty()) ret.push_back(fname);
    }
    return ret;
  }

  //----------------------------------------------------------------------
  std::map<std::string, std::string> GetCAFMetadata(TDirectory* dir)
  {
    std::map<std::string, std::string> ret;

    TTree* tr = (TTree*)dir->Get("metatree");
    if(!tr){
      std::cout << "Failed to find metadata tree in input CAF. Metadata will be blank." << std::endl;
      return ret;
    }

    std::string key, value;
    std::string* pkey = &key;
    std::string* pvalue = &value;
    tr->SetBranchAddress("key", &pkey);
    tr->SetBranchAddress("value", &pvalue);

    for(int i = 0; i < tr->GetEntries(); ++i){
      tr->GetEntry(i);
      ret[key] = value;
    }

    return ret;
  }

  //----------------------------------------------------------------------
  void CombineMetadata(std::map<std::string, std::string>& base,
                       const std::map<std::string, std::string>& add,
                       std::set<std::string>& mask)
  {
    for(auto it: add){
      const std::string& key = it.first;

      // Needs special handling anyway, skip it
      if(key == "parents") continue;

      // Accumulate the runs list
      if(key == "runs"){
        const std::string& r1 = base[key];
        const std::string& r2 = it.second;

        assert(!r2.empty());

        // Strip the outermost enclosing []
        const std::string& r2tmp = r2.substr(1, r2.length()-2);

        // Check that the run/subrun is not already present
        if (r1.find(r2tmp) != std::string::npos){
          std::cout << "Found files with duplicate run/subrun metadata: " << r2tmp << std::endl;

          abort();
        } 

        if(r1.empty()){
          base[key] = r2;
          continue;
        }

        // "[foo]" + "[bar]"
        std::string sum = r1+&r2[1]; // "[foo]bar]"
        sum[r1.size()-1] = ',';      // "[foo,bar]"
        base[key] = sum;
        continue;
      }

      if(base.find(key) == base.end()){
        // If it's new, add it
        base[key] = it.second;
      }
      else{
        if(key == "simulated.number_of_spills" ||
           key == "event_count" ||
           key == "online.totalevents"){
          // These two fields should be accumulated
          base[key] = TString::Format("%d",
                                      atoi(base[key].c_str()) +
                                      atoi(it.second.c_str())).Data();
        }
        else{
          // If it's a clash, record it
          if(base[key] != it.second) mask.insert(key);
        }
      }
    }
  }


  //----------------------------------------------------------------------
  void WriteCAFMetadata(TDirectory* dir,
                        const std::map<std::string, std::string>& meta)
  {
    TDirectory* tmp = gDirectory;
    dir->cd();

    TTree* trmeta = new TTree("metatree", "metatree");
    std::string key, value;
    trmeta->Branch("key", &key);
    trmeta->Branch("value", &value);
    for(const auto& keyval: meta){
      key = keyval.first;
      value = keyval.second;
      trmeta->Fill();
    }
    trmeta->Write();

    dir->Save();

    tmp->cd();
  }

  //----------------------------------------------------------------------
  bool RunningOnGrid()
  {
    static bool cache;
    static bool cache_set = false;
    if(!cache_set){
      cache = (getenv("_CONDOR_SCRATCH_DIR") != 0);
      cache_set = true;
    }

    return cache;
  }

  //----------------------------------------------------------------------
  size_t Stride(bool allow_default)
  {
    static int cache = -1;

    if(cache < 0){
      char* env = getenv("CAFANA_STRIDE");
      if(env){
        cache = std::atoi(env);
      }
      else{
        if(allow_default){
          cache = 1;
        }
        else{
          std::cout << "Stride() called, but CAFANA_STRIDE is not set (--stride not passed?)" << std::endl;
          abort();
        }
      }
    }

    return cache;
  }

  //----------------------------------------------------------------------
  size_t Offset(bool allow_default)
  {
    static int cache = -1;

    if(cache < 0){
      char* env = getenv("CAFANA_OFFSET");
      if(env){
        cache = std::atoi(env);
      }
      else{
        if(allow_default){
          cache = 0;
        }
        else{
          std::cout << "Offset() called, but CAFANA_OFFSET is not set (--offset not passed?)" << std::endl;
          abort();
        }
      }
    }

    return cache;
  }

  //----------------------------------------------------------------------
  int Limit()
  {
    static int cache = 0;

    if(cache == 0){
      char* env = getenv("CAFANA_LIMIT");
      if(env){
        cache = std::atoi(env);
      }
      else{
        cache = -1;
      }
    }

    return cache;
  }

  //----------------------------------------------------------------------
  size_t JobNumber()
  {
    if(!RunningOnGrid()){
      std::cout << "JobNumber() called, but we are not running on the grid" << std::endl;
      abort();
    }

    return Offset(false);
  }

  //----------------------------------------------------------------------
  size_t NumJobs()
  {
    if(!RunningOnGrid()){
      std::cout << "NumJobs() called, but we are not running on the grid" << std::endl;
      abort();
    }

    return Stride(false);
  }


  //----------------------------------------------------------------------
  bool AlmostEqual(double a, double b)
  {
    if(a == 0 && b == 0) return true;

    return fabs(a-b)/std::max(a, b) < .0001; // allow 0.01% error
  }


  //----------------------------------------------------------------------
  std::string pnfs2xrootd(std::string loc, bool unauth)
  {
    static bool first = true;
    static bool onsite = false;

    if (first && unauth) {
      first = false;
      char chostname[255];
      gethostname(chostname, 255);
      std::string hostname = chostname;

      if ( hostname.find("fnal.gov") != std::string::npos ){
        onsite = true;
        std::cout << "Using unauthenticated xrootd access (port 1095) while on-site, hostname: " << hostname << std::endl;
      }
      else {
        onsite = false;
        std::cout << "Using authenticated xrootd access (port 1094) access while off-site, hostname: " << hostname << std::endl;
      }
    }

    if(loc.rfind("/pnfs/", 0) == 0){ // ie begins with
      if ( onsite && unauth )
        loc = std::string("root://fndcagpvm01.fnal.gov:1095//pnfs/fnal.gov/usr/")+&loc.c_str()[6];
      else
        loc = std::string("root://fndcagpvm01.fnal.gov:1094//pnfs/fnal.gov/usr/")+&loc.c_str()[6];
    }
    return loc;
  }

  //----------------------------------------------------------------------
  FitToFourier::FitToFourier(TH1* h, double xlo, double xhi, int NOsc)
    : fHist(h), fxlo(xlo), fxhi(xhi), fNOsc(NOsc)
  {
  }

  //----------------------------------------------------------------------
  FitToFourier::~FitToFourier()
  {
  }

  //----------------------------------------------------------------------
  double FitToFourier::operator()(double *x, double *par) const
  {
    double x0 = x[0];
    double val = par[0];
    for (int i = 1; i <= fNOsc; i++)
      val += par[2*i-1]*sin(i*M_PI*x0) + par[2*i]*cos(i*M_PI*x0);
    return val;
  }

  //----------------------------------------------------------------------
  TF1* FitToFourier::Fit() const
  {
    //double s[fNOsc] = {0};
    //double c[fNOsc] = {0};

    std::vector<double> s(fNOsc, 0.0);
    std::vector<double> c(fNOsc, 0.0);

    int nBins = 0;
    for(int i = 1; i <= fHist->GetNbinsX(); ++i){
      const double x = M_PI * fHist->GetXaxis()->GetBinCenter(i);
      const double y = fHist->GetBinContent(i);

      if(y == 0) continue;
      ++nBins;

      for(int n = 0; n <= fNOsc; ++n){
        s[n] += y * sin(n*x);
        c[n] += y * cos(n*x);
      }
    }

    for(int n = 0; n <= fNOsc; ++n){
      s[n] *= 2./nBins;
      c[n] *= 2./nBins;
    }

    TF1* f = new TF1(UniqueName().c_str(), this, fxlo, fxhi, 2*fNOsc+1);

    f->SetParameter(0, c[0]/2);
    for(int n = 1; n <= fNOsc; ++n){
      f->SetParameter(n*2-1, s[n]);
      f->SetParameter(n*2,   c[n]);
    }

    // Because ROOT is having problems drawing f if I don't
    double min = fHist->GetMinimum();
    double max = fHist->GetMaximum();
    f->GetYaxis()->SetRangeUser(0.8*min, 1.2*max);
    return f;
  }

  //----------------------------------------------------------------------
  void EnsurePositiveDefinite(TH2* mat)
  {
    // Convert histogram to a proper matrix
    assert(mat->GetNbinsX() == mat->GetNbinsY());
    const int N = mat->GetNbinsX();
    TMatrixD m(N, N);
    for(int i = 0; i < N; ++i)
      for(int j = 0; j < N; ++j)
        m(i, j) = mat->GetBinContent(i+1, j+1);

    // Decompose it
    TVectorD evals;
    TMatrixD evecs = m.EigenVectors(evals);
    TMatrixD evalmat(N, N);
    // Force any negative eigenvalues slightly positive (floating point errors)
    for(int i = 0; i < N; ++i) evalmat(i, i) = std::max(1e-14, evals[i]);

    // Put the original matrix back together
    const TMatrixD evecs_inv(TMatrixD::kTransposed, evecs);
    m = evecs*evalmat*evecs_inv;

    // Decompose again to check for floating point problems
    m.EigenVectors(evals);
    for(int i = 0; i < N; ++i) assert(evals[i] > 0);

    // Copy the new matrix contents back into the histogram
    for(int i = 0; i < N; ++i)
      for(int j = 0; j < N; ++j)
        mat->SetBinContent(i+1, j+1, m(i, j));
  }

  //----------------------------------------------------------------------
  // Note that this does not work for 3D!
  TH1* GetMaskHist(const Spectrum& s, double xmin, double xmax, double ymin, double ymax)
  {
    if (s.GetBinnings().size() > 2){
      std::cout << "Error: unable to apply a mask in " << s.GetBinnings().size() << " dimensions" << std::endl;
      abort();
    }

    // The exposure isn't important here
    TH1* fMaskND  = s.ToTHX(s.POT());
    TH1D* fMask1D = s.ToTH1(s.POT());

    int ybins = fMaskND->GetNbinsY();

    for(int i = 0; i < fMask1D->GetNbinsX()+2; ++i){

      int ix = i / ybins;
      int iy = i % ybins;

      bool isMask = false;

      if (xmin < xmax){
	if (fMaskND->GetXaxis()->GetBinLowEdge(ix+1) < xmin) isMask=true;
	if (fMaskND->GetXaxis()->GetBinUpEdge(ix+1) > xmax) isMask=true;
      }

      if (ymin < ymax){
	if (fMaskND->GetYaxis()->GetBinLowEdge(iy+1) < ymin) isMask=true;
	if (fMaskND->GetYaxis()->GetBinUpEdge(iy+1) > ymax) isMask=true;
      }

      if (isMask) fMask1D->SetBinContent(i+1, 0);
      else fMask1D->SetBinContent(i+1, 1);

    }
    return fMask1D;
  }

  //----------------------------------------------------------------------
  double FindQuantile(double frac, std::vector<double>& xs)
  {
    // This turns out to be a much more fraught issue than you would naively
    // expect. This algorithm is equivalent to R-6 here:
    // https://en.wikipedia.org/wiki/Quantile#Estimating_quantiles_from_a_sample

    // In principle we could use std::nth_element(). Probably doesn't matter
    // much in practice since this is only for plotting.
    std::sort(xs.begin(), xs.end());

    const int N = xs.size();
    // The index we would ideally be sampling at
    const double h = frac*(N+1);
    // The indices on either side where we have to actually evaluate
    const unsigned int h0 = std::floor(h);
    const unsigned int h1 = std::ceil(h);
    if(h0 == 0) return xs[0]; // Don't underflow indexing
    if(h1 > xs.size()) return xs.back(); // Don't overflow indexing
    // The values at those indices
    const double x0 = xs[h0-1]; // wikipedia is using 1-based indexing
    const double x1 = xs[h1-1];

    if(h0 == h1) return x0;

    // Linear interpolation
    return (h1-h)*x0 + (h-h0)*x1;
  }
}
