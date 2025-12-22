// -*- mode: c++; c-basic-offset: 2; -*-
void load(std::string lib)
{
  std::cout << "." << std::flush;
  int ret = gSystem->Load(("lib"+lib).c_str());
  // In case of error, exit immediately with the error clearly showing, instead
  // of with a confusing secondary error about a page of output later.
  if(ret != 0){
    std::cout << std::endl << "gSystem->Load(lib"+lib+") failed with code " << ret << std::endl;
    exit(ret);
  }
}

void load_cafana_libs()
{
  // All the CINT exception handler does is obfuscate the stack. With this,
  // uncaught exceptions immediately show a useful backtrace under gdb.
  //  G__SetCatchException(0);

  TString qmrb = gSystem->Getenv("MRB_QUALS");
  // Mirror the optimization settings in use elsewhere
  if(qmrb.Contains("debug")) {
    gSystem->SetAclicMode(TSystem::kDebug);
  }
  else{
    gSystem->SetAclicMode(TSystem::kOpt);
  }

  // This magic incantation prevents ROOT doing slow cleanup work in
  // TList::RecursiveRemove() under ~TH1(). It also tends to lead to shutdown
  // crashes. This seems like a good compromise: go fast in batch mode
  // (probably fitting) and slow but not crashy in interactive move (probably
  // want to see the plots).
  if(gROOT->IsBatch()) gROOT->SetMustClean(false);

  // Colorize error messages. Would be better if we could just pick up the
  // flags novasoft uses, but they don't seem to be in any env var.
  gSystem->SetFlagsDebug(TString(gSystem->GetFlagsDebug())+" -fdiagnostics-color=auto");
  gSystem->SetFlagsOpt(TString(gSystem->GetFlagsOpt())+" -fdiagnostics-color=auto -UNDEBUG"); // match gcc's maxopt behaviour of retaining assert()

  std::string incdir;
  if(getenv("SBNANA_INC")){
    // ups product
    incdir = std::string(getenv("SBNANA_INC"));
  }
  else{
    // local mrb checkout
    char* mrbi = getenv("MRB_INSTALL");
    if(!mrbi){
      std::cout << "$MRB_INSTALL is not set" << std::endl;
      exit(1);
    }
    char* sbnv = getenv("SBNANA_VERSION");
    if(!sbnv){
      std::cout << "$SBNANA_VERSION is not set" << std::endl;
      exit(1);
    }

    incdir = std::string(mrbi)+"/sbnana/"+std::string(sbnv)+"/include/";
  }

  TString includes = "-I"+incdir+" -I$ROOTSYS/include -I$SRPROXY_INC -I$OSCLIB_INC -I$EIGEN_INC";

  // List of libraries to load. Dependency order.
  const std::vector<std::string> libs =
    {
      "Minuit2",
      //      "StandardRecord",
      "sbnanaobj_StandardRecordProxy",
      // "StandardRecord_dict",
      "OscLib",
      "CAFAnaCore",
      "CAFAnaUnfold",
      "CAFAnaVars",
      "CAFAnaCuts",
      "SBNAnaVars",
      "SBNAnaCuts",
      "CAFAnaSysts",
      "CAFAnaExtrap",
      "CAFAnaPrediction",
      "CAFAnaExperiment",
      "CAFAnaAnalysis",
    };

  // Actually load the libraries
  std::cout << "Loading libraries";
  for(const std::string& lib: libs) load(lib);
  std::cout << std::endl;


  gSystem->SetIncludePath(includes);

  // Doesn't seem to work
  //  gSystem->Setenv("IFDH_DEBUG", "0"); // shut ifdh up

  // Pick up standard style
  if(getenv("SBNANA_FQ_DIR")){
    gROOT->Macro("${SBNANA_FQ_DIR}/bin/rootlogon.C");
  }
  else{
    gROOT->Macro("${SBNANA_DIR}/sbnana/CAFAna/rootlogon.C");
  }
  gROOT->ForceStyle();

  TRint* rint = dynamic_cast<TRint*>(gApplication);
  if(rint) rint->SetPrompt("cafe [%d] ");


  // Do this last so that the profiler library is unloaded before our main
  // libraries, meaning that the code in ProfilerSupport runs at the right time
  // (after the profile file is made). This does mean that we aren't able to
  // profile library loading/startup. The ideal solution would be to make a
  // mini-library for just ProfilerSupport.
  if(gSystem->Getenv("CPUPROFILE")){
    // We were passed the --prof option
    const std::string cpuprof = gSystem->Getenv("CPUPROFILE");

    std::cout << "Profiling enabled." << std::endl;

    if(!qmrb.Contains("debug")){
      std::cout << "Note: profiling works much better in debug mode." << std::endl;
    }

    const char* pd = getenv("GPERFTOOLS_DIR");
    if(pd){
      gSystem->Load((std::string(pd)+"/lib/libprofiler.so").c_str());
      // Apparently the library load manages to corrupt the env var? Put it
      // back.
      gSystem->Setenv("CPUPROFILE", cpuprof.c_str());
    }
    else{
      std::cout << "Couldn't find gperftools library" << std::endl;
    }
  }
}
