#include "sbnana/CAFAna/Core/FileReducer.h"

#include "sbnana/CAFAna/Core/Progress.h"
#include "sbnana/CAFAna/Core/Utilities.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include "sbnanaobj/StandardRecord/StandardRecord.h"
#include "sbnanaobj/StandardRecord/SRGlobal.h"

#include <cassert>
#include <iostream>

#include <fenv.h>

#include "TChain.h"
#include "TFile.h"
#include "TH2.h"
#include "TTree.h"


namespace ana
{
  //----------------------------------------------------------------------
  void ClearTrueParticles(caf::StandardRecord* sr)
  {
    sr->true_particles.clear();
    sr->ntrue_particles = 0;
  }

  //----------------------------------------------------------------------
  FileReducer::FileReducer(const std::string& wildcard,
                           const std::string& outfile)
    : SpectrumLoaderBase(wildcard),
      fOutfile(outfile),
      fSpillCut(nullptr), fSliceCut(nullptr)
  {
  }

  //----------------------------------------------------------------------
  FileReducer::FileReducer(const std::vector<std::string>& fnames,
                           const std::string& outfile)
    : SpectrumLoaderBase(fnames),
      fOutfile(outfile),
      fSpillCut(nullptr),
      fSliceCut(nullptr)
  {
  }

  //----------------------------------------------------------------------
  FileReducer::~FileReducer()
  {
    delete fSpillCut;
    delete fSliceCut;
  }

  //----------------------------------------------------------------------
  void FileReducer::AddSpillCut(const SpillCut& cut)
  {
    if(fSpillCut){
      *fSpillCut = *fSpillCut && cut;
    }
    else{
      fSpillCut = new SpillCut(cut);
    }
  }

  //----------------------------------------------------------------------
  void FileReducer::AddSliceCut(const SliceCut& cut)
  {
    if(fSliceCut){
      *fSliceCut = *fSliceCut && cut;
    }
    else{
      fSliceCut = new SliceCut(cut);
    }
  }

  //----------------------------------------------------------------------
  void FileReducer::SetEventList(const std::string& fname)
  {
    FILE* f = fopen(fname.c_str(), "r");
    assert(f);

    while(!feof(f)){
      int run, subrun, event;
      fscanf(f, "%d %d %d", &run, &subrun, &event);
      fEventList.emplace(run, subrun, event);
    }

    fclose(f);
  }

  //----------------------------------------------------------------------
  void FileReducer::Go()
  {
    //    FloatingExceptionOnNaN fpnan;

    // Don't want overflow to happen. Set to 1 petabyte: effectively infinite.
    TTree::SetMaxTreeSize(1e15);

    if(fGone){
      std::cerr << "Error: can only call Go() once on a FileReducer" << std::endl;
      abort();
    }
    fGone = true;

    const int Nfiles = NFiles();

    Progress prog(TString::Format("Filling from %d files matching '%s'", Nfiles, fWildcard.c_str()).Data());

    TFile fout(fOutfile.c_str(), "RECREATE");
    TTree* trOut = 0; // created in first call to HandleFile()

    TH1* hPOTOut = new TH1F("TotalPOT", "", 1, 0, 1);
    TH1* hEventsOut = new TH1F("TotalEvents", "", 1, 0, 1);

    std::vector<std::string> fnames;

    std::map<std::string, std::string> meta;
    std::set<std::string> meta_mask;

    long nRecSeen = 0;
    long nRecPassed = 0;

    int fileIdx = -1;
    while(TFile* fin = GetNextFile()){
      ++fileIdx;

      assert(!fin->IsZombie());

      fnames.push_back(fin->GetName());

      TH1* hPOT = (TH1*)fin->Get("TotalPOT");
      assert(hPOT);
      hPOTOut->Add(hPOT);

      TDirectory* meta_dir = (TDirectory*)fin->Get("metadata");
      assert(meta_dir);
      CombineMetadata(meta, GetCAFMetadata(meta_dir), meta_mask);

      TH1* hEvents = (TH1*)fin->Get("TotalEvents");
      assert(hEvents);
      hEventsOut->Add(hEvents);

      HandleFile(fin, &fout, trOut,
                 Nfiles == 1 ? &prog : 0,
                 nRecSeen, nRecPassed);

      if(fileIdx == 0) CopyGlobalTree(fin, &fout);

      prog.SetProgress((fileIdx+1.)/Nfiles);
    } // end while GetNextFile

    fout.cd();
    trOut->Write();
    hPOTOut->Write();
    hEventsOut->Write();

    UpdateMetadata(meta, meta_mask, fnames);
    WriteCAFMetadata(fout.mkdir("metadata"), meta);

    fout.Close();

    prog.Done();

    std::cout << "Passed " << nRecPassed << " / " << nRecSeen << " records";
    std::cout << std::endl;
  }

  //----------------------------------------------------------------------
  void FileReducer::HandleFile(TFile* fin, TFile* fout, TTree*& trOut,
                               Progress* prog,
                               long& nRecSeen, long& nRecPassed)
  {
    TTree* recTree = (TTree*)fin->Get("recTree");
    assert(recTree);

    const caf::CAFType type = caf::GetCAFType(0, recTree);

    if(trOut){
      const caf::CAFType outtype = caf::GetCAFType(0, trOut);
      if(type != outtype){
        std::cerr << "FileReducer: Error: dataset contains mixed CAF types (flat vs nested)" << std::endl;
        abort();
      }
    }

    if(type == caf::kNested) HandleNestedTree(fout, recTree, trOut,
                                              prog,
                                              nRecSeen, nRecPassed);
    else if(type == caf::kFlatSingleTree) HandleFlatTree(fout, recTree, trOut,
                                                         prog,
                                                         nRecSeen, nRecPassed);

    else{
      std::cerr << "FileReducer: Error: Unrecognized file type: "
                << fin->GetName() << std::endl;
      abort();
    }
  }

  //----------------------------------------------------------------------
  void FileReducer::HandleNestedTree(TFile* fout,
                                     TTree* recTree, TTree*& trOut,
                                     Progress* prog,
                                     long& nRecSeen, long& nRecPassed)
  {
    if(!trOut){
      fout->cd();
      trOut = new TTree("recTree", "recTree");
      //      FloatingExceptionOnNaN fpnan(false);
      caf::StandardRecord dummy;
      trOut->Branch("rec", &dummy);
    }

    // Use this one for assessing cuts etc
    caf::Proxy<caf::StandardRecord> srProxy(0, recTree, "rec", 0, 0);

    // And if they pass load into this one for writing out
    caf::StandardRecord* sr = 0;
    recTree->SetBranchAddress("rec", &sr);

    caf::StandardRecord* oldsr = 0;

    const int Nentries = recTree->GetEntries();
    for(int n = 0; n < Nentries; ++n){
      ++nRecSeen;
      recTree->LoadTree(n);

      // Apply EventList cut if it's been enabled
      if(!fEventList.empty() &&
         !fEventList.count(std::make_tuple(srProxy.hdr.run,
                                           srProxy.hdr.subrun,
                                           srProxy.hdr.evt))) continue;

      /// Do we need to include the event? Either based on the selection...
      const bool passesSpillCut = !fSpillCut || (*fSpillCut)(&srProxy);
      // Or if it carries POT or spill counting information
      if(passesSpillCut || srProxy.hdr.pot > 0 || srProxy.hdr.ngenevt > 0){
        recTree->GetEntry(n);

        if(sr != oldsr){
          trOut->SetBranchAddress("rec", &sr);
          oldsr = sr;
        }

        if(!passesSpillCut){
          // We're only keeping this one for the exposure info
          Huskify(sr);
          trOut->Fill();
          continue;
        }
        // Otherwise, need to proceed to keep/reject individual slices

        std::vector<int> tocut;
        for(unsigned int i = 0; i < srProxy.slc.size(); ++i){
          if(fSliceCut && !(*fSliceCut)(&srProxy.slc[i])) tocut.push_back(i);
        }

        // Remove slices in reverse order so that the indices remain valid
        for(auto it = tocut.rbegin(); it != tocut.rend(); ++it){
          sr->slc.erase(sr->slc.begin() + *it);
        }

        // Apply any additional reduction steps
        for(const auto& f: fReductionFuncs) f(sr);
        // This is kind of problematic since the proxy and actual record
        // could be out of sync. Let's just disable this option for now.
        //          for(const auto & f: fReductionFuncsWithProxy) f(sr, &srProxy);

        ++nRecPassed;
        trOut->Fill();
      }

      // prog won't be passed if the caller doesn't want per-event updates
      if(n%100 == 0 && prog) prog->SetProgress(double(n)/Nentries);
    } // end for n
  }

  //----------------------------------------------------------------------
  void FileReducer::HandleFlatTree(TFile* fout, TTree* recTree, TTree*& trOut,
                                   Progress* prog,
                                   long& nRecSeen, long& nRecPassed)
  {
    if(!fEventList.empty()){
      std::cerr << "FileReducer: Event list not supported for FlatCAFs (yet)" << std::endl;
      abort();
    }

    if(fSpillCut){
      std::cerr << "FileReducer: Spill cuts not supported for FlatCAFs" << std::endl;
      abort();
    }

    if(fSliceCut){
      std::cerr << "FileReducer: Slice cuts not supported for FlatCAFs" << std::endl;
      abort();
    }

    if(!fReductionFuncs.empty()){
      std::cerr << "FileReducer: Reduction functions not supported for FlatCAFs" << std::endl;
      abort();
    }

    // Do the actual copy
    if(!trOut){
      fout->cd();
      trOut = recTree->CloneTree();
    }
    else{
      fout->cd();
      recTree->CopyAddresses(trOut);
      trOut->CopyEntries(recTree);
    }

    nRecSeen += recTree->GetEntries();
    nRecPassed += recTree->GetEntries();
  }

  //----------------------------------------------------------------------
  void FileReducer::CopyGlobalTree(TFile* fin, TFile* fout)
  {
    TTree* globalIn = (TTree*)fin->Get("globalTree");
    if(!globalIn) return;

    // Copy globalTree verbatim from input to output
    caf::SRGlobal global;
    caf::SRGlobal* pglobal = &global;
    globalIn->SetBranchAddress("global", &pglobal);
    fout->cd();
    TTree globalOut("globalTree", "globalTree");
    globalOut.Branch("global", "caf::SRGlobal", &global);
    assert(globalIn->GetEntries() == 1);
    // TODO check that the globalTree is the same in every file
    globalIn->GetEntry(0);
    globalOut.Fill();
    globalOut.Write();
  }

  //----------------------------------------------------------------------
  void FileReducer::UpdateMetadata(std::map<std::string, std::string>& meta,
                                   const std::set<std::string>& mask,
                                   const std::vector<std::string>& fnames) const
  {
    for(const std::string& m: mask){
      std::cerr << "Warning: metadata parameter '" << m << "' differs between input files and has been dropped from the output." << std::endl;
      meta.erase(m);
    }

    /*
    // change caf -> decaf in the metadata field, if parents have data_tier
    if(meta.find("data_tier") != meta.end()){
      std::string decaf_tier = meta["data_tier"];
      assert(decaf_tier.size() >= 3);
      // For indexing: note that all of this is surrounded by quotes
      assert(decaf_tier.substr(decaf_tier.size()-4,3) == "caf");
      // don't make 'decaf' into 'dedecaf', however
      if (decaf_tier.size() < 6 || decaf_tier.substr(decaf_tier.size()-6,5) != "decaf")
        decaf_tier.replace(decaf_tier.size()-4,3,"decaf");
      meta["data_tier"] = decaf_tier;
    }
    */

    std::string parents = "[";
    for(const std::string& f: fnames){
      parents += "{\"file_name\": \""+std::string(basename((char *)f.c_str()))+"\"}, ";
    }
    if(!fnames.empty()) parents.resize(parents.size()-2); // drop trailing ", "
    parents += "]";

    meta["parents"] = parents;

    // if there's one more than one parent, this is a concat
    if(fnames.size() > 1 && meta.find("data_tier") != meta.end()){
      meta["data_tier"] = "\"concat_caf\"";
    }

    // Allow user to override any metadata
    for(auto it: fMetaMap) meta[it.first] = it.second;
  }

  //----------------------------------------------------------------------
  void FileReducer::Huskify(caf::StandardRecord* sr) const
  {
    sr->hdr.husk = true;

    // It's not actually necessary for correctness to clear all these fields
    // out, so not a disaster if we happen to miss one that gets added in
    // future, for example. But it does help to minimize the size of the
    // resulting CAF.

    sr->slc.clear();
    sr->nslc = 0;
    sr->fake_reco.clear();
    sr->nfake_reco = 0;
    sr->true_particles.clear();
    sr->ntrue_particles = 0;
    sr->crt_hits.clear();
    sr->ncrt_hits = 0;
    sr->crt_tracks.clear();
    sr->ncrt_tracks = 0;

    sr->reco.trk.clear();
    sr->reco.ntrk = 0;
    sr->reco.shw.clear();
    sr->reco.nshw = 0;
    sr->reco.stub.clear();
    sr->reco.nstub = 0;

    sr->mc.nu.clear();
    sr->mc.nnu = 0;
    sr->mc.prtl.clear();
    sr->mc.nprtl = 0;
  }
} // namespace
