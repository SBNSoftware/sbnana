#include "sbnana/CAFAna/Core/FileReducer.h"
using namespace ana;

#include <iostream>
#include <string>
#include <vector>

template<class... T> void concat_cafs(T... args)
{
  std::vector<std::string> infiles = {args...};

  if(infiles.size() < 2){
    std::cout << "usage: cafe -bq concat_cafs INFILES    OUTFILE\n"
              << "       cafe -bq concat_cafs INFILES... OUTFILE\n"
              << "       cafe -bq concat_cafs 'WILDCARD' OUTFILE\n"
              << "       cafe -bq concat_cafs SAMDEF     OUTFILE\n"
              << "Warning: this argument order is different to hadd"
              << std::endl;
    exit(1);
  }

  const std::string outfile = infiles.back();
  infiles.pop_back();

  FileReducer* reducer = 0;
  if(infiles.size() == 1){
    // Could be an (escaped) wildcard, or SAM definition
    reducer = new FileReducer(infiles[0], outfile);
  }
  else{
    reducer = new FileReducer(infiles, outfile);
  }

  reducer->SetMetadata("data_tier", "concat_caf");
  reducer->SetMetadata("file_format", "concat_caf");

  reducer->Go();

  delete reducer;

  std::cout << "Created " << outfile << std::endl;
}
