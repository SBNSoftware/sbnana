#include "sbnana/CAFAna/Core/FileReducer.h"
using namespace ana;

#include <iostream>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
  if(argc < 3){
    std::cout << "usage: concat_cafs INFILES    OUTFILE\n"
              << "       concat_cafs INFILES... OUTFILE\n"
              << "       concat_cafs 'WILDCARD' OUTFILE\n"
              << "       concat_cafs SAMDEF     OUTFILE\n"
              << "Warning: this argument order is different to hadd"
              << std::endl;
    exit(1);
  }

  std::vector<std::string> infiles;
  for(int i = 1; i < argc-1; ++i) infiles.push_back(argv[i]);
  const std::string outfile = argv[argc-1];

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
