


#include "sbnana/CAFAna/Core/Binning.h"
#include "sbnana/CAFAna/Core/Spectrum.h"
#include "sbnana/CAFAna/Core/SpectrumLoader.h"

#include "training_vars.h"

using namespace ana;

#include "TCanvas.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TLegend.h"
#include <numeric>
#include <fstream>
#include <TROOT.h>
#include <TStyle.h>


// ---- VARS -----
std::vector<float> starttoend,endtostart,starttostart,endtoend,angles,matches,energies,slices,intersections;

const Var kFindSliceInfo([](const caf::SRSliceProxy* slc) -> int {
  std::vector<double> vStarttoend =  kStarttoend(slc);
  starttoend.insert(starttoend.end(), vStarttoend.begin(), vStarttoend.end());
  std::vector<double> vEndtostart =  kEndtostart(slc);
  endtostart.insert(endtostart.end(), vEndtostart.begin(), vEndtostart.end());
  std::vector<double> vStarttostart =  kStarttostart(slc);
  starttostart.insert(starttostart.end(), vStarttostart.begin(), vStarttostart.end());
  std::vector<double> vEndtoend =  kEndtoend(slc);
  endtoend.insert(endtoend.end(), vEndtoend.begin(), vEndtoend.end());
  std::vector<double> vAngles =  kAngles(slc);
  angles.insert(angles.end(), vAngles.begin(), vAngles.end());
  std::vector<double> vMatches =  kMatches(slc);
  matches.insert(matches.end(), vMatches.begin(), vMatches.end());
  std::vector<double> vEnergies =  kEnergies(slc);
  energies.insert(energies.end(), vEnergies.begin(), vEnergies.end());
  std::vector<double> vSlices =  kSlices(slc);
  slices.insert(slices.end(), vSlices.begin(), vSlices.end());
  std::vector<double> vIntersections =  kIntersections(slc);
  intersections.insert(intersections.end(), vIntersections.begin(), vIntersections.end());

  return 42;
});

void make_training_texts(const std::string inputName="nue_cosmic_recursiveshower.root")
{

  SpectrumLoader loader(inputName);

  const Binning bins = Binning::Simple(100, 0, 5000);
  Spectrum sFindSlice("", bins, loader, kFindSliceInfo, kNoSpillCut, kNoCut);

  loader.Go();
  ofstream output;
  output.open("training_text.txt");
  output << "Event,Energy,Nu Slice,Start to End,End to Start,Start to Start,End to End,Angle,Overlap,Match";
  output << "\n";
  for (unsigned int i = 0; i<matches.size(); i++) {
    output << i << "," << energies[i] << "," << slices[i] << "," << starttoend[i] << "," << endtostart[i] << "," << starttostart[i] << "," << endtoend[i] << "," << angles[i] << "," <<intersections[i]<<","<< matches[i];
    output << "\n";
  }
  output.close(); 
}
