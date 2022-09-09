#include "sbnana/CAFAna/Core/Multiverse.h"

#include "sbnana/CAFAna/Systs/SBNWeightSysts.h"

#include "TFile.h"

#include <iostream>

using namespace ana;

void test_multiverse(bool remake = true)
{
  if(remake){
    const std::vector<const ISyst*> systs = GetSBNGenieWeightSysts();

    const Multiverse& cross = Multiverse::Hypercross(systs, 3);
    const Multiverse& gas = Multiverse::RandomGas(systs, 100, 42);//Multiverse::kTrulyRandom);

    std::cout << cross.LatexName() << std::endl << cross.ShortName() << std::endl;
    std::cout << gas.LatexName() << std::endl << gas.ShortName() << std::endl;

    TFile fout("test_multiverse.root", "RECREATE");
    cross.SaveTo(&fout, "cross");
    gas.SaveTo(&fout, "gas");
  }

  std::cout << std::endl;

  GetSBNGenieWeightSysts(); // magic them into existence

  TFile fin("test_multiverse.root");
  const FitMultiverse& cross = *Multiverse::LoadFrom(&fin, "cross");
  const FitMultiverse& gas = *Multiverse::LoadFrom(&fin, "gas");

  std::cout << cross.LatexName() << std::endl << cross.ShortName() << std::endl;
  std::cout << gas.LatexName() << std::endl << gas.ShortName() << std::endl;
}
