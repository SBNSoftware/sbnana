#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"
#include "sbnana/CAFAna/Core/Utilities.h"
#include "sbnana/CAFAna/Core/MultiVar.h"
#include "sbnana/CAFAna/Core/Var.h"

#include "sbnanaobj/StandardRecord/Proxy/SRProxy.h"

#include <string>
#include <vector>

class TH1;

namespace ana {

    class DetectorSyst : public ISyst {
    public:
        void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override;
        void Shift(double sigma, caf::SRTrueInteractionProxy* slc, double& weight) const override {
          std::cout << "Error: Detector systematics not implmented for true interactions." << std::endl;
          abort();
        }
    protected:
        friend const DetectorSyst* GetDetectorSyst(const std::string&,
                                               std::string&,
                                               const std::string&,
                                               std::string&,
                                               Var,
                                               int,
                                               const std::string&
                                               );
        DetectorSyst(const std::string& dir,
                     std::string& prefix,
                     const std::string& name,
                     std::string& variable,
                     Var var,
                     int fIdx,
                     const std::string& systFile = "");
        std::string fHistName, fName, fVariable, fSystFilePath;
        Var fVar;
        int fDualSided;
    };

    const DetectorSyst* GetDetectorSyst(const std::string& dir,
                                        std::string& prefix,
                                        const std::string& name,
                                        std::string& variable,
                                        Var var,
                                        int fIdx,
                                        const std::string& systFile = "");

    std::vector<const ISyst*> GetDetectorSysts(std::string name_in, Var var_in);
}
