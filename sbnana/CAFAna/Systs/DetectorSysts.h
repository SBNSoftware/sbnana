#pragma once

#include "sbnana/CAFAna/Core/ISyst.h"

#include <string>
#include <vector>

class TH1;

namespace ana {

    class DetectorSyst : public ISyst {
    public:
        void Shift(double sigma, caf::SRSliceProxy* slc, double& weight) const override;
    protected:
        friend const DetectorSyst* GetDetectorSyst(const std::string&,
                                               const std::string&,
                                               const std::string&,
                                               const std::string&);
        DetectorSyst(const std::string& dir,
                     const std::string& prefix,
                     const std::string& name,
                     const std::string& systFile = "");
        std::string fHistName, fName, fSystFilePath;
    };

    const DetectorSyst* GetDetectorSyst(const std::string& dir,
                                        const std::string& prefix,
                                        const std::string& name,
                                        const std::string& systFile = "");

    std::vector<const ISyst*> GetDetectorSysts();
}