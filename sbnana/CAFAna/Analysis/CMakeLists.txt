
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME CAFAnaAnalysis
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore CAFAnaExperiment
                   ROOT::Minuit2
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(CAFAnaAnalysis PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_include_directories(CAFAnaAnalysis PRIVATE $ENV{OSCLIB_INC})

install_headers()
install_source()
