
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME CAFAnaExperiment
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore CAFAnaVars CAFAnaPrediction
                   -L$ENV{OSCLIB_LIB}
                   OscLib
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(CAFAnaExperiment PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_include_directories(CAFAnaExperiment PRIVATE $ENV{OSCLIB_INC})

install_headers()
install_source()
