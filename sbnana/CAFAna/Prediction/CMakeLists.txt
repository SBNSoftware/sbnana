
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME CAFAnaPrediction
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore CAFAnaSysts CAFAnaExtrap
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(CAFAnaPrediction PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_include_directories(CAFAnaPrediction PRIVATE $ENV{OSCLIB_INC})

install_headers()
install_source()
