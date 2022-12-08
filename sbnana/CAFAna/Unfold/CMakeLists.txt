
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME CAFAnaUnfold
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(CAFAnaUnfold PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

install_headers()
install_source()
