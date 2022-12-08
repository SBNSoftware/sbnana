
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME CAFAnaExtrap
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore CAFAnaCuts
        )

# Allow "StandardRecord/" to find the right directory
include_directories(CAFAnaExtrap PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

install_headers()
install_source()
