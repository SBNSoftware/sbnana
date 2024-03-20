
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME CAFAnaSysts
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
                   CAFAnaCuts
                   sbnanaobj::StandardRecordProxy
                   sbnanaobj::StandardRecord
                   cetlib::cetlib
                   cetlib_except::cetlib_except
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(CAFAnaSysts PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_include_directories(CAFAnaSysts PRIVATE $ENV{SRPROXY_INC})

install_headers()
install_source()
