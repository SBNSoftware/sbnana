
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME SBNAnaVars
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
                   sbnanaobj::StandardRecordProxy
                   -L$ENV{OSCLIB_LIB}
                   cetlib::cetlib
                   cetlib_except::cetlib_except
                   OscLib
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(SBNAnaVars PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_include_directories(SBNAnaVars PRIVATE $ENV{OSCLIB_INC})

target_include_directories(SBNAnaVars PRIVATE $ENV{SRPROXY_INC})

install_headers()
install_source()
