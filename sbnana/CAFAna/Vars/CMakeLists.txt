
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME CAFAnaVars
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
                   sbnanaobj::StandardRecordProxy
                   -L$ENV{OSCLIB_LIB}
                   OscLib
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories( CAFAnaVars PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..> )

target_include_directories( CAFAnaVars PRIVATE $ENV{OSCLIB_INC} )

target_include_directories( CAFAnaVars PRIVATE $ENV{SRPROXY_INC} )

install_headers()
install_source()
