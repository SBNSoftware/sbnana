
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME SBNAnaCuts
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
                   CAFAnaSysts
                   SBNAnaVars
                   sbnanaobj::StandardRecordProxy
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(SBNAnaCuts PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_include_directories(SBNAnaCuts PRIVATE $ENV{SRPROXY_INC})

install_headers()
install_source()
