file(GLOB lib_srcs *.cxx)

cet_make_library(LIBRARY_NAME CRTPMTMatching
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
                   CAFAnaSysts
                   SBNAnaVars
                   SBNAnaCuts
                   sbnanaobj::StandardRecordProxy
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(CRTPMTMatching PRIVATE $ENV{SRPROXY_INC})

install_headers()
install_source()
