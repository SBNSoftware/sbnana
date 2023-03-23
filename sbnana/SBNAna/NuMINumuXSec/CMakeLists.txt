file(GLOB lib_srcs *.cxx)

find_package(cetlib)
find_package(cetlib_except)

cet_make_library(LIBRARY_NAME NuMINumuXSec
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
                   CAFAnaSysts
                   SBNAnaVars
                   SBNAnaCuts
                   CRTPMTMatching
                   sbnanaobj::StandardRecordProxy
                   cetlib::cetlib
                   cetlib_except::cetlib_except
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(NuMINumuXSec PRIVATE $ENV{SRPROXY_INC})

install_headers()
install_source()
