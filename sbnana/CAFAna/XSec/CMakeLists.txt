
file(GLOB lib_srcs *.cxx)
cet_make_library(LIBRARY_NAME CAFAnaXSec
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
        )

target_include_directories( CAFAnaXSec PRIVATE $ENV{SRPROXY_INC} )

install_headers()
install_source()
