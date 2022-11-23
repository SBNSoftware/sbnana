file(GLOB lib_srcs *.cxx)

cet_make_library(LIBRARY_NAME CAFAnaCuts
         SOURCE ${lib_srcs}
         LIBRARIES CAFAnaCore
                   sbnanaobj::StandardRecordProxy
        )

# Allow "StandardRecord/" to find the right directory
target_include_directories(CAFAnaCuts PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_include_directories(CAFAnaCuts PRIVATE $ENV{OSCLIB_INC})

target_include_directories(CAFAnaCuts PRIVATE $ENV{SRPROXY_INC})

install_headers()
install_source()
