
file(GLOB lib_srcs *.cxx)

cet_make_library(LIBRARY_NAME CAFAnaCore
         SOURCE ${lib_srcs}
         LIBRARIES sbnanaobj::StandardRecordProxy
                   sbnanaobj::StandardRecord
                   # TODO how am I supposed to link to an external library?
#                  ${OSCLIB}
                   -L$ENV{OSCLIB_LIB}
                   OscLib
                   Threads::Threads
                   ROOT::Core
                   ifdh::ifdh
                   Eigen3::Eigen
                   Boost::system
        )


# Allow "StandardRecord/" to find the right directory
target_include_directories(CAFAnaCore PRIVATE  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>)

target_include_directories(CAFAnaCore PRIVATE $ENV{OSCLIB_INC})

target_include_directories(CAFAnaCore PRIVATE $ENV{SRPROXY_INC})

install_headers()
install_source()
