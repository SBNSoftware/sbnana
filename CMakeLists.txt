# ======================================================================
#  larsoft main build file
#
#  cd .../path/to/build/directory
#  source .../path/to/larsoft/ups/setup_for_development <-d|-p>
#  cmake [-DCMAKE_INSTALL_PREFIX=/install/path]
#        -DCMAKE_BUILD_TYPE=$CETPKG_TYPE
#        $CETPKG_SOURCE
#  make
#  make test
#  make install
#  make package (builds distribution tarfile)
# ======================================================================


cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
# cmake_policy(VERSION 3.19)


find_package(cetmodules )
project(sbnana VERSION 09.78.06 LANGUAGES CXX)

# cetbuildtools contains our cmake modules

message(STATUS "\n\n  ==========================   ${PROJECT_NAME}   ==========================")

include(CetCMakeEnv)
cet_cmake_env()

# The specification of the C++ standard is required because CET macros
# set it to C++98 by default, and then override it depending on the
# compiler version, as determined by the compiler qualifier.
# Removing that qualifier, we need to supply that information manually,
# and maintain it too.
cet_set_compiler_flags(DIAGS CAUTIOUS
  WERROR
  NO_UNDEFINED
  ALLOW_DEPRECATIONS
  EXTRA_FLAGS -pedantic -Wno-unused-local-typedefs -Wno-undefined-var-template -I $ENV{CLANG_FQ_DIR}/include/c++/v1
)

cet_report_compiler_flags()

# these are minimum required versions, not the actual product versions
find_package( sbnanaobj REQUIRED )
find_package( Boost COMPONENTS system filesystem REQUIRED )
find_package( ROOT REQUIRED )
find_package( ifdhc REQUIRED)
find_package( srproxy )
find_package( osclib  )
find_package( pygccxml  )
find_package( Threads REQUIRED )
find_package( Eigen3 REQUIRED )
find_package( cetlib REQUIRED EXPORT )
find_package( cetlib_except REQUIRED EXPORT )

#ifdhc
cet_find_library( IFDH    NAMES ifdh    PATHS ENV IFDHC_LIB   NO_DEFAULT_PATH )

add_subdirectory(sbnana)

# packaging utility
cet_cmake_config()
