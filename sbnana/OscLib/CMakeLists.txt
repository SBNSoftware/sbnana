include_directories( $ENV{GSL_INC} $ENV{BOOST_INC} )

# Want to be able to find files using just an OscLib/ prefix
include_directories(..)

cet_find_library( GSL      NAMES gsl      PATHS ENV GSL_LIB NO_DEFAULT_PATH )
cet_find_library( GSLCBLAS NAMES gslcblas PATHS ENV GSL_LIB NO_DEFAULT_PATH )

cet_make( LIBRARY_NAME  OscLib
          LIBRARIES ${ROOT_BASIC_LIB_LIST}
                    ${GSL}
	            ${GSLCBLAS}
        )

install_headers()
install_source()
