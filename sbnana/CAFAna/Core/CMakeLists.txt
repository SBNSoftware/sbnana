include_directories( $ENV{IFDHC_FQ_DIR}/inc )

# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME CAFAnaCore
         LIBRARIES StandardRecordProxy
                   StandardRecord
                   OscLib
                   pthread
                   ${ROOT_BASIC_LIB_LIST}
                   ${IFDH}
        )

install_headers()
install_source()
