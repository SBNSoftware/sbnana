include_directories( $ENV{IFDHC_FQ_DIR}/inc )

# Allow "StandardRecord/" to find the right directory
include_directories(..)

include_directories($ENV{OSCLIB_INC})

cet_make(LIBRARY_NAME CAFAnaCore
         LIBRARIES StandardRecordProxy
                   sbnanaobj_StandardRecord
                   # TODO how am I supposed to link to an external library?
#                  ${OSCLIB}
                   -L$ENV{OSCLIB_LIB}
                   OscLib
                   pthread
                   ${ROOT_BASIC_LIB_LIST}
                   ${IFDH}
        )

install_headers()
install_source()
