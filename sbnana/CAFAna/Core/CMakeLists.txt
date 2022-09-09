include_directories( $ENV{IFDHC_FQ_DIR}/inc )

include_directories($ENV{OSCLIB_INC})

cet_make(LIBRARY_NAME CAFAnaCore
         LIBRARIES sbnanaobj_StandardRecordProxy
                   sbnanaobj_StandardRecord
                   # TODO how am I supposed to link to an external library?
#                  ${OSCLIB}
                   -L$ENV{OSCLIB_LIB}
                   OscLib
                  -L$ENV{CAFANACORE_LIB}
                   CAFAnaCoreExt
                   pthread
                   ${ROOT_BASIC_LIB_LIST}
                   ifdh::ifdh
        )

install_headers()
install_source()
