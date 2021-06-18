# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME CAFAnaSysts
         LIBRARIES CAFAnaCore
                   CAFAnaCuts
                   StandardRecordProxy
        )

install_headers()
install_source()
