# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME CAFAnaCuts
         LIBRARIES CAFAnaCore
                   StandardRecordProxy
        )

install_headers()
install_source()
