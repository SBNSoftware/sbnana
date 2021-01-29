# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME CAFAnaVars
         LIBRARIES CAFAnaCore
                   StandardRecordProxy
                   OscLib
        )

install_headers()
install_source()
