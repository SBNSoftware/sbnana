# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME SBNAnaCuts
         LIBRARIES CAFAnaCore
                   SBNAnaVars
                   StandardRecordProxy
        )

install_headers()
install_source()
