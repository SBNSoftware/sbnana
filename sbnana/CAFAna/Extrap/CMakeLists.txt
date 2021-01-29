# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME CAFAnaExtrap
         LIBRARIES CAFAnaCore
        )

install_headers()
install_source()
