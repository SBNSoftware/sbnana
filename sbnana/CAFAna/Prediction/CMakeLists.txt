# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME CAFAnaPrediction
         LIBRARIES CAFAnaCore CAFAnaSysts CAFAnaExtrap
        )

install_headers()
install_source()
