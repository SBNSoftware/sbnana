# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME CAFAnaExperiment
         LIBRARIES CAFAnaCore CAFAnaVars CAFAnaPrediction
                   OscLib
        )

install_headers()
install_source()
