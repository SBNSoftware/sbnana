# Allow "StandardRecord/" to find the right directory
include_directories(..)

cet_make(LIBRARY_NAME CAFAnaAnalysis
         LIBRARIES CAFAnaCore CAFAnaExperiment
                   ${ROOT_MINUIT2}
        )

install_headers()
install_source()
