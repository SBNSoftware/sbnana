cet_make(LIBRARY_NAME CAFAnaAnalysis
         LIBRARIES CAFAnaCore CAFAnaExperiment
                   ${ROOT_MINUIT2}
        )

install_headers()
install_source()
