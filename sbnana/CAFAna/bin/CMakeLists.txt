cet_make_exec( diff_spectra
               SOURCE diff_spectra.cc
               LIBRARIES ${ROOT_BASIC_LIB_LIST}
               )

cet_make_exec( concat_cafs
               SOURCE concat_cafs.cc
               LIBRARIES ${ROOT_BASIC_LIB_LIST} CAFAnaCore
               )

install_headers()
install_source()
