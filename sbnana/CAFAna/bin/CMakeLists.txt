cet_make_exec( NAME diff_spectra
               SOURCE diff_spectra.cc
               LIBRARIES ROOT::Core ROOT::TreePlayer
               )

cet_make_exec( NAME concat_cafs
               SOURCE concat_cafs.cc
               LIBRARIES ROOT::Core CAFAnaCore
               )

target_include_directories( concat_cafs PRIVATE $ENV{SRPROXY_INC} )

install_headers()
install_source()
