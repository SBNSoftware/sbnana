cet_make_exec(NAME flatten_caf
              SOURCE flatten_caf.cc
              LIBRARIES sbnanaobj::StandardRecordFlat CAFAnaCore
                        ROOT::Core)

target_include_directories( flatten_caf PRIVATE $ENV{SRPROXY_INC} )
