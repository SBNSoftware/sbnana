FILE(GLOB SR_DEPENDENCIES $ENV{SBNANA_DIR}/sbnana/StandardRecord/*.h)

add_custom_command(# Rebuild if anything in StandardRecord/ changes
                   DEPENDS ${SR_DEPENDENCIES}
                   OUTPUT FlatRecord.cxx FlatRecord.h FwdDeclare.h
                   COMMAND gen_srproxy --flat -i sbnanaobj/StandardRecord/StandardRecord.h -o FlatRecord --target caf::StandardRecord --include-path $ENV{SBNANAOBJ_DIR}:$ENV{ROOT_INC} -op sbnana/FlatMaker --prolog $ENV{SBNANA_DIR}/sbnana/FlatMaker/Prolog.h --extra-cflags ' -D_Float16=short'
                   # Couldn't figure out how to make install() install the
                   # headers so do it directly here
                   COMMAND mkdir -p $ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}/include/FlatMaker/
                   COMMAND cp FlatRecord.h FwdDeclare.h $ENV{MRB_INSTALL}/sbnana/$ENV{SBNANA_VERSION}/include/FlatMaker/
  )

include_directories($ENV{SRPROXY_INC})

cet_make_library(LIBRARY_NAME FlatMaker
                 SOURCE       FlatRecord.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST} ${ROOT_TREEPLAYER} # CAFAnaCore
                 )

cet_make_exec(flatten_caf
              SOURCE flatten_caf.cc
              LIBRARIES FlatMaker CAFAnaCore
                        ${ROOT_BASIC_LIB_LIST})
