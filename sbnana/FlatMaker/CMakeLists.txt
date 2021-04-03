SET(INC_DIR $ENV{SBNANAOBJ_INC})
IF(DEFINED ENV{SBNANAOBJ_INC})
  # We are building with sbnanaobj as an external UPS product
  SET(INC_DIR $ENV{SBNANAOBJ_INC})
ELSE()
  # sbnanaobj is checked out in the mrb area with us
  SET(INC_DIR $ENV{SBNANAOBJ_DIR})
ENDIF()

FILE(GLOB SR_DEPENDENCIES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${INC_DIR}/sbnanaobj/StandardRecord/*.h)

add_custom_command(# Rebuild if anything in StandardRecord/ changes
                   DEPENDS ${SR_DEPENDENCIES}
                   OUTPUT FlatRecord.cxx FlatRecord.h FwdDeclare.h
                   COMMAND gen_srproxy --flat -i sbnanaobj/StandardRecord/StandardRecord.h -o FlatRecord --target caf::StandardRecord --include-path ${INC_DIR}:$ENV{ROOT_INC} -op sbnana/FlatMaker --prolog $ENV{SBNANA_DIR}/sbnana/FlatMaker/Prolog.h --extra-cflags ' -D_Float16=short'
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

install_headers(EXTRAS $ENV{MRB_BUILDDIR}/sbnana/sbnana/FlatMaker/FlatRecord.h $ENV{MRB_BUILDDIR}/sbnana/sbnana/FlatMaker/FwdDeclare.h)
