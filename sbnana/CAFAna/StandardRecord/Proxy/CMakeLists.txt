SET(INC_DIR $ENV{SBNANAOBJ_INC})
IF(DEFINED ENV{SBNANAOBJ_INC})
  # We are building with sbnanaobj as an external UPS product
  SET(INC_DIR $ENV{SBNANAOBJ_INC})
ELSE()
  # sbnanaobj is checked out in the mrb area with us
  SET(INC_DIR $ENV{SBNANAOBJ_DIR})
ENDIF()

FILE(GLOB SR_DEPENDENCIES $INC_DIR/sbnanaobj/StandardRecord/*.h)

add_custom_command(# Rebuild if anything in StandardRecord/ changes
                   DEPENDS ${SR_DEPENDENCIES}
                   OUTPUT SRProxy.cxx SRProxy.h FwdDeclare.h
                   COMMAND gen_srproxy -i sbnanaobj/StandardRecord/StandardRecord.h -o SRProxy --target caf::StandardRecord --include-path ${INC_DIR}:$ENV{ROOT_INC} -op sbnana/CAFAna/StandardRecord/Proxy/ --prolog $ENV{SBNANA_DIR}/sbnana/CAFAna/StandardRecord/Proxy/Prolog.h --epilog-fwd $ENV{SBNANA_DIR}/sbnana/CAFAna/StandardRecord/Proxy/EpilogFwd.h  --extra-cflags ' -D_Float16=short'
                   # Couldn't figure out how to make install() install the
                   # headers so do it directly here
                   COMMAND mkdir -p $ENV{MRB_INSTALL}/${sbnana_inc_dir}/sbnana/CAFAna/StandardRecord/Proxy/
                   COMMAND cp SRProxy.h FwdDeclare.h $ENV{MRB_INSTALL}/${sbnana_inc_dir}/sbnana/CAFAna/StandardRecord/Proxy/
  )

include_directories($ENV{SRPROXY_INC})

# This is a very picky error buried inside template instantiations
#add_definitions(-Wno-int-in-bool-context)

cet_make_library(LIBRARY_NAME StandardRecordProxy
                 SOURCE       SRProxy.cxx Instantiations.cxx
                 LIBRARIES    ${ROOT_BASIC_LIB_LIST} ${ROOT_TREEPLAYER})

#install(FILES SRProxy.h FwdDeclare.h
#         DESTINATION ${${product}_inc_dir}/CAFAna/StandardRecord/Proxy/)
