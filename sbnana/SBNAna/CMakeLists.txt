# SBNAna starts all its includes with just SBNAna/ -- allow that
#include_directories(..)

# But we also want to find the generated SRProxy.h under "StandardRecord/Proxy"
# include_directories($ENV{MRB_BUILDDIR}/sbnana/sbnana/CAFAna/)

#include_directories($ENV{SRPROXY_INC})

# And expects assert() to still function even in the prof build
add_definitions(-UNDEBUG)

# All variable length arrays
add_definitions(-Wno-vla)

#add_subdirectory(StandardRecord)

 add_subdirectory(Vars)
 add_subdirectory(Cuts)
# #add_subdirectory(Decomp)
# add_subdirectory(Prediction)
# add_subdirectory(Experiment)
# add_subdirectory(Extrap)
# add_subdirectory(Analysis)
# add_subdirectory(Systs)

cet_script(cafe)
# same place as cafe script so we can find it
cet_script(load_cafana_libs.C)
#cet_script(rootlogon.C)

install_headers()
install_source()
