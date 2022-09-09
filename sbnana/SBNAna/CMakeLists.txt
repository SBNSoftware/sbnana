include_directories($ENV{SRPROXY_INC})

include_directories($ENV{EIGEN_INC})

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
#cet_script(rootlogon.C)

install_headers()
install_source()
