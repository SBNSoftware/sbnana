
# CAFAna expects assert() to still function even in the prof build
add_definitions(-UNDEBUG)

# All variable length arrays
add_definitions(-Wno-vla)

add_subdirectory(Core)
add_subdirectory(Unfold)
add_subdirectory(Vars)
add_subdirectory(Cuts)
add_subdirectory(Systs)
# #add_subdirectory(Decomp)
add_subdirectory(Extrap)
add_subdirectory(Prediction)
add_subdirectory(Experiment)
add_subdirectory(XSec)
add_subdirectory(Analysis)
add_subdirectory(bin)

# These are in SBNAna for now
#cet_script(cafe)
# same place as cafe script so we can find it
#cet_script(load_cafana_libs.C)
cet_script(cafana.py)
cet_script(rootlogon.C)

install_headers()
install_source()
