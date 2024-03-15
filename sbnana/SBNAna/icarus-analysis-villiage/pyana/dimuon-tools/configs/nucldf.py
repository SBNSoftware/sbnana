def make_slc_trkdf_nucl(f):
    return make_slc_trkdf(f, mcs=False)
DFS = [make_slc_trkdf_nucl, make_nuclhitdf, make_hdrdf, make_mcnudf]
NAMES = ["trk", "hit", "hdr", "mcnu"]
