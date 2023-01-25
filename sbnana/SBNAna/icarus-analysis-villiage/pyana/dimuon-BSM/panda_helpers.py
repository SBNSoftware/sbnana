import pandas as pd
import numpy as np

def detect_vectors(tree, branch):
    ret = []
    hierarchy = branch.split(".")
    for i in range(len(hierarchy)):
        subbranch = ".".join(hierarchy[:i+1])
        lenbranch = subbranch + "..length"
        if lenbranch in tree.keys():
            ret.append(subbranch)
    return ret

def idarray(ids, lens):
    return np.repeat(ids.values, lens.values)

def loadbranches(tree, branches, **uprargs):
    vectors = []
    for i,branch in enumerate(branches):
        this_vectors = detect_vectors(tree, branch)
        if i == 0:
            vectors = this_vectors
        elif len(this_vectors) == 0: # This case is ok since it will automatically broadcast
            pass
        # All the branches must have the same vector structure for this to work
        elif vectors != this_vectors:
            raise ValueError("Branches %s and %s have different vector structures in the CAF." % (branches[0], branch))

    lengths = [tree.arrays([v+"..length"], library="pd", **uprargs) for v in vectors]
    data = tree.arrays(branches, library="pd", **uprargs)

    # If there's no vectors, we can just return the top guy
    if len(lengths) == 0:
        data.index.name = "entry"
        df = data
    else:
        tomerge = lengths + [data]
        # Otherwise, iteratively merge the branches
        df = tomerge[0]
        df.index.name = "entry"

        # handle the rest
        for i in range(1, len(tomerge)):
            thismerge = tomerge[i]
            v_ind = i - 1

            # Build the information in the right-hand table needed to do the join
            # The "upidx" will be matched to the index vector-by-vector
            for i in range(v_ind):
                thismerge[vectors[v_ind] + "..upidx" + str(i)] = idarray(df[vectors[i]+ "..index"], df[vectors[v_ind] + "..length"])

            # Inner join! Throw away rows in the right-hand with no match in the left-hand
            df = pd.merge(df, thismerge, how="inner",
                         left_on = ["entry"] + [v+"..index" for v in vectors[:v_ind]],
                         right_on = ["entry"] + [vectors[v_ind] + "..upidx" + str(i) for i in range(v_ind)],
                         validate="one_to_many")

            # Make sure no rows in the right-hand were dropped
            assert(df.shape[0] == thismerge.shape[0])

            # postprocess: build the index
            df[vectors[v_ind] + "..index"] = df.groupby(["entry"] + [v+"..index" for v in vectors[:v_ind]]).cumcount()

        # Set the index
        df.set_index([v+"..index" for v in vectors], append=True, verify_integrity=True, inplace=True)

        # Drop all the metadata info we don't need anymore
        df = df[branches]

    # Setup branch names so df reflects structure of CAF file
    bsplit = [b.split(".") for b in branches]
    # Replace any reserved names
    def unreserve(s):
        if s == "index":
            return "idx"
        if s[0].isdigit(): # make the name a legal field 
            return "I" + s
        return s

    bsplit = [[unreserve(s) for s in b] for b in bsplit]

    depth = max([len(b) for b in bsplit])

    def pad(b):
        return tuple(b + [""]*(depth - len(b)))

    df.columns = pd.MultiIndex.from_tuples([pad(b) for b in bsplit])

    return df


