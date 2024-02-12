import pandas as pd
import numpy as np

def broadcast(v, df):
    for vi, ii in zip(v.index.names, df.index.names):
        if vi != ii:
            raise ValueError("Value index (%s) does not match index (%s)." % (str(vi), str(ii)))
    if len(v.index.names) > len(df.index.names):
        raise ValueError("Value index too long.")
    if len(v.index.names) == len(df.index.names):
        return v

    rpt = df.groupby(level=list(range(v.index.nlevels))).size()
    has_value = v.index.intersection(rpt.index)
    v_rpt = np.repeat(v.loc[has_value].values, rpt)

    return pd.Series(v_rpt, df.index).rename(v.name) 

def multicol_add(df, s):
    if isinstance(s.name, str):
        s.name = (s.name,)

    nlevel = max(df.columns.nlevels, len(s.name))
    def pad(c):
       return tuple(list(c) + [""]*(nlevel - len(c))) 

    if df.columns.nlevels < nlevel:
        df.columns = pd.MultiIndex.from_tuples([pad(c) for c in df.columns])
    if len(s.name) < nlevel:
        s.name = pad(s.name)

    return df.join(s)

def multicol_addkey(df,s,fill=np.nan,inplace=False):
    #can't handle s larger then length of df
    if isinstance(s,str):
        s = [s]
    nlevel = df.columns.nlevels
    s_depth = max([len(k.split('.')) for k in s]) #depth of s
    # Pad DataFrame's keys if s has larger depth
    if s_depth > nlevel:
        padding = [''] * (s_depth - nlevel)
        df.columns = pd.MultiIndex.from_tuples([list(col) + padding for col in df.columns])
    col = getcolumns(s,depth=nlevel)
    data = np.full((len(df),len(s)),fill)
    new = pd.MultiIndex.from_tuples(col)
    if inplace:
        df[new] = pd.DataFrame(data, index=df.index, columns=new)
    else:
        return pd.concat([df, pd.DataFrame(data, index=df.index, columns=new)], axis=1)
    
    

def multicol_merge(lhs, rhs, **panda_kwargs):
    # Fix the columns
    lhs_col = lhs.columns
    rhs_col = rhs.columns

    nlevel = max(lhs_col.nlevels, rhs_col.nlevels)

    def pad(c):
       return tuple(list(c) + [""]*(nlevel - len(c))) 

    lhs.columns = pd.MultiIndex.from_tuples([pad(c) for c in lhs_col])
    rhs.columns = pd.MultiIndex.from_tuples([pad(c) for c in rhs_col])

    return lhs.merge(rhs, **panda_kwargs)

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
def unreserve(s):
    #change reserved names
    if s == "index":
        return "idx"
    if s[0].isdigit(): # make the name a legal field 
        return "I" + s
    return s
def pad(b,depth):
    return tuple(b + [""]*(depth - len(b)))
def getcolumns(branches,depth=None):
    # Setup branch names so df reflects structure of CAF file
    bsplit = [b.split(".") for b in branches]
    # Replace any reserved names
    bsplit = [[unreserve(s) for s in b] for b in bsplit]
    if depth is None:
        depth = max([len(b) for b in bsplit])
    return [pad(b,depth) for b in bsplit]
    

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
    #Convert branches to caf format for columns
    columns = getcolumns(branches)
    df.columns = pd.MultiIndex.from_tuples(columns)
    return df
