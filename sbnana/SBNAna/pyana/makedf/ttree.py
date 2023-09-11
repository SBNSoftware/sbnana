def get_tree_name(f,ttree_root):
  """
  return tree name for uproot file
  """
  trees = f.keys()
  matched_trees = [tree for tree in trees if ttree_root in tree]
  return matched_trees[0] #Return first matched instance

