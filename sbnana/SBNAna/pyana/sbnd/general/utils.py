"""
General helper functions
"""

import numpy as np
import shutil
from pathlib import Path
import os

def move_file(fname,folder_name,overwrite=True):
  # Move the file to the folder
  src_file = f'{fname}'
  dst_folder = Path(folder_name)
  
  # Make sure the folder exists
  os.makedirs(folder_name, exist_ok=True)

  # Check if the file exists in the source path
  if os.path.isfile(src_file):

    # Generate the destination path by appending the filename to the destination directory
    dst_file = os.path.join(folder_name, os.path.basename(src_file))

    if overwrite:
      shutil.move(src_file, dst_file)
    else:
      # If the file already exists in the destination directory, make a copy with a numerical suffix
      if os.path.exists(dst_file):
        suffix = 1
        while True:
          new_dst_file = os.path.join(folder_name, os.path.splitext(os.path.basename(src_file))[0] + "_" + str(suffix) + os.path.splitext(os.path.basename(src_file))[1])
          if not os.path.exists(new_dst_file):
            shutil.copy2(src_file, new_dst_file)
            dst_file = new_dst_file
            break
          suffix += 1
      shutil.move(src_file, dst_file)

def flatten_list(l):
  """
  Flatten a list of lists into a list.
  """
  ret_l = []
  for i,sublist in enumerate(l):
    ret_l.extend(sublist)
    
  return ret_l
  #return [item for sublist in l for item in sublist]

def sqrt_sum_of_squares(numbers):
  """
  Takes a list of numbers and returns the sum of their squares.
  """
  return np.sqrt(np.sum([x**2 for x in numbers]))

def find_indices_in_common(list1,list2):
  """
  Return indices in common between list1 and list2
  """
  set1 = set(list1)
  set2 = set(list2)
  return list(set1.intersection(set2))

def find_indices_not_in_common(list1,list2):
  """
  Find list of indices not in common between two lists.
  Return the values not in list 1
  """
  set1 = set(list1)
  set2 = set(list2)
  return list(set2.difference(set1))

def common_indices(*lists):
  """
  Return a list of common values from a list of lists.
  """
  common = []
  for elem in lists[0]:
    if all(elem in lst for lst in lists[1:]):
      common.append(elem)
  return common

