from pandas import DataFrame
import numpy as np
from pyanalib import panda_helpers
from sbnd.volume import *
from sbnd.constants import *

class CAF(DataFrame):
    def __init__(self,df):
      super().__init__(df)
      if not self.check_key_structure():
        raise ValueError("key structure not correct")
    def key_length(self):
      return len(self.keys()[0])
    def check_key_structure(self):
      #Keys should each be a tuple of the same size
      length = self.key_length()
      for i,key in enumerate(self.keys()[1:]):
        if len(key) != length: return False
      return True
    def add_key(self,keys,fill=np.nan):
      """
      Add key to dataframe
      """
      updated_df = panda_helpers.multicol_addkey(self, keys,fill=fill)
      # Update the current PFP object with the new DataFrame
      for col in updated_df.columns.difference(self.columns):
          self[col] = updated_df[col]
      return self
    def postprocess(self):
      """
      Run all post processing
      """
      pass
    def clean(self,dummy_vals=[-9999,-999,999,9999],fill=np.nan):
      """
      set all dummy vals to nan
      """
      for i,val in enumerate(dummy_vals):
        self.replace(val,fill,inplace=True)
      return self

