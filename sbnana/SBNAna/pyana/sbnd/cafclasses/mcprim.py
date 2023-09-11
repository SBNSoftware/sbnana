from pandas import DataFrame
import numpy as np
from pyanalib import panda_helpers
from sbnd.volume import *
from sbnd.constants import *
from .parent import CAF

class MCPRIM(CAF):
    def __init__(self,df):
      super().__init__(df)
    def postprocess(self):
      """
      Run all post processing
      """
      self.clean()
      self = self.drop_noninteracting()
      #self = self.drop_neutrinos() - this is taken care of by drop noninteracting
      self.add_fv()
      return MCPRIM(self) #Have to return since we're dropping rows
    def add_fv(self):
      """
      Add containment for start and end of particle
      """
      keys = [
        'in_tpc',
      ]
      self.add_key(keys)
      cols = panda_helpers.getcolumns(keys,depth=self.key_length())
      self.loc[:,cols[0]] = involume(self.start) | involume(self.end)
      return None
    def get_true_parts(self,remove_nan=True,**dropna_args):
      """
      return true particles from track and shower matching
      """
      particles = self.copy()
      if remove_nan:
        particles = particles.dropna(**dropna_args)
      return particles
    
    def get_true_parts_from_pdg(self,pdg,remove_nan=True,**dropna_args):
      """
      Return particles from pdg
      """
      particles = self.get_true_parts(remove_nan=remove_nan,**dropna_args)
      particles[particles.pdg == pdg]
      
      return particles
    def drop_neutrinos(self):
      """
      Drop rows with neutrinos
      """
      has_nu = (abs(self.pdg) == 14) | (abs(self.pdg) == 12)
      self = MCPRIM(self[~has_nu])
      return self
    def drop_noninteracting(self):
      """
      Drop rows where the visible energy is zero
      """
      visible = (((self.plane.I0.I0.nhit + self.plane.I0.I1.nhit + self.plane.I0.I2.nhit) > 0) &
                 [val is not np.nan for val in self.plane.I0.I0.nhit])
      self = MCPRIM(self[visible])
      return self
      
      