import numpy as np
from pandas import DataFrame
import pandas as pd
from pyanalib import panda_helpers
from sbnd.cafclasses.object_calc import *
from sbnd.volume import *
from sbnd.constants import *
from .parent import CAF

ntrk_key = ('ntrk','','','','')
nshw_key = ('nshw','','','','')
trkscore_key = ('trackScore','','','','')
pdg_key = ('shw','truth','p','pdg','')

class PFP(CAF):
  def __init__(self,df):
    super().__init__(df)
  def postprocess(self,fill=-5,dummy=-5,is_cheat=False):
    """
    Do all the post processing in the correct order
    """
    self.shw_energy_fix(fill=fill,dummy=dummy)
    #fill dummy's with nan before calculations
    self.clean() 
    self.count_tracks_showers(is_cheat=is_cheat)
    self.add_reco_containment()
    self.add_trk_bestenergy()
    self.add_neutrino_dir()
    self.add_theta()
    self.add_Etheta()
    self.clean(dummy_vals=[-5])
    #self.add_errs()
    
  def stats(self):
    self.add_stat(get_err,'mae',normalize=True)
    self.add_stat(get_err,'dif',normalize=False)
    
    
  def get_particles(self,pdgs,remove_nan=True,**dropna_args):
    """
    Return list of particles from list of pdgs
    """
    parts = [None]*len(pdgs)
    for i,pdg in enumerate(pdgs):
      parts[i] = self.get_true_parts_from_pdg(pdg,remove_nan=remove_nan,**dropna_args)
    return parts

  def shw_energy_fix(self,fill=-5,dummy=-5):
      ## correct the reconstructed shower energy column
      nhits2 = ((self.shw.plane.I2.nHits >= self.shw.plane.I1.nHits) & (self.shw.plane.I2.nHits>= self.shw.plane.I0.nHits))
      nhits1 = ((self.shw.plane.I1.nHits >= self.shw.plane.I2.nHits) & (self.shw.plane.I1.nHits>= self.shw.plane.I0.nHits))
      nhits0 = ((self.shw.plane.I0.nHits >= self.shw.plane.I2.nHits) & (self.shw.plane.I0.nHits>= self.shw.plane.I1.nHits))
      # if energy[plane] is positive
      energy2 = (self.shw.plane.I2.energy > 0 )
      energy1 = (self.shw.plane.I1.energy > 0 )
      energy0 = (self.shw.plane.I0.energy > 0 )
      conditions = [(nhits2 & energy2),
                  (nhits1 & energy1),
                  (nhits0 & energy0),
                  (((nhits2 & energy2)== False) & (energy1) & (self.shw.plane.I1.nHits>= self.shw.plane.I0.nHits)), # if 2 is invalid, and 1 is positive and 1>0, go with 1 
                  (((nhits2 & energy2)== False) & (energy0) & (self.shw.plane.I0.nHits>= self.shw.plane.I1.nHits)), # if 2 is invalid, and 0 is positive and 0>1, go with 0
                  (((nhits1 & energy1)== False) & (energy2) & (self.shw.plane.I2.nHits>= self.shw.plane.I0.nHits)), # if 1 is invalid, and 2 is positive and 2>0, go with 2 
                  (((nhits1 & energy1)== False) & (energy0) & (self.shw.plane.I0.nHits>= self.shw.plane.I2.nHits)), # if 1 is invalid, and 0 is positive and 0>2, go with 0
                  (((nhits0 & energy0)== False) & (energy2) & (self.shw.plane.I2.nHits>= self.shw.plane.I1.nHits)), # if 0 is invalid, and 2 is positive and 2>1, go with 2              
                  (((nhits0 & energy0)== False) & (energy1) & (self.shw.plane.I1.nHits>= self.shw.plane.I2.nHits)), # if 0 is invalid, and 1 is positive and 1>2, go with 1 
                  ((self.shw.plane.I2.nHits==dummy) & (self.shw.plane.I1.nHits==dummy) & (self.shw.plane.I0.nHits==dummy))]
      shw_choices = [ self.shw.plane.I2.energy,
                      self.shw.plane.I1.energy,
                      self.shw.plane.I0.energy,
                      self.shw.plane.I1.energy,
                      self.shw.plane.I0.energy,
                      self.shw.plane.I2.energy,
                      self.shw.plane.I0.energy,
                      self.shw.plane.I2.energy,
                      self.shw.plane.I1.energy,
                      -1]
      self.loc[:,panda_helpers.getcolumns(['shw.bestplane_energy'],depth=self.key_length())] = np.select(conditions, shw_choices, default = fill)
      return None
  #----This can and should be optimized in the future
  def count_tracks_showers(self,is_cheat=False):
    """
    Store number of tracks and showers in pfp df
    """
    keys = [
      'ntrk',
      'nshw',
      'true_ntrk',
      'true_nshw',
    ]
    cols = panda_helpers.getcolumns(keys,depth=self.key_length())
    if is_cheat:
      self.add_key('ntrk',fill=0)
      self.add_key('nshw',fill=1)
      self.add_key('true_ntrk',fill=0)
      self.add_key('true_nshw',fill=1)
    if not is_cheat:
      self.add_key(keys,fill=0)
      particles = self.get_true_parts(remove_nan=False)
      for ix,row in self.iterrows():
        self.loc[ix,cols[0]] = ((self.loc[ix,trkscore_key] > 0.5).values & (self.loc[ix,trkscore_key] <= 1).values).sum()
        self.loc[ix,cols[1]] = ((self.loc[ix,trkscore_key] <= 0.5).values & (self.loc[ix,trkscore_key] >=0).values).sum()
      for ind in particles.index.drop_duplicates():
        true_ntrk,true_nshw = count_true_tracks_showers(particles.pdg.loc[ind])
        self.loc[ind,cols[2]] = true_ntrk
        self.loc[ind,cols[3]] = true_nshw
    return None
  def get_total_reco_energy(self):
    """
    Return total reco energy from pfp
    """
    return get_row_vals(self.shw.bestplane_energy) #Add tracks later
  def get_min_theta(self):
    """
    Return minimum theta from pfp
    """
    return get_row_vals(self.shw.theta,mode='min')
  def remove_dummy_track_score(self):
    """
    Make a column that determines if the pfp is valid or not. 
    For now this just means it has a track score E [0,1]
    """
    mask_lower = self.loc[:,trkscore_key] >= 0
    mask_upper = self.loc[:,trkscore_key] <= 1
    mask = np.logical_and(mask_lower,mask_upper)
    self = self[mask]
    return None
  
  def get_true_parts(self,remove_nan=True,**dropna_args):
    """
    return true particles from track and shower matching
    """
    particles = self.copy()
    particles = pd.concat([particles.trk.truth.p,particles.shw.truth.p],axis=0).drop_duplicates()
    if remove_nan:
      particles = particles.dropna(**dropna_args)
    return particles
  
  def get_true_parts_from_pdg(self,pdg,remove_nan=True,**dropna_args):
    """
    Return particles from pdg
    """
    particles = self.get_true_parts(remove_nan=remove_nan,**dropna_args)
    
    return particles[particles.pdg == pdg]
  def add_true_tracks_showers(self):
    """
    count truth level tracks showers
    """
    pass
  
  def add_reco_containment(self):
    """
    Add containment 1 or 0 for each pfp
    """
    keys = [
      'shw.cont_tpc',
      'trk.cont_tpc',
    ]
    self.add_key(keys)
    cols = panda_helpers.getcolumns(keys,depth=self.key_length())
    self.loc[:,cols[0]] = involume(self.shw.start) & involume(self.shw.end)
    self.loc[:,cols[1]] = involume(self.trk.start) & involume(self.trk.end)
    return None
  
  def add_trk_bestenergy(self):
    """
    Get best energy for track
    """
    keys = [
      'trk.bestenergy',
    ]
    self.add_key(keys)
    cols = panda_helpers.getcolumns(keys,depth=self.key_length())
    cont = self.trk.cont_tpc #containment boolean
    #----Assume muon, for now
    self.loc[:,cols[0]][cont] = np.sqrt(self.trk.rangeP.p_muon**2+mu**2)
    self.loc[:,cols[0]][~cont] = np.sqrt(self.trk.mcsP.fwdP_muon**2+mu**2)
    return None
    
  def add_neutrino_dir(self):
    """
    add neutrino direction to df
    """
    keys = [
      'shw.nudir.x','shw.nudir.y','shw.nudir.z',
      'shw.truth.p.nudir.x','shw.truth.p.nudir.y','shw.truth.p.nudir.z',
      'trk.nudir.x','trk.nudir.y','trk.nudir.z',
      'trk.truth.p.nudir.x','trk.truth.p.nudir.y','trk.truth.p.nudir.z',
    ]
    self.add_key(keys)
    cols = panda_helpers.getcolumns(keys,depth=self.key_length())
    self.loc[:,cols[0:3]] = get_neutrino_dir(self.shw.start)
    self.loc[:,cols[3:6]]= get_neutrino_dir(self.shw.truth.p.start)
    self.loc[:,cols[6:9]] = get_neutrino_dir(self.trk.start)
    self.loc[:,cols[9:12]] = get_neutrino_dir(self.trk.truth.p.start)
    return None
  
  def add_theta(self):
    """
    add theta wrt nu direction
    """
    keys = [
      'shw.theta',
      'shw.truth.p.theta',
      'trk.theta',
      'trk.truth.p.theta',
    ]
    self.add_key(keys)
    cols = panda_helpers.getcolumns(keys,depth=self.key_length())
    self.loc[:,cols[0]] = get_theta(self.shw.dir.values,self.shw.nudir.values)
    self.loc[:,cols[1]]= get_theta(self.shw.truth.p.genp.values,self.shw.truth.p.nudir.values)
    self.loc[:,cols[2]] = get_theta(self.trk.dir.values,self.trk.nudir.values)
    self.loc[:,cols[3]] = get_theta(self.trk.truth.p.genp.values,self.trk.truth.p.nudir.values)
    return None
  def add_Etheta(self):
    """
    add theta wrt nu direction
    """
    keys = [
      'shw.Etheta2',
      'shw.truth.p.Etheta2',
      'trk.Etheta2',
      'trk.truth.p.Etheta2',
    ]
    self.add_key(keys)
    cols = panda_helpers.getcolumns(keys,depth=self.key_length())
    self.loc[:,cols[0]] = self.shw.theta**2*self.shw.bestplane_energy
    self.loc[:,cols[1]]= self.shw.truth.p.theta**2*self.shw.truth.p.genE
    self.loc[:,cols[2]] = self.trk.theta**2*self.trk.bestenergy
    self.loc[:,cols[3]]= self.trk.truth.p.theta**2*self.trk.truth.p.genE
    return None
  def add_visenergy(self):
    """
    Best plane true visible energy
    """
    keys = [
      'shw.truth.p.bestplane_energy',
      'trk.truth.p.bestplane_energy',
    ]
    self.add_key(keys)
    cols = panda_helpers.getcolumns(keys,depth=self.key_length())
    self.loc[:,cols[0]] = self.shw.truth
    
  def add_stat(self,func,key,**funkwargs):
    """
    add relevant errors to df
    """
    keys = [
      # Key names for comparing true and reco
      f'shw.stat.{key}.Etheta2',
      f'shw.stat.{key}.energy',
      f'shw.stat.{key}.theta',
      f'shw.stat.{key}.start.x',
      f'shw.stat.{key}.start.y',
      f'shw.stat.{key}.start.z',
      f'shw.stat.{key}.vtx',
      f'trk.stat.{key}.Etheta2',
      f'trk.stat.{key}.energy',
      f'trk.stat.{key}.theta',
      f'trk.stat.{key}.start.x',
      f'trk.stat.{key}.start.y',
      f'trk.stat.{key}.start.z',
      f'trk.stat.{key}.vtx',
    ]
    
    self.add_key(keys)
    cols = panda_helpers.getcolumns(keys,depth=self.key_length())
    
    #showers
    self.loc[:,cols[0]] = func(self.shw.Etheta2,self.shw.truth.p.Etheta2,**funkwargs)
    self.loc[:,cols[1]] = func(self.shw.bestplane_energy,self.shw.truth.p.genE,**funkwargs)
    self.loc[:,cols[2]] = func(self.shw.theta,self.shw.truth.p.theta,**funkwargs)
    self.loc[:,cols[3]] = func(self.shw.start.x,self.shw.truth.p.start.x,**funkwargs)
    self.loc[:,cols[4]] = func(self.shw.start.y,self.shw.truth.p.start.y,**funkwargs)
    self.loc[:,cols[5]] = func(self.shw.start.z,self.shw.truth.p.start.z,**funkwargs)
    self.loc[:,cols[6]] = func(self.shw.start.values,self.shw.truth.p.start.values,**funkwargs)
    #tracks
    self.loc[:,cols[7]] = func(self.trk.Etheta2,self.trk.truth.p.Etheta2,**funkwargs)
    self.loc[:,cols[8]] = func(self.trk.bestenergy,self.trk.truth.p.genE,**funkwargs)
    self.loc[:,cols[9]] = func(self.trk.theta,self.trk.truth.p.theta,**funkwargs)
    self.loc[:,cols[10]] = func(self.trk.start.x,self.trk.truth.p.start.x,**funkwargs)
    self.loc[:,cols[11]] = func(self.trk.start.y,self.trk.truth.p.start.y,**funkwargs)
    self.loc[:,cols[12]] = func(self.trk.start.z,self.trk.truth.p.start.z,**funkwargs)
    self.loc[:,cols[13]] = func(self.trk.start.values,self.trk.truth.p.start.values,**funkwargs)

    #If any stat is na, the others are also na
    self[self.trk.stat.isna().any(axis=1)].trk.stat = np.nan
    self[self.shw.stat.isna().any(axis=1)].shw.stat = np.nan

    return None
