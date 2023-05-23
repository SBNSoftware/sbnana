import numpy as np
import pandas as pd
from sbnd.general import utils
from sbnd.constants import *
from sbnd.prism import *
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, mean_squared_log_error
from scipy.stats import pearsonr
import copy
#from pyanalib import panda_helpers

def get_neutrino_dir(start):
  """
  start is start location is either shower, pfp, or other 
  """
  
  #Get the direction
  neutrino_dir = np.array([start.x - prism_centroid[0],
                          start.y - prism_centroid[1],
                          start.z + distance_from_bnb])
  neutrino_dir /= np.linalg.norm(neutrino_dir,axis=0)
  
  return neutrino_dir.T
def get_theta(pfpdir,nudir):
  """
  Return angle wrt beam modified by reconstucted vertex
  """
  if (pfpdir.shape != nudir.shape):
    print('neutrino and pfp direction are not the same shape')
  
  #Normalize
  pfpdir /= np.linalg.norm(pfpdir,axis=1).reshape((len(pfpdir),1))
  nudir /= np.linalg.norm(nudir,axis=1).reshape((len(nudir),1))
  
  #Calc theta
  theta = np.zeros(len(nudir))
  for i,(pfdir,ndir) in enumerate(zip(pfpdir,nudir)):
    if all(np.isnan(pfdir)) or all(np.isnan(ndir)): 
      theta[i] = np.nan
    else:
      theta[i] = np.arccos(np.dot(pfdir,ndir))
  return theta
def get_err(reco,true,normalize=True):
  """
  Get error between two series
  """
  if len(reco.shape) != len(true.shape):
    print('reco and true are different shapes')
  if len(reco.shape) == 2:
    err = np.linalg.norm(reco-true,axis=1)
    if normalize: err /= np.linalg.norm(true,axis=1)
  elif len(reco.shape) == 1:
    err = reco-true
    if normalize: err /= true
  return err

def get_row_vals(series,indices=None,mode='sum'):
  """
  returns a series of length indices where the series is the series to act on
  """
  if indices is None:
    indices = series.index.drop_duplicates()
  series_subset = series.loc[indices]
  if mode == 'sum':
    return series_subset.groupby(series_subset.index).sum()
  elif mode == 'min':
    return series_subset.groupby(series_subset.index).min()
  elif mode == 'mean':
    return series_subset.groupby(series_subset.index).mean()
  else: #assume mode is a key now - doesn't work
    max_idx = series_subset.groupby(series_subset.index)[mode].idxmax()
    return series_subset.loc[max_idx,key]
    

def get_r2(reco,true):
  """
  Get r2 between true and reco
  """
  if len(reco.shape) != len(true.shape):
    print('reco and true are different shapes')
  #2d not supported yet
  if len(reco.shape) == 2:
    return np.full(reco.shape,np.nan)
  elif len(reco.shape) == 1:
    corr_matrix = np.corrcoef(true,reco)
    corr = corr_matrix[0,1]
    R_sq = corr**2
  pass

def count_true_tracks_showers(pdgs):
  """
  Return true tracks and showers as two integers
  """
  nshw = 0
  ntrk = 0
  #if isinstance(pdgs,np.int32) or isinstance(pdgs,np.int64):
  #  return ntrk,nshw
  counts_dict = pdgs.value_counts().to_dict()
  for pdg,count in counts_dict.items():
    if abs(pdg) in [11,22]: nshw += count 
    if pdg in [111]: nshw += 2*count 
    if abs(pdg) in [13,2212,211]: ntrk += count
  return ntrk,nshw

def get_nonna(x,y):
  """
  REturn non na in common for x and y series
  """
  nainds = x.isna() | y.isna()
  return x[~nainds],y[~nainds]
  
def merge_objs(objs,keys,clas):
  """
  Merge two dataframe like objects and convert to clas
  """  
  objs_copy = [copy.deepcopy(obj) for obj in objs]
  #modify keys
  for i,obj in enumerate(objs_copy):
    cols = pd.MultiIndex.from_tuples([(keys[i],) + key for key in obj.columns])
    obj.columns = cols
    objs_copy[i] = obj #update
  return clas(pd.concat(objs_copy,axis=0))
  
  
  
  

  
  
