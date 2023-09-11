import numpy as np
import pandas as pd

av = [
  [-200,200],#x
  [-200,200],#y
  [0,500],#z
  ]

def involume(coords,volume=av):
  xb,yb,zb = volume
  if isinstance(coords, pd.DataFrame):
    return ((xb[0] < coords.x) & (coords.x < xb[1]) &
        (yb[0] < coords.y) & (coords.y < yb[1]) &
        (zb[0] < coords.z) & (coords.z < zb[1]))
  else:
    print('Make a dataframe or GET OUT OF MY FACE')
    raise Exception('deal with it')
