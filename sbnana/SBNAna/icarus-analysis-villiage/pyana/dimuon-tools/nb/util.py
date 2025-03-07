import pandas as pd
import numpy as np

def InFV(df, inzback, inx=10, iny=10, inzfront=15):
    xmin_C0 = -358.49
    xmax_C0 = -61.94

    xmax_C1 = -xmin_C0
    xmin_C1 = -xmax_C0

    ymin = -181.85999999999999
    ymax = 134.96

    zmin = -894.950652270838
    zmax = 894.950652270838

    xmin_C0 = xmin_C0 + inx
    xmax_C0 = xmax_C0 - inx
    xmin_C1 = xmin_C1 + inx
    xmax_C1 = xmax_C1 - inx
    
    ymin = ymin + iny
    ymax = ymax - iny

    zmin = zmin + inzfront
    zmax = zmax - inzback

    fid = (((df.x < xmax_C0) & (df.x > xmin_C0)) | ((df.x < xmax_C1) & (df.x > xmin_C1))) &\
        (df.y < ymax) & (df.y > ymin) & (df.z < zmax) & (df.z > zmin)

    # further cut bad regions in parts of the detector
    bad_WW = (df.x > 210.29) & (df.x < 358.49) & (df.y > 70) & (df.z > 0) # cable
    bad_EE = (df.x > -358.49) & (df.x < -210.29) & ((df.y > 115) | (df.y < -161.86)) # field cage

    return fid & ~bad_WW & ~bad_EE


def TrkInFV(df):
    return InFV(df, 15.)

def SlcInFV(df):
    return InFV(df, 100.)

def mag(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def magdf(df):
    return mag(df.x, df.y, df.z)

def dmagdf(df1, df2):
    return mag(df1.x - df2.x, df1.y - df2.y, df1.z - df2.z)

def dotdf(df1, df2):
    return df1.x*df2.x + df1.y*df2.y + df1.z*df2.z 

def unitdf(df):
    return df.divide(magdf(df), axis=0)

