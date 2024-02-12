import numpy as np

def mag(x, y, z):
    return np.sqrt(x**2 + y**2 + z**2)

def magdf(df):
    return mag(df.x, df.y, df.z)

def dmagdf(df1, df2):
    return mag(df1.x - df2.x, df1.y - df2.y, df1.z - df2.z)

