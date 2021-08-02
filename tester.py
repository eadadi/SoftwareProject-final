import sklearn.datasets as dt
import numpy as np
import pandas as pd
from sys import argv
import os
import matplotlib.pyplot as plt

np.random.seed(0)

t0 = 'ea15'
t1 = 'cmd /c "python spkmeans.py'
t2 = '.txt'

os.system('cmd /c "IF NOT EXIST '+t0+' mkdir '+t0+'"')
def merge(x):
    s = x[0]
    for a in range(1,len(x)):
        s += ' '+str(x[a])
    return s

def test():
    for _ in range(1):
        n = str(_)
        data,vector_to_cluster,centers = dt.make_blobs(100,2,return_centers=True, random_state=1)
        np.savetxt(t0+"\\in_"+n+t2,data,fmt="%.4f",delimiter=',')
        cmnd = merge([t1,0,'spk',t0+'\\in_'+n+t2,'>> '+t0+'\\out_'+str(_)+t2+'"'])
        os.system(cmnd)
        df = pd.read_csv('.\\'+t0+'\\'+'in_'+n+t2, skiprows=0, header=None)
        print(df)


#data = dt.make_blobs(n1,n2)[0]
#print(data)

test()

