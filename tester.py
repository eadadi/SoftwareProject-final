import sklearn.datasets as dt
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
from sys import argv
import os
import matplotlib.pyplot as plt

#np.random.seed(0)

t0 = 'ea15'
t1 = 'cmd /c "python spkmeans.py'
t10 = 'cmd /c ".\\spkmeans.exe'
t2 = '.txt'

os.system('cmd /c "IF NOT EXIST '+t0+' mkdir '+t0+'"')
def merge(x):
    s=str(x[0])
    for a in range(1,len(x)):
        s+=' '+str(x[a])
    return s

p1 = 1000
p2 = 2
if len(argv)>1:
    k = argv[1]
    sk = int(k)
else:
    k = 0
    sk = 1
def test():
    for _ in range(1):
        n = str(_)
        #, random_state=1)
        data,vector_to_cluster,centers = dt.make_blobs(p1,p2,return_centers=True) 
        np.savetxt(t0+"\\in_"+n+t2,data,fmt="%.4f",delimiter=',')
        np.savetxt(t0+"\\true_"+n+t2,centers,fmt="%.4f",delimiter=',')
        cmnd0 = merge([t10,k,'spk',t0+'\\in_'+n+t2,'>> '+t0+'\\out_c_'+n+t2+'"'])
        cmnd1 = merge([t1,k,'spk',t0+'\\in_'+n+t2,'>> '+t0+'\\out_py_'+n+t2+'"'])
        print((cmnd0,cmnd1))
        os.system(cmnd0)
        os.system(cmnd1)
        scikit_kmeans = pd.DataFrame(KMeans(n_clusters =sk).fit(data).cluster_centers_)
        print(scikit_kmeans)
        df1 = pd.read_csv('.\\'+t0+'\\'+'in_'+n+t2, header=None)
        df2 = pd.read_csv('.\\'+t0+'\\'+'out_c_'+n+t2, header=None)
        df3 = pd.read_csv('.\\'+t0+'\\'+'out_py_'+n+t2, skiprows=[0], header=None)
        df4 = pd.read_csv('.\\'+t0+'\\'+'true_'+n+t2, header = None)
        if p2==2:
            plt.scatter(df1[0],df1[1],s=100, c = 'black', alpha = 0.2, label="datapoints")
            plt.scatter(df2[0],df2[1],s=200, color = 'green', marker="+", alpha=0.7 , label="c_centroids")
            plt.scatter(df3[0],df3[1],s=200, c = 'red', marker ="x", alpha=0.7, label="py_centroids")
            plt.scatter(df4[0],df4[1],s=200, c='purple', marker="*", alpha=0.7, label="true_centers")
            plt.scatter(scikit_kmeans[0],scikit_kmeans[1],s=200, c='blue', marker="*", alpha=0.7, label="scikit_kmeans_result")
            plt.grid()
            plt.legend()
            plt.show()

#data = dt.make_blobs(n1,n2)[0]
#print(data)

test()

