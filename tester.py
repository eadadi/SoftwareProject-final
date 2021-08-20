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

p1 = 250
p2 = 2
del_flag = -1
k = 0
sk = 1
C = 3
output = 0
if len(argv)>1:
    if "-p" in argv:
        output = 1
    if "-1cluster" in argv:
        k = sk = 1
    if "-2clusters" in argv:
        k = sk = 2
    if "-3clusters" in argv:
        k = sk = 3
    if "-4clusters" in argv:
        k = sk = 4
    if "-5clusters" in argv:
        k = sk = 5
    if "-r" in argv:
        del_flag = 1
    if "-2d" in argv:
        p2 = 2
    if "-3d" in argv:
        p2 = 3
    if "-4d" in argv:
        p2 = 4
    if "-1000pts" in argv:
        p1 = 1000
    if "-750pts" in argv:
        p1 = 750
    if "-500pts" in argv:
        p1 = 500
    if "-250pts" in argv:
        p1 = 250
    if "-1C" in argv:
        C = 1
    if "-2C" in argv:
        C = 2
    if "-3C" in argv:
        C = 3
    if "-4C" in argv:
        C = 4
    if "-5C" in argv:
        C = 5
def test():
    for _ in range(1):
        n = str(_)
        #, random_state=1)
        data,vector_to_cluster,centers = dt.make_blobs(p1,p2,centers = C, return_centers=True )
        np.savetxt(t0+"/in_"+n+t2,data,fmt="%.4f",delimiter=',')
        np.savetxt(t0+"/true_"+n+t2,centers,fmt="%.4f",delimiter=',')
        cmnd0 = merge([t10,k,'spk',t0+'/in_'+n+t2,'>> '+t0+'/out_c_'+n+t2+'"'])
        cmnd1 = merge([t1,k,'spk',t0+'/in_'+n+t2,'>> '+t0+'/out_py_'+n+t2+'"'])
        print((cmnd0,cmnd1))
        os.system(cmnd0)
        os.system(cmnd1)
        scikit_kmeans = pd.DataFrame(KMeans(n_clusters =sk).fit(data).cluster_centers_)
        print("scikit_kmeans outputs:")
        print(scikit_kmeans)
        print()
        df1 = pd.read_csv('./'+t0+'/'+'in_'+n+t2, header=None)
        #print()
        df2 = pd.read_csv('./'+t0+'/'+'out_c_'+n+t2, header=None)
        #print()
        df3 = pd.read_csv('./'+t0+'/'+'out_py_'+n+t2, skiprows=[0], header=None)
        #print()
        df4 = pd.read_csv('./'+t0+'/'+'true_'+n+t2, header = None)
        #print()
        if output == 1 or p2>2:
            df2 = df2.sort_values(by=0)
            df3 = df3.sort_values(by=0)
            print("c outputs:")
            print(df2)
            print("py outputs:")
            print(df3)

        if p2==2:
            plt.scatter(df1[0],df1[1],s=100, c = 'black', alpha = 0.2, label="datapoints")
            plt.scatter(df2[0],df2[1],s=300, color = 'green', marker="+", alpha=0.7 , label="c_centroids")
            plt.scatter(df3[0],df3[1],s=300, c = 'red', marker ="x", alpha=0.7, label="py_centroids")
            plt.scatter(df4[0],df4[1],s=200, c='purple', marker="*", alpha=0.7, label="true_centers")
            plt.scatter(scikit_kmeans[0],scikit_kmeans[1],s=200, c='blue', marker="*", alpha=0.7, label="scikit_kmeans_result")
            plt.grid()
            plt.legend()
            plt.show()

        if del_flag==1:
            os.system('cmd /c "del ./'+t0+'/in_'+n+t2+'"')
            os.system('cmd /c "del ./'+t0+'/out_c_'+n+t2+'"')
            os.system('cmd /c "del ./'+t0+'/out_py_'+n+t2+'"')
            os.system('cmd /c "del ./'+t0+'/true_'+n+t2+'"')
#data = dt.make_blobs(n1,n2)[0]
#print(data)

test()

