import sklearn.datasets as dt
from sklearn.cluster import KMeans
import numpy as np
import pandas as pd
from sys import argv
import os
import matplotlib.pyplot as plt
import time


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
compile_first = 0
if os.name =='nt':
    nova = 0
else:
    nova = 1
rand = 0
sort_results = 0
h = 0
if len(argv)>1:
    if "-h" in argv:
        h = 1
    if "-compile" in argv:
        compile_first = 1
    if "-just_compile" in argv:
        compile_first = 2 
    if "-print" in argv:
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
    if "-5d" in argv:
        p2 = 5
    if "-6d" in argv:
        p2 = 6
    if "-7d" in argv:
        p2 = 7
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
    if "-not_random" in argv:
        rand = 1
        np.random.seed(0)
    if "-sort_results" in argv:
        sort_results = 1
def test():
    if h == 1:
        print("(order of different flags does not matter)")
        print("available options")
        print(" -h shows this message")
        print(" -{x}pts is data length (250,500,750,1000)")
        print(" -{x}C is number of true centers (2-5)")
        print(" -{x}d is feature dimensions (2-7)")
        print(" -print (results to cmd)")
        print(" -sort_results to sort printed results by first index")
        print(" -r to remove files after execution")
        print(" -not_random")
        print(" -compile")
        print(" -just_compile to exit after compilation")
        return

    t0 = 'ea15'
    t1 = "spkmeans.py"
    t10 = "spkmeans"
    if nova == 0:
        t1 = 'cmd /c "python '+t1
        t10 = 'cmd /c "'+t10+".exe"
        os.system('cmd /c "IF NOT EXIST '+t0+' mkdir '+t0+'"')
    else:
        t1 = "python3 "+t1
        os.system('mkdir -p '+t0)
    t2 = '.txt'

    if compile_first > 0:
        name0 = os.path.join(t0,"c_compile_output.txt")
        name1 = os.path.join(t0,"py_compile_output.txt")
        str0 = 'gcc -ansi -Wall -Wextra -Werror -pedantic-errors spkmeans.c -lm -o spkmeans > ' +name0
        str1 = 'setup.py build_ext --inplace > '+name1
        if nova == 0:
            str0 = 'cmd /c "'+str0+'"'
            str1 = 'cmd /c "python ' +str1+'"'
            str00 = 'cmd /c "del '+name0+'"'
            str11 = 'cmd /c "del '+name1+'"'
        else:
            str1 = "python3 "+str1
            str00 = 'rm '+name0
            str11 = 'rm '+name1
        print(str0)
        print(str1)

        os.system(str0)
        os.system(str1) 
        print(open(name0).read())
        print(open(name1).read())
        os.system(str00)
        os.system(str11)
        if compile_first > 1:
            return


    for _ in range(1):
        n = str(_)
        if rand == 1:
            data,vector_to_cluster,centers = dt.make_blobs(p1,p2,centers = C, return_centers=True, random_state=1)
        else:
            data,vector_to_cluster,centers = dt.make_blobs(p1,p2,centers = C, return_centers=True )

        np.savetxt(os.path.join(t0,'in_'+n+t2),data,fmt="%.4f",delimiter=',')
        np.savetxt(os.path.join(t0,'true_'+n+t2),centers,fmt="%.4f",delimiter=',')
        cmnd0 = merge([t10,k,'spk',os.path.join(t0,'in_'+n+t2),'> '+os.path.join(t0,"out_c_"+n+t2)])
        cmnd1 = merge([t1,k,'spk',os.path.join(t0,'in_'+n+t2),'> '+os.path.join(t0,"out_py_"+n+t2)])

        if nova == 0:
            cmnd0 += '"'
            cmnd1 += '"'

        print("Commands to run:")
        print("A." + cmnd0)
        print("B." + cmnd1)

        c_t_s = time.perf_counter()
        os.system(cmnd0)
        c_t_e = time.perf_counter()
        p_t_s = time.perf_counter()
        os.system(cmnd1)
        p_t_e = time.perf_counter()

        c_time = c_t_e-c_t_s
        py_time = p_t_e-p_t_s

        scikit_kmeans = pd.DataFrame(KMeans(n_clusters =sk).fit(data).cluster_centers_)
        if sort_results == 1:
            scikit_kmeans = scikit_kmeans.sort_values(by=0)
        print("scikit_kmeans outputs:")
        print(scikit_kmeans)
        print()

        df1 = pd.read_csv(os.path.join(t0,"in_"+n+t2), header=None)
        df2 = pd.read_csv(os.path.join(t0,"out_c_"+n+t2), header=None)
        df3 = pd.read_csv(os.path.join(t0,"out_py_"+n+t2), skiprows=[0], header=None)
        df4 = pd.read_csv(os.path.join(t0,"true_"+n+t2), header = None)
        if output == 1 or p2>2:
            if sort_results == 1:
                df2 = df2.sort_values(by=0)
                df3 = df3.sort_values(by=0)
            print("c outputs:")
            print(df2)
            print(f'runtime: {c_time:.3f}\n')
            print("py outputs:")
            print(df3)
            print(f'runtime: {py_time:.3f}\n')

        if p2==2 and nova == 0:
            plt.scatter(df1[0],df1[1],s=100, c = 'black', alpha = 0.2, label="datapoints")
            plt.scatter(df2[0],df2[1],s=300, color = 'green', marker="+", alpha=0.7 , label="c_centroids")
            plt.scatter(df3[0],df3[1],s=300, c = 'red', marker ="x", alpha=0.7, label="py_centroids")
            plt.scatter(df4[0],df4[1],s=200, c='purple', marker="*", alpha=0.7, label="true_centers")
            plt.scatter(scikit_kmeans[0],scikit_kmeans[1],s=200, c='blue', marker="*", alpha=0.7, label="scikit_kmeans_result")
            plt.grid()
            plt.legend()
            plt.title(f'py_runtime:{py_time:.3f}\nc_runtime:{c_time:.3f}')
            plt.show()

        if del_flag==1:
            pre0='cmd /c "del '
            end0='"'
            pre1='rm '
            f1 = os.path.join(t0,'in_'+n+t2)
            f2 = os.path.join(t0,'out_c_'+n+t2)
            f3 = os.path.join(t0,'out_py_'+n+t2)
            f4 = os.path.join(t0,'true_'+n+t2)
            for _f in [f1,f2,f3,f4]:
                if nova == 1:
                    os.system(pre1+_f)
                else:
                    os.system(pre0+_f+end0)

test()

