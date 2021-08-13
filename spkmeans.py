#imports:
import sys
import enum
import numpy as np
import pandas as pd
import spkmeans as sp
import matplotlib.pyplot as plt

#"goal" enum:
class goalEnum (enum.Enum):
    wam = "wam"
    ddg = "ddg"
    lnorm = "lnorm"
    jacobi = "jacobi"
    spk = "spk"

#methods:
def getInput():
    args = sys.argv
    assert(len(args)==4)
    k,goal,file_name = (args[1],args[2],args[3])
    k = int(k)
    goal = goalEnum(goal)
    file_name = file_name
    return(k,goal,file_name)

#kmeanspp logic:
def build_probabilities (D, p):
    s = sum(D)
    for i in range(len(D)):
        if s!=0:
            p[i] = D[i] /s
    return p

def squared_distance (u, v):
    res = 0
    for i in range(len(u)):
        res += (u[i]-v[i])**2
    return res

def updateDi (datapoints, i, Z, centroids, D):
    tmp = squared_distance(datapoints[i], centroids[Z-1])
    D[i] = min(D[i], tmp)

def first_selection(datapoints, selected_indexes):
    index = np.random.choice(len(datapoints))
    u1 = np.copy(datapoints[index])
    selected_indexes[0] = index
    return u1

def kmeanspp(datapoints, number_of_clusters):
    k = number_of_clusters
    #print("k is "+str(k))
    initial_centroids = [0]*k
    selected_indexes = [0]*k
    initial_centroids[0] = first_selection(datapoints, selected_indexes)
    Z = 1
    D = [squared_distance(datapoints[i],initial_centroids[0]) for i in range(len(datapoints))]
    probabilities = [0] * len(D)
    while Z < k:
        #print("iteration "+str(Z)+"...")
        for i in range(len(datapoints)):
            updateDi (datapoints, i, Z, initial_centroids, D)
        build_probabilities(D, probabilities)
        j = np.random.choice(range(0,len(datapoints)), p = probabilities)
        selected_indexes[Z] = j
        initial_centroids[Z] = np.copy(datapoints[j])
        Z += 1
    return (selected_indexes, [array.tolist() for array in initial_centroids])

#calculate actual centers based on the map: vector_i belongs to cluster_j
def calcCentroidsBasedOnMap (data, data_to_clusters_map, clustersNumber):
    sums = [[0]*len(data[0]) for i in range(clustersNumber)] 
    amounts = [0]*clustersNumber
    for i in range(len(data)):
        for j in range(len(data[0])):
            sums[data_to_clusters_map[i]][j] += data[i][j]
        amounts[data_to_clusters_map[i]] += 1
    for i in range(clustersNumber):
        for j in range(len(data[0])):
            if amounts[i]!=0:
                sums[i][j] /= amounts[i]
    return sums

#main method:
def main():
    np.random.seed(0)
    try:
        k, goal, file_name = getInput()
        pd_data = pd.read_csv(file_name, sep=",", header=None)
        data = pd.read_csv(file_name, sep=",", header=None).to_numpy().tolist()
        n = len (data)
        d = len(data[0])
        if goal == goalEnum.wam:
            result = sp.wam(data, n, d, k)
        elif goal == goalEnum.ddg:
            result = sp.ddg(data, n, d, k)
        elif goal == goalEnum.lnorm:
            result = sp.lnorm(data, n, d, k)
        elif goal == goalEnum.jacobi:
            result = sp.jacobi(data, n, d, k)
        elif goal == goalEnum.spk:
            Rnk,Rnk_to_Data = sp.spk(data, n, d, k)
            numpyData = np.array(Rnk)
            clustersNumber = len(Rnk[0])
            #print("strting choose initials...")
            initial_indexes,initial_vectors = kmeanspp(numpyData,clustersNumber)
            #print("going back to c....")
            vectors_to_clusters_map = sp.kmeans(Rnk , initial_vectors, 1000)
            result = calcCentroidsBasedOnMap(data, vectors_to_clusters_map, clustersNumber) 
            result = {'initial_centroids_indexes':initial_indexes, 'final_centroids':result}
    except ValueError:
        print("Invalid Input!")
        return
    #except:
    #    print("An Error Has Occured")
    #    return
    else:
        if goal ==goalEnum.spk:
            a = [str(y) for y in result['initial_centroids_indexes']]
            b = result['final_centroids']
            print(','.join(a))
            for x in b:
                x0 = ['%.4f'%y for y in x]
                print(','.join(x0))
        else:
            for x in result:
                x0 = ['%.4f'%y for y in x]
                print(','.join(x0))



#RUN:
main()
