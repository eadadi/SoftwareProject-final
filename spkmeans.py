#imports:
import sys
import enum
import numpy as np
import pandas as pd
import spkmeans as sp

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


#main method:
def main():
    np.random.seed(0)
    try:
        k, goal, file_name = getInput()
        data = pd.read_csv(file_name, sep=",", header=None).to_numpy().tolist()
        n = len (data)
        d = len(data[0]) 

        if goal == goalEnum.wam:
            result = sp.wam(data, n, d)
        if goal == goalEnum.ddg:
            result = sp.ddg(data, n, d)
        if goal == goalEnum.lnorm:
            result = sp.lnorm(data, n, d)
        if goal == goalEnum.jacobi:
            result = sp.jacobi(data, n, d)
        if goal == goalEnum.spk:
            result = sp.spk(data, n, d)

    except ValueError:
        print("Invalid Input!")
        return
    except:
        print("An Error Has Occured")
        return
    else:
        print(result)
        


#RUN:
main()
