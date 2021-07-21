#imports:
import sys
import enum
import numpy as np

#"goal" enum:
class goalEnum (enum.Enum):
    spk = "spk"
    wam = "wam"
    ddg = "ddg"
    lnorm = "lnorm"
    jacobi = "jacobi"

#methods:
def input():
    args = sys.argv
    if len(args)!=4:
        raise("Invalid Input!")
    k,goal,file_name = (args[1],args[2],args[3])
    try:
        k = int(k)
        goal = goalEnum(goal)
        file_name = file_name
    except ValueError:
        raise("Invalid Input!")
    else:
        return(k,goal,file_name)


#main method:
def main():
    np.random.seed(0)
    k, goal, file_name = input()
    print (k,goal,file_name)
#RUN:
main()
