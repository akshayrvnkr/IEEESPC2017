import numpy as np
def getcost(x,bpm):
    L = int(round(len(x)/bpm))
    thresh=bpm/5;
    C=np.zeros(bpm+thresh,1)





def assigncost(x,C,bpm):
    mincost=float('inf')
    for i in range(len(C)-bpm,0,-1):
        cost=C[i]+sum(pow(x[i:len(x)],2))
        if(cost<mincost):
            mincost=cost










