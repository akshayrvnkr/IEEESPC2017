import numpy as np

def getperiod(acf=None, wv=None, step=None, pmin=None, pmax=None):
    rcf = np.zeros(len(acf))
    numelem = 4
    for i in range(int(pmin-1),int(pmax - 1-1)):        # maximum beat period
        for a in range(0,numelem-1):            # number of comb elements
            for b in range(1 - a,a - 1):                # gs using normalization of comb elements
                rcf[i] = rcf[i] + (acf[a * i + b] * wv[i]) / (2 * a - 1)
    return rcf