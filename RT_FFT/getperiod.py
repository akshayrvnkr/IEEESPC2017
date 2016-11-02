import numpy as np


def getperiod(acf=None, wv=None, timesig=None, step=None, pmin=None, pmax=None):
    rcf = np.zeros(step)

    if (not timesig):    # timesig unknown, must be general state
        numelem = 4

        for i in range(pmin-1,pmax - 1-1):        # maximum beat period
            for a in range(0,numelem-1):            # number of comb elements
                for b in range(1 - a,a - 1):                # gs using normalization of comb elements
                    rcf[i] = rcf[i] + (acf[a * i + b] * wv[i]) / (2 * a - 1)

    else:
        numelem = timesig    # timesig known must be context dependent state
        for i in range(pmin-1,pmax - 1-1):        # maximum beat period
            for a in range(0,numelem-1):            # number of comb elements
                for b in range(1 - a,a - 1):                # cds not normalizing comb elements
                    rcf[i]= rcf[i] + acf[a * i + b] * wv[i]