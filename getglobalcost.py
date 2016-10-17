import numpy as np
def getglobalcost(x,bpm):
	T0 = bpm
	L = int(round(len(x)/T0))
	thresh = T0/5

	C = float('inf') + np.zeros((L,len(x)))
	B = float('inf') + np.zeros((L,len(x)))
	
	ll = [0] * L
	ul = [0] * L
	ll[0] = 1
	ul[0] = int((float(len(x))/L) + thresh)
	for i in range(1,L):
			ll[i] = int(round(max(1,ll[i-1]+T0-thresh)))
			ul[i] = len(x)-1
	for i in range(1,L):
		ll[i] = int(round(max(1.0,ll[i-1]+T0-thresh)))
		ul[i] = len(x)-1
		for K in range(0,ul[0]+1):
			n = range(0,K+1)
			t = np.tile(np.array(x[K]),len(n))
			n1 = n == np.tile(np.array([K]),len(n))
			t2 = n1.astype(int)
			C[0][K] = sum(pow((x[0:K+1]-(t*t2)),2))
			if L==1:
					n = range(0,len(x))
					t = np.tile(np.array(x[K]),len(n))
					n1 = n == np.tile(np.array([K]),len(n))
					t2 = n1.astype(int)
					C[0][K] = sum(pow(x[0:K+1]-(t*t2),2))
			B[0][K] = 0
		for l in range(1,L):
			for K in range(ll[l],ul[l]+1):
				for m in range(int(round(max(1.0,K-T0-thresh))),int(round(min(len(x),K-T0+thresh)))+1):
					if m!=K:
						locost = C[l-1][m] + getlocalcost(x,m,K,K)
						if l == L-1:
							locost = C[l-1][m] + getlocalcost(x,m,K,len(x))
						if locost < C[l][K]:
							C[l][K] = locost
							B[l][K] = m
	pos = [0] * L
	pos[0] = np.argmin((C[L-1][:]))
	ctr = 1
	for i in range(L-1,0,-1):
		pos[ctr] = B[i][pos[ctr-1]]
		print (pos[ctr])
	ctr = ctr + 1
	pos = pos[::-1]

def getlocalcost(x,gciloc1,gciloc2,End):
	n = range(gciloc1+1,End)
	t1 = np.tile(np.array(x[gciloc2]),len(n))
	n1 = n == np.tile(np.array([gciloc2]),len(n))
	t2 = n1.astype(int)
	locost = sum(pow((x[gciloc1:End]-(t1*t2)),2))
	return locost

x = [1,2,3,4,5]
getglobalcost(x,1)
