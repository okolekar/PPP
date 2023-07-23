#Heat_Transfer_verification
import numpy as np
import inputs as ip

class Verification():
    def __init__(self,Tinitial, n,nlist,alpha):
        self.T_analytical = np.zeros(ip.n).reshape(n,1)
        for i in range(n):
            x = nlist[i]
            self.T_analytical[i] = 2*(Tinitial[i])*(np.cos(np.pi*x))*np.exp((-np.pi**2*alpha*ip.delt))
