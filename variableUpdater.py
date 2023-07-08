import numpy as np

def Update(D,lambf,T=None,H=None):
    if T is not None:
        return np.where(T < 0, D*T, np.where(T == 0,D*lambf, D*T+D*lambf))
    elif H is not None:
        return np.where(H < 0, H/lambf, np.where(H <= lambf*D,0, (H-D*lambf)/D))
    else:
        print("Problem in the program detected. Terminating further execution")