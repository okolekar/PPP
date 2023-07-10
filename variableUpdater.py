import numpy as np

def Update(D,lambf,T=None,H=None):
    if T is not None:
        return np.where(T < 0, D*T, np.where(T == 0,D*lambf, D*T+D*lambf))
    elif H is not None:
        return np.where(H < 0, H/D, np.where(H <= lambf*D,0, (H-D*lambf)/D))
    else:
        print("Problem in the program detected. Terminating further execution")

def Update_amorphus(D,lambf,Tliq,T=None,H=None,Hliq=None):
    if T is not None:
        return np.where(T <= 0, D*T, np.where(np.logical_and(T >= 0, T<= Tliq),D*lambf*T/Tliq, D*T+D*lambf))
    elif H is not None:
        if Hliq is None:
            print("Hliq is not provided. Value override")
        else:
            return np.where(H <= 0, H/D, np.where(H <= Hliq,H*Tliq/(D*Tliq+D*lambf), (H-D*lambf)/D))
    else:
        print("Problem in the program detected. Terminating further execution")
