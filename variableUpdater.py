import numpy as np

def check_arguments(func,*args,**kwargs):
    def wrapper(*args,**kwargs):
        if func.__name__ == 'Update_Crys' and len(kwargs) != 4:
            raise ValueError("Incorrect number of input arguments to Update_Crys.")
        elif func.__name__ == 'Update_amorphus' and len(kwargs) != 6:
            raise ValueError("Incorrect number of input arguments to Update_amorphus.")
        return func(*args,**kwargs)
    return wrapper

@check_arguments
def Update_Crys(D,lambf,T=None,H=None):
    if T is not None:
        return np.where(T < 0, D*T, np.where(T == 0,D*lambf, D*T+D*lambf))
    elif H is not None:
        return np.where(H < 0, H/D, np.where(H <= lambf*D,0, (H-D*lambf)/D))
    else:
        print("Update_Crys() failed to execute")
        raise ValueError("Enthalpy and/or Temperature not provided. Hence, terminating further execution")
    
@check_arguments
def Update_amorphus(D,lambf,Tliq,T=None,H=None,Hliq=None):
    if T is not None:
        return np.where(T <= 0, D*T, np.where(np.logical_and(T >= 0, T<= Tliq),D*lambf*T/Tliq, D*T+D*lambf))
    elif H is not None:
        return np.where(H <= 0, H/D, np.where(H <= Hliq,H*Tliq/(D*Tliq+D*lambf), (H-D*lambf)/D))
    else:
        print("Update_amorphus() failed to execute")
        raise ValueError("Enthalpy and/or Temperature not provided. Hence, terminating further execution")
