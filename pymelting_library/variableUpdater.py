'''Topic: PPP
Library: - Variable Updater
#############################################################################################################################
Importing the required standard libraries   
-----------------------------------------------------------------------------------------------------------------------------
inputs -> Script where all the inputs are defined                                                                         '''
#############################################################################################################################
import numpy as np
import inputs as ip
'''
#############################################################################################################################
Check arguments Decorator: -
-----------------------------------------------------------------------------------------------------------------------------
Checks if atleast 1 argument is passed to this Class                                                                         
#############################################################################################################################
'''
def check_arguments(func,*args,**kwargs):
    def wrapper(*args,**kwargs):
        if func.__name__ == 'Update_Crys' and len(kwargs) > 2:
            raise ValueError("Incorrect number of input arguments to Update_Crys.")

        elif func.__name__ == 'Update_amorphus' and (len(kwargs) > 3):
            raise ValueError("Incorrect number of input arguments to Update_amorphus.")
        return func(*args,**kwargs)
    return wrapper
'''
#############################################################################################################################
Function Update_Crys: -
=============================================================================================================================
Based on the inputs, if Temperature is provided then Enthalpy is updated based on the Temperature, as follows: -
-----------------------------------------------------------------------------------------------------------------------------
If T < 0 --> H = D * T                  -> The material is still solid 
If T = 0 --> H = 0                      -> The material just melted
If T > 0 --> H = D * T + D * lambda     -> The material is completed melted
and vice versa if Enthalpy is provided.
#############################################################################################################################
'''
@check_arguments
def Update_Crys(T=None,H=None):
    if T is not None:
        return np.where(T < 0, ip.D*T, np.where(T == 0,ip.D*ip.lambf, ip.D*T+ip.D*ip.lambf))
    elif H is not None:
        return np.where(H < 0, H/ip.D, np.where(H <= ip.lambf*ip.D,0, (H-ip.D*ip.lambf)/ip.D))
    else:
        print("Update_Crys() failed to execute")
        raise ValueError("Enthalpy and/or Temperature not provided. Hence, terminating further execution")
'''
#############################################################################################################################
Function Update_amorphus: -
=============================================================================================================================
Based on the inputs, if Temperature is provided then Enthalpy is updated based on the Temperature, as follows: -
-----------------------------------------------------------------------------------------------------------------------------
If T < 0 --> H = D * T                          -> The material is still solid 
If T = 0 --> H = D * T + D * lambda * T/Tliq    -> The material just melted
If T > 0 --> H = D * T + D * lambda             -> The material is completed melted
and vice versa if Enthalpy is provided.
#############################################################################################################################
'''   
@check_arguments
def Update_amorphus(T=None,H=None,Hliq=None):
    if T is not None:
        return np.where(T <= 0, ip.D*T, np.where(np.logical_and(T >= 0, T <= ip.Tliq),ip.D*ip.lambf*T/ip.Tliq, ip.D*T+ip.D*ip.lambf))
    elif H is not None:
        return np.where(H <= 0, H/ip.D, np.where(H <= Hliq,H*ip.Tliq/(ip.D*ip.Tliq+ip.D*ip.lambf), (H-ip.D*ip.lambf)/ip.D))
    else:
        print("Update_amorphus() failed to execute")
        raise ValueError("Enthalpy and/or Temperature not provided. Hence, terminating further execution")
