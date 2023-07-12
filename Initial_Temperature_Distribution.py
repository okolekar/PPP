import numpy as np
import math as m 

def check_arguments(func):
        def wrapper(self, *args,**kwargs):
            if len(args) != 4:
                raise ValueError("Incorrect number of input arguments to Initial_Temp.")
            return func(self,*args,**kwargs)
        return wrapper

class Initial_Temp():
    'The class generates the initial temperature distribution'

    @check_arguments
    def __init__(self, ratioI, Ta, delt, n):
        if ratioI < 0:
            print("Task suspended at Initial_Temperature_Dist __init__")
            raise ValueError("ratioI cannot be negative")
        F = ratioI
        
        try:
            if len(n) < 2:
                print("Task suspended at Initial_Temperature_Dist __init__")
                raise ValueError("Mesh list cannot have just one node")
            list = n
        except TypeError:
            print("Mesh list cannot be scalar")
        
        try:
            self.T = np.zeros(len(list))
            self.Tdist_Enthalpy(F,Ta,delt,list)
        except UnboundLocalError:
            print("Scalar node list detected")

    def check_T(self,Ta):
        if self.T[0]<Ta:
            print(self.T)
            print("Temperature distribution might be invalid.")
            print("As the Temperature at node 1 is lower than the Ta.")

    def Tdist_Enthalpy(self,F,Ta,dt,list):
        for i in range(len(self.T)):
            if i == len(self.T)-1:
                self.T[i] = Ta #This takes care of the boundary condition
            else:
                self.T[i] = (F*(2*np.sqrt(dt/np.pi)
                               *np.exp(-list[i]**2/(4*dt))
                                -list[i]*m.erfc(list[i]/(2*np.sqrt(dt))))
                                  + Ta)
        self.check_T(Ta)
