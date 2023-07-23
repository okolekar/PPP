import numpy as np
import math as m
import inputs as ip 

def check_arguments(func):
    def wrapper(self, *args,**kwargs):
        if len(args) != 1:
            raise ValueError("Incorrect number of input arguments to Initial_Temp.")
        return func(self,*args,**kwargs)
    return wrapper

class Initial_Temp():
    'The class generates the initial temperature distribution'

    @check_arguments
    def __init__(self,list):
        self.T = np.zeros(ip.n)
        self.Tdist_Enthalpy(list)

    def check_T(self,Ta):
        if self.T[0]<Ta:
            print(self.T)
            print("Temperature distribution might be invalid.")
            print("As the Temperature at node 1 is lower than the Ta.")

    def Tdist_Enthalpy(self,list):
        F = ip.ratioI
        for i in range(len(self.T)):
            if i == len(self.T)-1:
                self.T[i] = ip.Ta #This takes care of the boundary condition
            else:
                self.T[i] = (F*(2*np.sqrt(ip.delt/np.pi)
                               *np.exp(-list[i]**2/(4*ip.delt))
                                -list[i]*m.erfc(list[i]/(2*np.sqrt(ip.delt))))
                                  + ip.Ta)
        self.T = self.T[:,np.newaxis]
        self.check_T(ip.Ta)
