import numpy as np
import math as m

class Initial_Temp():
    'This class throws back a vector of temperature. The input required are I/Iref, Ta, delta_t and node list'
    def __init__(self, ratioI, Ta, delt, n):
        self.F = ratioI
        self.Ta = Ta
        self.dt = delt
        self.list = n

    def Tdist_Enthalpy(self):
        T = np.zeros(len(self.list))
        for i in range(len(T)):
            if i == len(T)-1:
                T[i] = self.Ta
            else:
                T[i] = self.F*(2*np.sqrt(self.dt/np.pi)*np.exp(-self.list[i]**2/(4*self.dt))-self.list[i]*m.erfc(self.list[i]/(2*np.sqrt(self.dt)))) + self.Ta
        return T
