import numpy as np
import math as m
class initial_condition():
    def __init__(self,ratioI, delt, Ta,s = 0.01):

        self.ratioI = ratioI
        self.delt = delt
        self.Ta = Ta
        self.s = s

    def Optimize(self):
        result = self.ratioI*(2*np.sqrt(self.delt/np.pi)*np.exp(-self.s**2/(4*self.delt))-self.s*m.erfc(self.s/(2*np.sqrt(self.t))))+self.Ta
        while result > 10**-5:
            #Still need to code
            print('work under process')
            break            
        return self.s