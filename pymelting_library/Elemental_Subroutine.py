# Elemental Subroutine

import numpy as np
import inputs
from Material_Subroutine import Amorphus_model as Amor_routine

#e_1 is the the is the subscript used for ***new NRS*** vector and 
# e is used for the previous ***timestep*** vector 

class Elemental_Subroutine():
    def __init__(self):
        self.A = None
        self.Gg = np.zeros(inputs.n).reshape(inputs.n,1)
        self.dGg = np.zeros(inputs.n).reshape(inputs.n,1)

    def update_Aly(self,i):             # here the n represents the total number of elements and p represents the current element
        self.A = np.zeros((inputs.n,2))   # Assembly matrix
        self.A[i,0] = 1
        self.A[i+1,1] = 1

    def get_param(self,Hk,Hk_1,Tk_1,Mesh):
        SS304L = Amor_routine(Mesh)
        for i in range(inputs.n-1):
            self.update_Aly(i)
            He = np.matmul(self.A.T,Hk)
            He_1 = np.matmul(self.A.T,Hk_1)
            Te_1 = np.matmul(self.A.T,Tk_1)
            SS304L.get_param(He_1,Te_1,He,i)
            self.Gg = self.Gg + np.matmul(self.A,SS304L.G)
            self.dGg = self.dGg + np.matmul(np.matmul(self.A,SS304L.dG),self.A.T)
