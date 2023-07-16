import numpy as np
from slope_matrix import temp_der as dT
import inputs as mm
from Mesh import Mesh 

#Here He_1 corresponds to the updated enthalpy from the NRS
#And the He corresponds to the enthalpy at the previous time step not NRS scheme.

class Material_model():
    def __init__(self):
        self.b = 0
        self.M = 0
        self.N = 0
        self.F = 0
        self.G = 0
        self.Ele = Mesh(mm.n)

    def get_parent_param(self,boundry_node):
        
        if boundry_node == 0:
            b = np.array([[mm.ratioI],
                            [0]])
        elif boundry_node == 1:
            b = - mm.Ta*self.Ele.dphi[1]*self.Ele.dphi
        else:
            b = np.array([[0],
                          [0]])
            
        M = np.matmul(self.Ele.phi,self.Ele.phi.reshape(1,2))
        N = np.matmul(self.Ele.dphi,self.Ele.dphi.reshape(1,2))
        return[b,M,N]

    def Quad_Integ(self,ele_no):
        gpts = [-0.577350269189626, 0.577350269189626]
        check_J = 0
        for zeta in gpts:
            self.Ele.update_phi(zeta)
            [b,M,N] = self.get_parent_param(ele_no)
            if ele_no != 0:
                self.b = self.b + b
            else:
                self.b = b
            self.M = self.M + M
            self.N = self.N + N
            check_J = check_J + self.Ele.J
        if check_J != self.Ele.nl[ele_no+1] - self.Ele.nl[ele_no]:
            raise ValueError("Jacobi after Gauss summation is not equal to element length")

class Amorphus_model(Material_model):
    def __init__(self):
        super().__init__()
        self.dG = None

    def get_param(self,He_1,Te_1,He,ele_no):
        self.Ele.update_element(ele_no)
        self.Quad_Integ(ele_no)
        slope = dT(mm.D, He_1,mm.Tliq,mm.lambf,mm.Hliq)
        self.F = mm.D*(self.b-np.matmul(self.N,Te_1))
        self.G = np.matmul(self.M,(He_1-He)) - mm.delt*self.F
        self.dG = self.M + mm.delt*mm.D*np.matmul(self.N,slope.dT_dH)

class Crystal_model(Material_model):
    def __init__(self):
        super().__init__()
        self.dG = None

    def get_param(self,D,He_1,delt):
        slope = dT.Cryst_slop(D, He_1)
        self.dG = self.M + delt*D*np.matmul(self.N,slope.dT_dH)