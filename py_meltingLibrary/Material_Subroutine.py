'''Topic: PPP
Library: - Material Subroutine
#############################################################################################################################
Importing the required libraries : -
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs -> Script where all the inputs are defined
slope_matrix -> Calculates the dT/dH matrix
#############################################################################################################################
'''
import numpy as np
from slope_matrix import temp_der as dT
import Processed_Inputs as mm
'''
#############################################################################################################################
Class Material_model: -
=============================================================================================================================
This is a parent class and hence the parameter related methods have parent word associated with it
The methods in the parent class are only called by the child classes 
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    b,M,N,F     ->  Material specific vectors
    G           ->  Residual type vector 
    Ele         ->  Instance of Mesh class stores mesh related parameters passed by the FEA class
    runcount    ->  Checks if the method was already called before
-----------------------------------------------------------------------------------------------------------------------------
Other Variables: -
------------------
    He_1 ->  corresponds to the current enthalpy from the NRS scheme
    He   ->  corresponds to the enthalpy at the previous time step
-----------------------------------------------------------------------------------------------------------------------------
Methods: -
=============================================================================================================================
    get_parent_param -> generates the b,M,N vector/materices
    Quad_Integ -> Performs the Gaussian Integration also tests if the Jacobi calculated is correct
#############################################################################################################################
'''
class Material_model():
    def __init__(self,Mesh,runcount):
        self.b = 0
        self.M = 0
        self.N = 0
        self.F = 0
        self.G = 0
        self.Ele = Mesh
        self.runcount = runcount

    def get_parent_param(self,boundry_node):
        
        if boundry_node == 0:
            b = np.array([[mm.ratioI],
                          [    0    ]])
        elif boundry_node == 1:
            b = -mm.Ta*self.Ele.dphi[1]*self.Ele.dphi*self.Ele.J
        else:
            b = np.array([[0],
                          [0]])
        M = np.matmul(self.Ele.phi,self.Ele.phi.reshape(1,2))*self.Ele.J
        N = np.matmul(self.Ele.dphi,self.Ele.dphi.reshape(1,2))*self.Ele.J
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
            check_J = check_J + self.Ele.J              #Testing Jacobi 
        if check_J != self.Ele.nl[ele_no+1] - self.Ele.nl[ele_no]:
            raise ValueError("Jacobi after Gauss summation is not equal to element length")
        
'''
#############################################################################################################################
Child Class Amorphus_model: -
=============================================================================================================================
Inherited by the Material_model
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    dG   ->  Tangent Stiffness type Matrix 
-----------------------------------------------------------------------------------------------------------------------------
Other Variables: - 
------------------
    He_1 ->  corresponds to the current enthalpy from the NRS scheme
    He   ->  corresponds to the enthalpy at the previous time step
-----------------------------------------------------------------------------------------------------------------------------
Methods: -
=============================================================================================================================
    get_param -> generates the G and dG vector/materix. It uses slope function which calculates dT/dH matrix to which 
                 an 1 or 2 integer argument along with He_1 is passed. I use 2 to tell the slope library to use the 
                 amorphous method to calculate the dT/dH matrix.
#############################################################################################################################
'''
class Amorphus_model(Material_model):
    def __init__(self,Mesh,runcount):
        super().__init__(Mesh,runcount)
        if self.runcount<1:
            print('Amorphous routine activated in Material subroutine')
            self.runcount +=1
        self.dG = None

    def get_param(self,He_1,Te_1,He,ele_no):
        self.Ele.update_element(ele_no)
        self.Quad_Integ(ele_no)
        slope = dT(mm.mat_type,He_1,self.runcount)
        self.runcount +=1
        self.F = mm.D*(self.b-np.matmul(self.N,Te_1))
        self.G = np.matmul(self.M,(He_1-He)) - mm.delt*self.F
        self.dG = self.M + mm.delt*mm.D*np.matmul(self.N,slope.dT_dH)
'''
#############################################################################################################################
Child Class Crystal_model: -
=============================================================================================================================
Inherited by the Material_model
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    dG   ->  Tangent Stiffness type Matrix 
-----------------------------------------------------------------------------------------------------------------------------
Other Variables: - 
    He_1 ->  corresponds to the current enthalpy from the NRS scheme
    He   ->  corresponds to the enthalpy at the previous time step
-----------------------------------------------------------------------------------------------------------------------------
Methods: -
=============================================================================================================================
    get_param -> generates the G and dG vector/materix. It uses slope function which calculates dT/dH matrix to which 
                 an 1 or 2 integer argument along with He_1 is passed. I use 1 to tell the slope library to use the 
                 crystalline method to calculate the dT/dH matrix.
#############################################################################################################################
'''
class Crystal_model(Material_model):
    def __init__(self,Mesh,runcount):
        super().__init__(Mesh,runcount)
        if self.runcount<1:
            print('Crystalline routine activated in Material subroutine')
            self.runcount +=1
        self.dG = None

    def get_param(self,He_1,Te_1,He,ele_no):
        self.Ele.update_element(ele_no)
        self.Quad_Integ(ele_no)
        slope = dT(mm.mat_type,He_1,self.runcount)
        self.runcount +=1
        self.F = mm.D*(self.b-np.matmul(self.N,Te_1))
        self.G = np.matmul(self.M,(He_1-He)) - mm.delt*self.F
        self.dG = self.M + mm.delt*mm.D*np.matmul(self.N,slope.dT_dH)