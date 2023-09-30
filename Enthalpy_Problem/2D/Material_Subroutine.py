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
import Processed_Inputs as mm
'''
#############################################################################################################################
Class Material_model: -
=============================================================================================================================
This is a parent class and hence the parameter related methods have parent word
The methods in the parent class are only called by the child classes 
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    b,M,N       ->  Material specific vectors 
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
        self.b = np.zeros(4).reshape(-1,1)
        self.M = 0
        self.N = 0
        self.Ele = Mesh
        self.runcount = runcount

    def get_parent_param(self,element):
        b_dN = self.Ele.dN[1].reshape(-1,1)
        #b_dN Used in place of dN only for the b vector as it should contain only entries of dN at drichilet boundaries
        index = 0
        for i in self.Ele.node_list[element]:          
            if i not in self.Ele.boundary_nodes['Bottom_Nodes'] and i not in self.Ele.boundary_nodes['Right_Nodes']:
                    b_dN[index] = 0
            index = index + 1
        sum = np.sum(b_dN)
        b = b_dN*sum 
        M = np.matmul(self.Ele.N,self.Ele.N.reshape(1,-1))*self.Ele.det_J
        N = np.matmul(self.Ele.dN.T,self.Ele.dN)*self.Ele.det_J
        return[b,M,N]

    def Quad_Integ(self,ele_no):
        gpts = [-0.577350269189626, 0.577350269189626]
        check_J = 0
        for zeta in gpts:
            for eta in gpts:
                self.Ele.update_element(zeta,eta,ele_no)
                [b,M,N] = self.get_parent_param(ele_no)
                self.b = self.b + b
                self.M = self.M + M
                self.N = self.N + N
                check_J = check_J + self.Ele.det_J                                 #Testing Jacobi 
        self.b = -mm.Ta*self.b
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

    def get_param(self,He_1,ele_no):
        self.Ele.update_element(ele_no)
        self.b = 0
        self.M = 0
        self.N = 0
        self.Quad_Integ(ele_no)
        self.dT_dH = np.array([[1/mm.D if He_1[0] <= 0 else mm.Tliq/(mm.D*mm.Tliq + mm.D*mm.lambf) if np.logical_and(0 <= He_1[0], He_1[0] <= mm.Hliq) else 1/mm.D, 0],
                      [0,1/mm.D if He_1[1] <= 0 else mm.Tliq/(mm.D*mm.Tliq + mm.D*mm.lambf) if np.logical_and(0 <= He_1[1], He_1[1] <= mm.Hliq) else 1/mm.D]])
        self.runcount +=1
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
-----------------------------------------------------------------------------------------------------------------------------
Methods: -
=============================================================================================================================
    get_param -> generates the material Parameters
#############################################################################################################################
'''
class Crystal_model(Material_model):
    def __init__(self,Mesh,runcount):
        super().__init__(Mesh,runcount)
        if self.runcount<1:
            print('Crystalline routine activated in Material subroutine')
            self.runcount +=1

    def get_param(self,He_1,ele_no):
        self.b = np.zeros(4).reshape(-1,1)
        self.M = 0
        self.N = 0
        self.Quad_Integ(ele_no)
        if He_1[0] > 0 and He_1[0] < mm.D*mm.lambf:
                k = 0
        else:
                k = 1/mm.D
        if He_1[1] > 0 and He_1[1] < mm.D*mm.lambf:
                p = 0
        else:
                p = 1/mm.D
        if He_1[2] > 0 and He_1[2] < mm.D*mm.lambf:
                q = 0
        else:
                q = 1/mm.D
        if He_1[3] > 0 and He_1[3] < mm.D*mm.lambf:
                r = 0
        else:
                r = 1/mm.D
        self.dT_dH = np.array([[k, 0, 0, 0],
                               [0, p, 0, 0],
                               [0, 0, q, 0],
                               [0, 0, 0, r]])
        self.runcount +=1