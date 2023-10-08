'''
Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
Library: - Material Subroutine
#############################################################################################################################
Importing the required libraries : -
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs -> Script where all the inputs are defined
#############################################################################################################################
'''
import numpy as np
import Processed_Inputs as mm

class Material_model():
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
    Quad_Integ       -> Performs the Gaussian Integration also tests if the Jacobi calculated is correct
#############################################################################################################################
    '''
    def __init__(self,Mesh,runcount):
        self.b = 0
        self.M = 0
        self.N = 0
        self.Ele = Mesh
        self.runcount = runcount

    def get_parent_param(self,element):
        if element == 0:
            b = np.array([[mm.ratioI],
                          [0]])
        elif element == mm.n-2:
            b = self.Ele.dphi[1]*self.Ele.dphi*self.Ele.J
        else:
            b = np.array([[0],
                          [0]])
        M = np.matmul(self.Ele.phi,self.Ele.phi.reshape(1,-1))*self.Ele.J
        N = np.matmul(self.Ele.dphi,self.Ele.dphi.reshape(1,-1))*self.Ele.J
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
            check_J = check_J + self.Ele.J                                 #Testing Jacobi 
            if ele_no == mm.n-2:
                    self.b = -mm.Ta*self.b
        if check_J != self.Ele.nl[ele_no+1] - self.Ele.nl[ele_no]:
            raise ValueError("Jacobi after Gauss summation is not equal to element length")

class Amorphus_model(Material_model):
    '''
#############################################################################################################################
Child Class Amorphus_model: -
=============================================================================================================================
Inherited by the Material_model
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    dT_dH   ->  Temperature Enthalpy slope matrix 
-----------------------------------------------------------------------------------------------------------------------------
Other Variables: - 
------------------
    He_1 ->  corresponds to the current enthalpy from the NRS scheme
-----------------------------------------------------------------------------------------------------------------------------
Methods: -
=============================================================================================================================
    get_param -> generates the material specific vector/materix. It calculates slope function dT/dH matrix.
#############################################################################################################################
    '''
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

class Crystal_model(Material_model):
    '''
#############################################################################################################################
Child Class Crystal_model: -
=============================================================================================================================
Inherited by the Material_model
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    dT_dH   ->  Temperature Enthalpy slope matrix
-----------------------------------------------------------------------------------------------------------------------------
Other Variables: - 
    He_1 ->  corresponds to the current enthalpy from the NRS scheme
-----------------------------------------------------------------------------------------------------------------------------
Methods: -
=============================================================================================================================
    get_param -> generates the material specific vector/materix. It calculates slope function dT/dH matrix.
#############################################################################################################################
    '''
    def __init__(self,Mesh,runcount):
        super().__init__(Mesh,runcount)
        if self.runcount<1:
            print('Crystalline routine activated in Material subroutine')
            self.runcount +=1

    def get_param(self,He_1,ele_no):
        self.Ele.update_element(ele_no)
        self.b = 0
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
        self.dT_dH = np.array([[k, 0],
                               [0, p]])
        self.runcount +=1