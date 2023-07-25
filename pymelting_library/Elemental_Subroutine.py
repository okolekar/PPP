'''Topic: PPP
Library: - Elemental Subroutine
#############################################################################################################################
Importing the required libraries: -
-----------------------------------------------------------------------------------------------------------------------------
inputs -> Script where all the inputs are defined
Material_Subroutine -> All material specific equations and parameters are processed here'''
#############################################################################################################################
import numpy as np
import inputs
from Material_Subroutine import Amorphus_model as Amor_routine
'''
#############################################################################################################################
Class Elemental_Subroutine: -
=============================================================================================================================
g denotes global vector, e denotes element specific vector
e_1 is the is the subscript used for current NRS vector and e is used for the previous timestep vector
(As a note, here we don't need the values of the vector at previous NRS step)
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    A       ->      Stores the assembly matrix for the specific element
    Gg      ->      The global G vector (Similar to the Residual Vector)
    dGg     ->      The global dG vector (Similar to the Tangent Stiffness Matrix)
-----------------------------------------------------------------------------------------------------------------------------
Methods : -
=============================================================================================================================
update_Aly  ->      Builds the Assembly matrix based on the element number
get_param   ->      This method extracts the element specific vectors and feeds it to the Material Subroutine.
                    It expects element specific G and dG vector for all the elements, and formualtes the global G and dG 
                    vectors to return it to the FEA Class.'''
#############################################################################################################################
class Elemental_Subroutine():
    def __init__(self):
        self.A = None
        self.Gg = np.zeros(inputs.n).reshape(inputs.n,1)
        self.dGg = np.zeros(inputs.n).reshape(inputs.n,1)

    def update_Aly(self,i):               # here the n represents the total number of elements and p represents the current element
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
