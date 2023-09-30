'''Topic: PPP
Library: - Elemental Subroutine
#############################################################################################################################
Importing the required libraries: -
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs    -> Script where all the inputs are defined
Material_Subroutine -> All material specific equations and parameters are processed here
#############################################################################################################################
'''
import numpy as np
import Processed_Inputs as inputs
if inputs.mat_type == 1:
    from Material_Subroutine import Crystal_model as Mat_routine
    print('Crystalline material detected @ Elemental Subroutine')
elif inputs.mat_type == 2:
    from Material_Subroutine import Amorphus_model as Mat_routine
    print('Amorphous material detected @ Elemental Subroutine')
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
    A           ->      Stores the assembly matrix for the specific element
    Gg          ->      The global G vector (Similar to the Residual Vector)
    dGg         ->      The global dG vector (Similar to the Tangent Stiffness Matrix)
    runcount    ->      Records the number of times the methods of this class was called
-----------------------------------------------------------------------------------------------------------------------------
Methods : -
=============================================================================================================================
update_Aly  ->      Builds the Assembly matrix based on the element number
get_param   ->      This method extracts the element specific vectors and feeds it to the Material Subroutine.
                    It expects element specific G and dG vector for all the elements, and formualtes the global G and dG 
                    vectors to return it to the FEA Class.
#############################################################################################################################
'''
class Elemental_Subroutine():
    def __init__(self):
        self.Gg = np.zeros(inputs.n).reshape(inputs.n,1)
        self.dGg = np.zeros(inputs.n).reshape(inputs.n,1)
        self.runcount = 0

    def get_param(self,Hk,Hk_1,Tk_1,Mesh):
        Mat = Mat_routine(Mesh,self.runcount)
        self.b = np.zeros(inputs.n).reshape(-1,1)
        self.M = np.zeros([inputs.n,inputs.n])
        self.N = np.zeros([inputs.n,inputs.n])
        self.dT_dH = np.zeros([inputs.n,inputs.n])
        for i in range(inputs.n-1):
            He_1 = Hk_1[i:i+2].reshape(-1,1)
            Mat.get_param(He_1,i)
            self.runcount += 1
#---------------------------------------------------------------------------------------------------------------------------#
            self.M[i,i] = self.M[i,i] + Mat.M[0,0]
            self.M[i,i+1] = self.M[i,i+1] + Mat.M[0,1]
            self.M[i+1,i] = self.M[i+1,i] + Mat.M[1,0]
            self.M[i+1,i+1] = self.M[i+1,i+1] + Mat.M[1,1]
#---------------------------------------------------------------------------------------------------------------------------#
            self.N[i,i] = self.N[i,i] + Mat.N[0,0]
            self.N[i,i+1] = self.N[i,i+1] + Mat.N[0,1]
            self.N[i+1,i] = self.N[i+1,i] + Mat.N[1,0]
            self.N[i+1,i+1] = self.N[i+1,i+1] + Mat.N[1,1]
#---------------------------------------------------------------------------------------------------------------------------#
            self.dT_dH[i,i] = self.dT_dH[i,i] + Mat.dT_dH[0,0]
            self.dT_dH[i,i+1] = self.dT_dH[i,i+1] + Mat.dT_dH[0,1]
            self.dT_dH[i+1,i] = self.dT_dH[i+1,i] + Mat.dT_dH[1,0]
            self.dT_dH[i+1,i+1] = self.dT_dH[i+1,i+1] + Mat.dT_dH[1,1]
#---------------------------------------------------------------------------------------------------------------------------#
        self.F = inputs.D*(self.b-np.matmul(self.N,Tk_1))
        self.Gg = np.matmul(self.M,(Hk_1-Hk)) - inputs.delt*self.F
        self.dGg = self.M + inputs.delt*inputs.D*np.matmul(self.N,self.dT_dH)
