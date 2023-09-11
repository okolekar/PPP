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
from Stefan_Material_model import Stefan_Material_Model as SMat_routine
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
    Enthalpy Approach: -
    ---------------------
        A           ->      Stores the assembly matrix for the specific element
        runcount    ->      Records the number of times the methods of this class was called
        Gg          ->      The global G vector (Similar to the Residual Vector)
        dGg         ->      The global dG vector (Similar to the Tangent Stiffness Matrix)
    Stefan Approach: -
    ------------------
        Mlg         ->      Global Mass matrix for the liquid domain
        Nlg         ->      Material specific matrix for the liquid domain 
        blg         ->      Material specific column vector for the liquid domain
        Msg         ->      Global Mass matrix for the solid domain
        Nsg         ->      Material specific matrix for the solid domain 
        bsg         ->      Material specific column vector for the solid domain
-----------------------------------------------------------------------------------------------------------------------------
Input Arguments: -
------------------
    approach    ->      Indicates if Stefan Approach or enthalpy Approach is to be selected
-----------------------------------------------------------------------------------------------------------------------------
Methods : -
=============================================================================================================================
update_Aly              ->  Builds the Assembly matrix based on the element number
get_param               ->  This method extracts the element specific vectors and feeds it to the Material Subroutine.
                            It expects element specific G and dG vector for all the elements, and formualtes the global  
                            G and dG vectors to return it to the FEA Class.
get_global_Stefan_param ->  The method Updates the Jacobian for each element and passes the element current number to the  
                            Stefan Material Model to formulate the element specific matrices and  assembles them to form the
                            global materices which are later accessed by the FEM class. 
#############################################################################################################################
'''
class Elemental_Subroutine():
    def __init__(self,approach = None):
        self.A          = None
        self.runcount   = 0
        if approach == None:
            self.Gg     = np.zeros(inputs.n).reshape(inputs.n,1)
            self.dGg    = np.zeros(inputs.n).reshape(inputs.n,1)
        elif approach == 's':
            self.Mlg = np.zeros([inputs.n,inputs.n])
            self.Nlg = np.zeros([inputs.n,inputs.n])
            self.blg = np.zeros(inputs.n).reshape(inputs.n,1)
            self.Msg = np.zeros([inputs.m,inputs.m])
            self.Nsg = np.zeros([inputs.m,inputs.m])
            self.bsg = np.zeros(inputs.m).reshape(inputs.m,1)
#===========================================================================================================================#
                                            #Assembly matrix updater.
#===========================================================================================================================#
    def update_Aly(self,i,approach=None):
        if approach == 's':
            n = inputs.m
        else:
            n = inputs.n
        self.A          = np.zeros((n,2))   
        self.A[i,0]     = 1
        self.A[i+1,1]   = 1
#===========================================================================================================================#
                                        #Enthalpy Approach parameters maker.
#===========================================================================================================================#
    def get_param(self,Hk,Hk_1,Tk_1,Mesh):
        Mat = Mat_routine(Mesh,self.runcount)
        for i in range(inputs.n-1):
            self.update_Aly(i)
            He      = np.matmul(self.A.T,Hk)
            He_1    = np.matmul(self.A.T,Hk_1)
            Te_1    = np.matmul(self.A.T,Tk_1)
            Mat.get_param(He_1,Te_1,He,i)
            self.runcount += 1
            self.Gg     = self.Gg + np.matmul(self.A,Mat.G)
            self.dGg    = self.dGg + np.matmul(np.matmul(self.A,Mat.dG),self.A.T)
#===========================================================================================================================#
                                        #Stefan Approach parameters maker.
#===========================================================================================================================#
    def get_global_Stefan_param(self,meshl,meshs,ds_dt):
        material = SMat_routine()
#---------------------------------------------------------------------------------------------------------------------------#
                                                #liquid regime
#---------------------------------------------------------------------------------------------------------------------------#
        for i in range(inputs.n-1):
            meshl.update_Jacobi(i)
            material.get_liq_param(meshl,ds_dt,i)
            self.update_Aly(i)
            self.Mlg = self.Mlg + np.matmul(np.matmul(self.A,material.Ml),self.A.T)
            self.Nlg = self.Nlg + np.matmul(np.matmul(self.A,material.Nl),self.A.T)
            self.blg = self.blg + np.matmul(self.A,material.bl)
#---------------------------------------------------------------------------------------------------------------------------#
                                                #solid regime
#---------------------------------------------------------------------------------------------------------------------------#
        for i in range(inputs.m-1):
            meshs.update_Jacobi(i)
            material.get_sol_param(meshs,i,ds_dt)
            self.update_Aly(i,'s')
            self.Msg = self.Msg + np.matmul(np.matmul(self.A,material.Ms),self.A.T)
            self.Nsg = self.Nsg + np.matmul(np.matmul(self.A,material.Ns),self.A.T)
            self.bsg = self.bsg + np.matmul(self.A,material.bs)