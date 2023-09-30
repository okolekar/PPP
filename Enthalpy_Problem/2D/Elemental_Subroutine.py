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
    Gg          ->      The global G vector (Similar to the Residual Vector)
    dGg         ->      The global dG vector (Similar to the Tangent Stiffness Matrix)
    runcount    ->      Records the number of times the methods of this class was called
-----------------------------------------------------------------------------------------------------------------------------
Methods : -
=============================================================================================================================
get_param   ->      This method extracts the element specific vectors and feeds it to the Material Subroutine.
                    It expects element specific G and dG vector for all the elements, and formualtes the global G and dG 
                    vectors to return it to the FEA Class.
#############################################################################################################################
'''
class Elemental_Subroutine():
    def __init__(self):
        self.M   = None
        self.N   = None
        self.F   = None
        self.b   = None
        self.Gg  = None
        self.dGg = None
        self.runcount = 0

    def get_param(self,Hk,Hk_1,Tk_1,Mesh):
        Mat = Mat_routine(Mesh,self.runcount)
        self.b      = np.zeros(inputs.n*inputs.m).reshape(-1,1)
        self.M      = np.zeros([inputs.n*inputs.m,inputs.n*inputs.m])
        self.N      = np.zeros([inputs.n*inputs.m,inputs.n*inputs.m])
        self.dT_dH  = np.zeros([inputs.n*inputs.m,inputs.n*inputs.m])
        He_1        = np.zeros(4).reshape(-1,1)

        for i in Mesh.elementList:
            He_1[0] = Hk_1[Mesh.node_list[i][0]]
            He_1[1] = Hk_1[Mesh.node_list[i][1]]
            He_1[2] = Hk_1[Mesh.node_list[i][2]]
            He_1[3] = Hk_1[Mesh.node_list[i][3]]
            Mat.get_param(He_1,i)
            self.runcount += 1
            self.assemble(Mat,Mesh.node_list[i])

#=========================================Application of the laser on the top nodes=========================================#

        self.b[Mesh.boundary_node['Top_Nodes'][0]] = input.ratioI
        self.b[Mesh.boundary_node['Top_Nodes'][1]] = input.ratioI
#---------------------------------------------------------------------------------------------------------------------------#

        self.F = inputs.D*(self.b-np.matmul(self.N,Tk_1))
        self.Gg = np.matmul(self.M,(Hk_1-Hk)) - inputs.delt*self.F
        self.dGg = self.M + inputs.delt*inputs.D*np.matmul(self.N,self.dT_dH)

#=====================================================End of the Method=====================================================#

    def assemble(self,Mat,nodelist):
        for i in range(4):
            for j in range(4):
                self.M[nodelist[i],nodelist[j]]     = self.M[nodelist[i],nodelist[j]] + Mat.M[i,j]
                self.N[nodelist[i],nodelist[j]]     = self.N[nodelist[i],nodelist[j]] + Mat.N[i,j]
                self.dT_dH[nodelist[i],nodelist[j]] = self.dT_dH[nodelist[i],nodelist[j]] + Mat.dT_dH[i,j]
            self.b[nodelist[i]]   = self.b[nodelist[i]] + Mat.b[i]         
#---------------------------------------------------------------------------------------------------------------------------#