'''Topic: PPP
Library: - Element2D
#############################################################################################################################
Importing the required standard libraries
-----------------------------------------------------------------------------------------------------------------------------
Material_model2D    ->  The material model library
#############################################################################################################################
'''
import numpy as np
from Material_model2D import material_model as model
'''
#############################################################################################################################
Class Element_Subroutine: -
=============================================================================================================================
The attributes of this class store all the Assembled material matrices and vectors
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    K   ->  Conditivity matrix
    N   ->  Heat Energy matrix
    Vx  ->  Velocity vector in x direction
    Vy  ->  Velocity vector in y direction
-----------------------------------------------------------------------------------------------------------------------------
Variables: -
------------
    M   ->  object of the material model 
_____________________________________________________________________________________________________________________________
Main Methods : -
----------------
=============================================================================================================================
formGlobalParm  ->  calls the material model and updates material specific matrices and parametes for every element and then 
                    assembles them to form the global vectors and matrices.
#############################################################################################################################
'''
class Element_Subroutine():
    def __init__(self,n):
        self.K  = np.zeros([n,n])
        self.N  = np.zeros([n,n])
        self.Vx = np.zeros(n)
        self.Vy = np.zeros(n)
    '''
#############################################################################################################################
Method formGlobalParm(self,mesh): -
=============================================================================================================================
Arguments: -
------------
    mesh    ->  object of the class Mesh
-----------------------------------------------------------------------------------------------------------------------------
The method forms the global parameters which is solved by the solver.
#############################################################################################################################
'''
    def formGlobalParm(self,mesh):
        M = model()
#=====================================================#Elemental Loop#======================================================#
        for i in mesh.elementList:
            M.update_attributes(mesh,i)     #updating the material attributes      
#=================================================#Assembling the matrices#=================================================#
            
            for q in range(4):
#=================================================#Assembling the M matrix#=================================================#
                self.K[mesh.node_list[i][q],mesh.node_list[i][0]] =self.K[mesh.node_list[i][q],mesh.node_list[i][0]]+M.k[q,0]
                self.K[mesh.node_list[i][q],mesh.node_list[i][1]] =self.K[mesh.node_list[i][q],mesh.node_list[i][1]]+M.k[q,1]
                self.K[mesh.node_list[i][q],mesh.node_list[i][2]] =self.K[mesh.node_list[i][q],mesh.node_list[i][2]]+M.k[q,2]
                self.K[mesh.node_list[i][q],mesh.node_list[i][3]] =self.K[mesh.node_list[i][q],mesh.node_list[i][3]]+M.k[q,3]
#=================================================#Assembling the N matrix#=================================================#
                self.N[mesh.node_list[i][q],mesh.node_list[i][0]] =self.N[mesh.node_list[i][q],mesh.node_list[i][0]]+M.n[q,0]
                self.N[mesh.node_list[i][q],mesh.node_list[i][1]] =self.N[mesh.node_list[i][q],mesh.node_list[i][1]]+M.n[q,1]
                self.N[mesh.node_list[i][q],mesh.node_list[i][2]] =self.N[mesh.node_list[i][q],mesh.node_list[i][2]]+M.n[q,2]
                self.N[mesh.node_list[i][q],mesh.node_list[i][3]] =self.N[mesh.node_list[i][q],mesh.node_list[i][3]]+M.n[q,3]
#=================================================#Assembling the Vx vector#================================================#
            self.Vx[mesh.node_list[i][0]] = M.Vx[0]
            self.Vx[mesh.node_list[i][1]] = M.Vx[1]
            self.Vx[mesh.node_list[i][2]] = M.Vx[2]
            self.Vx[mesh.node_list[i][3]] = M.Vx[3]
#=================================================#Assembling the Vy vector#================================================#
            self.Vy[mesh.node_list[i][0]] = M.Vy[0]
            self.Vy[mesh.node_list[i][1]] = M.Vy[1]
            self.Vy[mesh.node_list[i][2]] = M.Vy[2]
            self.Vy[mesh.node_list[i][3]] = M.Vy[3]
#================================================#End of the Elemental Loop#================================================#
'''
=============================================================================================================================
--------------------------------------------------------End of Class---------------------------------------------------------
=============================================================================================================================
'''