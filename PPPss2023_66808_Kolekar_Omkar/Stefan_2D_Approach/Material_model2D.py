'''Topic: PPP
Library: - Material_model2D
#############################################################################################################################
Importing the required standard libraries
-----------------------------------------------------------------------------------------------------------------------------
Crystalline_inputs    ->  The non dimentionalised Inputs library
#############################################################################################################################
'''
import numpy as np
import Crystalline_inputs as ip
'''
#############################################################################################################################
Class material_model: -
=============================================================================================================================
The attributes of this class form material matrices and vectors
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    k   ->  element specific Conditivity matrix
    n   ->  element specific Heat Energy matrix
    Vx  ->  element specific Velocity vector in x direction
    Vy  ->  element specific Velocity vector in y direction
-----------------------------------------------------------------------------------------------------------------------------
Variables: -
------------
    mesh       ->  Object of the mesh class
    elem_no    ->  element number
    gausspts   ->  gaussian summation points
    wtpts      ->  weight points for the gaussian points
_____________________________________________________________________________________________________________________________
Main Methods : -
----------------
=============================================================================================================================
update_attributes  ->   The method builds the material matrices based on the element number. 
                        The method has a special power to update the shape function, its derivatives and the Jacobian.
                        It also performs the gaussian integration.
#############################################################################################################################
'''
class material_model():
    def __init__(self):
        self.k  = None
        self.n  = None
        self.Vx = None
        self.Vy = None
#===========================================================================================================================#
#-----------------------------------------------#Method: -Update_Attributes#------------------------------------------------#
#===========================================================================================================================#
    def update_attributes(self,mesh,elem_no):
        self.k = np.zeros(4)
        self.n = np.zeros(4)
#========================================================Gauss loop=========================================================#
        gausspts = [-0.86113,-0.33998,0.33998,0.86113]
        wtpts    = [0.34785,0.65214,0.65214,0.34785]
        #gauss loop only one perform twice
        for r in range(4):
            for i in range(4):
                mesh.update_element(gausspts[i],gausspts[i],elem_no)
                A = (1/mesh.J[0,0])*mesh.dN[0,:].reshape(1,-1)
                B = (1/mesh.J[1,1])*mesh.dN[1,:].reshape(1,-1)
                Q = A.reshape(-1,1)
                P = B.reshape(-1,1)
                k = ip.k*(np.matmul(Q,A)\
                    + np.matmul(P,B))*mesh.det_J
                n = np.matmul(mesh.N,mesh.N.reshape(1,-1))*mesh.det_J
                self.k = self.k + wtpts[i]*k
                self.n = self.n + wtpts[i]*n
        self.n = ip.rho*ip.c*self.n
        mesh.update_element(0,0,elem_no)
        beta = ip.c*ip.theta_m/ip.Lf
        self.Vx = (beta/mesh.det_J)*(mesh.dN[0,:]*mesh.J[1,1]-mesh.dN[1,:]*mesh.J[0,1])
        self.Vy = (beta/mesh.det_J)*(-mesh.dN[0,:]*mesh.J[1,0]+mesh.dN[1,:]*mesh.J[0,0])
'''
=============================================================================================================================
--------------------------------------------------------End of Class---------------------------------------------------------
=============================================================================================================================
'''