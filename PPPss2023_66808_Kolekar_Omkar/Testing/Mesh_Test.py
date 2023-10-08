'''Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
Testing 1D Mesh
#############################################################################################################################
Importing the required standard libraries                                                                
#############################################################################################################################
'''
from Mesh import Mesh as mesh
import numpy as np
#---------------------------------------------------------------------------------------------------------------------------#
def test_dphi():
    M = mesh()
    M.update_Jacobi(2)
    M.update_dphi()
#---------------------------------------------Test Shape Function Derivative------------------------------------------------#
    T = np.array([2,1]).reshape(-1,1)
    k = round(np.matmul(M.dphi.reshape(1,-1),T).item(),5)
    assert k == -100, "The Shape function derivative failed for T = [2,1]"
    T = np.array([1,3]).reshape(-1,1)
    k = round(np.matmul(M.dphi.reshape(1,-1),T).item(),5)
    assert k == 200, "The Shape function derivative failed for T = [1,3]"
    T = np.array([4,4]).reshape(-1,1)
    k = round(np.matmul(M.dphi.reshape(1,-1),T).item(),5)
    assert k == 0, "The Shape function derivative failed for T = [4,4]"
#------------------------------------------------Test Partition of Unity----------------------------------------------------#    
    M.update_phi(0.5)
    assert np.sum(M.phi) == 1,"Partition of unity test failed"
#-----------------------------------------------------Test Jacobian---------------------------------------------------------#
    expected_J = 0.5*(M.nl[3]-M.nl[2])
    assert M.J == expected_J, "Test Jacobian failed"
#-----------------------------------------------------End of Testing--------------------------------------------------------#