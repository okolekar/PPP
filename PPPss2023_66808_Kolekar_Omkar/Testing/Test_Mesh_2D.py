'''Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
Testing
#############################################################################################################################
Importing the required standard libraries                                                                
#############################################################################################################################
'''
import numpy as np
from Gridify import  Mesh
#============================================Form object of Mesh Class======================================================#
M = Mesh(5,5,1,1)
#=============================================test_partition_of_unity=======================================================#
def test_partition_of_unity():
    gpts = [-0.577350269189626, 0.577350269189626]
    for zeta in gpts:
        for eta in gpts:
            M.update_element(zeta,eta,3)
            assert np.sum(M.N) == 1, "Partition of unity failed sum of the shape function not = 1"
            assert np.sum(M.dN) == 0, "Partition of unity failed sum of the derivative of shape function not = 0"
#===================================================test_Jacobian===========================================================#
def test_Jacobian():
    gpts = [-0.577350269189626, 0.577350269189626]
    sum = 0
    for zeta in gpts:
        for eta in gpts:
            M.update_element(zeta,eta,3)
            assert M.det_J > 0, "test_Jacobian failed, determinant of the Jacobian = 0"
            sum = sum + M.det_J
    x = M.xy_list[M.node_list[3][2]][0]-M.xy_list[M.node_list[3][0]][0]
    y = M.xy_list[M.node_list[3][2]][1]-M.xy_list[M.node_list[3][0]][1]
    assert sum == x*y, "test_Jacobian failed, sum of determinant of the Jacobian after integration not equal to area of element"
#==================================================test_derivative==========================================================#
def test_derivative():
    M.update_element(1,1,3)
    del_x = M.xy_list[M.node_list[3][2]][0]-M.xy_list[M.node_list[3][0]][0]
    del_y = M.xy_list[M.node_list[3][2]][1]-M.xy_list[M.node_list[3][0]][1]
    T = [1,2,3,4]
    expect1 = (T[2]-T[3])/del_x
    expect2 = (T[2]-T[1])/del_y
    answer = M.dN@T
    assert answer[0] == expect1,"The test_derivative function failed to calculate the temperature derivatives"
    assert answer[1] == expect2,"The test_derivative function failed to calculate the temperature derivatives"
#===================================================End of Tests============================================================#