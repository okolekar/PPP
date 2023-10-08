'''Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
Testing the Crystalline material slope
#############################################################################################################################
Importing the required standard libraries                                                                
#############################################################################################################################
'''
import numpy as np
from slope_matrix import temp_der
import Processed_Inputs as ip
#--------------------------------------------For the Crystalline Settings---------------------------------------------------#
def test_temp_derivative():
    Hk = np.array([-1,ip.D*ip.lambf]).reshape(-1,1)
    dT_dH = temp_der(1,Hk,3).dT_dH
    answer = np.array([[1/ip.D,0],
                       [0,0]])
    for i in range(2):
        assert dT_dH[i,0] == answer[i,0], "Temperature Derivative wrt Enthlapy failed for '[-1,ip.D*lambf]' enthalpy inputs"
        assert dT_dH[i,1] == answer[i,1], "Temperature Derivative wrt Enthlapy failed for '[-1,ip.D*lambf]' enthalpy inputs"  
#-------------------------------------------------------Case2---------------------------------------------------------------#
def test_temp_derivative2():
    Hk = np.array([0,ip.D*ip.lambf+0.5]).reshape(-1,1)
    dT_dH = temp_der(1,Hk,3).dT_dH
    answer = np.array([[0,0],
                       [0,1/ip.D]])
    for i in range(2):
        assert dT_dH[i,0] == answer[i,0],"Temperature Derivative wrt Enthlapy failed for '[0,ip.Dip.lambf+0.5]' enthalpy inputs"
        assert dT_dH[i,1] == answer[i,1],"Temperature Derivative wrt Enthlapy failed for '[0,ip.Dip.lambf+0.5]' enthalpy inputs"
#-------------------------------------------------------Case3---------------------------------------------------------------#
def test_temp_derivative3():
    Hk = np.array([ip.lambf,-1]).reshape(-1,1)
    dT_dH = temp_der(1,Hk,3).dT_dH
    answer = np.array([[1/ip.D,0],
                       [0,1/ip.D]])
    for i in range(2):
        assert dT_dH[i,0] == answer[i,0],"Temperature Derivative wrt Enthlapy failed for '[ip.lambf,-1]' enthalpy inputs"
        assert dT_dH[i,1] == answer[i,1],"Temperature Derivative wrt Enthlapy failed for '[ip.lambf,-1]' enthalpy inputs"
#-------------------------------------------------------Case4---------------------------------------------------------------#
def test_temp_derivative4():
    Hk = np.array([0,ip.lambf*ip.D]).reshape(-1,1)
    dT_dH = temp_der(1,Hk,3).dT_dH
    answer = np.array([[0,0],
                       [0,0]])
    for i in range(2):
        assert dT_dH[i,0] == answer[i,0],"Temperature Derivative wrt Enthlapy failed for '[0,ip.lambf*ip.D]' enthalpy inputs"
        assert dT_dH[i,1] == answer[i,1],"Temperature Derivative wrt Enthlapy failed for '[0,ip.lambf*ip.D]' enthalpy inputs"
#---------------------------------------------------End of Testing----------------------------------------------------------#