'''Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
Testing Amorphous material slope
#############################################################################################################################
Importing the required standard libraries                                                                
#############################################################################################################################
'''
import numpy as np
from Slope_MatrixA import temp_der
import Processed_InputsA as ip
#--------------------------------------------For the Amorphous Settings-----------------------------------------------------#
def test_A_temp_derivative():
    Hk = np.array([-1,ip.Hliq+0.5]).reshape(-1,1)
    dT_dH = temp_der(2,Hk,3).dT_dH
    answer = np.array([[1/ip.D,0],
                       [0,1/ip.D]])
    for i in range(2):
        assert dT_dH[i,0] == answer[i,0],"Temperature Derivative wrt Enthlapy failed for '[-1,ip.Hliq+0.5]' enthalpy inputs"
        assert dT_dH[i,1] == answer[i,1],"Temperature Derivative wrt Enthlapy failed for '[-1,ip.Hliq+0.5]' enthalpy inputs"
#-------------------------------------------------------Case2---------------------------------------------------------------#
def test_A_temp_derivative2():
    Hk = np.array([0,ip.Hliq]).reshape(-1,1)
    dT_dH = temp_der(2,Hk,3).dT_dH
    answer = np.array([[1/ip.D,0],
                       [0,ip.Tliq/(ip.D*ip.Tliq + ip.D*ip.lambf)]])
    for i in range(2):
        assert dT_dH[i,0] == answer[i,0],"Temperature Derivative wrt Enthlapy failed for '[0,ip.Hliq]' enthalpy inputs"
        assert dT_dH[i,1] == answer[i,1],"Temperature Derivative wrt Enthlapy failed for '[0,ip.Hliq]' enthalpy inputs"
#-------------------------------------------------------Case3---------------------------------------------------------------#
def test_A_temp_derivative3():
    Hk = np.array([0.000001,ip.Hliq+0.000001]).reshape(-1,1)
    dT_dH = temp_der(2,Hk,3).dT_dH
    answer = np.array([[ip.Tliq/(ip.D*ip.Tliq + ip.D*ip.lambf),0],
                       [0,1/ip.D]])
    for i in range(2):
        assert dT_dH[i,0] == answer[i,0],"Temperature Derivative wrt Enthlapy failed for '[0.000001,ip.Hliq+0.000001]' enthalpy inputs"
        assert dT_dH[i,1] == answer[i,1],"Temperature Derivative wrt Enthlapy failed for '[0.000001,ip.Hliq+0.000001]' enthalpy inputs"
#-------------------------------------------------------Case4---------------------------------------------------------------#
def test_A_temp_derivative4():
    Hk = np.array([0,ip.Hliq-0.000001]).reshape(-1,1)
    dT_dH = temp_der(2,Hk,3).dT_dH
    answer = np.array([[1/ip.D,0],
                       [0,ip.Tliq/(ip.D*ip.Tliq + ip.D*ip.lambf)]])
    for i in range(2):
        assert dT_dH[i,0] == answer[i,0],"Temperature Derivative wrt Enthlapy failed for '[0,ip.Hliq-0.000001]' enthalpy inputs"
        assert dT_dH[i,1] == answer[i,1],"Temperature Derivative wrt Enthlapy failed for '[0,ip.Hliq-0.000001]' enthalpy inputs"
#---------------------------------------------------------------------------------------------------------------------------#
