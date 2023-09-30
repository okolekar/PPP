'''Topic: PPP
Library: - Heat_Transfer_verification 
#############################################################################################################################
Importing the required standard libraries 
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs -> Script where all the inputs are defined                                                                  
#############################################################################################################################
'''
import numpy as np
import Processed_Inputs as ip
'''
#############################################################################################################################
Class Verification: -
=============================================================================================================================
Generates the Analytical solution for the Heat Transfer.
-----------------------------------------------------------------------------------------------------------------------------
Attributes:-
------------
T_analytical    ->  Stores the Analytical Temperature distribution
-----------------------------------------------------------------------------------------------------------------------------
The Boundary conditions assumed are as follows: -
=============================================================================================================================
@ x = L -> T = Ta = 0. In dimentionaless form, this means that the material is completely melted and  
                       the temperature @ x = L is equal to the Melting Point.
@ x = 0 -> T > 0    
-----------------------------------------------------------------------------------------------------------------------------
*For this dimensionless form the Analytical Solution has assumed that Total length L = 1.                                                                                                          
#############################################################################################################################
'''
class Verification():
    def __init__(self,Tinitial, n,nlist,alpha):
        self.T_analytical = np.zeros(ip.n).reshape(n,1)
        for i in range(n):
            x = nlist[i]
            self.T_analytical[i] = (Tinitial[i])*(np.cos(np.pi*x))*np.exp((-np.pi**2*alpha*ip.delt))
