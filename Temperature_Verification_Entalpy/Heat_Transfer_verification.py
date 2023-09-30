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
The Boundary conditions assumed are as follows : -
=============================================================================================================================
Boundary Condition type 1: -
----------------------------
@ x = 0 and x = L -> T = Ta = 0. In dimentionaless form, this means that the material is completely melted and  
                                 the temperature @ x = L is equal to the Melting Point.
-----------------------------------------------------------------------------------------------------------------------------
Boundary Condition type 2: -
----------------------------
@ x = L -> T = Ta = 0. In dimentionaless form, this means that the material is completely melted and  
                       the temperature @ x = L is equal to the Melting Point.
@ x = 0 -> dT/dx = 0       Material in liquid state. 
-----------------------------------------------------------------------------------------------------------------------------
*For this dimensionless form the Analytical Solution has assumed that Total length L = 1.                                                                                                          
#############################################################################################################################
'''
class Verification():
    def __init__(self,n,nlist,t):
        self.T_analytical = np.zeros(ip.n).reshape(n,1)
        for i in range(n):
            x = nlist[i]
            self.T_analytical[i] = (3/2)*np.sin(np.pi*x)*np.exp(-np.pi**2*t) - (1/2)*np.sin(3*np.pi*x)*np.exp(-9*np.pi**2*t)
            #Both ends @ T = 0
            """P = 0
            for q in range(1,5):
                P = P + (((-1)**(q-1)/(2*q-1)))*np.cos((2*q-1)*np.pi*x/(2*ip.length))*np.exp(-((2*q-1)**2)*(np.pi**2)*t/(4*ip.length**2))
                #x = L end @ T = 0
            self.T_analytical[i] = (1*4/np.pi)*P"""