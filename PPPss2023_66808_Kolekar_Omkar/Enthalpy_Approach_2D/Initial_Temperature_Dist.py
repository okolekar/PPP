'''
Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
Library Initial_Temperature_Dist
#############################################################################################################################
Importing the required standard libraries 
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs -> Script where all the inputs are defined
#############################################################################################################################
'''
import numpy as np
import math as m
import Processed_Inputs as ip 
'''
#############################################################################################################################
Class Initial Temp: -
=============================================================================================================================
As soon as the init method is called the method calculating the Temperature Distribution is called.
Attributes:-
------------
T                -> Stores the Initial Temperature distribution
-----------------------------------------------------------------------------------------------------------------------------
Tests: -
=============================================================================================================================
check_T          -> Checks if the first entry in the Temperature is greater than the ambient Temperature.
-----------------------------------------------------------------------------------------------------------------------------
Methods: - 
=============================================================================================================================
Tdist_Enthalpy   -> Calculates the Temperature Distribution                                                               '''
#############################################################################################################################
class Initial_Temp():

    def __init__(self,list,Mesh):
        self.T = np.zeros(ip.n*ip.m)
        self.hist_t = None
        self.Tdist_Enthalpy(list,Mesh)

    def check_T(self,Ta):
        if self.T[0]<Ta:
            print(self.T)
            print("Temperature distribution might be invalid.")
            print("As the Temperature at node 1 is lower than the Ta.")
    '''
#############################################################################################################################
Method Tdist_Enthalpy()
-----------------------------------------------------------------------------------------------------------------------------
This method has the try and except as for the Test Heat Transfer Test Case, to capture any valueErrors that may arrise 
during the simualtion                                                                                                       
#############################################################################################################################
'''
    def Tdist_Enthalpy(self,list,Mesh):
        F = ip.ratioI
        try:
            for i in range(len(self.T)):
                if i in Mesh.boundary_nodes["Bottom_Nodes"]:
                    self.T[i] = ip.Ta #This takes care of the boundary condition
                elif i in Mesh.boundary_nodes["Right_Nodes"]:
                    self.T[i] = ip.Ta #This takes care of the boundary condition
                else:
                    self.T[i] = (F*(2*np.sqrt(ip.t_start/np.pi)
                                    *np.exp(-list[i]**2/(4*ip.t_start))
                                    -list[i]*m.erfc(list[i]/(2*np.sqrt(ip.t_start))))
                                 + ip.Ta)
        except ZeroDivisionError:
            for i in range(len(self.T)):
                if i == len(self.T)-1:
                    self.T[i] = ip.Ta #This takes care of the boundary condition
                else:
                    print("Zero Division error in Initial Temperature. Changing t_start to delt")
                    self.T[i] = (F*(2*np.sqrt(ip.t_start/np.pi)
                                    *np.exp(-list[i]**2/(4*ip.delt))
                                    -list[i]*m.erfc(list[i]/(2*np.sqrt(ip.delt))))
                                + ip.Ta)        
        self.T = self.T[:,np.newaxis]
        self.check_T(ip.Ta)
