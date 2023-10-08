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
import numpy as np 
'''
#############################################################################################################################
Check arguments Decorator: -
-----------------------------------------------------------------------------------------------------------------------------
Checks if atleast 1 argument is passed to this Class                                                                         
#############################################################################################################################
'''
def check_arguments(func):
    def wrapper(self, *args,**kwargs):
        if len(args) != 1:
            raise ValueError("Incorrect number of input arguments to Initial_Temp.")
        return func(self,*args,**kwargs)
    return wrapper
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
Tdist_Enthalpy   -> Calculates the Temperature Distribution                                                               
#############################################################################################################################
'''
class Initial_Temp():

    @check_arguments
    def __init__(self,list):
        self.T = np.zeros(ip.n)
        self.Tdist_Enthalpy(list)

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
    def Tdist_Enthalpy(self,list):
        for i in range(len(self.T)):
            x = list[i]
            self.T[i] = -0.5*np.sin(3*np.pi*x)+(3/2)*np.sin(np.pi*x)# for BCs Type 1 (sine)
            #self.T[i] = 1 # For BCs type 2
        self.T[-1] = 0
        self.T[0] = 0 ## For BCs type 1 sine    
        self.T = self.T[:,np.newaxis]
        self.check_T(ip.Ta)