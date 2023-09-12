'''Topic: PPP
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
import matplotlib.pyplot as plt
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
Tdist_Enthalpy   -> Calculates the Temperature Distribution                                                               '''
#############################################################################################################################
class Initial_Temp():

    @check_arguments
    def __init__(self,list):
        self.T = np.zeros(ip.n)
        self.s = 0
        self.s_array = [0]
        self.ds_dt = 0
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
        F = ip.ratioI
        try:
            for i in range(len(self.T)):
                if i == len(self.T)-1:
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


    def formulate_s(self):
        hist_t = [0]
        t = 0.01293
        historyds_dt = [0]
        delt    = ip.delt
        R1      = ip.ratioI*(2*np.sqrt(t/np.pi)*np.exp(-(self.s**2/(4*t)))-self.s*m.erfc(self.s/(2*t**0.5)))+ip.Ta 
        for i in range(ip.t_end):
            nrs = 0
            hist_t.append(t)
            while nrs < 20:
                R       =   R1
                dR      =  -ip.ratioI*m.erfc(self.s/(2*m.sqrt(t)))
                delt_s  =  -R/dR
                self.s  =   self.s + delt_s
                nrs     +=  1
                R1      =   ip.ratioI*(2*np.sqrt(t/np.pi)*np.exp(-(self.s**2/(4*t)))-self.s*m.erfc(self.s/(2*t**0.5)))+ip.Ta
                if abs(R1-R) < np.exp(-5):
                    self.s_array.append(self.s)
                    break
            self.ds_dt = (self.s)/t
            historyds_dt.append(self.ds_dt)        
            t += delt
        self.ds_dt = (self.s)/t