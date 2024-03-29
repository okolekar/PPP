'''
Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Stefan Problem               Matriculation Number: 66808
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
This class generates the Initial Temperature Distribution for Enthalpy Approach
As soon as the init method is called, the method calculating the Temperature Distribution is called.
-----------------------------------------------------------------------------------------------------------------------------
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
    
'''
#############################################################################################################################
Class Stefan_Initial_Conditions: -
=============================================================================================================================
This class generates the Initial Temperature Distribution for Stefan Approach
-----------------------------------------------------------------------------------------------------------------------------
Attributes:-
------------
Tl                -> Stores the Initial Temperature distribution for the liquid domain
Ts                -> Stores the Initial Temperature distribution for the solid domain
-----------------------------------------------------------------------------------------------------------------------------
Methods: - 
=============================================================================================================================
Tdist_Stefan   -> Calculates the Temperature Distribution for the liquid and the solid domain                                                               
#############################################################################################################################
'''
class Stefan_Initial_Conditions():
    def __init__(self):
        self.Tl     = []
        self.Ts     = []
        self.s      = 0
        self.s_array = []
        self.ds_dt  = 0
        self.t      = 0
        self.hist_t = []
        self.formulate_s()

    def Tdist_Stefan(self,list_l,list_s,ds_dt,s):
        F = ip.ratioI-ip.lambf*ds_dt
#---------------------------------------------------------------------------------------------------------------------------#
                                    #Liquid domain initial temperature distribution
#---------------------------------------------------------------------------------------------------------------------------#
        for z in list_l:
            self.Tl.append(ip.ratioI*(s-z))
        self.Tl[-1] = 0
#---------------------------------------------------------------------------------------------------------------------------#
                                    #Solid domain initial temperature distribution
#---------------------------------------------------------------------------------------------------------------------------#
        for z in list_s:
            self.Ts.append(np.abs(ip.Ta)*np.exp(-(z-s)**2*F**2/(ip.Ta**2*np.pi))-F*(z-s)\
                      *m.erfc((z-s)*F/(np.abs(ip.Ta)*m.sqrt(np.pi)))+ip.Ta)
#---------------------------------------------------------------------------------------------------------------------------#
                                    #Application of the Dirichlet Boundary condition 
#---------------------------------------------------------------------------------------------------------------------------#
        self.Ts[0] = 0
        self.Ts[-1] = ip.Ta

    def formulate_s(self):
        t = 0.1293
        historyds_dt = [0]
        delt    = 0.00001
        R       = 0
        R1      = ip.ratioI*(2*np.sqrt(t/np.pi)*np.exp(-(self.s**2/(4*t)))-self.s*m.erfc(self.s/(2*t**0.5)))+ip.Ta 
        nrs     = 0
        for i in range(500):
            nrs = 0
            self.hist_t.append(t)
            while nrs < 20:
                R       =   R1
                dR      =  -ip.ratioI*m.erfc(self.s/(2*m.sqrt(t)))
                delt_s  =  -R/dR
                self.s  =   self.s + delt_s
                nrs     +=  1
                R1      =   ip.ratioI*(2*np.sqrt(t/np.pi)*np.exp(-(self.s**2/(4*t)))-self.s*m.erfc(self.s/(2*t**0.5)))+ip.Ta
                if abs(R1-R) < np.exp(-5):
                    if self.s<0:
                        self.s = 0
                    self.s_array.append(self.s)
                    break
            self.ds_dt = (self.s)/t
            historyds_dt.append(self.ds_dt)        
            t += delt
        self.ds_dt = historyds_dt[0]
        for i in range(len(self.s_array)):
            if self.s_array[i] > 0:
                self.s = self.s_array[i]
                self.ds_dt = historyds_dt
                self.t = self.hist_t[i]
                break

        self.t = t

'''
=============================================================================================================================
----------------------------------------------------****End of Class****-----------------------------------------------------
=============================================================================================================================
'''