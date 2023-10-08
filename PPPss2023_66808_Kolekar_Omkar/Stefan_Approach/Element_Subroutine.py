"""
Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Stefan Problem               Matriculation Number: 66808
Library: Elemental Subroutine
#############################################################################################################################
Importing the required libraries: -
=============================================================================================================================
#############################################################################################################################
"""
import Processed_Inputs as inputs
import numpy as np

class Elemental_Subroutine():
    '''
#############################################################################################################################
Class Elemental_Subroutine: -
=============================================================================================================================
Attributes: -
-------------
        self.Mlg = Material specific mass matrix for liquid domain
        self.Nlg = Material specific derivative matrix for liquid domain
        self.blg = Material specific force type vector for liquid domain
        self.Msg = Material specific mass matrix for solid domain
        self.Nsg = Material specific derivative matrix for solid domain
        self.bsg = Material specific force type vector for solid domain
-----------------------------------------------------------------------------------------------------------------------------
Description and Methods: -
--------------------------
get_global_Stefan_param: - Forms the global material parameters.
-----------------------------------------------------------------------------------------------------------------------------          
#############################################################################################################################
'''
    def __init__(self,s=None):
        self.Mlg = np.zeros([inputs.n,inputs.n])
        self.Nlg = np.zeros([inputs.n,inputs.n])
        self.blg = np.zeros(inputs.n).reshape(-1,1)
        self.Msg = np.zeros([inputs.m,inputs.m])
        self.Nsg = np.zeros([inputs.m,inputs.m])
        self.bsg = np.zeros(inputs.m).reshape(-1,1)

    def get_global_Stefan_param(self,meshl,meshs,ds_dt):
        self.Mlg = np.zeros([inputs.n,inputs.n])
        self.Nlg = np.zeros([inputs.n,inputs.n])
        self.blg = np.zeros(inputs.n).reshape(-1,1)
        self.Msg = np.zeros([inputs.m,inputs.m])
        self.Nsg = np.zeros([inputs.m,inputs.m])
        self.bsg = np.zeros(inputs.m).reshape(-1,1)
        np.fill_diagonal(self.Mlg,meshl.h)
        self.Mlg[0,0]   = meshl.h/2
        self.Mlg[-1,-1] = 1
        self.blg[0]     = inputs.ratioI
        dia = ds_dt/(3*inputs.n) + 2/(meshl.h)
        np.fill_diagonal(self.Nlg,dia)
        for i in range(1,inputs.n-1):
            self.Nlg[i-1,i] = -(1/6)*((3*i-2)/(inputs.n))*ds_dt - 1/meshl.h
            self.Nlg[i+1,i] = (1/6)*((3*i+2)/(inputs.n))*ds_dt - 1/meshl.h
        self.Mlg[-1,-1] = 1
        self.Nlg[0,0] = (1/(6*inputs.n))*ds_dt+1/meshl.h
        self.Nlg[-1,-1] = 1
#---------------------------------------------------------------------------------------------------------------------------#        
        np.fill_diagonal(self.Msg,meshs.h)
        self.Msg[-1,-1] = 1
        self.Msg[0,0] = 1
        self.bsg[-2] = inputs.Ta*((1/(3*inputs.m))*ds_dt+1/(meshs.h))
        dia = -(1/(3*inputs.m))*ds_dt + 2/(meshs.h)
        np.fill_diagonal(self.Nsg,dia)
        for i in range(1,inputs.m-1):
            self.Nsg[i-1,i] = -(1/6)*((3*inputs.m-3*i+2)/(inputs.m))*ds_dt - 1/meshs.h
            self.Nsg[i+1,i] = (1/6)*((3*inputs.m-3*i-2)/(inputs.m))*ds_dt - 1/meshs.h
        self.Nsg[0,0] = 1
        self.Nsg[-1,-1] = 1