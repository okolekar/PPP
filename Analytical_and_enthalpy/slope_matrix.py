'''Topic: PPP
Library: - Slope Matrix
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
Check arguments Decorator: -
-----------------------------------------------------------------------------------------------------------------------------
Checks if exactly 2 argument are passed to this Class                                                                        
#############################################################################################################################
'''
def check_arguments(func):
    def wrapper(self, *args,**kwargs):
        if len(args) == 3:
            return func(self,*args,**kwargs)
        else:
            raise ValueError(f"{len(args)} are incorrect number of input arguments")
    return wrapper
'''
#############################################################################################################################
Class temp_der: -
=============================================================================================================================
This class has only one job to calculate dT/dH matrix
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    dT_dH       ->  Slope matrix
-----------------------------------------------------------------------------------------------------------------------------
Methods: -
=============================================================================================================================
    verify      -> Checks if the slope matrix was created correctly by checking the off diagonal elements to be zero 
    Amorph_slop -> Calculates the slope matrix for the amorphous materials
    Cryst_slop  -> Calculates the slope matrix for the crystalline materials                                                 
#############################################################################################################################
'''
class temp_der():

    @check_arguments
    def __init__(self,*args):
        self.dT_dH = np.zeros([1,1])
        runcount = args[2]
        if args[0] == 1:
            if runcount<2 and runcount > 0:
                print('Crystalline material detected in slope_matrix module')
            Hk_1 = args[1]
            self.Cryst_slop(Hk_1)

        if args[0] == 2:
            if runcount<2 and runcount > 0:
                print('Amorphous material detected in slope_matrix module')
            Hk_1 = args[1]
            self.Amorph_slop(Hk_1)
#===========================================================================================================================#
                            #Verify: - Checks if the off diagonal elements are zero
#===========================================================================================================================# 
    def verify(self,Hk_1,c):
        if self.dT_dH[0][1] != 0 or self.dT_dH[1][0] != 0:
            raise ValueError("The 2*2 dT/dH Matrix is not constructed correctly")
        if c == 1:
            if Hk_1[0] >= 0 and Hk_1[0] <= ip.lambf*ip.D and self.dT_dH[0][0] != 0:
                raise ValueError(f"Problem in dT/dH the Enthalpy at node 1 {Hk_1[0]} and corresponding dT/dH {self.dT_dH[0][0]}")
            if Hk_1[1] >= 0 and Hk_1[1] <= ip.lambf*ip.D and self.dT_dH[1][1] != 0:
                raise ValueError(f"Problem in dT/dH the Enthalpy at node 1 {Hk_1[0]} and corresponding dT/dH {self.dT_dH[0][0]}")
            '''if Hk_1[1]<0 or Hk_1[1]> ip.D*ip.lambf and self.dT_dH[1][1] != 1/ip.D:
                raise ValueError(f"The 2*2 dT/dH Matrix is not constructed correctly the matrix is {self.dT_dH} nad it was expected to be {1/ip.D}")
            if Hk_1[0]<0 or Hk_1[0]> ip.D*ip.lambf and self.dT_dH[0][0] != 1/ip.D:
                raise ValueError("The 2*2 dT/dH Matrix is not constructed correctly the matrix is {self.dT_dH}")'''
#===========================================================================================================================#
                        #Amorph_slop: - Defines the dT_dH slope matrix for Amorphous Material
#===========================================================================================================================#
    def Amorph_slop(self, Hk_1):
        try:
            self.dT_dH = [[1/ip.D if Hk_1[0] <= 0 else ip.Tliq/(ip.D*ip.Tliq + ip.D*ip.lambf) if np.logical_and(0 <= Hk_1[0], Hk_1[0] <= ip.Hliq) else 1/ip.D, 0],
                          [0,1/ip.D if Hk_1[1] <= 0 else ip.Tliq/(ip.D*ip.Tliq + ip.D*ip.lambf) if np.logical_and(0 <= Hk_1[1], Hk_1[1] <= ip.Hliq) else 1/ip.D]]
            self.verify(Hk_1,2)
        except ValueError as e:
            print(f"An error in dT_dH: {str(e)} following were the inputs")
            print(f"Hk+1 = {Hk_1}")
            raise ValueError("Check the variables passed")   
#===========================================================================================================================#
                        #Cryst_slop: - Defines the dT_dH slope matrix for Crystaline Material
#===========================================================================================================================#
    def Cryst_slop(self, Hk_1):
        try:
            self.dT_dH = np.array([[1/ip.D if Hk_1[0] < 0 else 0 if np.logical_and(Hk_1[0] >= 0, Hk_1[0] <= ip.D*ip.lambf)else 1/ip.D, 0],
                                   [0,1/ip.D if Hk_1[1] < 0 else 0 if np.logical_and(Hk_1[1] >= 0, Hk_1[1] <= ip.D*ip.lambf) else 1/ip.D]])
            self.verify(Hk_1,1)
            
        except ValueError as e:
            print(f"An error in dT_dH: {str(e)} following were the inputs")
            print(f"Hk+1 = {Hk_1}")
            raise ValueError("Check the variables passed")