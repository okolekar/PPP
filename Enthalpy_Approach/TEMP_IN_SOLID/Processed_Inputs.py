# Processed_Inputs
'''
#############################################################################################################################
Importing the standard library: -
-----------------------------------------------------------------------------------------------------------------------------
Crystalline_inputs  -> Script where all the inputs for crystalline material are defined
Amorphous_inputs    -> Script where all the inputs for amorphous material are defined
#############################################################################################################################
'''
import numpy as np
'''
#############################################################################################################################
get_input function: -
-----------------------------------------------------------------------------------------------------------------------------
Takes care of the scenario when no inputs are passed and runs the default test
#############################################################################################################################
'''
#---------------------------------------------------------------------------------------------------------------------------#
mat_type  = 2
test_case = 0

#===========================================================================================================================#
                         #Check if the material type is entered correctly
#===========================================================================================================================#
if mat_type != 1 and mat_type != 2:
    raise ValueError("Incorrect material type. Please enter 1 or 2 for the respective element.")

if mat_type == 1:
     import Crystalline_inputs as ip

if mat_type == 2:
     import Amorphous_inputs as ip
'''
#############################################################################################################################
Script: -
=============================================================================================================================
Variables: -
mat_type  ->     Indicates the type of material
test_case ->     Check whether to run the test case
ratioI    ->     Dimensionless input Laser Ratio
alpha     ->     Conductivity
delt      ->     Time step 
n         ->     Number of nodes
Ta        ->     Dimensionless Ambient Temperature
D         ->     Dimensionless Material Constant
length    ->     Dimensionless Length
lambf     ->     Dimensionless Material Constant
t_start   ->     Dimensionless start time of the simulation
tm        ->     Dimensionless time when the material starts to melt
t_end     ->     End time for the simulation
-----------------------------------------------------------------------------------------------------------------------------
Tests Performed: -
=============================================================================================================================
1) Checks that latent heat of fusion, vapour point, thermal conductivity, density,specific heat capacity 
   and time step values are greater than zero.
2) Check that at least there are two nodes present.
3) Checks if all the dimensionless material specific parameters are positive.                                                                     
#############################################################################################################################
'''
if ip.Lf < 0 or ip.vapourpt<0 or ip.k<0 or ip.delt<0 or ip.n < 2 or ip.rho < 0 or ip.c<0 or ip.I < 0 or ip.Iref <= 0\
       or ip.t_end < 0:
        raise ValueError("Incorrect material input arguments")
#===========================================================================================================================#
                                   #Defining the dimensionless Parameters
#===========================================================================================================================#
ratioI = ip.I/ip.Iref
alpha = ip.k/(ip.rho*ip.c)
delt = ip.delt
n = ip.n
t_start = delt
'''
#############################################################################################################################
Crystaline Material Specific Dimensionless Parameters: -
#############################################################################################################################
'''
if mat_type == 1:
    Ta = (293-ip.meltpt)/(ip.vapourpt-ip.meltpt)             
    D = ip.rho*ip.c*(ip.vapourpt-ip.meltpt)/((ip.rho*ip.c*ip.vapourpt+ip.Lf*ip.rho)-(ip.rho*ip.c*ip.meltpt))
    length = (ip.depth*ip.Iref)/(ip.k*(ip.vapourpt-ip.meltpt))
    lambf = ip.Lf/(ip.c*(ip.vapourpt-ip.meltpt))
'''
#############################################################################################################################
Amorphous Material Specific Dimensionless Parameters: -
=============================================================================================================================
Variables: -
Tliq      ->     Dimensionless Liquidus Temperature
Hliq      ->     Dimensionless Liquidus Enthalpy
-----------------------------------------------------------------------------------------------------------------------------
#############################################################################################################################
'''
if mat_type == 2:
     Tliq = (ip.liquidouspt-ip.soliduspt)/(ip.vapourpt-ip.soliduspt)                
     D = ip.rho*ip.c*(ip.vapourpt-ip.soliduspt)/((ip.rho*ip.c*ip.vapourpt+ip.Lf*ip.rho)-(ip.rho*ip.c*ip.soliduspt))
     Ta = (293-ip.soliduspt)/(ip.vapourpt-ip.soliduspt)                             
     lambf = ip.Lf/(ip.c*(ip.vapourpt-ip.soliduspt)) 	  	                        
     Hliq = D*Tliq + D*lambf
     alpha = ip.k/(ip.rho*ip.c)
     length = 100#ip.Iref/(ip.k*(ip.vapourpt-ip.soliduspt))
#===========================================================================================================================#
                         #Testing if the dimensionless parameters are correctly calculated
#===========================================================================================================================# 	  	     
if D < 0 or lambf < 0 or alpha < 0 or ratioI<0 or length < 0 :
     raise ValueError("Incorrect material input arguments")
'''
#===========================================================================================================================#
Checking the test case scenario.
For test case, setting the length to unity and Ta = 0 as per the assumption
#===========================================================================================================================#
'''
if test_case == 1:
     Ta = 0
     length = 1
else:
     tm = (ip.Iref**2*Ta**2*np.pi)/(4*ip.I**2)
     t_end = ip.t_end
