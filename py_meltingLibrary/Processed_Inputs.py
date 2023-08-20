'''Topic: PPP
Library: - Processed_Inputs
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
def get_input(prompt, default_value):
    user_input = input(f"{prompt} else default case {default_value}: - ")
    return user_input if user_input else default_value
#---------------------------------------------------------------------------------------------------------------------------#
try:
     mat_type  = get_input("Enter 1 for crystalline material and 2 for amorphous material,",1)
     test_case = get_input("Enter 1 to run the Heat Transfer Verification case or 0 to disable,",0)
     mat_type = int(mat_type)
     test_case = int(test_case)
except ValueError as e:
     raise e
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
test_case ->     Check whether to run the Verification case
ratioI    ->     Dimensionless input Laser Ratio
alpha     ->     Conductivity
delt      ->     Time step 
n         ->     Number of nodes (in case of Stefan approach its number of nodes in the liquid domain)
m         ->     Number of nodes in the solid domain for Stefan approach
t_start   ->     Dimensionless start time of the simulation
Ta        ->     Dimensionless Ambient Temperature
D         ->     Dimensionless Material Constant
length    ->     Dimensionless Length
lambf     ->     Dimensionless Material Constant
s         ->     Initial position of the solid-liquid interface
ds_dt     ->     Initial solid-liquid interface  speed of travel
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
       or ip.t_end < 0 or ip.m < 2:
        raise ValueError("Incorrect material input arguments")
#===========================================================================================================================#
                                   #Defining the dimensionless Parameters
#===========================================================================================================================#
ratioI    = ip.I/ip.Iref
alpha     = ip.k/(ip.rho*ip.c)
delt      = ip.delt
n         = ip.n
m         = ip.m
t_start   = delt
'''
#############################################################################################################################
Crystaline Material Specific Dimensionless Parameters: -
#############################################################################################################################
'''
if mat_type == 1:
    Ta         = (293-ip.meltpt)/(ip.vapourpt-ip.meltpt)             
    D          = ip.rho*ip.c*(ip.vapourpt-ip.meltpt)/((ip.rho*ip.c*ip.vapourpt+ip.Lf*ip.rho)-(ip.rho*ip.c*ip.meltpt))
    length     = (ip.depth*ip.Iref)/(ip.k*(ip.vapourpt-ip.meltpt))
    lambf      = ip.Lf/(ip.c*(ip.vapourpt-ip.meltpt))
    s          = 0.5
    ds_dt      = ip.ds_dt
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
     Tliq      = (ip.liquidouspt-ip.soliduspt)/(ip.vapourpt-ip.soliduspt)                
     D         = ip.rho*ip.c*(ip.vapourpt-ip.soliduspt)/((ip.rho*ip.c*ip.vapourpt+ip.Lf*ip.rho)-(ip.rho*ip.c*ip.soliduspt))
     Ta        = (293-ip.soliduspt)/(ip.vapourpt-ip.soliduspt)                             
     lambf     = ip.Lf/(ip.c*(ip.vapourpt-ip.soliduspt)) 	  	                        
     Hliq      = D*Tliq + D*lambf
     alpha     = ip.k/(ip.rho*ip.c)
     length    = ip.Iref/(ip.k*(ip.vapourpt-ip.soliduspt))
     s         = 0.5
     ds_dt     = ip.ds_dt
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