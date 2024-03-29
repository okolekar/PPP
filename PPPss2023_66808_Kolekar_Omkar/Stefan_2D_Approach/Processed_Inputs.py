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
import Crystalline_inputs as ip
'''
#############################################################################################################################
get_input function: -
-----------------------------------------------------------------------------------------------------------------------------
Takes care of the scenario when no inputs are passed and runs the default test
#############################################################################################################################
'''
T = ip.theta
alpha = ip.k/(ip.rho*ip.c)
x = ip.X/ip.L
y = ip.Y/ip.L
c = ip.c
theta_m
#---------------------------------------------------------------------------------------------------------------------------#
#===========================================================================================================================#
                         #Check if the material type is entered correctly
#===========================================================================================================================#

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
#===========================================================================================================================#
                                   #Defining the dimensionless Parameters
#===========================================================================================================================#
'''
#############################################################################################################################
Crystaline Material Specific Dimensionless Parameters: -
#############################################################################################################################
'''
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

#===========================================================================================================================#
                         #Testing if the dimensionless parameters are correctly calculated
#===========================================================================================================================# 	  	     

'''
#===========================================================================================================================#
Checking the test case scenario.
For test case, setting the length to unity and Ta = 0 as per the assumption
#===========================================================================================================================#
'''