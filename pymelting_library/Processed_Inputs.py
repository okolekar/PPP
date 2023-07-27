# Processed_Inputs
'''
#############################################################################################################################
Importing the standard library: -
-----------------------------------------------------------------------------------------------------------------------------
Crystalline_inputs  -> Script where all the inputs for crystalline material are defined
Amorphous_inputs    -> Script where all the inputs for amorphous material are defined
Material_Subroutine -> All material specific equations and parameters are processed here'''
#############################################################################################################################
import numpy as np
try:
     mat_type = int(input("Enter 1 for crystalline material and 2 for amorphous material: -  "))
except ValueError as e:
     raise ValueError("Invalid input detected. Please enter 1 or 2 for the respective element.")
if mat_type != 1 and mat_type != 2:
    raise ValueError("Incorrect material type. Please enter 1 or 2 for the respective element.")

if mat_type == 1:
     import inputs_Al_Material as ip

if mat_type == 2:
     import inputs as ip

'''
#############################################################################################################################
Script: -
=============================================================================================================================
Variables: -
mat_type  ->     Indicates the type of material
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
-----------------------------------------------------------------------------------------------------------------------------
Tests Performed: -
=============================================================================================================================
1) Checks that latent heat of fusion, vapour point, thermal conductivity, density,specific heat capacity 
   and time step values are greater than zero.
2) Check that at least there are two nodes present.
3) Checks if all the dimensionless material specific parameters are positive.                                                                     
#############################################################################################################################
'''
if ip.Lf < 0 or ip.vapourpt<0 or ip.k<0 or ip.delt<0 or ip.n < 2 or ip.rho < 0 or ip.c<0:
        raise ValueError("Incorrect material input arguments")
ratioI = ip.I/ip.Iref
alpha = ip.k/(ip.rho*ip.c)
delt = ip.delt
n = ip.n
if mat_type == 1:
    Ta = (293-ip.meltpt)/(ip.vapourpt-ip.meltpt)             # ambient temperature in Dimensionless.
    D = ip.rho*ip.c*(ip.vapourpt-ip.meltpt)/((ip.rho*ip.c*ip.vapourpt+ip.Lf*ip.rho)-(ip.rho*ip.c*ip.meltpt))
    length = (ip.depth*ip.Iref)/(ip.k*(ip.vapourpt-ip.meltpt))
    lambf = ip.Lf/(ip.c*(ip.vapourpt-ip.meltpt)) 	  	     # The lambda constant
    if D < 0 or lambf < 0 or alpha < 0 or ratioI<0 or ip.meltpt<0 :
        raise ValueError("Incorrect material input arguments")
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
     Tliq = (ip.liquidouspt-ip.soliduspt)/(ip.vapourpt-ip.soliduspt)                #Dimensionless 	  	         #Liq. Temp.
     D = ip.rho*ip.c*(ip.vapourpt-ip.soliduspt)/((ip.rho*ip.c*ip.vapourpt+ip.Lf*ip.rho)-(ip.rho*ip.c*ip.soliduspt))
     Ta = (293-ip.soliduspt)/(ip.vapourpt-ip.soliduspt)                             # ambient temperature in Dimensionless.
     lambf = ip.Lf/(ip.c*(ip.vapourpt-ip.soliduspt)) 	  	                                # The lambda constant
     Hliq = D*Tliq + D*lambf
     alpha = ip.k/(ip.rho*ip.c)
     length = ip.Iref/(ip.k*(ip.vapourpt-ip.soliduspt))

tm = (ip.Iref**2*Ta**2*np.pi)/(4*ip.I**2)
t_start = delt