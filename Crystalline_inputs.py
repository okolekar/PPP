# Crystalline_Material_Constants
# Default Material:- Aluminium
meltpt      = 9.3*10**2          #Melting Point                             # Kelvin
vapourpt    = 2.5*10**3          #Vapour Point                              # Kelvin
Iref        = 1.5*10**10         #Reference Intensity                       # W/m^2
I           = 1.5*10**10         #Used Intensity                            # W/m^2
k           = 2.3*10**2          #Conductivity                              # W/m-K
delt        = 1              #time step size                            # dimensionless
#===========================================================================================================================#
#For Stephan Problem: - 
# n         =###                 number of elements in liquid subdomain
# m         =###                 number of elements in solid subdomain
#===========================================================================================================================#                             
n           = 11                 #number of nodes
m           = 41                 #number of nodes for solid subdomain
#---------------------------------------------------------------------------------------------------------------------------#                       
c           = 9.0*10**2          #specific heat capacity                    # in J/kg/K
rho         = 2.7*10**3          #density of the material                   # kg/m**3  
Lf          = 3.6*10**5          #latent  heat of fusion                    # J/kg     
depth       = 1                  #total depth to be simulated               # meters
#ds_dt       = (0.5/0.53)         #Initial velocity of the the melt front    # dimensionless
t_end       = 1                 #End time for simulation                   # dimensionless