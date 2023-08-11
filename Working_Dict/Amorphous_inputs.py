# Amorphous_Material_Constants
# Default Material: - SS304L
soliduspt   = 1670               #Solidus Point                             # Kelvin
liquidouspt = 1727               #Liqudus Point                             # Kelvin
vapourpt    = 3375               #Vapour Point                              # Kelvin
Iref        = 1.5*10**10         #Reference Intensity                       # W/m^2
I           = 1.0*10**10         #Used Intensity                            # W/m^2
k           = 22                 #Conductivity                              # W/m-K 
delt        = 0.001              #time step size       
#===========================================================================================================================#
#For Stephan Problem: - n =  number of elements in liquid subdomain
#                       m =  number of elements in solid subdomain
#===========================================================================================================================#                             
n           = 21                  #number of nodes/elements
m           = 5                  #number of elements for solid subdomain
#---------------------------------------------------------------------------------------------------------------------------#
c           = 700                #specific heat capacity                    # in J/kg/K
rho         = 7200               #density of the material                   # kg/m**3  
Lf          = 261000             #latent  heat of fusion                    # J/kg   
depth       = 1                  #total depth to be simulated               # meters  
t_end       = 10                 #End time for simulation                   # dimensionless