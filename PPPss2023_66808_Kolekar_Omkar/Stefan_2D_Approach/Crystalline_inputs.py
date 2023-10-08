# Crystalline_Material_Constants
# Default Material:- Aluminium
theta_m     = 931              #Melting Point                        # Kelvin
k           = 225.5            #Conductivity                         # W/m-K
delt        = 0.0005           #time step size                                                     
c           = 1.0162*10**3     #specific heat capacity               # in J/kg/K
rho         = 2545             #density of the material              # kg/m**3  
Lf          = 396*10**3        #latent  heat of fusion               # J/kg     
L           = 0.5              #total depth to be simulated          # meters
X           = 0.5
Y           = 0.5

alpha = k/(rho*c)
x = X/L
y = Y/L