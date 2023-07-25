import numpy as np

soliduspt = 1670                                        #in K
liquidouspt = 1727                                      #in K
vapourpt  = 3375                                        #in K
Iref = 1.5*10**10                                       # W/m^2
I =    1.5*10**10                                       # W/m^2
k = 22                                                  #Conductivity W/mK
delt = 0.01                                             
ratioI = I/Iref                                         
n = 21                                                  #number of nodes  
c = 700                                                 #specific heat capacity         # in J/kg/K 
rho = 7200                                              #density of the material        #kg/m**3
Tliq = (liquidouspt-soliduspt)/(vapourpt-soliduspt)     #Dimensionless 	  	            #Liq. Temp.          
Lf = 261000                                             #J/kg
D = rho*c*(vapourpt-soliduspt)/((rho*c*vapourpt+Lf*rho)-(rho*c*soliduspt))  #0.7970
Ta = 0#(293-soliduspt)/(vapourpt-soliduspt)             # ambient temperature in Dimensionless.
lambf = Lf/(c*(vapourpt-soliduspt)) 	  	  # The lambda constant 
Hliq = D*Tliq + D*lambf
alpha = k/(rho*c)
length = 400*10**3                          
#tm is the time when the Temperature reaches the melting point
#t_start = time when we consider that the simulation starts
tm = (Iref**2*Ta**2*np.pi)/(4*I**2)
t_start = delt
