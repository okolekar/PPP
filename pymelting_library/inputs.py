soliduspt = 1670  #in K
liquidouspt = 1727
vapourpt  = 3375
k = 22 #Conductivity W/mK
delt = 0.01
ratioI = 1
n = 21                 #number of nodes  
c = 700              #specific heat capacity         # in J/kg/K 
rho = 7200              #density of the material        #kg/m**3
Tliq = (liquidouspt-soliduspt)/(vapourpt-soliduspt)               #Dimensionless 	  	            #Liq. Temp.          
Lf = 261000             #J/kg
D = rho*c*(vapourpt-soliduspt)/((rho*c*vapourpt+Lf*rho)-(rho*c*soliduspt))  #0.7970
Ta = 0 #(293-soliduspt)/(vapourpt-soliduspt)             # ambient temperature in Dimensionless.
lambf = Lf/(c*(vapourpt-soliduspt)) 	  	  # The lambda constant 
Hliq = D*Tliq + D*lambf
alpha = k/(rho*c)
