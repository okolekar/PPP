#Heat_Transfer_verification

from Initial_Temperature_Dist import Initial_Temp as IT
import numpy as np
import inputs as ip
from Mesh import Mesh
import postProcessor as pp

#T_bound = -ip.soliduspt/(ip.vapourpt - ip.soliduspt) #boundary temperature assumed to be 0K and dedimensionalized
nlist = Mesh(ip.n).nl
alpha = ip.k/(ip.rho*ip.c)

Tinitial = IT(ip.ratioI, ip.Ta, ip.delt, nlist).T
T_analytical = np.zeros(ip.n).reshape(ip.n,1)
for i in range(ip.n):
    x = nlist[i]
    #T_analytical[i] = (Tinitial[i]/(np.cos(np.pi*x)/2))*(np.cos(np.pi*x)/2)*np.exp((-np.pi**2*alpha**2)/(2))
    #T_analytical[i] = (Tinitial[i]/(np.cos(np.pi*x)/2))*(np.cos(np.pi*x)/2)*np.exp((-np.pi**2*alpha**2*1000)/(2))
    #T_analytical[i] = (Tinitial[0])*(np.sin(np.pi*x/2))*np.exp((-np.pi**2*alpha**2)/(2))
    #T_analytical[i] = (Tinitial[i])*(np.sin(np.pi*x/2))*np.exp((-np.pi**2*alpha**2)/(2))
    T_analytical[i] = 2*(Tinitial[i])*(np.cos(np.pi*x))*np.exp((-np.pi**2*alpha*ip.delt))



T = np.column_stack((Tinitial,T_analytical))


print(T)
pp.plotTemprature(Tinitial,T_analytical,1)