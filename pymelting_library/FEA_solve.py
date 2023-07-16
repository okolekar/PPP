import numpy as np
import inputs as ip
import postProcessor as pp
from Initial_Temperature_Dist import Initial_Temp as IT
from variableUpdater import  Update_amorphus as update
from Elemental_Subroutine import Elemental_Subroutine
from Mesh import Mesh

node_list = Mesh(ip.n).nl
class FEA():
    def __init__(self):
        self.Tg = IT(ip.ratioI, ip.Ta, ip.delt, node_list).T #The dot T shows that I access the Temperature variable.
        self.Hg = update(D=ip.D,lambf= ip.lambf,Tliq= ip.Tliq,T = self.Tg,H=None,Hliq=None) #subscript g shows its a global vector
    def solve(self):
        E = Elemental_Subroutine()
        for t in range(1):
            nrs = 0
            Hk = self.Hg[:,t][:,np.newaxis].copy()
            Hk_1 = Hk.copy()
            if not np.array_equal(Hk,self.Hg[:,-1][:,np.newaxis]): #Testing if the extraction was a success
                raise ValueError("Data extraction failed")
            Tk_1 = self.Tg[:,t][:,np.newaxis].copy()
            if not np.array_equal(Tk_1,self.Tg[:,-1][:,np.newaxis]): #Testing if the extraction was a success
                raise ValueError("Data extraction failed")
            while (nrs<10):
                Ht = Hk_1.copy()
                E.get_param(Hk,Hk_1,Tk_1)
                # Boundry Condition
                E.Gg[ip.n-1] = 0
                E.dGg[ip.n-1] = 0
                E.dGg[:,ip.n-1] = 0
                E.dGg[ip.n-1,ip.n-1] = 1
                Hk_1 = Hk_1 - np.matmul(np.linalg.inv(E.dGg),E.Gg) 
                Tk_1 = update(D=ip.D,lambf=ip.lambf,Tliq=ip.Tliq,H = Hk_1,Hliq=ip.Hliq,T=None)
                nrs = nrs + 1
                if(np.abs(Ht-Hk_1).all()<np.exp(-5)):
                    print(f"NRS converged at {nrs} step")
                    break
                elif(nrs == 9):
                    print('Convergence failed to achieve in 10 nrs steps')
            
            self.Tg = np.column_stack((self.Tg,Tk_1))
            self.Hg = np.column_stack((self.Hg,Hk_1)) 

EntlpyApp = FEA()
EntlpyApp.solve()
pp.plotTemprature(EntlpyApp.Tg[:,0],EntlpyApp.Tg[:,-1],1)
