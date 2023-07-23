import numpy as np
import inputs as ip
import postProcessor as pp
from Heat_Transfer_verification import Verification as HTv
from Initial_Temperature_Dist import Initial_Temp as IT
from variableUpdater import  Update_amorphus as update
from Elemental_Subroutine import Elemental_Subroutine
from Mesh import Mesh
from Check_Inputs import check_inputs

np.set_printoptions(edgeitems=11, threshold=np.inf, linewidth=np.inf, suppress=True, formatter={'int': lambda x: str(x)})

class FEA():
    @check_inputs
    def __init__(self):
        self.Mesh = Mesh()
        self.Tg = IT(self.Mesh.nl).T                    #The dot T shows that I access the Temperature variable. #it was 0.51
        self.Hg = update(T = self.Tg,H=None,Hliq=None)  #subscript g shows its a global vector
    
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
                E.get_param(Hk,Hk_1,Tk_1,self.Mesh)
                # Boundry Condition
                E.Gg[ip.n-1] = 0
                E.dGg[ip.n-1] = 0
                E.dGg[:,ip.n-1] = 0
                E.dGg[ip.n-1,ip.n-1] = 1
                Hk_1 = Hk_1 - np.matmul(np.linalg.inv(E.dGg),E.Gg) 
                Tk_1 = update(H = Hk_1,Hliq=ip.Hliq,T=None)
                nrs = nrs + 1
                if(np.all(np.abs(Ht-Hk_1)<np.exp(-5))):
                    print(f"NRS converged at {nrs} step")
                    break
                elif(nrs > 9):
                    raise ValueError('Convergence failed to achieve in 10 nrs steps')
            
            self.Tg = np.column_stack((self.Tg,Tk_1))
            self.Hg = np.column_stack((self.Hg,Hk_1)) 


EntlpyApp = FEA()
EntlpyApp.solve()

Tverify = HTv(EntlpyApp.Tg[:,0],ip.n,EntlpyApp.Mesh.nl,ip.alpha).T_analytical
#print(EntlpyApp.Tg)
#print(EntlpyApp.Hg)
#pp.plotTemprature(EntlpyApp.Tg[:,0],EntlpyApp.Tg[:,-1],1)
#pp.Simultaneous_Plotter(EntlpyApp.Tg,'Temperature',EntlpyApp.Hg,'Enthalpy')
pp.plotVerify(EntlpyApp.Tg[:,0],'Initial_Temperature',EntlpyApp.Tg[:,-1],'Final Temperature from the numerical scheme',Tverify,'Analytical Temperature')
