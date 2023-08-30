'''Topic: PPP Name title and the start date

#############################################################################################################################
Importing the required libraries and their functions: -
=============================================================================================================================
Processed_Inputs: ->              Script where all the inputs are defined
postProcessor: ->                 Graphs and visualization library, 
Heat_Transfer_verification: ->    Provides an analytical solution to verify the application Heat transfer
Initial_Temperature_Dist: ->      Provides the Initial Temperature Distribution 
variableUpdater: ->               Updates the Temperature and the Enthalpy 
Elemental_Subroutine: ->          Extracts the Temperatures and Enthalpy form the Global vector for each element
Mesh: ->                          Generates the Mesh specific parameters '''
#############################################################################################################################
import numpy as np
import Processed_Inputs as ip
import postProcessor as pp
from Heat_Transfer_verification import Verification as HTv
from Initial_Temperature_Dist import Initial_Temp as IT
from Elemental_Subroutine import Elemental_Subroutine
from Mesh import Mesh
if ip.mat_type == 1:
    from variableUpdater import  Update_Crys as update
elif ip.mat_type ==2:
    from variableUpdater import  Update_amorphus as update
else:
    raise ValueError('Incorrect input value for the mateial type')

np.set_printoptions(edgeitems=11, threshold=np.inf, linewidth=np.inf, suppress=True, formatter={'int': lambda x: str(x)})
'''
#############################################################################################################################
Class FEA: -
=============================================================================================================================
Attributes: - 
-------------
    Mesh   ->   Instance of Mesh class stores mesh related parameters
      Tg   ->   Global temperature vector stores the temperature distribution for each time step
      Hg   ->   Global enthalpy vector stores the enthalpy distribution for each time step
        *subscript g shows its a global vector
-----------------------------------------------------------------------------------------------------------------------------
1) The FEA class is the starting point of the numerical scheme.
2) The main goal of this class is to calculate Tg and Hg vectors for given duration of time using the solve method.          
#############################################################################################################################
'''
class FEA():
    def __init__(self):
        self.Mesh = Mesh()
        self.Tg = IT(self.Mesh.nl).T                    #The dot T shows that I access the Temperature variable.
        self.Hg = update(T = self.Tg,H=None,run = 0)                                   
    '''
#############################################################################################################################
Method solve(): -
=============================================================================================================================
Variables: -
nrs         ->     Newton Raphson Steps recorder (Breaks the NRS loop if steps are more than 9)
t           ->     varible for time iteration
End_time    ->     Total time for which Simulation runs
E           ->     Instance of Elemental subroutine class
Hk          ->     Enthalpy calculated in the previous time Step
Hk+1        ->     Enthalpy calculated in the current Newton Raphson Step
Ht          ->     Enthalpy in the previous NRS step, acts as a condition to break the NRS loop
Tk          ->     Temperature calculated in the previous Newton Raphson Step
Tk-1        ->     Temperature calculated in the current Newton Raphson Step
-----------------------------------------------------------------------------------------------------------------------------
The method works as follows: -
1) Extraction of the latest temperature and enthalpy values from the global temperature Tg and enthalpy Hg vectors 
   respectively by Tk and Hk. Setting Hk_1 = Hk for the first NRS step
2) Extracted variables are passed to the Elemental_Subroutine and Gg and dGg vectors (Residual and Tangent Stiffness Matrix
   like vectors) are obtained
3) Application of the boundary condition and finding the new Hk_1 and Tk_1
4) Comparing for convergence and increasing the time                                                                      '''
#############################################################################################################################   
    def solve(self):
        E = Elemental_Subroutine()
        End_time = (1 if ip.test_case ==1 else ip.t_end)
        for t in range(End_time):
            nrs = 0
            Hk = self.Hg[:,t][:,np.newaxis].copy()
            Hk_1 = Hk.copy()
            if not np.array_equal(Hk,self.Hg[:,-1][:,np.newaxis]):      #Testing if the extraction was a success
                raise ValueError("Data extraction failed")
            Tk_1 = self.Tg[:,t][:,np.newaxis].copy()
            if not np.array_equal(Tk_1,self.Tg[:,-1][:,np.newaxis]):    #Testing if the extraction was a success
                raise ValueError("Data extraction failed")
            while (nrs<10):
                Ht = Hk_1.copy()
                E.get_param(Hk,Hk_1,Tk_1,self.Mesh)
#===========================================================================================================================#
                # Application of the Dirichlet Boundary Condition at the last node.
#===========================================================================================================================#
                E.Gg[ip.n-1] = 0
                E.dGg[ip.n-1] = 0
                E.dGg[:,ip.n-1] = 0
                E.dGg[ip.n-1,ip.n-1] = 1
                Hk_1 = Hk_1 - np.matmul(np.linalg.inv(E.dGg),E.Gg)
                #Hk_1[1] = (Hk_1[0] + Hk_1[2])*0.5 
                Tk_1 = update(H = Hk_1,T=None,run = t)
                nrs = nrs + 1
                if nrs == 1 and t == 0:
                    print("\nAll tests were successfully performed and no abnormalities detected\n")
                if(np.all(np.abs(Ht-Hk_1)<np.exp(-5))):
                    print(f"Solution converged in {nrs} step, for t = {t}")
                    break
                elif(nrs > 9):
                    raise ValueError('Convergence failed to achieve in 10 nrs steps')
            self.Tg = np.column_stack((self.Tg,Tk_1))
            self.Hg = np.column_stack((self.Hg,Hk_1)) 

EntlpyApp = FEA()
EntlpyApp.solve()
signs = np.sign(EntlpyApp.Tg)
signs = signs.T
mushy = []
for i in range(signs.shape[0]):
    change = np.where(np.diff(signs[i])<0)[0]
    if len(change) == 0:
        mushy.append(0)
    else :
        mushy.append(change.item()+1)
pp.plotxy(mushy,'node number','Evolution of Mushy Zone over time','Time step')


if ip.test_case == 1:
    Tverify = HTv(EntlpyApp.Tg[:,1],ip.n,EntlpyApp.Mesh.nl,ip.alpha).T_analytical
    pp.plotVerify(EntlpyApp.Tg[:,-1],
              'Final Temperature from the numerical scheme',Tverify,'Analytical Temperature')
else:
    pp.plotTemprature(EntlpyApp.Tg[:,0],EntlpyApp.Tg[:,-1],1)
    pp.plotxy(EntlpyApp.Tg[:,-1],'Dimentionless temperature','Final Temperature Distribution','Node Number')
