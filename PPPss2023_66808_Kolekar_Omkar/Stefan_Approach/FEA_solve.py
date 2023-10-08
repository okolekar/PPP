'''
Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Stefan Problem               Matriculation Number: 66808
Library: - Main
#############################################################################################################################
Importing the required libraries and their functions: -
=============================================================================================================================
Processed_Inputs:           ->      Script where all the inputs are defined
postProcessor:              ->      Graphs and visualization library
variableUpdater:            ->      Updates the Temperature, Enthalpy, rate of melting and the melt distance
Initial_Conditions:         ->      Provides the Initial Conditions for both the approaches  
Elemental_Subroutine:       ->      Extracts the Temperatures and Enthalpy form the Global vector for each element
Mesh:                       ->      Generates the Mesh specific parameters
#############################################################################################################################
'''
import numpy as np
import Processed_Inputs as ip
import postProcessor as pp
import variableUpdater as vu
from Initial_Conditions import Initial_Temp as IC
from Initial_Conditions import Stefan_Initial_Conditions as SIC
from Element_Subroutine import Elemental_Subroutine
import matplotlib.pyplot as plt
from Mesh import Mesh
import time
#===========================================================================================================================#
                #The below section selects the material specific updater method for the enthalpy approach
#===========================================================================================================================#
if ip.mat_type   == 1:
    from variableUpdater import  Update_Crys as update
elif ip.mat_type == 2:
    from variableUpdater import  Update_amorphus as update
else:
    raise ValueError('Incorrect input value for the mateial type.')

np.set_printoptions(edgeitems=11, threshold=np.inf, linewidth=np.inf, suppress=True, formatter={'int': lambda x: str(x)})

#_____________________________________________*End of Import Section*_______________________________________________________#
'''
#############################################################################################################################
Class FEA: -
=============================================================================================================================
Attributes: -
-------------
    Enthalpy Approach: -
    ---------------------
    Mesh   ->   Instance of Mesh class stores mesh related parameters
      Tg   ->   Global temperature vector stores the temperature distribution for each time step
      Hg   ->   Global enthalpy vector stores the enthalpy distribution for each time step
        *subscript g shows it's a global vector

    Stefan Approach: -
    ------------------
    list_s -> Stores the distance of phase change for the each time step
    s      -> The current position of the solid liquid interface
    ds_dt  -> Velocity of the phase interface
    meshl  -> Stores the mesh related parameters for the liquid regime
    meshs  -> Stores the mesh related parameters for the solid regime
    Tl     -> The temperature distribution in the liquid regime
    Ts     -> The temperature distribution in the solid regime
-----------------------------------------------------------------------------------------------------------------------------
Description and Methods: -
--------------------------
1) The FEA class is the starting point of the numerical scheme.
    The Attributes are selected on the basis of the Approach which is passed as an argument to the class.
    If approach is s, then Stefan Approach is selected and the Stefan Approach specific parameters are selected else Enthalpy  
    Approach specific  parameters are selected.

2) The two main goals of this class are achieved by two methods respectively which are as follows: -
   a) solve:                    -> to calculate Tg and Hg vectors for given duration of time using the solve method for the
      ------                       Enthalpy Method and (Redudant hsere)

   b) Stefan_problem_solver:    -> to calculate the phase front position in the Stefan Approach
      ----------------------
-----------------------------------------------------------------------------------------------------------------------------          
#############################################################################################################################
'''
class FEA():
    def __init__(self,approach = None):
        if approach == None:
            self.Mesh   = Mesh()
            self.Tg     = IC(self.Mesh.nl).T                    #The dot T shows that I access the Temperature variable.
            self.Hg     = update(T = self.Tg,H=None,run = 0)
        elif approach == 's':
            self.list_s = []
            Initial     = SIC()
            self.s      = Initial.s_array[1]
            self.ds_dt  = Initial.ds_dt[1]
            self.t      = Initial.hist_t[1]
            self.histt  = [self.t]
            self.meshl  = Mesh(s = 0,nodes = ip.n,length = self.s)
            self.meshs  = Mesh(s = self.s,nodes = ip.m,length = ip.length)
            Initial.Tdist_Stefan(self.meshl.nl,self.meshs.nl,self.ds_dt,self.s)
            self.Tl     = np.array(Initial.Tl)[:, np.newaxis]
            self.Ts     = np.array(Initial.Ts)[:, np.newaxis]
        else:
            raise TypeError("Invalide Approach. Please enter s for Stefan Approach")                                  

    '''
#############################################################################################################################
Method Stefan_problem_solver(): -
=============================================================================================================================
Variables: -
Tl_History  ->  Stores the history of the Temperature distribution in the liquid phase
Ts_History  ->  Stores the history of the Temperature distribution in the solid phase
End_time    ->  Total time for which Simulation runs
E           ->  Instance of Elemental subroutine class
Ml,Nl,bl    ->  Liquid phase specific parameters at the previous time step
Ms,Ns,bs    ->  phase specific parameters at the previous time step
A,B,C       ->  Temporary variables used only for calculation purpose
LHS, RHS    ->  Left Hand Side and Right Hand Side of the liquid and solid Temperature distribution equation
t           ->  Time Step in the Euler Forward Time Integration Scheme
-----------------------------------------------------------------------------------------------------------------------------
The method works as follows: -
1) Based on the Initial conditions M,N,b values for the previous time step are calculated for both the phases and stored
   in the respective variables and the Temperature distribution is also stored in the History variable
2) The new melt front distance is updated and the phase change travel velocity is updated
3) Based on the new position of the melt front, the mesh is updated
4) Based on the new mesh, new Temperature distribution is calculated in both the phases
5) Process is repeated untill the end time is reached                                                                      
#############################################################################################################################
'''  
    def Stefan_problem_solver(self):
        node_l_history = self.meshl.nl
        node_s_history = self.meshs.nl
        history_tl = []
        history_ts = []
        dt_ds_tracker = []
        dt_ds_tracker.append(self.ds_dt)
        Tl_History = self.Tl.copy()
        Ts_History = self.Ts.copy()
        self.list_s.append(self.s)
        E = Elemental_Subroutine('s')
        E.get_global_Stefan_param(self.meshl,self.meshs,self.ds_dt)
#---------------------------------------------------------------------------------------------------------------------------#
        for tempo in range(500): 
            print(tempo)  
            Ml = E.Mlg.copy()
            Nl = E.Nlg.copy()
            bl = E.blg.copy()
            Ms = E.Msg.copy()
            Ns = E.Nsg.copy()
            bs = E.bsg.copy()
            ip.delt    = 0.00001
            history_tl.append(self.Tl[-2].item())
            history_ts.append(self.Ts[1].item())
            self.s = vu.s_update(self.s,self.Tl[-2].item(),self.Ts[1].item(), self.meshl.h,self.meshs.h)
            if self.s<0:
                pp.plotxy(self.list_s[:-2],"melt surface dist","Melt vs Time step","Time step")
                pp.plotxy(dt_ds_tracker[:-2],"ds/dt","evolution of speed","time step")
            self.meshl = Mesh(s = 0,nodes = ip.n,length = self.s)
            node_l_history = np.column_stack((node_l_history,self.meshl.nl))
            self.meshs = Mesh(s = self.s,nodes = ip.m,length = ip.length)
            node_s_history = np.column_stack((node_s_history,self.meshs.nl))
            E.get_global_Stefan_param(self.meshl,self.meshs,self.ds_dt)
#---------------------------------------------------------------------------------------------------------------------------#
            A = np.eye(ip.n)-0.5*ip.delt*np.matmul(np.linalg.inv(Ml),Nl)
            B = np.matmul(A,self.Tl)
            C = ip.delt*np.matmul((0.5*np.linalg.inv(E.Mlg)+0.5*np.linalg.inv(Ml)),bl)
            RHS = B + C
            A = np.matmul(np.linalg.inv(E.Mlg),E.Nlg)
            LHS = np.linalg.inv(np.eye(ip.n)+0.5*ip.delt*A)
            self.Tl = np.matmul(LHS,RHS)
#---------------------------------------------------------------------------------------------------------------------------#
        #Application of the Dirichlet boundary condition on the temperature distribution in the liquid regime
#---------------------------------------------------------------------------------------------------------------------------#
            self.Tl[-1] = 0            
            Tl_History = np.column_stack((Tl_History,self.Tl))
#---------------------------------------------------------------------------------------------------------------------------#
            A = np.eye(ip.m)-0.5*ip.delt*np.matmul(np.linalg.inv(Ms),Ns)
            B = ip.delt*(0.5*np.matmul(np.linalg.inv(E.Msg),E.bsg)+0.5*np.matmul(np.linalg.inv(Ms),bs))
            LHS = np.eye(ip.m)+0.5*ip.delt*np.matmul(np.linalg.inv(E.Msg),E.Nsg)
            RHS = np.matmul(A,self.Ts)+B
            self.Ts = np.matmul(np.linalg.inv(LHS),RHS)
#---------------------------------------------------------------------------------------------------------------------------#
        #Application of the Dirichlet boundary condition on the temperature distribution in the solid regime
#---------------------------------------------------------------------------------------------------------------------------#
            self.Ts[0] = 0   
            self.Ts[-1] = ip.Ta           
            Ts_History = np.column_stack((Ts_History,self.Ts))
#---------------------------------------------------------------------------------------------------------------------------#
            self.list_s.append(self.s)
            dt_ds_tracker.append(self.ds_dt)
            self.t = self.t + ip.delt
            self.histt.append(self.t)
        pp.plotTemprature(Tl_History[:,0],Tl_History[:,-1])
        pp.plotTemprature(Ts_History[:,0],Ts_History[:,-1])
#------------------------------------------------End of Stefan_problem_solver-----------------------------------------------#
'''
=============================================================================================================================
--------------------------------------------------------End of Class---------------------------------------------------------
=============================================================================================================================
'''
Initial     = SIC()
start = time.time()
Stefan = FEA(approach = 's')
Stefan.Stefan_problem_solver()
end = time.time()
print(f"The time taken to run the program is {end - start}")
print("\n")
fig, ax = plt.subplots()
ax.plot(Initial.hist_t,Initial.s_array,marker='o', linestyle='-',label='Analytical Scheme')
ax.plot(Stefan.histt,Stefan.list_s,marker='o', linestyle='-',label='Numerical Scheme')
ax.set_xlabel('Dimensionless Time')
ax.set_ylabel('Dimensionless Depth')
ax.set_title('Evolution of the phase front depth')
ax.legend()
plt.show()