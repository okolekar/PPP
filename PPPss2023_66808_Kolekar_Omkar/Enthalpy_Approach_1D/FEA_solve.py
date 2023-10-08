'''
Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
Library: Main Program
#############################################################################################################################
Importing the required libraries and their functions: -
=============================================================================================================================
(Standard libraries like numpy and others are not described here.)
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs:           ->    Source script for non dimensionalised constants
postProcessor:              ->    Graphs and visualization library, 
Initial_Temperature_Dist:   ->    Provides the Initial Temperature Distribution 
Elemental_Subroutine:       ->    Performs the element subroutine
Mesh:                       ->    Generates the Mesh specific parameters
variableUpdater:            ->    Updates the Temperature and the Enthalpy 
#############################################################################################################################
'''
import numpy as np
import time
import matplotlib.pyplot as plt
import Processed_Inputs as ip
import postProcessor as pp
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
#___________________________________________________________________________________________________________________________#

class FEA():
    '''
#############################################################################################################################
Class FEA: -
=============================================================================================================================
Attributes: - 
-------------
    Mesh              ->   Instance of Mesh class stores mesh related parameters
    Tg                ->   Global temperature vector stores the temperature distribution for each time step
    Hg                ->   Global enthalpy vector stores the enthalpy distribution for each time step
    analytical_s      ->   Stores the analytical phase front value
    Analytical_time   ->   Stores the analytical time for the evolution of the phase front
    time              ->   Stores the numerical time for the evolution of the phase front
        *subscript g shows its a global vector
-----------------------------------------------------------------------------------------------------------------------------
Variables: - 
------------
    K  ->   Instance of the Initial Temperature Distribution class
-----------------------------------------------------------------------------------------------------------------------------
1) The FEA class is the starting point of the numerical scheme.
2) The main goal of this class is to calculate Tg and Hg vectors for given duration of time using the solve method.          
#############################################################################################################################
'''
    def __init__(self):
        #-------------------------------------------#Creating Mesh#--------------------------------------------#
        Mstart      = time.time()
        self.Mesh   = Mesh()
        Mend        = time.time()
        print(f' The time taken to mesh is found to be {Mend - Mstart}')

        #----------------------------------#Defining the Initial Temperature#----------------------------------#
        ITstart     = time.time()
        K           = IT(self.Mesh.nl)                                                                                    
        self.Tg     = K.T                    #The dot T shows that I access the Temperature variable.
        ITend       = time.time()
        print(f' The time taken to find the Initial Temperature Distribution is found to be {ITend - ITstart}')

        #---------------------------------------#Updating the Enthalpy#----------------------------------------#
        updatestart         = time.time()
        self.Hg             = update(T = self.Tg,H=None,run = 0)
        updateend           = time.time()
        print(f' The time taken to Update is found to be {updateend - updatestart}')

        #-----------------------------------#Finding the Analytical Solution#----------------------------------#
        K.formulate_s()
        self.analytical_s    = K.s_array
        self.Analytical_time = K.hist_t
        self.time            = [0]
#===========================================================================================================================#    
    def solve(self):
        '''
#############################################################################################################################
Method solve(): -
=============================================================================================================================
Variables: -
E           ->     Instance of Elemental subroutine class
End_time    ->     Total time for which Simulation runs
tempt       ->     Dummy time variable to track time
t           ->     varible for time iteration and new Hg and Tg extraction
nrs         ->     Newton Raphson Steps recorder (Breaks the NRS loop if steps are more than 9)
Hk          ->     Enthalpy calculated in the previous time Step
Hk_1        ->     Enthalpy calculated in the current Newton Raphson Step
Ht          ->     Enthalpy in the previous NRS step, acts as a condition to break the NRS loop
Tk          ->     Temperature calculated in the previous Newton Raphson Step
Tk_1        ->     Temperature calculated in the current Newton Raphson Step
-----------------------------------------------------------------------------------------------------------------------------
The method works as follows: -
1) Extraction of the latest temperature and enthalpy values from the global temperature Tg and enthalpy Hg vectors 
   respectively by Tk and Hk. Setting Hk_1 = Hk for the first NRS step
2) Extracted variables are passed to the Elemental_Subroutine and Gg and dGg vectors (Residual and Tangent Stiffness Matrix
   like vectors) are obtained
3) Application of the boundary condition and finding the new Hk_1 and Tk_1
4) Comparing for convergence and increasing the time                                                                         
#############################################################################################################################
'''
        E = Elemental_Subroutine()
        End_time = (1 if ip.test_case ==1 else ip.t_end)
        tempt = 0       #Dummy variable to track the time. 
#===========================================================================================================================#
                                            # Time Discretisation Scheme
#===========================================================================================================================#
        for t in range(End_time):
            self.time.append(tempt)
            nrs = 0
            Hk = self.Hg[:,t][:,np.newaxis].copy()
            Hk_1 = Hk.copy()

            if not np.array_equal(Hk,self.Hg[:,-1][:,np.newaxis]):      #Testing if the extraction was a success
                raise ValueError("Data extraction failed")
            
            Tk_1 = self.Tg[:,t][:,np.newaxis].copy()
            if not np.array_equal(Tk_1,self.Tg[:,-1][:,np.newaxis]):    #Testing if the extraction was a success
                raise ValueError("Data extraction failed")
#===========================================================================================================================#
                                              # Newton Raphson Scheme
#===========================================================================================================================#
            while (nrs<20):
                print(f"Solution is at {nrs} nrs step, for t = {t}")
                Ht = Hk_1.copy()
                if t == 0:
                    Estart = time.time()
                    E.get_param(Hk,Hk_1,Tk_1,self.Mesh)
                    Eend = time.time()
                    print(f'The time to update the Elemental Subroutine is found to be {Eend - Estart}')
                else:
                    E.get_param(Hk,Hk_1,Tk_1,self.Mesh)
#===========================================================================================================================#
                            # Application of the Dirichlet Boundary Condition at the last node.
#===========================================================================================================================#
                E.Gg[ip.n-1] = 0
                E.dGg[ip.n-1,:] = 0
                E.dGg[:,ip.n-1] = 0
                E.dGg[ip.n-1,ip.n-1] = 1        #Diagonal entry at the last node of the tangent stiffness matrix = 1
#---------------------------------------------------------------------------------------------------------------------------#
                Hk_1 = Hk_1 - np.matmul(np.linalg.inv(E.dGg),E.Gg) 
                Tk_1 = update(H = Hk_1,T=None,run = t)
                nrs = nrs + 1
                if nrs == 1 and t == 0:
                    print("\nAll tests were successfully performed and no abnormalities detected\n")
                if(np.all(np.abs(Ht-Hk_1)<np.exp(-5))):
                    print(f"Solution converged in {nrs} step, for t = {t}")
                    tempt = tempt + ip.delt
                    break
                elif(nrs>19):
                    print(f'The nrs step was found to be {nrs}')  
                    signs = np.sign(self.Hg)
                    signs = signs.T
                    mushy = []
                    for i in range(signs.shape[0]):
                        change = np.where(np.diff(signs[i])<0)[0]
                        if len(change) == 0:
                            mushy.append(0)
                        else :
                            mushy.append(change[-1].item()+1)
                    mushy = np.array(mushy)
                    mushy = mushy.astype('float')
                    fig, ax = plt.subplots()
                    ax.plot(mushy,marker='o',linestyle='-',label='Enthalpy Problem Approach')
                    ax.plot(EntlpyApp.analytical_s,marker='o',linestyle='-',label='Analytical Solution')
                    ax.set_xlabel('Dimensionless Time')
                    ax.set_ylabel('Dimensionless depth')
                    ax.set_title('Evolution of the position of the solid liquid interface')
                    ax.legend()
                    plt.show()
                    raise ValueError('Convergence failed to achieve in 20 nrs steps')
                        
            self.Tg = np.column_stack((self.Tg,Tk_1))
            self.Hg = np.column_stack((self.Hg,Hk_1))
#===========================================================================================================================#
#---------------------------------------------------------#End Of Class#----------------------------------------------------#
#===========================================================================================================================#

EntlpyApp = FEA()
start = time.time()
EntlpyApp.solve()
end = time.time()
print(f'The total time is{end-start}')

signs = np.sign(EntlpyApp.Tg)
signs = signs.T
mushy = []
for i in range(signs.shape[0]):
    change = np.where(np.diff(signs[i])<0)[0]
    if len(change) == 0:
        mushy.append(0)
    else :
        mushy.append(change[-1].item()+1)
mushy = np.array(mushy)
mushy = mushy.astype('float')
mushy = mushy*(EntlpyApp.Mesh.nl[1]-EntlpyApp.Mesh.nl[0])
poststart = time.perf_counter()
postend = time.perf_counter()
print(f'The post processer takes around {postend-poststart} time')
fig, ax = plt.subplots()
ax.plot(EntlpyApp.time,mushy,label='Position of Interface from Enthalpy Problem Approach')
ax.plot(EntlpyApp.Analytical_time,EntlpyApp.analytical_s,label='Position of Interface from Analytical Solution')
ax.set_xlabel('Dimensionless Time')
ax.set_ylabel('Dimensionless depth')
ax.set_title('Evolution of the position of the solid liquid interface')
ax.legend()
plt.show()

pp.plotTemprature(EntlpyApp.Mesh.nl,EntlpyApp.Tg[:,0],EntlpyApp.Tg[:,-1],ip.t_end)
fig, ax = plt.subplots()
ax.scatter(EntlpyApp.Tg[1,:],EntlpyApp.Hg[1,:],marker = '*',label='Evolution of the Enthalpy at node 1')
ax.set_xlabel('Dimensionless Temperature')
ax.set_ylabel('Dimensionless Enthalpy')
ax.set_title('Evolution of the Enthalpy vs Temperature over time')
ax.legend()
plt.show()