'''
Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
#############################################################################################################################
Importing the required libraries and their functions: -
=============================================================================================================================
(Standard libraries like numpy and others are not described here.)
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs:           ->    Source script for non dimensionalised constants
postProcessor:              ->    Graphs and visualization library, 
Heat_Transfer_verification: ->    Provides an analytical solution to verify the application Heat transfer and the melt front
Initial_Temperature_Dist:   ->    Provides the Initial Temperature Distribution 
variableUpdater:            ->    Updates the Temperature and the Enthalpy 
Elemental_Subroutine:       ->    Performs the element subroutine
Mesh:                       ->    Generates the Mesh specific parameters
#############################################################################################################################
'''
import numpy as np
import time
import matplotlib.pyplot as plt
import Processed_Inputs as ip
from Initial_Temperature_Dist import Initial_Temp as IT
from Elemental_Subroutine import Elemental_Subroutine
from Gridify import Mesh
if ip.mat_type == 1:
    from variableUpdater import  Update_Crys as update
elif ip.mat_type ==2:
    from variableUpdater import  Update_amorphus as update
else:
    raise ValueError('Incorrect input value for the mateial type')
import matplotlib.pyplot as plt
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
        #------------------------Meshing--------------------------------#
        Mstart = time.time()
        self.Mesh = Mesh(ip.n,ip.m,1,1)
        Mend = time.time()
        print(f' The time taken to mesh is found to be {Mend - Mstart}')
        #---------------------Initial_Condition--------------------------#
        ITstart = time.time()
        K = IT(np.array(self.Mesh.xy_list)[:,1],self.Mesh)
        self.Tg = K.T                    #The dot T shows that I access the Temperature variable.
        ITend = time.time()
        print(f' The time taken to find the Initial Temperature Distribution is found to be {ITend - ITstart}')
        #---------------------------------------------------------------#
        updatestart = time.time()
        self.Hg = update(T = self.Tg,H=None,run = 0)
        updateend = time.time()
        print(f' The time taken to Update is found to be {updateend - updatestart}')
#===========================================================================================================================#    
    def solve(self):
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
                is_symmetric = np.array_equal(E.dGg, E.dGg.T)
                print(f"The global Tangent stiffness matrix is symmetric: =  {is_symmetric}")
                if not is_symmetric:
                    print("Unsymmetric Global Tangent Stiffness Matrix Symmetric")
                for i in self.Mesh.boundary_nodes["Bottom_Nodes"]:
                    E.Gg[i] = 0
                    E.dGg[i,:] = 0
                    E.dGg[:,i] = 0
                for i in self.Mesh.boundary_nodes["Right_Nodes"]:
                    E.Gg[i] = 0
                    E.dGg[i,:] = 0
                    E.dGg[:,i] = 0
                for i in self.Mesh.boundary_nodes["Bottom_Nodes"]:
                    E.dGg[i,i] = 1
                for i in self.Mesh.boundary_nodes["Right_Nodes"]:
                    E.dGg[i,i] = 1
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
                            Vect_Temp = []
                            vectory = []
                            for i in self.Mesh.boundary_nodes["Left_Nodes"]:
                                Vect_Temp.append(self.Tg[i][-1])
                                vectory.append(self.Mesh.xy_list[i][1])

                            plt.plot(vectory,Vect_Temp,label  = "Temperature after time 't' ")
                            plt.title("Evolution of temperature just below the heat source")
                            plt.xlabel("dimensionless depth")
                            plt.ylabel("Dimensionless temperature")
                            plt.show()
                            raise ValueError('Convergence failed to achieve in 20 nrs steps')
                        
            self.Tg = np.column_stack((self.Tg,Tk_1))
            self.Hg = np.column_stack((self.Hg,Hk_1))

EntlpyApp = FEA()
start = time.time()
EntlpyApp.solve()
end = time.time()
print(f'The total time is{end-start}')
vectorx = np.array(EntlpyApp.Mesh.xy_list)[:,0]
vectory = np.array(EntlpyApp.Mesh.xy_list)[:,1]
vectorT = EntlpyApp.Tg[:,-1]
plt.scatter(vectorx,vectory,c=vectorT,cmap='coolwarm', s=100)
plt.xlabel('X Axis Label')
plt.ylabel('Y Axis Label')
plt.title('Temperature Scatter Plot')

# Add a colorbar to indicate temperature values
plt.colorbar(label='Temperature (°C)')

# Show the plot
plt.show()

Vect_Temp = []
for i in EntlpyApp.Mesh.boundary_nodes["Left_Nodes"]:
    Vect_Temp.append(EntlpyApp.Tg[i][-1])

plt.plot(Vect_Temp,label  = "Temperature after time 't' ")
plt.title("Evolution of temperature just below the heat source")
plt.xlabel("Node")
plt.ylabel("Dimensionless temperature")
plt.show()