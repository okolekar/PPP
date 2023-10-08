'''
Topic:   Programming of the 2D FE Implementation of Stefan Problem and the Enthalpy Problem.
Program: The Enthalpy Problem               Matriculation Number: 66808
Library: - Gridify
#############################################################################################################################
Importing the required standard libraries                                                                
#############################################################################################################################
'''
import numpy as np
'''
#############################################################################################################################
Class Mesh: -
=============================================================================================================================
The attributes of this class store all the mesh related information
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    xy_list             ->      x,y co-ordinate list with index as the node number
    boundary_nodes      ->      Stores the boundary nodes
    elementList         ->      Stores the list of the elements 
    N                   ->      The bi-Linear shape function
    dN_zeta             ->      The derivative of bi-Linear shape function wrt local co-ordinates
    dN                  ->      The derivative of bi-Linear shape function wrt global co-ordinates
    J                   ->      The Jacobian
    det_J               ->      The determinant of Jacobian
    node_list           ->      Stores the node number for each element
    elementCordiante    ->      Stores the x,y co-ordinate for each element 
-----------------------------------------------------------------------------------------------------------------------------
Main Methods : -
----------------
=============================================================================================================================
__mesh_list     ->  Generates mesh data (Hidden)
update_N        ->  Updates the Shape function as per the gauss points
update_Jacobi   ->  Updates the Jacobi and the determinant of the Jacobi as per the element number
update_dN       ->  Updates the derivative of the shape function as per the element number.
update_element  ->  Runs update_Jacobi, update_N and update_dN with a condition that the element number is non negative.
-----------------------------------------------------------------------------------------------------------------------------
Private Methods : -
-------------------
=============================================================================================================================
__private_co_ornidate_list -> Generates the x and y co ordinates used internally for mesh generation
#############################################################################################################################
'''
class Mesh():
    def __init__(self,n,m,length_x,length_y):
        self.xy_list          = []
        self.boundary_nodes   = {}
        self.elementList      = []
        self.N                = None
        self.dN_zeta          = None
        self.dN               = None
        self.J                = None
        self.Jinv             = None
        self.det_J            = None
        self.node_list        = {}
        self.elementCordiante = {}
        self.__mesh_list(n,m,length_x,length_y)
    '''
#############################################################################################################################
Method __private_co_ornidate_list(): -
=============================================================================================================================
1) Generates the co ordinates for x and y direction
2) Checks if the list was properly created by ensuring that last entry of the mesh list is equal to the total length
#############################################################################################################################
'''
    def __private_co_ornidate_list(self,nodes,length):
        l = []
        for q in range(nodes):
            l.append(q*length/(nodes-1))
        if round(l[-1],5) != round(length,5) or len(l) == 0:
            print(f"\nError in Mesh module, the mesh list {l} is not generated correctly for length = {length}")
        return l
#===========================================================================================================================#
    def __mesh_list(self,n,m,length_x,length_y):
        x = self.__private_co_ornidate_list(n,length_x)
        y = self.__private_co_ornidate_list(m,length_y)
        k = 0
        for j in range(len(y)-1):
            for i in range(len(x)-1):
                self.node_list[k] = [j*n+i, j*n+i+1, (j+1)*n+i+1, (j+1)*n+i]
                self.elementCordiante[k] = [[x[i],y[j]],[x[i+1],y[j]],[x[i+1],y[j+1]],[x[i],y[j+1]]]
                self.elementList.append(k)
                k = k+1
        for j in range(m):
            for i in range(n):
                self.xy_list.append([x[i],y[j]])
        self.boundary_nodes['Bottom_Nodes'] = [q for q in range(n)]
        self.boundary_nodes['Top_Nodes'] = [(q + n*(m-1)) for q in range(n)]
        self.boundary_nodes["Right_Nodes"] = [n*q+n-1 for q in range(m)]
        self.boundary_nodes["Left_Nodes"] = [n*q for q in range(m)]
    '''
#############################################################################################################################
Method update_Jacobi(): -
=============================================================================================================================
1) Updates the  Jacobi
2) Updates the determinant of the Jacobian
#############################################################################################################################
'''
    def update_Jacobi(self,i):
        self.J = np.matmul(self.dN_zeta,self.elementCordiante[i])
        self.det_J = np.linalg.det(self.J)
        self.Jinv = np.linalg.inv(self.J)
    '''
#############################################################################################################################
Method update_dN(self,zeta,eta): -
=============================================================================================================================
1) Updates the derivative of Shape function for the given zeta and eta values
#############################################################################################################################
'''
    def update_dN(self,zeta,eta,i): 
        self.dN_zeta = (1/4)*np.array([[-(1-eta ), (1-eta ),(1+eta ),-(1+ eta)],
                                       [-(1-zeta),-(1+zeta),(1+zeta), (1-zeta)]])
        if np.sum(self.dN_zeta[0,:])!=0:
            raise ValueError(f"The sum of the dN is incorrect {self.dN_zeta}")
        if np.sum(self.dN_zeta[1,:])!=0:
            raise ValueError(f"The sum of the dN is incorrect {self.dN_zeta}")
        self.update_Jacobi(i)
        self.dN = self.Jinv @ self.dN_zeta
        if round(np.sum(self.dN[0,:]),2) != 0.0:
            raise ValueError(f"The sum of the dN is incorrect {np.sum(self.dN[0,:])}")
        if round(np.sum(self.dN[1,:]),2) != 0.0:
            raise ValueError(f"The sum of the dN is incorrect {np.sum(self.dN[0,:])}")
    '''
#############################################################################################################################
Method update_N(self,zeta,eta): -
=============================================================================================================================
1) Updates the  Shape function for the given zeta and eta values
#############################################################################################################################
'''
    def update_N(self,zeta,eta):
        if zeta <=1 and zeta >= -1 and eta <=1 and eta >= -1:
            self.N=np.array([1/4*(1-zeta)*(1-eta),1/4*(1+zeta)*(1-eta),1/4*(1+zeta)*(1+eta),1/4*(1-zeta)*(1+eta)])\
                                                                               .reshape(-1,1) # The shape function is Linear.
        if np.sum(self.N)!=1:
            raise ValueError(f"The sum of the shape function is incorrect {self.N}")
    '''
#############################################################################################################################
Method update_element(self,zeta,eta,i): -
=============================================================================================================================
updates the Shape function, derivative of the Shape function and the Jacobian
#############################################################################################################################
'''
    def update_element(self,zeta,eta,i):
        self.update_N(zeta,eta)
        self.update_dN(zeta,eta,i)
'''
=============================================================================================================================
--------------------------------------------------------End of Class---------------------------------------------------------
=============================================================================================================================
'''