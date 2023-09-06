'''Topic: PPP
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
    N                   ->      The bi-Linear shape function
    dN                  ->      The derivative of bi-Linear shape function
    J                   ->      The Jacobian
    det_J               ->      The determinant of Jacobian
    node_list           ->      Stores the node number for each element
    elementCordiante    ->      Stores the x,y co-ordinate for each element 
-----------------------------------------------------------------------------------------------------------------------------
Main Methods : -
=============================================================================================================================
__mesh_list     ->  Generates mesh data
update_N        ->  Updates the Shape function as per the gauss points
update_Jacobi   ->  Updates the Jacobi and the determinant of the Jacobi as per the element number
update_dN       ->  Updates the derivative of the shape function as per the element number.
update_element  ->  Runs update_Jacobi, update_N and update_dN with a condition that the element number is non negative.
-----------------------------------------------------------------------------------------------------------------------------
Private Methods : -
=============================================================================================================================
__private_co_ornidate_list -> Generates the x and y co ordinates used internally for mesh generation
#############################################################################################################################
'''
class Mesh():
    def __init__(self,n,m):
        self.xy_list          = []
        self.boundary_nodes   = {}
        self.N                = None
        self.dN               = None
        self.J                = None
        self.det_J            = None
        self.node_list        = {}
        self.elementCordiante = {}
        self.__mesh_list(n,m)
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
    def __mesh_list(self,n,m):
        x = self.__private_co_ornidate_list(n,4)
        y = self.__private_co_ornidate_list(m,4)
        k = 0
        for j in range(len(y)-1):
            for i in range(len(x)-1):
                self.node_list[k] = [j*n+i, j*n+i+1, (j+1)*n+i+1, (j+1)*n+i]
                self.elementCordiante[k] = [[x[i],y[j]],[x[i+1],y[j]],[x[i+1],y[j+1]],[x[i],y[j+1]]]
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
        self.J = np.matmul(self.dN,self.elementCordiante[i])
        self.det_J = np.linalg.det(self.J)
    '''
#############################################################################################################################
Method update_dN(self,zeta,eta): -
=============================================================================================================================
1) Updates the derivative of Shape function for the given zeta and eta values
#############################################################################################################################
'''
    def update_dN(self,zeta,eta): 
        self.dN = (1/4)*np.array([[-(1-eta),(1-eta),(1+eta),-(1+eta)],
                                    [-(1-zeta),-(1+zeta),(1+zeta),(1-zeta)]])
    '''
#############################################################################################################################
Method update_N(self,zeta,eta): -
=============================================================================================================================
1) Updates the  Shape function for the given zeta and eta values
#############################################################################################################################
'''
    def update_N(self,zeta,eta):
        if zeta <=1 and zeta >= -1 and eta <=1 and eta >= -1:
            self.N=np.array([1/4*(1-zeta)*(1-eta),1/4*(1+zeta)*(1-eta),1/4*(1+zeta)*(1+eta),1/4*(1-zeta)*(1+eta)]).reshape(-1,1) # The shape function is Linear.
        else: 
            raise ValueError("\nError in Mesh module, value of zeta and eta should be within [-1,1]")
    '''
#############################################################################################################################
Method update_element(self,zeta,eta,i): -
=============================================================================================================================
updates the Shape function, derivative of the Shape function and the Jacobian
#############################################################################################################################
'''
    def update_element(self,zeta,eta,i):
        self.update_N(zeta,eta)
        self.update_dN(zeta,eta)
        self.update_Jacobi(i)
