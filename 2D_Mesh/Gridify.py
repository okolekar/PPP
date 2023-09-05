'''Topic: PPP
Library: - Mesh 
#############################################################################################################################
Importing the required standard libraries 
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs -> Script where all the inputs are defined                                                                  
#############################################################################################################################
'''
import math
import numpy as np
#import Processed_Inputs as ip
'''
#############################################################################################################################
Class Mesh: -
=============================================================================================================================
The attributes of this class store all the mesh related information
-----------------------------------------------------------------------------------------------------------------------------
Attributes: -
-------------
    nl          ->      Node list
    phi         ->      The Linear shape function
    dphi        ->      The derivative of Linear shape function
    J           ->      The Jacobian
-----------------------------------------------------------------------------------------------------------------------------
Methods : -
=============================================================================================================================
test_phi        ->      Tests if the shape function was correctly formed formulated
test_dphi       ->      Tests if the derivative of the shape function was correctly formed formulated
mesh_list       ->      This method provides the node number and corresponding z co-ordinate value
update_phi      ->      Updates the Shape function as per the gauss points
update_Jacobi   ->      Updates the Jacobi as per the element number
update_dphi     ->      Updates the derivative of the shape function as per the element number.
update_element  ->      Runs update_Jacobi and update_dphi with a condition that the element number is non negative.
#############################################################################################################################
'''
class Mesh():
    def __init__(self,n,m):
        self.xy_list          = []
        self.boundary_nodes   = None
        self.N                = None
        self.dN               = None
        self.J                = None
        self.node_list        = {}
        self.elementCordiante = {}
        self.__mesh_list(n,m)
    '''
#############################################################################################################################
Method test_phi(): -
=============================================================================================================================
The following tests were performed: -
1) Interpolation Properity: - checks whether the value of the Shape function at it's node is 1 and at other node 0.
2) Partition of unity property: - Checks if the sum of all the shape functions is 1 for any value passed
#############################################################################################################################
'''
    def test_phi(self,zeta):
        if zeta <=1 and zeta >= -1:
#===========================================================================================================================#
                                            #Interpolation property Check
#===========================================================================================================================#
            if self.N.any() > 1 or self.N.any() < 0:
                raise ValueError(f"Incorrect Shape function, value exceeds 1 or is negative {self.N}")
            elif zeta == 1 and (self.N[0] != 0 or self.N[1] != 1):
                raise ValueError(f"Incorrect Shape function for the given zeta {self.N}")
            elif zeta == -1 and (self.N[0] != 1 or self.N[1] != 0):
                raise ValueError(f"Incorrect Shape function for the given zeta {self.N}")
#===========================================================================================================================#
                                           #Partition of unity property Check
#===========================================================================================================================#
            if self.N.sum()!=1:
                raise ValueError(f"Incorrect Shape function for the given zeta {self.N}")
            else:
                print("And following tests were performed on Shape function")
                print("1) Interpolation property check,")
                print("2) Partition of unity property check,")
                print("Shape function passed all the tests")
                print("\n")
                self.run_count += 1
        else: 
            raise ValueError("Value of zeta should be within [-1,1] and i should be int")
        '''
#############################################################################################################################
Method test_dphi(): -
=============================================================================================================================
Checks if the sum of the derivative of the shape function at the constant open coefficient is zero
#############################################################################################################################
'''
    def test_dphi(self):
        ans = np.sum(self.dN)
        if ans != 0:
            print("test_dphi failed")
            raise ValueError(f"Shape function evaluated to be {ans}")
        else:
            print("\nIn the Mesh module")
            print("The derivative of Shape function passed the test.")
            '''
#############################################################################################################################
Method __private_co_ornidate_list(): -
=============================================================================================================================
1) Generates the co ordinates for x and y co ordinates
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
    '''
#############################################################################################################################
Method update_Jacobi(): -
=============================================================================================================================
1) Updates the  Jacobi
2) Checks if the determinant of the Jacobi is positive and non zero
#############################################################################################################################
'''
    def update_Jacobi(self,i):
        self.J = np.matmul(self.dN,self.elementCordiante[i])
        """if self.J <= 0 or self.J == 0:
            raise ValueError("\nError in Mesh module,the determinant of Jacobi is negative or zero")"""
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
