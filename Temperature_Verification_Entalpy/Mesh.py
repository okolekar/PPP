'''Topic: PPP
Library: - Mesh 
#############################################################################################################################
Importing the required standard libraries 
-----------------------------------------------------------------------------------------------------------------------------
Processed_Inputs -> Script where all the inputs are defined                                                                  
#############################################################################################################################
'''
import numpy as np
import Processed_Inputs as ip
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
    runcount    ->      Checks if the method was already called before
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
    def __init__(self):
        self.nl = self.mesh_list()
        self.phi = None
        self.dphi = None
        self.J = None
        self.run_count = 0
    '''
#############################################################################################################################
Method test_phi(): -
=============================================================================================================================
The following tests were performed: -
1) Interpolation Properity: - checks whether the value of the Shape function at it's node is 1 and at other node 0.
2) Partition of unity property: - Checks if the sum of all the shape functions is 1 for any value passed'''
#############################################################################################################################
    def test_phi(self,zeta):
        if zeta <=1 and zeta >= -1:
            # Check Interpolation property
            if self.phi.any() > 1 or self.phi.any() < 0:
                raise ValueError(f"Incorrect Shape function, value exceeds 1 or is negative {self.phi}")
            elif zeta == 1 and (self.phi[0] != 0 or self.phi[1] != 1):
                raise ValueError(f"Incorrect Shape function for the given zeta {self.phi}")
            elif zeta == -1 and (self.phi[0] != 1 or self.phi[1] != 0):
                raise ValueError(f"Incorrect Shape function for the given zeta {self.phi}")
            #Check partition of unity property
            if self.phi.sum()!=1:
                raise ValueError(f"Incorrect Shape function for the given zeta {self.phi}")
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
Checks if the sum of the derivative of the shape function at the constant open coefficient is zero'''
#############################################################################################################################
    def test_dphi(self):
        ans = np.sum(self.dphi)
        if ans != 0:
            print("test_dphi failed")
            print("Ensure that the temperature field was constant at both nodes")
            raise ValueError(f"Shape function evaluated to be {ans}")
        else:
            print("\nIn the Mesh module")
            print("The derivative of Shape function passed the test.")
            '''
#############################################################################################################################
Method mesh_list(): -
=============================================================================================================================
1) Generates the mesh list 
2) Checks if the mesh list was properly created by ensuring that last entry of the mesh list is equal to the total length'''
#############################################################################################################################
    def mesh_list(self):
        l = []
        for q in range(ip.n):
            l.append(q*ip.length/(ip.n-1))
        if l[-1] != ip.length or len(l) == 0:
            print("\nError in Mesh module, the mesh list is not generated correctly")           
        return l

    def update_phi(self,zeta):
        if zeta <=1 and zeta >= -1:
            self.phi=np.array([0.5*(1-zeta),0.5*(1+zeta)]).reshape(2,1) # The shape function is Linear.
        else: 
            raise ValueError("\nError in Mesh module, value of zeta should be within [-1,1]")
        if self.run_count < 1:
            self.test_phi(zeta)
    '''
#############################################################################################################################
Method update_Jacobi(): -
=============================================================================================================================
1) Updates the  Jacobi
2) Checks if the determinant of the Jacobi is positive and non zero'''
#############################################################################################################################
    def update_Jacobi(self,i):
        x1 = self.nl[i]
        x2 = self.nl[i+1]
        self.J = 0.5*(x2-x1)
        if self.J <= 0 or self.J == 0:
            raise ValueError("\nError in Mesh module,the determinant of Jacobi is negative or zero")
        
    def update_dphi(self): # The shape function's first and second order partial deravatives.
        self.dphi = (self.J**-1)*np.array([-0.5,0.5]).reshape(2,1)
        if self.run_count < 1:
            self.test_dphi()
    
    def update_element(self,i):
        if  i<0 or type(i) != int:
            raise ValueError("\nError in Mesh module, Incorrect element number")
        self.update_Jacobi(i)
        self.update_dphi()