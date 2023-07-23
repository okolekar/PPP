import numpy as np
import inputs as ip

class Mesh():
    def __init__(self):
        self.nl = self.mesh_list() #node list
        self.phi = None
        self.dphi = None
        self.J = None

    def test_phi(self,zeta):
        if zeta <=1 and zeta >= -1:
            self.update_phi(zeta)
            # Check Interpolation property
            if self.phi.any() > 1 or self.phi.any() < 0:
                raise ValueError(f"Incorrect Shape function value exceeds 1 or is negative {self.phi}")
            elif zeta == 1 and (self.phi[0] != 0 or self.phi[1] != 1):
                raise ValueError(f"Incorrect Shape function for the given zeta {self.phi}")
            elif zeta == -1 and (self.phi[0] != 1 or self.phi[1] != 0):
                raise ValueError(f"Incorrect Shape function for the given zeta {self.phi}")
            #Check partition of unity property
            if self.phi.sum()!=1:
                raise ValueError(f"Incorrect Shape function for the given zeta {self.phi}")
            if self.dphi.sum() != 0:
                raise ValueError(f"Incorrect derivative of Shape function")
            else:
                print("Following tests were performed on Shape function")
                print("1) Interpolation property check,")
                print("2) Partition of unity property check,")
                print("Shape function passed all the tests")
        else: 
            raise ValueError("Value of zeta should be within [-1,1] and i should be int")

    def test_dphi(self,T):
        self.update_dphi()
        ans = np.matmul(self.dphi.reshape(1,2),T.reshape(2,1))
        if ans != 0:
            print("test_dphi failed")
            print("Ensure that the temperature field was constant at both nodes")
            raise ValueError(f"Shape function evaluated to be {ans}")
        else:
            print("dphi passed the test.")

    def mesh_list(self):
        l = []
        for q in range(ip.n):
            l.append(q/(ip.n-1))
        if l[-1] != 1 or len(l) == 0:
            print("The mesh list is not generated correctly")           
        return l

    def update_phi(self,zeta):
        if zeta <=1 and zeta >= -1:
            self.phi=np.array([0.5*(1-zeta),0.5*(1+zeta)]).reshape(2,1) # The shape function is Linear.
        else: 
            raise ValueError("Value of zeta should be within [-1,1]")

    def update_Jacobi(self,i):
        x1 = self.nl[i]
        x2 = self.nl[i+1]
        self.J = 0.5*(x2-x1)
        if self.J < 0:
            raise ValueError("The determinant of Jacobi is negative")
        
    def update_dphi(self): # The shape function's first and second order partial deravatives.
        self.dphi = (self.J**-1)*np.array([-0.5,0.5]).reshape(2,1)
    
    def update_element(self,i):
        if  i<0 or type(i) != int:
            raise ValueError("Incorrect element number")
        self.update_Jacobi(i)
        self.update_dphi()
