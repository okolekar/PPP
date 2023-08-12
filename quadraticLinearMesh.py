#3noded_Mesh
import numpy as np
import Processed_Inputs as ip
from Mesh_List_Generator import gen_Meshlist
class Mesh():
    def __init__(self,s):
        self.phi        = 0
        self.dphi       = 0
        self.J          = 0
        self.node_list  = []
        self.generate_node_list(s)

    def test_phi(self):
        if np.sum(self.phi) != 1:
            raise ValueError("The sum of the derivate of the shape function is not zero")
        else:
            print("The shape function passed the test")

    def test_dphi(self):
        if np.sum(self.dphi) != 0:
            raise ValueError("The sum of the derivate of the shape function is not zero")
        else:
            print("Derivative of the shape function passed the test")

    def generate_node_list(self,s):
        generate = gen_Meshlist(ip.n,ip.m,s,ip.length) 
        l1 = next(generate)
        l2 = next(generate)
        self.node_list = [l1,l2]

    def update_phi(self,xi):
        self.phi = np.array([[-(1-xi)*xi*0.5],
                             [   1-xi**2    ],
                             [ (1+xi)*xi*0.5]])
        
    def update_dphi(self,xi):
        self.dphi = np.array([[-(1-2*xi)*0.5],
                              [    -2*xi    ],
                              [ (1+2*xi)*0.5]])
        
    def update_Jacobian(self,elem_no):
        if elem_no < ip.n:
            self.J = np.matmul(self.dphi.reshape(1,3),np.array(list(self.node_list[0][elem_no].values())).reshape(-1,1))
        else:
            self.J = np.matmul(self.dphi.reshape(1,3),np.array(list(self.node_list[1][elem_no].values())).reshape(-1,1))

    def update_param(self,xi,elem_no):
        self.update_phi(xi)
        self.update_dphi(xi)
        self.update_Jacobian(elem_no)