import numpy as np
import Crystalline_inputs as ip
import matplotlib.pyplot as plt 
from Element2D import Element_Subroutine
from Gridify import Mesh
class Stefan2D():
    def __init__(self):
        self.X      = []
        self.Y      = []
        self.T      = np.zeros(50*50).reshape(-1,1)

    def solver(self):
        mesh = Mesh(50,50,ip.x,ip.y)
        for node in mesh.boundary_nodes['Bottom_Nodes']:
            x = mesh.xy_list[node][0]
            self.T[node] = -0.8480+0.9119*x**2-0.6079*x**3
        E = Element_Subroutine(50*50)
        E.formGlobalParm(mesh)
        A = E.K + (1/ip.delt)*E.N
        B = np.matmul((1/ip.delt)*E.N,self.T)
        self.T = np.linalg.solve(A,B)
        print(f"Temperature at firt node is {self.T[0]}")
        print(f"The temp at last node is {self.T[mesh.boundary_nodes['Bottom_Nodes'][-1]]}")
        plt.plot(self.T[mesh.boundary_nodes['Bottom_Nodes']],'-o')
        plt.xlabel("Dimensionless Distance")
        plt.ylabel("Dimensionless Temperature")
        plt.title("Temperature vs distance in Stefan Problem Approach 2D")
        plt.show()
        Vx = np.matmul(E.Vx,self.T)
        Vy = np.matmul(E.Vy,self.T)
        print("The velocity in the y direction is")
        print(Vy*ip.delt)


S = Stefan2D()
S.solver()