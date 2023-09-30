from Mesh import Mesh as mesh
from Material_Subroutine import Crystal_model as Mat
import numpy as np
M = mesh()
M.update_Jacobi(0)
M.update_dphi()
#Material = Mat(M,5)
T = np.array([2,1]).reshape(-1,1)
H = np.array([1,2]).reshape(-1,1)

#Material.get_param(H,T,H,2)
k = np.sign(np.matmul(M.dphi.reshape(1,-1),T).item())
print(f"The slope output is {k}")
print(f"The value of the slopw is found to be {np.matmul(M.dphi.reshape(1,-1),T).item()}")
M.update_phi(0.5)
print(np.sum(M.phi))