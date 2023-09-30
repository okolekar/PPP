import numpy as np
from slope_matrix import temp_der
import Processed_Inputs as ip
from Mesh import Mesh as mesh
from Material_Subroutine import Crystal_model as Mat

def test_temp_derivative():
    Hk = np.array([-1,ip.lambf]).reshape(-1,1)
    dT_dH = temp_der(1,Hk,3).dT_dH
    signs = np.sign(dT_dH)
    expected = [1,1]
    for i in range(2):
        assert signs[i,i] == expected[i]

def test_temp_derivative2():
    Hk = np.array([0,ip.lambf+0.5]).reshape(-1,1)
    dT_dH = temp_der(1,Hk,3).dT_dH
    signs = np.sign(dT_dH)
    expected = [0,1]
    for i in range(2):
        assert signs[i,i] == expected[i]

def test_temp_derivative3():
    Hk = np.array([ip.lambf,-1]).reshape(-1,1)
    dT_dH = temp_der(1,Hk,3).dT_dH
    signs = np.sign(dT_dH)
    expected = [1,1]
    for i in range(2):
        assert signs[i,i] == expected[i]

def test_temp_derivative4():
    Hk = np.array([0,ip.lambf-0.005]).reshape(-1,1)
    dT_dH = temp_der(1,Hk,3).dT_dH
    signs = np.sign(dT_dH)
    expected = [0,0]
    for i in range(2):
        assert signs[i,i] == expected[i]

def test_dphi():
    M = mesh()
    M.update_Jacobi(2)
    M.update_dphi()
    T = np.array([2,1]).reshape(-1,1)
    H = np.array([1,0]).reshape(-1,1)
    k = np.sign(np.matmul(M.dphi.reshape(1,-1),T).item())
    assert k == -1
    T = np.array([1,2]).reshape(-1,1)
    k = np.sign(np.matmul(M.dphi.reshape(1,-1),T).item())
    assert k == 1
    T = np.array([2,2]).reshape(-1,1)
    k = np.sign(np.matmul(M.dphi.reshape(1,-1),T).item())
    assert k == 0
    M.update_phi(0.5)
    assert np.sum(M.phi) == 1
