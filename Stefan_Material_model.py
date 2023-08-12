#Stefan Model 2 noded element
#Stefan_Material_Model
import numpy as np
import Processed_Inputs as mm
class Stefan_Material_Model():
    def __init__(self) -> None:
        self.Ml = 0
        self.Nl = 0
        self.bl = 0
        self.Ms = 0
        self.Ns = 0
        self.bs = 0
#===========================================================================================================================#
                                        #Get liquid subdomain parameters
#===========================================================================================================================#
    def get_liq_param(self,mesh,ds_dt,ele_no):
        if ele_no == 0:
            print(type(mesh.J))
            self.Ml = np.array([[mesh.J,           0],
                                [0     ,    mesh.J*2]])
            
            self.Nl = np.array([[(1/(6*mm.n))*ds_dt+1/(mesh.J*2),   -(1/(6*mm.n))*ds_dt-1/(mesh.J*2)],
                                [              0                ,    (1/(3*mm.n))*ds_dt+2/(mesh.J*2)]])

            self.bl = np.array([[mm.ratioI],
                                [    0    ]])
        else:
            self.Ml = np.array([[mesh.J*2,          0],
                                [0         , mesh.J*2]])
            
            self.Nl = np.array([[(1/(3*mm.n))*ds_dt+2/(mesh.J*2)              ,   -((3*(ele_no+1)-2)/(6*mm.n))*ds_dt-1/(mesh.J*2)],
                                [((3*ele_no+2)/(6*mm.n))*ds_dt-(1/(mesh.J*2)) ,    (1/(3*mm.n))*ds_dt+2/(mesh.J*2)]])

            self.bl = np.array([[0],
                                [0]])

    def get_sol_param(self,mesh,ele_no,ds_dt):
        self.Ms = np.array([[mesh.J*2,          0],
                            [0         , mesh.J*2]])
        
        self.Ns = np.array([[-(ds_dt/(mm.m*3))+2/(2*mesh.J)                  ,   -((3*mm.m-3*(ele_no+1) + 2)*ds_dt/(6*mm.m))-(1/mesh.J*2)],
                            [(3*mm.m-3*ele_no-2)*ds_dt/(6*mm.m)-1/(mesh.J*2) ,                               -(ds_dt/(mm.m*3))+2/(2*mesh.J)]])
        
        if ele_no == 0:
            b = mm.Ta*((ds_dt/(3*mm.m))+1/(mesh.J*2))
            self.bs = np.array([[b],
                                [0]])
        else:
            b = np.array([[0],
                          [0]])