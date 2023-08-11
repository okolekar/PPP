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
                                        #Gauss Summation
#===========================================================================================================================#
    def Gauss_Sum(self):
        q = [-0.774596669241483,    0,   0.774596669241483]
        w = [0.666666666665,        1,  0.6666666666666665]

#===========================================================================================================================#
                                        #Get liquid subdomain parameters
#===========================================================================================================================#
    def get_liq_param(self,mesh,ds_dt,ele_no):
        if ele_no == 0:

            global_coordinate = list(mesh.node_list[0][ele_no].values())
            he = global_coordinate[2]-global_coordinate[0]

            self.Ml = np.array([[0.5*he,       0     ,    0 ],
                                [0     ,       he    ,    0 ],
                                [0     ,       0     ,    he]])
            
            self.Nl = np.array([[(1/(6*mm.n))*ds_dt+1/(he) , -(1/(6*mm.n))*ds_dt-1/(he),  0],
                                [              0           ,  (1/(3*mm.n))*ds_dt+2/(he),  0],
                                []])

            self.bl = np.array([[mm.ratioI],
                                [    0    ],
                                [    0    ]])
        else:

            global_coordinate = list(mesh.node_list[0][ele_no].values())
            he = global_coordinate[2]-global_coordinate[0]

            self.Ml = np.array([[he  ,     0   ,   0 ],
                                [0   ,    he   ,   0 ],
                                [0   ,     0   ,   he]])
            
            self.bl = np.array([[0],
                                [0],
                                [0]])

    def get_sol_param(self,mesh):
        self.Ms = np.array([[mesh.J*2,          0],
                            [0         , mesh.J*2]])