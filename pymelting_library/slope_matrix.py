import numpy as np
import inputs as ip

def check_arguments(func):
    def wrapper(self, *args,**kwargs):
        if len(args) == 2:
            return func(self,*args,**kwargs)
        else:
            raise ValueError(f"{len(args)} are incorrect number of input arguments")
    return wrapper

class temp_der():
    @check_arguments
    def __init__(self,*args):
        
        self.dT_dH = np.zeros([1,1])

        if args[0] == 1:
            Hk_1 = args[1]
            self.Amorph_slop(Hk_1)

        if args[0] == 2:
            Hk_1 = args[1]
            self.Cryst_slop(Hk_1)


        
    def verify(self):
        if self.dT_dH[0][1] != 0 or self.dT_dH[1][0] != 0:
            raise ValueError("The 2*2 dT/dH Matrix is not constructed correctly")

    def Amorph_slop(self, Hk_1):
        try:
            self.dT_dH = [[1/ip.D if Hk_1[0] <= 0 else ip.Tliq/(ip.D*ip.Tliq + ip.D*ip.lambf) if np.logical_and(0 <= Hk_1[0], Hk_1[0] <= ip.Hliq) else 1/ip.D, 0],
                          [0,1/ip.D if Hk_1[1] <= 0 else ip.Tliq/(ip.D*ip.Tliq + ip.D*ip.lambf) if np.logical_and(0 <= Hk_1[1], Hk_1[1] <= ip.Hliq) else 1/ip.D]]
            self.verify()

        except ValueError as e:
            print(f"An error in dT_dH: {str(e)} following were the inputs")
            print(f"Hk+1 = {Hk_1}")
            raise ValueError("Check the variables passed")   
  
    def Cryst_slop(self, Hk_1):
        try:
            self.dT_dH = [[1/ip.D if Hk_1[0] < 0 else 0 if Hk_1[0] == 0 else 1/ip.D, 0],
                        [0,1/ip.D if Hk_1[1] < 0 else 0 if Hk_1[1] == 0 else 1/ip.D]]
            self.verify()
            
        except ValueError as e:
            print(f"An error in dT_dH: {str(e)} following were the inputs")
            print(f"Hk+1 = {Hk_1}")
            raise ValueError("Check the variables passed")
