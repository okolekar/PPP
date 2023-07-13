import numpy as np

def check_arguments(func):
    def wrapper(self, *args,**kwargs):
        if len(args) == 5:
            return func(self,*args,**kwargs)
        elif len(args) == 2:
            return func(self,*args,**kwargs)
        else:
            raise ValueError(f"{len(args)} are incorrect number of input arguments")
    return wrapper

class temp_der():
    @check_arguments
    def __init__(self,*args):
        if args[0] < 0:
            raise ValueError("D cannot be less than zero")
        
        self.dT_dH = np.zeros([1,1])
        D = args[0]

        if len(args) == 5:
            Hk_1 = args[1]
            Tliq = args[2]
            if Tliq < 0:
                raise ValueError("Tliq cannot be less than zero")
            lambf = args[3]
            if lambf < 0:
                raise ValueError("lambf cannot be less than zero")
            Hliq = args[4]
            if Hliq < 0:
                raise ValueError("Hliq cannot be less than zero")
            self.Amorph_slop(D, Hk_1,Tliq,lambf,Hliq)

        if len(args) == 2:
            Hk_1 = args[1]
            self.Cryst_slop(D, Hk_1)


        
    def verify(self):
        if self.dT_dH[0][1] != 0 or self.dT_dH[1][0] != 0:
            raise ValueError("The 2*2 dT/dH Matrix is not constructed correctly")

    def Amorph_slop(self,D, Hk_1,Tliq,lambf,Hliq):
        try:
            self.dT_dH = [[1/D if Hk_1[0] <= 0 else Tliq/(D*Tliq + D*lambf) if np.logical_and(0 <= Hk_1[0], Hk_1[0] <= Hliq) else 1/D, 0],
                          [0,1/D if Hk_1[1] <= 0 else Tliq/(D*Tliq + D*lambf) if np.logical_and(0 <= Hk_1[1], Hk_1[1] <= Hliq) else 1/D]]
            self.verify()

        except ValueError as e:
            print(f"An error in dT_dH: {str(e)} following were the inputs")
            print(f"D = {D}") 
            print(f"Hk+1 = {Hk_1}")
            print(f"Tliq = {Tliq}")
            print(f"lambf = {lambf}")
            print(f"Hliq = {Hliq}")
            raise ValueError("Check the variables passed")   
  
    def Cryst_slop(self,D, Hk_1):
        try:
            self.dT_dH = [[1/D if Hk_1[0] < 0 else 0 if Hk_1[0] == 0 else 1/D, 0],
                        [0,1/D if Hk_1[1] < 0 else 0 if Hk_1[1] == 0 else 1/D]]
            self.verify()
            
        except ValueError as e:
            print(f"An error in dT_dH: {str(e)} following were the inputs")
            print(f"D = {D}") 
            print(f"Hk+1 = {Hk_1}")
            raise ValueError("Check the variables passed")
