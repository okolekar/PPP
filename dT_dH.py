import numpy as np

def slope(D, Hk_1,Tliq,lambf,Hliq):
    try:
        temporary = [[1/D if Hk_1[0] < 0 else Tliq/(D*Tliq + D*lambf) if np.logical_and(0 <= Hk_1[0], Hk_1[0] < Hliq) else 1/D, 0],
                    [0,1/D if Hk_1[1] < 0 else Tliq/(D*Tliq + D*lambf) if np.logical_and(0 <= Hk_1[1], Hk_1[1] < Hliq) else 1/D]]
        return temporary
    except ValueError as e:
        print(f"An error in dT_dH: {str(e)} following were the inputs")
        print(f"D = {D}") 
        print(f"Hk+1 = {Hk_1}")
        print(f"Tliq = {Tliq}")
        print(f"lambf = {lambf}")
        print(f"Hliq = {Hliq}")   
    