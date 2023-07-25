'''Topic: PPP
Library: - Check_Inputs
#############################################################################################################################
Importing the required standard libraries 
-----------------------------------------------------------------------------------------------------------------------------
inputs -> Script where all the inputs are defined                                                                         '''
#############################################################################################################################
import inputs as ip
'''
#############################################################################################################################
Function check_inputs: -
=============================================================================================================================
Checks if the inputs from the inputs script is physically admissible.                                                                                                         
#############################################################################################################################
'''
def check_inputs(func):
    def Wrapper(self):
        if ip.soliduspt<0 or ip.liquidouspt <0 or ip.vapourpt<0 or ip.k<0 or ip.delt<0 or ip.ratioI<0 or ip.n < 0:
            raise ValueError("Incorrect material input arguments")
        elif ip.rho < 0 or ip.Lf < 0 or ip.D < 0 or ip.lambf < 0 or ip.alpha < 0:
            raise ValueError("Incorrect material input arguments")
        if (ip.n) < 2:
            raise ValueError("Mesh list cannot have just one node")
        return func(self)
    return Wrapper
