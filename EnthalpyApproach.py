import numpy as np

def N(zeta):
    shape_fun=np.array([0.5*(1-zeta),0.5*(1+zeta)]) # The shape function is a parabola.
    return shape_fun

def N_derivative(): # The shape function's first and second order partial deravatives.
    shape_fun = (0.5,0.5)
    return shape_fun

def M(zeta,ratioI,Ta): # The matrix element.
    b = ratioI*N(zeta)[2] - Ta*N_derivative()[2]*N_derivative()[2] # The term in the numerator.
    return np.array([[N(zeta)[1]*N(zeta)[1],N(zeta)[1]*N(zeta)[2]],
                     [N(zeta)[2]*N(zeta)[1],N(zeta)[2]*N(zeta)[2]]])