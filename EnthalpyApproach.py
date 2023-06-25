import numpy as np

#The length of the material is assumed to be 1m
#The number of ELEMENTS decided are  10
#below defined are the global variables
#As of now I have taken the materials as Al

np.set_printoptions(precision=2, suppress=True)

delt = 0.01
ratioI = 0.0123
n = 11                #number of nodes  
c = 9.0*10**2         #specific heat capacity           # in J/kgK 
rho = 2.7 * 10**3     #density of the material      #kg/m**3
Tliq = 933.3          #Kelvin 	  	  #Liq. Temp. [K] 	  	  #In the unit ??????????????????????
D = 139.64            
Ta = -0.4             # ambient temperature in Kelvin.
T = Ta*np.ones(11).reshape(11,1)    # The temperature array for the initial start. 	  	  #In the unit K.
lambf = 0.25 	  	  # The lambda constant 
Hk = rho*c*T          # The enthalpy was calculated as taking into account the below temperature array


def phi(zeta):
    shape_fun=np.array([0.5*(1-zeta),0.5*(1+zeta)]) # The shape function is Linear.
    return shape_fun.reshape(2,1)

def phi_derivative(): # The shape function's first and second order partial deravatives.
    shape_fun = np.array([0.5,0.5])
    return (J()**-1)*shape_fun.reshape(2,1)

def J():
    x1 = mesh_list()[0]
    x2 = mesh_list()[1]
    return x2-x1

def mesh_list():
    l = []
    ele_list = {}
    for i in range(10):
        l.append(i/10)        
    return l

def Matrial_model(zeta,i,Te,Hk_1,Hk):
    '''Physics subroutine or material subroutine'''

    global delt,ratioI,n,c,rho,D,Ta,lambf,Tliq

    if i == 0:
        b = ratioI*phi(zeta) - Ta*phi_derivative()[1]*phi_derivative() #Column vector
    else:
        b = - Ta*phi_derivative()[1]*phi_derivative()

    if Te.all()<0:    
        dT_dH = D
    elif Te.all() <= 0:
        dT_dH = D+D*lambf/Tliq
    elif Te.all() >= 0: 
        dT_dH = D 
    

    M = np.matmul(phi(zeta),phi(zeta).reshape(1,2))
    N = np.matmul(phi_derivative(),phi_derivative().reshape(1,2))
    F = D*(b-np.matmul(N,Te))
    G = np.matmul(M,(Hk_1-Hk)) - delt*F
    dG = M + delt*D*dT_dH*N
    return G,dG

def Aly(p):               # here the n represents the total number of elements and p represents the current element
    A = np.zeros((n,2)) # Assembly matrix
    A[p,0] = 1
    A[p+1,1] = 1 
    return A

def GaussLoop(He_1,He,Te,i):    
    Gauss_points = [-0.577350269189626, 0.577350269189626] # Gauss points for the integration
    Ge = np.zeros([2,1])
    dGe = np.zeros([2,2])
    for zeta in Gauss_points: # zeta is the Gauss coordinate value of the single element.
                ''' Material model below'''
                MaterialCall = Matrial_model(zeta,i,Te,He_1,He)
                Gt = MaterialCall[0]
                dGt = MaterialCall[1]
                Ge = Ge + Gt
                dGe = dGe +dGt
    return Ge,dGe

def main():
    global Hk,delt,ratioI,n,c,rho,D,Ta,lambf,T,Tliq
    for t in np.arange(0,1,delt):
         Ht= Hk.copy() # at the start we assume that the new and previous enthalpies are same
         #Ht = new enthalpy after the newton Raphson scheme          #Hk is the previous Enthalpy      #H_1e is the elemental new enthalpy     #He is the elemental old enthalpy 
         G = np.zeros([n,1])
         dG = np.zeros([n,n])
         while True:
            for i in range(n-1): #n is the number of elements in the mesh. i is the element number.
                A = Aly(i) 
                He = np.matmul(A.T,Hk) # Assembly matrix for element i. He is a 2x1 vector.
                He_1 = np.matmul(A.T,Ht)
                Te = np.matmul(A.T,T) #Elemental T extracted from universal T 
                [Ge,dGe] = GaussLoop(He_1,He,Te,i)
                G = G + np.matmul(A,Ge)   #G is a column vector
                dG = dG + np.matmul(np.matmul(A,dGe),A.T)    #(n,2)*(2,2)*(2,n) 

            #Activating the boundry condition
            G[n-1] = 0
            dG[n-1] = 0
            dG[:,n-1] = 0    

            np.set_printoptions(edgeitems=11, threshold=np.inf, linewidth=np.inf, suppress=True, formatter={'int': lambda x: str(x)})
            print(np.array2string(dG,separator = ',')) ####DELETE THIS LINE
                
            Ht = Hk - np.matmul(np.linalg.inv(dG),G)
            if (abs(Ht - Hk)<np.exp(-5)):
                break #If the difference between the new and previous enthalpy is lower than tolerance, then the NRS Scheme is over.
            else:
                Hk = Ht.copy() #Otherwise, the NRS Scheme continues to optimize the Enthalpy.
                T = Hk/(rho*c)
    return Ht,T

[Ht,T] = main()

print(Ht)
print("Printing T")
print(T)
