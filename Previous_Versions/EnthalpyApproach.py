import numpy as np
import Initial_Temperature_Dist as IT
import postProcessor as pp
import variableUpdater as vu

#The length of the material is assumed to be 1m
#The number of ELEMENTS decided are  10
#below defined are the global variables
#As of now I have taken the materials as Al

np.set_printoptions(precision=2, suppress=True)

delt = 4*10**-3
ratioI = 0.3333   #I/Iref  Iref/I = 3
n = 101               #number of nodes  
c = 9.0*10**2         #specific heat capacity           # in J/kgK 
rho = 2.7 * 10**3     #density of the material      #kg/m**3
Tliq = 933.3          #Kelvin 	  	  #Liq. Temp. [K] 	  	  #In the unit ??????????????????????
D = 0.7970           
Ta = -0.4             # ambient temperature in Kelvin.
#T = Ta*np.ones(11).reshape(n,1)    # The temperature array for the initial start. 	  	  #In the unit K.
lambf = 0.25 	  	  # The lambda constant 
np.set_printoptions(edgeitems=11, threshold=np.inf, linewidth=np.inf, suppress=True, formatter={'int': lambda x: str(x)})
# Hm is considered to be zero as it is dimensionless and at the start of melting the Hm = 0 

def mesh_list():
    global n
    l = []
    ele_list = {}
    for i in range(n):
        l.append(i/(n-1))           
    return l

T = IT.Initial_Temp(ratioI,Ta,2,mesh_list())\
    .Tdist_Enthalpy().copy().reshape(n,1) #The temperature array for the initial start.  #Dimensionless.
#print(f"The initial temprature distribution looks like\n {T} ") #1.14

H = vu.Update(D,lambf,T)          # The enthalpy was calculated as taking into account the above temperature array

def phi(zeta):
    shape_fun=np.array([0.5*(1-zeta),0.5*(1+zeta)]) # The shape function is Linear.
    return shape_fun.reshape(2,1)

def phi_derivative(i): # The shape function's first and second order partial deravatives.
    shape_fun = np.array([0.5,0.5])
    return (J(i)**-1)*shape_fun.reshape(2,1)

def J(i):
    node_list = mesh_list()
    x1 = node_list[i]
    x2 = node_list[i+1]
    return x2-x1

def Matrial_model(zeta,i,Te,Hk_1,Het_1):
    '''Physics subroutine or material subroutine'''

    global delt,ratioI,n,D,Ta,lambf

    if i == 0 and zeta<0: #This condition shows that we are on the first node
        b = ratioI*phi(zeta) - Ta*phi_derivative(i)[1]*phi_derivative(i) #Column vector
    else:
        b = - Ta*phi_derivative(i)[1]*phi_derivative(i)

    
    dT_dH = [[1/D if Hk_1[0] < 0 else 0 if Hk_1[0] ==  0 else 1/D, 0],
             [0, 1/D if Hk_1[1] < 0 else 0 if Hk_1[1] ==  0 else 1/D]]
    

    M = np.matmul(phi(zeta),phi(zeta).reshape(1,2))
    N = np.matmul(phi_derivative(i),phi_derivative(i).reshape(1,2))
    F = D*(b-np.matmul(N,Te))
    G = np.matmul(M,(Hk_1-Het_1)) - delt*F
    dG = M + delt*D*np.matmul(N,dT_dH)
    return G,dG

def Aly(p):               # here the n represents the total number of elements and p represents the current element
    A = np.zeros((n,2)) # Assembly matrix
    A[p,0] = 1
    A[p+1,1] = 1 
    return A

def GaussLoop(He_1,Te,i,Het_1):    
    Gauss_points = [-0.577350269189626, 0.577350269189626] # Gauss points for the integration
    Ge = np.zeros([2,1])
    dGe = np.zeros([2,2])
    for zeta in Gauss_points: # zeta is the Gauss coordinate value of the single element.
                ''' Material model below'''
                MaterialCall = Matrial_model(zeta,i,Te,He_1,Het_1)
                Gt = MaterialCall[0]
                dGt = MaterialCall[1]
                Ge = Ge + Gt
                dGe = dGe +dGt
    return Ge,dGe



def main():
    global H,delt,ratioI,n,D,Ta,lambf,T
    for t in range(1): 
        #print("The current time is {}".format(t))
        Hk= H[:,t][:,np.newaxis].copy() #Variable used for NRS scheme
        if not np.array_equal(Hk,H[:,-1][:,np.newaxis]): #Testing if the extraction was a success
            print("Enthalpy at previous time step failed to extract.")
            print("Program Terminated")
            break                        
        Ht = Hk.copy() # at the start we assume that the new and previous enthalpies are same
        Ht_1 = Hk.copy()
        Tnrs = T[:,t][:,np.newaxis].copy() #Tnrs is the netwon raphson scheme vector
        if not np.array_equal(Tnrs,T[:,-1][:,np.newaxis]): #Testing if the extraction was a success
            print("Temperature at previous time step failed to extract.")
            print("Program Terminated")
            break
        #Ht = new enthalpy after the newton Raphson scheme          #Ht_1 is the time-1 Enthalpy      #H_1e is the elemental new enthalpy     #He is the elemental old enthalpy 
        G = np.zeros([n,1])
        dG = np.zeros([n,n])
        while True:
            for i in range(n-1): #n is the number of elements in the mesh. i is the element number.
                A = Aly(i) 
                #He = np.matmul(A.T,Hk) # Assembly matrix for element i. He is a 2x1 vector.
                He_1 = np.matmul(A.T,Ht)
                Het_1 = np.matmul(A.T,Ht_1) 
                Te = np.matmul(A.T,Tnrs) #Elemental T extracted from universal T 
                [Ge,dGe] = GaussLoop(He_1,Te,i,Het_1)
                G = G + np.matmul(A,Ge)   #G is a column vector
                dG = dG + np.matmul(np.matmul(A,dGe),A.T)    #(n,2)*(2,2)*(2,n) 

            #Activating the boundry condition
            G[n-1] = 0
            dG[n-1] = 0
            dG[:,n-1] = 0
            dG[n-1,n-1] = 1     

            '''np.set_printoptions(edgeitems=11, threshold=np.inf, linewidth=np.inf, suppress=True, formatter={'int': lambda x: str(x)})
            print(np.array2string(dG,separator = ',')) ####DELETE THIS LINE'''
            Ht = Hk - np.matmul(np.linalg.inv(dG),G)
            #Condition = TCR.condition(Ht,Hk)
            diff = np.abs(Ht-Hk)
            if np.all(diff<np.exp(-5)):
                Tnrs = vu.Update(D,lambf,None,Hk)
                break #If the difference between the new and previous enthalpy is lower than tolerance, then the NRS Scheme is over.
            else:
                Hk = Ht.copy() #Otherwise, the NRS Scheme continues to optimize the Enthalpy.
                Tnrs = vu.Update(D,lambf,None,Hk)
        T = np.column_stack((T,Tnrs))
        H = np.column_stack((H,Ht))
    return Ht,T,t
try:
    [Ht,T,t] = main()
except ValueError:
    print("Matrix Multiplication Failed")
print("The enthalpy and Temperature distribution is shown below")
print("the first two columns correspond to the Enthalpy and the latter two ")
print("columns corresponds to Temperature")
print("[Ht Ht+1 Tt Tt+1]")
print(np.column_stack((H,T)))

pp.plotTemprature(T,t)
