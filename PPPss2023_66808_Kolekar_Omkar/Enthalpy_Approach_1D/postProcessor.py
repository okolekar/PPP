'''Topic: PPP
Library: - postProcessor
#############################################################################################################################
Importing the required standard libraries                                                                  
#############################################################################################################################
'''
import matplotlib
import matplotlib.pyplot as plt

def plotTemprature(x,T0,T1,t):
    '''
#############################################################################################################################
plotTemprature
-----------------------------------------------------------------------------------------------------------------------------
Takes Initial Temperature, Final Temperature and time and plots it                                                               
#############################################################################################################################
'''
    plt.plot(x,T0,label='Initial Temperature Distribution')
    plt.plot(x,T1, label='Final Temperature Distribution')
    plt.xlabel('Dimensionless Depth')
    plt.ylabel('Dimensionless Temperature')
    plt.title(f'Dimensionless Temperature Vs Dimensionless Depth at Dimensionless time {t}')
    plt.legend()
    plt.show()

def Simultaneous_Plotter(P,name1,Q,name2):
    '''
#############################################################################################################################
Simultaneous_Plotter
-----------------------------------------------------------------------------------------------------------------------------
Takes N number of Q and P inputs from a common family like temperature and plots it simultaneously
Note Common family means all vectors in P should belong to the temperature family or Enthalpy or etc                                             
#############################################################################################################################
'''
    matplotlib.use('TkAgg')
    for i in range(P.shape[1]):
        plt.plot(P[:,i], label=f't{i}')
    plt.xlabel('Node number')
    plt.ylabel(f'Dimensionless {name1}')
    plt.title(f'Variation of {name1} Vs Depth over time')
    plt.legend()
    plt.show(block=False)

    plt.figure()
    
    for i in range(Q.shape[1]):
        plt.plot(Q[:,i], label=f't{i}')
    plt.xlabel('Node number')
    plt.ylabel(f'Dimensionless {name2}')
    plt.title(f'Variation of {name2} Vs Depth over time isoparametric Scheme')
    plt.legend()
    plt.show()

def plotVerify(*args):
    '''
#############################################################################################################################
plotVerify
-----------------------------------------------------------------------------------------------------------------------------
Takes N number of Q and P inputs from a common family like temperature and plots it simultaneously
Plots the Temperature data for the verification case of the Enthalpy Approach                                        
#############################################################################################################################
'''
    for i in range(0,int(len(args)),2):
        plt.plot(args[i],label=args[i+1])
    #plt.xticks(range(len(args[0])))
    plt.xlabel('Node number')
    plt.ylabel('Dimensionless Temperature')
    plt.title('Temperature Verification')
    plt.legend()
    plt.show()

def plotxy(Y1,labelY1,Title,Xtitle):
    '''
#############################################################################################################################
plotxy
-----------------------------------------------------------------------------------------------------------------------------
This method can plot one vector of any type when the labels for the Y axis, Titel and the X axis                                         
#############################################################################################################################
'''
    plt.plot(Y1,label=labelY1)
    plt.xlabel(Xtitle)
    plt.ylabel(labelY1)
    plt.title(Title)
    plt.legend()
    plt.show()