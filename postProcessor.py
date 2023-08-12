import matplotlib
import matplotlib.pyplot as plt

def plotTemprature(T0,T1,t):
    plt.plot(T0,label='Initial Temperature Distribution')
    plt.plot(T1, label='Final Temperature Distribution')
    plt.xticks(range(len(T0)))
    plt.xlabel('Node number')
    plt.ylabel('Dimensionless Temperature')
    plt.title(f'Temperature Vs Depth at time {t}')
    plt.legend()
    plt.show()

def Simultaneous_Plotter(X,name1,Y,name2):
    matplotlib.use('TkAgg')
    for i in range(X.shape[1]):
        plt.plot(X[:,i], label=f't{i}')
    plt.xlabel('Node number')
    plt.ylabel(f'Dimensionless {name1}')
    plt.title(f'Variation of {name1} Vs Depth over time')
    plt.legend()
    plt.show(block=False)

    plt.figure()

    for i in range(Y.shape[1]):
        plt.plot(Y[:,i], label=f't{i}')
    plt.xlabel('Node number')
    plt.ylabel(f'Dimensionless {name2}')
    plt.title(f'Variation of {name2} Vs Depth over time isoparametric Scheme')
    plt.legend()
    plt.show()

def plotVerify(*args):
    for i in range(0,int(len(args)),2):
        plt.plot(args[i],label=args[i+1])
    plt.xticks(range(len(args[0])))
    plt.xlabel('Node number')
    plt.ylabel('Dimensionless Temperature')
    plt.title('Temperature Verification')
    plt.legend()
    plt.show()
