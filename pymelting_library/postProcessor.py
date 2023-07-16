import matplotlib.pyplot as plt
import numpy as np

def plotTemprature(T0,T1,t):
    plt.plot(T0,label='Initial Temperature Distribution')
    plt.plot(T1, label='Final Temperature Distribution')
    plt.xticks(range(len(T0)))
    plt.xlabel('Node number')
    plt.ylabel('Dimensionless Temperature')
    plt.title(f'Temperature Vs Depth at time {t}')
    plt.legend()
    plt.show()