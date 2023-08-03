import matplotlib.pyplot as plt

def plotTemprature(T,t):
    plt.plot(T)
    plt.xlabel('Node number')
    plt.ylabel('Dimensionless Temperature')
    plt.title(f'Temperature Vs Depth at time {t}')
    plt.show()