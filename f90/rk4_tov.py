import matplotlib.pyplot as plt
import numpy as np




def plot(R, M, file_name):
    plt.figure(figsize=(9, 7))
    plt.title('Radius vs. Mass', fontsize=20)
    plt.scatter(R, M, linewidths=0.5, color='r', label='EoS with crust')
    
    plt.xlabel('Radius in km', fontsize=20)
    plt.ylabel(r'Mass in $M_\odot$', fontsize=22)
    plt.xlim(4.0, 17.0)
    plt.ylim(0.0, 3.0)
    plt.grid()
    
    # plt.text(4, 1.5, f'Maximum mass with crust {np.max(M):.2f}', fontsize=17)
    
    plt.suptitle(file_name, fontsize=20)
    plt.savefig('plot')
    plt.legend(fontsize=20, loc='upper left')
    plt.show()

    return 0

data = np.loadtxt('output.txt', skiprows=1)

plot(data[:, 1], data[:, 0],"output.png")