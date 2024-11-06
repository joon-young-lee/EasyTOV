import numpy as np
import matplotlib.pyplot as plt

file_name0 = './results/result_1025/relativistic.dat'
file_name1 = './results/result_1025/relativistic_nocrust.dat'

# file_name1 = f'Result_20240815/{file_name}' # without metric

#file_name2 = f'Result_metric/{file_name}' # metric h, 1.0MeV/fm^3 ~ 100.0MeV/fm^3


data0 = np.loadtxt(file_name0)
data1 = np.loadtxt(file_name1)


M0 = data0[:, 1]
R0 = data0[:, 2]
M1 = data1[:, 1]
R1 = data1[:, 2]

# plt.subplot(1, 2, 2)
plt.figure(figsize=(9, 7))
# plt.title('Radius vs. Mass', fontsize=20)
plt.plot(R0, M0, '-', markersize=2, color='r', label='EOS with crust')
plt.plot(R1, M1, '-', markersize=2, color='b', label='EOS without crust')
plt.xlabel('Radius in km', fontsize=20)
plt.ylabel(r'Mass in $M_\odot$', fontsize=22)
plt.xlim(4.0, 25.0)
plt.ylim(0.0, 0.8)
plt.grid()
fontsize = 15
# ----------------------------------------------------------------
# maximum initial pressure
#plt.text(8.0, 1.5, r'initial pressure $10^5 MeV/fm^3$', fontsize=fontsize)
# plt.arrow(9.6, 1.6, R[-1] - 9.6, M[-1] - 1.6)
# ----------------------------------------------------------------
# maximum mass
#plt.text(11.5, 2.0, r'initial pressure $667 MeV/fm^3$', fontsize=fontsize)
# plt.arrow(12.9, 2.1, R[np.argmax(M)] - 12.9, np.max(M) - 2.1)
# ----------------------------------------------------------------
# maximum radius
#plt.text(11.0, 0.5, r'initial pressure $0.08 MeV/fm^3$', fontsize=fontsize)
# plt.arrow(13.2, 0.45, np.max(R) - 13.2, M[np.argmax(R)] - 0.45)
# plt.quiver(13.2, 0.45, np.max(R) - 13.2, M[np.argmax(R)] - 0.45)
# ----------------------------------------------------------------
#x = np.array([9.6, 12.9, 13.2])
#y = np.array([1.6, 2.1, 0.45])
# X, Y = np.meshgrid(x, y)
#x_direction = np.array([9.8-9.6, R[np.argmax(M)] - 12.9, np.max(R) - 13.2])
#y_direction = np.array([2.01-1.6, np.max(M) - 2.1, M[np.argmax(R)] - 0.45])
#plt.quiver(x, y, x_direction, y_direction, scale=1, units='xy')

# 1.0 ~ 5000MeV/fm^3 central energy density
plt.text(11, 0.7, f'Maximum mass with crust {np.max(M0):.4f}\nMaximum mass witout crust {np.max(M1):.4f}', fontsize=17)
plt.text(13.7, 0.6, r'central density up to $ \sim 5000MeV/fm^3$', fontsize = 15)
plt.suptitle('Pure Fermi Gas core EoS', fontsize=20)
plt.savefig(f'./results/result_1025/plot.png')
plt.legend(fontsize=20)
plt.show()