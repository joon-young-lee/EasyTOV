import numpy as np
import matplotlib.pyplot as plt

# Result_1 - with same metric
# Result_metric - with new metric
# p_i = 1.602179e+33 # 1.0 MeV/fm^3
# p_f = 1.602179e+36 # 1000.0 MeV/fm^3
# step = 1000
file_name = './result_0821/relativistic.dat'
# file_name1 = f'Result_20240815/{file_name}' # without metric

#file_name2 = f'Result_metric/{file_name}' # metric h, 1.0MeV/fm^3 ~ 100.0MeV/fm^3


data1 = np.loadtxt(file_name)
#data2 = np.loadtxt(file_name2)


M = data1[:, 1]
R = data1[:, 2]
#M = data2[:, 1]
#R = data2[:, 2]




# plt.figure(figsize=(25, 8))
# plt.subplot(1, 2, 1)
# plt.title('Radius vs. Mass, new metric h', fontsize=20)
# plt.plot(R, M, '.', markersize=2, color='b')
# plt.xlabel('Radius in km', fontsize=20)
# plt.ylabel(r'Mass in $M_0$', fontsize=20)
# plt.xlim(8.0, 30.0)
# plt.ylim(0, 2.5)
# plt.grid()

# plt.subplot(1, 2, 2)
plt.figure(figsize=(10, 8))
plt.title('Radius vs. Mass', fontsize=20)
plt.plot(R, M, '.', markersize=2, color='r')
plt.xlabel('Radius in km', fontsize=20)
plt.ylabel(r'Mass in $M_\odot$', fontsize=22)
plt.xlim(3.0, 26.0)
plt.ylim(0, 0.8)
plt.grid()
fontsize = 15
# ----------------------------------------------------------------
# maximum initial pressure
plt.text(4.9, 0.3, r'initial pressure $10^8 MeV/fm^3$', fontsize=fontsize)
plt.arrow(7.0, 0.32, R[-1] - 7.0, M[-1] - 0.32)
# ----------------------------------------------------------------
# maximum mass
plt.text(7.5, 0.6, r'initial pressure $225 MeV/fm^3$', fontsize=fontsize)
plt.arrow(9.5, 0.622, R[np.argmax(M)] - 9.5, np.max(M) - 0.622)
# ----------------------------------------------------------------
# maximum radius
plt.text(17.5, 0.1, r'initial pressure $0.08 MeV/fm^3$', fontsize=fontsize)
plt.arrow(23.3, 0.13, np.max(R) - 23.3, M[np.argmax(R)] - 0.13)
# ----------------------------------------------------------------
plt.text(16, 0.7, rf'Maximum mass: {np.max(M):.3f}$M_\odot$'+'\n'+f'Radius: {R[np.argmax(M)]:.3f}km', fontsize=17)
plt.suptitle('Pure neutron, relativistic case', fontsize=20)
plt.savefig('relativistic_pure.png')
plt.show()