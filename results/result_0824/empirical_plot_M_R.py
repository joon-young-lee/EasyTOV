import numpy as np
import matplotlib.pyplot as plt

# Result_1 - with same metric
# Result_metric - with new metric
# p_i = 1.602179e+33 # 1.0 MeV/fm^3
# p_f = 1.602179e+36 # 1000.0 MeV/fm^3
# step = 1000
file_name = './result_0824/empirical.dat'
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
plt.plot(R, M, '-', linewidth=2.5, color='r')
plt.xlabel('Radius in km', fontsize=20)
plt.ylabel(r'Mass in $M_\odot$', fontsize=22)
plt.xlim(7.0, 14.0)
plt.ylim(0.0, 3.0)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.grid()
fontsize = 20
# ----------------------------------------------------------------
# maximum initial pressure
plt.text(8.0, 1.5, r'initial pressure $10^5 MeV/fm^3$', fontsize=fontsize)
plt.arrow(9.6, 1.6, R[-1] - 9.6, M[-1] - 1.6)
# ----------------------------------------------------------------
# maximum mass
plt.text(10.9, 2.0, r'initial pressure $665 MeV/fm^3$', fontsize=fontsize)
plt.arrow(12.9, 2.1, R[np.argmax(M)] - 12.9, np.max(M) - 2.1)
# ----------------------------------------------------------------
# maximum radius
plt.text(10.7, 0.5, r'initial pressure $0.190 MeV/fm^3$', fontsize=fontsize)
plt.arrow(13.2, 0.45, 1.409216E+01 - 13.2, 3.149288E-02 - 0.45)
# ----------------------------------------------------------------


plt.text(7.0, 2.5, rf'Maximum mass: {np.max(M):.3f}$M_\odot$'+'\n'+f'Radius: {R[np.argmax(M)]:.2f}km', fontsize=20)
plt.suptitle('Empirical EoS', fontsize=25)
plt.savefig(f'{file_name}_plot.png')
plt.show()