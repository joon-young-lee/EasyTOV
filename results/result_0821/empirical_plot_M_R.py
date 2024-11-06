import numpy as np
import matplotlib.pyplot as plt

# Result_1 - with same metric
# Result_metric - with new metric
# p_i = 1.602179e+33 # 1.0 MeV/fm^3
# p_f = 1.602179e+36 # 1000.0 MeV/fm^3
# step = 1000
file_name = './result_0821/empirical.dat'
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
plt.xlim(7.0, 14.0)
plt.ylim(0.0, 3.0)
plt.grid()
fontsize = 15
# ----------------------------------------------------------------
# maximum initial pressure
plt.text(8.0, 1.5, r'initial pressure $10^5 MeV/fm^3$', fontsize=fontsize)
# plt.arrow(9.6, 1.6, R[-1] - 9.6, M[-1] - 1.6)
# ----------------------------------------------------------------
# maximum mass
plt.text(11.5, 2.0, r'initial pressure $667 MeV/fm^3$', fontsize=fontsize)
# plt.arrow(12.9, 2.1, R[np.argmax(M)] - 12.9, np.max(M) - 2.1)
# ----------------------------------------------------------------
# maximum radius
plt.text(11.0, 0.5, r'initial pressure $0.08 MeV/fm^3$', fontsize=fontsize)
# plt.arrow(13.2, 0.45, np.max(R) - 13.2, M[np.argmax(R)] - 0.45)
# plt.quiver(13.2, 0.45, np.max(R) - 13.2, M[np.argmax(R)] - 0.45)
# ----------------------------------------------------------------
x = np.array([9.6, 12.9, 13.2])
y = np.array([1.6, 2.1, 0.45])
# X, Y = np.meshgrid(x, y)
x_direction = np.array([9.8-9.6, R[np.argmax(M)] - 12.9, np.max(R) - 13.2])
y_direction = np.array([2.01-1.6, np.max(M) - 2.1, M[np.argmax(R)] - 0.45])
plt.quiver(x, y, x_direction, y_direction, scale=1, units='xy')


plt.suptitle('Empirical EoS', fontsize=20)
plt.savefig(f'{file_name}_plot.png')
plt.show()