import numpy as np
import matplotlib.pyplot as plt

# Result_1 - with same metric
# Result_metric - with new metric
# p_i = 1.602179e+33 # 1.0 MeV/fm^3
# p_f = 1.602179e+36 # 1000.0 MeV/fm^3
# step = 1000
file_name = 'eft_pnm32_000001_ldm_result.dat'
file_name1 = f'Result_20240815/{file_name}' # without metric
#file_name2 = f'Result_metric/{file_name}' # metric h, 1.0MeV/fm^3 ~ 100.0MeV/fm^3


data1 = np.loadtxt(file_name1)
#data2 = np.loadtxt(file_name2)


M_1 = data1[:, 1]
R_1 = data1[:, 2]
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
plt.title('Radius vs. Mass, without metric change', fontsize=20)
plt.plot(R_1, M_1, '.', markersize=2, color='r')
plt.xlabel('Radius in km', fontsize=20)
plt.ylabel(r'Mass in $M_0$', fontsize=20)
plt.xlim(7.0, 14.0)
plt.ylim(0, 3.0)
plt.grid()

plt.suptitle(f'{file_name}', fontsize=20)
plt.savefig(f'{file_name}_plot.png')
plt.show()