import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
import units_cgs as cgs
# all in cgs

e_0 = cgs.m_n ** 4 * cgs.c ** 5 / (cgs.pi ** 2 * cgs.hbar ** 3)
p = lambda x: e_0/24.0 * ((2 * x ** 3 - 3 * x) * np.sqrt(1 + x ** 2) 
                                            + 3 * np.arcsinh(x))
e = lambda x: e_0/8.0 * ((2 * x ** 3 + x) * np.sqrt(1 + x ** 2) 
                                            - np.arcsinh(x))

p_p0 = lambda x: p(x) - 1.e33
root = newton(p_p0, 1.0, fprime = None, maxiter = 100000, tol = 1.e-10)

print(root)
print(e(root))
print(p(root))

# x0 = np.linspace(1.e-10, 0.5, num=int(1e3), endpoint=False)
# x1 = np.linspace(0.5, 1.0, num=int(1e4), endpoint=False)
# x2 = np.linspace(1.0, 5.0, num=int(1e4), endpoint=False)
# x3 = np.linspace(5.0, 10.0, num=int(1e4), endpoint=False)
# x4 = np.linspace(10.0, 25.0, num=int(1e4))

# x = np.concatenate((x0, x1))
# x = np.concatenate((x, x2))
# x = np.concatenate((x, x3))
# x = np.concatenate((x, x4))
x = np.logspace(np.log10(1.e-2), np.log10(25.0), num=int(1e5), 
                endpoint=True, dtype=np.float64)
p_arr = p(x)
e_arr = e(x)
print(np.all(np.diff(p_arr) > 0))
print(x)
# for i in range(len(p_arr)-1):
#     if p_arr[i] >= p_arr[i+1]:
#         print(i, 'not increasing index')
    

# plt.plot(p_arr, e_arr)
# plt.show()


p_arr *= cgs.erg_cm_to_MeV_fm
e_arr *= cgs.erg_cm_to_MeV_fm
data = np.column_stack((p_arr, e_arr))
np.savetxt(f'./eos/relativistic.dat', data, fmt='%.6E', delimiter=' ', 
                header='', comments='')

