import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton
import units_cgs as cgs


F = lambda u: np.sqrt(u)
deriv_F = lambda u: 1/(2.0 * np.sqrt(u))
p = lambda n: cgs.n_0 * (2./3. * cgs.E_0 * (n/cgs.n_0) ** (5/3)
                        + cgs.A * (n/cgs.n_0) ** 2 / 2.
                        + cgs.B * cgs.sigma * ((n/cgs.n_0) ** (1 + cgs.sigma)) / (1 + cgs.sigma))

# a = neutron to proton ratio
p_a = lambda a, n: p(n) + cgs.n_0 * a ** 2 * ((2. ** (2/3) - 1) * \
                cgs.E_0 * (2./3. * (n/cgs.n_0) ** (5./3.)
                - (n/cgs.n_0) ** 2 * deriv_F(n/cgs.n_0))
                + cgs.S_0 * deriv_F(n/cgs.n_0) * (n/cgs.n_0) ** 2)
# p_pure_neutron = lambda n: p_a(1.0, n) - 1.e32

e = lambda n: n * cgs.GeV_m_n + cgs.E_0 * n * (n/cgs.n_0) ** (2/3) \
    + cgs.A / 2 * n ** 2 / cgs.n_0 + n * cgs.B * ((n/cgs.n_0) ** cgs.sigma)/(1 + cgs.sigma)

e_a = lambda a, n: e(n) + n * a ** 2 * ((2 ** (2/3) - 1) * cgs.E_0 
                                     * ((n/cgs.n_0) ** (2/3) - F(n)) + cgs.S_0 * F(n))


# root = newton(p_pure_neutron, 1.e10, fprime = None, maxiter = 100000, tol = 1.e-15)

# print(root)
# print(p_pure_neutron(root))
# print(e_pure_neutron(root))
# print(root / cgs.fm ** 3)


# x = np.array(np.linspace(1.e-8 * cgs.n_0, 1.e-3 * cgs.n_0, num=int(1e4), endpoint=False))
# y = np.array(np.linspace(1.e-3 * cgs.n_0, 1.e1 * cgs.n_0, num=int(1e4), endpoint=False))
# z = np.array(np.linspace(1.e1 * cgs.n_0, 4 * 1.e1 * cgs.n_0, num=int(1e4), endpoint=True))
# # print(len(x))
# x = np.concatenate((x, y))
# x = np.concatenate((x, z))
x = np.logspace(np.log10(1.e-8 * cgs.n_0), np.log10(4 * 1.e1 * cgs.n_0),
                num=int(1e5), endpoint=True)
# 1.e-8n_0 ~ 40n_0


p_arr = p_a(1.0, x)
e_arr = e_a(1.0, x)
print(np.all(np.diff(p_arr) > 0))

data = np.column_stack((p_arr, e_arr))
np.savetxt(f'./eos/empirical3.dat', data, fmt='%.6E', delimiter=' ', 
                header='', comments='')

# convert to cgs
p_arr *= cgs.MeV_fm_to_cgs
e_arr *= cgs. MeV_fm_to_cgs
plt.figure(figsize=(10,8))

plt.plot(p_arr, e_arr)
plt.xlabel(r'Pressure in $dyne/cm^2$')
plt.ylabel(r'Energy density in $dyne/cm^2$')
# plt.xlim(1.e25, 1.e38)
# plt.ylim(1.e33, 1.e38)
plt.xscale('log')
plt.yscale('log')

plt.show()


