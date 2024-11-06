import numpy as np
import matplotlib.pyplot as plt
import units_cgs as cgs

t0 = 1024.1 # MeV fm^3
t3 = 14600.8 # MeV fm^6

e_n = lambda n: cgs.MeV_m_n * n + 3.0 / (10.0 * cgs.MeV_m_n) \
                * (3 * cgs.pi ** 2 * cgs.MeV_hbar ** 3) ** (2.0/3.0) * n ** (5.0/3.0) \
                + t3 / 24.0 * n ** 3 - t0 / 4 * n ** 2

p_n = lambda n: 2 / (10 * cgs.MeV_m_n) * (3 * cgs.pi ** 2 * cgs.MeV_hbar ** 3) ** (2.0/3.0) * n ** (5.0/3.0) \
                + t3 / 12 * n ** 3 - t0 / 4 * n ** 2

x = np.logspace(np.log10(1.e-3 * cgs.n_0), np.log10(4 * 1.e1 * cgs.n_0),
                num=int(1e5), endpoint=True)

print(e_n(0.16))


p_arr = p_n(x)
e_arr = e_n(x)
data = np.column_stack((p_arr, e_arr))
np.savetxt(f'./eos/skyrme/eos.dat', data, fmt='%.6E', delimiter=' ', 
                header='', comments='')
p_arr *= cgs.MeV_fm_to_cgs
e_arr *= cgs.MeV_fm_to_cgs
plt.figure(figsize=(10,8))

plt.plot(p_arr, e_arr)
plt.xlabel(r'Pressure in $dyne/cm^2$')
plt.ylabel(r'Energy density in $erg/cm^3$')
plt.xlim(1.e30, 1.e38)
plt.ylim(1.e33, 1.e38)
plt.xscale('log')
plt.yscale('log')

plt.show()
