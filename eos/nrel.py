import numpy as np
import units_cgs as cgs

erg_cm_to_MeV_fm = 6.2414999e-34

K_nrel = cgs.hbar ** 2 / (15 * cgs.pi ** 2 * cgs.m_n) \
    * (3 * cgs.pi ** 2 / (cgs.m_n * cgs.c2)) ** (5/3)

e = lambda p: (p / K_nrel) ** (3.0/5.0)

x = np.logspace(np.log10(1.e24), np.log10(1.e34), num=int(1e6), endpoint=True)

# EoS 6.2415e-10MeV/fm^3 to 6.2415MeV/fm^3 pressure
p_arr = x * erg_cm_to_MeV_fm
e_arr = e(x) * erg_cm_to_MeV_fm

data = np.column_stack((p_arr, e_arr))
np.savetxt(f'./eos/nrel.dat', data, fmt='%.6E', delimiter=' ', 
                header='', comments='')
