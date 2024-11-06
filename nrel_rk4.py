import numpy as np
from scipy.interpolate import CubicSpline
import units_cgs as cgs

# Using Runge-Kutta 4th order
K_nrel = cgs.hbar ** 2 / (15 * cgs.pi ** 2 * cgs.m_n) \
    * (3 * cgs.pi ** 2 / (cgs.m_n * cgs.c2)) ** (5/3)
print(K_nrel)

# p = lambda e: K_nrel * e ** (5/3) * cgs.erg_cm_to_MeV_fm
 
e = lambda p: (p / cgs.MeV_fm_to_km / cgs.erg_cm_to_MeV_fm /K_nrel) ** (3.0/5.0)\
            * cgs.erg_cm_to_MeV_fm * cgs.MeV_fm_to_km


def TOV(EoS, p0, del_h, allocate): # Initial pressure unit MeV/fm^3
    
    dp_dh = lambda p: EoS(p) + p # eos in units km(c = G = 1)
    dr2_dh = lambda r2, p, m: -2 * r2 * (np.sqrt(r2) - 2 * m)\
                            /(m + 4 * np.pi * p * r2 ** (3.0/2.0))
    dm_dh = lambda r2, p, m: -4 * np.pi * EoS(p) * r2 ** (3.0/2.0) * \
                            (np.sqrt(r2) - 2 * m)/(m + 4 * np.pi * p * r2 ** (3.0/2.0))
    
    # convert to km
    p = p0 * cgs.MeV_fm_to_km
    
    r2 = -3 / (2 * np.pi * (EoS(p) + 3 * p)) * del_h
    m = 1.e-20
    # Runge-Kutta 4th order
    # print(r2)
    pressure = np.zeros(allocate)
    mass = np.zeros(allocate)
    radius = np.zeros(allocate)
    for i in range(allocate):
            
        
        k1_p = dp_dh(p)
        k1_m = dm_dh(r2, p, m)
        k1_r2 = dr2_dh(r2, p, m)
        
        
        p_2 = p + del_h * k1_p / 2
        m_2 = m + del_h * k1_m / 2
        r_2 = r2 + del_h * k1_r2 / 2
        k2_p = dp_dh(p_2)
        k2_m = dm_dh(r_2, p_2, m_2)
        k2_r2 = dr2_dh(r_2, p_2, m_2)
        
        
        p_3 = p + del_h * k2_p / 2
        m_3 = m + del_h * k2_m / 2
        r_3 = r2 + del_h * k2_r2 / 2
        k3_p = dp_dh(p_3)
        k3_m = dm_dh(r_3, p_3, m_3)
        k3_r2 = dr2_dh(r_3, p_3, m_3)
        
        
        p_4 = p + del_h * k3_p
        m_4 = m + del_h * k3_m
        r_4 = r2 + del_h * k3_r2
        k4_p = dp_dh(p_4)
        k4_m = dm_dh(r_4, p_4, m_4)
        k4_r2 = dr2_dh(r_4, p_4, m_4)
        

        p += (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_h
        # if i  % 100 == 0:
        #     print(p, 'pressure')
        m += (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_h
        # print(m, 'mass')
        r2 += (k1_r2 + 2 * k2_r2 + 2 * k3_r2 + k4_r2) / 6 * del_h
        
        
        if p / cgs.MeV_fm_to_km  < 1.e-10 or np.isnan(p) == True \
            or np.isnan(m) == True \
            or np.isnan(r2) == True \
            or allocate - 1 <= i :
            print(i, 'iterations, ended with NaN')
            print(f'Final pressure: {pressure[i-1] / cgs.MeV_fm_to_km}')
            
            return mass[i-1] * cgs.km_to_M0, np.sqrt(radius[i-1])
        
        else:
            pressure[i] = p
            mass[i] = m
            radius[i] = r2