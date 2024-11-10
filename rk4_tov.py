import numpy as np
from scipy.interpolate import CubicSpline
import units_cgs as cgs

# Using Runge-Kutta 4th order
def eos(file_name):
    data = np.loadtxt(file_name)
    BPS = np.loadtxt('./crust/BPS.txt', skiprows=1)
    BPS_e = BPS[:, 0] * cgs.c2 * cgs.erg_cm_to_MeV_fm
    BPS_p = BPS[:, 1] * cgs.erg_cm_to_MeV_fm

    e_arr = data[:, 1] #* cgs.MeV_fm_to_km
    p_arr = data[:, 0] #* cgs.MeV_fm_to_km
    
    if np.all(np.diff(p_arr) > 0) == False:
        p_arr = np.flip(p_arr)
        e_arr = np.flip(e_arr)
    

    i = 0
    while i < len(e_arr):
        if e_arr[i] >= np.max(BPS_e):
            m = i
            #print(m)
            break
        elif e_arr[i] < np.max(BPS_e):
            i += 1
    
    
    p = np.concatenate((BPS_p, p_arr[m:]))
    e = np.concatenate((BPS_e, e_arr[m:]))
    
    # print(p)
    # print(e)
    
    p *= cgs.MeV_fm_to_km
    e *= cgs.MeV_fm_to_km
    

    eos = CubicSpline(p, e, bc_type='natural', extrapolate=False)
    
    
    return eos # eos in units km(c = G = 1)


def TOV(EoS, p0, del_h, num): # Initial pressure unit MeV/fm^3
    
    dp_dh = lambda p: EoS(p) + p # eos in units km(c = G = 1)
    dr2_dh = lambda r2, p, m: -2 * r2 * (np.sqrt(r2) - 2 * m)\
                            /(m + 4 * np.pi * p * r2 ** (3.0/2.0))
    dm_dh = lambda r2, p, m: -4 * np.pi * EoS(p) * r2 ** (3.0/2.0) * \
                            (np.sqrt(r2) - 2 * m)/(m + 4 * np.pi * p * r2 ** (3.0/2.0))
    
    p = p0 * cgs.MeV_fm_to_km
    
    r2 = -3 / (2 * np.pi * (EoS(p) + 3 * p)) * del_h
    m = 0.0
    # n_b = 0.0
    # Runge-Kutta 4th order
    outer = 1
    inner = 1
    for i in range(num):
            
        
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
        
        dp = (k1_p + 2 * k2_p + 2 * k3_p + k4_p) / 6 * del_h
        dm = (k1_m + 2 * k2_m + 2 * k3_m + k4_m) / 6 * del_h
        dr2 = (k1_r2 + 2 * k2_r2 + 2 * k3_r2 + k4_r2) / 6 * del_h
        
        if np.isnan(dp) == True or np.isnan(dm) == True:
            mass = m * cgs.km_to_M0
            radius = np.sqrt(r2)
            iterations = i
            final_pressure = p / cgs.MeV_fm_to_km
            inner_crust = np.sqrt(r2)-inner_start
            outer_crust = np.sqrt(r2)-outer_start
            end = 'Ended with NaN (dp or dm)'
            return mass, radius, iterations, final_pressure, inner_crust, outer_crust, end
        
        elif (np.abs(dm) / (m+1.e-15) < 1.e-15) and m != 0.0:
            mass = m * cgs.km_to_M0
            radius = np.sqrt(r2)
            iterations = i
            final_pressure = p / cgs.MeV_fm_to_km
            inner_crust = np.sqrt(r2)-inner_start
            outer_crust = np.sqrt(r2)-outer_start
            end = 'Ended with mass difference'

            return mass, radius, iterations, final_pressure, inner_crust, outer_crust, end
        # outer crust, inner crust
        
        elif (4.460e11 * cgs.c2 * cgs.erg_cm_to_MeV_fm * cgs.MeV_fm_to_km > EoS(p)) and outer == 1:
            outer = 0
            outer_start = np.sqrt(r2)
            
        
        elif (5.094e14 * cgs.c2 * cgs.erg_cm_to_MeV_fm * cgs.MeV_fm_to_km > EoS(p)) and (inner == 1):
            inner = 0
            inner_start = np.sqrt(r2)
        else:
            p += dp
            m += dm
            r2 += dr2


print(5.094e14 * cgs.c2 * cgs.erg_cm_to_MeV_fm)