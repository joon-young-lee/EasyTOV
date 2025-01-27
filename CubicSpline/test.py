import numpy as np
from scipy.interpolate import CubicSpline
import units_cgs as cgs

def eos(file_name):
    data = np.loadtxt(file_name)
    BPS = np.loadtxt('BPS.txt', skiprows=1)
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
            print(m)
            break
        elif e_arr[i] < np.max(BPS_e):
            i += 1
    
    
    p = np.concatenate((BPS_p, p_arr[m:]))
    e = np.concatenate((BPS_e, e_arr[m:]))
    
    
    p *= cgs.MeV_fm_to_km
    e *= cgs.MeV_fm_to_km
    

    eos = CubicSpline(p, e, bc_type='natural', extrapolate=False)
    print(len(p))
    return eos # eos in units km(c = G = 1)


EoS = eos('eft_pnm32_000001_ldm_eos_s.dat')

print(EoS(100 * cgs.MeV_fm_to_km))