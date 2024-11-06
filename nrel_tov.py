import numpy as np
import units_cgs as cgs
import time
import nrel_rk4 as rk4




start_time = time.time()
del_h = -1.e-4 # step size
EoS = rk4.e
step = 800
M = np.zeros(step)
R = np.zeros(step)
P = np.zeros(step)  
p_i = 0.0001 # units in MeV/fm^3
p_f = 10 # units in MeV/fm^3

P = np.logspace(np.log10(p_i), np.log10(p_f), num = step, endpoint = True)

for l in range(1, step + 1):
    
    allocate = 40000
    #cut = 0
        
    # P[l-1] = p_i + (l-1) * (p_f-p_i)/step
    M[l-1], R[l-1] = rk4.TOV(EoS, P[l-1], del_h, allocate)
    print('-------------------------------------------')
    print('Polytropic', l/step * 100, '%')
    print(f'Initial pressure: {P[l-1] / cgs.MeV_fm_to_km} MeV/fm^3')
    print(M[l-1], 'M_0')
    print(R[l-1], 'km')
    print('-------------------------------------------')


print(f'Initial Pressure {p_i}MeV/fm^3 ~ {p_f}MeV/fm^3')
end_time = time.time()
time_diff = end_time - start_time
    
data = np.column_stack((P, M, R))
    
np.savetxt('./result_0824/Polytropic_P_M_R.dat', data, fmt='%.6E', delimiter=' ', 
                header='', comments='')
    
    
print('Polytropic', 'took',time_diff, 'sec\n')
print(fr'Final: {R[-1]} km, {M[-1]}M_0, Initial pressure {P[-1]}')
print('Maximum Mass:', np.max(M), 'Initial pressure:', P[np.argmax(M)])
print('Maximum Radius', np.max(R), 'Initial pressure:', P[np.argmax(R)])