import numpy as np
import units_cgs as cgs
import time
import rk4_tov as rk4
# All units in km ( c = G = 1 ) 

del_h = -1.e-5
def main():
    start_time = time.time()
    # file_name = './eos/relativistic.dat'
    file_name = './eos/empirical.dat'
    # file_name = './eos/empirical2.dat'
    # file_name = './eos/empirical3.dat'
    # file_name = './eos/nrel.dat'
    EoS = rk4.eos(file_name)
    # step = 800
    print('------Input for start of central pressure in units MeV/fm^3------')
    p_i = float(input('Initial pressure: '))
    p_f = float(input('Final pressure: '))
    step = int(input(f'Number of steps between {p_i}~{p_f}: '))
    M = np.zeros(step)
    R = np.zeros(step)
    P_c = np.zeros(step)
    # p_i = 0.006 # units in MeV/fm^3
    # p_f = 6.0 # units in MeV/fm^3
    # Empirical limit 10^5MeV/fm^3
    P_c = np.logspace(np.log10(p_i), np.log10(p_f), num = step, endpoint = True) # central pressure
    # print(len(P))
    # print(P[0])
    # print(P[-1])
    for l in range(1, step + 1):
        allocate = 1000000

        M[l-1], R[l-1] = rk4.TOV(EoS, P_c[l-1], del_h, allocate)
        print('------------------------------------------------')
        print(file_name, l/step * 100, '%')
        print(f'Central pressure: {P_c[l-1]} MeV/fm^3')
        print(M[l-1], 'M_0')
        print(R[l-1], 'km')
        print('------------------------------------------------')


    print(f'Central Pressure {p_i}MeV/fm^3 ~ {p_f}MeV/fm^3')
    end_time = time.time()
    time_diff = end_time - start_time
        
    data = np.column_stack((P_c, M, R))
        
    np.savetxt(f'./results/result_1025/{file_name[5:]}', data, fmt='%.6E', delimiter=' ', 
                    header='', comments='')
    
    
    print(file_name, 'took',time_diff, 'sec\n')
    print(fr'Final: {R[-1]} km, {M[-1]}M_o, Initial pressure {P_c[-1]}MeV/fm^3')
    print(f'Maximum Mass: {np.max(M):.2f} M_o, Initial pressure: {P_c[np.argmax(M)]:.2f} MeV/fm^3')
    print(f'Maximum Radius {np.max(R):.2f} km, Initial pressure: {P_c[np.argmax(R)]} MeV/fm^3')
    return 0


if __name__ == "__main__":
    main()