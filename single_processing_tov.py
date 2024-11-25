import numpy as np
import time
import rk4_tov as rk4
import os
from datetime import datetime
from units_cgs import start_inner_crust
# All units in km ( c = G = 1 ) 


def main(del_h = -1.e-5):
    start_time = time.time()
    current_time = datetime.now()
    file_name = str(input('Select EoS in \'eos\' dir: '))
    file_name = './eos/' + file_name
    EoS = rk4.eos(file_name)
    
    print('------Input for start of central pressure in units MeV/fm^3------')
    while True:
        p_i = float(input('Initial pressure: '))
        if p_i < start_inner_crust:
            print(f'Initial pressure must be greater than {start_inner_crust:.3f} MeV/fm^3')
            continue
        
        else:
            p_f = float(input('Final pressure: '))
            break
        
    

    step = int(input(f'Number of steps between {p_i}~{p_f}: '))
    M = np.zeros(step)
    R = np.zeros(step)
    P_c = np.zeros(step)
    P_c = np.logspace(np.log10(p_i), np.log10(p_f), num = step, endpoint = True) # array of central pressure
    
    for l in range(1, step + 1):
        allocate = 1000000

        M[l-1], R[l-1], iterations, final_pressure, inner_crust, outer_crust, end = rk4.TOV(EoS, P_c[l-1], del_h, allocate)
        print('------------------------------------------------')
        print(file_name, l/step * 100, '%')
        print(f'Central pressure: {P_c[l-1]} MeV/fm^3')
        print(f'{M[l-1]:.5f} M_0')
        print(f'{R[l-1]:.5f} km')
        print(f'Final pressure: {final_pressure}MeV/fm^3')
        print(f'inner crust length {inner_crust:.5f} km')
        print(f'outer crust length {outer_crust:.5f} km')
        print(f'{iterations} iterations') # iterations taken
        print(end) # end method in rk4
        print('------------------------------------------------\n')


    print(f'Central Pressure {p_i} MeV/fm^3 ~ {p_f} MeV/fm^3')
    # make new directory in results using current date and time
    
    
    dir_name = f'results/{current_time.date()}'
    os.makedirs(dir_name, exist_ok=True)
    
    end_time = time.time()
    time_diff = end_time - start_time
    
    data = np.column_stack((P_c, M, R))
    time_of_file_name = str(current_time.time())[:8]
    np.savetxt(f'./{dir_name}/{file_name[5:]}_{time_of_file_name}.txt', data, fmt='%.6E', delimiter=' ', 
                    header='', comments='')
    
    # output of final, maximum results
    print(f'{file_name} took {time_diff:.4f} s')
    print(fr'Final: {R[-1]:.3f} km, {M[-1]:.3f} M_0, central pressure {P_c[-1]:.5f} MeV/fm^3')
    print(f'Maximum Mass: {np.max(M):.3f} M_0, central pressure: {P_c[np.argmax(M)]:.5f} MeV/fm^3')
    print(f'Maximum Radius {np.max(R):.3f} km, central pressure: {P_c[np.argmax(R)]:.5f} MeV/fm^3')
    
    rk4.plot(R, M, f'{dir_name}/plot_{time_of_file_name}')
    
    return 0


if __name__ == "__main__":
    main()