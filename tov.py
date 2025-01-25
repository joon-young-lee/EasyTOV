import os
import time
import numpy as np
import rk4_tov as rk4
from units_cgs import start_inner_crust
from datetime import datetime
from multiprocessing import Pool, cpu_count
# Uses multiprocessing
# Worker function to handle the computation for each pressure
def compute_star(args):
    EoS, P_c, del_h, allocate, step_index, total_steps, file_name = args
    M, R, iterations, final_pressure, inner_crust, outer_crust, end = rk4.TOV(EoS, P_c, del_h, allocate)
    
    print('------------------------------------------------')
    print(f'{file_name} {step_index/total_steps * 100:.2f} %')
    print(f'Central pressure: {P_c} MeV/fm^3')
    print(f'{M:.5f} M_0')
    print(f'{R:.5f} km')
    print(f'Final pressure: {final_pressure} MeV/fm^3')
    print(f'Inner crust length: {inner_crust:.5f} km')
    print(f'Outer crust length: {outer_crust:.5f} km')
    print(f'{iterations} iterations')  # iterations taken
    print(end)  # end method in rk4
    print('------------------------------------------------\n')
    
    return M, R

# Main function
def main(del_h=-1.e-5):
    start_time = time.time()
    current_time = datetime.now()
    print(f'Number of cores available: {cpu_count()}')
    file_name = str(input("Select EoS in 'eos' dir: "))
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
    P_c = np.logspace(np.log10(p_i), np.log10(p_f), num=step, endpoint=True)  # Array of central pressures

    # Multiprocessing setup
    allocate = 1000000
    pool_args = [(EoS, P_c[l], del_h, allocate, l + 1, step, file_name) for l in range(step)]
    with Pool(processes=cpu_count()) as pool:
        results = pool.map(compute_star, pool_args)

    # Collect results
    M, R = zip(*results)  # Separate masses and radii
    M = np.array(M)
    R = np.array(R)

    print(f'Central Pressure {p_i} MeV/fm^3 ~ {p_f} MeV/fm^3')
    
    # Create results directory
    dir_name = f'results/{current_time.date()}'
    os.makedirs(dir_name, exist_ok=True)
    
    end_time = time.time()
    time_diff = end_time - start_time
    
    # Save results to file
    data = np.column_stack((P_c, M, R))
    time_of_file_name = str(current_time.time())[:8]
    np.savetxt(f'./{dir_name}/{file_name[5:]}_{time_of_file_name}.txt', data, fmt='%.6E', delimiter=' ', 
               header='', comments='')

    # Output final and maximum results
    print(f'{file_name} took {time_diff:.4f} s')
    print(fr'Final: {R[-1]:.3f} km, {M[-1]:.3f} M_0, central pressure {P_c[-1]:.5f} MeV/fm^3')
    print(f'Maximum Mass: {np.max(M):.3f} M_0, central pressure: {P_c[np.argmax(M)]:.5f} MeV/fm^3')
    print(f'Maximum Radius {np.max(R):.3f} km, central pressure: {P_c[np.argmax(R)]:.5f} MeV/fm^3')
    
    rk4.plot(R, M, f'{dir_name}/plot_{time_of_file_name}')
    return 0

if __name__ == "__main__":
    main()
 