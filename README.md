Easy Tolman-Oppenheimer-Volkoff equation solver using RK4 method.

How to run:
  1. Make a EoS data sheet inside eos directory.
  2. $ python -u tov.py
  3. Put in the inputs as desried.

Inputs:
  1. Name of EoS data file in 'eos' directory
  2. Initial pressure and Final pressure
  3. Number of neutron stars you want to calcualte in the interval between Initial pressure and Final pressure

Output:
  1. M-R data in text file inside 'results' directory.
   First column: Central pressure
   Second column: Total mass of neutron star in M_0
   Third column: Total radius of neutron star in km

Terminal Outputs:
  1. Plot of M-R
  2. Mass, radius, outer crust, inner crust length, and number of iterations of the rk4 method.
  3. Final, maximum mass, and maximum radius of neutron stars in interval of Initial pressure and Final pressure.
