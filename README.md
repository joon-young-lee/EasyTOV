# Easy Tolman-Oppenheimer-Volkoff equation solver using RK4 method.

# How to run multi-core:
  ```console
  $ python -u tov.py
  ```
# How to run single-core:
  ```console
  $ python -u single_processing_tov.py
  ```
# Inputs:
  1. Name of EoS data file in 'eos' directory. (First column pressure, second column energy density both in MeV/fm^3 units.)
  2. Initial pressure and Final pressure for central pressure.
  3. Number of neutron stars you want to calcualte in the interval between Initial pressure and Final pressure.

# Output (makes automatic directory of the current date):
  1. Plot of M-R
  2. M-R data in text file inside 'results' directory. <br />
   First column: Central pressure in MeV/fm^3 <br />
   Second column: Total mass of neutron star in M_0 <br />
   Third column: Total radius of neutron star in km
  
# Output Log:
  1. Mass, radius, outer crust, inner crust length, and number of iterations of the rk4 method.
  2. Final, maximum mass, and maximum radius of neutron stars in interval of Initial pressure and Final pressure.

# Updates:
  1. Change functions to classes.
  2. Conversion to Fortran 90 is expected!

# Reference
  1. Chagne of metric: tovh.pdf, James Lattimer
  2. Crust EoS: Baym, 1971, inlcusion of inner and outer crust is default. (unchangable for now!)
