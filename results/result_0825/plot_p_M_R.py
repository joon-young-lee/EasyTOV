import numpy as np
import matplotlib.pyplot as plt


# Read the data from the file
data = np.loadtxt('./result_0825/nrel.dat')

# Separate the columns into variables
MeV_fm_to_cgs = 1.602179e+33
p0 = data[:, 0] * MeV_fm_to_cgs
M = data[:, 1]
R = data[:, 2]

# Create plots
# plt.figure(figsize=(12, 6))

# Plot initial pressure vs mass
# plt.subplot(1, 2, 1)
# plt.plot(p0, M, label=r'Mass in $M_\odot$')
# plt.xlabel(r'Initial Pressure($p_0$)')
# plt.ylabel(r'Total Mass($M_\odot$)')
# plt.xscale('log')
# plt.title('Initial pressure vs Mass')
# plt.legend()

# Plot mass vs radius
# plt.subplot(1, 2, 2)
# plt.plot(p0, R, label='Radius in km', color='red')
# plt.xlabel(r'Initial Pressure($p_0$)')
# plt.ylabel('Total Radius(km)')
# plt.xscale('log')
# plt.title('Initial pressure vs Radius')
# plt.legend()

# # Show plots
# plt.tight_layout()
# plt.savefig('p0_M_R.png')
# plt.show()
plt.figure(figsize=(12, 6))
plt.plot(R, M, color='k')
plt.title('Polytropic EoS M-R', fontsize=20)
plt.xlabel('Radius in km', fontsize=20)
plt.ylabel(r'Mass in $M_\odot$', fontsize=20)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.savefig('./result_0825/nrel_M_R.png')
# 9.244553E-03 1.015523E-01 3.248423E+01
plt.text(25.0, 0.25, r'initial pressure $9.24\times 10^{-3}MeV/fm^3$', 
         fontsize=18)
plt.arrow(30.62, 0.23, 3.248423e1-30.62, 1.015523e-1-0.23, color='b')
# 2.001268E+00 4.130981E-01 1.856932E+01
plt.text(20.5, 0.45, r'initial pressure $2.00 MeV/fm^3$'+'\n'+'becoming relativistic', fontsize=18)
plt.arrow(23.3, 0.43, 1.856932e1-23.3, 4.130981e-1-0.43, color='b')

plt.grid()
plt.savefig('./result_0825/nrel_M_R.png')
# plt.show()


# Create the First Plot
fig, ax1 = plt.subplots(figsize=(10, 6))  # Adjust the figure size here

# Create Twin Axes
ax2 = ax1.twinx()

# Plot the Second Data on Twin Axes
ax1.plot(p0, R, '--',color='k', label='Radius in km')
ax1.set_xlabel(r'Initial Pressure($dyne/cm^2$)', fontsize=20)
ax1.set_ylabel('Radius in km', fontsize=20)
ax1.tick_params(axis='y', labelsize='15')
ax1.set_xscale('log')
ax1.set_title(label='Polytropic EoS', fontsize=20)
ax1.legend(fontsize=15, loc='center right')

ax2.plot(p0, M, color='k', label=r'Mass in $M_\odot$')
ax2.set_ylabel(r'Mass in $M_\odot$', fontsize=20)
ax2.tick_params(axis='y', labelsize='15')
ax2.set_xscale('log')

ax2.legend(fontsize=15, loc='center left')

plt.savefig('./result_0825/nrel_p0_M_R.png')




# plt.grid()
plt.show()