import numpy as np
import matplotlib.pyplot as plt
import units_cgs as cgs

# Read the data from the file
data = np.loadtxt('./tov_data_.txt')

# Separate the columns into variables

p0 = data[:, 0]
M = data[:, 1]
R = data[:, 2]

# Create plots
plt.figure(figsize=(12, 6))

# Plot initial pressure vs mass
plt.subplot(1, 2, 1)
plt.plot(p0, M, label=r'Mass in $M_\odot$')
plt.xlabel(r'Initial Pressure($p_0$)')
plt.ylabel(r'Total Mass($M_\odot$)')
plt.xscale('log')
plt.title('Initial pressure vs Mass')
plt.legend()

# Plot mass vs radius
plt.subplot(1, 2, 2)
plt.plot(p0, R, label='Radius in km', color='red')
plt.xlabel(r'Initial Pressure($p_0$)')
plt.ylabel('Total Radius(km)')
plt.xscale('log')
plt.title('Initial pressure vs Radius')
plt.legend()

# Show plots
plt.tight_layout()
plt.savefig('p0_M_R.png')
# plt.show()


# Create the First Plot
fig, ax1 = plt.subplots(figsize=(10, 6))  # Adjust the figure size here

# Create Twin Axes
ax2 = ax1.twinx()

# Plot the Second Data on Twin Axes
ax1.plot(p0, R, '--',color='k', label='Radius in km')
ax1.set_xlabel(r'Initial Pressure($dyne/cm^2$)')
ax1.set_ylabel('Radius in km', fontsize=15)
ax1.tick_params(axis='y')
ax1.set_xscale('log')

ax1.legend(fontsize=15, loc='center right')

ax2.plot(p0, M, color='k', label=r'Mass in $M_\odot$')
ax2.set_ylabel(r'Mass in $M_\odot$', fontsize=15)
ax2.tick_params(axis='y')
ax2.set_xscale('log')

ax2.legend(fontsize=15, loc='center left')

plt.savefig('p0_M_R.png')




# plt.grid()
plt.show()