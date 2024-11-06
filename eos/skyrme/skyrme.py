import numpy as np
import units_cgs as cgs
import matplotlib.pyplot as plt
import scipy as sci

# zero-range Skyrme force, the density functional, symmetric matter
# E. Chabanat et al./Nuclear Physics A 627 (1997) 710-746
# Asymmetric infinite nuclear matter

# coefficients
rho_0 = 0.16 / cgs.fm ** 3
# SIII, SGII, SkM*, RATP, SkP, T6
sigma = np.array([1.0, 1.0/6.0, 1.0/6.0, 1.0/5.0, 1.0/6.0, 1.0/3.0])
W0 = np.array([120.0, 105.0, 130.0, 120.0, 100.0, 107.0])
t0 = np.array([-1128.75, -2645.0, -2645.0, -2160.0, -2931.7, -1794.2])
t1 = np.array([395.0, 340.0, 410.0, 513.0, 320.62, 294.0]) 
t2 = np.array([-95.0, -41.90, -135.0, 121.0, -337.41, -294.0])
t3 = np.array([1.4e4, 15595.0, 15595.0, 11600.0, 18708.97, 12817.0])
x0 = np.array([0.45, 0.09, 0.09, 0.418, 0.29215, 0.392])
x1 = np.array([0.0, -0.0588, 0.0, -0.36, 0.65318, -0.5])
x2 = np.array([0.0, 1.425, 0.0, -2.29, -0.53732, -0.5])
x3 = np.array([1.0, 0.06044, 0.0, 0.586, 0.18103, 0.5])

i = 1
sigma = sigma[i]
W0 = W0[i] * cgs.MeV * cgs.fm ** 5
t0 = t0[i] * cgs.MeV * cgs.fm ** 3
t1 = t1[i] * cgs.MeV * cgs.fm ** 5
t2 = t2[i] * cgs.MeV * cgs.fm ** 5
t3 = t3[i] * cgs.MeV * cgs.fm ** (3 + 3 * sigma)
x0 = x0[i]
x1 = x1[i]
x2 = x2[i]
x3 = x3[i]
theta_s = 3 * t1 + (5 + 4 * x2) * t2

def F(m, Y_p):
    val = 2 ** (m-1) * (Y_p ** m + (1 - Y_p) ** m)
    return val

E_A = lambda rho, Y_p: 3 * cgs.hbar ** 2 / (10 * cgs.m_n) * (3.0 * cgs.pi ** 2 * rho / 2.0) ** (2.0/3.0) * F(5.0/3.0, Y_p) \
                + 1.0/8.0 * t0 * rho * (2 * (x0 + 2) - (2 * x0 + 1) * F(2.0, Y_p)) \
                + 1.0/48.0 * t3 * rho ** (sigma + 1) * (2 * (x3 + 2) - (2 * x3 + 1) * F(2.0, Y_p)) \
                + 3.0/40.0 * (3.0 * cgs.pi ** 2 / 2.0) ** (2.0/3.0) * rho ** (5.0/3.0) * (
                    (t1*(x1+2) + t2*(x2+2)) * F(5.0/3.0, Y_p) 
                    + 1.0/2.0 * (t2*(2*x2+1) - t1*(2*x1+1)) * F(8.0/3.0, Y_p))
E_A0 = lambda rho: E_A(rho, 0.1) + 16.0 * cgs.MeV

e = lambda rho, Y_p: E_A(rho, Y_p) * rho

# finding rho_0
rho_0 = sci.optimize.root(E_A0, 0.04 / cgs.fm ** 3, method = 'anderson', tol = 0.01).x
print(rho_0 * cgs.fm ** 3)
print(E_A(rho_0, 0.1) / cgs.MeV)

