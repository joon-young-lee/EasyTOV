import math


c = 2.99792458e10 #cm/s
c2 = c ** 2 
hbar = 1.05457266e-27 # erg * s
G = 6.67259e-8 # cm ** 3 * g ** (-1) * s ** (-2)
M0 = 1.9884e33 # g
pi = math.pi
MeV = 1.602179e-6 # erg
fm = 1.0e-13 # cm
MeV_fm_to_cgs = 1.602179e+33
# MeV_fm_to_cgs = 1.0
erg_cm_to_MeV_fm = 6.2414999e-34
hbar3 = hbar ** 3
m_n = 1.6749286e-24 # g
m_p = 1.67262192e-24 # g
m_e = 9.1093837e-28 # g
 
# print(5.38125462909e23 * erg_cm_to_MeV_fm)

# units in MeV and fm^3
# c = G = 1
g_to_km = G / c2 * 1.e-5  # grams in km
MeV_to_km = 1.60218 * 1.e-6 * g_to_km / c2 # MeV in km
fm_to_km = 1.e-18
MeV_fm_to_km = MeV_to_km / fm_to_km ** 3
km_to_M0 = 1 / (M0 * g_to_km) # 1M0 
km_to_gram = 1 / g_to_km


n_0 = 0.16 # (fm ** 3)
# Empirical formula units in MeV
E_0 = 22.1
A = -118.2
B = 65.39
sigma = 2.112
S_0 = 30.0
GeV_m_n = 939.57
def main():
    # print(5.38125462909e23 * erg_cm_to_MeV_fm)
    # print(3.3587e-10 * MeV_fm_to_km)
    # print(4.445229102499874e-16 / MeV_fm_to_km)
    # print(1.3234969191948892e-6 / MeV_fm_to_km)
    # print(g)
    # print(MeV_to_km)
    # print(km_to_M0)
    # print(g_to_km)
    # print(MeV_to_km)
    print(1.3234969191948892e-3 / MeV_fm_to_km)
    print(MeV_fm_to_km)
    return 0

if __name__ == '__main__':
    main()

