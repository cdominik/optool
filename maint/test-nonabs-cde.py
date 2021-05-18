# Compage the CDE results for the scattering cross section
# close to the non-absorbing case.  In optool, we use the
# analytical limit of k -> 0 when |k/n| <1e-6.  This test shows
# that this is OK.  1e-4 would be OK as well.


import numpy as np
import matplotlib.pyplot as plt

nvals = 100
n = 2.
k = n*np.logspace(-6,0,nvals)

m      = np.zeros(nvals, dtype=np.complex)
cscat1 = np.zeros(nvals, dtype=np.float)
cscat2 = np.zeros(nvals, dtype=np.float)
for i in range(nvals):
    m[i] = complex(n,k[i])
    wvno = 1.
    V    = 1.
    pi   = 3.1415926
    print(m[i],abs(m[i]**2),abs(m[i]**2-1.))
    cscat1[i] =  wvno**4*V**2*(abs(m[i]**2-1e0))**2/(3e0*pi*(m[i]**2).imag) * (m[i]**2/(m[i]**2-1e0)*np.log(m[i]**2)).imag
    m1 = complex(m[i].real,0e0)
    cscat2[i] = wvno**4*V**2/(3e0*pi) * (m1**2-1e0-np.log(m1**2))
line1, = plt.semilogx(k/n,cscat2)
line2, = plt.semilogx(k/n,cscat1)
plt.xlabel('k/n')
plt.ylabel('c_scat')
plt.show(block=False)
