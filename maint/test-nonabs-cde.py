
import numpy as np
import matplotlib.pyplot as plt

n = 1.3
k = np.logspace(-10,1,10)

m=np.zeros(10, dtype=np.complex)
cscat1 = np.zeros(10, dtype=np.float)
cscat2 = np.zeros(10, dtype=np.float)
for i in range(10):
    m[i] = complex(n,k[i])
    wvno = 1.
    V    = 1.
    pi   = 3.1415926
    print("ii",m[i],m[i]**2,(m[i]**2).imag)
    cscat1[i] =  wvno**4*V**2*(abs(m[i]**2-1e0))**2/(3e0*pi*(m[i]**2).imag) * (m[i]**2/(m[i]**2-1e0)*np.log(m[i]**2)).imag
    m1 = complex(m[i].real,0e0)
    cscat2[i] = wvno**4*V**2/(3e0*pi) * (m1**2-1e0-np.log(m1**2))

