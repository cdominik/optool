import optool
import numpy as np

# Reading a file without scattering data
file1 = 'dustkappa.dat'
header,lam,kabs,ksca,g = optool.readoutputfile(file1,scat=False)
nlam = lam.shape[0]

# Reading a file with scattering data
file = 'dustkapscatmat.dat'
header,lam,kabs,ksca,phase_g,scatang,f11,f12,f22,f33,f34,f44 = \
    optool.readoutputfile(file,scat=True)
nlam = lam.shape[0]
nang = f11.shape[1]
ang  = np.arange(nang) + 0.5*nang/180.


# The shape of the fxx arrays is (nlam,nang)
