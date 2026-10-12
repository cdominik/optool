# This is a Python code snippet showing how to read a fits file
# written by the write_fits method of optool.particle in optool.py

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

file = 'dustkappa.fits'
hdul = fits.open(file)

opac = hdul[0]
nlam = opac.shape[1]
lam  = opac.data[0,:]
kext = opac.data[1,:]
kabs = opac.data[2,:]
ksca = opac.data[3,:]

scat = hdul[1]
nang = scat.shape[0]
ang  = (np.arange(nang) + 0.5) * 180./nang
f11  = scat.data[:,0,:]
f12  = scat.data[:,1,:]
f22  = scat.data[:,2,:]
f33  = scat.data[:,3,:]
f34  = scat.data[:,4,:]
f44  = scat.data[:,5,:]
# The shape of the fxx arrays is (nang,nlam)
