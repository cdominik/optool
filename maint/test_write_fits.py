#!/usr/bin/env python3
"""Test optool.particle.write_fits.

Run from the optool directory, after compiling optool:

    python3 maint/test_write_fits.py

This needs numpy, matplotlib and astropy.  Without astropy, the test is
skipped.  The exit code is 1 if any check fails.
"""

import os
import sys
import tempfile

sys.path.insert(0, os.getcwd())
import numpy as np

try:
    from astropy.io import fits
except ImportError:
    print("test_write_fits: SKIPPED, astropy is not installed")
    sys.exit(0)

import optool

failures = 0

def check(ok, what):
    global failures
    print("%-60s %s" % (what, "ok" if ok else "FAILED"))
    if not ok:
        failures += 1

def expect_error(error, func, what):
    try:
        func()
    except error:
        check(True, what)
    else:
        check(False, what)

with tempfile.TemporaryDirectory() as d:

    # A particle with two materials and the scattering matrix
    p = optool.particle('./optool pyr-mg70 0.87 c-z 0.13 -p 0.25 -na 10 -nl 30 -s -q',
                        silent=True)
    file = os.path.join(d, 'a.fits')
    p.write_fits(file)
    h = fits.open(file)

    # Layout as documented in the UserGuide and read by code4reading/readfits.py
    check(len(h) == 2, "two HDUs")
    check(h[0].data.shape == (4, p.nlam), "HDU 0 has shape (4,nlam)")
    check(h[1].data.shape == (p.nang, 6, p.nlam), "HDU 1 has shape (nang,6,nlam)")

    o, sc = h[0].data, h[1].data
    check(np.array_equal(o[0], p.lam), "HDU 0: wavelength")
    check(np.array_equal(o[1], p.kext[0]), "HDU 0: kext")
    check(np.array_equal(o[2], p.kabs[0]), "HDU 0: kabs")
    check(np.array_equal(o[3], p.ksca[0]), "HDU 0: ksca")
    for k, f in enumerate([p.f11, p.f12, p.f22, p.f33, p.f34, p.f44]):
        check(np.array_equal(sc[:, k, :], f[0].T), "HDU 1: matrix element %d" % (k+1))

    hdr = h[0].header
    check(hdr['N_MAT'] == 2, "header: N_MAT")
    check(hdr['N_RADII'] == 10, "header: N_RADII")
    check(hdr['FILE01'] == 'pyr-mg70' and hdr['FILE02'] == 'c-z', "header: FILE01, FILE02")
    check(abs(hdr['FRAC01']-0.87) < 1e-6 and abs(hdr['FRAC02']-0.13) < 1e-6,
          "header: FRAC01, FRAC02")
    check(abs(hdr['POROSITY']-0.25) < 1e-6, "header: POROSITY")
    check(abs(hdr['R_MIN']-0.05) < 1e-6 and abs(hdr['R_MAX']-3000) < 1e-6,
          "header: R_MIN, R_MAX")
    h.close()

    expect_error(OSError, lambda: p.write_fits(file), "refuse to overwrite an existing file")
    p.write_fits(file, overwrite=True)
    check(True, "overwrite=True replaces the file")

    # Multi-particle objects need select()
    q = optool.particle('./optool pyr-mg70 -na 10 -nl 20 -d 3 -s -q', silent=True)
    expect_error(TypeError, lambda: q.write_fits(os.path.join(d, 'b.fits')),
                 "multi-particle object is refused")
    q.select(1).write_fits(os.path.join(d, 'b.fits'))
    check(abs(fits.getheader(os.path.join(d, 'b.fits'))['A1'] - q.a1[1]) < 1e-12,
          "select(1) writes the second particle")

    # Without -s there is no scattering matrix
    r = optool.particle('./optool pyr-mg70 -na 10 -nl 20 -q', silent=True)
    expect_error(TypeError, lambda: r.write_fits(os.path.join(d, 'c.fits')),
                 "particle without scattering matrix is refused")

print("test_write_fits: %s" % ("all checks passed" if failures == 0 else "%d FAILED" % failures))
sys.exit(1 if failures else 0)
