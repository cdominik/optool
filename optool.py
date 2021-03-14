"""
NAME
    optool

DESCRIPTION

    This module provides an interface to the optool program (available
    at https://github.com/cdominik/optool), and tools to plot and
    convert the results.
    It also provides tools to prepare refractive index data for use
    with the tool.
"""

import numpy as np
import matplotlib.pyplot as plt
import math as m
import re
import os
import subprocess
from distutils.spawn import find_executable
import random

class particle:
    """
NAME
    optool

DESCRIPTION

    Provides an interface to the optool program for computing dust opacities.

    The optool program can be found at https://github.com/cdominik/optool .

Arguments
---------

    cmd : string 
         A shell command to run optool.  The output produced by this command
         will be read in and stored in an instance of the optool.particle
         class.

Keywords
--------

    keep : Boolean, default False
           When True, do not clean up the directory with the output from
           running optool.
    """
    def __init__(self,cmd,keep=False):
        "Create a new optool.particle opject."
        self.cmd = cmd
        
        # Convert command string into list if necessary
        if (isinstance(cmd, str)):
            cmd = cmd.split()

        if cmd[0].startswith("~"):
            cmd[0] = os.path.expanduser(cmd[0])
            
        # Find the optool executable
        try:
            bin = find_executable(cmd[0])
        except:
            print('ERROR: executable not found:',cmd[0])
            return -1
            
        if (not bin):
            print('ERROR: executable not found:',cmd[0])
            return -1

        # Wrap the main part into try - finally to make sure we clean up
        try:
            # create a directory for the output and make sure it is empty
            random.seed(a=None)
            tmpdir = 'optool_tmp_output_dir_'+str(int(random.random()*1e6))
            os.system('rm -rf '+tmpdir)
            os.system('mkdir '+tmpdir)
            cmd.append('-o'); cmd.append(tmpdir)
    
            # Run optool to produce the opacities
            cmd[0] = bin; subprocess.Popen(cmd).wait()
            
            # Check if there is output we can use
            scat,ext = check_for_output(tmpdir)
            self.scat = scat
            self.massscale = 1.
    
            kabs=[]; ksca=[]; kext=[]; gg=[]
            f11=[]; f12=[]; f22=[]; f33=[]; f34=[]; f44=[]
            nfiles=0; header=[];
            materials = []
            rho = []
            
            for i in range(500):
                if scat:
                    file = ("%s/dustkapscatmat_%03d.%s") % (tmpdir,(i+1),ext)
                else:
                    file = ("%s/dustkappa_%03d.%s") % (tmpdir,(i+1),ext)
                if (not os.path.exists(file)): break
                nfiles = nfiles+1
                x = readoutputfile(file,scat)
                header.append(x[0])
                lam = x[1]
                kabs.append(x[2])
                ksca.append(x[3])
                kext.append(x[2]+x[3])
                gg.append(x[4])
                if scat:
                    scatang = x[5]
                    f11.append(x[6])
                    f12.append(x[7])
                    f22.append(x[8])
                    f33.append(x[9])
                    f34.append(x[10])
                    f44.append(x[11])
                    self.scat = scat
                self = parse_headers(header,self)
                self.nlam = len(lam)
                self.kabs = np.array(kabs)
                self.ksca = np.array(ksca)
                self.kext = np.array(kext)
                self.gsca = np.array(gg)
                self.lam  = lam
                if scat:
                    self.nang = len(scatang)
                    self.scatang = scatang
                    self.f11  = np.array(f11)
                    self.f12  = np.array(f12)
                    self.f22  = np.array(f22)
                    self.f33  = np.array(f33)
                    self.f34  = np.array(f34)
                    self.f44  = np.array(f44)
                else:
                    self.nang = 0
            self.nsize = nfiles
            self.np = nfiles
        finally:
            if keep:
                print("Keeping the temporary directory for inspection: "+tmpdir)
            else:
                print("Cleaning up temporary directory "+tmpdir)
                os.system('rm -rf '+tmpdir)

    def plot(self):
        """Create interactive plots of the opacities in SELF."""

        # Check if mean opacities have been computed
        if hasattr(self, 'kplanck'):
            # llamfmt = np.round(np.log10(self.lam),decimals=3)
            kplanck = self.kplanck
            kross   = self.kross
            temp = self.temp
            viewarr([kplanck,kross],index=1,ylabel=['kplanck','kross'],
                    idxnames=['grain index','log lambda [um]'],
                    idxvals=[np.array(range(self.nsize))+1,temp])

        # Extract the kappas and g
        kabs   = np.copy(self.kabs)
        ksca   = np.copy(self.ksca)
        kext   = kabs+ksca
        gg     = np.copy(self.gsca)
    
        # limit the kappa plotting range
        minkap = 1e0
        kabs   = np.maximum(kabs,minkap)
        ksca   = np.maximum(ksca,minkap)
        kext   = np.maximum(kext,minkap)
    
        # We will plot the logarithms of the Kappa values
        kabs   = np.log10(kabs)
        ksca   = np.log10(ksca)
        kext   = np.log10(kext)
    
        # Scale g such that it will fill the y range of the kappa plot
        kmin   = np.amin(np.array([np.amin(kabs),np.amin(ksca),np.amin(kext)]))
        kmax   = np.amax(np.array([np.amax(kabs),np.amax(ksca),np.amax(kext)]))
        ggscal = gg*(kmax-kmin)+kmin
        
        # Extract and plot the scattering matrix elements
        if self.scat:
            bottom = 1e-2
            f11 = logscale_with_sign(np.copy(self.f11),bottom)
            f12 = logscale_with_sign(np.copy(self.f12),bottom)
            f22 = logscale_with_sign(np.copy(self.f22),bottom)
            f33 = logscale_with_sign(np.copy(self.f33),bottom)
            f34 = logscale_with_sign(np.copy(self.f34),bottom)
            f44 = logscale_with_sign(np.copy(self.f44),bottom)
            f00 = f11*0.
    
            # Make version of grid variables with fewer digits
            lamfmt  = np.round(self.lam,decimals=3)
            angfmt  = np.round(self.scatang,decimals=3)

            # interactive plot of the scattering matric elements
            viewarr([f00,f00+2,f00-2,f00+4,f00-4,f11,f12,f22,f33,f34,f44],
                    index=2,ylabel=['<1e-2','±1','','±1e2','','f11','f12',
                                    'f22','f33','f34','f44'],
                    idxnames=['grain index','lambda [um]','angle'],
                    idxvals=[np.array(range(self.nsize))+1,lamfmt,angfmt])

        # interactive plot of kabs, ksca, kext, and g
        llamfmt = np.round(np.log10(self.lam),decimals=3)
        viewarr([ggscal,kext,ksca,kabs],index=1,
                ylabel=['gg','kext','ksca','kabs'],
                idxnames=['grain index','log lambda [um]'],
                idxvals=[np.array(range(self.nsize))+1,llamfmt])

    def scatnorm(self,norm=""):
        """Check or change the normalization of the scattering matrix.

        Without an argument, check the current normalization of the
        scattering matrix.

            p = optool.particle('./optool -s')
            p.scatnorm()

        Calling the method with an argument will change the normalization
        to one of the following conventions

            'b'  Bohren & Huffman
            'm'  Mishchenko
            'r'  RADMC-3D
            'h'  Hovenier
        """

        # analyze the NORM parameter
        if (norm == ""):
            renorm = False
            conv = self.norm
        else:
            renorm = True
            conv = norm
            
        conv = conv.lower()
        if (conv in ['h','hovenier']):
            self.norm = "hovenier"
            name = "Hovenier"
            normalization = "4 pi"
            units = "sr^-1"
        elif (conv in ['b','bh','bohren','bohrenhuffman']):
            self.norm = "bohrenhuffman"
            name = "Bohren & Huffman"
            normalization = "kappa_scat m_grain (2pi/lambda)^2"
            units = "sr^-1"
        elif (conv in ['m','mish','mishchenko']):
            self.norm = "mishchenko"
            name = "Mishchenko"
            normalization = "kappa_scat m_grain"
            units = "cm^2 sr^-1"
        elif (conv in ['r','radmc','radmc3d']):
            self.norm = "radmc3d"
            name = "RADMC-3D"
            normalization = "kappa_scat"
            units = "cm^2 g^-1 sr^-1"
        else:
            print("ERROR: Unknown normalization ",conv)
        
        ang   = self.scatang
        lam   = self.lam
        wav   = 2.*np.pi/(lam*1e-4)      # need cm here, not micrometer
        ratio = np.zeros([self.np,self.nlam])

        # Compute values and weights for the integration
        if (self.gridtype == "boundary"):
            # Matrix values are on cell boundaries
            if (ang[0] != 0):
                print("Error: inconsistency between gridtype \"boundary\" and angle values")
                return -1
            thetab = ang*np.pi/180.
            mub = np.cos(thetab)
            dmu = mub[:-1]-mub[1:]   # Defined negatively for mu integral
            fc = 0.5*(self.f11[:,:,1:]+self.f11[:,:,:-1]) 
        else:
            # This is the standard grid with values on cell midpoints
            if (ang[0] == 0):
                print("Error: inconsistency between gridtype \"center\" and angle values")
                return -1
            th1 = (ang-0.5)*np.pi/self.nang; mu1 = np.cos(th1)
            th2 = (ang+0.5)*np.pi/self.nang; mu2 = np.cos(th2)
            dmu = mu1-mu2  # Defined negatively for the mu integral
            fc  = self.f11

        for ip in (range(self.np)):
            for il in (range(self.nlam)):
                integ = 2.*np.pi*np.sum(fc[ip,il,:]*dmu)
                if (self.norm == "radmc3d"):
                    nn = self.ksca[ip,il]
                elif (self.norm == "hovenier"):
                    nn = 4.*np.pi
                elif (self.norm == "bohrenhuffman"):
                    mgrain = (4./3.)*np.pi * self.a3[ip]**3 * self.rho[ip]
                    nn = self.ksca[ip,il] * wav[il]**2 * mgrain
                elif (self.norm == "mishchenko"):
                    mgrain = (4./3.)*np.pi * self.a3[ip]**3 * self.rho[ip]
                    nn = self.ksca[ip,il] * mgrain
                if (norm):
                    self.f11[ip,il,:] = self.f11[ip,il,:] * nn/integ
                    self.f12[ip,il,:] = self.f12[ip,il,:] * nn/integ
                    self.f22[ip,il,:] = self.f22[ip,il,:] * nn/integ
                    self.f33[ip,il,:] = self.f33[ip,il,:] * nn/integ
                    self.f34[ip,il,:] = self.f34[ip,il,:] * nn/integ
                    self.f44[ip,il,:] = self.f44[ip,il,:] * nn/integ
                    ratio[ip,il] = 1.
                else:
                    ratio[ip,il] = integ / nn
        
        if (norm):
            print("New     nomalization is       ",name," convention")
        else:
            print("Current nomalization is       ",name," convention")
        print("Units of matrix elements are  ",units)
        print("Integral F_11 d Omega =       ",normalization)
        if (not norm):
            maxerr = np.amax(np.abs(ratio-1.))
            print("Maximum deviation              %7.2e" % maxerr)

    def computemean(self, tmin=10., tmax=1500., ntemp=100):
        """Compupte mean opacities from the opacities in self.

        Keyword parameters are
        tmin - minimum temperature for which to compute mean opacities
        tmax - maximum temperature for which to compute mean opacities
        ntemp - number of temperature steps between tmin and tmax
        """
        self.tmin    = tmin
        self.tmax    = tmax
        self.ntemp   = ntemp
        self.temp    = np.logspace(np.log10(tmin),np.log10(tmax),ntemp)
        self.kross   = np.zeros([self.nsize,self.ntemp])
        self.kplanck = np.zeros([self.nsize,self.ntemp])

        cl = 2.99792458e10          # Speed of light [cgs]
        nu = 1e4*cl/self.lam        # 10^4 because lam is in um - we need cm
        dnu = -1. * np.hstack([nu[1]-nu[0],0.5 * (nu[2:]-nu[:-2]), nu[-1] - nu[-2]  ])

        for it in range(self.ntemp):
            bnu    = bplanck(self.temp[it],nu)
            bnudt  = bplanckdt(self.temp[it],nu)
            dumbnu = np.sum(bnu*dnu)
            dumdb  = np.sum(bnudt*dnu)
            for ip in range(self.nsize):
                kap_p  = np.sum(self.kabs[ip,:]*bnu*dnu) / dumbnu
                kap_r  = dumdb / np.sum(bnudt * dnu / ( self.kabs[ip,:] + self.ksca[ip,:]*(1.-self.gsca[ip,:])))
                self.kplanck[ip,it] = kap_p
                self.kross[ip,it]   = kap_r


    def __add__(s,o):
        """Addition of optool.particle objects.
        This can be used to mix different grain types together
        into a dust model.
        
        # Make a silicate grain and a carbonatieous grain
        p1 = optool.particle('./optool -a 0.01 0.3 pyr-mg70')
        p2 = optool.particle('./optool -1 0.03 0.1 c-z')

        # Mix the particles with a mass ration 0.75 : 0.25
        # Make sure abundances add up to 1, or the opacities will
        # not be per g of dust!
        p = 0.75*p1 + 0.25*p2

        # Apply a dust-to-gas ratio, so that the opacities will be
        # per unit of GAS mass
        dtg = 0.01
        p   = dtg * p

        # Plot the opacities
        p.plot()

        """
        #
        # First, check if the particles are compatible
        #
        if ((s.np > 1) or (o.np>1)):
            raise NameError('Cannot add multi-particles')
        if ((s.nlam != o.nlam) or (np.abs((s.lam-o.lam)/s.lam).any()>1e-4)):
            raise NameError('Wavelength grids differ')
        if (s.scat):
            if ((s.nang != o.nang) or (np.abs((s.scatang[1:]-o.scatang[1:])/s.scatang[1:]).any()>1e-4)):
                # We don't check the first value, could be 0
                raise NameError('Angular grids differ')
            if (s.norm != o.norm):
                raise NameError('Scattering normalizations differ')
        #
        # Now do the adding
        #
        import copy
        x = copy.deepcopy(s)
        x.kabs = x.kabs+o.kabs
        x.ksca = x.ksca+o.ksca
        x.kext = x.kext+o.kext
        # F11 is linear in the integral for the computation of g.
        # So we can just take the weighted mean for g.
        x.gsca = (x.ksca*x.gsca + o.ksca*o.gsca) / (x.ksca+o.ksca)
        x.massscale = s.massscale + o.massscale
        if s.scat:
            # There is a scattering matrix
            for ip in range(s.np):
                for il in range(s.nlam):
                    if s.norm == 'hovenier':
                        # average, weighted by kappa_scat
                        ws, wo, wn = s.ksca[ip,il],o.ksca[ip,il],s.ksca[ip,il]+o.ksca[ip,il]
                    else:
                        # Add the values
                        ws, wo, wn = 1.,1.,1.
                    x.f11[ip,il,:] = (s.f11[ip,il,:]*ws + o.f11[ip,il,:]*wo) / wn
                    x.f12[ip,il,:] = (s.f12[ip,il,:]*ws + o.f12[ip,il,:]*wo) / wn
                    x.f22[ip,il,:] = (s.f22[ip,il,:]*ws + o.f22[ip,il,:]*wo) / wn
                    x.f33[ip,il,:] = (s.f33[ip,il,:]*ws + o.f33[ip,il,:]*wo) / wn
                    x.f34[ip,il,:] = (s.f34[ip,il,:]*ws + o.f34[ip,il,:]*wo) / wn
                    x.f44[ip,il,:] = (s.f44[ip,il,:]*ws + o.f44[ip,il,:]*wo) / wn
        #
        # Invalidate variables tha no longer make sense.
        #
        x.materials = np.hstack((x.materials,o.materials))
        if (x.fmax    != o.fmax   ): x.fmax    = -1
        if (x.pcore   != o.pcore  ): x.pcore   = -1
        if (x.pmantle != o.pmantle): x.pmantle = -1
        if (x.amin    != o.amin   ): x.amin    = -1
        if (x.amax    != o.amax   ): x.amax    = -1
        if (x.nsub    != o.nsub   ): x.nsub    = -1
        if (x.apow    != o.apow   ): x.apow    = -1
        if (x.rho     != o.rho    ): x.rho     = -1
        if (x.chop    != o.chop   ): x.chop    = -1
        x.a1,x.a2,x.a3 = -1,-1,-1

        return x
        
    def __mul__(s,o):
        """Multiplication for optool.particle objects.
        
        This is indended for the multiplication of such an object with
        a number.  The way to think about it is like this.  Such an
        contains opacities in units cm^2/g.  Multiplying it with a
        number means that the opacities are now per a different mass.
        This sounds strange, but it makes sens together with addition
        of particles - which see.
        """
        import copy
        x = copy.deepcopy(s)
        x.kabs = x.kabs*o
        x.ksca = x.ksca*o
        x.kext = x.kext*o
        x.massscale = x.massscale*o
        if s.norm != 'hovenier':
            # We need to change the matrix as well,
            # because is it normalized to ksca
            x.f11 = x.f11*o
            x.f12 = x.f12*o
            x.f22 = x.f22*o
            x.f33 = x.f33*o
            x.f34 = x.f34*o
            x.f44 = x.f44*o
        return x

    def __rmul__(s,o):
        """Rightsided multiplication of optool.particle object by a number."""
        return s*o
    def __div__(s,o):
        """Division of optool.particle object by a number."""
        return s * (1./o)
    def __truediv__(s,o):
        """Division of optool.particle object by a number."""
        return s * (1./o)

class lnktable:
    """NAME
    
    optool.lnktable

DESCRIPTION

    The is a clall to work with lnk files. lnk stands for lambda, n,
    and k, where and and k are the real and imaginary components of
    the refractive index of a material.
    

Arguments
---------
        
   file : string
          the file name from which to read the lnk data

Keywords
--------

   i_lnk : numpy array
           the column numbers where to find lambda, the real part of the
           refractive index and the imaginary part of it, respectively.
           The default is [1,2,3] .

   nskip : integer
        Number of lines to skil at the beginning.  Lines starting with
        `#', `!' or `*` are stored as header lines and ar skipped in
        this way. So this parameter is for dealing with files that are
        not yet formatted in the standard way for optool.  The default
        is 0.

   nlam_rho : Boolean
        True means, the first unskipped line contains the number of
        wavelengths points and the specific density of the material.
        False means no such line exists, and the lines have to be
        counted.  Rho will be se to 0 then, to indicate that the value
        is not know at this point.

Conversion
----------

    The standard format of these files is described in the optool user
    guide.  The class can also read files that are formatted
    differently, in order to create properly formatted version.  For
    example, if you have a file starting with 4 unimportant lines, and
    then data columns where n an k are in column 1 and 2,
    respectively, and the wavelength is given in units of cm^-1 in
    column 3, you can do the conversion in this way:

    new = optool.lnktable('x.dat',i_lnk=[3,1,2], nskip=4)
    new.lam = 10000./new.lam   # convert cm^-1 -> micrometer
    new.sort()                 # sort arrays according to lambda
    new.rho = 3.2              # set density in g/cm^3
    new.header = "# This is a silicate from Dorschner+1995)"
    new.write('sil-Dorschner1995.lnk')

    """
    def __init__(self,file,i_lnk=[1,2,3],nskip=0,nlam_rho=True):
        self.filename = file
        try:
            rfile = open(file, 'r')
        except:
            print('ERROR: file not found:',file)
            return -1
        print('Reading lnk file ',file,'...')

        # Skip lines that are irrelevant
        for i in range (nskip): dum = rfile.readline()

        # Read the header/comment field
        header = ''
        dum = rfile.readline()
        while ((dum.strip()[0]=='#') or (dum.strip()[0]=='*') or (dum.strip()[0]=='!')):
            header = header + dum
            dum = rfile.readline()
        self.header = header

        # Extract the number of wavelengths points, and the material density
        if (nlam_rho):
            dum = dum.split()
            self.nlam = int(dum[0])
            self.rho  = float(dum[1])
            dum = rfile.readline()
        else:
            self.nlam = 1
            self.rho  = 0.0
            print("Warning: density rho is nt known! Make sure to set it by hand.")

        # Prepare the arrays
        self.lam = []
        self.n   = []
        self.k   = []
        ilam     = 0
        
        # Fill the arrays
        while True:
            dum  = dum.split()
            ilam = ilam+1
            self.lam.append(float(dum[i_lnk[0]-1]))
            self.n.append(  float(dum[i_lnk[1]-1]))
            self.k.append(  float(dum[i_lnk[2]-1]))
            dum = rfile.readline()
            if ((len(dum) == 0) or dum.isspace()):
                # No more data. Truncate the arrays and stop reading
                if (not (self.nlam == ilam)):
                    print("WARNING: found %d lines of data, not %d" % (ilam,self.nlam))
                # Convert to numpy arrays and exit
                self.nlam = ilam
                self.lam = np.array(self.lam)
                self.n   = np.array(self.n)
                self.k   = np.array(self.k)
                break
        rfile.close()

    def sort(self):
        """Sort lam, n, and k according to lambda array."""
        sortinds = self.lam.argsort()
        self.lam = self.lam[sortinds]
        self.n   = self.n[sortinds]
        self.k   = self.k[sortinds]

    def smooth(self,size=10):
        """Smooth n and k with a medium filter of SIZE bins."""
        from scipy.ndimage import median_filter
        self.n = median_filter(self.n,size)
        self.k = median_filter(self.k,size)

    def decimate(self,step=2,size=0):
        """Decimate the arrays by a factor STEP.
        When SIZE is given instead, decimate to that size."""
        # FIXME: should we force to keep the first and last values?
        from math import floor
        if (size > 0):
            nlam = self.lam.size        
            step = floor(nlam/size)
            print("Decimating in steps of ",step," to reach size ",size)
        self.lam = self.lam[:-step:step]
        self.n   = self.n[:-step:step]
        self.k   = self.k[:-step:step]
        self.nlam = self.lam.size        

    def klimit(self,limit=0.):
        """Make sure imaginary part k is never smaller than LIMIT"""
        self.k[self.k<limit] = limit

    def fromwav(self):
        """Convert from wavenumbers and sort.
        Assuming that the self.lam array is actually wavenumbers,
        convert them to microns, and then sort the arraus so that
        lambda is increasing.
        """
        self.lam = 10000./self.lam
        self.sort()

    def compute_absorbance(self,dlayer=5):
        """Compute the absorbance for a thin layer of the material.

        dlayer   is the thickness of the layer to be used for the
                 computation, in micrometer units.  The default is 5um.

        Simple assumptions: infinite vacuum "substrate".
        R. Swaneloel 1983, J.Phys. E 16,1214
        We use the expressions (A2), for an infinite substrate."""
        n = self.n
        k = self.k
        s = 1.    # assume a vacuum subtrate
        lam = self.lam/10000.      # units are cm now
        d = dlayer /1e4            # Units: cm
        phi   = 4.*np.pi*n*d/lam
        alpha = 4.*np.pi*k/lam
        x  = np.exp(-alpha*d)
        A  = 16.*s*(n**2+k**2)
        B  =  ( (n+1.)**2 + k**2 ) * ( (n+s)**2 + k**22 )
        C1 = ( (n**2-1+k**2) * (n**2-s**2+k**2) + 4.*k**2*s  )    *2.*np.cos(phi)
        C2 = k * ( 2.*(n**2-s**2-k**2 ) + 2.* s *(n**2-1+k**2 ) ) *2.*np.sin(phi)
        C  = C1-C2
        D  = ( (n-1)**2 + k**2 ) * ( (n-s)**2 + k**2 )
        self.transmission  = A*x / (B-C*x-D*x**2)
        self.absorptivity = -np.log(self.trans)
        self.d = dlayer    # Record the density that was used.
            
    def plot(self):
        """Plot the refractive index aas a function of wavelength."""
        fig,ax = plt.subplots()
        ax.semilogx(self.lam,self.n,label='n',color="blue")
        ax.set_title(self.filename)
        ax.set_xlabel(r"log $\lambda$ [$\mu$m]")
        ax.set_ylabel(r'real part: $n$',color="blue")
        ax2=ax.twinx()
        ax2.loglog(self.lam,self.k,label='k',color="orange")
        ax2.set_ylabel(r'imaginary part: log $k$',color="orange")
        plt.show(block=False)

    def write(self,file):
        """Write the table to a file."""
        try:
            wfile = open(file, 'w')
        except:
            print('ERROR: Cannot write to file: ',file)
            return -1
        wfile.write(self.header)
        wfile.write("  %d  %g\n" % (self.nlam,self.rho))
        for i in range(self.nlam):
            wfile.write("  %16.6e %16.6e %16.6e\n" %
                        (self.lam[i],self.n[i],self.k[i]))
        wfile.close()

def logscale_with_sign(array,bottom):
    # Take the log10 of the absolute value of ARRAY, but transfer the
    # sign back onto the result.  Compress the region between
    # -BOTTOM and +BOTTOM into zero, smoothly.
    # This is a clever way to make a logarithmic plot of a variable
    # that has positive and negative values covering more then
    # one order of magnitude.
    lb = np.log10(bottom)
    a  =  np.where(array>0)
    b  =  np.where(array<=0)
    array[a] =  np.log10(array[a]+bottom)  - lb
    array[b] = -np.log10(-array[b]+bottom) + lb
    return array

def check_for_output(tmpdir):
    # Check for and if necessary rename input files
    for ext in['dat','inp']:
        if (os.path.exists(tmpdir+'/dustkapscatmat_001.'+ext)):
            return True, ext
        elif (os.path.exists(tmpdir+'/dustkappa_001.'+ext)):
            return False, ext
        elif (os.path.exists(tmpdir+'/dustkapscatmat.'+ext)):
            os.system('mv '+tmpdir+'/dustkapscatmat.'+ext+' '+tmpdir+'/dustkapscatmat_001.'+ext)
            return True, ext
        elif (os.path.exists(tmpdir+'/dustkappa.'+ext)):
            os.system('mv '+tmpdir+'/dustkappa.'+ext+' '+tmpdir+'/dustkappa_001.'+ext)
            return False, ext
    raise RuntimeError('No valid OpTool output files found')

def parse_headers(headers,b):
    # Extract information on run parameters from headers
    n = len(headers)
    b.amin  = np.zeros(n); b.amax = np.zeros(n); b.apow  = np.zeros(n)
    b.a1    = np.zeros(n); b.a2 = np.zeros(n); b.a3 = np.zeros(n);
    b.nsub  = np.zeros(n,dtype=np.int8)
    b.pcore = np.zeros(n); b.pmantle = np.zeros(n); b.fmax = np.zeros(n);
    b.chop  = np.zeros(n);
    b.materials = []
    b.rho = []

    for i in range(n):
        mat = []
        m = re.search(r" amin \[um\]\s*=\s*(-?[0-9.]+)",headers[i])
        b.amin[i]=float(m.group(1))
        m = re.search(r" amax \[um\]\s*=\s*(-?[0-9.]+)",headers[i])
        b.amax[i]=float(m.group(1))
        m = re.search(r" na\s*=\s*(-?[0-9.]+)",headers[i])
        b.nsub[i]=int(m.group(1))
        m = re.search(r" <a\^n>\s*=\s*([-+0-9.eE]+)\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)",headers[i])
        b.a1[i]=float(m.group(1))
        b.a2[i]=float(m.group(2))
        b.a3[i]=float(m.group(3))

        for m in re.finditer(r"^#\s+(core|mantle|grain)\s+([.0-9]+)\s+([.0-9]+)\s*(\S.*?)$",headers[0],re.MULTILINE):
            mat.append([m.group(1),float(m.group(2)),float(m.group(3)),m.group(4)])
            if m.group(1) == "grain": b.rho.append(float(m.group(3)))
        b.materials.append(mat)
    b.rho = np.array(b.rho)
    b.materials = np.array(b.materials)
    m = re.search(r" apow\s*=\s*(-?[0-9.]+)",headers[0])
    b.apow[i]=float(m.group(1))
    m = re.search(r" porosity\s*=\s*([0-9.]+)",headers[0])
    b.pcore[i]=float(m.group(1))
    m = re.search(r" p_mantle\s*=\s*(-?[0-9.]+)",headers[0])
    b.pmantle[i]=float(m.group(1))
    m = re.search(r" fmax\s*=\s*([0-9.]+)",headers[0])
    b.fmax[i]=float(m.group(1))
    m = re.search(r" chop\s*=\s*([0-9.]+)",headers[0])
    b.chop[i]=float(m.group(1))

    m = re.search(r" RADMC-3D",headers[0])
    if m:
        b.radmc = True
        b.gridtype = "boundary"
        b.norm = "radmc"
    else:
        b.radmc = False
        b.gridtype = "center"
        b.norm = "hovenier"

    return b

def readoutputfile(file,scat):
    # Read OpTool output file FILE.
    # scat=True marks if the file contains a scattering matrix
    try:
        rfile = open(file, 'r')
    except:
        print('ERROR: file not found:',file)
        return -1
    print('Reading ',file,'...')

    # Read the header/comment field
    header = ''
    dum = rfile.readline()
    while dum.strip()[0]=='#':
        header = header + dum
        dum = rfile.readline()

    # Read the file format
    iformat = int(dum)

    # Read the number of wavelengths in the file
    nlam = int(rfile.readline())

    if scat:
        # Read the scattering angular grid size
        dum = rfile.readline()
        while len(dum.strip())<2: dum = rfile.readline()
        nang = int(dum)

    # Prepare a few arrays
    lam=np.zeros(nlam); kabs=np.zeros(nlam); ksca=np.zeros(nlam); phase_g=np.zeros(nlam)
    if scat:
        scatang = np.zeros(nang)
        f11=np.zeros([nlam,nang]); f12=np.zeros([nlam,nang]); f22=np.zeros([nlam,nang])
        f33=np.zeros([nlam,nang]); f34=np.zeros([nlam,nang]); f44=np.zeros([nlam,nang])

    # Read the opacities
    for ilam in range(nlam):
        dum = rfile.readline()
        while len(dum.strip())<2: dum = rfile.readline()
        dum           = dum.split()
        lam[ilam]     = float(dum[0])
        kabs[ilam]    = float(dum[1])
        ksca[ilam]    = float(dum[2])
        phase_g[ilam] = float(dum[3])

    if scat:
        # Read the angular grid
        for iang in range(nang):
            dum        = rfile.readline()
            while len(dum.strip())<2: dum = rfile.readline()
            scatang[iang] = float(dum)

        # Read the scattering matrix
        for ilam in range(nlam):
            for iang in range(nang):
                dum = rfile.readline()
                while len(dum.strip())<2: dum = rfile.readline()
                dum = dum.split()
                f11[ilam,iang] = float(dum[0])
                f12[ilam,iang] = float(dum[1])
                f22[ilam,iang] = float(dum[2])
                f33[ilam,iang] = float(dum[3])
                f34[ilam,iang] = float(dum[4])
                f44[ilam,iang] = float(dum[5])

    rfile.close()
    if scat:
        return [header,lam,kabs,ksca,phase_g,scatang,f11,f12,f22,f33,f34,f44]
    else:
        return [header,lam,kabs,ksca,phase_g]

def viewarr(data,index=0,x=None,ymin=None,ymax=None,ylabel=None,idxnames=None,idxvals=None,idxformat=''):
    """
    For details about this function see https://github.com/dullemond/interactive_plot
    """
    if type(data)==list:
        shape =  data[0].shape
        ndim  = len(shape)
    else:
        shape = data.shape
        ndim  = len(shape)
    assert index<ndim, "Index out of range"
    idxorder  = list(range(ndim))
    idxorder.pop(index)
    idxorder.append(index)
    if type(data)==list:
        datatrans = []
        for d in data:
            datatrans.append(d.transpose(idxorder))
        shapetrans = datatrans[0].shape
    else:
        datatrans  = data.transpose(idxorder)
        shapetrans = datatrans.shape
    def func(x,param,fixedpar={"datatrans":datatrans}):
        datatrans = fixedpar["datatrans"]
        if type(datatrans)==list:
            answer = []
            for dslice in datatrans:
                for i in range(len(param)):
                    dslice = dslice[param[i]]
                answer.append(dslice)
        else:
            dslice = datatrans
            for i in range(len(param)):
                dslice = dslice[param[i]]
            answer = dslice
        answer = np.array(answer)
        return answer
    params=[]
    for i in range(ndim-1):
        params.append(np.arange(shapetrans[i]))
    if x is None:
        if idxvals is None:
            x = np.arange(shapetrans[-1])
        else:
            x = np.array(idxvals[index])
    if ymin is None:
        if type(data)==list:
            ymin = []
            for d in data:
                ymin.append(d.min())
            ymin = np.array(ymin).min()
        else:
            ymin = data.min()
    if ymax is None:
        if type(data)==list:
            ymax = []
            for d in data:
                ymax.append(d.max())
            ymax = np.array(ymax).max()
        else:
            ymax = data.max()
    if idxvals is not None:
        paramsalt = []
        for i in range(ndim-1):
            paramsalt.append(idxvals[idxorder[i]])
    else:
        paramsalt = None
    fig = None
    ax  = None
    if idxnames is None:
        parnames = []
        for i in range(ndim-1):
            s = 'Parameter {}'.format(idxorder[i])
            parnames.append(s)
        xname    = 'Parameter {}'.format(index)
    else:
        parnames = []
        for i in range(ndim-1):
            parnames.append(idxnames[idxorder[i]]+" =")
        xname    = idxnames[index]
    fig = plt.figure()
    ax  = plt.axes(xlim=(x.min(),x.max()),ylim=(ymin,ymax))
    ax.set_xlabel(xname)
    if ylabel is not None:
        if type(ylabel)==list:
            label = r''
            glue  = ''
            for l in ylabel:
                label += glue+l
                glue = ', '
            ax.set_ylabel(label)
        else:
            ax.set_ylabel(ylabel)
    if type(data)==list:
        axmodel = []
        if ylabel is None:
            for i in range(len(datatrans)):
                axm0,  = ax.plot(x,x,label='{}'.format(i))
                axmodel.append(axm0)
        else:
            for i in range(len(datatrans)):
                if (len(datatrans)>4 and i==0):
                    axm0,  = ax.plot(x,x,'-',color='0.9',linewidth=5,label=ylabel[i])
                elif (len(datatrans)>4 and i==1): 
                    axm0,  = ax.plot(x,x,'--',color='0.8',linewidth=2,label=ylabel[i])
                elif (len(datatrans)>4 and i==2): 
                    axm0,  = ax.plot(x,x,'--',color='0.8',linewidth=2,label=ylabel[i])
                elif (len(datatrans)>4 and i==3): 
                    axm0,  = ax.plot(x,x,':',color='0.6',linewidth=2,label=ylabel[i])
                elif (len(datatrans)>4 and i==4): 
                    axm0,  = ax.plot(x,x,':',color='0.6',linewidth=2,label=ylabel[i])
                else:
                    axm0,  = ax.plot(x,x,label=ylabel[i])
                axmodel.append(axm0)
        ax.legend()
    else:
        axmodel = None
    interactive_plot(x, func, params, ymin=ymin, ymax=ymax, parnames=parnames, parunits=None, fig=fig, ax=ax, axmodel=axmodel, parstart=None, iparstart=None, plotbutton=False, fixedpar=None, returnipar=False, block=False, paramsalt=paramsalt, altformat=idxformat)

def interactive_plot(x, func, params, ymin=None, ymax=None, parnames=None, parunits=None, fig=None, ax=None, axmodel=None, parstart=None, iparstart=None, plotbutton=False, fixedpar=None, returnipar=False, block=False, paramsalt=None, altformat='', **kwargs):
    """
    For details about this function see https://github.com/dullemond/interactive_plot
    """
    from matplotlib.widgets import Slider, Button, RadioButtons

    # Compute spacing of plot, sliders and button
    hslider  = 0.03
    nslidrscl= 6
    if(len(params)>nslidrscl):
        hslider *= float(nslidrscl)/len(params)
    dyslider = hslider*(4./3.)
    xslider  = 0.3
    wslider  = 0.3
    hbutton  = 0.06
    wbutton  = 0.15
    xbutton  = 0.3
    dybutton = hbutton+0.01
    panelbot = 0.0
    controlh = panelbot + len(params)*dyslider
    if plotbutton: controlh += dybutton
    controltop = panelbot + controlh
    bmargin  = 0.15
    
    # generate figure
    if fig is None: fig = plt.figure()
    fig.subplots_adjust(top=0.95,bottom=controltop+bmargin)

    # Set the initial values
    indexinit = np.zeros(len(params),dtype=int)
    if parstart is not None:
        for i in range(len(params)):
            if parstart[i] in params[i]:
                idx = np.where(np.array(params[i])==parstart[i])[0]
                if len(idx)>0:
                    indexinit[i] = idx[0]
            else:
                if params[i][-1]>params[i][0]:
                    idx = np.where(np.array(params[i])<parstart[i])[0]
                    if len(idx)>0:
                        indexinit[i] = idx[-1]
                else:
                    idx = np.where(np.array(params[i])>parstart[i])[0]
                    if len(idx)>0:
                        indexinit[i] = idx[0]
    if iparstart is not None:
        indexinit[:] = iparstart[:]

    # select first image
    par = []
    for i in range(len(params)):
        par.append(params[i][indexinit[i]])
    if fixedpar is not None:
        f = func(x,par,fixedpar=fixedpar)
    else:
        f = func(x,par)

    # set range
    if ymin is None: ymin = f.min()
    if ymax is None: ymax = f.max()
    
    # display function(s)
    if ax is None:      ax       = plt.axes(xlim=(x.min(),x.max()),ylim=(ymin,ymax))
    if axmodel is None:
        if len(f.shape)==1:
            # Normal case: a single model function
            axmodel, = ax.plot(x,f,**kwargs)
        else:
            # Special case: multiple model functions: f[imodel,:]
            assert len(f.shape)==2, 'Model returns array with more than 2 dimensions. No idea what to do.'
            axmodel = []
            for i in range(f.shape[0]):
                axm, = ax.plot(x,f[i,:],**kwargs)
                axmodel.append(axm)
            
    sliders = []
    for i in range(len(params)):
    
        # define slider
        axcolor = 'lightgoldenrodyellow'
        axs = fig.add_axes([xslider, controltop-i*dyslider, xslider+wslider, hslider], facecolor=axcolor)

        if parnames is not None:
            name = parnames[i]
        else:
            name = 'Parameter {0:d}'.format(i)

        slider = Slider(axs, name, 0, len(params[i]) - 1,
                    valinit=indexinit[i], valfmt='%i')
        sliders.append(slider)

    if plotbutton:
        axb = fig.add_axes([xbutton, panelbot+0.2*hbutton, xbutton+wbutton, hbutton])
        pbutton = Button(axb,'Plot')
    else:
        pbutton = None

    class callbackplot(object):
        def __init__(self,x,func,params,sliders,pbutton=None,fixedpar=None,ipar=None):
            self.x        = x
            self.func     = func
            self.params   = params
            self.sliders  = sliders
            self.pbutton  = pbutton
            self.fixedpar = fixedpar
            self.parunits = parunits
            self.paramsalt= paramsalt
            self.altformat= altformat
            self.closed   = False
            if ipar is None:
                self.ipar = np.zeros(len(sliders),dtype=int)
            else:
                self.ipar = ipar
        def handle_close(self,event):
            self.closed   = True
        def myreadsliders(self):
            for isl in range(len(self.sliders)):
                ind = int(self.sliders[isl].val)
                self.ipar[isl]=ind
            par = []
            for i in range(len(self.ipar)):
                ip = self.ipar[i]
                value = self.params[i][ip]
                par.append(value)
                name = self.sliders[i].label.get_text()
                if '=' in name:
                    namebase = name.split('=')[0]
                    if self.paramsalt is not None:
                        vls  = "{0:" + self.altformat + "}"
                        name = namebase + "= " + vls.format(self.paramsalt[i][ip])
                    else:
                        if self.parunits is not None:
                            valunit = self.parunits[i]
                        else:
                            valunit = 1.0
                        name = namebase + "= {0:13.6e}".format(value/valunit)
                    self.sliders[i].label.set_text(name)
            return par
        def myreplot(self,par):
            x = self.x
            if self.fixedpar is not None:
                f = self.func(x,par,fixedpar=self.fixedpar)
            else:
                f = self.func(x,par)
            if len(f.shape)==1:
                axmodel.set_data(x,f)
            else:
                for i in range(f.shape[0]):
                    axmodel[i].set_data(x,f[i,:])
            plt.draw()
        def mysupdate(self,event):
            par = self.myreadsliders()
            if self.pbutton is None: self.myreplot(par)
        def mybupdate(self,event):
            par = self.myreadsliders()
            if self.pbutton is not None: self.pbutton.label.set_text('Computing...')
            plt.pause(0.01)
            self.myreplot(par)
            if self.pbutton is not None: self.pbutton.label.set_text('Plot')

    mcb = callbackplot(x,func,params,sliders,pbutton=pbutton,fixedpar=fixedpar,ipar=indexinit)

    mcb.mybupdate(0)

    if plotbutton:
        pbutton.on_clicked(mcb.mybupdate)
    for s in sliders:
        s.on_changed(mcb.mysupdate)

    fig._mycallback    = mcb

    if block:
        plt.show(block=True)
    if returnipar:
        return mcb.ipar
        

def interactive_curve(t, func, params, xmin=None, xmax=None, ymin=None, ymax=None, parnames=None, parunits=None, fig=None, ax=None, axmodel=None, parstart=None, iparstart=None, plotbutton=False, fixedpar=None, returnipar=False, block=False, **kwargs):
    """
    For details about this function see https://github.com/dullemond/interactive_plot
    """
    from matplotlib.widgets import Slider, Button, RadioButtons

    # Compute spacing of plot, sliders and button
    hslider  = 0.03
    nslidrscl= 6
    if(len(params)>nslidrscl):
        hslider *= float(nslidrscl)/len(params)
    dyslider = hslider*(4./3.)
    xslider  = 0.3
    wslider  = 0.3
    hbutton  = 0.06
    wbutton  = 0.15
    xbutton  = 0.3
    dybutton = hbutton+0.01
    panelbot = 0.0
    controlh = panelbot + len(params)*dyslider
    if plotbutton: controlh += dybutton
    controltop = panelbot + controlh
    bmargin  = 0.15
    
    # generate figure
    if fig is None: fig = plt.figure()
    fig.subplots_adjust(top=0.95,bottom=controltop+bmargin)

    # Set the initial values
    indexinit = np.zeros(len(params),dtype=int)
    if parstart is not None:
        for i in range(len(params)):
            if parstart[i] in params[i]:
                idx = np.where(np.array(params[i])==parstart[i])[0]
                if len(idx)>0:
                    indexinit[i] = idx[0]
            else:
                if params[i][-1]>params[i][0]:
                    idx = np.where(np.array(params[i])<parstart[i])[0]
                    if len(idx)>0:
                        indexinit[i] = idx[-1]
                else:
                    idx = np.where(np.array(params[i])>parstart[i])[0]
                    if len(idx)>0:
                        indexinit[i] = idx[0]
    if iparstart is not None:
        indexinit[:] = iparstart[:]

    # select first image
    par = []
    for i in range(len(params)):
        par.append(params[i][indexinit[i]])
    if fixedpar is not None:
        x, y = func(t,par,fixedpar=fixedpar)
    else:
        x, y = func(t,par)

    # set range
    if xmin is None: xmin = x.min()
    if xmax is None: xmax = x.max()
    if ymin is None: ymin = y.min()
    if ymax is None: ymax = y.max()
    
    # display function
    if ax is None: ax   = plt.axes(xlim=(xmin,xmax),ylim=(ymin,ymax))
    if axmodel is None:
        if len(x.shape)==1:
            # Normal case: a single model function
            assert len(x.shape)==1, 'Cannot have multiple y and single x'
            axmodel, = ax.plot(x,y,**kwargs)
        else:
            # Special case: multiple model functions: f[imodel,:]
            assert len(x.shape)==2, 'Model returns array with more than 2 dimensions. No idea what to do.'
            assert len(y.shape)==2, 'Cannot have multiple x and single y'
            axmodel = []
            for i in range(x.shape[0]):
                axm, = ax.plot(x[i,:],y[i,:],**kwargs)
                axmodel.append(axm)
    
    sliders = []
    for i in range(len(params)):
    
        # define slider
        axcolor = 'lightgoldenrodyellow'
        axs = fig.add_axes([xslider, controltop-i*dyslider, xslider+wslider, hslider], facecolor=axcolor)

        if parnames is not None:
            name = parnames[i]
        else:
            name = 'Parameter {0:d}'.format(i)
            
        slider = Slider(axs, name, 0, len(params[i]) - 1,
                    valinit=indexinit[i], valfmt='%i')
        sliders.append(slider)

    if plotbutton:
        axb = fig.add_axes([xbutton, panelbot+0.2*hbutton, xbutton+wbutton, hbutton])
        pbutton = Button(axb,'Plot')
    else:
        pbutton = None

    class callbackcurve(object):
        def __init__(self,t,func,params,sliders,pbutton=None,fixedpar=None,ipar=None):
            self.t        = t
            self.func     = func
            self.params   = params
            self.sliders  = sliders
            self.pbutton  = pbutton
            self.fixedpar = fixedpar
            self.parunits = parunits
            self.closed   = False
            if ipar is None:
                self.ipar = np.zeros(len(sliders),dtype=int)
            else:
                self.ipar = ipar
        def handle_close(self,event):
            self.closed   = True
        def myreadsliders(self):
            for isl in range(len(self.sliders)):
                ind = int(self.sliders[isl].val)
                self.ipar[isl]=ind
            par = []
            for i in range(len(self.ipar)):
                ip = self.ipar[i]
                value = self.params[i][ip]
                par.append(value)
                name = self.sliders[i].label.get_text()
                if '=' in name:
                    namebase = name.split('=')[0]
                    if self.parunits is not None:
                        valunit = self.parunits[i]
                    else:
                        valunit = 1.0
                    name = namebase + "= {0:13.6e}".format(value/valunit)
                    self.sliders[i].label.set_text(name)
            return par
        def myreplot(self,par):
            t = self.t
            if self.fixedpar is not None:
                x,y = self.func(t,par,fixedpar=self.fixedpar)
            else:
                x,y = self.func(t,par)
            if len(x.shape)==1:
                axmodel.set_data(x,y)
            else:
                for i in range(x.shape[0]):
                    axmodel[i].set_data(x[i,:],y[i,:])
            plt.draw()
        def mysupdate(self,event):
            par = self.myreadsliders()
            if self.pbutton is None: self.myreplot(par)
        def mybupdate(self,event):
            par = self.myreadsliders()
            if self.pbutton is not None: self.pbutton.label.set_text('Computing...')
            plt.pause(0.01)
            self.myreplot(par)
            if self.pbutton is not None: self.pbutton.label.set_text('Plot')

    mcb = callbackcurve(t,func,params,sliders,pbutton=pbutton,fixedpar=fixedpar,ipar=indexinit)
            
    mcb.mybupdate(0)
        
    if plotbutton:
        pbutton.on_clicked(mcb.mybupdate)
    for s in sliders:
        s.on_changed(mcb.mysupdate)

    fig._mycallback    = mcb
    
    if block:
        plt.show(block=True)
    if returnipar:
        return mcb.ipar

def bplanck(temp,nu):
    """
----------------------------------------------------------------------------
                THE BLACKBODY PLANCK FUNCTION B_nu(T)

     This function computes the Blackbody function 

                    2 h nu^3 / c^2
        B_nu(T)  = ------------------    [ erg / cm^2 s ster Hz ]
                   exp(h nu / kT) - 1

     ARGUMENTS:
        nu    [Hz]            = Frequency (may be an array)
        temp  [K]             = Temperature
----------------------------------------------------------------------------
    """
    if (temp == 0.e0): return nu*0.e0
    bplanck = 1.47455e-47 * nu**3 /  (np.exp(4.7989e-11 * nu / temp)-1.e0) + 1.e-290
    return bplanck

def bplanckdt(temp,nu):
    """
----------------------------------------------------------------------------
           THE TEMPERATURE DERIVATIVE OF PLANCK FUNCTION 
     
      This function computes the temperature derivative of the
      Blackbody function 
      
         dB_nu(T)     2 h^2 nu^4      exp(h nu / kT)        1 
         --------   = ---------- ------------------------  ---
            dT          k c^2    [ exp(h nu / kT) - 1 ]^2  T^2
     
      ARGUMENTS:
         nu    [Hz]            = Frequency (may be an array)
         temp  [K]             = Temperature
----------------------------------------------------------------------------
    """
    bplanckdt = np.zeros(len(nu))
    exponent = 4.7989e-11*nu/temp
    mask = (exponent <= 76.)
    bplanckdt[mask] = 7.07661334104e-58 * nu[mask]**4 * np.exp(exponent[mask]) /  \
        ( (np.exp(exponent[mask])-1.e0)**2 * temp**2 ) + 1.e-290
    mask = (exponent > 76.)
    bplanckdt[mask] = 7.07661334104e-58 * nu[mask]**4 /  \
            ( np.exp(exponent[mask]) * temp**2 ) + 1.e-290
    return bplanckdt


def plotall():
    import matplotlib.pyplot as plt
    import numpy as np
    import optool

    files = [
        'lnk_data/pyr-mg100-Dorschner1995.lnk',
        'lnk_data/pyr-mg95-Dorschner1995.lnk',
        'lnk_data/pyr-mg80-Dorschner1995.lnk',
        'lnk_data/pyr-mg70-Dorschner1995.lnk',
        'lnk_data/pyr-mg60-Dorschner1995.lnk',
        'lnk_data/pyr-mg50-Dorschner1995.lnk',
        'lnk_data/pyr-mg40-Dorschner1995.lnk',
        'lnk_data/pyr-c-mg96-Jaeger1998.lnk',
        'lnk_data/ol-mg50-Dorschner1995.lnk',
        'lnk_data/ol-mg40-Dorschner1995.lnk',
        'lnk_data/ol-c-mg100-Steyer1974.lnk',
        'lnk_data/astrosil-Draine2003.lnk',

        'lnk_data/c-z-Zubko1996.lnk',
        'lnk_data/c-p-Preibisch1993.lnk',
        'lnk_data/c-gra-Draine2003.lnk',
        'lnk_data/c-org-Henning1996.lnk',
        'lnk_data/c-nano-Mutschke2004.lnk',

        'lnk_data/fe-c-Henning1996.lnk',
        'lnk_data/fes-Henning1996.lnk',
        'lnk_data/sic-Draine1993.lnk',

        'lnk_data/cor-c-Koike1995.lnk',

        'lnk_data/h2o-w-Warren2008.lnk',
        'lnk_data/co2-w-Warren1986.lnk',
        'lnk_data/nh3-m-Martonchik1983.lnk',
        
        'lnk_data/co-a-Palumbo2006.lnk',
        'lnk_data/co2-a-Gerakines2020.lnk',
        'lnk_data/co2-c-Gerakines2020.lnk',
        'lnk_data/ch4-a-Gerakines2020.lnk',
        'lnk_data/ch4-c-Gerakines2020.lnk',
        'lnk_data/ch3oh-a-Gerakines2020.lnk',
        'lnk_data/ch3oh-c-Gerakines2020.lnk'
       ]
    
    # Some example data to display
    x = np.linspace(0, 2 * np.pi, 400)
    y = np.sin(x ** 2)
    nx = 6
    ny = 6
    fig = plt.figure(figsize=(10,9))
    gs = fig.add_gridspec(nx,ny, hspace=0, wspace=0)
    gs1 = gs.subplots(sharex='all',sharey='all')
    print(gs1.shape)
    for iy in range(ny):
        for ix in range(nx):
            nn = ix+iy*nx
            print(ix,iy,nn)
            if (nn >= len(files)):
                break
            file = files[nn]
            p=optool.lnktable(file)
            ax = gs1[ix,iy]
            ax.loglog(p.lam,p.k+1e-5)
            ax.set_xlim(1,300)
            ax.set_ylim(1e-4,1e3)
            ax.text(1.2,100.,file[9:-4],fontsize='xx-small')
    fig.show()
    fig.savefig("maint/all_k.pdf", bbox_inches='tight')
        
    
