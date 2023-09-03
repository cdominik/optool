"""
NAME
    optool

DESCRIPTION

    This module provides an interface to the optool program (available
    at https://github.com/cdominik/optool), and tools to plot, convert,
    and to compute with the results.
    It also provides tools to prepare refractive index data for use
    with the tool.

EXAMPLES

Compute pyroxene with an ice mantle in 24 different grain sizes
and plot the results

  import optool
  p = optool.particle(’~/bin/optool pyr 0.8 -m ice 0.2 -na 24 -d’)
  p.plot()

Compute the opacities of 100 olivine silicate grain sizes and of 50
carbon grain sizes, and store the opacities in cache directories. This
works by specifying the directory as the second argument. In a new
session, if the directories still exist and were produced using the
same commands, the opacities are simply read back in.

  import optool
  import numpy as np
  sil  = optool.particle('optool -d -a 0.001 100 0 100 ol-mg50',cache='sil')
  carb = optool.particle('optool -d -a 0.001 3.0 0 50  c',cache='carb')

Apply powerlaw size distributions, and limit the size of the
contributing grains.  Note that a power law f(a)\propto a^{-3.5}
implies using a power a^{-2.5} when computing the number of particles
per size bin on a logarithmic size grid. No normalization is
necessary - the =sizedist= method will take care of that.

  nsil = sil.a1**(-2.5)             # power law, no normalization required
  nsil[sil.a1<0.01] = 0             # no grains smaller than 0.01um
  nsil[sil.a1>0.3]  = 0             # no grains larger  than 0.3um
  sil_pl = sil.sizedist(nsil)       # pass the relative number for each size

  nc = carb.a1**(-2.5)              # power law, no normalization required
  nc[carb.a1>0.3]=0                 # no grains larger than 0.3um
  carb_pl = carb.sizedist(nc)       # pass the relative number for each size

sil_pl and carb_pl are now objects with a single opacity each,
obtained by adding opacities with the weights of the size
distribution. The opacities are still per g of total grain mass.
Let's add these two opacities with mass weights, to get something
resembling an interstellar dust opacity produced by a mixture of
silicate and carbon grains:

  ptot = 0.7*sil_pl + 0.3*carb_pl   # weights should add up to 1
  ptot.plot()                       # plot the resulting opacity

Now let's assume we are looking at an interstellar cloud, where the
dust is just one percent of the total mass.  We want to have the
opacity per unit of /gas mass/ instead, and we need Planck and
Rosseland mean opacities:

  p_ism = ptot * 0.01               # dilute the opacity
  p_ism.computemean(tmax=1300)      # Compute mean opacities
  p_ism.plot()                      # Plot the results

Other size distributions can be made just as easily.  Here is a
log-normal size distribution for the silicate grains, with a
peak abundance at a size of a_m=1.3 microns, and a logarithmic width
of \sigma=1.2:

  sil_ln = sil.sizedist( np.exp( -0.5*(np.log(sil.a1/1.3)/1.2)**2) )
  sil_ln.write('dkap_ln.dat')       # write opacity to a file

"""
import copy
import numpy as np
import matplotlib.pyplot as plt
import math as m
import re
import os
import shutil
import subprocess
from distutils.spawn import find_executable
import tempfile

class particle:
    """Run optool and turn output into a python object.

        Provides an interface to the optool program for computing dust
        opacities. The optool program can be found on GitHub, at this address:
        https://github.com/cdominik/optool .

        Attributes
        ----------

        cmd : str
             The full command given in the particle() call
        radmc : boolean
             Output follows RADMC conventions
        scat : boolean
             Scattering matrix is available
        nlam : int
             Number of wavelength points
        lam : float[nlam]
             The wavelength grid
        nang : int
             Number of scattering angles
        scatang : float[nang]
             The angular grid
        materials : [[[...]...]... ]
             Lists with [location,m_{frac},\rho,material
        np : int
             Number of particles, either 1 or (with -d) n_a
        fmax : float[np]
             Maximum volume fraction of vacuum for DHS
        pcore, pmantle : float[np]
             Porosity of the core/mantle material
        amin : float[np]
             min grain size used for each particle
        amax : float[np]
             max grain size used for each particle
        nsub : int[np]
             Number of sizes averaged for each particle
        apow : float[np]
             Negative size distribution power law (e.g. 3.5)
        amean : float[np]
             Mean size for (log-)nornal size distributions
        asig : float[np]
             Standard deviation for (log-)normal distribution
        a1 : float[np]
            Mean grain radius
        a2 : float[np]
            Radius of the grain with mean surface area
        a3 : float[np]
            Radius of the grain with mean volume
        rho : float[np]
             Specific density of grains
        kabs : float[np,nlam]
             Absorption cross section
        ksca : float[np,nlam]
             Scattering cross section
        kext : float[np,nlam]
             Extinction cross section
        gsca : float[np,nlam]
             Asymmetry parameter
        f11, ..., f44 : float[np,nlam,nang]
             Scattering matrix element F_11, ... ,F_44
        chop : float[np]
             Degrees chopped off forward scattering
        tmin : float
             Minimum temperature for mean opacities
        tmax : float
             Maximum temperature for mean opacities
        ntemp : int
             Number of temperatures for mean opacities
        temp : float[ntemp]
             Temperatures used for mean opacities
        kplanck : float[np,ntemp]
             Planck mean opacities, after calling computemean()
        kross : float[np,ntemp]
             Rosseland mean opacities, after calling computemean()
        norm : string
             Current scattering matrix normalization

        Methods
        -------

        plot()
             Plot the opacities and the scattering matrix

        computemean(tmin=10,tmax=1500,ntemp=100)
             Compute Planck and Rosseland mean opacities

        scatnorm(norm='')
             Check or change the normalization of the scattering matrix

        sizedist(N_of_a)
             Compute opacity of a size distribution of elements of SELF
        """
    def __init__(self,cmd,cache='',silent=False):
        """Create a new optool.particle opject.

        Parameters
        ---------=

        cmd  : str or False
               A shell command to run optool. The output produced by this
               command will be read in and stored in an instance of the
               optool.particle class.
               If this is False or the empty string, just read what an
               earlier run of optool has put into the directory given by
               the second parameter CACHE.

        cache  : str, optional
               The diretory to cache the optool output files in, so that
               they can be read instead of recomputed the next time
               the same command is used. The cache is automatically
               cleared when CMD changes between runs.
               If CMD was False or empty we do not check what command
               made the directory. Instead we simply read what is there.
        
        silent : boolean, optional
               If True no messages or warnings will be printed on screen.
        """
        if (not cmd):
            # no command, only read the cache directory
            cmd = ''
        elif (type(cmd)==list):
            self.cmd = " ".join(cmd)
        elif (type(cmd)==str):
            self.cmd = cmd
        else:
            raise RuntimeError("First argument CMD needs to be string or list")

        if (cache and not cmd):
            # No command, just read directory
            if not silent:
                print("Reading files in directory:",cache,"...")
            # Set cmd to the emty string, to signal not to run a command
            cmd = ''
        elif (cache and checkcmd(cache,self.cmd)):
            # Directory was created by the exact same command - just read
            if not silent:
                print("Using result cache in directory:",cache,"...")
            # Set cmd to the emty string, to signal not to run a command
            cmd = ''
        else:
            # Convert command string into list if necessary
            if (isinstance(cmd, str)):
                cmd = cmd.split()

                if cmd[0].startswith("~"):
                    cmd[0] = os.path.expanduser(cmd[0])

            # Find the optool executable
            bin = find_executable(cmd[0])
            if (not bin):
                raise RuntimeError("Executable not found: "+cmd[0])

        # Wrap the main part into try - finally to make sure we clean up
        try:
            if (cache):
                dir = cache
            else:
                # create temporary directory in /tmp/
                dir = tempfile.mkdtemp(prefix="optool_")
            if cmd:
                if cache:
                    # make sure directory is new and empty
                    shutil.rmtree(dir,ignore_errors=True)
                    os.mkdir(dir)
                # Store the command line we are using.  We store the
                # string version of the command, not the list version.
                writecmd(dir,self.cmd)
                # tell optool to use the directory as writing desination
                cmd.append('-o'); cmd.append(dir)
    
                # Run optool to produce the opacities
                stdout = subprocess.DEVNULL if silent else None
                stderr = subprocess.DEVNULL if silent else None
                cmd[0] = bin; subprocess.Popen(cmd, stdout=stdout, stderr=stderr).wait()
            
            # Check if there is output we can use
            scat,ext,translate = check_for_output(dir)
            self.scat = scat
            self.massscale = 1.

            kabs=[]; ksca=[]; kext=[]; gg=[]
            f11=[]; f12=[]; f22=[]; f33=[]; f34=[]; f44=[]
            nfiles=0; header=[];
            materials = []
            rho = []
            
            for i in range(5000):
                if scat:
                    file = ("%s/dustkapscatmat_%03d.%s") % (dir,(i+1),ext)
                else:
                    file = ("%s/dustkappa_%03d.%s") % (dir,(i+1),ext)
                file = translate.get(file,file)
                if (not os.path.exists(file)): break
                nfiles = nfiles+1
                x = readoutputfile(file,scat,silent=silent)
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
            self.np = nfiles
        finally:
            if cache:
                if not silent:
                    print("Files remain available in directory: "+dir)
            else:
                if not silent:
                    print("Cleaning up temporary directory "+dir)
                shutil.rmtree(dir)

    def plot(self,minkap=1e0):
        """Create interactive plots of the opacities in SELF.

        Furthermore, a plot for the scattering matric elements and, if the
        computemean() method has been called, a plot of the mean opacities
        are produces as well.
        """

        # Check if mean opacities have been computed
        if hasattr(self, 'kplanck'):
            # llamfmt = np.round(np.log10(self.lam),decimals=3)
            kplanck = self.kplanck
            kross   = self.kross
            temp = self.temp
            viewarr([kplanck,kross],index=1,ylabel=['kplanck','kross'],
                    idxnames=['grain index','log lambda [um]'],
                    idxvals=[np.array(range(self.np))+1,temp])

        # Extract the kappas and g
        kabs   = np.copy(self.kabs)
        ksca   = np.copy(self.ksca)
        kext   = kabs+ksca
        gg     = np.copy(self.gsca)
        maxkap = np.amax(kabs)
        if (maxkap<minkap*100):
            print('WARNING: you may want to change minkap to ',maxkap/100,' or smaller')
    
        # limit the kappa plotting range
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
                    idxvals=[np.array(range(self.np))+1,lamfmt,angfmt])

        # interactive plot of kabs, ksca, kext, and g
        llamfmt = np.round(np.log10(self.lam),decimals=3)
        viewarr([ggscal,kext,ksca,kabs],index=1,
                ylabel=['gg','kext','ksca','kabs'],
                idxnames=['grain index','log lambda [um]'],
                idxvals=[np.array(range(self.np))+1,llamfmt])


    def plotpi(self,ymin=None,ymax=None):
        """Create interactive plots of I,P,p.
        I is the intensity, f11
        P it the polarization, abs(f12)
        p is the degree of polarization, P/I

        The y-axis scaling is based on the values of P and p, but ignores I since
        it can become so large. ymax and ymin are keyword parameters to change
        the y-axis scaling.    
        """

        # Extract and plot the scattering matrix elements
        if self.scat:
            bottom = 1e-2
            int_I = self.f11
            int_P = np.sqrt(self.f12**2)  # abs would be good enough, but to remind of general case
            deg_P = abs(int_P/int_I)
            if (not ymin):
                ymin = min(np.min(np.matrix.flatten(int_P)),np.min(np.matrix.flatten(deg_P)))
            if (not ymax):
                ymax = max(np.max(np.matrix.flatten(int_P)),np.max(np.matrix.flatten(deg_P)))
    
            # Make version of grid variables with fewer digits
            lamfmt  = np.round(self.lam,decimals=3)
            angfmt  = np.round(self.scatang,decimals=3)

            # interactive plot of the scattering matric elements
            viewarr([int_I,int_P,deg_P],
                    index=2,ylabel=['f11','-f12',
                                    'f12/f11'],
                    idxnames=['grain index','lambda [um]','angle'],
                    idxvals=[np.array(range(self.np))+1,lamfmt,angfmt],
                    ymin=ymin,ymax=ymax)

    def select(self,i):
        """Select just one bin from a multi-particle object.
        A multi-particle opject is produced when running optool with
        a -d switch.

        This is useful for doing particle arithmetic, which only works for
        single particle objects.
        """
        x = copy.deepcopy(self)

        x.np = 1
        j = i+1
        
        x.fmax    = x.fmax[i:j]
        x.pcore   = x.pcore[i:j]
        x.pmantle = x.pmantle[i:j]

        x.amin    = x.amin[i:j]
        x.amax    = x.amax[i:j]
        x.nsub    = x.nsub[i:j]
        x.apow    = x.apow[i:j]
        x.amean   = x.amean[i:j]
        x.asig    = x.asig[i:j]
        x.a1      = x.a1[i:j]
        x.a2      = x.a2[i:j]
        x.a3      = x.a3[i:j]
        x.rho     = x.rho[i:j]
        x.chop    = x.chop[i:j]
        
        x.kabs    = x.kabs[i:j,:]
        x.ksca    = x.ksca[i:j,:]
        x.kext    = x.kext[i:j,:]
        x.gsca    = x.gsca[i:j,:]

        if x.scat:
            x.f11     = x.f11[i:j,:,:]
            x.f12     = x.f12[i:j,:,:]
            x.f22     = x.f22[i:j,:,:]
            x.f33     = x.f33[i:j,:,:]
            x.f34     = x.f34[i:j,:,:]
            x.f44     = x.f44[i:j,:,:]

        if (hasattr(x,'kross')):
            x.kplanck = x.kplanck[i:j,:]
            x.kross   = x.kross[i:j,:]

        return x

    def sizedist(self,N_of_a):
        """Compute opacity of a size distribution of elements of SELF.

        Arguments
        ---------
                
        N_of_a : numpy array containing the sumber of partiles of each size
                 available in SELF (as given by self.a1)
        """
        # Check if N_of_a is compatible with self.a1
        if (len(N_of_a) != len(self.a1)):
            raise RuntimeError('N_of_a and a1 arrays differ in length')
            
        # create a particle object to return
        x = copy.deepcopy(self)

        # Fill all attributes that make sense
        x.np = 1

        x.cmd     = ''

        x.materials = self.materials[0:1]
        
        x.fmax    = x.fmax[0:1]
        x.pcore   = x.pcore[0:1]
        x.pmantle = x.pmantle[0:1]

        x.amin    = x.a1[0:1]
        x.amax    = x.a1[-1:]
        x.nsub    = self.nsub[0]*self.np
        x.apow    = x.apow[0:1]
        x.amean   = x.amean[0:1]
        x.asig    = x.asig[0:1]
        x.a1 = x.a2 = x.a3 = -1;
        x.rho     = x.rho[0:1]
        x.chop    = x.chop[0:1]

        # Turn N_of_a into mass fractions, normalized to 1
        mass   = (4./3.) * np.pi * (self.a1*1e-4)**3 * self.rho
        m_of_a = N_of_a*mass
        mtot = np.sum(m_of_a)
        mfrac  = m_of_a/mtot
        x.massscale = 1

        # add up the opacities
        x.kabs = np.sum(self.kabs*mfrac[:,None],axis=0)
        x.ksca = np.sum(self.ksca*mfrac[:,None],axis=0)
        x.kabs = x.kabs[None,:]; x.ksca = x.ksca[None,:] # add particle size axis
        x.kext = x.kabs+x.ksca

        # compute gsca
        x.gsca = np.sum(self.ksca*self.gsca*mfrac[:,None],axis=0) / x.ksca[0]
        x.gsca = x.gsca[None,:]  # add particle size axis

        # compute the scattering matrix elements
        if x.scat:
            if self.norm == 'hovenier':
                w  = (self.ksca*mfrac[:,None])[:,:,None]
                wn = x.ksca[0,:,None]
                x.f11 = np.sum(self.f11*w,axis=0)/wn
                x.f12 = np.sum(self.f12*w,axis=0)/wn
                x.f22 = np.sum(self.f22*w,axis=0)/wn
                x.f33 = np.sum(self.f33*w,axis=0)/wn
                x.f34 = np.sum(self.f34*w,axis=0)/wn
                x.f44 = np.sum(self.f44*w,axis=0)/wn
            else:
                w  = mfrac[:,None,None]
                x.f11 = np.sum(self.f11*w,axis=0)
                x.f12 = np.sum(self.f12*w,axis=0)
                x.f22 = np.sum(self.f22*w,axis=0)
                x.f33 = np.sum(self.f33*w,axis=0)
                x.f34 = np.sum(self.f34*w,axis=0)
                x.f44 = np.sum(self.f44*w,axis=0)
            # Add the particle size axis
            x.f11 = x.f11[None,:]; x.f12 = x.f12[None,:]; x.f22 = x.f22[None,:]; 
            x.f33 = x.f33[None,:]; x.f34 = x.f34[None,:]; x.f44 = x.f44[None,:]; 

        # Return the new object
        return x
    
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
            return -1
        
        ang   = self.scatang
        lam   = self.lam
        wav   = 2.*np.pi/(lam*1e-4)      # need cm here, not micrometer
        ratio = np.zeros([self.np,self.nlam])

        # Compute values and weights for the integration
        if (self.gridtype == "boundary"):
            # Matrix values are on cell boundaries
            if (ang[0] != 0):
                raise RuntimeError("Inconsistency between gridtype \"boundary\" and angle values")
            thetab = ang*np.pi/180.
            mub = np.cos(thetab)
            dmu = mub[:-1]-mub[1:]   # Defined negatively for mu integral
            fc = 0.5*(self.f11[:,:,1:]+self.f11[:,:,:-1]) 
        else:
            # This is the standard grid with values on cell midpoints
            if (ang[0] == 0):
                raise RuntimeError("Inconsistency between gridtype \"center\" and angle values")
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

        Parameters
        ----------

        tmin : float
             minimum temperature for which to compute mean opacities
        tmax : float
             maximum temperature for which to compute mean opacities
        ntemp : int
             number of temperature steps between tmin and tmax
        """
        self.tmin    = tmin
        self.tmax    = tmax
        self.ntemp   = ntemp
        self.temp    = np.logspace(np.log10(tmin),np.log10(tmax),ntemp)
        self.kross   = np.zeros([self.np,self.ntemp])
        self.kplanck = np.zeros([self.np,self.ntemp])

        cl = 2.99792458e10          # Speed of light [cgs]
        nu = 1e4*cl/self.lam        # 10^4 because lam is in um - we need cm
        dnu = -1. * np.hstack([nu[1]-nu[0],0.5 * (nu[2:]-nu[:-2]), nu[-1] - nu[-2]  ])

        for it in range(self.ntemp):
            bnu    = bplanck(self.temp[it],nu)
            bnudt  = bplanckdt(self.temp[it],nu)
            dumbnu = np.sum(bnu*dnu)
            dumdb  = np.sum(bnudt*dnu)
            for ip in range(self.np):
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
            raise TypeError('Cannot add multi-particle objects')
        if ((s.nlam != o.nlam) or (np.abs((s.lam-o.lam)/s.lam).any()>1e-4)):
            raise RuntimeError('Wavelength grids differ')
        if (s.scat):
            if ((s.nang != o.nang) or
                (np.abs((s.scatang[1:]-o.scatang[1:])/s.scatang[1:]).any()>1e-4)):
                # We don't check the first value, could be 0
                raise RuntimeError('Angular grids differ')
            if (s.norm != o.norm):
                raise RuntimeError('Scattering normalizations differ')
        #
        # Now do the adding
        #
        x = copy.deepcopy(s)
        x.kabs = x.kabs+o.kabs
        x.ksca = x.ksca+o.ksca
        x.kext = x.kext+o.kext
        # F11 is linear in the integral for the computation of g.
        # So we can just take the weighted mean for g.
        x.gsca = (x.ksca*x.gsca + o.ksca*o.gsca) / (x.ksca+o.ksca)
        x.massscale = s.massscale + o.massscale

        if s.scat:
            # There is a scattering matrix.
            if s.norm == 'hovenier':
                # Add, weighted by kappa_scat
                ws = s.ksca[:,:,None]
                wo = o.ksca[:,:,None]
                wn = ws+wo
            else:
                # Just add the values
                ws, wo, wn = 1.,1.,1.
            x.f11 = (s.f11*ws + o.f11*wo) / wn
            x.f12 = (s.f12*ws + o.f12*wo) / wn
            x.f22 = (s.f22*ws + o.f22*wo) / wn
            x.f33 = (s.f33*ws + o.f33*wo) / wn
            x.f34 = (s.f34*ws + o.f34*wo) / wn
            x.f44 = (s.f44*ws + o.f44*wo) / wn
        #
        # Invalidate attributes that no longer make sense.
        #
        x.materials = np.hstack((x.materials,o.materials))
        if (x.fmax    != o.fmax   ): x.fmax    = -1
        if (x.pcore   != o.pcore  ): x.pcore   = -1
        if (x.pmantle != o.pmantle): x.pmantle = -1
        if (x.amin    != o.amin   ): x.amin    = -1
        if (x.amax    != o.amax   ): x.amax    = -1
        if (x.nsub    != o.nsub   ): x.nsub    = -1
        if (x.apow    != o.apow   ): x.apow    = -1
        if (x.amean   != o.amean  ): x.amean   = -1
        if (x.asig    != o.asig   ): x.asig    = -1
        if (x.rho     != o.rho    ): x.rho     = -1
        if (x.chop    != o.chop   ): x.chop    = -1
        x.a1,x.a2,x.a3 = -1,-1,-1

        if hasattr(s, 'kplanck'):
            kplanck = -1
            kross   = -1 
            temp    = -1

        return x
        
    def __mul__(s,o):
        """Multiplication for optool.particle objects.
        
        This is intended for the multiplication of such an object with
        a number.  The way to think about it is like this.  Such an
        contains opacities in units cm^2/g.  Multiplying it with a
        number means that the opacities are now per a different mass.
        This sounds strange, but it makes sense together with addition
        of particles - which see.
        """
        if (not (isinstance(o,int) or isinstance(o,float))):
            raise TypeError('optool.particle object can only be multiplied by a number')
        x = copy.deepcopy(s)
        x.kabs = x.kabs*o; x.ksca = x.ksca*o; x.kext = x.kext*o
        x.massscale = x.massscale*o
        if (s.scat and (s.norm != 'hovenier')):
            # We need to change the matrix as well, it's normalized to ksca
            x.f11 = x.f11*o; x.f12 = x.f12*o; x.f22 = x.f22*o
            x.f33 = x.f33*o; x.f34 = x.f34*o; x.f44 = x.f44*o
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

    def write(s,filename,header="Opacity file written by optool.particle.write"):
        """Write a single particle object to a file.
        
        The format of the file will be similar to the dustkappa.dat and
        dustkapscatmat.dat files produced by the optool FORTRAN program,
        with the difference that the header will not contain the detailed
        information about the computation.  But the file would be readable
        with the `readoutputfile' function.

        Arguments
        =========

        filename:  String, pointing the file name to which output should
                   be written.
        
        header:    A string that should be put at the beginning of the
                   file, as a commend describing the dataset.  The string
                   may have several lines, the # comment character will
                   automatically be added to the beginning of every line.
        """

        if (s.np>1):
            raise TypeError('Writing is not supported for multi-particle objects')
        try:
            wfile = open(filename, 'w')
        except:
            raise RuntimeError('Cannot write to file: '+filename)

        headerlines = header.splitlines()
        for i in range(len(headerlines)):
            wfile.write("# %s\n" % headerlines[i])
        if s.scat:
            wfile.write('  0\n')
            wfile.write('  %d\n' % s.nlam)
            wfile.write('  %d\n' % s.nang)
            wfile.write('\n')
        else:
            wfile.write('  3\n')
            wfile.write('  %d\n' % s.nlam)
            
        for i in range(s.nlam):
            # write the lambda grid and the opacities
            wfile.write(' %15.5e %15.5e %15.5e %15.5e\n' % (s.lam[i],s.kabs[0,i],s.ksca[0,i],s.gsca[0,i]))
            
        if s.scat:
            # we have a scattering matrix
            wfile.write('\n')
            # Write the angular grid
            for i in range(s.nang):
                wfile.write("%9.2f\n" % s.scatang[i])
            wfile.write('\n')
            # Write the scattering matrix
            for il in range(s.nlam):
                for ia in range(s.nang):
                    wfile.write('  %15.5e %15.5e %15.5e %15.5e %15.5e %15.5e\n' %
                                (s.f11[0,il,ia],s.f12[0,il,ia],s.f22[0,il,ia],
                                 s.f33[0,il,ia],s.f34[0,il,ia],s.f44[0,il,ia]))
        wfile.close()

class lnktable:
    """Class to work with lnk files.

lnk stands for lambda, n, and k, where and and k are the real and
imaginary components of the refractive index of a material.
    


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
    new.header = "# This is a silicate from Dorschner+1995"
    new.write('sil-Dorschner1995.lnk')

    """
    def __init__(self,file,i_lnk=[1,2,3],nskip=0,nlam_rho=True):
        """Create a new optool.lnktable object

        Parameters
        ----------
        
        file : str
             the file name from which to read the lnk data

        i_lnk : numpy array, optional
             the column numbers where to find lambda, the real part of the
             refractive index and the imaginary part of it, respectively.
             The default is [1,2,3] .

        nskip : int, optional
             Number of lines to skil at the beginning.  Lines starting with
            `#', `!' or `*` are stored as header lines and ar skipped in
             this way. So this parameter is for dealing with files that are
             not yet formatted in the standard way for optool.  The default
             is 0.

        nlam_rho : boolean, optional
             True means, the first unskipped line contains the number of
             wavelengths points and the specific density of the material.
             False means no such line exists, and the lines have to be
             counted.  Rho will be se to 0 then, to indicate that the value
             is not know at this point.
        """
        self.filename = file
        try:
            rfile = open(file, 'r')
        except:
            print('ERROR: File not found:',file)
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
            print("Warning: density rho is not known! Make sure to set it by hand.")

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
            raise RuntimeError('Cannot write to file: '+file)
        wfile.write(self.header)
        wfile.write("  %d  %g\n" % (self.nlam,self.rho))
        for i in range(self.nlam):
            wfile.write("  %16.6e %16.6e %16.6e\n" %
                        (self.lam[i],self.n[i],self.k[i]))
        wfile.close()
    def powerlaws(self):
        """Compute the extrapolation powerlaws."""
        print("n: ",(np.log(self.n[-1])-np.log(self.n[-2]))/(np.log(self.lam[-1])-np.log(self.lam[-2])))
        print("k: ",(np.log(self.k[-1])-np.log(self.k[-2]))/(np.log(self.lam[-1])-np.log(self.lam[-2])))


def logscale_with_sign(array,bottom):
    # Take the log10 of the absolute value of ARRAY, but transfer the
    # sign back onto the result.  Compress the region between
    # -BOTTOM and +BOTTOM into zero, smoothly.
    # This is a clever way to make a logarithmic plot of a variable
    # that has positive and negative values covering more than
    # one order of magnitude.
    lb = np.log10(bottom)
    a  =  np.where(array>0)
    b  =  np.where(array<=0)
    array[a] =  np.log10(array[a]+bottom)  - lb
    array[b] = -np.log10(-array[b]+bottom) + lb
    return array

def check_for_output(dir):
    # Check for and if necessary rename input files
    for ext in['dat','inp']:
        if (os.path.exists(dir+'/dustkapscatmat_001.'+ext)):
            return True, ext, dict()
        elif (os.path.exists(dir+'/dustkappa_001.'+ext)):
            return False, ext, dict()
        elif (os.path.exists(dir+'/dustkapscatmat.'+ext)):
            return True, ext, { dir+'/dustkapscatmat_001.'+ext : dir+'/dustkapscatmat.'+ext }
        elif (os.path.exists(dir+'/dustkappa.'+ext)):
            return False, ext, { dir+'/dustkappa_001.'+ext : dir+'/dustkappa.'+ext }
    raise RuntimeError('No valid OpTool output files found')

def parse_headers(headers,b):
    # Extract information on run parameters from headers
    n = len(headers)
    b.amin  = np.zeros(n); b.amax = np.zeros(n);
    b.apow  = np.zeros(n); b.amean = np.zeros(n); b.asig = np.zeros(n);
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
    # apow may or may not be present, so we need to test
    if m:
        b.apow[i]=float(m.group(1))
    # lgnm may or may not be present, so we need to test
    m = re.search(r" (lgnm|norm)\s*=\s*(-?[-0-9.eE]+):(-?[-0-9.eE]+)",headers[0])
    if m:
        b.amean[i]=float(m.group(2)); b.asig[i]=float(m.group(3))
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

def readoutputfile(file,scat,silent=False):
    """Read OpTool output file FILE.

    Parameters
    ----------

    file : str
         The file name to read
    scat : bool
         When True, the file contains a scattering matrix
    silent : boolean
         If True, no message or warning will be printed on screen

    Returns
    -------

    Depending on the SCAT flag, Returns a list with these elements
   
    [header,lam,kabs,ksca,phase_g]  or
    [header,lam,kabs,ksca,phase_g,scatang,f11,f12,f22,f33,f34,f44]
    """
    try:
        rfile = open(file, 'r')
    except:
        raise RuntimeError('File not found: '+file)
    if not silent:
        print('Reading',file,'...')

    # Read the header/comment field
    header = ''
    dum = rfile.readline()
    while dum.strip()[0]=='#':
        header = header + dum
        dum = rfile.readline()

    # Read the file format
    while len(dum.strip())<1: dum = rfile.readline() # skip any empty lines
    iformat = int(dum)

    # Read the number of wavelengths in the file and prepare arrays
    nlam = int(rfile.readline())
    lam=np.zeros(nlam); kabs=np.zeros(nlam); ksca=np.zeros(nlam); phase_g=np.zeros(nlam)

    if scat:
        # Read the scattering angular grid size and prepare arrays
        nang = int(rfile.readline())
        scatang = np.zeros(nang)
        f11=np.zeros([nlam,nang]); f12=np.zeros([nlam,nang]); f22=np.zeros([nlam,nang])
        f33=np.zeros([nlam,nang]); f34=np.zeros([nlam,nang]); f44=np.zeros([nlam,nang])

    # Read the opacities
    dum = rfile.readline()
    while len(dum.strip())<1: dum = rfile.readline() # skip any empty lines
    for ilam in range(nlam):
        dum           = dum.split()
        lam[ilam]     = float(dum[0])
        kabs[ilam]    = float(dum[1])
        ksca[ilam]    = float(dum[2])
        phase_g[ilam] = float(dum[3])
        dum = rfile.readline()

    if scat:
        # Read the angular grid
        while len(dum.strip())<1: dum = rfile.readline() # skip any empty lines
        for iang in range(nang):
            scatang[iang] = float(dum)
            dum = rfile.readline()

        # Read the scattering matrix
        while len(dum.strip())<1: dum = rfile.readline()
        dums = rfile.readlines()
        dums.insert(0,dum)
        data = np.fromstring("".join(dums),sep=' ')
        data = np.reshape(data,(nlam,nang,6),'C')
        f11[:,:]=data[:,:,0]; f12[:,:]=data[:,:,1]; f22[:,:]=data[:,:,2]
        f33[:,:]=data[:,:,3]; f34[:,:]=data[:,:,4]; f44[:,:]=data[:,:,5]

        rfile.close()
    if scat:
        return [header,lam,kabs,ksca,phase_g,scatang,f11,f12,f22,f33,f34,f44]
    else:
        return [header,lam,kabs,ksca,phase_g]

def writecmd(dir,cmd):
    """Store the CMD string in file DIR/cmd.
    """
    if (os.path.isdir(dir)):
        # Directory does not exist
        dir = dir.rstrip('/')
        filename = dir+"/cmd"
        try:
            wfile = open(filename, 'w')
        except:
            print('ERROR: Cannot write to file: ',filename)
            return False
        newcmd=cmd.strip()
        wfile.write(cmd+"\n")
        wfile.close()
        return True
    else:
        return False
    
def checkcmd(dir,cmd):
    """Check if new command line is the same as the old one.

    This functions checks if the directory DIR contains a file
    called CMD, and if the first line in thie directory is the
    same as the string passed with the DIR parameter.
    """
    if (not os.path.isdir(dir)):
        # Directory does not exist
        return False
    if (len(os.listdir(dir))<=1):
        # There are less than one file in the directory. So either the cmd
        # file does not exist, or no output files are present.
        return False
    dir = dir.rstrip('/')
    filename = dir+"/cmd"
    if (not os.path.exists(dir+"/cmd")):
        # The command file does not exist
        return False
    try:
        rfile = open(filename, 'r')
    except:
        print('ERROR: Cannot read file: ',filename)
        return False
    dum = rfile.readline()
    rfile.close()
    cached_cmd = dum.strip()
    new_cmd    = cmd.strip()
    return (cached_cmd == new_cmd)
    
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

p = particle  # Make p() an alias for the particle() method
