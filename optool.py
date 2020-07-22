import numpy as np
import matplotlib.pyplot as plt
import re
import os
import subprocess
from distutils.spawn import find_executable
import random

class particle:
    def __init__(self,cmd,keep=False):
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
    
            kabs=[]; ksca=[]; kext=[]; gg=[]
            f11=[]; f12=[]; f22=[]; f33=[]; f34=[]; f44=[]
            nfiles=0; header=[];
            
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
        finally:
            if keep:
                print("Keeping the temporary directory for inspection: "+tmpdir)
            else:
                print("Cleaning up temporary directory "+tmpdir)
                os.system('rm -rf '+tmpdir)

    def plot(self):
        # Create interactive plots of the opacities in SELF.
        
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
                    index=2,ylabel=['<1e-2','±1','','±1e2','','f11','f12','f22','f33','f34','f44'],
                    idxnames=['grain index','lambda [um]','angle'],
                    idxvals=[np.array(range(self.nsize))+1,lamfmt,angfmt])

        # interactive plot of kabs, ksca, kext, and g
        llamfmt = np.round(np.log10(self.lam),decimals=3)
        viewarr([ggscal,kext,ksca,kabs],index=1,ylabel=['gg','kext','ksca','kabs'],
                idxnames=['grain index','log lambda [um]'],
                idxvals=[np.array(range(self.nsize))+1,llamfmt])

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

    for i in range(n):
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
        b.materials.append([m.group(1),float(m.group(2)),float(m.group(3)),m.group(4)])
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
    else:
        b.radmc = False

    # FIXME: parse the composition as well

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

    
