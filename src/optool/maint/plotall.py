import copy
import numpy as np
import matplotlib.pyplot as plt
import math as m
import re
import os
import subprocess
from distutils.spawn import find_executable
import random
all_lnk_files = [
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
    'lnk_data/ol-c-mg100-Suto2006.lnk',
    'lnk_data/ol-c-mg95-Fabian2001.lnk',
    'lnk_data/ol-c-mg00-Fabian2001.lnk',
    'lnk_data/astrosil-Draine2003.lnk',

    'lnk_data/c-z-Zubko1996.lnk',
    'lnk_data/c-p-Preibisch1993.lnk',
    'lnk_data/c-gra-Draine2003.lnk',
    'lnk_data/c-org-Henning1996.lnk',
    'lnk_data/c-nano-Mutschke2004.lnk',

    'lnk_data/fe-c-Henning1996.lnk',
    'lnk_data/fes-Henning1996.lnk',
    'lnk_data/sic-Draine1993.lnk',
    
    'lnk_data/sio2-Kitamura2007.lnk',
    'lnk_data/cor-c-Koike1995.lnk',
    
    'lnk_data/h2o-w-Warren2008.lnk',
    'lnk_data/h2o-a-Hudgins1993.lnk',
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

def plotallk():
    import matplotlib.pyplot as plt
    import numpy as np
    import optool

    files = all_lnk_files

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
            #print(ix,iy,nn)
            if (nn >= len(files)):
                break
            file = files[nn]
            p=optool.lnktable(file)
            ax = gs1[ix,iy]
            ax.loglog(p.lam,p.k+1e-5)
            ax.set_xlim(0.05,300)
            ax.set_ylim(1e-4,1e3)
            ax.text(0.1,100.,file[9:-4],fontsize='xx-small')
    fig.show()
    fig.savefig("maint/all_k.pdf", bbox_inches='tight')

def plotalln():
    import matplotlib.pyplot as plt
    import numpy as np
    import optool

    files = all_lnk_files

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
            #print(ix,iy,nn)
            if (nn >= len(files)):
                break
            file = files[nn]
            p=optool.lnktable(file)
            ax = gs1[ix,iy]
            ax.semilogx(p.lam,p.n)
            ax.set_xlim(0.05,300)
            ax.set_ylim(0,10)
            ax.text(0.1,9.,file[9:-4],fontsize='xx-small')
    fig.show()

def extrapolateallk():
    import matplotlib.pyplot as plt
    import numpy as np
    import optool

    files = all_lnk_files

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
            #print(ix,iy,nn)
            if (nn >= len(files)):
                break
            file = files[nn]

            cmd = './optool -b -l 0.001 1e6 2000 '+file
            cmd = cmd.split()
            subprocess.Popen(cmd).wait()
            p=optool.lnktable('blended.lnk')
            ax = gs1[ix,iy]
            ax.loglog(p.lam,p.k+1e-5)
            ax.set_xlim(0.001,1e6)
            ax.set_ylim(1e-4,1e3)
            ax.text(0.1,100.,file[9:-4],fontsize='xx-small')
    fig.show()

def extrapolatealln():
    import matplotlib.pyplot as plt
    import numpy as np
    import optool

    files = all_lnk_files

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
            #print(ix,iy,nn)
            if (nn >= len(files)):
                break
            file = files[nn]
            cmd = './optool -b -l 0.001 1e6 2000 '+file
            cmd = cmd.split()
            subprocess.Popen(cmd).wait()
            p=optool.lnktable('blended.lnk')
            ax = gs1[ix,iy]
            ax.semilogx(p.lam,p.n)
            ax.set_xlim(0.001,1e6)
            ax.set_ylim(0,10)
            ax.text(0.1,9.,file[9:-4],fontsize='xx-small')
    fig.show()
