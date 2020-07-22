import glob
import re
import matplotlib.pyplot as plt
import numpy as np

def read_lnk(file):
    try:
        rfile = open(file, 'r')
    except:
        print('ERROR: file not found:',file)
        return -1
    print('Reading ',file,'...')
    dum = rfile.readline()
    while dum.strip()[0]=='#':
        dum = rfile.readline()
    dum = dum.split()
    nlam = int(dum[0])
    lam = np.zeros(nlam, dtype = np.float64)
    n = np.zeros(nlam, dtype = np.float64)
    k = np.zeros(nlam, dtype = np.float64)
    for ilam in range(nlam):
        dum = rfile.readline()
        dum = dum.split()
        lam[ilam]  = float(dum[0])
        n[ilam]    = float(dum[1])
        k[ilam]    = float(dum[2])
    rfile.close()
    return [lam,n,k]

files = np.array(glob.glob('lnk_data/*.lnk'))
files = np.sort(files)
#files = np.array(['lnk_data/c-org-Henning1900'])
nfiles = len(files)
keys = []
for x in files:
    y = re.match(r'lnk_data/(.*)-.*.lnk',x)
    keys.append(y.group(1))

if (nfiles == 1):
    fig, axs = plt.subplots(1, 1, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
elif (nfiles <=4):
    fig, axs = plt.subplots(2, 2, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
elif (nfiles <=6):
    fig, axs = plt.subplots(2, 3, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
elif (nfiles <=9):
    fig, axs = plt.subplots(3, 3, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
elif (nfiles <=12):
    fig, axs = plt.subplots(3, 4, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
elif (nfiles <=16):
    fig, axs = plt.subplots(4, 4, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
elif (nfiles <=25):
    fig, axs = plt.subplots(5, 5, sharex='col', sharey='row', gridspec_kw={'hspace': 0, 'wspace': 0})
else:
    print('too many lots')
    exit()
    
for i in range(len(files)):
    ax = axs.flat[i]
    data = read_lnk(files[i])
    lam = data[0]
    n = data[1]
    k = data[2]
    x = np.where(lam<1000)
    n = n[x]
    k = k[x]
    lam = lam[x]
    x=np.where(n<5)
    ax.plot(np.log10(lam[x]),n[x])
    x=np.where(k<5)
    ax.plot(np.log10(lam[x]),k[x])
    ax.plot([-2,3],[-0.5,-0.5])
    ax.annotate(keys[i], (0.1, 0.5), textcoords='axes fraction', size=7)
    
fig.show()


    
