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
keys = []
for x in files:
    y = re.match(r'lnk_data/(.*)-.*.lnk',x)
    keys.append(y.group(1))
print(files)
print(keys)

fig, axs = plt.subplots(4, 4, sharex='col', sharey='row',
                        gridspec_kw={'hspace': 0, 'wspace': 0})
#(ax1, ax2, ax3, ax4),(ax5, ax6, ax7, ax8),(ax9, ax10, ax11, ax12),(ax13, ax14, ax15, ax16) = axs

for i in range(16):
    ax = axs.flat[i]
    data = read_lnk(files[i])
    lam = data[0]
    n = data[1]
    k = data[2]
    ax.plot(np.log10(lam),n)
    ax.plot(np.log10(lam),np.log10(k))
    ax.annotate(keys[i], (0.1, 0.5), textcoords='axes fraction', size=7)
    
fig.show()


    
