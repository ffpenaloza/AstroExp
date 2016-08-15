from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.io import ascii
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os

# Se extraen las extensiones
f555w = fits.open('hst_10396_a4_acs_hrc_f555w_drz.fits')
f814w = fits.open('hst_10396_a4_acs_hrc_f814w_drz.fits')

f555w[1].writeto('sci_f555w_n121.fits',clobber=True)
f555w[2].writeto('invar_f555w_n121.fits',clobber=True)

f814w[1].writeto('sci_f814w_n121.fits',clobber=True)
f814w[2].writeto('invar_f814w_n121.fits',clobber=True)

f555w.close()
f814w.close()

# Corre SExtractor
os.system('sextractor sci_f555w_n121.fits -c f555w.sex')
os.system('sextractor sci_f814w_n121.fits -c f814w.sex')

#######################
# Se realiza el match #
#######################

# Se cargan listas con RA y DEC para cada imagen
RAf555w = np.loadtxt('f555w.cat',usecols=(3,))
DEf555w = np.loadtxt('f555w.cat',usecols=(4,))
RAf814w = np.loadtxt('f814w.cat',usecols=(3,))
DEf814w = np.loadtxt('f814w.cat',usecols=(4,))

# Match por parte de astropy
catalog = SkyCoord(ra=RAf555w*u.degree, dec=DEf555w*u.degree)  
c = SkyCoord(ra=RAf814w*u.degree, dec=DEf814w*u.degree)  
idx = c.match_to_catalog_sky(catalog,nthneighbor=1)

# Del catalogo f555w.cat se extraen las filas que indica el match
matches = list(idx[0])
f555w = np.loadtxt('f555w.cat')
f814w = np.loadtxt('f814w.cat')
out = []

j = 0
for i in matches:
	out.append(np.concatenate([f555w[i],f814w[j]]))
	j = j+1

# Salida a archivo
np.savetxt('n121_match.cat',out,
	fmt='%d\t%.4f\t%.4f\t%.7f\t%.7f\t%d\t%.4f\t%.4f\t%.7f\t%.7f',
	header='f555wN\tf555wMAG\tf555wMAGERR\tf555wALPHA\tf555wDELTA\tf814wN\tf814wMAG\tf814wMAGERR\tf814wALPHA\tf814wDELTA',
	footer='END')

# Genera plot
tbl = ascii.read('n121_match.cat')
plt.scatter((tbl["f555wMAG"]+0.24) - (tbl["f814wMAG"]+0.452) +  (3.180-1.803), (tbl["f555wMAG"]+0.24-3.180),c='black',s=1)
plt.xlim(0, 3)
plt.ylim(reversed(plt.ylim(14,22.5)))
plt.xlabel("$F555W-F814W$", fontsize=20)
plt.ylabel("$F555W$", fontsize=20)
plt.savefig('cmd_n121', dpi=None, facecolor='w', edgecolor='w',
        orientation='portrait', papertype=None, format=None,
        transparent=False, bbox_inches=None, pad_inches=0.1,
        frameon=None)

