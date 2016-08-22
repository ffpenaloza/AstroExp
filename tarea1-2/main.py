from astropy.io import ascii
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from isochrones.dartmouth import Dartmouth_Isochrone
from isochrones import StarModel
import HessCMD

tbl = ascii.read('n121_match.cat')

def transform(v,i):
	c1f555 = [-0.09,-0.124]
	c2f555 = [0.034,0.018]
	c1f814 = [0.06,0.001]
	c2f814 = [-0.099,0.013]
	for j in range(8):
		tcol = v-i
		v = np.where(tcol<0.6, (tbl['f555wMAG']+c1f555[0]*tcol+c2f555[0]*tcol*tcol), (tbl['f555wMAG']+c1f555[1]*tcol+c2f555[1]*tcol*tcol))
		i = np.where(tcol<0.1, (tbl['f814wMAG']+c1f814[0]*tcol+c2f814[0]*tcol*tcol), (tbl['f814wMAG']+c1f814[1]*tcol+c2f814[1]*tcol*tcol))
	return v,i

tbl['V'], tbl['I'] = transform(tbl['f555wMAG'],tbl['f814wMAG'])

# Genera catalogo
data = [tbl['V'],tbl['I'],tbl['f814wALPHA'],tbl['f814wDELTA']]
ascii.write(data,'n121_ubvri.cat',delimiter='\t')

# Genera plots
HessCMD.plotHess(tbl['V'], tbl['V']-tbl['I'], cbarrtitle='Density', saveas='HessCMD_n121.png', cbarr='Yes')

fig, axis = plt.subplots(figsize=(8,6))
axis.plot(tbl['V']-tbl['I'],tbl['V'],'k.',alpha=0.4,ms=3)
axis.set_xlim(-0.5,2)
axis.set_ylim(16,25)
axis.invert_yaxis()
axis.set_xlabel('$V-I$',fontsize=20)
axis.set_ylabel('$V$',fontsize=20)
plt.savefig('cmd_n121_ubvri', dpi=300)

# Ajusta isocrona
iso = Dartmouth_Isochrone(bands=['V','I'])
model = iso.isochrone(age=10.021189299,feh=-1.56,distance=61000,AV=0.1272)
axis.plot(model.V_mag - model.I_mag,model.V_mag,'g',lw=2)
plt.savefig('iso_n121', dpi=300)

# Calcula ridge line
s = 0.22
x = []
y1 = np.arange(16.5,19.4,step=s)
y2 = np.arange(20,24.5,step=s)
y = np.append(y1,y2)
for i in y:
	a = np.where((i<tbl['V']) & (tbl['V']<i+s))
	x.append(np.median(tbl['V'][a]-tbl['I'][a]))

p = np.poly1d(np.polyfit(y,x,6))
#axis.plot(x,y,'r',lw=2)
axis.plot(p(y),y,'b',lw=2)
plt.savefig('rline_n121', dpi=300)

# Se busca iscocrona fit a mano
model2 = iso.isochrone(age=10.005,feh=-1.2,distance=61000,AV=0.1272)
axis.plot(model2.V_mag - model2.I_mag,model2.V_mag,'r',lw=2)
plt.savefig('iso_fit_n121', dpi=300)

# mueve el cmd a la derecha
plt.close()
fig, axis = plt.subplots(figsize=(8,6))
axis.set_xlim(-0.5,2)
axis.set_ylim(16,25)
axis.invert_yaxis()
axis.set_xlabel('$V-I$',fontsize=20)
axis.set_ylabel('$V$',fontsize=20)

axis.plot(tbl['V']-tbl['I']+0.1,tbl['V']+0.1,'k.',alpha=0.4,ms=3)
axis.plot(p(y)+0.1,y+0.1,'b',lw=2)
axis.plot(model2.V_mag - model2.I_mag,model2.V_mag,'r',lw=2)
plt.savefig('iso_fit2_n121', dpi=300)

